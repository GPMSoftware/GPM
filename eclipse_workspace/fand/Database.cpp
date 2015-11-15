/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Database.cpp
 *
 *  Created on: Jun 29, 2010
 *      Author: gareth
 */

// TODO: We assume the loci are given in enum Locus. There are num_loci of them.
// Probably these should be parameters.

#include "Database.h"
#include "loci.h"
#include "Assert.h"
#include "util.h"

#include "cuda_accel/GPUDevices.h"

#include <UnitTest++/UnitTest++.h>

#include "fand/MessageStream.h"
INIT_MESSAGES("Database");
#include "fand/messages.h"

#include <my_global.h>
#include <my_sys.h>

using namespace std;

extern char *progname;

pthread_mutex_t init_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t conn_mutex = PTHREAD_MUTEX_INITIALIZER;


// separators within the blob
const string allele_sep   = ";";
const string locus_sep    = ",";

static unsigned int  port_num    = 0;       /* port number (0 means use built-in value) */
static char         *socket_name = NULL;    /* socket name (use built-in value) */
static char         *table_name  = "Profiles";  /* database name (default=none) */

static unsigned int  flags       = 0;       /* connection flags (none) */
//static bool          ask_password = false;  /* whether to solicit password */

Database::Database(std::string const &db_name, char *meta_fields_for_unit_tests)
: m_name(db_name)
, m_num_data_fields(10)
, m_host_name("localhost")
, m_port(port_num)
, m_user_name("")
, m_password("")
{
    // read environment variables
    m_host_name = getStringEnv("FAND_SERVER", m_host_name.c_str());
    m_port      = getIntEnv("FAND_PORT", m_port);
    m_user_name = getStringEnv("FAND_USER", m_user_name.c_str());
    m_password  = getStringEnv("FAND_PASS", m_password.c_str());

    // read the metadata fields
    bool meta_ok = false;

    if (meta_fields_for_unit_tests)
    {
        // in a unit test
        istringstream iss(meta_fields_for_unit_tests);
        meta_ok = readMetaFields(iss, m_meta_fields);
    }
    else
    {
        // in real life
        meta_ok = readMetaFields(m_meta_fields);
    }

    if (!meta_ok)
    {
        error << startl << "Error reading Metadata field definitions" << endl;
        m_meta_fields.clear();
    }

    m_num_fields = m_num_data_fields + m_meta_fields.size();
}

Database::~Database()
{
    disconnect();
}

bool
Database::clear()
{
    (void) erase();
    return create();
}

MYSQL *
Database::getConnection()
{
    int task = GPUDevices::getTask();

    MYSQL *ret = NULL;

    pthread_mutex_lock( &conn_mutex );

    if (m_conn.find(task) != m_conn.end())
    {
        ret = (MYSQL*)m_conn[task];
    }

    pthread_mutex_unlock( &conn_mutex );

    return ret;
}

void
Database::closeConnection()
{
    int task = GPUDevices::getTask();

    pthread_mutex_lock( &conn_mutex );

    Connection_map::iterator it = m_conn.find(task);
    if (it != m_conn.end() && it->second != 0)
    {
        mysql_close ((MYSQL*)it->second);

     // do not terminate client library since other threads may still use it
     // mysql_library_end();

        it->second = 0;
    }

    pthread_mutex_unlock( &conn_mutex );
}

void
Database::closeAllConnections()
{
    // ! must be called from the main thread when it is the only thread
    pthread_mutex_lock( &conn_mutex );

    Connection_map::iterator it;
    for (it = m_conn.begin(); it != m_conn.end(); ++it)
    {

        if (it->second != NULL)
        {
            mysql_close ((MYSQL*)it->second);

         // do not terminate client library since other threads may still use it
         // and we might want to re-connect.
         // mysql_library_end();

            it->second = NULL;
        }
    }
    pthread_mutex_unlock( &conn_mutex );
}

void
Database::setHost(string host_name)
{
    m_host_name = host_name;
}

void
Database::setPort(unsigned int port)
{
    m_port = port;
}

void
Database::setUser(string user_name)
{
    m_user_name = user_name;
}

void
Database::setPass(string password)
{
    m_password = password;
}

string
Database::getHost()
{
    return m_host_name;
}

unsigned int
Database::getPort()
{
    return m_port;
}

string
Database::getUser()
{
    return m_user_name;
}

string
Database::getPass()
{
    return m_password;
}

bool
Database::connect()
{
    // initialize client library (once only - in the main thread)
    static bool mysql_library_init_done = false;

    if ( ! mysql_library_init_done )
    {
        Assert2(isMainThread(), "Database::connect(): mysql_library_init should only be called from main thread");

        pthread_mutex_lock( &init_mutex );

        MY_INIT (progname);

        if (mysql_library_init (0, NULL, NULL))
        {
            print_error ("mysql_library_init() failed");
            return false;
        }
        mysql_library_init_done = true;

        pthread_mutex_unlock( &init_mutex );
    }

    // initialize connection handler
    MYSQL *conn = getConnection();
    if (conn != NULL)
    {
        // already connected
        return true;
    }

    conn = mysql_init (NULL);

    if (conn == NULL)
    {
        print_error ("mysql_init() failed (probably out of memory)");
        return false;
    }

//    m_host_name = "127.0.0.1";
    // connect to server
    if (mysql_real_connect (conn, m_host_name.c_str(), m_user_name.c_str(), m_password.c_str(),
            m_name.c_str(), m_port, socket_name, flags) == NULL)
    {
        print_error ("mysql_real_connect() failed");
        mysql_close (conn);
        return false;
    }

    // connection OK - remember it
    int task = GPUDevices::getTask();
    pthread_mutex_lock( &conn_mutex );
    m_conn[task] = (MYSQL_ptr)conn;
    pthread_mutex_unlock( &conn_mutex );

    // If "Profiles" table does not exist in database, create it
    return create();
}

bool
Database::disconnect()
{
    // disconnect from server
    closeConnection();
    return true;
}

bool
Database::disconnectAll()
{
    closeAllConnections();
    return true;
}

bool
Database::insert(DBProfile const &p)
{
    return write("INSERT", p);
}

bool
Database::overwrite(DBProfile const &p)
{
    return write("REPLACE", p);
}

bool isValid(DBProfile const &p)
{
    if (p.m_data.m_id.size() == 0 ||
        p.m_data.m_sample_id.size() == 0 ||
        p.m_data.m_profile_id.size() == 0 ||
        p.m_data.m_dataset.size() == 0)
    {
        return false;
    }
    return true;
}

void
Database::outputDBProfile(MYSQL_ROW const &row, int cols, MYSQL_FIELD *fields, DBProfile &p)
{
    // store data fields
    p.m_data.m_id               = row[0]==0 ? "" : row[0];
    p.m_data.m_dataset          = row[1]==0 ? "" : row[1];
    p.m_data.m_sample_id        = row[2]==0 ? "" : row[2];
    p.m_data.m_profile_id       = row[3]==0 ? "" : row[3];
    p.m_data.m_kit_type         = (ProfileType)(row[4]==0 ? prof_type_unknown : atoi(row[4]));
    p.m_data.m_evidence_type    = (EvidenceType)(row[5]==0 ? ev_type_unknown : atoi(row[5]));
    // date_added = row[6]
    p.m_data.m_num_contributors = row[7]==0 ? 1 : atoi(row[7]);
    p.m_data.m_error_rate       = row[8]==0 ? 0 : atof(row[8]);
    // profile_blob = row[9]

    // parse the blob
    string profile_blob = row[9]==0 ? "" : row[9];
    vector<string> locus_blob = split2(profile_blob, locus_sep);

    Assert(locus_blob.size() <= num_loci);

    for (size_t i=0; i< locus_blob.size(); ++i)
    {
        p.m_b11text[Locus(i)] = split(locus_blob[i], allele_sep);
    }

    // store metadata fields
    map<string, DBField>::const_iterator it;

    // step through the remaining fields, which should contain metadata, but in unknown order
    for (int i = m_num_data_fields; i < cols; ++i)
    {
        char *field_name = fields[i].name;

        // check it exists in the configuration file
        // TODO use a map for efficiency
        for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
        {
            if (it->second.field_name == field_name)
            {
                p.m_metadata[it->second.field_name] = row[i]==0 ? "" : row[i];
            }
        }
    }
}

std::string
Database::profileQuery(std::string const &profile_key)
{
    return string("SELECT * FROM ") + table_name + " WHERE profile_key = '" + profile_key + "'";
}

bool
Database::readProfile(std::string const &profile_key, DBProfile &p)
{
    // read the profile
    string sql_stmt = profileQuery(profile_key);

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    MYSQL *conn = getConnection();
    if (conn == NULL)
    {
        print_error ("readProfile(): no connection");
        return false;
    }

    MYSQL_RES *res_set = mysql_store_result(conn);

    if (res_set == NULL)
    {
        print_error ("readProfile(): mysql_store_result() failed");
        return false;
    }

    size_t ncols = mysql_num_fields(res_set);
    if (ncols != m_num_fields)
    {
        warn << startl << "expected " << m_num_fields << " columns but found " << ncols << endl;
        warn << alignl << "Metadata definitions file may be missing or out of date" << endl;
    }

    MYSQL_FIELD *fields;
    if ((fields = mysql_fetch_fields(res_set)) == NULL)
    {
        print_error("readProfile(): mysql_fetch_fields() failed");
        return false;
    }

    MYSQL_ROW row = mysql_fetch_row(res_set);
    if (row == NULL || mysql_errno(conn) != 0)
    {
        // this may not be an error - we may be just checking if the Profile already exists
//        print_error("read(): mysql_fetch_row() failed");
        return false;
    }

    outputDBProfile(row, ncols, fields, p);

    mysql_free_result (res_set);

    if (!isValid(p))
    {
        error << startl << "Profile read from database is invalid" << endl;
        return false;
    }

    return true;
}

std::string
Database::allProfilesQuery()
{
    return string("SELECT * FROM ") + table_name;
}

bool
Database::readAll(std::vector<DBProfile> &dbpvec)
{
    // read the dataset
    string sql_stmt = allProfilesQuery();

    return readQuery(sql_stmt, dbpvec);
}

std::string
Database::datasetQuery(std::string const &dataset_name)
{
    return string("SELECT * FROM ") + table_name + " WHERE dataset = '" + dataset_name + "'";
}

bool
Database::readDataset(std::string const &dataset_name, std::vector<DBProfile> &dbpvec)
{
    // read the dataset
    string sql_stmt = datasetQuery(dataset_name);

    return readQuery(sql_stmt, dbpvec);
}

// NB the query must return entire profiles, i.e. "SELECT * FROM Profiles WHERE ..."
bool
Database::readQuery(std::string const &sql_query, std::vector<DBProfile> &dbpvec)
{
//    cout << "Database::readQuery(): " << sql_query << endl;

    if (!do_sql(sql_query))
    {
        return false;
    }

    return outputDBProfiles(dbpvec);
}

// Query must be of the form "SELECT * FROM ..."
// We change this to "SELECT COUNT(*) FROM ...
int
Database::countQuery(std::string const &sql_query)
{
    std::string tmp = sql_query;
    to_upper(tmp);
    int pos = tmp.find("SELECT * FROM ");
    Assert(pos == 0);

    string count_query = "SELECT COUNT(*) FROM " + sql_query.substr(13);

    if (!do_sql(count_query))
    {
        return -1;
    }

    // read result
    vector<string> rc;
    if (! returnResultColumn(rc, 0))
    {
        return false;
    }

    return atoi(rc.begin()->c_str());
}

bool
Database::outputDBProfiles(std::vector<DBProfile> &dbpvec)
{
    MYSQL *conn = getConnection();
    if (conn == NULL)
    {
        print_error ("outputDBProfiles(): no connection");
        return false;
    }

    MYSQL_RES *res_set = mysql_store_result(conn);

    if (res_set == NULL)
    {
        print_error ("outputDBProfiles(): mysql_store_result() failed");
        return false;
    }

    size_t ncols = mysql_num_fields(res_set);
    if (ncols != m_num_fields)
    {
        warn << startl << "expected " << m_num_fields << " columns but found " << ncols << endl;
        warn << alignl << "Metadata definitions file may be missing or out of date" << endl;
    }

    MYSQL_FIELD *fields;
    if ((fields = mysql_fetch_fields(res_set)) == NULL)
    {
        print_error("readDataset(): mysql_fetch_fields() failed");
        return false;
    }

    MYSQL_ROW row;
    int count = 0;
    while ((row = mysql_fetch_row (res_set)) != NULL)
    {
        count++;

        if (mysql_errno(conn) != 0)
        {
            print_error("readDataset(): mysql_fetch_row() failed");
            return false;
        }

        ProfileData data;
        DBProfile dbprof(data);

        outputDBProfile(row, ncols, fields, dbprof);

        if (isValid(dbprof))
        {
            dbpvec.push_back(dbprof);
        }
        else
        {
            error << startl << "Profile read from database is invalid" << endl;
        }
    }

    if (count == 0)
    {
        // no records found
        return false;
    }

    mysql_free_result (res_set);

    return true;
}

bool
Database::listValues(string field, set<string> &values, int nmax)
{
    // read the column given by 'field'
    string sql_stmt = string("SELECT DISTINCT " + field + " FROM ") + table_name;

    if (nmax > 0)
    {
        std::ostringstream oss;
        oss << " ORDER BY " << field << " LIMIT " << nmax;
        sql_stmt += oss.str();
    }

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    vector<string> values_vec;
    if (!returnResultColumn(values_vec, 0)) return false;

    vec2set(values_vec, values);
    return true;
}

bool
Database::listDatasets(set<string> &datasets)
{
    return listValues("dataset", datasets);
}

// return rows affected by last operation
int
Database::rowsAffected()
{
    // call row_count
    string sql_stmt = string("SELECT ROW_COUNT()");

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    // read result
    vector<string> rc;
    if (! returnResultColumn(rc, 0))
    {
        return false;
    }

    return atoi(rc.begin()->c_str());
}

bool
Database::deleteProfile(std::string const &key)
{
    // read the dataset
    string sql_stmt = string("DELETE FROM ") + table_name + " WHERE profile_key = '" + key + "'";

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    // NB call to DELETE succeeds if record does not exist, but
    // rowsAffected will return 0
    return rowsAffected() > 0;
}

bool
Database::returnResultColumn(std::vector<std::string> &results, int col)
{
    results.clear();

    MYSQL *conn = getConnection();
    if (conn == NULL)
    {
        print_error ("returnResultColumn(): no connection");
        return false;
    }

    MYSQL_RES *res_set = mysql_store_result(conn);

    if (res_set == NULL)
    {
        print_error ("returnResultColumn(): mysql_store_result() failed");
        return false;
    }

    MYSQL_ROW row;
    while ((row = mysql_fetch_row (res_set)) != NULL)
    {
        if (mysql_errno(conn) != 0)
        {
            print_error("returnResultColumn(): mysql_fetch_row() failed");
            return false;
        }

        string s = row[col]==0 ? "" : row[col];
        results.push_back(s);
    }

    mysql_free_result (res_set);

    return true;
}

bool
Database::listProfileIDs(std::set<std::string> &profile_keys, std::string sql_stmt)
{
    if ( sql_stmt.empty())
    {
        // if no query given, select all profile IDs
        sql_stmt = string("SELECT profile_key FROM ") + table_name;
    }
    else
    {
        // if a query is given change the * to profile_key
        size_t pos = sql_stmt.find('*');
        Assert2(pos != sql_stmt.npos, "Database::listProfileIDs: Bad query");

        sql_stmt.replace(pos, 1, "profile_key");
    }

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    vector<string> profile_keys_vec;
    if (!returnResultColumn(profile_keys_vec, 0)) return false;

    vec2set(profile_keys_vec, profile_keys);
    return true;
}

bool
Database::moveProfiles(string dataset, string where_clause)
{
    // SQL: UPDATE Profiles SET dataset = whatever WHERE where_clause

    if ( where_clause.empty())
    {
        // nothing to do
        return true;
    }

    if ( dataset.empty())
    {
        return false;
    }

    string sql_stmt = string("UPDATE ") + table_name + " SET dataset = '" + dataset + "' WHERE " + where_clause;

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    return true;
}

bool
Database::deleteProfiles(std::string sql_query)
{
    // Delete all profiles matching sql_query.
    // NB query is of the form "SELECT * FROM profiles WHERE ..."
    // Edit this to read "DELETE FROM profiles WHERE ..."

    if ( sql_query.empty())
    {
        // nothing to do
        return true;
    }

    // if a query is given change the SELECT * to DELETE
    size_t pos = sql_query.find('*');
    Assert2(pos != sql_query.npos, "Database::deleteProfiles: Bad query");

    sql_query.replace(0, pos+1, "DELETE");

    if (!do_sql(sql_query))
    {
        return false;
    }

    return true;
}

void
Database::print_error (string const &message)
{
    MYSQL *conn = getConnection();

    error << startl << message << endl;
    if (conn != NULL)
    {
        error << alignl << "Error " << mysql_errno(conn) <<
                           " ("     << mysql_sqlstate(conn) <<
                           "): "    << mysql_error(conn) << endl;
    }
}

// quote a string for inclusion in SQL
string quote(const string &s)
{
    if (s.size() == 0 || s.find('\'') != string::npos)
    {
        return "NULL"; // unquoted, i.e. an SQL NULL
    }

    return "'" + s + "'";
}

// a comma separated list of metadata field names
string
Database::metaNames()
{
    string ret;

    map<string, DBField>::const_iterator it;

    for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
    {
        ret += ", ";
        ret += it->second.field_name; // field name
    }
    return ret;
}

// Field names and values formatted for an SQL CREATE statement
string
Database::createMetaFields()
{
    string ret;

    map<string, DBField>::const_iterator it;

    for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
    {
        ret += it->second.field_name + " " +   // field name
               it->second.field_type + ",\n";  // type
    }
    return ret;
}

// a comma separated list of the metadata values held in p
string
Database::metaValues(DBProfile const &p)
{
    string ret;

    map<string, DBField>::const_iterator it;

    for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
    {
        string field_name = it->second.field_name;

        ret += ", ";

        MetaData::const_iterator mit;
        if ((mit = p.m_metadata.find(field_name)) != p.m_metadata.end())
        {
            ret += "'" + mit->second + "'";
        }
        else
        {
            ret += "NULL";
        }
    }
    return ret;
}

bool
Database::write(char *op, DBProfile const &p)
{
    if (!isValid(p))
    {
        error << startl << "error: Attempt to write invalid Profile" << endl;
        return false;
    }

    string profile_key         = quote(p.m_data.m_id);
    string dataset             = quote(p.m_data.m_dataset);
    string sample_id           = quote(p.m_data.m_sample_id);
    string profile_id          = quote(p.m_data.m_profile_id);
    ProfileType kit_type       = p.m_data.m_kit_type;
    EvidenceType evidence_type = p.m_data.m_evidence_type;
    int num_donors             = p.m_data.m_num_contributors;
    float error_rate           = p.m_data.m_error_rate;

    std::ostringstream blob;

    for (int i=0; i<num_loci; ++i)
    {
        DBProfile::LocusMap::const_iterator it = p.m_b11text.find((Locus)i);
        if (it != p.m_b11text.end())
        {
//            Locus locus = it->first;
            const vector<string> text = it->second;

            for (size_t j=0; j<text.size(); ++j)
            {
                blob << text[j];

                if (j+1 < text.size())
                {
                    blob << allele_sep;
                }
            }
        }

        if (i+1 < num_loci)
        {
            blob << locus_sep;
        }
    }
    blob << flush;

    string blob_str = quote(blob.str());

//    info << startl << "blob_str.c_str() = [" << blob_str.c_str() << "]" << endl;

    // TODO: unsafe - can overflow. rewrite with string?
    char sql_stmt[16192];

    // NB date_added is a timestamp - it is given as NULL automatically
    // given the current tie by MySQL
    size_t n = snprintf(sql_stmt, sizeof(sql_stmt),
             "%s INTO %s ("
             "profile_key, dataset, sample_id, profile_id, kit_type, "
             "evidence_type, date_added, num_donors, error_rate, profile_blob%s) "
             "VALUES(%s, %s, %s, %s, '%d', '%d', NULL, '%d', '%f', %s%s)",
             op, table_name, metaNames().c_str(),
             profile_key.c_str(), dataset.c_str(), sample_id.c_str(), profile_id.c_str(), kit_type,
             evidence_type, /*date_added,*/ num_donors, error_rate, blob_str.c_str(), metaValues(p).c_str());

    if (n>=sizeof(sql_stmt))
    {
        error << startl << "sql_stmt: buffer overflow" << endl;
        return false;
    }

    return do_sql(sql_stmt);
}

bool
Database::do_sql(string const &sql_stmt)
{
    MYSQL *conn = getConnection();
    if (conn == NULL)
    {
        error << startl << "No database connection" << endl;
        return false;
    }

    if (mysql_query (conn, sql_stmt.c_str()) != 0)
    {
        print_error("mysql_query failed: " + sql_stmt);
        return false;
    }

    return true;
}

bool
Database::erase()
{
    string sql_stmt = string("DROP TABLE IF EXISTS ") + table_name;

    return do_sql(sql_stmt);
}

// create table, or update it if it already exists
bool
Database::create()
{
    bool ret = true;

    string sql_stmt = string("CREATE TABLE IF NOT EXISTS ") + table_name + "\n"
    "("
        "profile_key   VARCHAR(32),\n"
        "dataset       VARCHAR(32),\n"
        "sample_id     VARCHAR(16),\n"
        "profile_id    VARCHAR(15),\n"
        "kit_type      INT,\n"
        "evidence_type INT,\n"
        "date_added    TIMESTAMP,\n"
        "num_donors    INT,\n"
        "error_rate    FLOAT,\n"
        "profile_blob  VARCHAR(8096),\n"

            + createMetaFields() +

         "PRIMARY KEY   (profile_key)\n"
    ")";

    ret &= do_sql(sql_stmt);

    // now check if (if the table already existed) we need to add any metadata fields
    ret &= (checkMetaFields() >= 0);

    return ret;
}

long
Database::size()
{
    string sql_stmt = string("SELECT * FROM ") + table_name;

    return countQuery(sql_stmt);
}

// read the MetaData info from file
bool
Database::readMetaFields(
        std::map<std::string, DBField> & meta_fields)
{
    string meta_file = getStringEnv("META_DEFS");

    if (meta_file.size() == 0)
    {
        meta_file = homeDirStr() + "/.fand_meta.txt";
    }

    ifstream ifs(meta_file.c_str());

    if (! ifs.good())
    {
        error << startl << "Can't find Metadata definitions file: " << meta_file << endl;
        return false;
    }

    return readMetaFields(ifs, meta_fields);
}

// read the MetaData info from stream
bool
Database::readMetaFields(
        std::istream &is,
        std::map<std::string, DBField> & meta_fields)
{
    string column_header, field_name, field_type;

    string line;
    while (getline(is, line))
    {

        vector<string> words = split2(line, ",");
        cleanup(words);

        if (words.size() != 3)
        {
            error << startl << "Bad line in Metadata definitions file: " << line << endl;
            return false;
        }

        column_header = words[0];

        DBField field;
        field.field_name    = words[1];
        field.field_type    = words[2];

        meta_fields.insert(make_pair(column_header, field));

    // TESTING
//        cout << column_header << " " << field_name << " " << field_type << endl;
    }

    return !is.bad();
}


// list columns in table
bool
Database::listColumns(vector<string> &cols)
{
    // read the dataset
    string sql_stmt = string("SHOW columns IN ") + table_name;

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    return returnResultColumn(cols, 0);
}

bool
Database::listMetaColumns(set<string> &metacols)
{
    vector<string> cols;
    if (!listColumns(cols))
    {
        return false;
    }

    vec2set(cols.begin() + m_num_data_fields, cols.end(), metacols);

    return true;
}

// return the names of columns that correspond to ProfileData fields
bool
Database::listDataHeadings(set<string> &datacols)
{
    datacols.clear();

    datacols.insert("KIT TYPE");
    datacols.insert("CRIME/REF");
    datacols.insert("DATASET");
    datacols.insert("ERROR RATE");

    return true;
}

// check (and if necessary add) MetaData fields to database
int // number of fields added
Database::checkMetaFields()
{
    set<string> metacols;
    if (!listMetaColumns(metacols))
    {
        return -1;
    }

    int added = 0;
    string prev_field_name = "profile_blob";
    map<string, DBField>::const_iterator it;
    for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
    {
        string curr_field_name = it->second.field_name;
        if (metacols.find(curr_field_name) != metacols.end())
        {
            //  TODO: check field definition
//            if ()
//            {
//                // warning
//            }
        }
        else
        {
            // add
            if (!addMetaFieldAfter(it->second, prev_field_name))
            {
                return -1;
            }
            ++added;
        }

        prev_field_name = curr_field_name;
    }

    return added;
}

// add MetaField to database
bool
Database::addMetaFieldAfter(DBField const &mfield, string const &prev_field)
{
    string sql_stmt = string("ALTER TABLE ") + table_name
    		+ " ADD COLUMN " + mfield.field_name + " " + mfield.field_type
            + " AFTER " + prev_field;

    if (!do_sql(sql_stmt))
    {
        return false;
    }

    return true;
}

string
locusNameList()
{
    ostringstream oss;

    for (int i=0; i< num_loci; ++i)
    {
        oss << locus_name[i] << ",";
    }
    return oss.str();
}

void
Database::writeB11Header(std::ostream &ofs)
{
    ofs << "SAMPLE,PROFILE,DONORS,ALLELE,";

    ofs << locusNameList();

    ofs << ",KIT TYPE,CRIME/REF,DATASET,ERROR RATE,"; // NB: empty column between loci and metadata

    map<string, DBField>::const_iterator it;
    for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
    {
        ofs << it->first << ",";
    }

    ofs << endl;
}

void
Database::writeB11Profile(std::ostream &ofs, const DBProfile& dbp)
{
    ofs << dbp.m_data.m_sample_id << ","
        << dbp.m_data.m_profile_id << ","
        << dbp.m_data.m_num_contributors << ",";

    for (int n=0; n<2*dbp.m_data.m_num_contributors; ++n) // allele number
    {
        if (n>0)
        {
            ofs << ",,,";
        }
        ofs << n+1 << ",";

        // write loci
        for (int loc=0; loc<num_loci; ++loc)
        {
            vector<string> allele_strings;
            DBProfile::LocusMap::const_iterator it = dbp.m_b11text.find((Locus)loc);

            if (it != dbp.m_b11text.end())
            {
                allele_strings = it->second;
            }

            string data = (size_t)n<allele_strings.size() ? allele_strings[n] : "F";
            ofs << data << ",";
        }

        // write ProfileData at end of first line only
        ofs << ","; // empty column between loci and metadata

        if (n == 0)
        {
            ofs << kitName(dbp.m_data.m_kit_type) << ",";
            ofs << evidenceName(dbp.m_data.m_evidence_type) << ",";
            ofs << dbp.m_data.m_dataset << ",";
            ofs << dbp.m_data.m_error_rate << ",";
        }
        else
        {
            ofs << ",,,,";
        }

        // write Metadata
        map<string, DBField>::const_iterator it;
        for (it = m_meta_fields.begin(); it != m_meta_fields.end(); ++it)
        {

            // write metadata (if present) at end of first line only
            if (n == 0)
            {
                // if we have this field, output its value
                string db_field_name = it->second.field_name;
                MetaData::const_iterator it2 = dbp.m_metadata.find(db_field_name);
                if (it2 != dbp.m_metadata.end())
                {
                    ofs << it2->second;
                }
            }

            ofs << ",";
        }

        ofs << endl;
    }
}

// export in B11 format
bool
Database::exportDB(
    std::ostream &ofs)
{
    vector<DBProfile> dbp;

    if (!readAll(dbp))
    {
        error << startl << "export(): database read failed" << endl;
        return false;
    }

    writeB11Header(ofs);

    vector<DBProfile>::const_iterator it;
    for (it = dbp.begin(); it != dbp.end(); ++it)
    {
        writeB11Profile(ofs, *it);
    }

    return true;
}

bool
Database::import(
    bool overwrite_flag,
    std::istream &ifs,
    ProfileFilter::FileType ftype,
    ProfileType prof_type,
    string const &dataset,
    EvidenceType evidence_type,
    float error_rate,
    int sample_col,
    int profile_col)
{
    bool ret = true;

    std::vector<DBProfile>             dbpvec;

    // set up default ProfileData. These will be overridden by if present in the file
    ProfileData pd_def;
    pd_def.m_kit_type       = prof_type;
    pd_def.m_evidence_type  = evidence_type;
    pd_def.m_dataset        = dataset;
    pd_def.m_error_rate     = error_rate;
    dbpvec.push_back(DBProfile(pd_def));

    // Read DBProfiles
    if (ftype == ProfileFilter::B11)
    {
        std::vector<Profile> pvec; // not used
        if( ! ProfileFilter::theProfileFilter().readB11Profiles(pvec, /*index,*/ ifs, &dbpvec, this) )
        {
            error << startl << "readB11Profiles() failed" << endl;
            return false;
        }
    }
    else
    {
        // NB Identifiler, Powerplex16 and other variants are now all handled by readAnyProfiles.

        if (ftype == ProfileFilter::Identifiler)
        {
            sample_col = 0; profile_col = 1;
        }
        else if (ftype == ProfileFilter::PowerPlex16)
        {
            sample_col = 0; profile_col = -1;
        }

        std::vector<Profile> pvec; // not used
        if( ! ProfileFilter::theProfileFilter().readAnyProfiles(*this, ifs, &dbpvec, sample_col, profile_col) )
        {
            error << startl << "readAnyProfiles() failed" << endl;
            return false;
        }
    }

    int nrec = dbpvec.size();

    ok << "Writing " << nrec << " profiles to the database" << endl;

    int size_before = size();

    for (int i=0; i< nrec; ++i)
    {
        DBProfile &dbp = dbpvec[i];

        if (overwrite_flag)
        {
            if (!overwrite(dbp))
            {
                error << startl << "Overwrite Failed for " << dbp.m_data.m_id << " (continuing)" << endl;
                ret = false;
            }
        }
        else if (!insert(dbp))
        {
            error << startl << "Insert Failed for " << dbp.m_data.m_id << " (continuing)" << endl;
            ret = false;
        }

    }

    int size_after = size();

    ok << size_after - size_before << " new profiles added to the database" << endl;

    return ret;
}

static char meta_defns[] =
//      "WIS   wis   VARCHAR(16)\n"
//      "CEXE  cexe  VARCHAR(16)\n"
//      "LAB   lab   VARCHAR(16)\n"
        "TAT,  tat,  VARCHAR(16)\n"
        "DEF,  def,  VARCHAR(16)\n";

static char meta_defns2[] =
        "WIS,  wis,  VARCHAR(16)\n"
        "CEXE, cexe, VARCHAR(16)\n"
        "LAB,  lab,  VARCHAR(16)\n"
        "TAT,  tat,  VARCHAR(16)\n"
        "DEF,  def,  VARCHAR(16)\n";

void test_database(Database &db)
{
    CHECK(db.connect());
    CHECK(db.clear());

    vector<string> cols;
    CHECK(db.listColumns(cols));
    CHECK_EQUAL(10+2, cols.size());

    int i=0;
    CHECK_EQUAL(cols[i++], "profile_key");
    CHECK_EQUAL(cols[i++], "dataset");
    CHECK_EQUAL(cols[i++], "sample_id");
    CHECK_EQUAL(cols[i++], "profile_id");
    CHECK_EQUAL(cols[i++], "kit_type");
    CHECK_EQUAL(cols[i++], "evidence_type");
    CHECK_EQUAL(cols[i++], "date_added");
    CHECK_EQUAL(cols[i++], "num_donors");
    CHECK_EQUAL(cols[i++], "error_rate");
    CHECK_EQUAL(cols[i++], "profile_blob");
    // Metadata in alphabetic order
    CHECK_EQUAL(cols[i++], "def");
    CHECK_EQUAL(cols[i++], "tat");

    CHECK(db.disconnect());

    Database db2("test", meta_defns2);
    CHECK(db2.connect()); // Fields should be added to database

    CHECK(db2.listColumns(cols));
    CHECK_EQUAL(10+5, cols.size());

    // Metadata in alphabetic order
    i = 10;
    CHECK_EQUAL(cols[i++], "cexe");
    CHECK_EQUAL(cols[i++], "def");
    CHECK_EQUAL(cols[i++], "lab");
    CHECK_EQUAL(cols[i++], "tat");
    CHECK_EQUAL(cols[i++], "wis");

    CHECK(db2.disconnect());

    Database db3("test", meta_defns2);
    CHECK(db3.connect()); // No fields should be added this time

    CHECK(db3.listColumns(cols));
    CHECK_EQUAL(10+5, cols.size());
}

TEST (DATABASE_LOCAL)
{
    Database db("test", meta_defns);
    // by default we connect to a local database

    test_database(db);
}

//#define DB_REMOTE
// This test works only if there actually is a remote database
// So it is not part of the regular test-suite
#ifdef DB_REMOTE
TEST (DATABASE_REMOTE)
{
    Database db("test", meta_defns);
    db.setHost("fandserver");
    db.setUser("fanduser");
    db.setPass("fandpass");

    test_database(db);
}
#endif

struct Db_test
{
    Db_test()
    : db("test", meta_defns2)
    , id(Identifiler,
            "SAMPLEXX-PROFILEYY",    // ID
            1,                       // contributors
            0.001,                   // delta
            0,                       // mutation rate
            reference,               // evidence_type
            "BADDIES",               // dataset
            "SAMPLEXX",              // sample_id
            "PROFILEYY")             // profile_id
    , p(id), pB(id)
    {
        p.m_metadata["wis"]  = "SUSPECT22";
        p.m_metadata["cexe"] = "C-TEST";
        p.m_metadata["def"]  = "D-TEST";
        p.m_metadata["lab"]  = "L-TEST";
        p.m_metadata["tat"]  = "T-TEST";

        for (int i=0; i<num_loci; ++i)
        {
            p.m_b11text[Locus(i)].push_back("F");
            p.m_b11text[Locus(i)].push_back("11.1@.9");
        }

        pB = p;
        pB.m_data.m_id = "SAMPLEVV-PROFILEWW";
        pB.m_data.m_sample_id = "SAMPLEVV";
        pB.m_data.m_profile_id = "PROFILEWW";
        pB.m_data.m_dataset = "TESTXXX";
        pB.m_metadata["wis"]  = "SUSPECT33";
    }

    Database db;
    ProfileData id;
    DBProfile p;
    DBProfile pB;
};

#if 1

void checkBlob(DBProfile &p2)
{
    for (int i=0; i<num_loci; ++i)
    {
        DBProfile::LocusMap::const_iterator it = p2.m_b11text.find((Locus)i);

        if (it != p2.m_b11text.end())
        {
            if (p2.m_b11text[Locus(i)].size()==2)
            {
                CHECK_EQUAL("F", p2.m_b11text[Locus(i)][0]);
                CHECK_EQUAL("11.1@.9", p2.m_b11text[Locus(i)][1]);
            }
            else
            {
                CHECK(0);
            }
        }
        else
        {
            CHECK(0);
        }
    }
}

// check adding, reading and deleting profiles
TEST_FIXTURE(Db_test, DATABASE1)
{

    CHECK_EQUAL(true, db.connect());
    CHECK_EQUAL(true, db.clear());
    CHECK_EQUAL(0, db.size());
    CHECK_EQUAL(true, db.insert(p));
    CHECK_EQUAL(1, db.size());

    // insert again - should fail (duplicate entry)
    CHECK_EQUAL(false, db.insert(p));

    // overwrite - OK
    CHECK_EQUAL(true, db.overwrite(p));

    // read back into p2
    ProfileData id2;
    DBProfile p2(id2);
    CHECK_EQUAL(true, db.readProfile("SAMPLEXX-PROFILEYY", p2));

    // delete
    CHECK_EQUAL(false, db.deleteProfile("NOT-THERE")); // fails because not there
    CHECK_EQUAL(true, db.deleteProfile("SAMPLEXX-PROFILEYY"));
    CHECK_EQUAL(0, db.size());
    CHECK_EQUAL(false, db.deleteProfile("SAMPLEXX-PROFILEYY")); // fails because already gone

    CHECK_EQUAL("SAMPLEXX-PROFILEYY", p2.m_data.m_id);
    CHECK_EQUAL("BADDIES", p2.m_data.m_dataset);
    CHECK_EQUAL("SAMPLEXX", p2.m_data.m_sample_id);
    CHECK_EQUAL("PROFILEYY", p2.m_data.m_profile_id);
    CHECK_EQUAL((ProfileType)Identifiler, p2.m_data.m_kit_type);
    CHECK_EQUAL((EvidenceType)reference, p2.m_data.m_evidence_type);
    CHECK_EQUAL(1, p2.m_data.m_num_contributors);
    CHECK_CLOSE(0.001, p2.m_data.m_error_rate, 1e-6);
    CHECK_EQUAL("SUSPECT22", p2.m_metadata["wis"]);

    // check blob
    checkBlob(p2);

    // try to read a non-existent record
    CHECK_EQUAL(false, db.readProfile("", p2));
    CHECK_EQUAL(false, db.readProfile("SAMPLEXX-PROFILEZZ", p2));

}

// check exporting and importing profiles
TEST_FIXTURE(Db_test, DATABASE2)
{
    CHECK_EQUAL(true, db.connect());
    CHECK_EQUAL(true, db.clear());
    CHECK_EQUAL(0, db.size());
    CHECK_EQUAL(true, db.insert(p));
    CHECK_EQUAL(true, db.insert(pB));
    CHECK_EQUAL(2, db.size());

    ostringstream oss;
    CHECK(db.exportDB(oss));

    CHECK_EQUAL(true, db.clear());
    CHECK_EQUAL(0, db.size());

    istringstream iss(oss.str());

    CHECK(db.import(
        true,               // bool overwrite_flag
        iss,                // istream ifs
        ProfileFilter::B11, // ProfileFilter::FileType ftype
//        Identifiler,        // ProfileType prof_type         // Default value - may be overwritten
//        "BADDIES",          // string dataset                //
//        reference,          // EvidenceType evidence_type    //
//        0.001,              // float error_rate              //
        (ProfileType)0,     // ProfileType prof_type         // Default value - may be overwritten
        "",                 // string dataset                //
        (EvidenceType)0,    // EvidenceType evidence_type    //
        0,                  // float error_rate              //
        0,                  // int sample_col
        1));                // int profile_col

    CHECK_EQUAL(2, db.size());

    // read p into p2
    ProfileData id2;
    DBProfile p2(id2);
    CHECK_EQUAL(true, db.readProfile("SAMPLEXX-PROFILEYY", p2));

    // check we read it correctly
    CHECK_EQUAL("SAMPLEXX-PROFILEYY", p2.m_data.m_id);
    CHECK_EQUAL("BADDIES", p2.m_data.m_dataset);
    CHECK_EQUAL("SAMPLEXX", p2.m_data.m_sample_id);
    CHECK_EQUAL("PROFILEYY", p2.m_data.m_profile_id);
    CHECK_EQUAL((ProfileType)Identifiler, p2.m_data.m_kit_type);
    CHECK_EQUAL((EvidenceType)reference, p2.m_data.m_evidence_type);
    CHECK_EQUAL(1, p2.m_data.m_num_contributors);
    CHECK_CLOSE(0.001, p2.m_data.m_error_rate, 1e-6);

    // check metadata
    CHECK_EQUAL("SUSPECT22", p2.m_metadata["wis"]);
    CHECK_EQUAL("C-TEST",    p2.m_metadata["cexe"]);
    CHECK_EQUAL("L-TEST",    p2.m_metadata["lab"]);
    CHECK_EQUAL("T-TEST",    p2.m_metadata["tat"]);
    CHECK_EQUAL("D-TEST",    p2.m_metadata["def"]);

    // check blob
    checkBlob(p2);

    // read pB into p3
    ProfileData id3;
    DBProfile p3(id3);
    CHECK_EQUAL(true, db.readProfile("SAMPLEVV-PROFILEWW", p3));

    // check we read it correctly
    CHECK_EQUAL("SAMPLEVV-PROFILEWW", p3.m_data.m_id);
    CHECK_EQUAL("TESTXXX", p3.m_data.m_dataset);
    CHECK_EQUAL("SAMPLEVV", p3.m_data.m_sample_id);
    CHECK_EQUAL("PROFILEWW", p3.m_data.m_profile_id);
    CHECK_EQUAL((ProfileType)Identifiler, p3.m_data.m_kit_type);
    CHECK_EQUAL((EvidenceType)reference, p3.m_data.m_evidence_type);
    CHECK_EQUAL(1, p3.m_data.m_num_contributors);
    CHECK_CLOSE(0.001, p3.m_data.m_error_rate, 1e-6);

    // check metadata
    CHECK_EQUAL("SUSPECT33", p3.m_metadata["wis"]);
    CHECK_EQUAL("C-TEST",    p3.m_metadata["cexe"]);
    CHECK_EQUAL("L-TEST",    p3.m_metadata["lab"]);
    CHECK_EQUAL("T-TEST",    p3.m_metadata["tat"]);
    CHECK_EQUAL("D-TEST",    p3.m_metadata["def"]);

    // check blob
    checkBlob(p3);

}

#endif

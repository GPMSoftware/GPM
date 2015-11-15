/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Database.h
 *
 *  Created on: Jun 29, 2010
 *      Author: gareth
 *
 *      NB: The database named in the constructor must exist.
 *      Create a MySQL Database with e.g:
 *      mysql> CREATE DATABASE fand;
 *      mysql> GRANT ALL ON fand.* TO ''@'localhost'
 */

#ifndef DATABASE_H_
#define DATABASE_H_

#include "ProfileData.h"
#include "ProfileFilter.h"

#include <string>
#include <map>
#include <vector>
#include <set>
#include <mysql.h>

extern const std::string locus_sep;
extern const std::string allele_sep;

// db_field_name --> value;
typedef std::map<std::string, std::string> MetaData;

struct DBField
{
    DBField() {}
    DBField(std::string const &name, std::string const &type)
    : field_name(name), field_type(type) {}
    std::string field_name;
    std::string field_type;
};

// file_header --> (db_field_name, db_field_type);
typedef std::map<std::string, DBField> MetaField;

class DBProfile
{
public:
    DBProfile(ProfileData const &data) : m_data(data) {}

    ProfileData m_data;
    MetaData    m_metadata;

    typedef std::map<Locus, std::vector<std::string> > LocusMap;
    LocusMap m_b11text; // the 'blob' text for each locus (from profile_blob in the database)

    std::vector<std::string>&
    operator[](int i) { return m_b11text[(Locus)i]; }

    LocusMap::iterator
    find(int loc) { return m_b11text.find((Locus)loc); }

    int
    numLoci() { return m_b11text.size(); }
};

class Database
{
public:
    static const int SAMPLE_ID_LEN  = 16;
    static const int PROFILE_ID_LEN = 15;

    Database(std::string const &db_name, char *meta_fields_for_unit_tests = 0);
    virtual
    ~Database();

    bool connect();
    bool disconnect();
    bool isConnected() { return getConnection() != NULL; }
    bool disconnectAll();
    bool clear();
    bool insert(DBProfile const &p);
    bool overwrite(DBProfile const &p);
    bool readProfile(std::string const &profile_key, DBProfile &p);
    bool deleteProfile(std::string const &key);
    bool deleteProfiles(std::string sql_query);
    bool moveProfiles(std::string dataset, std::string where_clause);
    bool readDataset(std::string const &dataset_name, std::vector<DBProfile> &dbpvec);
    bool readQuery(std::string const &sql_query, std::vector<DBProfile> &dbpvec);
    int  countQuery(std::string const &sql_query);
    bool readAll(std::vector<DBProfile> &dbpvec);
    bool listDatasets(std::set<std::string> &datasets);
    bool listValues(std::string field, std::set<std::string> &values, int nmax = 0);
    bool listProfileIDs(std::set<std::string> &profile_keys, std::string sql_stmt = "");
    long size();
    bool listColumns(std::vector<std::string> &cols);      // list columns in table
    bool listMetaColumns(std::set<std::string> &metacols); // list metadata columns in table
    bool listMetaFields(MetaField &meta_fields) const      // list metadata fields
    { meta_fields = m_meta_fields; return true; }

    bool exportDB(std::ostream &ofs);
    bool import(
        bool overwrite_flag,
        std::istream &ifs,
        ProfileFilter::FileType ftype,
        ProfileType prof_type,
        std::string const &dataset,
        EvidenceType evidence_type,
        float error_rate,
        int sample_col,
        int profile_col);

    void setHost(std::string host_name);
    void setPort(unsigned int port);
    void setUser(std::string user_name);
    void setPass(std::string password);

    std::string getHost();
    unsigned int getPort();
    std::string getUser();
    std::string getPass();

    static std::string profileQuery(std::string const &profile_key);
    static std::string allProfilesQuery();
    static std::string datasetQuery(std::string const &dataset_name);

    // list ProfileData headings in files
    static bool listDataHeadings(std::set<std::string> &datacols);

private:
    Database(Database const &);
    Database &operator=(Database const &);

    std::string metaNames();                          // a comma separated list of metadata field names
    std::string createMetaFields();                   // Field names and values formatted for an SQL CREATE statement
    std::string metaValues(DBProfile const &p);       // a comma separated list of the metadata values held in p

    int  checkMetaFields();                           // check (and if necessary add) MetaData fields to database
    bool addMetaFieldAfter(DBField const &mfield,     // add MetaField to database
            std::string const &prev_field);

    void print_error (std::string const &message);
    bool do_sql(std::string const &sql_stmt);
    bool write(char *op, DBProfile const &p);
    bool erase();
    bool create();
    bool returnResultColumn(std::vector<std::string> &results, int col);
    void outputDBProfile(MYSQL_ROW const &row, int cols, MYSQL_FIELD *fields, DBProfile &p);
    bool outputDBProfiles(std::vector<DBProfile> &dbpvec);
    int  rowsAffected();

    static bool readMetaFields(            // read the MetaData info from file
            std::istream &is,
            MetaField & meta_fields);

    // read metadata field definitions from file
    static bool readMetaFields(
            MetaField & meta_fields);

    void writeB11Header(std::ostream &ofs);
    void writeB11Profile(std::ostream &ofs, const DBProfile& dbp);

    MYSQL *getConnection();   // get connection for this task (thread)
    void   closeConnection(); // close connection for this task (thread)
    void closeAllConnections(); // close all connections for the database

    std::string const                m_name;        // database name
    MetaField                        m_meta_fields; // definition of metadata fields

    typedef void* MYSQL_ptr;
    typedef std::map<int, MYSQL_ptr> Connection_map;
    Connection_map m_conn;                          // task (thread) -> connection handler

    size_t m_num_fields, m_num_data_fields;

    std::string m_host_name; // server host (default=localhost)
    unsigned int m_port;     // server port number
    std::string m_user_name; // username (default=login name)
    std::string m_password;  // password (default=none)

    friend class TestDATABASE0; // unit tests

//    friend class SelectProfilesPanel; // uses do_sql to define functions
    friend void addDBFunctions(Database & db);

};

#endif /* DATABASE_H_ */

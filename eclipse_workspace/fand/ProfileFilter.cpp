/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ProfileFilter.cpp
 *
 *  Created on: Mar 11, 2010
 *      Author: gareth
 */

#include "ProfileFilter.h"
#include "Parser.h"
#include "Database.h"
#include "loci.h"
#include "dbmatch.h"
#include "util.h"
#include "Assert.h"
#include "cuda_accel/cuda_accel.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <algorithm>
#include <cctype>
#include <UnitTest++/UnitTest++.h>

#include "MessageStream.h"
INIT_MESSAGES("ProfileFilter")
#include "messages.h"

using namespace std;

// certainty indicated by '(12)' notation
static double BRACKETS_FREQ  = BRACKETS_FREQ_DEFAULT;

// certainty indicated by '?(12)' notation
static double QBRACKETS_FREQ = QBRACKETS_FREQ_DEFAULT;

// B11 format
static const string id_sep       = "-";
static const string term_sep     = "/";
static const string prob_sep     = "@";
static const string unknown_char = "F";
static const string sample_char  = "D"; // "S" better ?

// Other formats
static const ProfileFilter::FileFormat identifiler_fmt = { ProfileFilter::Identifiler, 1, 16, 2, " " };
static const ProfileFilter::FileFormat powerplex16_fmt = { ProfileFilter::PowerPlex16, 0, 16, 3, "-" };

ProfileFilter::ProfileFilter()
{
    // do locus name comparisons in upper case
	for (int i=0; i<locus_name_size; ++i)
	{
	    string s(locus_name[i]);
	    to_upper(s);
		m_locus_name_map[s] = (Locus)locus_index[i];
	}

	for (int i=0; i<num_loci; ++i)
	{
        string s(locus_name[i]);
//        to_upper(s); // no
		m_locus_names[(Locus)locus_index[i]] = s;
	}

	// Read environment variables
	float f = getFloatEnv("BRACKETS_FREQ", BRACKETS_FREQ_DEFAULT);
	if (f>=0 && f <=1)
	{
		BRACKETS_FREQ = f;
	}

	f = getFloatEnv("QBRACKETS_FREQ", QBRACKETS_FREQ_DEFAULT);
	if (f>=0 && f <=1)
	{
		QBRACKETS_FREQ = f;
	}
}

ProfileFilter::~ProfileFilter() {
}

// only one object exists - and that is only needed to ensure the c'tor is called
ProfileFilter const &
ProfileFilter::theProfileFilter()
{
	static ProfileFilter pf;
	return pf;
}

Profile
makeProfile(DBProfile const &dbp, PopulationData *pop)
{
    // copy data fields
    Profile ret(dbp.m_data);

    // parse the blob
    DBProfile::LocusMap::const_iterator it;

    for (it = dbp.m_b11text.begin(); it != dbp.m_b11text.end(); ++it)
    {
        Locus loc = it->first;
        vector<string> allele_strings = it->second;

        if (!pop || pop->hasLocus(loc))
        {
            for (size_t i=0; i<allele_strings.size(); ++i)
            {
                // parse the data
                PMF<Allele> pmf;
                if (Parser::theParser().parse(allele_strings[i], pmf, pop, loc) == Parser::yacc_error)
                {
                    Assert2(false, "Bad blob in database");
                }

                // add allele to profile
                ret[loc].setAllele(i, pmf);
            }
        }
        else
        {
            // ignore locus. No error message - it is normal for some loci to be absent
        }
    }

    return ret;
}

ostream &
operator<<(ostream &os, const PMF<Allele> &pmf)
{
	if (pmf.empty())
	{
		os << unknown_char;
	}
	else
	{
		// To get the formatting right (such as field widths)
		// there should be a single write to os,
		// so accumulate the output in a string first.
		ostringstream oss;
		PMF<Allele>::const_iterator it = pmf.begin();

		Assert (it != pmf.end());

		while (true)
		{
			oss << it->first;
			if (it->second != 1)
			{
				char buf[10];
				sprintf(buf, "%.2g", it->second);
				char *p = buf;
				if (buf[0] == '0' && buf[1] == '.') ++p; // remove leading 0
				oss << prob_sep << p;
			}

			if (++it == pmf.end()) break;
			oss << term_sep;
		}
		os << oss.str();
	}

	return os;
}

#if 0
// output a Profile in human-readable form (on one line)
ostream &
operator<<(ostream &os, Profile const &p)
{
	os << setw(16) << p.data().m_name << ": ";

	for (int locus = 0; locus < num_loci; ++locus)
	{
		os << setw(10) << p[locus].getAllele(0);
		os << setw(10) << p[locus].getAllele(1);
	}
	return os;
}
#else
// output a Profile in human-readable form (on multiple lines)
ostream &
operator<<(ostream &os, Profile const &p)
{
	for (int i=0; i<p.data().m_num_contributors*2; ++i)
	{
		if (i==0)
		{
			os << setw(16) << p.data().m_id  << ": " << setw(6) << i+1;
		}
		else
		{
			os << setw(24)  << i+1;
		}

		for (int locus = 0; locus < num_loci; ++locus)
		{
			os << " " << setw(9) << p[locus].getAllele(i);
		}
		os << endl;
	}

	return os;
}
#endif

Allele
decodeAllele(string const &s)
{
	if (s == "X")
		return Allele::X;

	if (s == "Y")
		return Allele::Y;

	istringstream iss(s);
	float allele_name; // e.g. 16.2
	if (! (iss >> allele_name))
	{
		// doesn't look like a float
		error << startl << "not an allele: " << s << endl;
		throw std::exception();
	}

	allele_name = fabs(allele_name); // ignore '-' (should not happen)

	int repeats = int(allele_name);
	int variant = int( (allele_name - repeats + 0.01) * 10); // nasty!

    return Allele(repeats, variant);
}

float
decodeProb(string const &s)
{
	istringstream iss(s);
	float allele_prob; // e.g. 16.2
	if (! (iss >> allele_prob))
	{
		// doesn't look like a float
		error << startl << "not a float: " << s << endl;
		throw std::exception();
	}
	return allele_prob;
}

int
decodeInt(string const &s)
{
	istringstream iss(s);
	int ret;
	if (! (iss >> ret))
	{
		// doesn't look like an int (may not matter)
		return 0;
	}
	return ret;
}

PMF<Allele>
decodeKitAllelePMF(string const &s)
{
	PMF<Allele> ret;

	// if it is 'F' (unknown) return an empty PMF
	if (s == "F") return ret;

	double freq = 1;
	Allele a;

	if (s[0] == '(')
	{
		freq = BRACKETS_FREQ;
		string ss = s.substr(1, s.size()-2);
		a = decodeAllele(ss);
	}
	else if (s[0] == '?' && s[1] == '(')
	{
		freq = QBRACKETS_FREQ;
		string ss = s.substr(2, s.size()-3);
		a = decodeAllele(ss);
	}
	else if (s[0] == '-') // treat - as () (spreadsheet magic)
	{
		freq = BRACKETS_FREQ;
		string ss = s.substr(1);
		a = decodeAllele(ss);
	}
	else
	{
		a = decodeAllele(s);
	}

	ret[a] = freq;

	if (ret.sum() > 1)
	{
		ret.normalize();
		warn << startl << "WARNING: allele frequencies sum to more than one (normalizing)" << endl;
	}

	return ret;
}

// add to PMF each allele in sample at background frequencies, normalize
// and then scale them all by freq.
void
oneOfTheAbove(PMF<Allele> &pmf, float freq, set<Allele> const &sample, PMF<Allele> background)
{
	set<Allele>::const_iterator it;
	for (it = sample.begin(); it != sample.end(); ++it)
	{
		pmf[*it] = background[*it];
	}
	pmf.normalize();
	pmf *= freq;
}

// add each Allele in pmf to sample
void
addToSample(set<Allele> &sample, PMF<Allele> const &pmf)
{
	PMF<Allele>::const_iterator it;
	for (it = pmf.begin(); it != pmf.end(); ++it)
	{
		sample.insert(it->first);
	}
}

// check for unknown alleles and if possible add them to the population database
//
// There are a number of other things we could do here
// * build a frequency table of rare alleles
// * build a frequency table of all alleles (updating those in the population database)
//
int checkAlleles(PMF<Allele> &pmf, PopulationData &popdata, Locus loc, const string& locus_name, Unknowns unknown, int pop_size)
{
	int added = 0;

	PMF<Allele>::iterator ia, current;
	for(ia = pmf.begin(); ia != pmf.end();)
	{
		current = ia++;

        // check allele is in population database
        double f = popdata.getFrequency(loc, current->first, &added);
        if (f == 0)
        {
            // Ignore allele by removing from PMF. That is equivalent to replacing it with 'F'.
            warn << startl << "erasing loc = " << (Locus)loc << " allele = " <<  current->first << endl;
            pmf.erase(current);
            warn << alignl << "pmf.size() = " << pmf.size() << endl;
        }
	}

	return added;
}

// de-duplicate a duplicated identifier
void
deDup(std::map<std::string, int> const &index, string & ident)
{
	// check for duplicated identifier
	if (index.find(ident) != index.end())
	{
        error << startl << "duplicate Profile ID: rejecting record: " << ident << endl;
		throw std::exception();

// The commented out code creates a unique identifier

//		string alt_ident;
//		for (char c = 'A'; c <= 'Z'; c++)
//		{
//			string s = ident + '-' + c;
//			if (index.find(s) == index.end())
//			{
//				alt_ident = s;
//				break;
//			}
//		}
//
//		if (alt_ident == "")
//		{
//			error << startl << "duplicate Profile ID: rejecting record: " << ident<< endl;
//			throw std::exception();
//		}
//		else
//		{
//			warn << startl << "duplicate Profile ID: replacing " << ident << " with " << alt_ident << endl;
//			ident = alt_ident;
//		}
	}
}

// NB this only produces 'simple' text i.e. no 'D' or 'B'
void
pmfsToText(vector<string> &ret, vector< PMF<Allele> > const &pmfs)
{
    ret.clear();
    ostringstream tstream;

    for (size_t i=0; i<pmfs.size(); ++i)
    {
        tstream.str("");
        PMF<Allele>::const_iterator current, it = pmfs[i].begin();
        while (it != pmfs[i].end())
        {
            current = it++;
            tstream << current->first;
            if (current->second < 1)
            {
                tstream << "@" << current->second;
            }

            if (it != pmfs[i].end())
            {
                tstream << '/';
            }
        }

        string text = tstream.str();
        if (text.size() == 0)
        {
            text = "F";
        }
        ret.push_back(text);
    }
}

// columns for a particular locus:
// m_col[0] is the column labeled '1'
// m_col[1] is the column labeled '2'
// m_col[2] is the column labeled '3'
struct Cols
{
    Cols() {m_col[0] = m_col[1] = m_col[2] = -1;}

    int &
    operator[](int i)
    {
        Assert(0<=i && i<3);
        return m_col[i];
    }

private:
    int m_col[3];
};

typedef map<int, Cols> LocusCols;

bool
ProfileFilter::readHeaders(
        Database        const &db,
        std::istream          &ifs,   // IN:  input stream
        map<int, int> &locus_index,   // OUT: column -> Locus
        map<int, string> &metad_index,// OUT: column -> metadata field name
        map<int, string> &profd_index // OUT: column -> ProfileData field name
        ) const
{
    // Determine if a column header is Locus or Metadata,
    // and create maps (column -> Locus) and (column -> Metadata field name)
    locus_index.clear();
    metad_index.clear();
    profd_index.clear();

    MetaField meta_fields;
    db.listMetaFields(meta_fields);

    std::set<std::string> data_fields;
    db.listDataHeadings(data_fields);

    // keep track of columns associated with each allele number for each locus
    LocusCols cols;

    string line;
    getline(ifs, line);

    vector<string> headers = split2(line, locus_sep);
    cleanup(headers);

    PopulationData pop = knownLoci();

    int n; // column index
    vector<string>::const_iterator it;
    for(n = 0, it = headers.begin();
        it != headers.end();
        ++n, ++it)
    {

        // If column header ends in "-1", "-2" "-3" or " 1", " 2", " 3"
        // assume it is supposed to be a Locus

        int len = it->size();
        char y  = (*it)[len-2];
        char z  = (*it)[len-1];

        if ( (y == '-' || y == ' ' || y == '_') && (z == '1' || z == '2' || z == '3'))
        {
            string locus_name = it->substr(0, len-2);
            int allele_num = z - '1'; // (0, 1, 2) for ('1', '2', '3')

            string s = locus_name; to_upper(s);
            map<string, Locus>::const_iterator map_it = m_locus_name_map.find(s);
            int locus = map_it->second;

            if (map_it == m_locus_name_map.end())
            {
                error << startl << "unrecognized locus: " << locus_name << endl;
                return false;
            }

            locus_index[n] = locus;
            cols[locus][allele_num] = n;

            if (! pop.hasLocus(locus))
            {
                // give warning now so we do so only once
                warn << startl << "locus not in population database (will not be used): " << locus_name << endl;
            }
        }
        else if (data_fields.find(*it) != data_fields.end())
        {
            profd_index[n] = *it;
        }
        else if (meta_fields.find(*it) != meta_fields.end())
        {
            string field_name = meta_fields[*it].field_name;
            metad_index[n] = field_name;
        }
    }

    // check we have a least a '1' and a '2' for every locus
    LocusCols::iterator cit;
    for (cit = cols.begin(); cit != cols.end(); ++cit)
    {
        if (cit->second[0] < 0 || cit->second[1] < 0)
        {
            error << startl << "need alleles '1' and '2' for locus " << m_locus_names.at((Locus)cit->first) << endl;
            return false;
        }
    }

    return true;
}

void
setProfileData(
        ProfileData  &pd,          // OUT
        string const &field_name,  // IN
        string const &val)         // IN
{
    if (field_name == "KIT TYPE")
    {
        pd.m_kit_type = kitID(val);
    }
    else if (field_name == "CRIME/REF")
    {
        pd.m_evidence_type = evidenceID(val);
    }
    else if (field_name == "DATASET")
    {
        pd.m_dataset = val;
    }
    else if (field_name == "ERROR RATE")
    {
        pd.m_error_rate = atof(val.c_str());
    }
    else
    {
        // Programming error - must be one of these
        Assert2(false, "ProfileData field not found");
    }
}

// Read profiles in "database mode". Return the profiles in the dbpvec array
//
// The first column is assumed to be the identifier unless another column is indicated.
//
// Any column header not recognized as a locus is checked to see if it is metadata.
//
bool
ProfileFilter::readAnyProfiles(
        Database                     const &db,
        std::istream                       &ifs,          // IN
        std::vector<DBProfile>             *dbpvec,       // OUT DBProfiles read
        int sample_col,                                   // IN column containing sample ID
        int profile_col                                   // IN column containing profile ID
        ) const
{
    Assert(dbpvec != 0);

    map<int, int> locus_index;
    map<int, string> metad_index; // map from column order to metadata headers
    map<int, string> profd_index; // map from column order to profile data headers

    std::map<std::string, int> index;

    if (!readHeaders(db, ifs, locus_index, metad_index, profd_index))
    {
        return false;
    }

    info2 << endl << startl << "locus_index: " << endl;
    map<int, int>::const_iterator it2;
    for (it2 = locus_index.begin(); it2 != locus_index.end(); ++it2)
    {
        info2 << it2->first << " --> " << it2->second << ", ";
    }
    info2 << endl;

    info2 << endl << startl << "metad_index: " << endl;
    map<int, string>::const_iterator it3;
    for (it3 = metad_index.begin(); it3 != metad_index.end(); ++it3)
    {
        info2 << it3->first << " --> " << it3->second << ", ";
    }
    info2 << endl;

    info2 << endl << startl << "profd_index: " << endl;
    map<int, string>::const_iterator it4;
    for (it4 = profd_index.begin(); it4 != profd_index.end(); ++it4)
    {
        info2 << it4->first << " --> " << it4->second << ", ";
    }
    info2 << endl;

    PopulationData pop = knownLoci();

    // default ProfileData
    ProfileData id(prof_type_unknown,
                   "No ID",
                   1,     // contributors
                   0);    // delta

    // if a dbpvec is passed in, use its m_data as the default
    if (dbpvec->size() > 0)
    {
        id = (*dbpvec)[0].m_data;
        dbpvec->clear();
    }

    // read profiles
    string line;
    while (getline(ifs, line))
    {
        try
        {
            vector<string> words = split2(line, locus_sep);
            cleanup(words);

            // skip blank lines (require at least one locus present)
            if (words.size() < 3) continue;

            map<int, vector<PMF<Allele> > > loci; // loci read on the current line
            DBProfile dbp(id);

            string sample_id;
            string profile_id;

            // Read all columns
            vector<string>::const_iterator it;
            int n; // index of locus in the file
            for(n = 0, it = words.begin();
                it != words.end();
                ++n, ++it)
            {

                if (n == profile_col)
                {
                    profile_id = *it;
                    continue;
                }

                if (n == sample_col)
                {
                    sample_id = *it;
                    continue;
                }

                if (locus_index.find(n) != locus_index.end())
                {
                    // handle locus fields

                    if (pop.hasLocus(locus_index[n]))
                    {
                        if (it->size() > 0) // non-empty column
                        {
                            loci[locus_index[n]].push_back(decodeKitAllelePMF(*it));
                        }
                    }
                    else
                    {
                        // Unknown locus - warning already issued. Skip this column.
                    }

                }
                else if (profd_index.find(n) != profd_index.end())
                {
                    // handle ProfileData fields

                    string val = *it; // prof data value
                    if (val.size() == 0) continue;

                    string col_header = profd_index[n];

                    setProfileData(dbp.m_data, col_header, val);
                }
                else if (metad_index.find(n) != metad_index.end())
                {
                    // handle metadata fields

                    string val = *it; // metadata value
                    if (val.size() == 0) continue;

                    string meta_field = metad_index[n];
                    dbp.m_metadata[meta_field] = val;
                }
            }

            // all columns read - store the profile
            if (loci.size() > 0)
            {
                // we need a profile ID!
                if (profile_id.empty())
                {
                    profile_id = "NA";
                }

                string ident = sample_id.substr(0, Database::SAMPLE_ID_LEN) + "-" + profile_id.substr(0, Database::PROFILE_ID_LEN);

                // check for duplicated identifier
                (void)deDup(index, ident); // throws if fails

                // Use the ProfileData in the DBProfile (which might have changed)
                ProfileData &id = dbp.m_data;
                id.m_id = ident;
                id.m_sample_id  = sample_id;
                id.m_profile_id = profile_id;

                // construct the locus data in a Profile
                Profile p(id);

                map<int, vector<PMF<Allele> > >::iterator il;
                for (il = loci.begin(); il != loci.end(); ++il)
                {
                    int loc = il->first;
                    vector<PMF<Allele> > &vec = il->second;

                    // tri-allele. We do not yet handle this, so if we have >2 alleles
                    // replace this locus with 'F' 'F'
                    if (vec.size() > 2 )
                    {
                        vec.clear();
                    }

                    // copy vec into profile
                    for (size_t i=0; i<vec.size(); ++i)
                    {
                        p[loc].setAllele(i, vec[i]);
                    }

                    // and make sure we have at last two alleles by padding with 'F'
                    for (size_t i=vec.size(); i<2; ++i)
                    {
                        p[loc].setAllele(i, PMF<Allele>());
                    }

                    // add text representation of the locus to DBProfile
                    pmfsToText( dbp[loc], p[loc].getAllAlleles() );
                }

                info << startl << "readAnyProfiles(): " << p << endl;

                // store DBProfile
                index[ident] = dbpvec->size();
                dbpvec->push_back(dbp);
            }
        }
        catch(std::exception e)
        {
            error << startl  << "error reading line (skipping):" << endl;
            error << alignl << line << endl;
            continue;
        };
    }

    return true;
}

// read profiles and merge unknown alleles into the background
//
// If dbpvec is non-NULL then we are in "database mode" and do not merge with the background.
// Instead put the profiles in the dbpvec array
//
bool
ProfileFilter::readKitProfiles(
        std::vector<Profile>               &pvec,         // OUT Profiles read
		std::istream                       &ifs,          // IN
		ProfileFilter::FileFormat    const &format,       // IN
        std::vector<DBProfile>             *dbpvec        // OUT DBProfiles read (optional)
        ) const
{
    bool database_mode = (dbpvec != 0); // read input stream for input into database

    std::map<std::string, int> index;   // check for duplicate identifiers

	Unknowns unknown = ignore;
	std::string unknown_alleles = getStringEnv("UNKNOWN_ALLELES");

	if (unknown_alleles == "IGNORE")
	{
		unknown = ignore;
	}
	else if (unknown_alleles == "ADD_AS_RARE")
	{
		unknown = add_as_rare;
	}
	else if (unknown_alleles != "")
	{
		warn << startl << "UNKNOWN_ALLELES: should be \"IGNORE\" or \"ADD_AS_RARE\" (treating as IGNORE)" << endl;
	}

	// Read column headers
	// * first column is identifier
	// * second column is profile type
	// * subsequent columns are in pairs: "LOCUSNAME 1" "LOCUSNAME 2"
	//   where LOCUSNAME is one of the loci in m_locus_name_map and m_locus_names
	string line;
	getline(ifs, line);

	vector<string> headers = split2(line, locus_sep);
	cleanup(headers);

	// skip columns
	vector<string>::const_iterator it = headers.begin();
	for (size_t i=0; i<format.skip_columns + 1 && it != headers.end(); ++i) ++it;

	map<int, int> locus_index; // map from column order to internal order of loci
	int n = 0; // index of locus in the file

	// prepare to add unknown alleles to the population database
	PopulationData pop;
	int pop_size = 0;

	if (database_mode)
	{
        pop = knownLoci();
	}
	else
	{
	    pop = populationData(); // copy
	    pop_size = pop.sampleSize();
	}

	int added = 0;

	while(it != headers.end() && n < format.num_loci)
	{
		vector< vector<string> > h;
		for (int i=0; i<format.cols_per_locus; ++i)
		{
			Assert2(it != headers.end(), "unexpected end of line");

			string col_header = *it++;
			h.push_back( split(col_header, format.header_sep) );

			if (h[i].size() != 2)
			{
				error << startl << "malformed header: '" << col_header << "'" << endl;
				return false;
			}

			char col_num[2] = { '0'+i+1, 0 };
			if ( h[i][1] != col_num )
			{
				error << startl << "misnumbered header: '" << col_header << "' ==> '" << h[i][0] << "' '" << h[i][1] << "' != " << col_num << endl;
				return false;
			}

			if ( (i>0 && ( h[i][0] != h[i-1][0] )))
			{
				error << startl << "mismatched header: '" << col_header << "'" << endl;
				return false;
			}
		}

        string s = h[0][0]; to_upper(s);

		map<string, Locus>::const_iterator map_it = m_locus_name_map.find(s);

		if (map_it == m_locus_name_map.end())
		{
			error << startl << "error: unrecognized locus: " << h[0][0] << endl;
			return false;
		}
		locus_index[n++] = map_it->second;

		if (! pop.hasLocus(map_it->second))
		{
			warn << startl << "locus not in population database (will not be used): " << m_locus_names.at(map_it->second) << endl;
		}
	}

	info2 << endl << endl << startl << "locus_index: " << endl;
	map<int, int>::const_iterator tit;
	for (tit = locus_index.begin(); tit != locus_index.end(); ++tit)
	{
		info2 << tit->first << "-->" << tit->second << ", ";
	}
	info2 << endl << endl;

	Assert2(n == format.num_loci, "error: wrong number of loci");

	// read profiles
	while (getline(ifs, line))
	{
		try
		{
			vector<string> words = split2(line, locus_sep);
			cleanup(words);

			// skip blank lines (require at least one locus present)
			if (words.size() < 1 + format.skip_columns + 2) continue;

			vector<string>::const_iterator it = words.begin();
			string sample_id = *it++;
			string profile_id;

			for (size_t i=0; i<format.skip_columns; ++i)
			{
			    if (i == 0)
			    {
                    // if there is anything in this column, it is the profile_id
                    profile_id = *it;
			    }
			    it++;
			}

			// we need a profile ID!
            if (profile_id.empty())
            {
                profile_id = "NA";
            }

            string ident = sample_id.substr(0, Database::SAMPLE_ID_LEN) + "-" + profile_id.substr(0, Database::PROFILE_ID_LEN);

			// check for duplicated identifier
			(void)deDup(index, ident); // throws if fails

			ProfileData id((ProfileType)format.kit_id, // TODO translate these enums properly
					       ident,
					       1,     // contributors
					       0);    // delta

			id.m_sample_id  = sample_id;
            id.m_profile_id = profile_id;

			Profile p(id);
            DBProfile dbp(id);
			int p_added = 0;

			for(n = 0; // index of locus in the file
			    it != words.end() && n < format.num_loci;
			    ++n)
			{
				if (pop.hasLocus(locus_index[n]))
				{
					vector<PMF<Allele> > pmfs;
					for (int i=0; i<format.cols_per_locus; ++i)
					{
						Assert2(it != words.end(), "unexpected end of line");
						const string &s = *it++;

						if (s.size() > 0) // non-empty column
						{
							pmfs.push_back(decodeKitAllelePMF(s));

							if (!database_mode)
							{
                                // check for alleles not in the population database
                                p_added += checkAlleles(pmfs.back(), pop, (Locus)locus_index[n], m_locus_names.at((Locus)locus_index[n]), unknown, pop_size);
							}
						}
						else if (i<2)
						{
							// To ensure there are at least two alleles recorded, add an 'F' (empty PMF)
							pmfs.push_back(PMF<Allele>());
						}
					}

					// tri-allele. We do not yet handle this, so if we have 3 alleles
					// replace this locus with 'F' 'F'
					if (pmfs.size() > 2 )
					{
						pmfs.clear();
						pmfs.push_back(PMF<Allele>());
						pmfs.push_back(PMF<Allele>());
					}

					p[locus_index[n]] = AlleleSet(pmfs);

                    if (database_mode)
                    {
                        // add text to DBProfile
                        pmfsToText( dbp[locus_index[n]], pmfs );
                    }
				}
				else
				{
					// warning already issued

					// SKIP THESE COLUMNS
					it += format.cols_per_locus;
				}
			}

			const int add_max = 2;
			if (p_added > add_max)
			{
				warn << startl << p_added << " rare alleles added: line may be corrupted (line accepted)"<< endl;
				warn << alignl << line << endl;
			}

			added += p_added;

			if (p.numLoci() > 0)
			{
				info << startl << "readKitProfiles(): " << p << endl;

				pvec.push_back(p);
				index[ident] = pvec.size() - 1;

				if (database_mode)
				{
				    dbpvec->push_back(dbp);
				}
			}
		}
		catch(std::exception e)
		{
			error << startl  << "error reading line (skipping):" << endl;
			error << alignl << line << endl;
			continue;
		};

	}

	// save the population database
	if (added)
	{
		popSet(pop); // NB ignore increase in sample size
	}

	return true;
}

// new version - full language
//
// This function handles reading the file line by line and error handling.
// The entry for each PMF is parsed with a YACC parser, and the resulting PMF entered in the database.
// On error, the entire *sample" is abandoned.
// On completion of a sample, b11PostProcessSample() is called to replace 'D' entries with the
// required PMF (i.e. all the other alleles in the sample, at relative background frequencies);
//
// If dbp is non-NULL then the vector of DBProfiles will be filled and the vector of Profiles should be ignored.
bool
ProfileFilter::readB11Profiles(
		std::vector<Profile>               &pvec,           // OUT Profiles read
		std::istream                       &ifs,            // IN  file data (.b11.csv)
        std::vector<DBProfile>             *dbpvec,         // OUT DBProfiles read (optional)
        Database                     const *db
		) const
{

    std::set<std::string> data_fields;
    MetaField meta_fields;

    bool database_mode = (dbpvec != 0); // read in database mode

    std::map<std::string, int> index;

    if (database_mode)
    {
        Assert2(db, "Database pointer is null");
        db->listMetaFields(meta_fields);
    }

    Database::listDataHeadings(data_fields);

	Unknowns unknown = ignore;
	std::string unknown_alleles = getStringEnv("UNKNOWN_ALLELES");

	if (unknown_alleles == "IGNORE")
	{
		unknown = ignore;
	}
	else if (unknown_alleles == "ADD_AS_RARE")
	{
		unknown = add_as_rare;
	}
	else if (unknown_alleles != "")
	{
		warn << startl << "UNKNOWN_ALLELES: should be \"IGNORE\" or \"ADD_AS_RARE\" (treating as IGNORE)" << endl;
	}

	//
	// Read column headers
	// * first column is identifier
	// * second column is profile type
	// * "Donors"
	// * "Allele"
	// * subsequent columns locus names: "LOCUSNAME"
	//   where LOCUSNAME is one of the loci in m_locus_name_map and m_locus_names
	//
	// after an empty column, there may be MetaData
	string line;
	getline(ifs, line);

	vector<string> headers = split2(line, locus_sep);
	cleanup(headers);

	// skip columns
	vector<string>::const_iterator it = headers.begin();
	size_t skip_columns = 4;
	for (size_t i=0; i<skip_columns && it != headers.end(); ++i) ++it;

	map<int, int> locus_index;       // map from column order to internal order of loci
	map<int, string> metadata_index; // map from column order to metadata headers
    map<int, string> profdata_index; // map from column order to metadata headers
	int n = 0;                       // column number

	// get a copy of the population database (may be modified)
	PopulationData pop = database_mode ? (PopulationData &)knownLoci() : populationData();

	// read locus/metadata headers and generate index
	bool meta_part = false;
    int num_loci = 0;
	while(it != headers.end())
	{
		string col_header = *it++;

		// the first empty column divides loci from metadata
		if (col_header == "")
        {
		    meta_part = true;
		    n++;
		    continue;
        }

		if (!meta_part)
		{
		    // Reading loci. Make sure we recognize the locus name and put it in the index
		    string s = col_header; to_upper(s);
            map<string, Locus>::const_iterator map_it = m_locus_name_map.find(s);

            if (map_it == m_locus_name_map.end())
            {
                error << startl << "readB11Profiles(): error: unrecognized locus: " << col_header << endl;

//                cout << populationData().getLocusInfo() << endl;
//                cout << "m_locus_name_map.size() = " << m_locus_name_map.size() << endl;
//
//                std::map<std::string, Locus>::const_iterator it;
//                for (it = m_locus_name_map.begin(); it != m_locus_name_map.end(); ++it)
//                {
//                    cout << it->first << endl;
//                }
                return false;
            }

            locus_index[n++] = map_it->second;
            num_loci++;

            if (! pop.hasLocus(map_it->second))
            {
                warn << startl << "locus not in population database (will not be used): " << m_locus_names.at(map_it->second) << endl;
            }
        }
		else if (data_fields.find(col_header) != data_fields.end()) // check if column is a ProfileData field
        {
            profdata_index[n++] = col_header;
        }
        else if (meta_fields.find(col_header) != meta_fields.end()) // check if column is a metadata field
        {
            string field_name = meta_fields[col_header].field_name;
            metadata_index[n++] = field_name;
        }
        else
        {
            // Unrecognized field: warn and then ignore
            if (database_mode)
            {
                warn << startl << "Column not Metadata or ProfileData (ignored): " << col_header << endl;
            }
            n++;
        }
	}

	info2 << startl << "METADATA:" << endl;
	map<int, string>::const_iterator itt;
	for (itt = metadata_index.begin(); itt != metadata_index.end(); ++itt)
	{
	    info2 << alignl << itt->first << " = " << itt->second << endl;
	}

	if (num_loci < 1)
	{
		error << startl << "readB11Profiles(): no locus column headers found in file: " << endl;
		return false;
	}

	info2 << endl << endl << startl << "locus_index: " << endl;
	map<int, int>::const_iterator tit;
	for (tit = locus_index.begin(); tit != locus_index.end(); ++tit)
	{
		info2 << tit->first << "-->" << tit->second << ", ";
	}
	info2 << endl << endl;

	//
	// read profiles
	//
	string tok, sample_id, profile_id, ident;
	int n_allele     = -1; // Allele index within profile (0, 1, ...)
	int n_prof       = 0;  // Number of profiles in sample (1, 2, ... )
	bool skip_sample = false;
    int added = 0;

	// default ProfileData
    ProfileData id((ProfileType)B11, // TODO translate these enums properly
                   "No ID",
                   1,     // contributors
                   0);    // delta

    // if a dbpvec is passed in, use its m_data as the default
	if (database_mode && dbpvec->size() > 0)
	{
	    id = (*dbpvec)[0].m_data;
        dbpvec->clear();
	}

	Profile const p(id);     // empty Profile
    DBProfile const dbp(id); // empty DBProfile

	// read lines
	while (getline(ifs, line))
	{
		try
		{
			vector<string> words = split2(line, locus_sep);
			cleanup(words);

			// skip blank lines (require at least one locus present)
			if (words.size() < skip_columns + 1) continue;

			vector<string>::const_iterator it = words.begin();
			tok = *it++; // sample name
			if (tok.size() != 0)
			{
				// new sample
				if (n_prof > 0)
				{
				    added += b11PostProcessSample(pvec, n_prof, locus_index, pop, dbpvec);
				}

				sample_id = tok;
				n_allele = -1;
				n_prof = 0;
				skip_sample = false;
			}

			tok = *it++; // profile name
			if (tok.size() != 0)
			{
				// new profile
				profile_id = tok;
                ident = sample_id.substr(0, Database::SAMPLE_ID_LEN) + id_sep + profile_id.substr(0, Database::PROFILE_ID_LEN);

				n_allele = 0;
				n_prof++;
			}
			else
			{
				n_allele++;
			}

			if (n_allele < 0)
			{
				// no profile ID
				throw std::exception();
			}

			if (skip_sample)
			{
//				error << startl << "readB11Profiles(): skipping sample due to error:" << endl;
				error << alignl << line << endl;
				continue;
			}

			int n_donors = decodeInt(*it++);
			int c_allele = decodeInt(*it++);

			Assert2(c_allele == n_allele + 1, "Alleles mis-numbered");

			if (n_allele == 0) // new profile
			{
				// check for duplicated identifier
				(void)deDup(index, ident); // throws if fails

//				info << startl << "Reading profile: " << ident << endl;

				if (pvec.size() > 0)
				{
					// pvec.back() is a complete profile - take a look
					info << startl << "Read profile with " << pvec.back().numLoci() << " loci: " << endl << pvec.back() << endl;
				}

				// push new empty profile onto DB
				pvec.push_back(p);

				// add data fields
				Assert(n_donors > 0);
				ProfileData d = pvec.back().data();
				d.m_id = ident;
                d.m_profile_id = profile_id;
                d.m_sample_id = sample_id;
				d.m_num_contributors = n_donors;
				pvec.back().setData(d);

				// add profile to index
				index[ident] = pvec.size() - 1;

				if (database_mode)
				{
				    // copy profile data to the DBProfile
				    dbpvec->push_back(dbp);
				    dbpvec->back().m_data = d;
				}
			}

			for(n = 0; // index of locus in the file
			    it != words.end() && n < num_loci;
			    ++n)
			{
				if (pop.hasLocus(locus_index[n]))
				{
					PMF<Allele> pmf; // interpret empty column as 'F' (empty PMF)

					Assert2(it != words.end(), "unexpected end of line");
					const string &s = *it++;

					if (s.size() > 0) // non-empty column
					{
				    	// parse the data
					    PopulationData *pop_ptr = (database_mode? 0 : &pop);
						if (Parser::theParser().parse(s, pmf, pop_ptr) == Parser::yacc_error)
						{
							throw std::exception(); // abandon profile
						}
					}

//					info << startl  << "n = " << n << endl;
//					info << alignl << "locus_index[n] = " << locus_index[n] << endl;
//					info << alignl << "pmf = " << pmf << endl;
//					info << alignl << "db.back()[" << locus_index[n] << "].pushBack(pmf)" << endl;

					// add allele to profile
					pvec.back()[locus_index[n]].setAllele(n_allele, pmf /*, s.c_str()*/ );

					if (database_mode)
					{
					    // add text to DBProfile
	                    Assert2(dbpvec->back().m_b11text[(Locus)locus_index[n]].size() == (size_t)n_allele, "I can't count");
                        dbpvec->back().m_b11text[(Locus)locus_index[n]].push_back(s.c_str());
					}
				}
				else
				{
					// warning already issued

					// SKIP THIS COLUMN
					it++;
				}
			}

            if (n_allele == 0) // ProfileData and metadata must appear on the 'PROFILE' line
            {
                for(; it != words.end(); ++n)
                {
                    string val = *it++; // metadata value
                    if (val.size() == 0) continue;

                    // see if this is ProfileData
                    if (profdata_index.find(n) != profdata_index.end())
                    {
                        // ProfileData
                        string field_name = profdata_index[n];
                        setProfileData(dbpvec->back().m_data, field_name, val);
                    }
                    // else see if it is metadata
                    else if (metadata_index.find(n) != metadata_index.end())
                    {
                        // metadata
                        if (database_mode)
                        {
                            string field_name = metadata_index[n];
                            dbpvec->back().m_metadata[field_name] = val;
                        }
                    }
                    else
                    {
                        // Unrecognized field - warning already given
                    }
                }
            }
		}
		catch(std::exception e) // NB Asserts in the try block throw AssertException which are caught here
		{
			// Some kind of syntax error - skip this entire sample
			error << startl << "error reading line (skipping sample):" << endl;
			error << alignl << line << endl;
			skip_sample = true;

			// clean up profiles already read in this sample
			while (n_prof--)
			{
			    if (pvec.size()) // should be true
			    {
                    Profile &p = pvec.back();
                    index.erase(p.data().m_id);
                    pvec.pop_back();
			    }

				if (database_mode && dbpvec->size() /* may not be true depending which exception was thrown */)
				    dbpvec->pop_back();
			}

			continue;
		};
	} // read lines

	if (!sample_id.empty())
	{
		added += b11PostProcessSample(pvec, n_prof, locus_index, pop, dbpvec);
	}

   // save the population database
    if (!database_mode && added)
    {
        popSet(pop); // NB ignore increase in sample size
    }

	if (pvec.size() > 0)
	{
		// complete profile - take a look
		info << startl << "Read profile with " << pvec.back().numLoci() << " loci: " << endl << pvec.back() << endl;
	}

	Assert(pvec.size() == index.size());

	return true;
}

// The current sample is the last nprof profiles in pvec
// Build the set of all alleles in the sample.
// Set each allele to its background frequency (This is 'D')
// Replace occurrences of 'D' with this set.
//
// If dbpvec is set, do a text substitution for 'D' here instead
// (This is for use when reading into the back-end database)
int
ProfileFilter::b11PostProcessSample(
		std::vector<Profile>                          &pvec,        // IN/OUT Profiles read
		int                                            nprof,       // IN Number of profiles in Sample
		std::map<int, int>                            &locus_index, // IN Locus index
		PopulationData                                &pop,         // IN Population database
        std::vector<DBProfile>                        *dbpvec       // OUT BDProfiles (optional)
		) const
{
    int added = 0;

    bool database_mode = (dbpvec != 0); // read input stream for input into database

	int num_loci = locus_index.size();
	Assert(num_loci > 0);

	for (int n = 0; n<num_loci; ++n)
	{
		int loc = locus_index[n];
		if (! pop.hasLocus(loc)) continue;

	    PMF<Allele> sample_all; // the set of all alleles in the sample, at this locus (ignore the p values)

	    // build sample_all
	    vector<Profile>::reverse_iterator pit = pvec.rbegin();
	    for (int p=0; p<nprof; ++p, ++pit)
	    {
			AlleleSet &s = (*pit)[loc];
			for (int i=0; i<s.size(); ++i)
			{
				const PMF<Allele> &pmf = s.getAllele(i);
				sample_all += pmf;
			}
	    }
	    // NB sample_all may still contain entries for 'D' and 'ignored' alleles

	    // For the DBProfiles, build a string representing 'D'
	    // NB include 'ignored' alleles
	    ostringstream dstream;
	    int dstr_count = 0;
	    dstream << "(";

	    // For the Profiles, set all alleles in sample_all to their background frequencies,
	    // and remove D's and ignored alleles.
        PMF<Allele>::iterator it, current;
        for(it = sample_all.begin(); it != sample_all.end();)
        {
			current = it++;

			if (!database_mode)
			{
			    current->second = pop.getFrequency(loc, current->first, &added);
			}

			if (current->first != Allele::unknown) // not 'D'
			{
			    dstream << current->first;
			    ++dstr_count;

			    if (it != sample_all.end())
			    {
			        dstream << "/";
			    }
			}

        	if (current->first == Allele::unknown || current->second == 0) // 'D' || ignored
        	{
        		sample_all.erase(current);
        	}
        }

        // and normalize
        sample_all.normalize();

        dstream << ")@B";

        // Either we are in 'live' mode (doing a match now) and we need the Profiles, or we are
        // in 'database mode' (writing into the back-end database) when we need only the DBProfiles
        if (database_mode)
        {
            // Go through the DBProfiles and substitute 'D' with dstr
            vector<DBProfile>::reverse_iterator db_it;
            db_it = dbpvec->rbegin();
            for (int p=0; p<nprof; ++p, ++db_it)
            {
                DBProfile::LocusMap::iterator it = db_it->find(loc);
                Assert(it != db_it->m_b11text.end());
                vector<string> &text = it->second;
                for (size_t i=0; i<text.size(); ++i)
                {
                    // search for a 'D' (should be first character if present)
                    size_t pos = text[i].find('D');
                    if (pos != string::npos)
                    {
                        Assert(pos == 0);
                        if (dstr_count > 0)
                        {
                            text[i].replace(0, 1, dstream.str());
                        }
                        else
                        {
                            warn << startl << "'D' used with no alleles in the sample (assuming 'F')" << endl;
                            text[i] = "F";
                        }
                    }
                }
            }
        }
        else // 'immediate' mode
        {
            // Go through the profiles and substitute 'D' with sample_all @ prob
            // Also, remove any ignored alleles
            pit = pvec.rbegin();
            for (int p=0; p<nprof; ++p, ++pit)
            {
                AlleleSet &as = (*pit)[loc];
                for (int i=0; i<as.size(); ++i)
                {
                    PMF<Allele> pmf = as.getAllele(i);

                    for(it = pmf.begin(); it != pmf.end();)
                    {
                        current = it++;

                        if (current->first == Allele::unknown) // 'D'
                        {
                            Assert(pmf.size()==1); // D should be on its own
                            double d_prob = current->second;
                            pmf = sample_all;
                            pmf *= d_prob;
                        }
                        else
                        {
                            // check allele is in population database
                            double f = pop.getFrequency(loc, current->first, &added);
                            if (f == 0)
                            {
                                // Ignore. If the PMF is left empty, that is an 'F'
                                warn << startl << "erasing loc = " << (Locus)loc << " allele = " <<  current->first << endl;
                                pmf.erase(current);
                                warn << alignl << "pmf.size() = " << pmf.size() << endl;
                            }
                        }
                    }

                    as.setAllele(i, pmf);
                }
            }
        }
	}

	return added;
}

enum ProfileFilter::FileType
fileType(string const & filename)
{
	// filename should be:
	// name.idf.csv (Identifiler)
	// name.p16.csv (PowerPlex16)
	// name.sgm.csv (SGM_Plus)

	size_t n = filename.size();
	if (n < 8 || filename.substr(n-4) != ".csv") return ProfileFilter::type_unknown;

	string type_extn = filename.substr(n-8, 4);

	if (type_extn == ".idf") return ProfileFilter::Identifiler;
	if (type_extn == ".p16") return ProfileFilter::PowerPlex16;
	if (type_extn == ".sgm") return ProfileFilter::SGM_Plus; // not handled yet
	if (type_extn == ".b11") return ProfileFilter::B11;

	return ProfileFilter::other_csv;
}

struct ProfileFilter::FileFormat const *
fileFormat(ProfileFilter::FileType ftype)
{
    ProfileFilter::FileFormat const *file_fmt = 0;

    switch (ftype)
    {
    case ProfileFilter::Identifiler:
        file_fmt = &identifiler_fmt;
        break;
    case ProfileFilter::PowerPlex16:
        file_fmt = &powerplex16_fmt;
        break;
    default:
        // file_fmt = 0;
        break;
    }

    return file_fmt;
}

struct ProfileFilter::FileFormat const *
fileFormat(string const &filename)
{
    return fileFormat(fileType(filename));
}

bool readProfilesFromFile(
        std::vector<Profile>               &pvec,           // OUT Profiles read
		std::string                  const &filename,     // IN
		std::vector<DBProfile>             *dbpvec,       // OUT (optional)
        Database                     const *db)
{

	// open file
	std::ifstream ifs(filename.c_str());
	if (! ifs.good()){
		warn << startl << "error opening file: " << filename << endl;
		return false;
	}

	ProfileFilter::FileFormat const *file_fmt = fileFormat(filename);

	if (file_fmt)
	{
		return ProfileFilter::theProfileFilter().readKitProfiles(pvec, ifs, *file_fmt, dbpvec);
	}
	else if (fileType(filename) == ProfileFilter::B11)
	{
		return ProfileFilter::theProfileFilter().readB11Profiles(pvec, ifs, dbpvec, db);
	}
	else
	{
		warn << startl << "unrecognised file format: " << filename << endl;
		return false;
	}
}

bool readProfilesFromDB(
        std::vector<Profile>               &pvec,         // OUT
        std::string                  const &dataset_name, // IN
        std::string                  const &db_name)      // IN
{
    Database db(db_name);
    if (!db.connect())
    {
        cout << "Failed to connect to database. Giving up" << endl;
        return false;
    }

    return dbReadDataset(db, dataset_name, pvec);
}

bool readProfiles(
        ProfileRange                       &pr,           // OUT
        std::string                  const &source,       // IN - file or dataset name (if database given)
        std::string                  const &db_name)      // IN - database name (optional)

{
    if (!db_name.empty())
    {
        // don't read - just check the size and return a ProfileRangeDB
        boost::shared_ptr<Database> db( new Database(db_name) );

        if (!db->connect())
        {
            cout << "Failed to connect to database. Giving up" << endl;
            return false;
        }

        string sql_stmt = Database::datasetQuery(source);
        ProfileRangeDB prdb(db, sql_stmt);
        pr = ProfileRange(prdb);
        return true;

    }
    else
    {
        // read from file
        std::vector<Profile> *pv = new std::vector<Profile>;
        if ( ! readProfilesFromFile(*pv, source) )
        {
            delete pv;
            return false;
        }

        ProfileRangeRefCtdVec pr_rcv(pv);
        pr = ProfileRange(pr_rcv);
        return true;
    }
}

// write profiles in Identifiler format
void writeProfiles(vector<Profile> const &db,
				   string const &filename,
				   map<Locus, string> const &locus_names /*, vector of loci in this file
                                                             (not necessarily all of them,
                                                              not necessarily in order) */ )
{
	// open file
	std::ofstream ofs(filename.c_str());
	if (! ofs.good()){
		error << startl << "could not open file: " << filename << endl;
		exit_debug (-2);
	}

	// Write column headers
	// * first column is identifier
	// * second column is profile type
	// * subsequent columns are in pairs: "LOCUSNAME 1" "LOCUSNAME 2"
	//   where LOCUSNAME is one of the identifiler loci in identifiler_name[]

	// first two column headers
	ofs << "XXXXXXEL-C Number" << locus_sep << "PROFILE" << locus_sep;

	// TODO: use vector of loci in this file.
	// For now just use the first profile, and take the Loci in sorted order
	LocusSet::const_iterator it;
	for(it = db[0].begin(); it != db[0].end(); ++it)
	{
	    map<Locus, string>::const_iterator loc_it = locus_names.find(it->first);
	    if (loc_it == locus_names.end())
	    {
	        error << startl << it->first << " not found while writing " << filename << endl;
	    }
	    else
	    {
	        string loc_name = loc_it->second;
            ofs << loc_name << " 1" << locus_sep << loc_name << " 2" << locus_sep;
	    }
	}
	ofs << endl;

	// write profiles
	vector<Profile>::const_iterator db_it;
	for (db_it = db.begin(); db_it != db.end(); ++db_it)
	{
		ofs << db_it->data().m_id << locus_sep;
		ofs << "NA" << locus_sep;

		LocusSet::const_iterator it;
		for(it = db_it->begin(); it != db_it->end(); ++it)
		{
		    if (locus_names.find(it->first) == locus_names.end())
		    {
                continue; // not an identifiler locus
		    }

			const AlleleSet& as = it->second;

			for (int i=0; i<2; ++i)
			{
				if (i < as.size())
				{
					const PMF<Allele> &pmf = as.getAllele(i);
					if (!pmf.empty())
					{
						Allele a = pmf.mode();
						PMF<Allele>::const_iterator ii = pmf.find(a);
						if (ii->second == 1)
						{
							ofs << a.string() << locus_sep;
						}
						else
						{
							ofs << "(" << a.string() << ")" << locus_sep;
						}

						continue;
					}
				}

				ofs << "F" << locus_sep;
			}
		}
		ofs << endl;
	}

	ofs.close();
}

static char meta_defns[] =
        "TAT,  tat,  VARCHAR(16)\n"
        "DEF,  def,  VARCHAR(16)\n";

// set up a test population database (and destroy it afterwards)
struct Pop_fixture
{
	Pop_fixture() : db("test", meta_defns)
	{
		// population database for FGA:
		PMF<Allele> fga;
		fga[Allele(18)]    = 0.02649;
		fga[Allele(19)]    = 0.05298;
		fga[Allele(20)]    = 0.12748;
		fga[Allele(21)]    = 0.18543;
		fga[Allele(21, 2)] = 0.00497;

		// population database for CSF1PO:
		PMF<Allele> csf1po;
		csf1po[Allele(8)]  = 0.00497;
		csf1po[Allele(9)]  = 0.01159;
		csf1po[Allele(10)] = 0.21689;
		csf1po[Allele(11)] = 0.30132;
		csf1po[Allele(12)] = 0.36093;
		csf1po[Allele(13)] = 0.09603;

		PopulationData testpop;
		testpop.setLocus(CSF1PO, csf1po);
		testpop.setLocus(FGA, fga);
		testpop.setSampleSize(302);

		popSet(testpop);
	}

	~Pop_fixture()
	{
		popClear();
	}

	Database db;
};

TEST_FIXTURE(Pop_fixture, ProfileFilter1)
{
    vector<Profile>  pvec;
    map<string, int> index;

    // generates error message "no locus column headers found in file"
    istringstream iss1("not a valid b11.csv file");

    CHECK( ! ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss1) );
}


static char file1[] = {
"SAMPLE      ,PROFILE ,DONORS ,ALLELE ,CSF1PO  ,FGA     ,        ,Notes   ,\n"
"EXAMPLE1    ,FULL    ,1      ,1      ,8       ,18      ,        ,A single contributor profile.\n"
"            ,        ,       ,2      ,9       ,19      ,        ,        ,\n"
};

void
makeIndex(
    vector<Profile>  const &pvec,
    map<string, int> &index)
{
    size_t n = pvec.size();
    for (size_t i=0; i<n; ++i)
    {
        index[pvec[i].data().m_id] = i;
    }
}

TEST_FIXTURE(Pop_fixture, ProfileFilter2)
{
	vector<Profile>  pvec;
	map<string, int> index;

	istringstream iss1(file1);

	CHECK( ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss1) );
	CHECK_EQUAL(1, pvec.size());       // one profile
	CHECK_EQUAL(2, pvec[0].numLoci()); // with two loci

	makeIndex(pvec, index);

	CHECK_EQUAL(1, index.size());
	CHECK(index.find("EXAMPLE1-FULL") != index.end());
	CHECK_EQUAL(0, index.find("EXAMPLE1-FULL")->second);

	// CSFPO is 8,9 (detailed check)
	CHECK_EQUAL(2, pvec[0][CSF1PO].size());                               // two alleles
	CHECK_EQUAL(1, pvec[0][CSF1PO].getAllele(0).size());                  // one component in PMF
	CHECK_EQUAL(1, pvec[0][CSF1PO].getAllele(0).find(Allele(8))->second); // 8 @ probability 1
	CHECK_EQUAL(1, pvec[0][CSF1PO].getAllele(1).size());                  // one component in PMF
	CHECK_EQUAL(1, pvec[0][CSF1PO].getAllele(1).find(Allele(9))->second); // 9 @ probability 1

	// (simple check)
	Allele a1, a2;
	CHECK(pvec[0][CSF1PO].simple(a1, a2));
	CHECK_EQUAL(Allele(8), a1);
	CHECK_EQUAL(Allele(9), a2);

	// FGA is 18, 19 (simple check)
	CHECK(pvec[0][FGA].simple(a1, a2));
	CHECK_EQUAL(Allele(18), a1);
	CHECK_EQUAL(Allele(19), a2);
}

static char file2[] = {
"SAMPLE      ,PROFILE ,DONORS ,ALLELE ,CSF1PO    ,FGA   ,        ,Notes   ,\n"
"EXAMPLE2    ,FULL    ,1      ,1      ,8@.9      ,F     ,        ,8@.9 implies background at 0.1,\n"
"            ,        ,       ,2      ,9@.5/8@.1 ,F     ,        ,9@.5/8@.1 implies background at 0.4.,\n"
};

TEST_FIXTURE(Pop_fixture, ProfileFilter4)
{
	vector<Profile>  pvec;
	map<string, int> index;

	istringstream iss2(file2);

	CHECK( ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss2) );
	CHECK_EQUAL(1, pvec.size());       // one profile
	CHECK_EQUAL(2, pvec[0].numLoci()); // with two loci

    makeIndex(pvec, index);

	CHECK_EQUAL(1, index.size());
	CHECK(index.find("EXAMPLE2-FULL") != index.end());
	CHECK_EQUAL(0, index.find("EXAMPLE2-FULL")->second);

	// CSFPO is 8@.9, 9@.5/8@.1
	CHECK_EQUAL(2,   pvec[0][CSF1PO].size());                                     // two alleles
	CHECK_EQUAL(1,   pvec[0][CSF1PO].getAllele(0).size());                        // one component in PMF
	CHECK_CLOSE(0.9, pvec[0][CSF1PO].getAllele(0).find(Allele(8))->second, 1e-6); // 8 @ 0.9
	CHECK_EQUAL(2,   pvec[0][CSF1PO].getAllele(1).size());                        // two components in PMF
	CHECK_CLOSE(0.1, pvec[0][CSF1PO].getAllele(1).find(Allele(8))->second, 1e-6); // 8 @ 0.1
	CHECK_CLOSE(0.5, pvec[0][CSF1PO].getAllele(1).find(Allele(9))->second, 1e-6); // 9 @ 0.5

	// FGA is F, F
	CHECK_EQUAL(2, pvec[0][FGA].size());                                          // two alleles
	CHECK_EQUAL(0, pvec[0][FGA].getAllele(0).size());                             // zero components in PMF
	CHECK_EQUAL(0, pvec[0][FGA].getAllele(1).size());                             // zero components in PMF
}

static char file3[] = {
"SAMPLE      ,PROFILE  ,DONORS  ,ALLELE  ,CSF1PO  ,FGA      ,        ,Notes   ,\n"
"EXAMPLE5    ,MAJ      ,1       ,1       ,8       ,18       ,        ,        ,\n"
"            ,         ,        ,2       ,9       ,18       ,        ,        ,\n"
"            ,MINX     ,2       ,1       ,10      ,19       ,        ,        ,\n"
"            ,         ,        ,2       ,11      ,20       ,        ,        ,\n"
"            ,         ,        ,3       ,12      ,21.2     ,        ,        ,\n"
"            ,         ,        ,4       ,13      ,D@0.8    ,        ,D@.8 means any of (18/19/20/21.2) @.8,\n"
};

TEST_FIXTURE(Pop_fixture, ProfileFilter5)
{
	vector<Profile>  pvec;
	map<string, int> index;

	istringstream iss3(file3);

	CHECK( ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss3) );
	CHECK_EQUAL(2, pvec.size());                             // two profiles

    makeIndex(pvec, index);

	CHECK(index.find("EXAMPLE5-MAJ") != index.end());
	CHECK_EQUAL(0, index.find("EXAMPLE5-MAJ")->second);    // called EXAMPLE5-MAJ
	CHECK(index.find("EXAMPLE5-MINX") != index.end());
	CHECK_EQUAL(1, index.find("EXAMPLE5-MINX")->second);   // and EXAMPLE5-MINX

	//
	// EXAMPLE5-MAJ
	//
	CHECK_EQUAL(1, pvec[0].data().m_num_contributors);

	// CSFPO is 8,9
	Allele a1, a2;
	CHECK(pvec[0][CSF1PO].simple(a1, a2));
	CHECK_EQUAL(Allele(8), a1);
	CHECK_EQUAL(Allele(9), a2);

	// FGA is 18, 19
	CHECK(pvec[0][FGA].simple(a1, a2));
	CHECK_EQUAL(Allele(18), a1);
	CHECK_EQUAL(Allele(18), a2);

	//
	// EXAMPLE5-MINX
	//
	CHECK_EQUAL(2, pvec[1].data().m_num_contributors);

	// CSFPO is 10, 11, 12, 13
	CHECK_EQUAL(4,   pvec[1][CSF1PO].size());                                      // four alleles
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(0).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(0).find(Allele(10))->second, 1e-6); // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(1).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(1).find(Allele(11))->second, 1e-6); // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(2).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(2).find(Allele(12))->second, 1e-6); // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(3).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(3).find(Allele(13))->second, 1e-6); // 10 @ 1

	// FGA is 19, 20, 21.2, D@0.8
	CHECK_EQUAL(4,   pvec[1][FGA].size());                                         // four alleles
	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(0).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(0).find(Allele(19))->second, 1e-6);    // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(1).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(1).find(Allele(20))->second, 1e-6);    // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(2).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(2).find(Allele(21, 2))->second, 1e-6); // 10 @ 1

	// population database for FGA:
	//
	//	18 0.02649
	//	19 0.05298
	//	20 0.12748
	//	21 0.18543
	//	21.2 0.00497
	//  ...

	CHECK_EQUAL(4,   pvec[1][FGA].getAllele(3).size());                            // four components in PMF

	float sum  = 0.02649 + 0.05298 + 0.12748 + 0.00497;
	float p18  = 0.8 * 0.02649 / sum;
	float p19  = 0.8 * 0.05298 / sum;
	float p20  = 0.8 * 0.12748 / sum;
	float p212 = 0.8 * 0.00497 / sum;

	CHECK_CLOSE(p18, pvec[1][FGA].getAllele(3).find(Allele(18))->second, 1e-6);
	CHECK_CLOSE(p19, pvec[1][FGA].getAllele(3).find(Allele(19))->second, 1e-6);
	CHECK_CLOSE(p20, pvec[1][FGA].getAllele(3).find(Allele(20))->second, 1e-6);
	CHECK_CLOSE(p212, pvec[1][FGA].getAllele(3).find(Allele(21, 2))->second, 1e-6);
}

// New B11 features: To test:
// * D not the last component
// * use of brackets
// * error recovery - skip entire sample

static char file4[] = {
"SAMPLE      ,PROFILE  ,DONORS  ,ALLELE  ,CSF1PO    ,FGA      ,        ,TAT ,Notes   ,\n"
"EXAMPLE2    ,BAD1     ,1       ,1       ,8@.9      ,F        ,        ,E1  ,        ,\n"
"            ,         ,        ,2       ,9@.5/     ,F        ,        ,EX  ,syntax error at CSFPO - reject sample,\n"
"EXAMPLE5    ,MAJ      ,1       ,1       ,8         ,18       ,        ,E5a ,        ,\n"
"            ,         ,        ,2       ,9         ,18       ,        ,    ,        ,\n"
"            ,MINX     ,2       ,1       ,10        ,19       ,        ,E5b ,        ,\n"
"            ,         ,        ,2       ,13        ,D@0.8    ,        ,    ,        ,\n"
"            ,         ,        ,3       ,11        ,20       ,        ,    ,        ,\n"
"            ,         ,        ,4       ,12        ,21.2     ,        ,    ,        ,\n"
"EXAMPLE3    ,FULL     ,1       ,1       ,8@.9      ,F        ,        ,E3  ,        ,\n"
"            ,         ,        ,2       ,9@.5/8@.1 ,FF       ,        ,    ,syntax error at FGA - reject sample,\n"
"EXAMPLE6    ,MAJ      ,1       ,1       ,8         ,18       ,        ,E6  ,        ,\n"
"            ,         ,        ,2       ,9         ,18       ,        ,    ,        ,\n"
"            ,MINX     ,2       ,1       ,10        ,19       ,        ,    ,        ,\n"
"            ,         ,        ,2       ,13        ,D@@0.8   ,        ,    ,syntax error at FGA - reject sample,\n"
"            ,         ,        ,3       ,11        ,20       ,        ,    ,        ,\n"
"            ,         ,        ,4       ,12        ,21.2     ,        ,    ,        ,\n"
"EXAMPLE4    ,FULL     ,1       ,1       ,8@.9      ,F        ,        ,E4  ,        ,\n"
"            ,         ,        ,2       ,9@.1/8@.1 ,F        ,        ,    ,        ,\n"
};

// batch (pmatch) mode
TEST_FIXTURE(Pop_fixture, ProfileFilter6)
{
	vector<Profile>  pvec;
	map<string, int> index;

	istringstream iss3(file4);

	CHECK( ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss3) );
	CHECK_EQUAL(3, pvec.size());                             // three profiles

    makeIndex(pvec, index);

	CHECK(index.find("EXAMPLE5-MAJ") != index.end());
	CHECK_EQUAL(0, index.find("EXAMPLE5-MAJ")->second);    // called EXAMPLE5-MAJ
	CHECK(index.find("EXAMPLE5-MINX") != index.end());
	CHECK_EQUAL(1, index.find("EXAMPLE5-MINX")->second);   // and EXAMPLE5-MINX
	CHECK(index.find("EXAMPLE4-FULL") != index.end());
	CHECK_EQUAL(2, index.find("EXAMPLE4-FULL")->second);   // and EXAMPLE4-FULL

	//
	// EXAMPLE5-MAJ
	//
	CHECK_EQUAL(1, pvec[0].data().m_num_contributors);

	// CSFPO is 8,9
	Allele a1, a2;
	CHECK(pvec[0][CSF1PO].simple(a1, a2));
	CHECK_EQUAL(Allele(8), a1);
	CHECK_EQUAL(Allele(9), a2);

	// FGA is 18, 18
	CHECK(pvec[0][FGA].simple(a1, a2));
	CHECK_EQUAL(Allele(18), a1);
	CHECK_EQUAL(Allele(18), a2);

	//
	// EXAMPLE5-MINX
	//
	CHECK_EQUAL(2, pvec[1].data().m_num_contributors);

	// CSFPO is 10, 13, 11, 12
	CHECK_EQUAL(4,   pvec[1][CSF1PO].size());                                      // four alleles
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(0).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(0).find(Allele(10))->second, 1e-6); // 10 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(1).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(1).find(Allele(13))->second, 1e-6); // 13 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(2).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(2).find(Allele(11))->second, 1e-6); // 11 @ 1
	CHECK_EQUAL(1,   pvec[1][CSF1PO].getAllele(3).size());                         // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(3).find(Allele(12))->second, 1e-6); // 12 @ 1

	// FGA is 19, 20, 21.2, D@0.8
	CHECK_EQUAL(4,   pvec[1][FGA].size());                                         // four alleles
	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(0).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(0).find(Allele(19))->second, 1e-6);    // 10 @ 1

	// see below for Allele 1

	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(2).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(2).find(Allele(20))->second, 1e-6);    // 10 @ 1

	CHECK_EQUAL(1,   pvec[1][FGA].getAllele(3).size());                            // one component in PMF
	CHECK_CLOSE(1.0, pvec[1][FGA].getAllele(3).find(Allele(21, 2))->second, 1e-6); // 10 @ 1

	// Allele 1 is D@0.8
	// population database for FGA:
	//
	//	18 0.02649
	//	19 0.05298
	//	20 0.12748
	//	21 0.18543
	//	21.2 0.00497
	//  ...

	CHECK_EQUAL(4,   pvec[1][FGA].getAllele(1).size());                            // four components in PMF

	float sum  = 0.02649 + 0.05298 + 0.12748 + 0.00497;
	float p18  = 0.8 * 0.02649 / sum;
	float p19  = 0.8 * 0.05298 / sum;
	float p20  = 0.8 * 0.12748 / sum;
	float p212 = 0.8 * 0.00497 / sum;

	CHECK_CLOSE(p18, pvec[1][FGA].getAllele(1).find(Allele(18))->second, 1e-6);
	CHECK_CLOSE(p19, pvec[1][FGA].getAllele(1).find(Allele(19))->second, 1e-6);
	CHECK_CLOSE(p20, pvec[1][FGA].getAllele(1).find(Allele(20))->second, 1e-6);
	CHECK_CLOSE(p212, pvec[1][FGA].getAllele(1).find(Allele(21, 2))->second, 1e-6);


	//
	// EXAMPLE4-FULL
	//

	CHECK_EQUAL(2, pvec[2].numLoci()); // two loci

	// CSFPO is 8@.9, (9/8)@.1
	CHECK_EQUAL(2,   pvec[2][CSF1PO].size());                                     // two alleles
	CHECK_EQUAL(1,   pvec[2][CSF1PO].getAllele(0).size());                        // one component in PMF
	CHECK_CLOSE(0.9, pvec[2][CSF1PO].getAllele(0).find(Allele(8))->second, 1e-6); // 8 @ 0.9
	CHECK_EQUAL(2,   pvec[2][CSF1PO].getAllele(1).size());                        // two components in PMF
	CHECK_CLOSE(0.1, pvec[2][CSF1PO].getAllele(1).find(Allele(8))->second, 1e-6); // 8 @ 0.1
	CHECK_CLOSE(0.1, pvec[2][CSF1PO].getAllele(1).find(Allele(9))->second, 1e-6); // 9 @ 0.1

	// FGA is F, F
	CHECK_EQUAL(2, pvec[2][FGA].size());                                          // two alleles
	CHECK_EQUAL(0, pvec[2][FGA].getAllele(0).size());                             // zero components in PMF
	CHECK_EQUAL(0, pvec[2][FGA].getAllele(1).size());                             // zero components in PMF

}

// Test database mode
TEST_FIXTURE(Pop_fixture, ProfileFilter7)
{
    vector<Profile>  pvec;
    vector<DBProfile>  dbvec;
    map<string, int> index;

    istringstream iss3(file4);

    CHECK( ProfileFilter::theProfileFilter().readB11Profiles(pvec, iss3, &dbvec, &db) );
    CHECK_EQUAL(3, dbvec.size());                          // three profiles

    makeIndex(pvec, index);

    CHECK(index.find("EXAMPLE5-MAJ") != index.end());
    CHECK_EQUAL(0, index.find("EXAMPLE5-MAJ")->second);    // called EXAMPLE5-MAJ
    CHECK(index.find("EXAMPLE5-MINX") != index.end());
    CHECK_EQUAL(1, index.find("EXAMPLE5-MINX")->second);   // and EXAMPLE5-MINX
    CHECK(index.find("EXAMPLE4-FULL") != index.end());
    CHECK_EQUAL(2, index.find("EXAMPLE4-FULL")->second);   // and EXAMPLE4-FULL

    //
    // EXAMPLE5-MAJ
    //
    CHECK_EQUAL(1, dbvec[0].m_data.m_num_contributors);

    // CSFPO is 8,9
    CHECK_EQUAL(2, dbvec[0][CSF1PO].size()); // two alleles
    CHECK_EQUAL("8", dbvec[0][CSF1PO][0]);
    CHECK_EQUAL("9", dbvec[0][CSF1PO][1]);

    // FGA is 18, 18
    CHECK_EQUAL(2, dbvec[0][FGA].size()); // two alleles
    CHECK_EQUAL("18", dbvec[0][FGA][0]);
    CHECK_EQUAL("18", dbvec[0][FGA][1]);

    //
    // EXAMPLE5-MINX
    //
    CHECK_EQUAL(2, dbvec[1].m_data.m_num_contributors);

    // CSFPO is 10, 13, 11, 12
    CHECK_EQUAL(4, dbvec[1][CSF1PO].size()); // four alleles
    CHECK_EQUAL("10", dbvec[1][CSF1PO][0]);
    CHECK_EQUAL("13", dbvec[1][CSF1PO][1]);
    CHECK_EQUAL("11", dbvec[1][CSF1PO][2]);
    CHECK_EQUAL("12", dbvec[1][CSF1PO][3]);

    // FGA is 19, D@0.8, 20, 21.2
    CHECK_EQUAL(4,                          dbvec[1][FGA].size()); // four alleles
    CHECK_EQUAL("19",                       dbvec[1][FGA][0]);
    CHECK_EQUAL("(18/19/20/21.2)@B@0.8",    dbvec[1][FGA][1]);
    CHECK_EQUAL("20",                       dbvec[1][FGA][2]);
    CHECK_EQUAL("21.2",                     dbvec[1][FGA][3]);

    //
    // EXAMPLE4-FULL
    //
    CHECK_EQUAL(2, dbvec[2].numLoci()); // two loci

    // CSFPO is 8@.9, 9@.1/8@.1
    CHECK_EQUAL(2,              dbvec[2][CSF1PO].size()); // two alleles
    CHECK_EQUAL("8@.9",         dbvec[2][CSF1PO][0]);
    CHECK_EQUAL("9@.1/8@.1",    dbvec[2][CSF1PO][1]);

    // FGA is F, F
    CHECK_EQUAL(2,              dbvec[2][FGA].size()); // two alleles
    CHECK_EQUAL("F",            dbvec[2][FGA][0]);
    CHECK_EQUAL("F",            dbvec[2][FGA][1]);

    // check metadata
    CHECK(dbvec[0].m_metadata.end() != dbvec[0].m_metadata.find("tat"));
    CHECK(dbvec[1].m_metadata.end() != dbvec[1].m_metadata.find("tat"));
    CHECK(dbvec[2].m_metadata.end() != dbvec[2].m_metadata.find("tat"));
    CHECK_EQUAL("E5a", dbvec[0].m_metadata["tat"]);
    CHECK_EQUAL("E5b", dbvec[1].m_metadata["tat"]);
    CHECK_EQUAL("E4",  dbvec[2].m_metadata["tat"]);

}

static char file5[] = {
"XXXXXXEL-C Number, PROFILE,    D8S1179 1,  D8S1179 2,  D21S11 1,   D21S11 2,   D7S820 1,   D7S820 2,   CSF1PO 1,   CSF1PO 2,   D3S1358 1,  D3S1358 2,  THO1 1,     THO1 2,     D13S317 1,  D13S317 2,  D16S539 1,  D16S539 2,  D2S1338 1,  D2S1338 2,  D19S433 1,  D19S433 2,  vWA 1,  vWA 2,  TPOX 1, TPOX 2, D18S51 1,   D18S51 2,   AMXXXXXXEL 1,AMXXXXXXEL 2,  D5S818 1,   D5S818 2,   FGA 1,  FGA 2,   TAT\n"
"C0,                FULL,       15,         11,         29,         29,         9,          10,         11,         F,          15,         15,         8,          6,          11,         12,         13,         11,         23,         16,         14.2,       13,         16,     15,     11,     11,     14,         14,         X,          X,              12,         11,         18,     19,      ONE\n"
"C1,                FULL,       11,         14,         28,         30,         7,          10,         12,         R,          16,         15,         6,          9.3,        12,         13,         12,         9,          20,         23,         15,         14,         19,     19,     8,      11,     14,         16,         X,          X,              13,         13,         18,     19,      TWO\n"
"C2,                FULL,       11,         13,         32.2,       30,         10,         9,          11,         (10),       18,         17,         7,          9.3,        11,         10,         11,         12,         23,         26,         13,         13.2,       18,     16,     11,     10,     15,         15,         X,          Y,              12,         11,         21.2,   21,      TH REE\n"
"C3,                FULL,       13,         10,         29,         29,         10,         11,         ?(12),      -12,        16,         16,         6,          6,          14,         11,         11,         11,         17,         20,         15,         14,         16,     14,     8,      8,      15,         15,         X,          Y,              10,         13,         21,     21.2,    FOUR\n"
};

// readKitProfiles - batch (pmatch) mode
TEST_FIXTURE(Pop_fixture, ProfileFilter8)
{
    vector<Profile>  pvec;
    map<string, int> index;

    istringstream iss3(file5);

    CHECK( ProfileFilter::theProfileFilter().readKitProfiles(pvec, iss3, identifiler_fmt) );
    CHECK_EQUAL(3, pvec.size());                     // 4 - 1 profiles

    makeIndex(pvec, index);

    CHECK(index.find("C0-FULL") != index.end());
    CHECK_EQUAL(0, index.find("C0-FULL")->second);   // called C0-FULL
    CHECK(index.find("C1-FULL") == index.end());     // and C1-FULL has a syntax error
//  CHECK_EQUAL(1, index.find("C1-FULL")->second);
    CHECK(index.find("C2-FULL") != index.end());
    CHECK_EQUAL(1, index.find("C2-FULL")->second);   // and C2-FULL
    CHECK(index.find("C3-FULL") != index.end());
    CHECK_EQUAL(2, index.find("C3-FULL")->second);   // and C3-FULL

    //
    // C0-FULL
    //
    CHECK_EQUAL(1, pvec[0].data().m_num_contributors);

    // CSFPO is 11,F
    CHECK_EQUAL(2, pvec[0][CSF1PO].size());              // number of alleles
    CHECK_EQUAL(1, pvec[0][CSF1PO].getAllele(0).size()); // components in PMF
    CHECK_CLOSE(1.0, pvec[0][CSF1PO].getAllele(0).find(Allele(11))->second, 1e-6); // 11@1
    CHECK_EQUAL(0, pvec[0][CSF1PO].getAllele(1).size()); // components in PMF (F)

    // FGA is 18, 19
    Allele a1, a2;
    CHECK(pvec[0][FGA].simple(a1, a2));
    CHECK_EQUAL(Allele(18), a1);
    CHECK_EQUAL(Allele(19), a2);

    //
    // C2-FULL
    //
    CHECK_EQUAL(1, pvec[1].data().m_num_contributors);

    // CSFPO is 11, (10)
    CHECK_EQUAL(2, pvec[1][CSF1PO].size());              // number of alleles
    CHECK_EQUAL(1, pvec[1][CSF1PO].getAllele(0).size()); // components in PMF
    CHECK_CLOSE(1.0, pvec[1][CSF1PO].getAllele(0).find(Allele(11))->second, 1e-6); // 11@1
    CHECK_EQUAL(1, pvec[1][CSF1PO].getAllele(1).size()); // components in PMF
    CHECK_CLOSE(BRACKETS_FREQ, pvec[1][CSF1PO].getAllele(1).find(Allele(10))->second, 1e-6); // (10)

    // FGA is 21.2, 21
    CHECK(pvec[1][FGA].simple(a1, a2));
    CHECK_EQUAL(Allele(21), a1);         // held numerically ordered
    CHECK_EQUAL(Allele(21, 2), a2);

    //
    // C3-FULL
    //
    CHECK_EQUAL(1, pvec[2].data().m_num_contributors);

    // CSFPO is ?(12), -12
    CHECK_EQUAL(2, pvec[2][CSF1PO].size());              // number of alleles
    CHECK_EQUAL(1, pvec[2][CSF1PO].getAllele(0).size()); // components in PMF
    CHECK_CLOSE(QBRACKETS_FREQ, pvec[2][CSF1PO].getAllele(0).find(Allele(12))->second, 1e-6); // ?(12)
    CHECK_EQUAL(1, pvec[2][CSF1PO].getAllele(1).size()); // components in PMF
    CHECK_CLOSE(BRACKETS_FREQ, pvec[2][CSF1PO].getAllele(1).find(Allele(12))->second, 1e-6); // -12

    // FGA is 21.2, 21
    CHECK(pvec[2][FGA].simple(a1, a2));
    CHECK_EQUAL(Allele(21), a1);         // held numerically ordered
    CHECK_EQUAL(Allele(21, 2), a2);
}

void pf9_test_index(map<string, int> &index)
{
    CHECK(index.find("C0-FULL") != index.end());
    CHECK_EQUAL(0, index.find("C0-FULL")->second);   // called C0-FULL
    CHECK(index.find("C1-FULL") == index.end());     // and C1-FULL has a syntax error
    CHECK(index.find("C2-FULL") != index.end());
    CHECK_EQUAL(1, index.find("C2-FULL")->second);   // and C2-FULL
    CHECK(index.find("C3-FULL") != index.end());
    CHECK_EQUAL(2, index.find("C3-FULL")->second);   // and C3-FULL
}

void pf9_test_vec(vector<DBProfile> &dbvec)
{
    CHECK_EQUAL(3, dbvec.size());                     // C0-FULL, C2-FULL, C3-FULL

    // C0-FULL
    CHECK_EQUAL(1, dbvec[0].m_data.m_num_contributors);

    // CSFPO is 11,F
    CHECK_EQUAL(2, dbvec[0][CSF1PO].size());              // number of alleles
    CHECK_EQUAL("11", dbvec[0][CSF1PO][0]);
    CHECK_EQUAL("F", dbvec[0][CSF1PO][1]);

    // FGA is 18, 19
    CHECK_EQUAL(2, dbvec[0][FGA].size());              // number of alleles
    CHECK_EQUAL("18", dbvec[0][FGA][0]);
    CHECK_EQUAL("19", dbvec[0][FGA][1]);

    //
    // C2-FULL
    //
    CHECK_EQUAL(1, dbvec[1].m_data.m_num_contributors);

    // CSFPO is 11, (10)
    CHECK_EQUAL(2, dbvec[1][CSF1PO].size());              // number of alleles

    char buf[32];
    CHECK_EQUAL("11", dbvec[1][CSF1PO][0]);
    sprintf(buf, "10@%g", BRACKETS_FREQ);
    CHECK_EQUAL(buf, dbvec[1][CSF1PO][1]);

    // FGA is 21.2, 21
    CHECK_EQUAL(2, dbvec[1][FGA].size());              // number of alleles
    CHECK_EQUAL("21.2", dbvec[1][FGA][0]);
    CHECK_EQUAL("21", dbvec[1][FGA][1]);

    //
    // C3-FULL
    //
    CHECK_EQUAL(1, dbvec[2].m_data.m_num_contributors);

    // CSFPO is ?(12), -12
    CHECK_EQUAL(2, dbvec[2][CSF1PO].size());              // number of alleles
    sprintf(buf, "12@%g", QBRACKETS_FREQ);
    CHECK_EQUAL(buf, dbvec[2][CSF1PO][0]);
    sprintf(buf, "12@%g", BRACKETS_FREQ);
    CHECK_EQUAL(buf, dbvec[2][CSF1PO][1]);

    // FGA is 21, 21.2
    CHECK_EQUAL(2, dbvec[2][FGA].size());              // number of alleles
    CHECK_EQUAL("21", dbvec[2][FGA][0]);
    CHECK_EQUAL("21.2", dbvec[2][FGA][1]);
}

// readKitProfiles - Database mode
TEST_FIXTURE(Pop_fixture, ProfileFilter9)
{
    vector<Profile>  pvec;
    vector<DBProfile>  dbvec;
    map<string, int> index;

    istringstream iss3(file5);

    CHECK( ProfileFilter::theProfileFilter().readKitProfiles(pvec, iss3, identifiler_fmt, &dbvec) );

    pf9_test_vec(dbvec);
    makeIndex(pvec, index);
    pf9_test_index(index);

    // check metadata (none is read by readKitProfiles)
    CHECK(dbvec[0].m_metadata.end() == dbvec[0].m_metadata.find("tat"));
    CHECK(dbvec[1].m_metadata.end() == dbvec[1].m_metadata.find("tat"));
    CHECK(dbvec[2].m_metadata.end() == dbvec[2].m_metadata.find("tat"));
}

// test readAnyProfiles
TEST_FIXTURE(Pop_fixture, ProfileFilter10)
{
    vector<DBProfile>  dbvec;

    istringstream iss3(file5);

    CHECK( ProfileFilter::theProfileFilter().readAnyProfiles(db, iss3, &dbvec, 0, 1) );

    // same tests as readKitProfiles
    pf9_test_vec(dbvec);

    // check ProfileData - should be the defaults set in readAnyProfiles
    CHECK_EQUAL("C0-FULL", dbvec[0].m_data.m_id);
    CHECK_EQUAL("", dbvec[0].m_data.m_dataset);
    CHECK_EQUAL("C0", dbvec[0].m_data.m_sample_id);
    CHECK_EQUAL("FULL", dbvec[0].m_data.m_profile_id);
    CHECK_EQUAL((ProfileType)0, dbvec[0].m_data.m_kit_type);
    CHECK_EQUAL((EvidenceType)0, dbvec[0].m_data.m_evidence_type);
    CHECK_EQUAL(1, dbvec[0].m_data.m_num_contributors);
    CHECK_CLOSE(0, dbvec[0].m_data.m_error_rate, 1e-6);

    // check metadata
    CHECK_EQUAL("ONE", dbvec[0].m_metadata["tat"]);
    CHECK_EQUAL("TH REE", dbvec[1].m_metadata["tat"]);
    CHECK_EQUAL("FOUR", dbvec[2].m_metadata["tat"]);
}

// repeat ProfileFilter10 but passing in default ProfileData
TEST_FIXTURE(Pop_fixture, ProfileFilter11)
{
    vector<DBProfile>  dbvec;

    // set up default ProfileData
    ProfileData pd_def;
    pd_def.m_kit_type       = Identifiler;
    pd_def.m_evidence_type  = crime;
    pd_def.m_dataset        = "test";
    pd_def.m_error_rate     = 0.002;
    dbvec.push_back(DBProfile(pd_def));

    istringstream iss3(file5);

    CHECK( ProfileFilter::theProfileFilter().readAnyProfiles(db, iss3, &dbvec, 0, 1) );

    // same tests as readKitProfiles
    pf9_test_vec(dbvec);

    // check ProfileData - should be the values set in pd_def
    CHECK_EQUAL("C0-FULL", dbvec[0].m_data.m_id);
    CHECK_EQUAL("test", dbvec[0].m_data.m_dataset);
    CHECK_EQUAL("C0", dbvec[0].m_data.m_sample_id);
    CHECK_EQUAL("FULL", dbvec[0].m_data.m_profile_id);
    CHECK_EQUAL(Identifiler, dbvec[0].m_data.m_kit_type);
    CHECK_EQUAL(crime, dbvec[0].m_data.m_evidence_type);
    CHECK_EQUAL(1, dbvec[0].m_data.m_num_contributors);
    CHECK_CLOSE(0.002, dbvec[0].m_data.m_error_rate, 1e-6);

    // check metadata
    CHECK_EQUAL("ONE", dbvec[0].m_metadata["tat"]);
    CHECK_EQUAL("TH REE", dbvec[1].m_metadata["tat"]);
    CHECK_EQUAL("FOUR", dbvec[2].m_metadata["tat"]);
}

static char file6[] = {
"XXXXXXEL-C Number, PROFILE,    DATASET, D8S1179 1,  D8S1179 2,  D21S11 1,   D21S11 2,   D7S820 1,   D7S820 2,   CSF1PO 1,   CSF1PO 2,   D3S1358 1,  D3S1358 2,  THO1 1,     THO1 2,     D13S317 1,  D13S317 2,  D16S539 1,  D16S539 2,  D2S1338 1,  D2S1338 2,  D19S433 1,  D19S433 2,  vWA 1,  vWA 2,  TPOX 1, TPOX 2, D18S51 1,   D18S51 2,   AMXXXXXXEL 1,AMXXXXXXEL 2,  D5S818 1,   D5S818 2,   DEF,  FGA 1,  FGA 2,    TAT\n"
"C0,                FULL,       TEST1,      15,         11,         29,         29,         9,          10,         11,         F,          15,         15,         8,          6,          11,         12,         13,         11,         23,         16,         14.2,       13,         16,     15,     11,     11,     14,         14,         X,          X,              12,         11,      D1,   18,     19,      ONE\n"
"C1,                FULL,       TEST2,      11,         14,         28,         30,         7,          10,         12,         R,          16,         15,         6,          9.3,        12,         13,         12,         9,          20,         23,         15,         14,         19,     19,     8,      11,     14,         16,         X,          X,              13,         13,      D2,   18,     19,      TWO\n" // syntax error 'R'
"C2,                FULL,       TEST3,      11,         13,         32.2,       30,         10,         9,          11,         (10),       18,         17,         7,          9.3,        11,         10,         11,         12,         23,         26,         13,         13.2,       18,     16,     11,     10,     15,         15,         X,          Y,              12,         11,      D3,   21.2,   21,      TH REE\n"
"C3,                FULL,       TEST4,      13,         10,         29,         29,         10,         11,         ?(12),      -12,        16,         16,         6,          6,          14,         11,         11,         11,         17,         20,         15,         14,         16,     14,     8,      8,      15,         15,         X,          Y,              10,         13,      D4,   21,     21.2,    FOUR\n"
};

// repeat ProfileFilter11 but with the "DATASET" column present in the file
TEST_FIXTURE(Pop_fixture, ProfileFilter12)
{
    vector<DBProfile>  dbvec;

    // set up default ProfileData
    ProfileData pd_def;
    pd_def.m_kit_type       = Identifiler;
    pd_def.m_evidence_type  = crime;
    pd_def.m_dataset        = "test";
    pd_def.m_error_rate     = 0.002;
    dbvec.push_back(DBProfile(pd_def));

    istringstream iss3(file6);

    CHECK( ProfileFilter::theProfileFilter().readAnyProfiles(db, iss3, &dbvec, 0, 1) );

    // same tests as readKitProfiles
    pf9_test_vec(dbvec);

    // Check C0-FULL
    // check ProfileData - should be the values set in pd_def
    CHECK_EQUAL("C0-FULL", dbvec[0].m_data.m_id);
    CHECK_EQUAL("TEST1", dbvec[0].m_data.m_dataset);
    CHECK_EQUAL("C0", dbvec[0].m_data.m_sample_id);
    CHECK_EQUAL("FULL", dbvec[0].m_data.m_profile_id);
    CHECK_EQUAL(Identifiler, dbvec[0].m_data.m_kit_type);
    CHECK_EQUAL(crime, dbvec[0].m_data.m_evidence_type);
    CHECK_EQUAL(1, dbvec[0].m_data.m_num_contributors);
    CHECK_CLOSE(0.002, dbvec[0].m_data.m_error_rate, 1e-6);

    CHECK_EQUAL("C2-FULL", dbvec[1].m_data.m_id);
    CHECK_EQUAL("TEST3", dbvec[1].m_data.m_dataset);

    // Check C3-FULL
    // check ProfileData - should be the values set in pd_def
    CHECK_EQUAL("C3-FULL", dbvec[2].m_data.m_id);
    CHECK_EQUAL("TEST4", dbvec[2].m_data.m_dataset);
    CHECK_EQUAL("C3", dbvec[2].m_data.m_sample_id);
    CHECK_EQUAL("FULL", dbvec[2].m_data.m_profile_id);
    CHECK_EQUAL(Identifiler, dbvec[2].m_data.m_kit_type);
    CHECK_EQUAL(crime, dbvec[2].m_data.m_evidence_type);
    CHECK_EQUAL(1, dbvec[2].m_data.m_num_contributors);
    CHECK_CLOSE(0.002, dbvec[2].m_data.m_error_rate, 1e-6);

    // check metadata
    CHECK_EQUAL("ONE", dbvec[0].m_metadata["tat"]);
    CHECK_EQUAL("TH REE", dbvec[1].m_metadata["tat"]);
    CHECK_EQUAL("FOUR", dbvec[2].m_metadata["tat"]);
    CHECK_EQUAL("D1", dbvec[0].m_metadata["def"]);
    CHECK_EQUAL("D3", dbvec[1].m_metadata["def"]);
    CHECK_EQUAL("D4", dbvec[2].m_metadata["def"]);

}

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ProfileFilter.h
 *
 *  Created on: Mar 11, 2010
 *      Author: gareth
 */

#ifndef PROFILEFILTER_H_
#define PROFILEFILTER_H_

#include <map>
#include <vector>
#include <string>
#include <map>

#include "Profile.h"
#include "ProfileRange.h"
struct DBProfile;
class Database;

const double BRACKETS_FREQ_DEFAULT  = 0.7;
const double QBRACKETS_FREQ_DEFAULT = 0.6;

class ProfileFilter {
public:
	enum FileType     // NB corresponds with ProfileType but they are not really the same!
	{
		type_unknown = 0,
		Identifiler,
		PowerPlex16,
		SGM_Plus,
		B11,
		other_csv,
	};

	struct FileFormat
	{
		FileType                            kit_id;         // E.g. Identifiler
		size_t                              skip_columns;   // no of columns after the ID and before the first data column
		int                                 num_loci;       // number of loci in file
		int                                 cols_per_locus; // number of columns for alleles at each locus
		char                  const * const header_sep;     // separator between locus name and column number (e.g. "-" for "AMEL-1 AMEL-2" format)
	};

	ProfileFilter();
	virtual ~ProfileFilter();

	static ProfileFilter const & theProfileFilter();

	bool
	readAnyProfiles(
	        Database                     const &db,
	        std::istream                       &ifs,          // IN
	        std::vector<DBProfile>             *dbpvec,       // OUT DBProfiles
	        int sample_col   =  0,                            // IN column containing sample ID
	        int profile_col  = -1                             // IN column containing profile ID
	        ) const;

	// read Identifiler/PowerPlex16 files
	bool
	readKitProfiles(
			std::vector<Profile>               &pvec,          // OUT
			std::istream                       &ifs,           // IN
			FileFormat                   const &format,        // IN
            std::vector<DBProfile>             *dbpvec = 0     // OUT DBProfiles read
			) const;

	// read B11 files
	// If dbp is non-NULL then the vector of DBProfiles will be filled and the vector of Profiles should be ignored.
	bool
	readB11Profiles(
	        std::vector<Profile>               &pvec,           // OUT Profiles read
	        std::istream                       &ifs,            // IN  file data (.b11.csv)
	        std::vector<DBProfile>             *dbpvec = 0,     // OUT DBProfiles read
            Database                     const *db = 0
	        ) const;

private:
	ProfileFilter(const ProfileFilter &);
	ProfileFilter & operator=(const ProfileFilter &);

	int                                                                 // RET alleles added
	b11PostProcessSample(
	        std::vector<Profile>                          &pvec,        // IN/OUT Profile read
	        int                                            nprof,       // IN Number of profiles
	        std::map<int, int>                            &locus_index, // IN Locus index
	        PopulationData                                &pop,         // IN Population database
	        std::vector<DBProfile>                        *dbpvec       // OUT Profile read
	        ) const;

    bool
    readHeaders(
            Database             const &db,
            std::istream               &ifs,         // IN:  input stream
            std::map<int, int>         &locus_index, // OUT: column -> Locus
            std::map<int, std::string> &metad_index, // OUT: column -> metadata field name
            std::map<int, std::string> &profd_index  // OUT: column -> ProfileData field name
            ) const;

	std::map<std::string, Locus> m_locus_name_map; // map from identifiler name to position in AlleleSet
	std::map<Locus, std::string> m_locus_names;    // map from position to name
};

// output a PMF<Allele> in human-readable form
std::ostream &
operator<<(std::ostream &os, const PMF<Allele> &pmf);

// output a Profile in human-readable form
std::ostream &
operator<<(std::ostream &os, const Profile &p);

// read profiles and merge unknown alleles into the background
bool readProfilesFromFile(
        std::vector<Profile>               &pvec,         // OUT Profiles read
        std::string                  const &filename,     // IN
        std::vector<DBProfile>             *dbpvec = 0,   // OUT (optional)
        Database                     const *db = 0);

bool readProfilesFromDB(
        std::vector<Profile>               &pvec,         // OUT
        std::string                  const &dataset_name, // IN
        std::string                  const &db_name);     // IN

bool readProfiles(
        ProfileRange                         &pr,       // OUT
        std::string                  const &source,       // IN - file or dataset name
        std::string                  const &db_name = "");// IN - database name ("" if file)

void writeProfiles(std::vector<Profile> const &db,
				   std::string const &filename,
				   std::map<Locus, std::string> const &locus_names);
					// TODO FileFormat format);

// if pop = 0 then accept any Locus/Allele combination
Profile
makeProfile(DBProfile const &dbp, PopulationData *pop = 0);

enum ProfileFilter::FileType
fileType(std::string const & filename);

struct ProfileFilter::FileFormat const *
fileFormat(ProfileFilter::FileType ftype);

struct ProfileFilter::FileFormat const *
fileFormat(std::string const &filename);

#endif /* PROFILEFILTER_H_ */

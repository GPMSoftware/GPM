/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * dbmatch.h
 *
 *  Created on: Jul 16, 2010
 *      Author: gareth
 */

#ifndef DBMATCH_H_
#define DBMATCH_H_

#include "Profile.h"
#include "Database.h"
#include "ProfileFilter.h"
#include "match.h"
#include "SubPopModel.h"
#include "MutModel.h"

#include "cuda_accel/GPUDevices.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

// merges with population database
bool dbReadProfile(
        Database                           &db,           // IN
        std::string                  const &profile_name, // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index = 0);   // OUT

// merges with population database
bool dbReadDataset(
        Database                           &db,           // IN
        std::string                  const &dataset_name, // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index = 0);   // OUT (if non-null)

// does not merge with population database
bool dbReadRawDataset(
        Database                           &db,           // IN
        std::string                  const &dataset_name, // IN
        std::vector<Profile>               &pvec);        // OUT

// read SQL query into vector of Profiles
// NB the query must return entire profiles, i.e. "SELECT * FROM Profiles WHERE ..."
bool dbReadQuery(
        boost::shared_ptr<Database>        &db,           // IN
        std::string                  const &sql_query,    // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index = 0);   // OUT

bool dbReadQueryNoBackground(
        boost::shared_ptr<Database>        &db,           // IN
        std::string                  const &sql_query,    // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index = 0);   // OUT

void
copyNtoNM(NMmatchResults &matches, NmatchResults const &results, bool one_to_many, int n);

// match profile against profile
// NB single match is always done on the CPU (general method)
double
single_match(
    Profile& p1,
    Profile& p2,
    MatchType const &match_type,
    SubPopModel const & spm,
    MutModel const &mut = MutModel());

int
n_match(NmatchResults& matches,
        ProfileRange & db,
        Profile& p,
        double lr_threshold,
        double delta,
        std::vector<MatchType> const &match_types,
        SubPopModel const & spm,
        bool n_cuda_accel);

int n2_match(NMmatchResults &matches,
             ProfileRange & db,
             double lr_threshold,
             double delta,
             std::vector<MatchType> const &match_types,
             SubPopModel const & spm,
             bool n2_cuda_accel,
             Split split);

int nm_match(NMmatchResults &matches,
             ProfileRange & db1,
             ProfileRange & db2,
             double lr_threshold,
             double delta1,
             double delta2,
             std::vector<MatchType> const &match_types,
             SubPopModel const & spm,
             bool n2_cuda_accel,
             Split split);

struct Env
{
    Env() {}

    bool run_unit_tests;
    std::string crime_datasets;
    std::string ref_datasets;
    std::string crime_files;
    std::string ref_files;
    std::string out_file;
    std::string crime_id;
    std::string ref_id;
    std::string relatives;
    double crime_delta;
    double ref_delta;
    double lr_threshold;
    std::string popdata;
    std::string unknown_alleles;
    std::string sub_pop;
    double theta;
    double brackets_freq;
    double qbrackets_freq;
    std::string processor;
    int  ngpus;
    bool output_profs;
    bool crime_from_db;
    bool ref_from_db;

    void read();
};

std::ostream & operator<<(std::ostream &os, Env const &env);

struct MatchParams
{
    MatchParams();
    MatchParams(Env const &env);

    std::vector<MatchType> match_types;
    double lr_threshold;
    double crime_delta;
    double ref_delta;
    MutModel mut;
    SubPopModel spm;
    bool n_cuda_accel;
    bool n2_cuda_accel;
};

//void readEnv(Env &env, bool &n_cuda_accel, bool &n2_cuda_accel, bool &output_profs);

void
matchTypesFromString(
        std::vector<MatchType>         &match_types, // OUT
        std::vector<std::string> const & rel_list);  // IN
bool
import(
    Database &db,
    bool overwrite,
    std::ifstream &ifs,
    ProfileFilter::FileType ftype,
    ProfileType prof_type,
    std::string const &dataset,
    EvidenceType evidence_type,
    float error_rate,
    int sample_col   =  0,                            // IN column containing sample ID
    int profile_col  = -1                             // IN column containing profile ID
    );

bool
exportDB(
    Database &db,
    std::ofstream &ofs);

// the match function used by the gui (and the regression tests)
template <class OUT> // OUT is anything supporting <<
void
do_gmatch(
    NMmatchResults          &matches,
    OUT                     &os,
    Split                   split,
    ProfileRange            &crime_db,
    ProfileRange            &ref_db,
    bool                    ref_is_crime,
    float                   crime_delta,
    float                   ref_delta,
    std::vector<MatchType> const &match_types,
    SubPopModel const       &spm,
    MutModel const          &mut,              // used for single match only
    float                   lr_threshold,
    bool                    n_cuda_accel,
    bool                    n2_cuda_accel)
{
    using namespace std;

    if (crime_db.size() == 1)
    {
        int crime_record = 0;
        Profile p1 = crime_db.get(crime_record);
        p1.setErrorRate(crime_delta);

        if (ref_db.size() == 1)
        {
            // one-to-one match
            int ref_record   = 0;
            Profile p2 = ref_db.get(ref_record);
            p2.setErrorRate(ref_delta);

            os << "Matching these profiles: " << p1.data().m_id << ", " << p2.data().m_id << endl << endl;

            matches.n_comparisons = 0; // NA
            for (vector<MatchType>::const_iterator it = match_types.begin(); it != match_types.end(); ++it)
            {
                double likelihood_ratio = single_match(p1, p2, *it, spm, mut);
                matches[ make_pair(crime_record, ref_record) ].push_back(Result(*it, likelihood_ratio));
            }
        }
        else
        {
            // N match: p1 against the reference database

            os << "Matching the following crime profile against the reference database:" << p1.data().m_id << endl << endl;

            NmatchResults results(ref_db.size());
            int n = n_match(results, ref_db, p1, lr_threshold, ref_delta, match_types, spm, n_cuda_accel);

            copyNtoNM(matches, results, true, crime_record);

            os << n << " matches found" << endl;
        }
    }
    else
    {
        if (ref_db.size() == 1)
        {
            // N match: p2 against the crime database

            int ref_record = 0;
            Profile p2 = ref_db.get(ref_record);
            p2.setErrorRate(ref_delta);

            os << "Matching the following reference profile against the crime database:" << p2.data().m_id << endl << endl;

            NmatchResults results(crime_db.size());
            int n = n_match(results, crime_db, p2, lr_threshold, crime_delta, match_types, spm, n_cuda_accel);

            copyNtoNM(matches, results, false, ref_record);

            os << n << " matches found" << endl;
        }
        else
        {
            // N2 match
            if (ref_is_crime)
            {
                os << "Matching the crime database against itself..." << endl;
                int s = crime_db.size();
                matches.n_comparisons = s*(s-1)/2;
                int n = n2_match(matches, crime_db, lr_threshold, crime_delta, match_types, spm, n2_cuda_accel, split);
                os << n << " matches found" << endl;
            }
            else
            {
                os << "Matching crime database against reference database..." << endl;
                matches.n_comparisons = crime_db.size() * ref_db.size();
                int n = nm_match(matches, crime_db, ref_db, lr_threshold, crime_delta, ref_delta, match_types, spm, n2_cuda_accel, split);
                os << n << " matches found" << endl;
            }
        }
    }
}

struct MatchResults
{
    std::vector<MatchType>     const &match_types;   // The match types that were tested
    NMmatchResults             const &results;       // Matches found
    std::map<int, std::string> const &c_ids;         // IDs of matching profiles in the C-set
    std::map<int, std::string> const &r_ids;         // IDs of matching profiles in the R-set
    double                            lr_thresh;     // threshold applied

    MatchResults(
            std::vector<MatchType>     const &t,
            NMmatchResults             const &res,
            std::map<int, std::string> const &c,
            std::map<int, std::string> const &r,
            double                            lr)
    : match_types(t)
    , results(res)
    , c_ids(c)
    , r_ids(r)
    , lr_thresh(lr)
    {}
};

void
outputMatches(
        std::ostream        &ofs,
        MatchResults  const &matches);

void getProfileNames(std::map<int, std::string> &index, ProfileRange const &db, int chunk_size);

#endif /* DBMATCH_H_ */

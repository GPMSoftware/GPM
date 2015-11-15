/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * dbmatch.cpp
 *
 *  Created on: Jul 16, 2010
 *      Author: gareth
 */

#include "dbmatch.h"
#include "cuda_accel/cuda_match.h" // TODO: remove?
#include "cuda_accel/GPUDevices.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("dbmatch");
#include "fand/messages.h"
using namespace std;

// convert to profiles and make index (if requested)
bool convertDBPstoProfiles(
        vector<DBProfile>            const &dbpvec,       // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index,        // OUT
        PopulationData                     *pop)          // IN
{
    for (size_t i=0; i<dbpvec.size(); ++i)
    {
        Profile p = makeProfile(dbpvec[i], pop);
        pvec.push_back(p);
        if (index)
        {
            (*index)[p.data().m_id] = pvec.size() - 1;
        }

        info2 << "readDBDataset(): " << p << endl;
    }

    return true;
}

bool dbReadDatasetImpl(
        Database                           &db,           // IN
        std::string                  const &dataset_name, // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index,        // OUT
        PopulationData                     *pop)          // IN
{
    // read DBProfiles
    vector<DBProfile> dbpvec;

    if (! db.readDataset(dataset_name, dbpvec))
    {
        return false;
    }

    return convertDBPstoProfiles(dbpvec, pvec, index, pop);
}

bool dbReadQueryImpl(
        boost::shared_ptr<Database>        &db,           // IN
        std::string                  const &sql_query,    // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index,        // OUT
        PopulationData                     *pop)          // IN
{
    // read DBProfiles
    vector<DBProfile> dbpvec;

    if (! db->readQuery(sql_query, dbpvec))
    {
        return false;
    }

    return convertDBPstoProfiles(dbpvec, pvec, index, pop);
}

// does not merge with population database
bool dbReadRawDataset(
        Database                           &db,           // IN
        std::string                  const &dataset_name, // IN
        std::vector<Profile>               &pvec)         // OUT
{
    return dbReadDatasetImpl(db, dataset_name, pvec, 0, 0);
}

bool dbReadDataset(
        Database                           &db,           // IN
        std::string                  const &dataset_name, // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index)        // OUT
{
    PopulationData pop = populationData();

    bool ret;
    if (ret = dbReadDatasetImpl(db, dataset_name, pvec, index, &pop) != 0)
    {
        // The populationData may have been modified
        popSet(pop);
    }

    return ret;
}

bool dbReadQuery(
        boost::shared_ptr<Database>        &db,           // IN
        std::string                  const &sql_query,    // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index)        // OUT
{
    PopulationData pop = populationData();

    bool ret;
    if (ret = dbReadQueryImpl(db, sql_query, pvec, index, &pop) != 0)
    {
        // The populationData may have been modified
        popSet(pop);
    }

    return ret;
}

bool dbReadQueryNoBackground(
        boost::shared_ptr<Database>        &db,           // IN
        std::string                  const &sql_query,    // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index)        // OUT
{
    return dbReadQueryImpl(db, sql_query, pvec, index, 0);
}

bool dbReadProfile(
        Database                           &db,           // IN
        std::string                  const &profile_name, // IN
        std::vector<Profile>               &pvec,         // OUT
        std::map<std::string, int>         *index)        // OUT
{
    // read DBProfile
    ProfileData data;
    DBProfile dbprof(data);

    if (! db.readProfile(profile_name, dbprof))
    {
        return false;
    }

    // convert to profile and add to index (if requested)
    PopulationData pop = populationData();

    Profile p = makeProfile(dbprof, &pop);
    pvec.push_back(p);
    if (index)
    {
        (*index)[p.data().m_id] = pvec.size() - 1;
    }

    info2 << "readDBProfile(): " << p << endl;

    // The populationData may have been modified
    popSet(pop);

    return true;
}

void
copyNtoNM(NMmatchResults &matches, NmatchResults const &results, bool one_to_many, int n)
{
    NmatchResults::const_iterator it;
    for(it = results.begin(); it != results.end(); ++it)
    {
        if (one_to_many)
        {
            // one profile in db1 vs all in db2
            matches[ make_pair(n, it->first) ] = it->second;
        }
        else
        {
            // all in db1 vs one profile in db2
            matches[ make_pair(it->first, n) ] = it->second;
        }
    }
    matches.n_comparisons = results.n_comparisons;
}

// calculate LR as a linear combination of basic LRs
// Balding p126 Eqn. 7.14 (Rp, Ri the other way up)
// NB you must do this *at each locus separately*
double
lr_linear(
    MatchType const &match_type,
    double lr_ident,
    double lr_d1)
{
    double k0=0, k1=0, k2=0; // ibd coefficients
    getIBD(match_type, k0, k1, k2);

    return k0 + k1*lr_d1 + k2*lr_ident;
}

// match profile against profile
// NB single match is always done on the CPU (general method)
double
single_match(
    Profile& p1,
    Profile& p2,
    MatchType const &match_type,
    SubPopModel const &spm,
    MutModel const &mut)
{
    // set mutation rate. The ancestor must be the first argument to relmatch

    p1.setMutRate(0);
    p2.setMutRate(0);

    switch (mut.type)
    {
    case MutModel::CRIME_IS_ANCESTOR:
        p1.setMutRate(mut.mutation_rate);
        // fall through

    case MutModel::NO_MUTATIONS:
        return relmatch(match_type, p1, p2, spm);
        break;

    case MutModel::REF_IS_ANCESTOR:
        p2.setMutRate(mut.mutation_rate);
        return relmatch(match_type, p2, p1, spm);
        break;

    default:
        Assert2(false, "Unknown mutation model");
        break;
    }
} // compiler warning: OK

// match Profile against Profile database (GPU)
int
n_match_cuda(NmatchResults& matches,
             ProfileRange &db,
             Profile& p,
             double lr_threshold,
             double delta,
             MatchType const &match_type,
             SubPopModel const & spm)
{
    return gpuDevices().n_match(matches, db, p, lr_threshold, delta, match_type, spm);
}

// match Profile database against itself (GPU)
int
n2_match_cuda(NMmatchResults &matches,
              ProfileRange &db,
              double lr_threshold,
              double delta,
              MatchType const &match_type,
              SubPopModel const& spm,
              Split split)
{
    return gpuDevices().n2_match(matches, db, lr_threshold, delta, match_type, split, spm);
}

int
nm_match_cuda(NMmatchResults &matches,
              ProfileRange & dbA,
              ProfileRange & dbB,
              double lr_threshold,
              double delta1,
              double delta2,
              MatchType const &match_type,
              SubPopModel const& spm,
              Split split)
{
    return gpuDevices().nm_match(matches, dbA, dbB, lr_threshold, delta1, delta2, match_type, split, spm);
}

// match Profile against Profile database (CPU)
int
n_match_cpu(NmatchResults& matches,
            ProfileRange &db,
            Profile& p,
            double lr_threshold,
            double delta,
            vector<MatchType> const &match_types,
            SubPopModel const & spm)
{
    int ntypes = match_types.size();
    vector<double> lr(ntypes);

    p.setErrorRate(delta);

    db.readMyData();

    for(size_t i=0; i<db.size(); ++i)
    {
        // calculate the LRs we have been asked for, for this pair of profiles
        double lr_max = 0;
        for (int k=0; k<ntypes; ++k)
        {
            lr[k] = single_match(p, db[i], match_types[k], spm);

            if (lr[k] > lr_max)
            {
                lr_max = lr[k];
            }
        }

        // if any of them pass the threshold, store them all
        if (lr_max > lr_threshold)
        {
            // user feedback
           ok << ".";
           ok.flush();

            for (int k=0; k<ntypes; ++k)
            {
                matches[i].push_back(Result(match_types[k], lr[k]));
            }
        }
    }
    ok << endl;

    return matches.size();
}

int
n_match(NmatchResults& matches,
        ProfileRange & db,
        Profile& p,
        double lr_threshold,
        double delta,
        vector<MatchType> const &match_types,
        SubPopModel const & spm,
        bool n_cuda_accel)
{
    Timer t;
    MemMeter mem;

	if (n_cuda_accel && gpuDevices().GPUsEnabled() > 0)
    {
        // on cuda we do each relationship match separately,
        // (so each one is separately thresholded)
        for (vector<MatchType>::const_iterator it = match_types.begin(); it != match_types.end(); ++it)
        {
            (void)n_match_cuda(matches, db, p, lr_threshold, delta, *it, spm);
        }
    }
    else
    {
        // on CPU we threshold on the highest reported LR for each pair of profiles
        (void)n_match_cpu(matches, db, p, lr_threshold, delta, match_types, spm);
    }

	t.stop();
    info << startl << "total elapsed time =  " << t << " seconds" << endl;
    info << startl << "memory in use " << mem << endl;

	return matches.size();
}

// match Profile database against itself (CPU)
int
n2_match_cpu(NMmatchResults &matches,
             ProfileRange &db,
             double lr_threshold,
             double delta,
             vector<MatchType> const &match_types,
             SubPopModel const & spm)
{
    Timer t;
    MemMeter mem;

    int ntypes = match_types.size();
    vector<double> lr(ntypes);

    db.readMyData();

    for(size_t i=0; i<db.size(); ++i)
    {
        db[i].setErrorRate(delta);

        for(size_t j=0; j<i; ++j)
        {
            // calculate the LRs we have been asked for, for this pair of profiles
            double lr_max = 0;
            for (int k=0; k<ntypes; ++k)
            {
                lr[k] = single_match(db[i], db[j], match_types[k], spm);

                if (lr[k] > lr_max)
                {
                    lr_max = lr[k];
                }
            }

            // if any of them pass the threshold, store them all
            if (lr_max > lr_threshold)
            {
                // user feedback
               ok << ".";
               ok.flush();

                for (int k=0; k<ntypes; ++k)
                {
                    // NB (j, i) to get the results in the right order
                    matches[ make_pair(j, i) ].push_back(Result(match_types[k], lr[k]));
                }
            }
        }
    }
    ok << endl;

    t.stop();
    info << startl << "total elapsed time =  " << t << " seconds" << endl;
    info << startl << "memory in use " << mem << endl;

    return matches.size();
}

int n2_match(NMmatchResults &matches,
             ProfileRange &db,
             double lr_threshold,
             double delta,
             vector<MatchType> const &match_types,
             SubPopModel const & spm,
             bool n2_cuda_accel,
             Split split)
{
    Timer t;
    MemMeter mem;

	if (n2_cuda_accel && gpuDevices().GPUsEnabled() > 0)
    {
        // on cuda we do each relationship match separately,
        // (so each one is separately thresholded)
        for (vector<MatchType>::const_iterator it = match_types.begin(); it != match_types.end(); ++it)
        {
            (void)n2_match_cuda(matches, db, lr_threshold, delta, *it, spm, split);
        }
    }
    else
    {
        // on CPU we threshold on the highest reported LR for each pair of profiles
        (void)n2_match_cpu(matches, db, lr_threshold, delta, match_types, spm);
    }

    t.stop();
    info << startl << "total elapsed time =  " << t << " seconds" << endl;
    info << startl << "memory in use " << mem << endl;

    return matches.size();
}

// match Profile database against a different Profile database
int
nm_match_cpu(NMmatchResults &matches,
             ProfileRange & db1,
             ProfileRange & db2,
             double lr_threshold,
             double delta1,
             double delta2,
             vector<MatchType> const &match_types,
             SubPopModel const & spm)
{
    int ntypes = match_types.size();
    vector<double> lr(ntypes);

    Timer t;
    db1.readMyData();
    db2.readMyData();
    t.stop();
    info << startl << "readMyData() took " << t << " seconds" << endl;

    for(size_t i=0; i<db1.size(); ++i)
    {
        db1[i].setErrorRate(delta1);

        for(size_t j=0; j<db2.size(); ++j)
        {
            db2[j].setErrorRate(delta2);

            // calculate the LRs we have been asked for, for this pair of profiles
            double lr_max = 0;
            for (int k=0; k<ntypes; ++k)
            {
                lr[k] = single_match(db1[i], db2[j], match_types[k], spm);

                if (lr[k] > lr_max)
                {
                    lr_max = lr[k];
                }
            }

            // if any of them pass the threshold, store them all
            if (lr_max > lr_threshold)
            {
                // user feedback
               ok << ".";
               ok.flush();

                for (int k=0; k<ntypes; ++k)
                {
                    matches[ make_pair(i, j) ].push_back(Result(match_types[k], lr[k]));
                }
            }
        }
    }
    ok << endl;

    return matches.size();
}

int nm_match(NMmatchResults &matches,
             ProfileRange & db1,
             ProfileRange & db2,
             double lr_threshold,
             double delta1,
             double delta2,
             vector<MatchType> const &match_types,
             SubPopModel const & spm,
             bool n2_cuda_accel,
             Split split)
{
    Timer t;
    MemMeter mem;

	if (n2_cuda_accel && gpuDevices().GPUsEnabled() > 0)
    {
        // on cuda we do each relationship match separately,
        // (so each one is separately thresholded)
        for (vector<MatchType>::const_iterator it = match_types.begin(); it != match_types.end(); ++it)
        {
            (void)nm_match_cuda(matches, db1, db2, lr_threshold, delta1, delta2, *it, spm, split);
        }
    }
    else
    {
        // on CPU we threshold on the highest reported LR for each pair of profiles
        (void)nm_match_cpu(matches, db1, db2, lr_threshold, delta1, delta2, match_types, spm);
    }

	t.stop();
    info << startl << "total elapsed time =  " << t << " seconds" << endl;
    info << startl << "memory in use " << mem << endl;

    return matches.size();
}

float fl(string s)
{
    vector<string> vs = split(s, "/");

    if (vs.size() == 2) // int / int
    {
        return atoi(vs[0].c_str()) / (float)atoi(vs[1].c_str());
    }
    else // float
    {
        return atof(s.c_str());
    }
}

void
matchTypesFromString(
        std::vector<MatchType>         &match_types, // OUT
        std::vector<std::string> const & rel_list)   // IN
{
    for (size_t i=0; i<rel_list.size(); ++i)
    {
        MatchType mt(ident_t);
        string const &relative = rel_list[i];
        try
        {
            if (relative == "IDENT" || relative == "IDENTITY")
            {
                mt = MatchType(ident_t);
            }
            else if (relative == "SIB" || relative == "SIBLING")
            {
                mt = MatchType(sibling_t);
            }
            else if (relative == "D1" || relative == "DEGREE_1")
            {
                mt = MatchType(degree_1_t);
            }
            else if (relative == "D2" || relative == "DEGREE_2")
            {
                mt = MatchType(degree_2_t);
            }
            else if (relative.size() == 3 && relative[0] == 'D' && isdigit(relative[1]) && isdigit(relative[2])) // DXY
            {
                int p = relative[1] - '0';
                int q = relative[2] - '0';

                // ensure p <= q
                if (p>q) swap(p, q);

                // interpret (0,X) as (X, INF)
                if (p==0 && q>0)
                {
                    p=q;
                    q = MatchType::INF;
                }
                mt = MatchType(degree_pq_t, p, q);
            }
            else if (relative.substr(0, 7) == "DEGREE_") // DEGREE_X_Y
            {
                vector<string> vs = split(relative, "_");
                Assert(vs[0] == "DEGREE");
                int p=-1, q=-1;
                istringstream issp(vs[1]);
                issp >> p;
                istringstream issq(vs[2]);
                issq >> q;

                if (p<0 || q<0)
                {
                    warn << startl << "RELATIVE: bad value" << endl;
                    throw std::exception();
                }
                else if (p==1 && q==1)
                {
                    warn << startl << "RELATIVE: DEGREE (1,1) not possible for humans" << endl;
                    throw std::exception();
                }

                // ensure p <= q
                if (p>q) swap(p, q);

                // interpret (0,X) as (X, INF)
                if (p==0 && q>0)
                {
                    p=q;
                    q = MatchType::INF;
                }
                mt = MatchType(degree_pq_t, p, q);
            }
            else if (relative.substr(0, 4) == "GEN[" && *(relative.end()-1) == ']') // general 4-number relationship
            {
                vector<string> vs = split(relative.substr(4, relative.size()-5), ",");
                Assert(vs.size() == 4);

                mt = MatchType(gen_t, fl(vs[0]), fl(vs[1]), fl(vs[2]), fl(vs[3]));
            }
            else if (relative.substr(0, 4) == "INV[" && *(relative.end()-1) == ']') // general 4-number inverse relationship
            {
                vector<string> vs = split(relative.substr(4, relative.size()-5), ",");
                Assert(vs.size() == 4);

                mt = MatchType(inv_t, fl(vs[0]), fl(vs[1]), fl(vs[2]), fl(vs[3]));
            }
            else if (relative != "")
            {
                // some other illegal value
                throw std::exception();
            }

            match_types.push_back(mt);
        }
        catch(...)
        {
            warn << startl << "RELATIVE ignored: should be \"\" or \"IDENT[ITY]\" or \"SIB[LING]\" or \"D1\" or \"D2\" or \"Dpq\" or \"DEGREE_p_q\"" << endl;
    //          warn << alignl << "assuming IDENT" << endl;
        }
    }

    if (match_types.empty())
    {
        match_types.push_back(MatchType(ident_t));
    }
}

void
Env::read()
{
    run_unit_tests   = getBoolEnv("RUN_UNIT_TESTS", true);

    crime_datasets   = getStringEnv("CRIME_DATASETS");    // database
    ref_datasets     = getStringEnv("REF_DATASETS");      // database

    crime_files      = getStringEnv("CRIME_FILES");      // file
    ref_files        = getStringEnv("REF_FILES");        // file
    out_file         = getStringEnv("OUT_FILE");         // file

    crime_id         = getStringEnv("CRIME_ID");
    ref_id           = getStringEnv("REF_ID");
    relatives        = getStringEnv("RELATIVE");
    crime_delta      = getFloatEnv("CRIME_DELTA");
    ref_delta        = getFloatEnv("REF_DELTA");
    lr_threshold     = getFloatEnv("LR_THRESHOLD");
    popdata          = getStringEnv("POPDATA");
    unknown_alleles  = getStringEnv("UNKNOWN_ALLELES");
    sub_pop          = getStringEnv("SUB_POP_MODEL");
    theta            = getFloatEnv("THETA");
    brackets_freq    = getFloatEnv("BRACKETS_FREQ", BRACKETS_FREQ_DEFAULT);
    qbrackets_freq   = getFloatEnv("QBRACKETS_FREQ", QBRACKETS_FREQ_DEFAULT);
    processor        = getStringEnv("PROC");
    ngpus            = getIntEnv("NGPUS", 0);
    output_profs     = getBoolEnv("PRINT_PROFS", true);

    crime_from_db = crime_files.empty();
    ref_from_db   = ref_files.empty();
}

MatchParams::MatchParams()
: lr_threshold(0)
, crime_delta(0)
, ref_delta(0)
, spm(SubPopModel::HW, 0)
, n_cuda_accel(false)
, n2_cuda_accel(true)
{
}

MatchParams::MatchParams(Env const &env)
: lr_threshold(0)
, crime_delta(0)
, ref_delta(0)
, spm(SubPopModel::HW, 0)
, n_cuda_accel(false)
, n2_cuda_accel(true)
{
    crime_delta    = env.crime_delta;
    ref_delta      = env.ref_delta;
    lr_threshold   = env.lr_threshold;

    if (crime_delta < 0 || crime_delta > 1)
    {
        warn << startl << "environment variable CRIME_DELTA out of range, using 0" << endl;
        crime_delta = 0;
    }

    if (ref_delta < 0 || ref_delta > 1)
    {
        warn << startl << "environment variable REF_DELTA out of range, using 0" << endl;
        ref_delta = 0;
    }

    if (lr_threshold < 0)
    {
        warn << startl << "environment variable LR_THRESHOLD out of range, using 0" << endl;
        lr_threshold = 0;
    }

    if (env.processor == "GPU")
    {
        n_cuda_accel = true;
        n2_cuda_accel = true;
    }
    else if (env.processor == "CPU")
    {
        n_cuda_accel = false;
        n2_cuda_accel = false;
    }
    else if (env.processor != "")
    {
        warn << startl << "environment variable PROC ignored (should be \"GPU\" or \"CPU\")" << endl;
    }

    vector<string> rel_list = split(env.relatives,   ":");
    cleanup(rel_list);

    matchTypesFromString(match_types, rel_list);

    if (env.sub_pop == "" || env.sub_pop == "HW")
    {
        // OK
    }
    else if (env.sub_pop == "NRC4.4")
    {
        spm = SubPopModel(SubPopModel::NRC4_4, env.theta);
    }
    else if (env.sub_pop == "NRC4.10")
    {
        spm = SubPopModel(SubPopModel::NRC4_10, env.theta);
    }
    else
    {
        warn << startl << "environment variable SUB_POP_MODEL ignored: should be \"\" or \"HW\" or \"NRC4.4\" or \"NRC4.10\"" << endl;
        warn << alignl << "assuming HW" << endl;
    }
}

std::ostream &
operator<<(std::ostream &os, Env const & env)
{
//    bool crime_from_db;
//    bool ref_from_db;

    os << "Environment variables:" << std::endl;
    os << "RUN_UNIT_TESTS  = " << (env.run_unit_tests ? "true" : "false") << std::endl;
    os << "CRIME_DATASETS  = " << env.crime_datasets << std::endl;
    os << "REF_DATASETS    = " << env.ref_datasets << std::endl;
    os << "CRIME_FILES     = " << env.crime_files << std::endl;
    os << "REF_FILES       = " << env.ref_files << std::endl;
    os << "OUT_FILE        = " << env.out_file << std::endl;
    os << "CRIME_ID        = " << env.crime_id << std::endl;
    os << "CRIME_DELTA     = " << env.crime_delta << std::endl;
    os << "REF_ID          = " << env.ref_id << std::endl;
    os << "REF_DELTA       = " << env.ref_delta << std::endl;
    os << "RELATIVE        = " << env.relatives << std::endl;
    os << "LR_THRESHOLD    = " << env.lr_threshold << std::endl;
    os << "POPDATA         = " << env.popdata << std::endl;
    os << "UNKNOWN_ALLELES = " << env.unknown_alleles << std::endl;
    os << "SUB_POP_MODEL   = " << env.sub_pop << std::endl;
    os << "THETA           = " << env.theta << std::endl;
    os << "PROC            = " << env.processor << std::endl;
    os << "NGPUS           = " << env.ngpus << std::endl;
    os << "PRINT_PROFS     = " << (env.output_profs ? "true" : "false") << std::endl;
    os << "BRACKETS_FREQ   = " << env.brackets_freq << std::endl;
    os << "QBRACKETS_FREQ  = " << env.qbrackets_freq << std::endl;

    return os;
}

void
outputMatches(
        ostream             &ofs,
        MatchResults  const &matches)
{
    string below_thresh = "0";

    if (matches.lr_thresh > 0)
    {
        ostringstream oss;
        oss << "<" << matches.lr_thresh;
        below_thresh = oss.str();
    }

    // write headers, and create a map of the columns
    ofs << "Profile Key 1,Profile Key 2";

    map<MatchType, int> cols;
    int i=0;
    vector<MatchType>::const_iterator mit;
    for(mit = matches.match_types.begin(); mit != matches.match_types.end(); ++mit)
    {
        ofs << "," << *mit;
        cols[*mit] = i++;
    }
    ofs << "\n";

    int n = matches.match_types.size();
    float results[n];

    NMmatchResults::const_iterator it;
    for(it = matches.results.begin(); it != matches.results.end(); ++it)
    {
        int c_index = it->first.first;
        string p1_id = matches.c_ids.at(c_index);

        int r_index = it->first.second;
        string p2_id = matches.r_ids.at(r_index);

        ofs << p1_id << "," << p2_id;

        for(int i=0; i<n; ++i)
        {
            results[i] = -1;
        }

        const vector<Result> &res = it->second;
        for(size_t i=0; i<res.size(); ++i)
        {
            int col = cols[res[i].match_type];
            results[col] = res[i].likelihood_ratio;
        }

        for(int i=0; i<n; ++i)
        {
            if (results[i] == -1)
            {
                ofs << "," << below_thresh;
            }
            else
            {
                ofs << "," << results[i];
            }
        }
        ofs << "\n";
    }
}

// fill in the ids of the profiles indexed in 'index', reading 'chunk_size' profiles at a time
// the index is already filled with the required integers and the strings are overwritten
void
getProfileNames(std::map<int, std::string> &index, ProfileRange const &db, int chunk_size)
{
    std::map<int, std::string>::iterator it;

    if (db.gotData())
    {
        // We have the data - do it all in one go
        for (it = index.begin(); it != index.end(); ++it)
        {
            it->second = db.get(it->first).data().m_id;
        }
    }
    else
    {
        // We must read the data - do it in chunks.
        // NB we iterate through the indices and the chunks in numeric order

        it = index.begin();
        for(size_t s = 0; s < db.size(); s += chunk_size)
        {
            int n = min(s+chunk_size, db.size());
            ProfileRange chunk(db, s, n);
            chunk.readMyData();

            while (it != index.end())
            {
                int i = it->first;
                Assert((size_t)i >= chunk.begin());
                if ((size_t)i < chunk.end())
                {
                    it->second = chunk.get(i-s).data().m_id;
                    ++it; // next index
                }
                else
                {
                    break; // next chunk
                }
            }
        }
    }
}

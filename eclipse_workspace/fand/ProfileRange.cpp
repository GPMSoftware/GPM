/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ProfileRange.cpp
 *
 *  Created on: Jul 1, 2011
 *      Author: gareth
 */

#include "ProfileRange.h"
#include "Database.h"
#include "dbmatch.h"

#include "cuda_accel/GPUDevices.h"

#include <strstream>
#include <UnitTest++/UnitTest++.h>

#include "MessageStream.h"
INIT_MESSAGES("ProfileRange")
#include "messages.h"


ProfileRangeDB::ProfileRangeDB(boost::shared_ptr<Database> &db, std::string const &query)
: ProfileRangeBase(0, 0), m_db(db), m_query(query)
{
    // find number of records returned by the query
    int c = m_db->countQuery(m_query);

    if (c < 0)
    {
        warn << "ProfileRangeDB::ProfileRangeDB(): Database query failed: " << query << endl;
        c = 0;
    }

    m_end = c;

    // TODO - if unknowns policy is ADD_AS_RARE we need to read all the profiles now
    // (without storing them) and update the frequency database.
    Assert(populationData().getPolicy() == ignore);
}

bool
ProfileRangeDB::readData(std::vector<Profile> &vec, size_t begin, size_t end)
{
    // read records from database from m_begin to m_end

    // make sure we have a connection to the database in this thread
    if (!m_db->connect())
    {
        error << startl << "Database connection failed" << endl;
        return false;
    }

    std::ostringstream oss;
    oss << m_query << " ORDER BY profile_key LIMIT " << begin << ", " << end - begin;

    vec.clear();
    bool ret = dbReadQuery(m_db, oss.str(), vec);

    if (!ret || (end - begin != vec.size()))
    {
        error << startl << "ProfileRangeDB::readData(): Failed to read profiles from database" << endl;
        return false;
    }

    return ret;
}

bool
ProfileRangeDB::readMyData()
{
    m_data_ok = readData(m_vec, m_begin, m_end);
//    cout << "ProfileRangeDB::readMyData(): " << this << " m_data_ok = " << m_data_ok << endl;
    return m_data_ok;
}

Profile &
ProfileRangeDB::operator[](int i)
{
    Assert(m_data_ok && !(i<0) && m_begin + i < m_end);
    return m_vec[i]; // index into local array
}

Profile
ProfileRangeDB::get(int i) const
{
    Assert( !(i<0) && m_begin + i < m_end);

    if (m_data_ok)
    {
        return const_cast<ProfileRangeDB*>(this)->operator[](i); // efficient read
    }
    else
    {
        // Read profile from database
        std::vector<Profile> vec;
        const_cast<ProfileRangeDB*>(this)->readData(vec, i, i+1);
        return vec[0];
    }
}

Profile
ProfileRangeDB::find(std::string id) const
{
    // find directly in database

    // make sure we have a connection to the database in this thread
    Assert2(m_db->connect(), "ProfileRangeDB::find(): Database connection failed" );

    std::ostringstream oss;
    oss << m_query << " AND profile_key = '" << id << "'";

    std::vector<Profile> vec;
//    bool ok = dbReadQuery(m_db, oss.str(), vec);
    bool b = dbReadQueryNoBackground(m_db, oss.str(), vec);

    Assert(b && vec.size() == 1);

    return vec[0];
}

bool
ProfileRangeVec::readMyData()
{
    m_data_ok = true;
//    cout << "ProfileRangeVec::readMyData(): " << this << " m_data_ok = " << m_data_ok << endl;
    return m_data_ok;
}

Profile
ProfileRangeVec::get(int i) const
{
    Assert(!(i<0) && m_begin + i < m_end);
    return m_vec[m_begin + i]; //index into external array
}

Profile &
ProfileRangeVec::operator[](int i)
{
    Assert(m_data_ok);
    Assert(!(i<0) && m_begin + i < m_end);
    return m_vec[m_begin + i]; //index into external array
}

Profile &
findInVec(std::vector<Profile> &vec, size_t begin, size_t end, std::string id)
{
    // linear search (slow but avoids using any more memory. If you
    // need to do lots of finds create your own index in a map)

    int size = end - begin;
    Assert (size >= 0);
    for (int i=0; i<size; ++i)
    {
        if (vec[begin + i].data().m_id == id)
        {
            return vec[begin + i];
        }
    }

    Assert2(false, "record not found in ProfileRange");
}

Profile
ProfileRangeVec::find(std::string id) const
{
    return findInVec(m_vec, m_begin, m_end, id);
}

bool
ProfileRangeRefCtdVec::readMyData()
{
    m_data_ok = true;
//    cout << "ProfileRangeRefCtdVec::readMyData(): " << this << " m_data_ok = " << m_data_ok << endl;
    return m_data_ok;
}

Profile
ProfileRangeRefCtdVec::get(int i) const
{
    Assert(!(i<0) && m_begin + i < m_end);
    return (*m_vec)[m_begin + i]; //index into external array
}

Profile &
ProfileRangeRefCtdVec::operator[](int i)
{
    Assert(m_data_ok);
    Assert(!(i<0) && m_begin + i < m_end);
    return (*m_vec)[m_begin + i]; //index into external array
}

Profile
ProfileRangeRefCtdVec::find(std::string id) const
{
    return findInVec(*m_vec, m_begin, m_end, id);
}


TEST(ProfileRangeVec)
{
    std::vector<Profile> pvec;

    for (int i=0; i<10; ++i)
    {
        std::ostringstream oss;
        oss << "P" << i;
        ProfileData pd(Identifiler, oss.str());
        Profile p(pd);
        pvec.push_back(p);
    }

    CHECK_EQUAL(10, pvec.size());
    CHECK_EQUAL("P0", pvec[0].data().m_id);
    CHECK_EQUAL("P9", pvec[9].data().m_id);

    ProfileRangeVec prvec(pvec);
    CHECK_EQUAL(10, prvec.size());

    ProfileRange pr(prvec);
    CHECK_EQUAL(10, pr.size());

    ProfileRange pr2(pr);
    CHECK_EQUAL(10, pr2.size());

    CHECK(!pr2.gotData());
    CHECK_EQUAL("P1", pr2.find("P1").data().m_id);
    CHECK(!pr2.gotData());

    pr2.readMyData();
    CHECK(pr2.gotData());

    CHECK_EQUAL("P0", pr2[0].data().m_id);
    CHECK_EQUAL("P9", pr2[9].data().m_id);

    ProfileRange pr3(pr, pr.begin(), pr.begin() + 5);
    CHECK_EQUAL(5, pr3.size());
    pr3.readMyData();
    CHECK_EQUAL("P0", pr3[0].data().m_id);
    CHECK_EQUAL("P4", pr3[4].data().m_id);

    ProfileRange pr4(pr, pr.begin() + 5, pr.begin() + 10);
    CHECK_EQUAL(5, pr4.size());
    pr4.readMyData();
    CHECK_EQUAL("P5", pr4[0].data().m_id);
    CHECK_EQUAL("P9", pr4[4].data().m_id);

}

TEST(ProfileRangeRefCtdVec)
{
    std::vector<Profile> *pvec = new std::vector<Profile>;

    for (int i=0; i<10; ++i)
    {
        std::ostringstream oss;
        oss << "P" << i;
        ProfileData pd(Identifiler, oss.str());
        Profile p(pd);
        pvec->push_back(p);
    }


    CHECK_EQUAL(10, pvec->size());
    CHECK_EQUAL("P0", (*pvec)[0].data().m_id);
    CHECK_EQUAL("P9", (*pvec)[9].data().m_id);

    ProfileRangeRefCtdVec prvec(pvec);
    CHECK_EQUAL(10, prvec.size());

    CHECK_EQUAL(1, prvec.m_vec.use_count());

    {

        ProfileRange pr(prvec);
        CHECK_EQUAL(10, pr.size());

        ProfileRange pr2(pr);
        CHECK_EQUAL(10, pr2.size());

        CHECK(!pr2.gotData());
        CHECK_EQUAL("P2", pr2.find("P2").data().m_id);
        CHECK(!pr2.gotData());

        pr2.readMyData();
        CHECK(pr2.gotData());

        CHECK_EQUAL("P0", pr2[0].data().m_id);
        CHECK_EQUAL("P9", pr2[9].data().m_id);

        ProfileRange pr3(pr, pr.begin(), pr.begin() + 5);
        CHECK_EQUAL(5, pr3.size());
        pr3.readMyData();
        CHECK_EQUAL("P0", pr3[0].data().m_id);
        CHECK_EQUAL("P4", pr3[4].data().m_id);

        ProfileRange pr4(pr, pr.begin() + 5, pr.begin() + 10);
        CHECK_EQUAL(5, pr4.size());
        pr4.readMyData();
        CHECK_EQUAL("P5", pr4[0].data().m_id);
        CHECK_EQUAL("P9", pr4[4].data().m_id);
    }

    // outside the above scope we should be back to 1 reference
    CHECK_EQUAL(1, prvec.m_vec.use_count());
}

void checkData(Profile &p)
{
    CHECK_EQUAL("BADDIES", p.data().m_dataset);
    CHECK_EQUAL("SAMPLEXX", p.data().m_sample_id);
    CHECK_EQUAL("PROFILEYY", p.data().m_profile_id);
    CHECK_EQUAL((ProfileType)Identifiler, p.data().m_kit_type);
    CHECK_EQUAL((EvidenceType)reference, p.data().m_evidence_type);
    CHECK_EQUAL(1, p.data().m_num_contributors);
    CHECK_CLOSE(0.001, p.data().m_error_rate, 1e-6);
}

void checkAlleles(Profile &p)
{
    CHECK_EQUAL(1, p.data().m_num_contributors);

    PMF<Allele> a0 = p[CSF1PO].getAllele(0);
    CHECK_EQUAL(0, a0.size()); // FF

    PMF<Allele> a1 = p[CSF1PO].getAllele(1);
    CHECK_EQUAL(1, a1.size()); // "8@.9"
    CHECK_CLOSE(0.9, a1(Allele(8)), 1e-6);

    PMF<Allele> a2 = p[FGA].getAllele(0);
    CHECK_EQUAL(0, a2.size()); // FF

    PMF<Allele> a3 = p[FGA].getAllele(1);
    CHECK_EQUAL(1, a3.size()); // "21.2@.9"
    CHECK_CLOSE(0.9, a3(Allele(21, 2)), 1e-6);

}

// set up a test population database (and destroy it afterwards)
struct Pop_fixture
{
    Pop_fixture()
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

};

TEST_FIXTURE(Pop_fixture, ProfileRangeDB)
{
    boost::shared_ptr<Database> db( new Database("test") );

    GPUDevices::initThread(ThreadData(MAIN_THREAD_NO, DEFAULT_CUDA_DEV)); // thread must be initialized before we can connect

    CHECK_EQUAL(true, db->connect());
    CHECK_EQUAL(true, db->clear());

    CHECK_EQUAL(0, db->size());

    for (int i=0; i<10; ++i)
    {
        std::ostringstream oss;
        oss << "P" << i;

        ProfileData id(Identifiler,
                       oss.str(),               // ID
                       1,                       // contributors
                       0.001,                   // delta
                       0,
                       reference,               // evidence_type
                       "BADDIES",               // dataset
                       "SAMPLEXX",              // sample_id
                       "PROFILEYY");            // profile_id

        DBProfile p(id);

        p.m_b11text[CSF1PO].push_back("F");
        p.m_b11text[CSF1PO].push_back("8@.9");

        p.m_b11text[FGA].push_back("F");
        p.m_b11text[FGA].push_back("21.2@.9");

        CHECK_EQUAL(true, db->insert(p));
    }

    CHECK_EQUAL(10, db->size());

    ProfileRangeDB prdb(db, "select * from Profiles WHERE dataset = 'BADDIES'");
    CHECK_EQUAL(10, prdb.size());

    ProfileRange pr(prdb);
    CHECK_EQUAL(10, pr.size());

    ProfileRange pr2(pr);
    CHECK_EQUAL(10, pr2.size());

    CHECK(!pr2.gotData());
    CHECK_EQUAL("P3", pr2.find("P3").data().m_id);
    CHECK(!pr2.gotData());

    pr2.readMyData();
    CHECK(pr2.gotData());

    CHECK_EQUAL("P0", pr2[0].data().m_id);
    CHECK_EQUAL("P9", pr2[9].data().m_id);

    for (int i=0; i<10; ++i)
    {
        checkData(pr2[i]);
        checkAlleles(pr2[i]);
    }

    ProfileRange pr3(pr, pr.begin(), pr.begin() + 5);
    CHECK_EQUAL(5, pr3.size());
    pr3.readMyData();
    CHECK_EQUAL("P0", pr3[0].data().m_id);
    CHECK_EQUAL("P4", pr3[4].data().m_id);

    for (int i=0; i<5; ++i)
    {
        checkData(pr3[i]);
        checkAlleles(pr3[i]);
    }

    ProfileRange pr4(pr, pr.begin() + 5, pr.begin() + 10);
    CHECK_EQUAL(5, pr4.size());

    pr4.readMyData();

    CHECK_EQUAL("BADDIES", pr4[0].data().m_dataset);

    CHECK_EQUAL("P5", pr4[0].data().m_id);
    CHECK_EQUAL("P9", pr4[4].data().m_id);

    for (int i=0; i<5; ++i)
    {
        checkData(pr4[i]);
        checkAlleles(pr4[i]);
    }

}

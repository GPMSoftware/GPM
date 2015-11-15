/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Profile.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: gareth
 */

#include "Profile.h"
#include "populationdata.h"

#include "cuda_accel/cuda_accel.h"
#include <UnitTest++/UnitTest++.h>
#include <numeric>
#include <exception>

#include "MessageStream.h"
INIT_MESSAGES("Profile")
#include "messages.h"

using namespace std;

const AlleleSet Profile::m_empty;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

Profile::Profile(ProfileData prof_data, const LocusSet &loci)
: m_data(prof_data)
, m_loci(loci)
, m_bin()
, m_hash(0)
{
}

Profile::Profile(Profile const &other)
: m_data(other.m_data)
, m_loci(other.m_loci)
, m_bin()
, m_hash(0)
{
    // NB when copying a Profile we make no attempt to copy its binary data.
    // (You could have a go at reference counting or something)
}

Profile &
Profile::operator=(Profile const &other)
{
    // see comment for copy ctor
    m_data = other.m_data;
    m_loci = other.m_loci;
    m_hash = 0;

    m_bin.clear();

    return *this;
}

std::vector<float> const&
Profile::bin(PopulationData const &popdata, SubPopModel const &spm) const
{
    pthread_mutex_lock( &mutex );

    int hash = popdata.getHash();

    if ( hash != m_hash || m_bin.empty() )
    {
        makeBin(popdata, spm);
        m_hash = hash;
    }

    Assert(! m_bin.empty());

    pthread_mutex_unlock( &mutex );

    return m_bin;
}

// create and cache a binary representation of the profile
// (this version contains all the loci in the background and sizes can vary)
//
void
Profile::makeBin(
        PopulationData const &popdata,
        SubPopModel const &spm) const
{
	info << "Profile::makeBin(): m_name = " << m_data.m_id << endl;

	LocusInfo const &loc_info = popdata.getLocusInfo();

	m_bin.clear();

	LocusSet::const_iterator ip1 = m_loci.begin();

	Assert2((int)popdata.numLoci() <= cuda_num_loci, "Profile::makeBin(): too many loci");

	for (int i=0; i<cuda_num_loci; ++i)
	{
        int locus_size = loc_info.locus_size[i];

		if (popdata.hasLocus(i))   // locus exists in the population database
		{
			Locus locus = locus_none;
			if(ip1 != m_loci.end())
			{
				locus = ip1->first;
			}

			if (locus == (Locus)i) // Profile contains locus i
			{
				info << "Profile::makeBin(): adding locus " << locus << endl;

				const AlleleSet &as = ip1->second;
				ip1++;

				if (locus != AMEL) // ignore Amelogenin (write -1)
				{
					HMatrix h;
					if (spm.type == SubPopModel::HW)
					{
                        h = as.getMatchHMatrix(popdata(locus), m_data, false); // sparse = false
					}
					else // non-HW
                    {
					    PMF<Allele> back = popdata(locus);
                        BackFreq hback(back, spm);

                        h = as.getMatchHMatrixNHW(hback, m_data, false); // sparse = false
                    }
					const std::vector<float>& hvec = h.vec();
					int hvec_size = hvec.size();

                    Assert(hvec_size == locus_size);

                    for (int j=0; j<locus_size; ++j)
					{
                        m_bin.push_back(hvec[j]);
					}

					continue; // we are done: don't write -1s

				} // else AMEL, write -1
				else
				{
					info << "Profile::makeBin(): AMEL is not used: " << i << endl;
				}

			} // else locus is missing from profile, write -1
			else
			{
				info << "Profile::makeBin(): locus is missing from profile: " << i << endl;
			}

		} // else locus is missing from population database, write -1
		else
		{
			info << "Profile::makeBin(): locus is missing from population database: " << i << endl;
		}

		info << "Profile::makeBin(): ignoring locus " << i << endl;

		// for whatever reason, ignore this locus (fill with -1)
		for (int j=0; j<locus_size; ++j)
		{
            m_bin.push_back(-1);
		}
	}

    Assert2((int)m_bin.size() == loc_info.profile_size, "error in Profile::makeBin()");
}

// reset all cached state
void
Profile::reset() const
{
    resetH();
    m_bin.clear();
}

void
Profile::resetH() const
{
	for (LocusSet::const_iterator ip1 = m_loci.begin();
		ip1 != m_loci.end();
		++ip1)
	{
		const AlleleSet &as = ip1->second;
        as.setData(m_data);
	}
}

AlleleSet&
Profile::operator[](int loc)
{
    // reset just the changed locus
	m_loci[(Locus)loc].reset();

	// m_bin will need to be recalculated
    m_bin.clear();

    return m_loci[(Locus)loc];
}

const AlleleSet&
Profile::operator[](int locus) const
{
    LocusSet::const_iterator it = m_loci.find((Locus)locus);
    if (it == m_loci.end())
    {
        // that's ok - not every locus need be present
        return m_empty;
    }
    else
    {
        return it->second;
    }
}

void
Profile::setData(ProfileData const &d)
{
    pthread_mutex_lock( &mutex );

    m_data = d;
    reset();

    pthread_mutex_unlock( &mutex );
}

void
Profile::setErrorRate(double d)
{
    if (m_data.m_error_rate == d)
    {
        return; //someone already set it
    }

    pthread_mutex_lock( &mutex );

    m_data.m_error_rate = d;
    reset();

    pthread_mutex_unlock( &mutex );
}

void
Profile::setMutRate(double d)
{
    if (m_data.m_mutation_rate == d)
    {
        return; //someone already set it
    }

    pthread_mutex_lock( &mutex );

    m_data.m_mutation_rate = d;
    reset();

    pthread_mutex_unlock( &mutex );
}

void
Profile::setNumContributors(int n)
{
    pthread_mutex_lock( &mutex );

    m_data.m_num_contributors = n;
    reset();

    pthread_mutex_unlock( &mutex );

}

TEST(Profile)
{
	ProfileData pd(Identifiler, "P01");
	Profile p(pd);
	CHECK_EQUAL(Identifiler, p.data().m_kit_type);
	CHECK_EQUAL("P01", p.data().m_id);
	CHECK(p.empty());

	p[FGA] = AlleleSet(Allele(26, 2), Allele(26, 2));
//	p.m_loci.insert(std::make_pair(FGA, AlleleSet(Allele(26, 2), Allele(26, 2))));
}

enum TestLoci
{
  LOCUS1 = 0
, LOCUS2
, LOCUS3
, LOCUS4
, test_size
};

enum TestAlleles
{
  F = 0 // unknown
, A = 1
, B
, C
, D
};

static PMF<Allele>::POD LOCUS1_FREQ[] =
{
  { Allele(A), 0.1 } // Allele() optional
, { B, 0.2 }
, { C, 0.3 }
, { D, 0.4 }
};

struct SpreadsheetFixtureProf
{
	SpreadsheetFixtureProf()
	: id1(test, "P01")
	, id2(test, "P02")
	, p1(id1)
	, p2(id2)
	, bg_locus1(LOCUS1_FREQ, sizeof(LOCUS1_FREQ)/sizeof(PMF<Allele>::POD))
	{
		popdata.setLocus(LOCUS1, bg_locus1); // (all the same)
		popdata.setLocus(LOCUS2, bg_locus1); //
		popdata.setLocus(LOCUS3, bg_locus1); //
		popdata.setLocus(LOCUS4, bg_locus1); //
	}

	ProfileData id1;
	ProfileData id2;
	Profile p1;
	Profile p2;
	PMF<Allele> bg_locus1;
	PopulationData popdata;
};


TEST_FIXTURE(SpreadsheetFixtureProf, profile_1)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	p1[LOCUS2] = AlleleSet(A, B);
	p2[LOCUS2] = AlleleSet(A, B);

	p1[LOCUS3] = AlleleSet(C, C);
	p2[LOCUS3] = AlleleSet(C, C);

	p1[LOCUS4] = AlleleSet(C, D);
	p2[LOCUS4] = AlleleSet(C, D);

//	const float *b = p1.bin(popdata);
//  CHECK((float*)0 != b);

//	CachedArray<float>::ConstTempRef b = p1.binRef(popdata);

	vector<float> const &b = p1.bin(popdata);

	LocusInfo const &loc_info = popdata.getLocusInfo();
//    CHECK_EQUAL(loc_info.profile_size, b.numValues());
    CHECK_EQUAL(loc_info.profile_size, b.size());

	// each allele has 10 entries. One is non-zero:
	// AA AB AC AD BB BC BD CC CD DD
	// 0  1  2  3  4  5  6  7  8  9
	//
	// p1 is {AA, AB, CC, CD} --> { 0, 1, 7, 8 }

//	LocusInfo const &loc_info = popdata.getLocusInfo();

//	CHECK_EQUAL(1.0, b[0  + 0]);
//	CHECK_EQUAL(1.0, b[1 * cuda_locus_size + 1]);
//	CHECK_EQUAL(1.0, b[2 * cuda_locus_size + 7]);
//	CHECK_EQUAL(1.0, b[3 * cuda_locus_size + 8]);
    CHECK_EQUAL(1.0, b[loc_info.locus_offset[0] + 0]);
    CHECK_EQUAL(1.0, b[loc_info.locus_offset[1] + 1]);
    CHECK_EQUAL(1.0, b[loc_info.locus_offset[2] + 7]);
    CHECK_EQUAL(1.0, b[loc_info.locus_offset[3] + 8]);

	// other entries are zero so sum is 4
//	float sum = std::accumulate(b, b + 4*cuda_locus_size, 0.0);
//  float sum = std::accumulate(b, b + loc_info.profile_size, 0.0);
//    float const *p = b.values();
    float const *p = &b[0];
    float sum = std::accumulate(p, p + loc_info.profile_size, 0.0);
	CHECK_EQUAL(4, sum);

}

#if 0
TEST_FIXTURE(SpreadsheetFixtureProf, when_is_memory_freed)
{
    // see if memory that is freed is visible in the system monitor.
    // Answer: For c-array of char it is visible as soon as it is deleted.
    //         For vector<Profile> the memory is not visible until the program exits
    //         Same for vector<char>: vectors do not give up their memory unless resized.

    p1[LOCUS1] = AlleleSet(A, A);
    p1[LOCUS2] = AlleleSet(A, B);
    p1[LOCUS3] = AlleleSet(C, C);
    p1[LOCUS4] = AlleleSet(C, D);

    long n = 1000000000;

//    char *p = new char[n]; // 1GB
//    std::vector<Profile> p;
    std::vector<char> p;

    cout << "Array allocated. filling it up with " << n << " objects ... " << std::endl;

    for (int i=0; i<n; ++i)
    {
//        p[i] = 'X';
//        p.push_back(p1);
//        p[i].makeBin();
        p.push_back('X');
    }

    // wait 10 seconds
    cout << "done. Now sleep for 10 seconds ... " << std::flush;
    sleep(10);
    cout << "awake!" << std::endl;

//    delete [] p;
    p.clear();
    p.reserve(0); // no effect seen in system monitor

    // wait 10 seconds
    cout << n << " objects deleted. Sleeping for 10 seconds ... " << std::flush;
    sleep(10);
    cout << "awake!" << std::endl;
}

TEST(vector_size)
{
    // how big is a vector of a million floats?

    Timer t;
    MemMeter m;

    long m0 = m.read();

    cout << "m0  = : " << m << endl;

    std::vector<float> vec;

    for (int i=0; i<1000000; ++i)
    {
        vec.push_back((float)i);
    }

    t.stop();
    long m1 = m.read();

    cout << "vector of a million floats: " << endl;
    cout << "time     = : " << t << endl;
    cout << "m1       = : " << m << endl;
    cout << "m1 - m0  = : " << m1 - m0 << endl;

}
#endif

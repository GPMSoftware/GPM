/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * population.cpp
 *
 *  Created on: Dec 7, 2009
 *      Author: gareth
 *
 *      Interface to population databases
 */

#include "populationdata.h"
#include "loci.h"
#include "Profile.h"
#include "ProfileFilter.h"
#include <UnitTest++/UnitTest++.h>

#include "MessageStream.h"
INIT_MESSAGES("populationdata")
#include "messages.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

using namespace std;

pthread_mutex_t mutex2 = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;

// Default population database
static PopulationData default_pop;


PopulationData::PopulationData(Unknowns policy)
: m_sample_size(0)
, m_unknowns_policy(policy)
{
    if (policy == read_from_env)
    {
        getPolicyFromEnv();
    }

    m_hash = 0; // all empty PopulationData have m_hash = 0
}

PopulationData::PopulationData(string const &path, Unknowns policy)
: m_sample_size(0)
, m_unknowns_policy(policy)
{
    if (policy == read_from_env)
    {
        getPolicyFromEnv();
    }

    read(path);
}

void
PopulationData::getPolicyFromEnv()
{
    m_unknowns_policy = ignore;
    std::string unknown_alleles = getStringEnv("UNKNOWN_ALLELES");

    if (unknown_alleles == "IGNORE")
    {
        m_unknowns_policy = ignore;
    }
    else if (unknown_alleles == "ADD_AS_RARE")
    {
        m_unknowns_policy = add_as_rare;
    }
    else if (unknown_alleles != "")
    {
        warn << startl << "UNKNOWN_ALLELES: should be \"IGNORE\" or \"ADD_AS_RARE\" (treating as IGNORE)" << endl;
    }
}

int
PopulationData::numLoci() const
{
	return m_data.size();
}

bool
PopulationData::normalize()
{
    pthread_mutex_lock( &mutex2 );

    bool ret = true;
    std::map< Locus, PMF<Allele> >::iterator it;
    for (it = m_data.begin(); it != m_data.end(); ++it)
    {
        ret &= it->second.normalize();
    }
    // NB we do NOT reHash() as a result of normalize(). A population database should
    // always be normalized, so it is not considered to have changed.
    pthread_mutex_unlock( &mutex2 );

    return ret;
}

bool
PopulationData::hasLocus(int const &loc) const
{
    pthread_mutex_lock( &mutex2 );

    bool ret = true;
	std::map< Locus, PMF<Allele> >::const_iterator it;
	if (loc<0 || loc>=num_loci || (it=m_data.find((Locus)loc)) == m_data.end())
	{
		ret = false;
	}
    pthread_mutex_unlock( &mutex2 );

	return ret;

}

const PMF<Allele>&
PopulationData::operator()(int const &loc) const
{
    pthread_mutex_lock( &mutex2 );

    const PMF<Allele> *ret = 0;

    std::map< Locus, PMF<Allele> >::const_iterator it;
    if (loc<0 || loc>=num_loci || (it=m_data.find((Locus)loc)) == m_data.end())
    {
        // throw std::out_of_range("PopulationData: locus not found");
        Assert("PopulationData: locus not found");
    }
    else
    {
        ret = &(it->second);
    }

    pthread_mutex_unlock( &mutex2 );

    if (ret == 0)
    {
        throw std::out_of_range("PopulationData: locus not found");
    }
	return *ret;
}

void
PopulationData::setLocus(int const &loc, PMF<Allele> const &pmf)
{
    pthread_mutex_lock( &mutex2 );

    m_data[(Locus)loc] = pmf;
    m_locus_info.update(*this);
    reHash();

    pthread_mutex_unlock( &mutex2 );
}

const LocusInfo&
PopulationData::getLocusInfo() const
{
    Assert2(! m_data.empty(), "PopulationData::getLocusInfo: no data");
    return m_locus_info;
}

std::string
PopulationData::path() const
{
    return m_path;
}

double
PopulationData::getFrequency(int const &loc, Allele const &a, int *added_total)
{
    Assert2(! m_data.empty(), "PopulationData::getFrequency: no data");

	if (a == Allele(Allele::unknown))
	{
		return 0;
	}

	double f;
	int added = 0;
	double prare = 0;

    pthread_mutex_lock( &mutex2 );

	if (m_sample_size > 0)
	{
		prare = 5.0 / (2.0 * m_sample_size);
	}

	PMF<Allele> &background = m_data[(Locus)loc];
	PMF<Allele>::iterator it = background.find(a);
	if (it != background.end())
	{
		f = it->second;
	}
	else
	{
		// allele is not in the database
		if (m_unknowns_policy == add_as_rare && prare > 0)
		{
			// add to background at prare
			warn << startl << "found unknown allele " << locus_name[loc] << ":" << a.string() << " (adding to background: prare = " << prare << ")" << endl;

			// the first time we need to change anything, grab the mutex
			if (!added)
			{
//			    pthread_mutex_lock( &mutex2 );
			}

			f = prare;
			background.insert(make_pair(a, f));
			background.normalize();
			++added;
		}
		else if (m_unknowns_policy == ignore)
		{
			warn << startl << "unknown allele frequency requested " << locus_name[loc] << ":" << a.string() << " (returning 0)" << endl;
			f = 0;
		}
        else
        {
            Assert2(false, "m_unknowns_policy has invalid value");
        }

		if ((int)background.size() > CUDA_MAX_ALLELES)
		{
			error << alignl << "Fatal: CUDA_MAX_ALLELES exceeded!" << endl;
            cout << "Fatal: CUDA_MAX_ALLELES exceeded!" << endl;
            exit(1);
		}
	}

	if (added)
	{
	    m_locus_info.update(*this);

	    reHash();
//	    pthread_mutex_unlock( &mutex2 );

	    if (added_total)
	    {
	        *added_total += added;
	    }
	}

    pthread_mutex_unlock( &mutex2 );

	return f;
}

void
PopulationData::clear()
{
    pthread_mutex_lock( &mutex2 );

	m_data.clear();
    m_locus_info.update(*this);

    m_hash = 0; // all empty PopulationData have m_hash = 0
    m_path = "";

    pthread_mutex_unlock( &mutex2 );
}

// read data from files
bool
PopulationData::read(string const &path)
{
	pthread_mutex_lock( &mutex2 );

	// Read data for each locus from the specified path.
	// NB we assume the files in path correspond one-to-one
	// with the elements of enum Locus

	m_data.clear();

	double lowest_freq = 1.0;

	for (int i=0; i< num_loci; ++i)
	{
		string filename = path + "/" + popdata_file[i];
		info << startl << "opening " << filename << endl;

		ifstream ifs_locus(filename.c_str());

		// if failed to open file skip it (locus may not be needed)
		if (! ifs_locus.good())
		{
//		    cout << "PopulationData::read(): path = " << path << " filename = " << filename << " : file not found" << endl;
			continue;
		}

		PMF<Allele> a;
		double freq;
		Allele x;

		while (ifs_locus >> x >> freq)
		{
            a.insert(make_pair(x, freq));
			if (lowest_freq > freq) lowest_freq = freq;
		}

		// make sure the data is normalized
		a.normalize();
		m_data[(Locus)i] = a;
	}

	// look for a file telling us the sample size
	string filename = path + "/sample_size";
	info << startl << "opening " << filename << endl;

	ifstream sampsize(filename.c_str());

	if (sampsize.good())
	{
		sampsize >> m_sample_size;
		info << alignl << "read m_sample_size from file: " << m_sample_size << endl;
	}
	else
	{
		// assume lowest_freq == 1/2n
//		m_sample_size = (int)(0.5/lowest_freq);
//		warn << startl << "No sample_size file: setting m_sample_size = 1/2*lowest_freq: " << m_sample_size << endl;

		// Assume 300. This is a common size for frequency databases.
		// (the above method may give too small an estimate for m_sample_size which
		// results in a rare allele frequency that is too large, bringing down the other frequencies and overestimating the LR)
		m_sample_size = 300;
		warn << startl << "No sample_size file: setting m_sample_size = " << m_sample_size << endl;
	}

	m_locus_info.update(*this);
	reHash();

	pthread_mutex_unlock( &mutex2 );

	// if we read any data at all we are happy
	m_path = path;
	return ! m_data.empty();
}

// write data to files
bool
PopulationData::write(string const &path) const
{
    pthread_mutex_lock( &mutex2 );

    // check path exists - if not try to create it
    struct stat stbuf;
    if (stat(path.c_str(), &stbuf) != 0) // 0 means success
    {
        if (mkdir(path.c_str(), 0777) != 0) // 0 means success
        {
            error << startl << "unable to create " << path << endl;
            return false;
        }
        else
        {
            info << startl << "created " << path << endl;
        }
    }
    else
    {
        // check it is a directory
        if (!S_ISDIR(stbuf.st_mode))
        {
            error << startl << path << " already exists and is not a directory" << endl;
            return false;
        }
    }

    // delete any files in the directory (ensure there is no stale data)
    DIR *dfd = opendir(path.c_str());
    struct dirent *dp;
    while ((dp = readdir(dfd)) != NULL)
//    while ((readdir(dfd, dp)) != NULL)
    {
        if (strcmp(dp->d_name, ".") == 0
                || strcmp(dp->d_name, "..") == 0 ) continue;

        string filename = path + "/" + dp->d_name;

        if (unlink(filename.c_str()) != 0) // 0 means success
        {
            warn << startl << "unable to delete stale file: " << filename << endl;
        }
    }

    map< Locus, PMF<Allele> >::const_iterator it;
    for (it = m_data.begin(); it != m_data.end(); ++it)
    {
        int loc = it->first;
        PMF<Allele> pmf = it->second;
        string filename = path + "/" + popdata_file[loc];

        if (pmf.empty()) continue; // no data for this locus - continue

        info << startl << "opening " << filename << endl;

        ofstream ofs_locus(filename.c_str());

        // if failed to open file issue warning
        if (! ofs_locus.good())
        {
            warn << startl << "unable to update stale file: " << filename << endl;
            continue;
        }

        PMF<Allele>::const_iterator pit;
        for (pit = pmf.begin(); pit != pmf.end(); ++pit)
        {
            ofs_locus << pit->first << " " << pit->second << "\n";
        }

        ofs_locus.close();
    }

    // write the sample size
    string filename = path + "/sample_size";
    info << startl << "opening " << filename << endl;

    ofstream sampsize(filename.c_str());

    if (sampsize.good())
    {
        sampsize << m_sample_size << "\n";
    }
    sampsize.close();

    pthread_mutex_unlock( &mutex2 );

    return true;
}

// return the next global PopulationData index.
// Call this function whenever the PopulationData changes
// ALWAYS call from within a mutex lock
void
PopulationData::reHash()
{
    static long int m_index = 1; // start at 1: value will never be 0
    m_hash = ++m_index;

    info << startl << "reHash(): " << m_hash << endl;
}

int
PopulationData::sampleSize() const
{
    return m_sample_size;
}

void
PopulationData::setSampleSize(int n)
{
    pthread_mutex_lock( &mutex2 );

    m_sample_size = n;
    reHash();

    pthread_mutex_unlock( &mutex2 );
}

long int
PopulationData::getHash() const
{
    return m_hash;
}

void
PopulationData::setPolicy(Unknowns policy)
{
    pthread_mutex_lock( &mutex2 );

    m_unknowns_policy = policy;
    reHash();

    pthread_mutex_unlock( &mutex2 );
}

Unknowns
PopulationData::getPolicy() const
{
    return m_unknowns_policy;
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
  A = 1 // not 0! that means unknown
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

static PMF<Allele>::POD LOCUS2_FREQ[] =
{
  { A, 0.1 }
, { B, 0.9 }
};

static PMF<Allele>::POD LOCUS3_FREQ[] =
{
  { C, 0.9 }
, { D, 0.1 }
};

static PMF<Allele>::POD LOCUS4_FREQ[] =
{
  { A, 0.5 }
, { D, 0.5 }
};


PopulationData &
testPopulationData()
{
	static PopulationData pop;

	if (pop.numLoci() == 0)
	{
		pop.setLocus(LOCUS1, PMF<Allele>(LOCUS1_FREQ, sizeof(LOCUS1_FREQ)/sizeof(PMF<Allele>::POD)));
		pop.setLocus(LOCUS2, PMF<Allele>(LOCUS2_FREQ, sizeof(LOCUS2_FREQ)/sizeof(PMF<Allele>::POD)));
		pop.setLocus(LOCUS3, PMF<Allele>(LOCUS3_FREQ, sizeof(LOCUS3_FREQ)/sizeof(PMF<Allele>::POD)));
		pop.setLocus(LOCUS4, PMF<Allele>(LOCUS4_FREQ, sizeof(LOCUS4_FREQ)/sizeof(PMF<Allele>::POD)));
	}

	return pop;
}

PopulationData const & knownLoci()
{
    static PopulationData pop;

    if (pop.numLoci() == 0)
    {
        for (int i=0; i<num_loci; ++i)
        {
            pop.setLocus((Locus)i, PMF<Allele>());
        }
    }

    return pop;
}

void popSet(const PopulationData &newpop)
{
    pthread_mutex_lock( &mutex2 ); // TODO might be better to put lock in operator=()

	default_pop = newpop;

    pthread_mutex_unlock( &mutex2 );

}

void popClear()
{
	default_pop.clear();
	default_pop.setSampleSize(0);
}

void printPop(const PopulationData &pop)
{
	static int n = 0;

	pthread_mutex_lock( &mutex2 );

	ostringstream filename;
	filename  << "pop" << ++n;

	ofstream ftmp(filename.str().c_str());

	for (int i=0; i< pop.numLoci(); ++i)
	{
		ftmp << i << ") " << pop(i) << endl;
	}

	pthread_mutex_unlock( &mutex2 );
}

PopulationData const &
populationData()
{
    static bool default_has_been_read = false;

	if (default_pop.numLoci() == 0 && ! default_has_been_read)
	{
		// no default population database has been set - read from files
		string path = getStringEnv("POPDATA");

		Assert2 (path.size(), "No population database set (environment variable POPDATA)");
		Assert2 (default_pop.read(path), "Population database (environment variable POPDATA) can't be read");

		info << startl << "populationData(): default_pop = " << default_pop.path() << " numLoci() == " << default_pop.numLoci() << endl;

		default_has_been_read = true;
	}

	Assert2(default_pop.numLoci() != 0, "populationData(): population database has been corrupted!");

	return default_pop;
}

bool
operator==(PopulationData const &p1, PopulationData const &p2)
{
    for (int i=0; i< num_loci; ++i)
    {
        if (p1.hasLocus(i) && p2.hasLocus(i))
        {
            if (p1(i) != p2(i))
            {
                cout << "p1("<< i << ") = " << p1(i) << endl;
                cout << "p2("<< i << ") = " << p2(i) << endl;
                return false;
            }
        }
        else if (! p1.hasLocus(i) && ! p2.hasLocus(i))
        {
            continue; // ok
        }
        else
        {
            return false;
        }
    }

    return true;
}

void
LocusInfo::update(PopulationData const &p)
{
    // everything is calculated from num_alleles
    for (int i=0; i<cuda_num_loci; ++i)
    {
        num_alleles[i] = p.hasLocus(i) ? p(i).size() : 0;
    }

    for (int i=0; i<cuda_num_loci; ++i)
    {
        locus_size[i]   = num_alleles[i] * (num_alleles[i] + 1) / 2;
        locus_offset[i] = (i == 0) ? 0 : locus_offset[i-1] + locus_size[i-1];
        back_offset[i]  = (i == 0) ? 0 : back_offset[i-1] + num_alleles[i-1];
    }

    int i_back = cuda_num_loci - 1;

    profile_size = locus_offset[i_back] + locus_size[i_back];
    back_size = back_offset[i_back] + num_alleles[i_back];
}

ostream &
operator<<(ostream &os, const CudaLocusInfo& loc_info)
{
    os << "num_alleles   : ";
    for (int i=0; i<cuda_num_loci; ++i)
    {
        os << setw(5) << loc_info.num_alleles[i] << " ";
    }

    os << "\nlocus_size  : ";
    for (int i=0; i<cuda_num_loci; ++i)
    {
        os << setw(5) << loc_info.locus_size[i] << " ";
    }

    os << "\nlocus_offset: ";
    for (int i=0; i<cuda_num_loci; ++i)
    {
        os << setw(5) << loc_info.locus_offset[i] << " ";
    }

    os << "\nback_offset : ";
    for (int i=0; i<cuda_num_loci; ++i)
    {
        os << setw(5) << loc_info.back_offset[i] << " ";
    }

    os << "\nprofile_size: " << loc_info.profile_size;
    os << "\nback_size   : " << loc_info.back_size << endl;

    return os;
}

TEST(population_data1)
{
	PMF<Allele> pmf(LOCUS1_FREQ, sizeof(LOCUS1_FREQ)/sizeof(PMF<Allele>::POD));

	pmf.normalize();

	CHECK_EQUAL(0.1, pmf[A]);
	CHECK_EQUAL(0.2, pmf[B]);
	CHECK_EQUAL(0.3, pmf[C]);
	CHECK_EQUAL(0.4, pmf[D]);
}

TEST(population_data2)
{
	PopulationData testKitData = testPopulationData();

	CHECK_EQUAL(4, testKitData.numLoci());
}

TEST(population_data3)
{
	PopulationData identifilerData = testPopulationData();

	PMF<Allele>::const_iterator ipmf;
	CHECK((ipmf = identifilerData(A).find(1)) != identifilerData(A).end());
	CHECK_EQUAL(ipmf->second, 0.1);
}

char *testdb = "test_pop_data";

TEST(population_data4) // read and write
{
    PopulationData testpop = testPopulationData();

    string cmd = string("rm -rf ") + testdb;
    system(cmd.c_str());

    // create directory and write files
    CHECK(testpop.write(testdb));

    // overwrite files
    CHECK(testpop.write(testdb));

    PopulationData mypop;

    // read
    CHECK(mypop.read(testdb));
    CHECK_EQUAL(4, mypop.numLoci());
    CHECK_EQUAL(true, mypop.hasLocus(0));
    CHECK_EQUAL(true, mypop.hasLocus(1));
    CHECK_EQUAL(true, mypop.hasLocus(2));
    CHECK_EQUAL(true, mypop.hasLocus(3));
    CHECK_EQUAL(false, mypop.hasLocus(4));

    CHECK(testpop == mypop);

    system(cmd.c_str());
}


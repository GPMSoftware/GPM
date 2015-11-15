/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * util.cpp
 *
 *  Created on: Dec 7, 2009
 *      Author: gareth
 */

#include "util.h"
#include "Assert.h"
#include "cuda_accel/cuda_accel.h"
#include "cuda_accel/GPUDevices.h"
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <UnitTest++/UnitTest++.h>

const pthread_t main_thread_id = pthread_self();

bool isMainThread()
{
    return pthread_self() == main_thread_id;
}

double drand()
{
	return double(rand()) / RAND_MAX;
}

int rand(int range)
{
	return int(drand() * range);
}

int rand(int start, int end)
{
	return start + rand(end - start);
}

double Rand::drand()
{
    return double(rand_r(&m_state)) / RAND_MAX;
}

int Rand::rand(int range)
{
    return int(drand() * range);
}

int Rand::rand(int start, int end)
{
    return start + rand(end - start);
}

TEST(Rand)
{
    Rand r1(0), r2(42);

    const int TSIZE = 100;
    double runA_1[TSIZE], runA_2[TSIZE];

    for (int i=0; i< TSIZE; ++i)
    {
        runA_1[i] = r1.drand();
    }

    for (int i=0; i< TSIZE; ++i)
    {
        runA_2[i] = r2.drand();
    }

// They are not the same (by inspection)
//    for (int i=0; i< TSIZE; ++i)
//    {
//        CHECK_EQUAL(runA_1[i], runA_2[i]);
//    }

    // check we get the same results again
    Rand s1(0), s2(42);
    double runB_1[TSIZE], runB_2[TSIZE];
    for (int i=0; i< TSIZE; ++i)
    {
        runB_1[i] = s1.drand();
    }

    for (int i=0; i< TSIZE; ++i)
    {
        runB_2[i] = s2.drand();
    }

    for (int i=0; i< TSIZE; ++i)
    {
        CHECK_EQUAL(runA_1[i], runB_1[i]);
        CHECK_EQUAL(runA_2[i], runB_2[i]);
    }

    // now interleave them and check again
    Rand t1(0), t2(42);
    for (int i=0; i< TSIZE; ++i)
    {
        runB_1[i] = t1.drand();
        runB_2[i] = t2.drand();
    }

    for (int i=0; i< TSIZE; ++i)
    {
        CHECK_EQUAL(runA_1[i], runB_1[i]);
        CHECK_EQUAL(runA_2[i], runB_2[i]);
    }
}

std::ostream& operator<<(std::ostream &os, const Timer& timer)
{
	os.setf(std::ios_base::fixed, std::ios_base::floatfield);
	os << std::setprecision(2) << timer.read();
	os.setf((std::_Ios_Fmtflags)0, std::ios_base::floatfield);
	return os;
}

std::ostream& operator<<(std::ostream &os, const MemMeter& m)
{
	int mem = m.read();
	if (mem >= 0)
	{
		os << m.read() << " KiB";
	}
	else
	{
		os << "error";
	}
	return os;
}

using namespace std;

// parse a list of tokens into a vector of strings
// repeated separators are treated as a single separator  ("A,,B" --> "A", "B")
vector<string>
split(string const &s, string separator)
{
	vector<string> ret;
	size_t i = 0, j = 0;

	i = s.find_first_not_of(separator, j);
	while(i != string::npos)
	{
		j = s.find_first_of(separator, i);
		ret.push_back(s.substr(i, j-i));
		i = s.find_first_not_of(separator, j);
	}

	return ret;
}

// parse a list of tokens into a vector of strings
// repeated separators denote empty strings ("A,,B" --> "A", "", "B")
vector<string>
split2(string const &s, string separator)
{
	vector<string> ret;
	size_t i = 0, j = 0, n = s.size();

	do
	{
		j = s.find_first_of(separator, i);
		if (j == string::npos) j = n;
		ret.push_back(s.substr(i, j-i));
		i = j+1;
	}
	while (i <= n);

	return ret;
}

// remove leading and trailing spaces from vector of strings
void cleanup(vector<string> &words)
{
	size_t n = words.size();

	for (size_t i=0; i<n; ++i)
	{
		cleanup(words[i]);
	}
}

// remove leading and trailing spaces from string
void cleanup(string &word)
{
    size_t start = word.find_first_not_of(' ');
    if (start == string::npos)
    {
        word.clear();
        return;
    }

    size_t end   = word.find_last_not_of(' ');
    Assert(end != string::npos);
    Assert(end >= start);

    size_t newlen = end-start+1;

    if (newlen<word.size())
    {
        word = word.substr(start, newlen);
    }
}

// suggestion: set a break point in this function
void exit_debug (int status)
{
	exit(status);
}

// get float environment variable (0 if not set)
float getFloatEnv(const char *env_name, float default_val)
{
	const char *env = getenv(env_name);
	if (!env) return default_val;

	float ret;
	istringstream iss(env);
	if ( !(iss >> ret) || !iss.eof())
	{
		ret = default_val;
	}
	return ret;
}

// get float environment variable (0 if not set)
int getIntEnv(const char *env_name, int default_val)
{
    const char *env = getenv(env_name);
    if (!env) return default_val;

    int ret;
    istringstream iss(env);
    if ( !(iss >> ret) || !iss.eof())
    {
        ret = default_val;
    }
    return ret;
}

// get bool environment variable
// returns the default value unless explicitly set to the other value
// (ensures that set to "" behaves the same as unset)
bool getBoolEnv(const char *env_name, bool default_val)
{
	const char *env = getenv(env_name);
	if (!env) return default_val;

	const char *non_default_string = default_val ? "false" : "true";

	if (strcmp(env, non_default_string) == 0)
	{
		return ! default_val;
	}

	return default_val;
}

// get string environment variable ("" if not set)
const char * getStringEnv(const char *env_name, const char *default_val)
{
	const char *env = getenv(env_name);
    return env ? env : default_val;
}

std::string
baseName(std::string const &filename)
{
    size_t p1 = filename.rfind('/');
    if (p1 == std::string::npos)
    {
        p1 = 0;
    }
    else
    {
        ++p1;
    }

    size_t p2 = filename.find('.', p1+1);

    return filename.substr(p1, p2 - p1);
}

std::string
path(std::string const &filename)
{
    size_t p1 = filename.rfind('/');
    if (p1 == std::string::npos)
    {
        p1 = 0;
    }
    else
    {
        ++p1;
    }

    return filename.substr(0, p1);
}

std::string
dateStr()
{
    time_t rawtime;
    struct tm *timeinfo;
    const int size = 80;
    char buffer[size];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime (buffer, size, "%d%m%y%H%M", timeinfo);

    return std::string(buffer);
}

const char *
homeDir()
{
    return getStringEnv("HOME");
}

std::string
homeDirStr()
{
    return homeDir();
}

TEST(path)
{
    string p = path("/home/me/stuff/x.txt");
    CHECK_EQUAL("/home/me/stuff/", p);
}

TEST(split2)
{
	vector<string> args;

	args = split2("");
	CHECK_EQUAL(1, args.size());
	CHECK_EQUAL("", args[0]);

	args = split2(",");
	CHECK_EQUAL(2, args.size());
	CHECK_EQUAL("", args[0]);
	CHECK_EQUAL("", args[1]);

	args = split2("a");
	CHECK_EQUAL(1, args.size());
	CHECK_EQUAL("a", args[0]);

	args = split2(",a");
	CHECK_EQUAL(2, args.size());
	CHECK_EQUAL("", args[0]);
	CHECK_EQUAL("a", args[1]);

	args = split2("a,");
	CHECK_EQUAL(2, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("", args[1]);

	args = split2(",a,");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("", args[0]);
	CHECK_EQUAL("a", args[1]);
	CHECK_EQUAL("", args[2]);

	args = split2("a,bb,ccc");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);

	args = split2("a,bb,,ccc");
	CHECK_EQUAL(4, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("", args[2]);
	CHECK_EQUAL("ccc", args[3]);

	args = split2("a,,,ccc,");
	CHECK_EQUAL(5, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("", args[1]);
	CHECK_EQUAL("", args[2]);
	CHECK_EQUAL("ccc", args[3]);
	CHECK_EQUAL("", args[4]);
}

TEST(split)
{
	vector<string> args;

	args = split("");
	CHECK_EQUAL(0, args.size());

	args = split(" ");
	CHECK_EQUAL(0, args.size());

	args = split("a");
	CHECK_EQUAL(1, args.size());
	CHECK_EQUAL("a", args[0]);

	args = split(" a ");
	CHECK_EQUAL(1, args.size());
	CHECK_EQUAL("a", args[0]);

	args = split("a bb ccc");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);

	args = split("a bb   ccc ");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);

	args = split(" a bb   ccc ");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);

	args = split(" a bb   ccc  ");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);

	args = split("a,bb,ccc", ",");
	CHECK_EQUAL(3, args.size());
	CHECK_EQUAL("a", args[0]);
	CHECK_EQUAL("bb", args[1]);
	CHECK_EQUAL("ccc", args[2]);
}

TEST(mapInv)
{
    map<int, string> m;
    m[0] = "dog";
    m[1] = "dog";
    m[3] = "cat";
    m[6] = "pig";
    multimap<string, int> mi;
    mapInv(m, mi);

    CHECK_EQUAL(2, mi.count("dog"));
    pair<multimap<string, int>::const_iterator, multimap<string, int>::const_iterator> its = mi.equal_range("dog");

    // NB repeated key values appear in insert order
    multimap<string, int>::const_iterator it = its.first;
    CHECK(it != its.second);
    CHECK_EQUAL(0, it->second); // dog --> 0
    CHECK(++it != its.second);
    CHECK_EQUAL(1, it->second); // dog --> 1
    CHECK(++it == its.second);

    CHECK_EQUAL(1, mi.count("cat"));
    CHECK_EQUAL(1, mi.count("pig"));
}

TEST(memmeter)
{
    MemMeter m;
    long mem_total = m.total();
    CHECK(mem_total > 0);
    cout << "System memory = " << mem_total << endl;
}

#if 0
TEST(speed_test)
{
    // NB if you repeat any test with the same memory it gets
    // much faster. Due to processor cache?

    // can set the device here - as long as no context has been created yet
//    setDevice(1);
//    Assert(getDevice() == 1);

    const long BUFSIZE = /*230*/ 100 *(1<<20); // 250m floats = 1GB

    Timer t1;

    for (int i=0; i<1; ++i)
    {
        // normal to normal
        {
            float *normal1 = new float[BUFSIZE];
            float *normal2 = new float[BUFSIZE];
            t1.start();
            memcpy(normal1, normal2, BUFSIZE * sizeof(float));
            t1.stop();

            std::cout << "memcpy " << BUFSIZE * sizeof(float) << " bytes, normal -> normal took "
                      << t1 << " seconds" << std::endl;

            delete normal1;
            delete normal2;
        }

        // normal to pinned
        {
            float *normal1 = new float[BUFSIZE];
            float *pinned1;
            cudaHostAlloc((void**)&pinned1, BUFSIZE * sizeof(float), cudaHostAllocDefault);
            t1.start();
            memcpy(pinned1, normal1, BUFSIZE * sizeof(float));
            t1.stop();

            std::cout << "memcpy " << BUFSIZE * sizeof(float) << " bytes, normal -> pinned took "
                      << t1 << " seconds" << std::endl;

            cudaFreeHost(pinned1);
            delete normal1;
        }

        // normal to device
        {
            float *normal1 = new float[BUFSIZE];

            float *device1;
            Assert2(cudaMalloc((void **) &device1, BUFSIZE * sizeof(float)) != cudaErrorMemoryAllocation,
                    "speed_test: cudaMalloc failed");
            Assert2(cudaThreadSynchronize() == cudaSuccess, "speed_test: cudaThreadSynchronize failed");
            Assert(device1);

            t1.start();
            cudaMemcpy(device1, normal1, BUFSIZE * sizeof(float), cudaMemcpyHostToDevice);
            t1.stop();

            std::cout << "cudaMemcpy " << BUFSIZE * sizeof(float) << " bytes, normal -> device took "
                      << t1 << " seconds" << std::endl;

            delete normal1;
            cudaFree(device1);
        }

        // pinned to device
        {
            float *pinned1;
            cudaHostAlloc((void**)&pinned1, BUFSIZE * sizeof(float), cudaHostAllocDefault);

            float *device1;
            Assert2(cudaMalloc((void **) &device1, BUFSIZE * sizeof(float)) != cudaErrorMemoryAllocation,
                    "speed_test: cudaMalloc failed");
            Assert2(cudaThreadSynchronize() == cudaSuccess, "speed_test: cudaThreadSynchronize failed");
            Assert(device1);

            cudaStream_t stream0;
            cudaError_t cts;
            cts = cudaStreamCreate(&stream0);
            Assert2(cts == cudaSuccess, "cudaStreamCreate failed");

            t1.start();
            cudaMemcpyAsync(device1, pinned1, BUFSIZE * sizeof(float), cudaMemcpyHostToDevice, stream0);
            Assert2(cudaThreadSynchronize() == cudaSuccess, "cudaMemcpyAsync failed");
            t1.stop();

            std::cout << "cudaMemcpyAsync " << BUFSIZE * sizeof(float) << " bytes, pinned -> device took "
                      << t1 << " seconds" << std::endl;

            cudaFreeHost(pinned1);
            cudaFree(device1);
        }
    }

    // can't do setDevice here - even if the same as previously. We already have a CUDA context.
//    setDevice(1);
//    Assert(getDevice() == 1);

}

#endif

typedef void *(*FUNC)(void *);

struct ThreadTestData
{
    int thread_index;
    bool result;
};

const int NFLOATS  = 10 * (1<<20);
//const int NDEVICES = 2;

// Allocate memory on a device
bool
cudaMallocTest(int n)
{
    float *data_d;
    CHECK(cudaMalloc((void **) &data_d, sizeof(float)*n) != cudaErrorMemoryAllocation);
    CHECK(cudaThreadSynchronize() == cudaSuccess);
    CHECK(data_d != 0);
    return data_d != 0;
}

void to_upper(std::string &s)
{
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))toupper);
}

#if 0

static void
threadFunc1(ThreadTestData *p)
{
    CHECK(!p->result);

    std::set<int> device_ids;
    int ndevices = gpuDevices().GPUsEnabled(&device_ids);

//    enableGPUs(2);

    int dev = *device_ids.find(p->thread_index % ndevices);
    setDevice(dev);

    for (int i=0; i<2; ++i)
    {
        Timer t;
        t.start();

        p->result = cudaMallocTest(NFLOATS);

        t.stop();
        cout << "pass=" << i << " thread=" << p->thread_index << " device=" << dev << ": cudaMalloc(" << NFLOATS << " floats) took " << t << " seconds" << endl;
    }

    cudaThreadExit();
}

TEST(cudaMalloc_thread_test)
{
    // What slows down cudaMallocs in threaded apps?
    // Answer: creating a CUDA context in each new thread.

    // run test in main thread
    for (int i=0; i<4; ++i)
    {
        Timer t;
        t.start();
        CHECK(cudaMallocTest(NFLOATS));
        t.stop();
        cout << "pass=" << i << " thread=(main) device=(default): cudaMalloc(" << NFLOATS << " floats) took " << t << " seconds" << endl;
    }

    // run in subthreads
    const int NTHREADS = 4;

    cout << "cudaMalloc_thread_test: NTHREADS = " << NTHREADS << endl;
    pthread_t threadID[NTHREADS];
    ThreadTestData data[NTHREADS];

    for(int i = 0; i < NTHREADS; i++)
    {
        data[i].thread_index = i;
        data[i].result = false;

        pthread_create(&threadID[i],
                       NULL,
                       (FUNC)threadFunc1,
                       (void*)(data+i));
    }

    for(int i = 0; i < NTHREADS; i++)
    {
        pthread_join(threadID[i], NULL);
    }

    bool ret = true;
    for(int i = 0; i < NTHREADS; i++)
    {
        ret = ret && data[i].result;
    }

    CHECK(ret);

    // shut down the context in the main thread - otherwise we can't call setDevice() again!
    cudaThreadExit();
}
#endif

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Cached.cpp
 *
 *  Created on: Mar 18, 2011
 *      Author: gareth
 */

#include "Cached.h"
#include "util.h"
#include "cuda_accel/cuda_accel.h"
#include "cuda_accel/GPUDevices.h"

#include <vector>

#include <UnitTest++/UnitTest++.h>

size_t             Cached::s_cached_mem_limit = 1000000000; // default 1GB
size_t             Cached::s_cached_mem = 0;
std::list<Cached*> Cached::s_cached_list;

pthread_mutex_t list_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;

#define PSIZE 10

Cached::Cached()
: m_size(0)
{
}

Cached::~Cached()
{
}

size_t
Cached::setLimit(size_t n) // bytes
{
    size_t old_limit = s_cached_mem_limit;
    s_cached_mem_limit = n;
    makeRoom(0);
    return old_limit;
}

size_t
Cached::getLimit() // bytes
{
    return s_cached_mem_limit;
}

size_t
Cached::getBytesUsed()
{
    return s_cached_mem;
}

void
Cached::clear()
{
    makeRoom(s_cached_mem_limit);
}

// NB: called by the derived class allocate()
void
Cached::allocate(size_t n) // bytes
{
    // The actual allocation is done by the derived class
    // This function must record the number of bytes allocated
    // and add an entry to the end of the list

// Mutex protection provided by the calling functions - but take care if modifying
//    pthread_mutex_lock(&list_mutex);
    m_size = n;
    s_cached_mem += m_size;
    s_cached_list.push_back(this);
    m_my_entry = s_cached_list.end();
    --m_my_entry;
//    pthread_mutex_unlock(&list_mutex);
}

// NB: called by the derived class deallocate()
void
Cached::deallocate()
{
    // Decrement the byte count and remove this entry from list
    if (m_size != 0)
    {
// Mutex protection provided by the calling functions - but take care if modifying
//        pthread_mutex_lock(&list_mutex);
        s_cached_list.erase(m_my_entry);
        s_cached_mem -= m_size;
        m_my_entry = s_cached_list.end();
        m_size = 0;
//        pthread_mutex_unlock(&list_mutex);
    }
}

// Make room for n bytes by deallocating from the start of the list
void
Cached::makeRoom(size_t n) // bytes
{
    // Other threads must not access the cached list
    pthread_mutex_lock(&list_mutex);

    list<Cached*>::iterator it = s_cached_list.begin();
    while ( s_cached_mem_limit < n + s_cached_mem && it !=  s_cached_list.end() )
    {
        (*it++)->deallocate(); // NB calls the derived class deallocate
    }

    pthread_mutex_unlock(&list_mutex);
}

// P1 caches data without using CachedArray
//
class P1 // caches the reciprocals of 1..10
{
public:
	P1() : count(0), m_f(0) {}
	~P1() { delete m_f; }

    float f(int i)
    {
    	if (m_f == 0)
    	{
    		make_f();
    	}

    	if (0<i && i<=PSIZE)
    	{
    	    return m_f[i-1];
    	}
    	else
    	{
    		return 0;
    	}
    }

    int count;

private:
    float *m_f;

    void make_f()
    {
    	m_f = new float[PSIZE];
    	for (int i=1; i<=PSIZE; ++i)
    	{
    		m_f[i-1] = 1/i;
    	}
    	++count;
    }
};

TEST(cached_P1)
{
	P1 p;

    CHECK_EQUAL(0, p.count);
	for (int i=1; i<=PSIZE; ++i)
	{
	    CHECK(p.f(i) == 1/i);
	}
	CHECK_EQUAL(1, p.count);
}


// P2 caches data using CachedArray
//
class P2 // caches the reciprocals of 1..10 using Cached
{
public:
    P2() : m_f() { /*cout << "default constructing" << endl;*/ }

    P2(P2 const & other) : m_f() // copy ctor does NOT copy the cached data
    { /*cout << "default constructing" << endl;*/ }

    ~P2() {}

    float f(int i)
    {
        CachedArray<float>::TempRef ref(m_f); // valid until ref goes out of scope

        if ( ! ref.hasValues())
        {
            make_f();
        }

        if (0<i && i<=PSIZE)
        {
            return ref[i-1];
        }
        else
        {
            return 0;
        }
    }

    CachedArray<float>::ConstTempRef getRef()
    {
        CachedArray<float>::TempRef ref(m_f);

        if ( ! ref.hasValues())
        {
            make_f();
        }

        Assert2(ref.hasValues(), "invalid reference");
        Assert2(ref.numValues() > 0, "invalid reference");
        Assert2(ref.values() > 0, "invalid reference");

        CachedArray<float>::ConstTempRef tref(ref);

        Assert2(tref.hasValues(), "invalid reference");
        Assert2(tref.numValues() > 0, "invalid reference");
        Assert2(tref.values() > 0, "invalid reference");

        return tref;

        // TempRef is converted to ConstTempRef
//        return CachedArray<float>::ConstTempRef(ref);
    }

    static int count;

private:
    CachedArray<float> m_f;

    void make_f()
    {
        CachedArray<float>::TempRef ref(m_f);

        ref.allocate(PSIZE);

        for (int i=1; i<=PSIZE; ++i)
        {
            ref[i-1] = 1/(float)i;
        }
        ++count;
    }
};

int P2::count = 0;

#if 1
// test use of copy ctors and assignment
TEST(cached_copy)
{
    CachedArray<float> m_f;
    CachedArray<float>::TempRef ref(m_f);
    ref.allocate(1);
    CHECK_EQUAL(1, ref.numValues());

    CachedArray<float>::ConstTempRef *refp1 = new CachedArray<float>::ConstTempRef(ref);
    CHECK_EQUAL(1, refp1->numValues());

    CachedArray<float>::TempRef *refp2 = new CachedArray<float>::TempRef(ref);
    CHECK_EQUAL(1, refp2->numValues());

//    *refp2 = *refp1; // TempRef = ConstTempRef : not allowed - won't compile

    *refp1 = *refp2;   // ConstTempRef = TempRef
    CHECK_EQUAL(1, refp1->numValues());

    *refp2 = ref;      // TempRef = TempRef
    CHECK_EQUAL(1, refp2->numValues());

    delete refp1;
    delete refp2;

}

TEST(cached_P2)
{
    P2 p;

    CHECK_EQUAL(0, p.count);
    for (int i=1; i<=PSIZE; ++i)
    {
        CHECK(p.f(i) == 1/(float)i);
    }
    CHECK_EQUAL(1, p.count);


    CachedArray<float>::ConstTempRef ref = p.getRef();
    CHECK_EQUAL(PSIZE*sizeof(float), Cached::getBytesUsed());

    CHECK_EQUAL(1, p.count);
    for (int i=1; i<=PSIZE; ++i)
    {
        CHECK(ref[i-1] == 1/(float)i);
    }
    CHECK_EQUAL(1, p.count);

}

TEST(cached_P2_copy)
{
    P2 p;
    CachedArray<float>::ConstTempRef *refp1 = 0, *refp2 = 0;

    {
        CachedArray<float>::ConstTempRef ref = p.getRef();
        refp1 = new CachedArray<float>::ConstTempRef(ref); // copy ctor
        refp2 = new CachedArray<float>::ConstTempRef(ref); // copy ctor
        *refp2 = *refp1; // operator=
    }

    // ref now out of scope, refp1, refp2 point to copies of it

    for (int i=1; i<=PSIZE; ++i)
    {
        CHECK((*refp1)[i-1] == 1/(float)i);
        CHECK((*refp2)[i-1] == 1/(float)i);
    }

    delete refp1;
    delete refp2;
}

TEST(cached_P2_limit)
{
    // Cached limit is not exceeded
    size_t old_limit = Cached::setLimit(1000);
    P2::count = 0;

    CHECK_EQUAL(1000, Cached::getLimit());
    CHECK_EQUAL(0, Cached::getBytesUsed());

    // Allocate 30 P2 objects
    const int NOBJ = 30;
    P2 p2array[NOBJ];

    // Since f() has not been called, nothing is cached yet
    CHECK_EQUAL(0, Cached::getBytesUsed());

    // Make each P2 cache its results. The first 25 will use up 1000 bytes
    // So the first 5 will be freed up to allow the last 5 to be cached

    for (int i=0; i<NOBJ; ++i)
    {
        CHECK_EQUAL(0.125, p2array[i].f(8));
    }

    // We are at the cached limit
    CHECK_EQUAL(1000, Cached::getBytesUsed());

    // Each P2 has calculated its data once
    CHECK_EQUAL(30, P2::count);

    // Do the last 10 again
    for (int i=20; i<NOBJ; ++i)
    {
        CHECK_EQUAL(0.125, p2array[i].f(8));
    }

    // No extra cost - they were all cached
    CHECK_EQUAL(30, P2::count);

    // Do the first 10 again
    for (int i=0; i<10; ++i)
    {
        CHECK_EQUAL(0.125, p2array[i].f(8));
    }

    // We had to calculate all those again
    // As we cached 0 we deleted 5
    // As we cached 1 we deleted 6, etc
    CHECK_EQUAL(40, P2::count);

    // return to default value
    Cached::setLimit(old_limit);
}

#endif

typedef void *(*FUNC)(void *);

struct ThreadTestData
{
    int thread_index;
    std::vector<P2> *p2s;
    bool result;
};

static void
threadFunc(ThreadTestData *p)
{

    int seed = p->thread_index * 42 + 2;
    srand(seed);

    // Set the device (just to see how it works):
    // 1) The default device in a new thread is 0, this is not affected by any set in the main thread
    // 2) You can set the device as many times as you like BEFORE a CUDA context is created in the thread
    // 3) Once there is a CUDA context context it is an error to call cudeSetDevice (even with the current value)
    // 4) It is OK for different threads to have the same device. (Kernels will be queued)
#if 0
    Assert(getDevice() == 0);
    int n_devices = gpuDevices().GPUsPresent();
    int device_id = p->thread_index % n_devices;
    setDevice(1);
    Assert(getDevice() == 1);
//    setDevice(device_id);
//    Assert(getDevice() == device_id);

    // create a CUDA context in this thread (we don't use it)
    // This seems to take a long time!
    float *pinned;
    cudaHostAlloc((void**)&pinned, 100 * sizeof(float), cudaHostAllocDefault);

//    setDevice(0); // can't set the device here

    // after calling cudaThreadExit() the context is gone and we can set the device again
    cudaThreadExit();
    setDevice(0);

    Assert(getDevice() == 0);
#endif

    int n = p->p2s->size();

    p->result = true;

    for (int i=0; i<10000; ++i)
    {
        int j = rand(n);
        P2 & p2 = (*(p->p2s))[j];
        CachedArray<float>::ConstTempRef ref = p2.getRef();

        Assert2(ref.hasValues(), "invalid reference in threadFunc");
        Assert2(ref.numValues() > 0, "invalid reference in threadFunc");
        Assert2(ref.values() > 0, "invalid reference in threadFunc");

        if (ref[8-1] != 0.125)
//        if ( p2.f(8) != 0.125 )
        {
            p->result = false;
        }

//        cout << p->thread_index << " : " << i << " : " << endl;

        // add some randomness
        struct timespec rqtp = { 0, 0 };
        struct timespec rmtp;
        if (rand(100)==0)
        {
            rqtp.tv_nsec = rand(50);
//            cout << p->thread_index << " : " << i << " : " << rqtp.tv_nsec << endl;
            nanosleep(&rqtp, &rmtp);
        }
    }
}

static bool
threadTest(std::vector<P2> *p2s)
{
    const int NTHREADS = 3;

    pthread_t threadID[NTHREADS];
    ThreadTestData data[NTHREADS];

    for(int i = 0; i < NTHREADS; i++)
    {
        data[i].thread_index = i;
        data[i].p2s = p2s;
        data[i].result = false;

        pthread_create(&threadID[i],
                       NULL,
                       (FUNC)threadFunc,
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

    return ret;
}

#if 1
TEST(cached_P2_threaded)
{
    // thread safe

    // Cached limit is not exceeded
    size_t old_limit = Cached::setLimit(1000);
    P2::count = 0;

    CHECK_EQUAL(1000, Cached::getLimit());
    CHECK_EQUAL(0, Cached::getBytesUsed());

    // allocate 25 P2 objects (all can be cached at once)
    std::vector<P2> p2a(25);

    // launch two threads that randomly call them
    CHECK_EQUAL(0, P2::count); // each value calculated only once
    CHECK(threadTest(&p2a));

    CHECK_EQUAL(25, P2::count); // each value calculated only once

    // allocate 50 P2 objects (will force lots of recalculation)
    std::vector<P2> p2b(50);

    // launch two threads that randomly call them
    CHECK(threadTest(&p2b));

//    CHECK_EQUAL(25, P2::count); // each value calculated lots of times

    // return to default value
    Cached::setLimit(old_limit);

}
#endif

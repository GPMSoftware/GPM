/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Cached.h
 *
 *  Created on: Mar 18, 2011
 *      Author: gareth
 *
 * A class to hold cached data.
 *
 * There is a global limit on the amount of cached data,
 * and when it is exceeded cache data is reclaimed.
 *
 * The owner of the cached data must test to see if a cached value
 * exists, and if not be prepared to recalculate the data.
 *
 * Thread safe.
 *
 * Currently there is a single cache of data. Could be extended to multiple caches.
 */

#ifndef CACHED_H_
#define CACHED_H_

#include "Assert.h"
#include "list.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <pthread.h>

// This mutex controls access to the list of cached objects
extern pthread_mutex_t list_mutex;

// Base class for CachedArray
class Cached
{
public:
    Cached();
    virtual ~Cached();

    static size_t setLimit(size_t n); // bytes

    static size_t getLimit(); // bytes

    static size_t getBytesUsed();

    static void clear();

protected:

    // NB: called by the derived class allocate()
    virtual void allocate(size_t n); // bytes

    // NB: called by the derived class deallocate()
    virtual void deallocate();

    // Make room for n bytes by deallocating from the start of the list
    static void makeRoom(size_t n); // bytes

private:
    size_t m_size; // bytes allocated by this object
    list<Cached*>::iterator m_my_entry;

    static size_t      s_cached_mem_limit; // max cached memory (bytes)
    static size_t      s_cached_mem;       // current cached memory (bytes)
    static std::list<Cached*> s_cached_list;      // list of cached objects
};


template <typename T>
class CachedArray : public Cached
{
public:

    // A thread safe reference to the cached array.
    // The data is mutex-locked until this object is destroyed.
    // NB: Keeping a permanent copy of it will result in deadlock!
    //
    // ConstTempRef: does not allow the data to be modified, and can be provided to clients
    //             : holds only the m_mutex for the CachedArray
    //
    // TempRef: adds the non-const functions needed by the owner of the CachedArray
    //        : holds m_mutex and list_mutex to allow (de)allocation
    //

    class TempRefBase
    {
    public:
        TempRefBase(CachedArray<T> *a) : m_a(a) { }
        TempRefBase(TempRefBase const &other) : m_a(other.m_a)  { }
        virtual ~TempRefBase()  { }

        bool hasValues() const            { return m_a->hasValues(); }
        T const * values() const          { return m_a->values(); }
        T const & operator[](int i) const { return (*m_a)[i]; }
        size_t numValues() const   { return m_a->numValues(); }

    protected:
        CachedArray<T> *m_a;

        // Surprisingly, we can't actually see the above though a TempRefBase reference without this
        friend class CachedArray<T>::ConstTempRef;
    };

    class ConstTempRef : public TempRefBase
    {
    public:
        ConstTempRef(CachedArray<T> &a) : TempRefBase(&a)
        {
            lockMutextes();
        }

        virtual ~ConstTempRef()
        {
            unlockMutextes();
        }

        // Can construct a ConstTempRef from either a ConstTempRef or a TempRef
        ConstTempRef(TempRefBase const &other) : TempRefBase(other.m_a)
        {
            // NB this is a recursive lock. The copy ctor
            // puts another lock on, so that the lock is in
            // force until all TempRefBase objects are destroyed.
            // (In particular this makes returning a copy of a local TempRefBase safe,
            // even without the return value optimization).

            // Note that if (other) is a TempRef then when
            // (other) goes out of scope it releases its list_mutex lock and
            // we are left with just the m_mutex lock

            lockMutextes();
        }

        // Can assign a ConstTempRef from either a ConstTempRef or a TempRef
        ConstTempRef &
        operator=(TempRefBase const &other)
        {
            // See comment for copy ctor
            unlockMutextes(); // decrement locks on currently referenced object
            ConstTempRef::m_a = other.m_a;
            lockMutextes();   // increment locks on newly referenced object

            return *this;
        }

    private:
        void lockMutextes()
        {
            pthread_mutex_lock(&ConstTempRef::m_a->m_mutex);
        }

        void unlockMutextes()
        {
            pthread_mutex_unlock(&ConstTempRef::m_a->m_mutex);
        }

    };

    class TempRef : public TempRefBase
    {
    public:

        TempRef(CachedArray<T> &a) : TempRefBase(&a)
        {
            lockMutextes();
        }

        virtual ~TempRef()
        {
            unlockMutextes();
        }

        // non-const functions
        void allocate(size_t n)              { TempRefBase::m_a->allocate(n); }
        void deallocate()                 { TempRefBase::m_a->deallocate(); }
        T* values()                       { return TempRefBase::m_a->values(); }
        T& operator[](int i)              { return (*TempRefBase::m_a)[i]; }

        // Can construct a TempRef from a TempRef but NOT from a ConstTempRef
        // (because if we are going to have list_mutex we must grab it before the m_mutex)
        TempRef(TempRef const &other) : TempRefBase(other.m_a)
        {
            // See comment for ConstTempRef
            lockMutextes();
        }

        TempRef &
        operator=(TempRef const &other)
        {
            // See comment for copy ctor
            unlockMutextes(); // decrement locks on currently referenced object
            TempRefBase::m_a = other.m_a;
            lockMutextes();   // increment locks on newly referenced object

            return *this;
        }

    private:
        // NB These are NOT allowed. Since they won't be generated it is best not event to define them.
//        TempRef(ConstTempRef const &other);
//        TempRef &operator=(ConstTempRef const &other);

        void lockMutextes()
        {
            // Acquire list_mutex first
            pthread_mutex_lock(&list_mutex);
            pthread_mutex_lock(&TempRefBase::m_a->m_mutex);
        }

        void unlockMutextes()
        {
            pthread_mutex_unlock(&list_mutex);
            pthread_mutex_unlock(&TempRefBase::m_a->m_mutex);
        }
    };

	CachedArray()
	: m_numvalues(0)
	, m_values(0)
	{
	     pthread_mutexattr_t attr;
	     pthread_mutexattr_init(&attr);
	     pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
	     pthread_mutex_init(&m_mutex, &attr);
	}

    virtual ~CachedArray()
	{
	    deallocate();
	}

private:
    // No copying. It is up to the client class (the class that
    // contains a CachedArray) to manage it. We do not want compiler-
    // generated functions for client classes to copy the CachedArray.
    // A client class should explicitly define its own copy and assignment
    // functions and decide what to do with its cached data (Options
    // would include creating an empty CachedArray, or some kind of
    // reference counting. We can't decide that.
    // Also, we don't want CachedArrays passed to functions by value.
    // They should be using ConstTempRef or TempRef.
    CachedArray(CachedArray const & other);
    CachedArray& operator=(CachedArray const &other); // { Assert(false); return *this; }

	void allocate(size_t n) // objects of type T
	{
	    // Other threads must not access this object OR the cached list
        // Always grab the list_mutex first to avoid deadlock
        pthread_mutex_lock(&list_mutex);
        pthread_mutex_lock(&m_mutex);

	    Cached::makeRoom(n * sizeof(T));

	    m_numvalues = n;

// Putting the cache in pinned memory does not work
// - we get lots of errors in the unit tests to do with threading.
// It appears that cudaHostAlloc is not thread safe
//#define PINNED
#ifndef PINNED
	    m_values = new T[n]; // allocate with new
#else
	    m_values = 0;
	    cudaError_t err;
        err = cudaHostAlloc((void**)&m_values, n * sizeof(T), cudaHostAllocPortable); // allocate pinned memory
        Assert2(err == cudaSuccess, "cudaHostAlloc failed");
        Assert(m_values != 0);
#endif
	    Cached::allocate(m_numvalues * sizeof(T));

        pthread_mutex_unlock(&m_mutex);
        pthread_mutex_unlock(&list_mutex);
	}

	void deallocate()
	{
        // Other threads must not access this object OR the cached list
        // Always grab the list_mutex first to avoid deadlock
        pthread_mutex_lock(&list_mutex);
        pthread_mutex_lock(&m_mutex);

#ifndef PINNED
        delete [] m_values;
#else
        if (m_values)
        {
            cudaError_t err;
            err = cudaFreeHost((void*)m_values);
            Assert2(err == cudaSuccess, "cudaFreeHost failed");
        }
#endif

	    m_values = 0;
	    m_numvalues = 0;
	    Cached::deallocate();

        pthread_mutex_unlock(&m_mutex);
        pthread_mutex_unlock(&list_mutex);
	}

	bool hasValues() const
	{
	    return m_values !=0;
	}

    size_t numValues() const
    {
        return m_numvalues;
    }

	T * values()
	{
	    return m_values;
	}

    T const * values() const
    {
        return m_values;
    }

    void checkIndex(long i) const
    {
        if (i < 0 || (unsigned)i >= m_numvalues)
        {
            Assert2(false, "CachedArray index out of range");
        }
    }

	T& operator[](long i)
	{
	    checkIndex(i);
	    return m_values[i];
	}

    T const & operator[](long i) const
    {
        checkIndex(i);
        return m_values[i];
    }

	size_t m_numvalues;
    T * m_values;
    pthread_mutex_t m_mutex; // This recursive mutex controls access to the individual object

};

#endif /* CACHED_H_ */

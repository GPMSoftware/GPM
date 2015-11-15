/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * util.h
 *
 *  Created on: Dec 7, 2009
 *      Author: gareth
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sys/times.h>
#include <sys/resource.h>

#include "fand/Assert.h"

// random double in range 0 to 1
double drand();

// random integer from 0 to range - 1
int rand(int range);

// random integer from start to end inclusive
int rand(int start, int end);

// NB the original implementation using clock() is fine for single-threaded apps.
// However in multi-threded apps it seems to include clock cycles for ALL threads (not just the calling thread)
// and so overestimates the time.
//
// The new implementation uses times() which returns an actual clock time.
//
class Timer
{
public:
	Timer() : clocks_per_sec(sysconf(_SC_CLK_TCK)) { start(); }
	Timer& start()  { m_start = times(0); return *this; }
	Timer& stop()   { m_stop  = times(0); return *this; }
	double read() const { return double(m_stop - m_start)/clocks_per_sec; }

private:
	long clocks_per_sec;
	clock_t m_start, m_stop;
};

std::ostream& operator<<(std::ostream &os, const Timer& timer);

class MemMeter
{
public:
	MemMeter() : k_per_page(getpagesize()/1024) {}

	// return the memory currently in use by this process (k)
	long read() const
	{
		// NB getrusage does not return memory data in linux!
		// use the proc filesystem
		std::ifstream f("/proc/self/statm");
		long total_mem_k;
		f >> total_mem_k;
		total_mem_k *= k_per_page;    // the memory in statm is in PAGES not k
		return total_mem_k;
	}

    // return the total physical memory on the system (k)
    long total()
    {
      return sysconf(_SC_PHYS_PAGES) * k_per_page;
    }

    int k_per_page;
};

// A class for counting other classes. Derive from it: class A : public Counted<A>
template<class T>
struct Counted
{
public:
	static int m_created, m_destroyed;

	Counted() { m_created++; }

	Counted(const Counted &) { m_created++; }

	virtual ~Counted() { m_destroyed++; Assert(m_created >= m_destroyed); }

	static int count() { return m_created - m_destroyed; }
};

template<class T> int Counted<T>::m_created = 0;
template<class T> int Counted<T>::m_destroyed = 0;

// a function for seeing if we are in the main thread
extern const pthread_t main_thread_id;
bool isMainThread();

std::ostream& operator<<(std::ostream &os, const MemMeter& m);

// parse a list of tokens into a vector of strings
// repeated separators are treated as a single separator  ("A,,B" --> "A", "B")
std::vector<std::string>
split(std::string const &s, std::string separator = " ");

// parse a list of tokens into a vector of strings
// repeated separators denote empty strings ("A,,B" --> "A", "", "B")
std::vector<std::string>
split2(std::string const &s, std::string separator = ",");

// remove leading and trailing spaces from vector of strings
void cleanup(std::string &word);
void cleanup(std::vector<std::string> &words);

// output a vector in human-readable form
template <class T>
std::ostream &
operator<<(std::ostream &os, std::vector<T> const &x)
{
	typename std::vector<T>::const_iterator it = x.begin();
	os << "(";
	if (it != x.end()) os << *it++;
	while(it != x.end())
	{
		os << ", " << *it++;
	}
	os << ")";
	return os;
}

template<class T>
void order(T &a, T&b)
{
    if (b < a)
    {
        std::swap(a, b);
    }
}

template<class T>
int delta(T const &a, T const &b)
{
    return (a == b);
}

void exit_debug (int status);

// get float environment variable
float getFloatEnv(const char *env_name, float default_val = 0);

// get int environment variable
int getIntEnv(const char *env_name, int default_val = 0);

// get bool environment variable
// returns the default value unless explicitly set to the other value
// (ensures that set to "" behaves the same as unset)
bool getBoolEnv(const char *env_name, bool default_val);

// get string environment variable ("" if not set)
const char * getStringEnv(const char *env_name, const char *default_val = "");

std::string
baseName(std::string const &filename);

std::string
path(std::string const &filename);

std::string
dateStr();

template<class AA, class BB>
void
mapInv(std::map<AA, BB> const &map_in, std::multimap <BB, AA> &mmap_out)
{
    mmap_out.clear();

    typename std::map<AA, BB>::const_iterator it;

    for (it = map_in.begin(); it != map_in.end(); ++it)
    {
        mmap_out.insert(make_pair(it->second, it->first));
    }
}

template<class T>
void
vec2set(std::vector<T> const &vec_in, std::set <T> &set_out)
{
    set_out.clear();

    size_t size = vec_in.size();

    for (size_t i = 0; i < size; ++i)
    {
        set_out.insert(vec_in[i]);
    }
}

template<class T>
void
vec2set(typename std::vector<T>::const_iterator it_begin,
        typename std::vector<T>::const_iterator it_end,
        std::set <T> &set_out)
{
    set_out.clear();

    typename std::vector<T>::const_iterator it;

    for (it = it_begin; it != it_end; ++it)
    {
        set_out.insert(*it);
    }
}

// users home directory ($HOME)
const char *homeDir();
std::string homeDirStr();

void to_upper(std::string &s);

const double pi_180 = 3.1415926535 / 180;
inline double radians(double deg) { return deg * pi_180; }
inline double degrees(double rad) { return rad / pi_180; }

// a random number generator.
// Each instance generates a repeatable random sequence,
// independent of other instances. Uses POSIX rand_r.
class Rand
{
public:
    Rand(unsigned int seed = 1) : m_state(seed) {}

    // random double in range 0 to 1
    double drand();

    // random integer from 0 to range - 1
    int rand(int range);

    // random integer from start to end inclusive
    int rand(int start, int end);

private:
    unsigned int m_state;
};

// hash functions (simple)
#if 0
//jenkins_one_at_a_time_hash
uint32_t hash(unsigned char *key, size_t key_len)
{
    uint32_t hash = 0;
    size_t i;

    for (i = 0; i < key_len; i++) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

template <class InputIterator, class T>
uint32_t hash ( InputIterator first, InputIterator last)
{
    uint32_t hash = 0;

    while ( first!=last )
    {
        hash += hash(*first++);
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}

template <class T>
uint32_t hash(const T&);

// specializations for built-in types
//template <> uint32_t hash(const int& i) { return i; }
// ...

template <class T>
uint32_t hash(const std::vector<T> &v)
{
	return hash(v.begin(), v.end());
}

template <class T, class U>
uint32_t hash(const T &t, const U &u)
{
	return 0; // what ???
}

template <class T, class U>
uint32_t hash(const std::pair<T, U> &p)
{
	return hash(p.first, p.second);
}

template <class T, class U>
uint32_t hash(const std::map<T, U> &m)
{
	return hash(m.begin(), m.end());
}
#endif


#endif /* UTIL_H_ */

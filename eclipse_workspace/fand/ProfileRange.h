/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ProfileRange.h
 *
 *  Created on: Jul 1, 2011
 *      Author: gareth
 */

#ifndef PROFILERANGE_H_
#define PROFILERANGE_H_

#include <stddef.h>
#include <vector>
#include <boost/smart_ptr.hpp>

#include "Profile.h"

class Database;

class ProfileRangeBase // base class
{
public:
	ProfileRangeBase(size_t begin=0, size_t end=0)
	: m_begin(begin), m_end(end<begin ? begin : end), m_data_ok(false)
	{}

	ProfileRangeBase(const ProfileRangeBase &other)
	: m_begin(other.m_begin), m_end(other.m_end), m_data_ok(other.m_data_ok)
    {
        Assert(m_data_ok == false); // no copying once the data is readable
    }

	virtual ~ProfileRangeBase() {}

	virtual ProfileRangeBase *clone() const = 0;

	virtual bool readMyData() = 0;

	// in the case of a database read, operator[] is an "efficient" read and can be used
	// only after readMyData() has explicitly been called to read all the data in the range.
	//
	// get() is a "slow" read - If the data has not been read the database will be queried for the
	// individual profile.
	//
	// find() reads from the database based on id
	//
	// (For other derived classes the behavior is the same).
	//
    virtual Profile & operator[](int i) = 0;

    virtual Profile   get(int i) const = 0;

    virtual Profile   find(std::string id) const = 0;

	size_t begin() const { return m_begin; }

	size_t end()   const { return m_end; }

    size_t size()  const { return m_end - m_begin; }

    bool gotData() const { return m_data_ok; }

protected:
	size_t m_begin, m_end;
	bool m_data_ok;

private:
    ProfileRangeBase &operator=(ProfileRangeBase const &);

	friend class ProfileRange;
};

// ProfileRangeVec holds a reference to an external vector
class ProfileRangeVec : public ProfileRangeBase
{
public:
    ProfileRangeVec(std::vector<Profile> &vec)
    : ProfileRangeBase(0, vec.size()), m_vec(vec) {}

    ProfileRangeVec(const ProfileRangeVec& other)
    : ProfileRangeBase((const ProfileRangeBase &)other), m_vec(other.m_vec)
    {
    }

	virtual ~ProfileRangeVec() {}

	virtual ProfileRangeVec *clone() const
	{
		return new ProfileRangeVec(*this);
	}

    virtual bool readMyData();

    virtual Profile & operator[](int i);

    virtual Profile   get(int i) const;

    virtual Profile   find(std::string id) const;

private:
    ProfileRangeVec &operator=(ProfileRangeVec const &);

	std::vector<Profile> &m_vec; // reference to an external vector
};


// ProfileRangeRefCtdVec holds a reference counted vector
class ProfileRangeRefCtdVec : public ProfileRangeBase
{
public:
    ProfileRangeRefCtdVec(std::vector<Profile> *pv) // must be allocated in the heap!
    : ProfileRangeBase(0, pv->size()), m_vec(pv)
    {
    }

    ProfileRangeRefCtdVec(const ProfileRangeRefCtdVec& other)
    : ProfileRangeBase((const ProfileRangeBase &)other), m_vec(other.m_vec)
    {
    }

    virtual ~ProfileRangeRefCtdVec()
    {
    }

    virtual ProfileRangeRefCtdVec *clone() const
    {
        return new ProfileRangeRefCtdVec(*this);
    }

    virtual bool readMyData();

    virtual Profile & operator[](int i);

    virtual Profile   get(int i) const;

    virtual Profile   find(std::string id) const;


private:
    ProfileRangeRefCtdVec &operator=(ProfileRangeRefCtdVec const &);

    boost::shared_ptr< std::vector<Profile> > m_vec;

    friend class TestProfileRangeRefCtdVec;
};

class ProfileRangeDB : public ProfileRangeBase
{
public:
    ProfileRangeDB(boost::shared_ptr<Database> &db, std::string const &query);

    ProfileRangeDB(const ProfileRangeDB& other)
    : ProfileRangeBase((const ProfileRangeBase &)other)
    , m_db(other.m_db)
    , m_query(other.m_query)
    , m_vec(other.m_vec)
    {}

    virtual ~ProfileRangeDB() {}

	virtual ProfileRangeDB *clone() const
	{
        return new ProfileRangeDB(*this);
	}

    virtual bool readMyData();

    virtual Profile & operator[](int i);          // efficient read (PRE: readMyData() has been called)

    virtual Profile   get(int i) const;           // possibly inefficient read (will read an individual element from source if necessary)

    virtual Profile   find(std::string id) const; // search external source (inefficient)

private:
    ProfileRangeDB &operator=(ProfileRangeDB const &);

    bool readData(std::vector<Profile> &vec, size_t begin, size_t end);

    mutable boost::shared_ptr<Database> m_db;   // Reference-counted pointer to MySQL database
	std::string m_query;                // SQL query
	std::vector<Profile> m_vec;         // vector read from database
};

class ProfileRange // wrapper for a ProfileRangeBase object
{
public:
    ProfileRange() : m_obj(0) {}

	ProfileRange(ProfileRangeBase const &obj)
	: m_obj(obj.clone())
	{
	}

    ProfileRange(ProfileRange const &pr)
    : m_obj(pr.m_obj->clone())
    {
    }

    ProfileRange(ProfileRange const &pr, size_t begin, size_t end)
    : m_obj(pr.m_obj->clone())
    {
        m_obj->m_begin = begin;
        m_obj->m_end   = end < begin ? begin : end;
    }

    // convenience c'tor
    ProfileRange(std::vector<Profile> &vec)
    : m_obj(new ProfileRangeVec(vec)) {}

	~ProfileRange() { delete m_obj; }

	void clear() { delete m_obj; m_obj = 0; }

	ProfileRange &operator=(ProfileRange const & other)
	{
	    delete m_obj;
	    m_obj = other.m_obj->clone();
	    return *this;
	}

    size_t size()  const { return m_obj == 0 ? 0 : m_obj->size(); }

    bool gotData() const { return m_obj == 0 ? 0 : m_obj->gotData(); }

    // The following functions require size() != 0 or gotData() == true

    bool readMyData()
    {
        return m_obj->readMyData();
    }

    Profile & operator[](int i) { return (*m_obj)[i]; } // requires readMyData() = true, or gotData() == true

    Profile   get(int i) const { return m_obj->get(i); }

    Profile   find(std::string id) const { return m_obj->find(id); }

    size_t begin() const { return m_obj->begin(); }

    size_t end()   const { return m_obj->end(); }

private:
	ProfileRangeBase *m_obj;
};


#endif /* PROFILERANGE_H_ */

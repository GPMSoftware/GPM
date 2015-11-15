/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * PMF.h
 *
 *  Created on: Dec 4, 2009
 *      Author: gareth
 */

#ifndef PMF_H_
#define PMF_H_

#include "util.h"
#include "Assert.h"

#include <map>
#include <vector>
#include <stdexcept>

// Probability Mass Function
//
// This class is a map from objects of type T to probabilities
// representing a Probability Mass Function, i.e.
//
// PMF<Obj> pmf;
// Obj A, B, C, D;
//
// pmf[A] = 0.1; // set the probability of the pmf yielding A, B, C, D
// pmf[B] = 0.2;
// pmf[C] = 0.3;
// pmf[D] = 0.4;
//
// pmf.normalize(); // just in case
//
// const Obj& random_obj = pmf.pick();
//
// ! Dangerous to derive from std::map (no virtual destructor!)
//


template<class T>
class PMF // : public Counted<PMF<T> >
{
public:
	typedef typename std::map<T, double>::const_iterator const_iterator;
	typedef typename std::map<T, double>::iterator       iterator;
	typedef typename std::map<T, double>::value_type     value_type;

	PMF() {}

	struct POD // used for 'static' construction
	{
		T m_key;
		double m_value;
	};

	PMF(const POD *data, size_t data_size)
	{
		for (size_t i=0; i<data_size; ++i)
		{
			m_map.insert(std::make_pair(data[i].m_key, data[i].m_value));
		}
	}

	PMF(std::vector<POD> const & data)
	{
		for (int i=0; i<data.size(); ++i)
		{
			m_map.insert(std::make_pair(data[i].m_key, data[i].m_value));
		}
	}

// Now OK if we need it
//	virtual ~PMF();

	inline iterator       begin()         { return m_map.begin();  }
	inline const_iterator begin()   const { return m_map.begin();  }
	inline iterator       end()           { return m_map.end();    }
	inline const_iterator end()     const { return m_map.end();    }
	inline iterator       rbegin()        { return m_map.rbegin(); }
	inline const_iterator rbegin()  const { return m_map.rbegin(); }
	inline size_t         size()    const { return m_map.size();   }
	inline bool           empty()   const { return m_map.empty();  }
	inline void           clear()         { return m_map.clear();  }
	inline double& operator[] (const T&x) { return m_map[x];       }
    inline double  operator() (const T&x) const { return m_map.find(x)->second; }
	inline iterator       find (const T&x) { return m_map.find(x); }
	inline const_iterator find (const T&x) const { return m_map.find(x); }
	inline std::pair<iterator, bool> insert (const value_type&x) { return m_map.insert(x); }
	inline void           erase(iterator pos) { m_map.erase(pos); }

	inline double val(const T&x) const // returning 0 if value not present
	{
	    const_iterator it = m_map.find(x);
	    if (it == m_map.end())
	    {
	        return 0;
	    }
	    else
	    {
	        return it->second;
	    }
	}

	double sum() const
	{
		double ret = 0;

		for (const_iterator it = m_map.begin(); it != m_map.end(); ++it)
		{
			ret += it->second;
		}
		return ret;
	}

	// ensure probabilities sum to val
	bool normalize(double val = 1.0)
	{
		double s = sum();

		if (s == 0 || s == val)
		{
			return (s == val);
		}

		double w = val/s;
		for (iterator it = m_map.begin(); it != m_map.end(); ++it)
		{
			it->second *= w;
		}

		return true;
	}

	// pick a T at random according to the probability distribution
	const T& pick() const
	{
		if (m_map.empty() )
			throw std::out_of_range("PMF::pick(): empty container");

		double r = drand();

		double sum = 0;
		for (const_iterator it = m_map.begin(); it != m_map.end(); ++it)
		{
			sum += it->second;
			if (r < sum)
			{
				return it->first;
			}
		}

		return m_map.rbegin()->first; // just in case
	}

	// return reference to value with largest probability
	const T& mode() const
	{
		if (m_map.empty() )
			throw std::out_of_range("PMF::mode(): empty container");

		const_iterator it = m_map.begin();
		const_iterator max_it = it;
		for (++it; it != m_map.end(); ++it)
		{
			if (it->second > max_it->second)
			{
				max_it = it;
			}
		}

		return max_it->first;
	}

	// pre: PMFs need NOT have corresponding entries to (this). Entries are added as needed.
	PMF<T> operator+=(PMF<T> const &pmf)
	{
		for (typename PMF<T>::const_iterator it = pmf.begin(); it != pmf.end(); ++it)
		{
			m_map[it->first] += it->second;
		}
		return *this;
	}

	// pre: PMFs need NOT have corresponding entries to (this). Entries are added as needed.
	PMF<T> operator-=(PMF<T> const &pmf)
	{
		for (typename PMF<T>::const_iterator it = pmf.begin(); it != pmf.end(); ++it)
		{
			m_map[it->first] -= it->second;
		}
		return *this;
	}

	bool operator==(PMF<T> const &rhs) const
	{
		return this->m_map == rhs.m_map;
	}

	bool operator!=(PMF<T> const &rhs) const
	{
		return ! (*this == rhs);
	}

	bool operator<(PMF<T> const &rhs) const
	{
		return this->m_map < rhs.m_map;
	}

	bool operator>(PMF<T> const &rhs) const
	{
		return rhs.m_map < this->m_map;
	}

private:
	// use default copy ctor and assign
	// PMF(PMF<T> const &other);
	// PMF<T>& operator=(PMF<T> const &other);

	std::map<T, double> m_map;
};

template<class T>
PMF<T>& operator*=(PMF<T> &pmf, double d)
{
	for (typename PMF<T>::iterator it = pmf.begin(); it != pmf.end(); ++it)
	{
		it->second *= d;
	}
	return pmf;
}

template<class T>
PMF<T>& operator+=(PMF<T> &pmf, double d)
{
    for (typename PMF<T>::iterator it = pmf.begin(); it != pmf.end(); ++it)
    {
        it->second += d;
    }

    return pmf;
}

template<class T>
PMF<T> operator*(PMF<T> const &pmf, double d)
{
	PMF<T> ret = pmf;
	ret *= d;
	return ret;
}

template<class T>
PMF<T> operator*(double d, PMF<T> const &pmf)
{
	return pmf * d;
}

// pre: PMFs have corresponding entries
//template<class T>
//PMF<T> mmin(PMF<T> const &pmf1, PMF<T> const &pmf2)
//{
//	PMF<T> ret;
//	for (typename PMF<T>::const_iterator it1 = pmf1.begin(); it1 != pmf1.end(); ++it1)
//	{
//		typename PMF<T>::const_iterator it2 = pmf2.find(it1->first);
//		Assert(it2 != pmf2.end());
//		ret.insert(make_pair(it1->first, min(it1->second, it2->second)));
//	}
//	return ret;
//}

// pre: PMFs have corresponding entries
template<class T>
PMF<T> operator+(PMF<T> const &pmf1, PMF<T> const &pmf2)
{
	PMF<T> ret = pmf1;
	return ret += pmf2;
}

// pre: PMFs have corresponding entries
template<class T>
PMF<T> operator-(PMF<T> const &pmf1, PMF<T> const &pmf2)
{
	PMF<T> ret = pmf1;
	return ret -= pmf2;
}

// sparse implementation: PMFs need NOT have corresponding entries
template<class T>
double dot(PMF<T> const &pmf1, PMF<T> const &pmf2)
{
	double ret = 0;
	for (typename PMF<T>::const_iterator it1 = pmf1.begin(); it1 != pmf1.end(); ++it1)
	{
		typename PMF<T>::const_iterator it2 = pmf2.find(it1->first);
		if(it2 != pmf2.end())
		{
			ret += it1->second * it2->second;
		}
	}
	return ret;
}

// sparse implementation: PMFs need NOT have corresponding entries
template<class T>
PMF<T> operator*(PMF<T> const &pmf1, PMF<T> const &pmf2)
{
	PMF<T> ret;
	for (typename PMF<T>::const_iterator it1 = pmf1.begin(); it1 != pmf1.end(); ++it1)
	{
		typename PMF<T>::const_iterator it2 = pmf2.find(it1->first);
		if(it2 != pmf2.end())
		{
			double p = it1->second * it2->second;
			if (p != 0)
			{
				ret[it1->first] = p;
			}
		}
	}
	return ret;
}

// sparse implementation: PMFs need NOT have corresponding entries
template<class T>
PMF<T> operator/(PMF<T> const &pmf1, PMF<T> const &pmf2)
{
	PMF<T> ret;
	for (typename PMF<T>::const_iterator it1 = pmf1.begin(); it1 != pmf1.end(); ++it1)
	{
		typename PMF<T>::const_iterator it2 = pmf2.find(it1->first);
		if(it2 != pmf2.end() && it2->second != 0)
		{
			double p = it1->second / it2->second;
			if (p != 0)
			{
				ret[it1->first] = p;
			}
		}
	}
	return ret;
}

template <class T, class U>
std::ostream &
operator<<(std::ostream &os, std::pair<T, U> const &x)
{
    os << "[" << x.first << ", " << x.second << "]";
    return os;
}

// output a vector in human-readable form
template <class T>
std::ostream &
operator<<(std::ostream &os, PMF<T> const &x)
{
    os << "(";

    typename PMF<T>::const_iterator it;
    for (it = x.begin(); it != x.end(); ++it)
    {
        if (it != x.begin()) os << ", ";
        os << "(" << it->first << " ==> " << it->second << ")";
    }
    os << ")";
    return os;
}

#endif /* PMF_H_ */

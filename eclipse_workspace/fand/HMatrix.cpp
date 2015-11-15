/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * HMatrix.cpp
 *
 *  Created on: Jan 7, 2010
 *      Author: gareth
 */

#include "HMatrix.h"

#include <UnitTest++/UnitTest++.h>

#include <math.h>

HMatrix::HMatrix()
{
}

// construct HMatrix from allele frequencies and a subpopulation model
//
HMatrix::HMatrix(const PMF<Allele> &pmf, const SubPopModel &spm)
{
	PMF<Allele>::const_iterator ipmf1, ipmf2;

	for (ipmf1 = pmf.begin(); ipmf1 != pmf.end(); ++ipmf1)
	{
		for (ipmf2 = ipmf1; ipmf2 != pmf.end(); ++ipmf2)
		{
			double val = spm.prob((ipmf1 == ipmf2), ipmf1->second, ipmf2->second);
			m_pmf[std::make_pair(ipmf1->first, ipmf2->first)] += val;
		}
	}

//	normalize();
}

HMatrix::HMatrix(const PMF<Allele> &pmf1, const PMF<Allele> &pmf2)
{
	PMF<Allele>::const_iterator ipmf1, ipmf2;

	for (ipmf1 = pmf1.begin(); ipmf1 != pmf1.end(); ++ipmf1)
	{
		for (ipmf2 = pmf2.begin(); ipmf2 != pmf2.end(); ++ipmf2)
		{
			// HW
			double val = ipmf1->second * ipmf2->second;

#if 0
			PMF< std::pair<Allele, Allele> >::const_iterator it;
			if ((it = find(ipmf1->first, ipmf2->first)) != m_pmf.end()) // already exists
			{
				val += it->second;
			}
			set(ipmf1->first, ipmf2->first, val);
#else
			// optimized - 30% faster

			if (ipmf1->first < ipmf2->first)
			{
				m_pmf[std::make_pair(ipmf1->first, ipmf2->first)] += val;
			}
			else
			{
				m_pmf[std::make_pair(ipmf2->first, ipmf1->first)] += val;
			}
#endif
		}
	}
}

PMF< std::pair<Allele, Allele> >::const_iterator
HMatrix::find(const Allele &a1, const Allele &a2) const
{
	PMF< std::pair<Allele, Allele> >::const_iterator it;
	if (a1 < a2)
	{
		it = m_pmf.find(std::make_pair(a1, a2));
	}
	else
	{
		it = m_pmf.find(std::make_pair(a2, a1));
	}
	return it;
}

double
HMatrix::operator()(const Allele &a1, const Allele &a2) const
{
	PMF< std::pair<Allele, Allele> >::const_iterator it = find(a1, a2);
//	Assert2(it != m_pmf.end(), "HMatrix::operator(): element does not exist");
	if (it == m_pmf.end()) return 0; // SPARSE

	return it->second;
}

double
HMatrix::margP(const Allele &a) const
{
    double p = 0;
    PMF< std::pair<Allele, Allele> >::const_iterator it;
    for(it = m_pmf.begin(); it != m_pmf.end(); ++it)
    {
//        cout << "(" << it->first.first << ", " <<
        p += ( int(it->first.first == a) + int(it->first.second == a) ) * it->second;
    }
    return p/2;
}

void
HMatrix::set(const Allele &a1, const Allele &a2, double d)
{
	// invalidates m_bin
	m_bin.clear();

	if (a1 < a2)
	{
		m_pmf[std::make_pair(a1, a2)] = d;
	}
	else
	{
		m_pmf[std::make_pair(a2, a1)] = d;
	}

	ive_changed();
}

bool
HMatrix::normalize()
{
	// invalidates m_bin
	m_bin.clear();

	double ret = m_pmf.normalize();

	ive_changed();

	return ret;
}

const std::vector<float> &
HMatrix::vec() const
{
	if (m_bin.empty())
	{
		makeBin();
	}
	return m_bin;
}

const float *
HMatrix::bin() const
{
	return &(vec()[0]);
}

void
HMatrix::makeBin() const
{
	m_bin.clear();

	PMF< std::pair<Allele, Allele> >::const_iterator it;
	for(it = m_pmf.begin(); it != m_pmf.end(); ++it)
	{
		m_bin.push_back(it->second);
	}
}

HMatrix&
HMatrix::operator=(const HMatrix &h2)
{
    this->m_pmf = h2.m_pmf;
    ive_changed();
    return *this;
}

HMatrix&
HMatrix::operator+=(const HMatrix &h2)
{
	this->m_pmf += h2.m_pmf;
	ive_changed();
	return *this;
}

double
dot(const HMatrix &h1, const HMatrix &h2)
{
	return dot(h1.m_pmf, h2.m_pmf);
}

// halve the off-diagonal elements of h (sparse)
HMatrix
halveOffDiag(const HMatrix &h)
{
    HMatrix ret;
    for (PMF< std::pair<Allele, Allele> >::const_iterator it = h.m_pmf.begin();
         it != h.m_pmf.end();
         ++it)    // for all i,j in hback
    {
        Allele i = it->first.first;
        Allele j = it->first.second;
        double p = it->second;
        ret.set(i, j, (i==j)? p : p/2);
    }
    return ret;
}

double
hdot(const HMatrix &h1, const HMatrix &h2)
{
    double ret = dot(h1, halveOffDiag(h2) );
    return ret;
}

HMatrix
operator+(const HMatrix &h1, const HMatrix &h2)
{
	HMatrix ret;
	ret.m_pmf = h1.m_pmf + h2.m_pmf;
	ret.ive_changed();
	return ret;
}

HMatrix
operator-(const HMatrix &h1, const HMatrix &h2)
{
	HMatrix ret;
	ret.m_pmf = h1.m_pmf - h2.m_pmf;
    ret.ive_changed();
	return ret;
}

HMatrix operator*(const HMatrix &h1, double d)
{
	HMatrix ret;
	ret.m_pmf = h1.m_pmf * d;
    ret.ive_changed();
	return ret;
}

HMatrix operator/(const HMatrix &h1, double d)
{
    return h1 * (1/d);
}

HMatrix operator*(double d, const HMatrix &h1)
{
	HMatrix ret;
	ret.m_pmf = d * h1.m_pmf;
    ret.ive_changed();
	return ret;
}

HMatrix
operator*(const HMatrix &h1, const HMatrix &h2)
{
	HMatrix ret;
	ret.m_pmf = h1.m_pmf * h2.m_pmf;
    ret.ive_changed();
	return ret;
}

HMatrix
operator/(const HMatrix &h1, const HMatrix &h2)
{
	HMatrix ret;
	ret.m_pmf = h1.m_pmf / h2.m_pmf;
    ret.ive_changed();
	return ret;
}

HMatrix &
operator*=(HMatrix & h, double d)
{
    h.m_pmf *= d;
    h.ive_changed();
    return h;
}

HMatrix &
operator/=(HMatrix & h, double d)
{
    h.m_pmf *= 1/d;
    h.ive_changed();
    return h;
}

HMatrix &
operator+=(HMatrix & h, double d)
{
    h.m_pmf += d;
    h.ive_changed();
    return h;
}

HMatrix &
operator-=(HMatrix & h, double d)
{
    h.m_pmf += -d;
    h.ive_changed();
    return h;
}

BackFreq::BackFreq(const PMF<Allele> &pmf, SubPopModel const &spm)
: HMatrix(pmf, spm)
{
    ive_changed();
}

// probability of a (unconditionally, same as original pmf)
double
BackFreq::operator()(Allele const &a) const
{
    PMF<Allele>::const_iterator it = m_allele_pmf.find(a);

    if (it != m_allele_pmf.end())
    {
        return it->second;
    }
    return 0;
}

// probability of g (allele pair, in either order)
double
BackFreq::operator()(std::pair<Allele, Allele> const &g) const
{
    PMF< std::pair<Allele, Allele> >::const_iterator it = find(g.first, g.second);

    if (it != m_pmf.end())
    {
        return it->second;
    }
    return 0;
}

// probability of g (allele pair, in specific order)
double
BackFreq::pOrdered(std::pair<Allele, Allele> const &g) const
{
    double p = (*this)(g);

    return (g.first == g.second) ? p : p/2;
}

// probability of a given b (conditional)
double
BackFreq::operator()(Allele const &p, Allele const &q) const
{
    return pOrdered(std::make_pair(p, q))/ (*this)(q);
}

PMF<Allele> const &
BackFreq::alleleFreqs() const
{
    return m_allele_pmf;
}

void
BackFreq::ive_changed()
{
    m_allele_pmf.clear();

    PMF< std::pair<Allele, Allele> >::const_iterator it;

    for (it = m_pmf.begin(); it != m_pmf.end(); ++it)
    {
        m_allele_pmf[it->first.first] += it->second / 2;
        m_allele_pmf[it->first.second] += it->second / 2;
    }
}

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
  { A, 0.1 }
, { B, 0.2 }
};

static PMF<Allele>::POD LOCUS2_FREQ[] =
{
  { A, 0.3 }
, { B, 0.4 }
};

TEST(HMatrix0)
{
	HMatrix h1(PMF<Allele>(LOCUS1_FREQ, sizeof(LOCUS1_FREQ)/sizeof(PMF<Allele>::POD)),
			   PMF<Allele>(LOCUS2_FREQ, sizeof(LOCUS2_FREQ)/sizeof(PMF<Allele>::POD)));

	CHECK_EQUAL(0.1 * 0.3,             h1(A, A));
	CHECK_EQUAL(0.1 * 0.4 + 0.2 * 0.3, h1(A, B));
	CHECK_EQUAL(0.1 * 0.4 + 0.2 * 0.3, h1(B, A)); // same element
	CHECK_EQUAL(0.2 * 0.4,             h1(B, B));
}

TEST(HMatrix1)
{
	HMatrix h1;

//	CHECK_THROW(h1(A, A), UnitTest::AssertException);
	// we now allow sparse access (non-existent elements return 0)
	CHECK_EQUAL(0, h1(A, A));

	h1.set(A, A, 0.5);
	CHECK_EQUAL(0.5, h1(A, A));

	h1.set(A, B, 0.5);
	h1.set(A, C, 0.5);
	h1.set(A, D, 0.5);

	// (A,B) and (B,A) are the same element
	CHECK_EQUAL(0.5, h1(A, B));
	CHECK_EQUAL(0.5, h1(B, A));

	h1.normalize();

	CHECK_EQUAL(0.25, h1(A, B));
	CHECK_EQUAL(0.25, h1(B, A));
}

struct HMatrix_fixture
{
	HMatrix h1, h2;

	HMatrix_fixture()
	{
		h1.set(A, A, 1);
		h1.set(A, B, 2);
		h1.set(B, B, 3);

		h2.set(A, A, 8);
		h2.set(A, B, 4);
		h2.set(B, B, 2);
	}
};

TEST_FIXTURE(HMatrix_fixture, HMatrix2)
{
	HMatrix p = h1 * h2;
	CHECK_EQUAL(1.0*8, p(A, A));
	CHECK_EQUAL(2.0*4, p(A, B));
	CHECK_EQUAL(3.0*2, p(B, B));
}

TEST_FIXTURE(HMatrix_fixture, HMatrix3)
{
	HMatrix q = h1 / h2;
	CHECK_EQUAL(1.0/8, q(A, A));
	CHECK_EQUAL(2.0/4, q(A, B));
	CHECK_EQUAL(3.0/2, q(B, B));
}

TEST_FIXTURE(HMatrix_fixture, HMatrix4)
{
	CHECK_EQUAL(1.0*8 + 2.0*4 + 3.0*2, dot(h1, h2));
}

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * AlleleSet.cpp
 *
 *  Created on: Dec 22, 2009
 *      Author: gareth
 */

#include "AlleleSet.h"
#include "loci.h"
#include "match.h"
#include <UnitTest++/UnitTest++.h>

#include "MessageStream.h"
INIT_MESSAGES("AlleleSet")
#include "messages.h"

using namespace std;

AlleleSet::AlleleSet()
: m_pmfs()
, m_error_rate(0)
, m_mutation_rate(0)
, m_theta(0)
, m_num_contributors(0)
{
	// treat as two unknown alleles ('F'): at this point we do not know the background so represent 'F's as empty PMFs
	PMF<Allele> pmf1;

	m_pmfs.push_back(pmf1); // copy
	m_pmfs.push_back(pmf1); // copy

	cacheSimple();
}

AlleleSet::AlleleSet(const Allele a1, const Allele a2)
: m_pmfs()
, m_error_rate(0)
, m_mutation_rate(0)
, m_theta(0)
, m_num_contributors(0)
{
	// Unknown alleles: at this point we do not know the background so represent 'F's as empty PMFs
	PMF<Allele> pmf1, pmf2;

	if (a1 != Allele::unknown)
	{
		pmf1[a1] = 1;
	}

	if (a2 != Allele::unknown)
	{
		pmf2[a2] = 1;
	}

	m_pmfs.push_back(pmf1);
	m_pmfs.push_back(pmf2);

	cacheSimple();
}

AlleleSet::AlleleSet(const Allele a1, const Allele a2, const Allele a3, const Allele a4)
: m_pmfs()
, m_error_rate(0)
, m_mutation_rate(0)
, m_theta(0)
, m_num_contributors(0)
{
	// Unknown alleles: at this point we do not know the background so represent 'F's as empty PMFs
	PMF<Allele> pmf1, pmf2, pmf3, pmf4;

	if (a1 != Allele::unknown)
		pmf1[a1] = 1;

	if (a2 != Allele::unknown)
		pmf2[a2] = 1;

	if (a3 != Allele::unknown)
		pmf3[a3] = 1;

	if (a4 != Allele::unknown)
		pmf4[a4] = 1;

	m_pmfs.push_back(pmf1);
	m_pmfs.push_back(pmf2);
	m_pmfs.push_back(pmf3);
	m_pmfs.push_back(pmf4);

	cacheSimple();
}

AlleleSet::AlleleSet(std::vector< PMF<Allele> > const &pmf_list)
: m_pmfs(pmf_list)
, m_error_rate(0)
, m_mutation_rate(0)
, m_theta(0)
, m_num_contributors(0)
{
	cacheSimple();
}

// check for alleles not in population database
bool
AlleleSet::checkBackground(PMF<Allele> const &background) const
{
    bool ret = true;

	std::vector< PMF<Allele> >::const_iterator ip;
	for(ip = m_pmfs.begin(); ip != m_pmfs.end(); ++ip)
	{
		PMF<Allele>::const_iterator ia;
		for(ia = ip->begin(); ia != ip->end(); ++ia)
		{
			if (background.find(ia->first) == background.end())
			{
				// Allele not in database.

				warn << startl << "allele not in population database: " << ia->first.string() << std::endl;
				ret = false;
			}
		}
	}
	return ret;
}

// return an allele PMF
const PMF<Allele>&
AlleleSet::getAllele(int n) const
{
	static const PMF<Allele> f;

	if ((int)m_pmfs.size() < n+1)
	{
		return f; // empty PMF == 'F'
	}
	return m_pmfs[n];
}

const vector< PMF<Allele> >&
AlleleSet::getAllAlleles() const
{
    return m_pmfs;
}


// set an allele PMF
void
AlleleSet::setAllele(int n, PMF<Allele> const &pmf)
{
	if ((int)m_pmfs.size() < n+1 )
	{
		m_pmfs.resize(n+1);
	}
	m_pmfs[n] = pmf;

	reset();
	cacheSimple();
}

// construct PMF with the given background and delta, from a PMF representing the observed profile.
// if pmf unknown (empty) return the background
// NB The input PMF need NOT be normalized. The returned PMF is normalized.
//
static PMF<Allele>
makePMF(const PMF<Allele> &pmf,
        const PMF<Allele> &background,
        float delta,
        bool sparse)
{
	// How much background to add in?
	// If the input PMF is not normalized then we make up the difference
	// with background. Otherwise we add in delta.
	double bglevel = std::max((double)delta, 1 - pmf.sum());

	PMF<Allele> ret;
	if (!sparse)
	{
		// This implementation ensures there is an element of the
		// returned PMF for each element of the background
		ret = background;
		if (! pmf.empty())
		{
			ret *= bglevel;
            ret += pmf;
			ret.normalize();
		}

		Assert2(ret.size() == background.size(), "PMF contains unknown alleles");
	}
	else
	{
		// This implementation allows the PMF (and therefore the HMatrix) to remain sparse
		if (pmf.empty())
		{
			ret = background;
		}
		else if (bglevel == 0)
		{
			ret = pmf;
			ret.normalize();
		}
		else
		{
			ret = background;
			ret *= bglevel;
            ret += pmf;
			ret.normalize();
		}
	}

//    cout << "pmf = " << pmf << endl;
//    cout << "background = " << background << endl;
//    cout << "bglevel = " << bglevel << endl;
//    cout << "ret = " << ret << endl << endl;
	return ret;
}

// term of non-HW HMatrix
inline double term_ij(
        double Bp,
        double Bq,
        const PMF<Allele> &p,
        const PMF<Allele> &q,
        const BackFreq &hback,
        Allele const &i,
        Allele const &j)
{
    // NB in this term (i, j) is ordered

    double ret    = p.val(i) * q.val(j);
    if (Bq>0)       ret += Bq * p.val(i) * hback(j, i);
    if (Bp>0)       ret += Bp * q.val(j) * hback(i, j);
    if (Bp*Bq>0)    ret += Bp * Bq * hback.pOrdered(make_pair(i,j));
    return ret;
}

// This gives non-HW treatment of the 'F' term
HMatrix
makeHMatrixNHW(
        const PMF<Allele> &p,
        const PMF<Allele> &q,
        const BackFreq &hback,
        float delta,
        bool sparse)
{
    // How much background to add in?
    // If the input HMatrix is not normalized then we make up the difference
    // with background. Otherwise we add in delta.
    double Bp = std::max((double)delta, 1 - p.sum()); // background for p
    double Bq = std::max((double)delta, 1 - q.sum()); // background for q

    // apply this formula:
    // H(ij) = p(i)q(j) + Bq p(i)b(j|i) + Bp (q(j)b(i|j) + Bp Bq Bij

    HMatrix ret;

    // loop over all elements in the background. This gives us the upper
    // triangular terms only. We need to sum over all terms.

    HMatrix &href = (HMatrix&)hback; // to use base class member
    PMF< std::pair<Allele, Allele> >::iterator it;
    for (it = href.m_pmf.begin(); it != href.m_pmf.end(); ++it)
    {
        Allele i = it->first.first;
        Allele j = it->first.second;

        double h_ij = term_ij(Bp, Bq, p, q, hback, i, j); // upper-triangular term

        if (i != j)
        {
            h_ij += term_ij(Bp, Bq, p, q, hback, j, i); // lower-triangular term
        }

        if (!sparse || h_ij > 0) ret.set(i, j, h_ij);
    }

    ret.normalize(); // just in case
    return ret;
}

// HW case
void
AlleleSet::makeMatchHMatrix(
        const PMF<Allele> &background,
        bool sparse) const
{
    Assert(m_num_contributors > 0);

    HMatrix& match = ( sparse ? m_match_sparse : m_match);
    match.clear();

    if (m_num_contributors == 1)
    {
        // single contributor - two alleles (we do not yet handle tri-allele)
        if (m_pmfs.size() != 2)
        {
            error << startl << "should be 2 alleles but there are " <<  m_pmfs.size() << endl;
            exit_debug(-1);
        }

        match = HMatrix(makePMF(m_pmfs[0], background, m_error_rate, sparse),
                        makePMF(m_pmfs[1], background, m_error_rate, sparse));

        // NB the product of two normalized PMFs must be normalized
        // match.normalize();
    }
    else if ((int)m_pmfs.size() == 2 * m_num_contributors)
    {
        // mixture - two alleles per contributor (some may be duplicated, but we know which and how many times)
        int n = m_pmfs.size();

        for (int i=0; i<n; ++i)
        {
            for (int j=i+1; j<n; ++j)
            {
                match += HMatrix(makePMF(m_pmfs[i], background, m_error_rate, sparse),
                                 makePMF(m_pmfs[j], background, m_error_rate, sparse));
            }
        }

        match.normalize();
    }
    else
    {
        // we know only which alleles are present - not how many copies of each
        // we need the general method  - like Fung & Hu

        // Assume the alleles we don't know are just a background mixture of the known alleles.
        // (NB this treatment is not quite equivalent to Fung & Hu. because the conditional probabilities
        // are not exactly equal to the background frequencies).

        // construct the mixture
        int n = m_pmfs.size();

        PMF<Allele> missing_allele;

        for (int i=0; i<n; ++i)
        {
            missing_allele += m_pmfs[i] * background;
        }
        missing_allele.normalize();

        // add in copies of it
        std::vector<PMF<Allele> > pmfs_copy(m_pmfs);
        for (int i = n; i < 2 * m_num_contributors; ++i)
        {
            pmfs_copy.push_back(missing_allele);
        }

        // proceed as before
        n = pmfs_copy.size();
        for (int i=0; i<n; ++i)
        {
            for (int j=i+1; j<n; ++j)
            {
                match += HMatrix(makePMF(pmfs_copy[i], background, m_error_rate, sparse),
                                 makePMF(pmfs_copy[j], background, m_error_rate, sparse));
            }
        }

        match.normalize();
    }
}

// Non-HW case
HMatrix
AlleleSet::identNHW(
        const BackFreq &hback,
        bool sparse) const
{
	Assert(m_num_contributors > 0);

	HMatrix match;

	if (m_num_contributors == 1)
	{
		// single contributor - two alleles (we do not yet handle tri-allele)
		if (m_pmfs.size() != 2)
		{
			error << startl << "should be 2 alleles but there are " <<  m_pmfs.size() << endl;
			exit_debug(-1);
		}

		match = makeHMatrixNHW(m_pmfs[0], m_pmfs[1], hback, m_error_rate, sparse);

		// NB the product of two normalized PMFs must be normalized
		// match.normalize();
	}
	else if ((int)m_pmfs.size() == 2 * m_num_contributors)
	{
		// mixture - two alleles per contributor (some may be duplicated, but we know which and how many times)
		int n = m_pmfs.size();

		for (int i=0; i<n; ++i)
		{
			for (int j=i+1; j<n; ++j)
			{
			    match += makeHMatrixNHW(m_pmfs[i], m_pmfs[j], hback, m_error_rate, sparse);
			}
		}

		match.normalize();
	}
	else
	{
		// we know only which alleles are present - not how many copies of each
		// we need the general method  - like Fung & Hu

		// Assume the alleles we don't know are just a background mixture of the known alleles.
		// (NB this treatment is not quite equivalent to Fung & Hu. because the conditional probabilities
		// are not exactly equal to the background frequencies).

		// construct the mixture
		int n = m_pmfs.size();

		PMF<Allele> missing_allele;

		for (int i=0; i<n; ++i)
		{
            missing_allele += m_pmfs[i] * hback.alleleFreqs(); // multiply corresponding entries, sparsely
		}
		missing_allele.normalize();

		// add in copies of it until we have the right number of alleles
		std::vector<PMF<Allele> > pmfs_copy(m_pmfs);
		for (int i = n; i < 2 * m_num_contributors; ++i)
		{
			pmfs_copy.push_back(missing_allele);
		}

		// proceed as before
		n = pmfs_copy.size();
		for (int i=0; i<n; ++i)
		{
			for (int j=i+1; j<n; ++j)
			{
                match += makeHMatrixNHW(pmfs_copy[i], pmfs_copy[j], hback, m_error_rate, sparse);
			}
		}

		match.normalize();
	}

	return match;
}

void
AlleleSet::makeSibHMatrix(
        HMatrix const &prof,
        const PMF<Allele> &background,
        bool sparse) const
{
	HMatrix& hsib = ( sparse ? m_sib_sparse : m_sib);
	hsib = sib(prof, background); // NB sib() is declared as a friend of HMatrix and defined in match.cpp
}

const HMatrix &
AlleleSet::getMatchHMatrixNHW(
        const BackFreq &hback,
        const ProfileData &data,
        bool sparse) const
{
    setData(data);

    if (m_match_nhw.empty())
    {
        m_match_nhw = identNHW(hback, sparse);
    }
    return m_match_nhw;
}

const HMatrix &
AlleleSet::getMatchHMatrix(
        const PMF<Allele> &background,
        const ProfileData &data,
        bool sparse) const
{
    setData(data);

    HMatrix& hmatch = ( sparse ? m_match_sparse : m_match);
    if (hmatch.empty())
    {
        makeMatchHMatrix(background, sparse);
    }
    return hmatch;
}

const HMatrix &
AlleleSet::getSibHMatrix(
        HMatrix const &prof,
        const PMF<Allele> &background,
        const ProfileData &data,
        bool sparse) const
{
    setData(data);

	HMatrix& hsib = ( sparse ? m_sib_sparse : m_sib);
	if (hsib.empty())
	{
		makeSibHMatrix(prof, background, sparse);
	}
	return hsib;
}

// HW case
HMatrix
AlleleSet::getRelHMatrix(
        const MatchType &match_type,
        const PMF<Allele> &background,
        const ProfileData &data,
        bool sparse) const
{
	// Todo caching ?

    // Ensure match matrix exists (and call setData)
    (void)getMatchHMatrix(background, data, sparse);

    HMatrix *hmatch = &( sparse ? m_match_sparse : m_match); // point to cached matrix

    switch(match_type.m_rel_type)
    {
    case sibling_t:
        return getSibHMatrix(*hmatch, background, data, sparse); // cached

    case degree_1_t:
    case degree_2_t:
    case degree_pq_t:
        return degree(match_type.m_path1steps, match_type.m_path2steps, *hmatch, background);

    case gen_t:
        return genRel(match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, *hmatch, background);

    case inv_t:
    {
        return invRel(match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, *hmatch, background);
    }
    default:
        Assert2(false, "getRelHMatrix called with unknown MatchType");
        break;
    }
} // Compiler warning: It's OK. Honest.

const HMatrix &
AlleleSet::getSpmcHMatrix(
        const ProfileData &data,
        const PMF<Allele> &background,
        SubPopModel const &spm) const
{
    Assert(spm.type == SubPopModel::B11);

    setTheta(spm.theta_bar);

    if (m_spmc.empty())
    {
        BackFreq hback(background, spm);
        HMatrix hI = getMatchHMatrixNHW(hback, data, true);
        m_spmc = FSTCorrection(hI, hback, m_theta);
    }
    return m_spmc;
}

void AlleleSet::cacheSimple()
{
    // An AlleleSet is simple if it has exactly two alleles,
    // known with certainty (no component of 'F').
	m_simple = (m_pmfs.size() == 2) &&
	           (m_pmfs[0].size() == 1) &&
   	           (m_pmfs[0].begin()->second == 1) &&
   	           (m_pmfs[1].size() == 1) &&
	           (m_pmfs[1].begin()->second == 1);

	if (m_simple)
	{
		m_a1 = m_pmfs[0].begin()->first;
		m_a2 = m_pmfs[1].begin()->first;
		if (m_a2<m_a1) std::swap(m_a1, m_a2);
	}
}

bool AlleleSet::simple(Allele &a1, Allele &a2) const
{
	if (m_simple)
	{
		a1 = m_a1;
		a2 = m_a2;
	}
	return m_simple;
}

void AlleleSet::setErrorRate(double val) const
{
    if (val != m_error_rate)
    {
        m_error_rate = val;

        // reset all cached matrices
        reset();
    }
}

void AlleleSet::setMutRate(double val) const
{
    if (val != m_mutation_rate)
    {
        m_mutation_rate = val;

        // reset cached relationship matrices
        m_sib.clear();
        m_sib_sparse.clear();
    }
}

void AlleleSet::setTheta(double val) const
{
    if (val != m_theta)
    {
        m_theta = val;

        // reset subpop correction matrix
        m_spmc.clear();
    }
}

void AlleleSet::setNumContribs(int val) const
{
    if (val != m_num_contributors)
    {
        m_num_contributors = val;

        // reset all cached matrices
        reset();
    }
}

void AlleleSet::setData(const ProfileData &data) const
{
    setErrorRate(data.m_error_rate);
    setMutRate(data.m_mutation_rate);
    setNumContribs(data.m_num_contributors);
}

void AlleleSet::reset() const
{
    m_match.clear();
    m_sib.clear();
    m_match_sparse.clear();
    m_sib_sparse.clear();
    m_spmc.clear();
    m_match_nhw.clear();
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
  { C, 0.3 }
, { D, 0.4 }
};

TEST(AlleleSet0)
{
	AlleleSet as1;
	CHECK(true);
}

TEST(AlleleSet1)
{
	AlleleSet as1(A, B);
	CHECK(true);
}

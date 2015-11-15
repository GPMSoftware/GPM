/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * match.cpp
 *
 *  Created on: Dec 9, 2009
 *      Author: gareth
 *
 *
 */

#include "match.h"
#include "loci.h"
#include "Assert.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("match")
#include "fand/messages.h"

#include <math.h>

#define SPARSE_OPTIMIZATION

int SIMPLE_OPTIMIZATION = 1; // global so we can turn it off for testing

// delta function
#define d(a, b) ((int)(a==b))

// Beta-binomial sampling formula
// P(allele a | y a among n alleles)
double BBSF(int y,
            int n,
            double fst,
            double f)   // background frequency P(a)
{
    return (y * fst + (1-fst) * f) / (1 + (n-1) * fst);
}

// probability a new allele is p given i, j observed
double
probBNp_giv_ij(
        Allele const &p,
        Allele const &i,
        Allele const &j,
        double fp,
        double fst)
{
    // P(p|ij)

    int y = delta(p,i) + delta(p,j);
    return BBSF(y, 2, fst, fp);
}

// probability p and a new allele is i, j
double
probBNpX_is_ij(
    Allele const &p,
    Allele const &q,
    Allele const &i,
    Allele const &j,
    double fi,
    double fj,
    double fst)
{
    // P(pX == ij) | pq where X is drawn from the background

    if (p==i)
    {
        return probBNp_giv_ij(j,p,q,fj,fst); // P(j|pq)
    }
    else if (p==j)
    {
        return probBNp_giv_ij(i,p,q,fi,fst); // P(i|pq)
    }
    else
    {
        return 0;
    }
}

// probability two new alleles are p,q given i,j observed
// (NB this is a generalized version of NRC4_10)
double
probBNpq_giv_ij(
        Allele const &p,
        Allele const &q,
        Allele const &i,
        Allele const &j,
        double fp,
        double fq,
        double fst)
{
    // P(pq|ij) = P(p|ij) * P(q|pij) : *2 if p!=q

    int y1 = delta(p,i) + delta(p,j);
    int y2 = delta(q,p) + delta(q,i) + delta(q,j);

    double x = BBSF(y1, 2, fst, fp) * BBSF(y2, 3, fst, fq);

    if (p != q) x *= 2;

    return x;
}

// forward declaration
double baldingNichols(
        MatchType const &mtype,
        Profile const &p1,
        Profile const &p2,
        PopulationData const &popdata,
        double theta);

char *cpu_options = "\"CPU running General algorithm\"";

//#include "CheckMacros.h" // hacked version of UnitTest++ header
#include <UnitTest++/UnitTest++.h>

namespace UnitTest
{

template< typename Expected, typename Actual >
void CheckGreater(TestResults& results, Expected const& expected, Actual const& actual, TestDetails const& details)
{
    if (!(actual > expected))
    {
        UnitTest::MemoryOutStream stream;
        stream << "Expected >" << expected << " but was " << actual;

        results.OnTestFailure(details, stream.GetText());
    }
}

template< typename Expected, typename Actual >
void CheckLess(TestResults& results, Expected const& expected, Actual const& actual, TestDetails const& details)
{
    if (!(actual < expected))
    {
        UnitTest::MemoryOutStream stream;
        stream << "Expected <" << expected << " but was " << actual;

        results.OnTestFailure(details, stream.GetText());
    }
}

#define CHECK_GREATER(expected, actual) \
    do \
    { \
        try { \
            UnitTest::CheckGreater(*UnitTest::CurrentTest::Results(), expected, actual, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        } \
        catch (...) { \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_EQUAL(" #expected ", " #actual ")"); \
        } \
    } while (0)

#define CHECK_LESS(expected, actual) \
    do \
    { \
        try { \
            UnitTest::CheckLess(*UnitTest::CurrentTest::Results(), expected, actual, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        } \
        catch (...) { \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_EQUAL(" #expected ", " #actual ")"); \
        } \
    } while (0)

} // namespace UnitTest

using namespace std;


static float tol = 1e-4;


MatchType::MatchType(RelType rel, int p1, int p2)
{
    Assert2(rel != gen_t, "wrong constructor called for gen_t");
    Assert2(rel != inv_t, "wrong constructor called for inv_t");
    Assert(p1>=0);
    Assert(p2>=p1);
    Assert(!(p1==1 && p2==1));

	m_rel_type = rel;
	m_path1steps = p1;
	m_path2steps = p2;

	// ensure consistent representation
	switch(m_rel_type)
	{
	case ident_t:
		m_path1steps = 0;
		m_path2steps = 0;
		break;
	case sibling_t:
		m_path1steps = -2; // sibling is different from other 2,2 relationships
		m_path2steps = -2;
		break;
	case degree_1_t:
		m_path1steps = 1;
		m_path2steps = INF;
		break;
	case degree_2_t:
		m_path1steps = 2;
		m_path2steps = INF;
		break;
	case degree_pq_t:
		if ((m_path1steps == 0) && (m_path2steps == 0))
		{
			m_rel_type = ident_t;
		}
		else if ((m_path1steps == 1) && (m_path2steps == INF))
		{
			m_rel_type = degree_1_t;
		}
		else if ((m_path1steps == 2) && (m_path2steps == INF))
		{
			m_rel_type = degree_2_t;
		}

		Assert(m_path1steps>=0);
		Assert(m_path2steps>=m_path1steps);
		Assert(!(m_path1steps==1 && m_path2steps==1));
		break;
	default:
		m_path1steps = INF;
		m_path2steps = INF;
	}
}

MatchType::MatchType(RelType rel, float a1, float b1, float a2, float b2)
{
    m_rel_type = rel;
    Assert(rel == gen_t || rel == inv_t);

    m_a1 = a1;
    m_b1 = b1;
    m_a2 = a2;
    m_b2 = b2;
}

// Dnm with fractional n, m
// use -1 for infinity
MatchType::MatchType(RelType rel, double n, double m)
{
    m_rel_type = rel;
    Assert(rel == gen_t || rel == inv_t);

    m_a1 = (n<0) ? 0 : pow(0.5, n);
    m_b1 = m_a1;
    m_a2 = (m<0) ? 0 : pow(0.5, m);
    m_b2 = m_a2;
}

std::string
MatchType::string() const
{
	switch(m_rel_type)
	{
	case ident_t:
		return "IDENT";
	case sibling_t:
		return "SIB";
	case degree_1_t:
		return "D1";
	case degree_2_t:
		return "D2";
	case degree_pq_t:
	{
		ostringstream oss;

		oss << "D(" << m_path1steps;
		if (m_path2steps != INF)
		{
			oss <<  "," << m_path2steps;
		}
		oss << ")";

		return oss.str();
	}
	case gen_t:
	{
        ostringstream oss;

        oss << "GEN[" << m_a1 << "," << m_b1 << "," << m_a2 << "," << m_b2 << "]";

        return oss.str();
	}
    case inv_t:
    {
        ostringstream oss;

        oss << "INV[" << m_a1 << "," << m_b1 << "," << m_a2 << "," << m_b2 << "]";

        return oss.str();
    }
    case none_t:
    {
        return "unrelated";
    }
	default:
		Assert2(false, "Unknown RelType");
	}
	return "ERROR";
}

int MatchType::generations() const
{
    // ensure consistent representation
    switch(m_rel_type)
    {
    case ident_t:
        return 1;
        break;
    case sibling_t:
        return 1;
        break;
    case degree_1_t:
        return 1;
        break;
    case degree_2_t:
        return 2;
        break;
    case degree_pq_t:
        return m_path1steps; // TODO: this is an approximation for bilinear relationships
        break;
    default:
        Assert2(false, "unknown RelType in MatchType::generations");
    }
}

ostream &
operator<<(ostream &os, const MatchType &a)
{
	os << a.string();
	return os;
}

bool
operator<(const MatchType &m1, const MatchType &m2)
{
    if (m1.m_rel_type != m2.m_rel_type)
    {
        return (m1.m_rel_type < m2.m_rel_type);
    }
    else if (m1.m_rel_type == degree_pq_t)
    {
        return (m1.m_path1steps != m2.m_path1steps) ? (m1.m_path1steps < m2.m_path1steps) :
               (m1.m_path2steps != m2.m_path2steps) ? (m1.m_path2steps < m2.m_path2steps) : false;
    }
    else if (m1.m_rel_type == gen_t || m1.m_rel_type == inv_t)
    {
        return (m1.m_a1 != m2.m_a1) ? (m1.m_a1 < m2.m_a1) :
               (m1.m_b1 != m2.m_b1) ? (m1.m_b1 < m2.m_b1) :
               (m1.m_a2 != m2.m_a2) ? (m1.m_a2 < m2.m_a2) :
               (m1.m_b2 != m2.m_b2) ? (m1.m_b2 < m2.m_b2) : false;
    }
    else
    {
        return false; // identical
    }
}

// construct PMF of degree-1 relative of prof (general method)
// i.e parent or child
HMatrix
degree1(HMatrix const &prof,
	const PMF<Allele> &back)
{
	HMatrix ret;

	// for each element in the HMatrix construct the degree1, and sum them

	for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
		 it != prof.m_pmf.end();
		 ++it)
	{
		Allele a1 = it->first.first;
		Allele a2 = it->first.second;

		PMF<Allele> p1;
		p1[a1] += 0.5;
		p1[a2] += 0.5;

		PMF<Allele> p2 = back;

		HMatrix h_bit(p1, p2);
		ret.m_pmf += h_bit.m_pmf * it->second;
	}

	return ret;
}

// construct PMF of degree-N relative of prof (general method)
// i.e (n-2)greats-grandchild[grandparent]
HMatrix
degree(
	int n,                    // degree of relationship
	HMatrix const &prof,
	const PMF<Allele> &back)
{
	Assert(n >= 0);

	if (n == 0)
	{
		return prof;
	}

#if 0
	// recursively from degree-1 : TEST ONLY (will be inefficient)
	return degree(n-1, degree1(prof, back), back);

#else

	// directly
	HMatrix ret;

	for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
		 it != prof.m_pmf.end();
		 ++it)
	{
		Allele a1 = it->first.first;
		Allele a2 = it->first.second;

		PMF<Allele> p1;
		p1[a1] += 1;
		p1[a2] += 1;

		double Cr = pow(2, -n);
		p1 *= Cr;

		p1 += back * (1 - 2*Cr);

		PMF<Allele> p2 = back;

		HMatrix h_bit(p1, p2);
		ret.m_pmf += h_bit.m_pmf * it->second;
	}

	return ret;
#endif
}

// calculate parental contribution as PMF<Allele>
void
pcontrib(
	HMatrix const &prof,   // parent
	PMF<Allele>& contrib)  // contribution to offspring
{
	contrib.clear();

	for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
		 it != prof.m_pmf.end();
		 ++it)
	{
		Allele a1 = it->first.first;
		Allele a2 = it->first.second;

		double p = 0.5 * it->second;
		contrib[a1] += p;
		contrib[a2] += p;
	}
}

// construct PMF of full sibling of prof (general method)
HMatrix
childOf(
	HMatrix const &prof1,
	HMatrix const &prof2)
{
	// combine contributions from each parent
	PMF<Allele> p1, p2;
	pcontrib(prof1, p1);
	pcontrib(prof2, p2);
	return HMatrix(p1, p2);
}

HMatrix
degree(
	int p, int q,        // degree of relationship
	HMatrix const &prof,
	const PMF<Allele> &back)
{
	Assert(p>=0 && q>=p);
	Assert(!(p==1 && q==1));

	if (q == MatchType::INF || q == 0)
	{
		return degree(p, prof, back);
	}
	else
	{
		// a (p,q) is calculated as the child of a (p-1, 0) and a (q-1, 0)
		HMatrix p_contrib = degree(p-1, prof, back);
		HMatrix q_contrib = degree(q-1, prof, back);
		return childOf(p_contrib, q_contrib);
	}
}

// construct SIB on Balding-Nichols model
double
sibBN(
      Allele const &p,
      Allele const &q,
      Allele const &r,
      Allele const &s,
      double fr,
      double fs,
      double fst)
{
    // return P(SIB(pq)==rs)

    return 0.25 * ( delta(p,r) * delta(q,s)            // P(pq == rs)
                  + probBNpX_is_ij(p,q,r,s,fr,fs,fst)  // P(p? == rs)
                  + probBNpX_is_ij(q,p,r,s,fr,fs,fst)  // P(q? == rs)
                  + probBNpq_giv_ij(r,s,p,q,fr,fs,fst) // P(?? == rs) == P(rs|pq)
                  );
}

// construct Dn on Balding-Nichols model
double
DnBN(int n,
     Allele const &p,
     Allele const &q,
     Allele const &r,
     Allele const &s,
     double fr,
     double fs,
     double fst)
{
    // return P(Dn(pq)==rs)

    double a = pow(2, n);

    return (1/a) * ( 0                                          // P(pq == rs)
                   + probBNpX_is_ij(p,q,r,s,fr,fs,fst)          // P(p? == rs)
                   + probBNpX_is_ij(q,p,r,s,fr,fs,fst)          // P(q? == rs)
                   + (a-2) * probBNpq_giv_ij(r,s,p,q,fr,fs,fst) // P(?? == rs) == P(rs|pq)
                   );
}

// Construct relative on Balding-Nichols model
// Return Hij such that Hij = P(Rel(pq) = ij)
HMatrix
relBN(MatchType type,
      Allele const &p,
      Allele const &q,
      const PMF<Allele> &back,
      double fst)
{
    HMatrix ret;

    PMF<Allele>::const_iterator it1, it2;
    for (it1 = back.begin(); it1 != back.end(); ++it1)
    {
        Allele i = it1->first;

        for (it2 = it1; it2 != back.end(); ++it2)
        {
            Allele j = it2->first;
            Assert (!(j < i));

            double val = 0;
            switch(type.m_rel_type)
            {
            case sibling_t:
            {
                val = sibBN(p, q, i, j, back.val(i), back.val(j), fst);
                break;
            }
            case degree_pq_t:
            {
                Assert(type.m_path2steps == MatchType::INF);
            }
            case degree_1_t:
            case degree_2_t:
            {
                val = DnBN(type.m_path1steps, p, q, i, j, back.val(i), back.val(j), fst);
                break;
            }
            default:
                Assert2(false, "Unknown match type in relBN()");
            }

            ret.m_pmf[std::make_pair(i, j)] = val;
        }
    }

    return ret;
}

// Construct relative on Balding-Nichols model
HMatrix
relBN(MatchType type,
      HMatrix const &prof,
      const PMF<Allele> &back,
      double fst)
{
    HMatrix ret;

    // for each element in the HMatrix construct a relative, and sum them

    for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
         it != prof.m_pmf.end();
         ++it)
    {
        Allele p = it->first.first;
        Allele q = it->first.second;
        Assert (!(q < p));

        double w = it->second;

        HMatrix Hpq = relBN(type, p, q, back, fst);
        ret.m_pmf += Hpq.m_pmf * w;
    }

    return ret;
}
template <typename T>
std::pair<T,T>
make_ordered_pair(T const &a, T const &b)
{
    return (a<b) ? std::pair<T,T>(a, b) : std::pair<T,T>(b, a);
}

// Stepwise mutation model
HMatrix
mutate_1s(
      HMatrix const &prof,
      const PMF<Allele> &background,
      double r)
{
    HMatrix ret;

    // perform a convolution with a kernel that looks like this
    //
    //     0    r/2   0
    //     r/2  1-2r  r/2
    //     0    r/2   0
    //
    // TODO: FFT and multiply?
    //

    for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
         it != prof.m_pmf.end();
         ++it)
    {
        Allele p = it->first.first;
        Allele q = it->first.second;
        Assert (!(q < p));

        Allele p_plus_1  = Allele(p.m_repeats + 1, p.m_variant);
        Allele p_minus_1 = Allele(p.m_repeats - 1, p.m_variant);
        Allele q_plus_1  = Allele(q.m_repeats + 1, q.m_variant);
        Allele q_minus_1 = Allele(q.m_repeats - 1, q.m_variant);

        // If allele is not in the pop database, ignore it
        // (there is nothing for it to match with)

        double w = it->second;
        double s = r*w/2, t = w;

        if (background.find(p_plus_1) != background.end())
        {
            ret.m_pmf[make_ordered_pair(p_plus_1, q)] += s;
            t -= s;
        }

        if (background.find(p_minus_1) != background.end())
        {
            ret.m_pmf[make_ordered_pair(p_minus_1, q)] += s;
            t -= s;
        }

        if (background.find(q_plus_1) != background.end())
        {
            ret.m_pmf[make_ordered_pair(p, q_plus_1)]  += s;
            t -= s;
        }

        if (background.find(q_minus_1) != background.end())
        {
            ret.m_pmf[make_ordered_pair(p, q_minus_1)] += s;
            t -= s;
        }

        ret.m_pmf[make_ordered_pair(p, q)] += t;
    }

    return ret;
}

// construct PMF of full sibling of prof (HW)
HMatrix
sib(HMatrix const &prof,
	const PMF<Allele> &back)
{
	HMatrix ret;

#if 0
	// for each element in the HMatrix construct a sibling, and sum them

	for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
		 it != prof.m_pmf.end();
		 ++it)
	{
		// contribution of one parent
		Allele a1 = it->first.first;
		PMF<Allele> p1 = back * 0.5;
		p1[a1] += 0.5;

		// contribution of other parent
		Allele a2 = it->first.second;
		PMF<Allele> p2 = back * 0.5;
		p2[a2] += 0.5;

		HMatrix sib_bit(p1, p2);
		ret.m_pmf += sib_bit.m_pmf * it->second;
	}
#else
	// An equivalent method that may be faster
	// NB this is based on algebraically reducing the expression for the PDF.
	// See Study Report (DGW Software) for details

	PMF<Allele> Pq; // Sum_q Pqj
	PMF<Allele> Pr; // Sum_r Pir

	for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
		 it != prof.m_pmf.end();
		 ++it)
	{
		Allele i = it->first.first;
		Allele j = it->first.second;
		double p = it->second;

		Pq[j] += p;
		Pr[i] += p;
	}

	// this needs to be a sum over all i, j, not just where Pij != 0
	for (PMF<Allele>::const_iterator iti = back.begin(); iti != back.end(); ++iti)
	{
		for (PMF<Allele>::const_iterator itj = back.begin(); itj != back.end(); ++itj)
		{

			double Sij = 0;

			Allele i = iti->first;
			Allele j = itj->first;

			// b(i)b(j)
			Sij += iti->second * itj->second;

			// b(i) * Sum_q Pqj
			Sij += iti->second * Pq[j];

			// b(j) * Sum_r Pir
			Sij += itj->second * Pr[i];

			// Pij. NB since this is an HMatrix take the value as 0 below the diagonal
			// - do not allow reflection by operator(i,j)!
			if (i<j || i == j)
				Sij += prof(i,j);

			Sij *= 0.25;

			if (ret.find(i,j) != ret.m_pmf.end())
			{
				Sij += ret(i,j);
			}
			ret.set(i, j, Sij);
		}
	}

#endif

	return ret;
}

PMF<Allele>
genVec(
    float wi,
    float wj,
    Allele const &i,
    Allele const &j,
    const PMF<Allele> &back)
{
    float wf = 1 - wi - wj;     // weight of background

    PMF<Allele> ret = back * wf;

    ret[i] += wi;
    ret[j] += wj;

    return ret;
}

HMatrix
genRel(
    float a1,
    float b1,
    float a2,
    float b2,
    HMatrix const &prof,
    const PMF<Allele> &back)
{
    HMatrix ret;

    for (PMF< std::pair<Allele, Allele> >::const_iterator it = prof.m_pmf.begin();
         it != prof.m_pmf.end();
         ++it)
    {
        Allele i = it->first.first;
        Allele j = it->first.second;

        PMF<Allele> p1 = genVec(a1, b1, i, j, back);
        PMF<Allele> p2 = genVec(a2, b2, i, j, back);

        PMF<Allele> p1_s = genVec(b1, a1, i, j, back);
        PMF<Allele> p2_s = genVec(b2, a2, i, j, back);

        ret += ( HMatrix(p1, p2) + HMatrix(p1_s, p2_s) ) * (it->second / 2.0);
    }

    return ret;
}

// Make the inverse. (Use Bayes theorem - see the maths paper).
HMatrix
invRel(
    float a1,
    float b1,
    float a2,
    float b2,
    HMatrix const &prof,
    const PMF<Allele> &back)
{
    HMatrix ret;

    for (PMF<Allele>::const_iterator i = back.begin();
         i != back.end();
         ++i)
    {
        for (PMF<Allele>::const_iterator j = i;
             j != back.end();
             ++j)
        {
            double b_ij = i->second * j->second * (i==j ? 1 : 2);
            double ret_ij = 0;

            HMatrix AiAj;
            AiAj.set(i->first, j->first, 1);

            HMatrix rel = genRel(a1, b1, a2, b2, AiAj, back);

            for (PMF< std::pair<Allele, Allele> >::const_iterator it = rel.m_pmf.begin();
                 it != rel.m_pmf.end();
                 ++it)
            {
                Allele p = it->first.first;
                Allele q = it->first.second;

                double b_pq = back(p) * back(q) * (p==q ? 1 : 2);

                ret_ij += prof(p,q) * rel(p,q) / b_pq;
            }

            ret.set(i->first, j->first, ret_ij * b_ij);
        }
    }

    return ret;
}

// used in unit tests only
// based on analysis
double
delta1_exact(
        int p,      // 1st allele, 1st profile (the one with the delta attached)
        int q,      // 2nd allele, 1st profile
        int r,      // 1st allele, 2nd profile
        int s,      // 2nd allele, 2nd profile
        double b_r, // background frequency, r
        double b_s, // background frequency, s
        double k)   // delta value (for 1st profile)
{
    double k_ = 1-k;

    return (k_ * k_) * ( d(r,p)*d(s,q) + d(s,p)*d(r,q) ) / (2 * b_r * b_s) +
           (k  * k_) * ( (d(r,p) + d(r,q)) / (2*b_r) + (d(s,p) + d(s,q)) / (2*b_s) ) +
           (k  * k);
}

// used in unit tests only
// Based on analysis
// Strange but true: we need the backgrounds for just one profile
// (they are used only when the same as one or other of the
double
delta2_exact(
        int p,      // 1st allele, 1st profile
        int q,      // 2nd allele, 1st profile
        int r,      // 1st allele, 2nd profile
        int s,      // 2nd allele, 2nd profile
        double b_p, // background frequency, p
        double b_q, // background frequency, q
        double k1,  // delta value, 1st profile
        double k2)  // delta value, 2nd profile
{
    double k1_ = 1-k1;
    double k2_ = 1-k2;
    double k1k2 = k1 * k2;
    double k1_k2_ = k1_ * k2_;

    double A = k1_k2_ * k1_k2_;
    double B = k1_k2_ * (1 - k1_k2_);
    double C = (k1*k1) + (k2*k2) + (2*k1k2*(1 - k1 - k2 + (k1k2/2)));

    return A/2 * ( (d(p,r)*d(q,s)) + (d(q,r)*d(p,s)) ) / (b_p * b_q) +
           B/2 * ( ((d(p,r) + d(p,s)) / (b_p)) + ((d(q,r) + d(q,s)) / (b_q))) +
           C;
}

// First order in delta - based on the above.
// used in unit tests only
double
delta2_approx(
        int p,      // 1st allele, 1st profile
        int q,      // 2nd allele, 1st profile
        int r,      // 1st allele, 2nd profile
        int s,      // 2nd allele, 2nd profile
        double b_p, // background frequency, p
        double b_q, // background frequency, q
        double k1,  // delta value, 1st profile
        double k2)  // delta value, 2nd profile
{
    double k = k1 + k2;

    double A = 1 - 2*k;
    double B = k;
    double C = k * k;

    return A/2 * ( d(p,r)*d(q,s) + d(q,r)*d(p,s) ) / (b_p * b_q) +
           B/2 * ( (d(p,r) + d(p,s)) / (b_p) + (d(q,r) + d(q,s)) / (b_q)) +
           C;
}

// simple match
static double
match_simple(
      const Allele &s1_a1,
	  const Allele &s1_a2,
	  const Allele &s2_a1,
	  const Allele &s2_a2,
	  PMF<Allele> &back,
	  double d1,
	  double d2,
	  SubPopModel const &spm)
{
    double d = (d1 + d2);

    double ret = 0;

    if ((s1_a1 == s2_a1) && (s1_a2 == s2_a2))
	{
		// a match
        // Use SPM and approximate effect of delta
		double bf = spm.prob((s1_a1 == s1_a2), back[s1_a1], back[s1_a2]);

        ret += (1 - d) / bf;
	}

    // The smaller terms assume HW and approximate effect of delta
    ret += (s1_a1 == s2_a1) * 0.5 * d / back[s1_a1];
    ret += (s1_a1 == s2_a2) * 0.5 * d / back[s1_a1];

    ret += (s1_a2 == s2_a1) * 0.5 * d / back[s1_a2];
    ret += (s1_a2 == s2_a2) * 0.5 * d / back[s1_a2];

    ret += d * d;

    return ret;
}

static double // likelihood ratio, general method, calculating LR separately for each case
match(HMatrix const &prof1,
      HMatrix const &prof2,
      BackFreq const &hback)
{
	return dot(prof1, prof2 / hback);
}

BackFreq
calcBackForAAAATEST(PMF<Allele> const &back, SubPopModel const &spm)
{
    BackFreq hback(back, spm);
    Allele A = Allele(1);

    for (PMF< std::pair<Allele, Allele> >::iterator it = hback.m_pmf.begin();
         it != hback.m_pmf.end();
         ++it)    // for all i,j in hback
    {
        Allele i = it->first.first;
        Allele j = it->first.second;
    //            double p = it->second;

        double fi = back.val(i);
        double fj = back.val(j);
        double fst = spm.theta_bar;

        // P(ij|AA) = P(i|AA) * P(j|iAA) : *2 if i!=j
        int y1 = 2 * delta(i,A);
        int y2 = 2 * delta(j,A) + delta(j,i);
        double p = (2-delta(i,j)) * BBSF(y1, 2, fst, fi) * BBSF(y2, 3, fst, fj);

        hback.set(i, j, p);
    }

    return hback;
}

//
// exact match - general method
//
static double // likelihood ratio (= 0 if no match)
match(const Locus &locus,
	  const AlleleSet &s1,
	  const AlleleSet &s2,
	  SubPopModel spm,
	  PopulationData const &popdata,
	  const ProfileData &d1,
	  const ProfileData &d2,
	  const MatchType &match_type)
{
	if (! popdata.hasLocus(locus))
	{
		// no population data - ignore this locus
		return 1;
	}

    // use NRC4_4 for IDENT only, otherwise use B11
    if (spm.type == SubPopModel::NRC4_4 && match_type.m_rel_type != ident_t)
    {
        spm.type = SubPopModel::B11;
    }

	PMF<Allele> back = popdata(locus); // copy - may be modified
//	cout << "back = " << back << endl;
//    cout << "back.sum() = " << back.sum() << endl;

	if (!s1.checkBackground(back) || !s2.checkBackground(back))
	{
        error << startl << "unknown alleles: Locus ignored at " << locus_name[locus] << std::endl;

        // Ignore this locus. NB the better way to deal with unknown alleles is to
        // either add them at the rare allele frequency, or ignore (replace them with 'F').
        // This should be done by ProfileFilter when profiles are read in. At this stage we are
        // handling an error condition.

	    return 1;
	}

    if (SIMPLE_OPTIMIZATION)
    {
        // simple case optimization
        // Todo: extend to relative-vs-simple case
        if ((match_type.m_rel_type == ident_t) && (spm.type != SubPopModel::B11))
        {
            Allele s1_a1, s1_a2, s2_a1, s2_a2;

            if (d1.m_num_contributors == 1 &&
                d2.m_num_contributors == 1 &&
                s1.simple(s1_a1, s1_a2) &&
                s2.simple(s2_a1, s2_a2))
            {
                return match_simple(
                             s1_a1, s1_a2,
                             s2_a1, s2_a2,
                             back,
                             d1.m_error_rate, d2.m_error_rate,
                             spm);
            }
        }
    }

	// non-simple case: use  HMatrix
#ifdef SPARSE_OPTIMIZATION
	bool sparse = true;
#else
	bool sparse = false;
#endif

    // construct background frequencies consistent with sub-population model
    BackFreq hback(back, spm);

    // If we have a Subpopulation model (non-HW), use a specific function
    // Otherwise use the old HW functions which may be faster

    HMatrix h1, h2, hR;

    if (spm.type == SubPopModel::HW)
    {
//        PMF<Allele> brel = back; // background to be used when calculating relatives

        h1 = s1.getMatchHMatrix(back, d1, sparse); // cached
        h2 = s2.getMatchHMatrix(back, d2, sparse);

        if (match_type.m_rel_type == ident_t)
        {
            hR = h1;
        }
        else
        {
            hR = s1.getRelHMatrix(match_type, back, d1, sparse); // not yet cached
        }
    }
    else if (spm.type == SubPopModel::B11) // B11 method with correction
    {
        // B11: CONDITIONAL subpopulation model. Equivalent to NRC4_10 for Identity and non-probabilstic)
        // This model is HW *within the subpopulation*
        // but uses NRC4_4 background frequencies.

        h1 = s1.getMatchHMatrixNHW(hback, d1, sparse); // need NHW treatment of 'F' (cached)
        h2 = s2.getMatchHMatrixNHW(hback, d2, sparse);

        if (match_type.m_rel_type == ident_t)
        {
            hR = h1;
        }
        else
        {

            hR = relBN(match_type, h1, back, spm.theta_bar);
        }

    }
    else // other non-HW cases: NRC4_4 identity only (NRC4_10 is depreciated)
	{
        Assert2(spm.type == SubPopModel::NRC4_4 && match_type.m_rel_type == ident_t,
                "Can't handle Subpopulation Model");

	    // Construct HMatrices for SPM
        h1 = s1.getMatchHMatrixNHW(hback, d1, sparse);
        h2 = s2.getMatchHMatrixNHW(hback, d2, sparse);
        hR = h1;
	}

    if (d1.m_mutation_rate > 0 && match_type.m_rel_type != ident_t)
    {
        // Apply mutation to hR. NB it is correct that we apply the
        // mutation after calculating hR (not the other way around).
        // This is because of the way conditional probabilities are handled.

        double r = d1.m_mutation_rate * match_type.generations();
        hR = mutate_1s(hR, back, r);
    }

    double lr = match(hR, h2, hback);

    if (spm.type == SubPopModel::B11)
    {
        // Calculate correction factor and divide
        // (It does not matter whether we calculate C for s2 and multiply by s1 or vice-versa)
        HMatrix C = s2.getSpmcHMatrix(d2, back, spm); // cached
        double f = hdot(h1, C);
        lr = lr / f;
    }

    return lr;
}

// likelihood ratio of exact match between p1, p2
// i.e.
//     (likelihood of p1 from a individual with profile p2)
//     ______________________________________________________
//     (likelihood of p1 from a randomly selected individual)
//
// NB simple case: assume numerator = 1 if there is an exact match, 0 otherwise
//

// This implementation does the match one locus at a time and is therefore not suitable for CUDA profile/profile acceleration
// However this is the fastest general solution on the CPU
double
genmatch(const Profile &p1,
	     const Profile &p2,
	     const MatchType &match_type,
	     SubPopModel const &spm,
		 PopulationData const &popdata)
{
	double likelihood_ratio = 1;

	// If all loci **that are present in both profiles** match return LR, else 0

	for (LocusSet::const_iterator ip1 = p1.begin();
		ip1 != p1.end();
		++ip1)
	{
		Locus locus = ip1->first;

		if (locus == AMEL)
		{
			// ignore amelogenin
			continue;
		}

		LocusSet::const_iterator ip2 = p2.find(locus);

		// if second profile contains this locus
		if (ip2 != p2.end())
		{
			// see if alleles match
			double lhr;
			if ((lhr = match(locus, ip1->second, ip2->second, spm, popdata, p1.data(), p2.data(), match_type)) == 0)
			{
				if ( match_type.m_rel_type != ident_t && match_type.m_path1steps != 1 )
				{
//					warn << startl << "Exclusion (p=0) seen for match of type: " << match_type << endl;
//					warn << alignl << "(Exclusions should be possible only for IDENTITY or DEGREE(1, X) relationships)" << endl;
				}

				// The ability to short-circuit a match considerably speeds up
				// match as compared with sibmatch (but not when we set the error rate)

				return 0; // no match
			}
			else
			{
                info << startl << "locus = " << locus << " LR = " << lhr << endl;
				likelihood_ratio *= lhr;
			}
//			cout << "genmatch(): lhr = " << lhr << endl;
		}
		else
		{
			// missing locus - ignore
		}
	}

	return likelihood_ratio;

}

// likelihood ratio that p1, p2 come from siblings
// i.e.
//     (likelihood of p1 from a sibling of individual with p2)
//     ______________________________________________________
//     (likelihood of p1 from an unrelated individual)
//
// return the product of the LRs over all loci. If we find an exclusion
// (i.e. a locus at which the numerator is 0) we can return 0 immediately.
// NB there are in fact no exclusions in a sibling match.
//

double
match(const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm)
{
	return match(p1, p2, spm, popdata);
}

double
match(const Profile &p1, const Profile &p2, SubPopModel const &spm, PopulationData const &popdata)
{
	return genmatch(p1, p2, MatchType(ident_t), spm, popdata);
}

double
sibmatch(const Profile &p1, const Profile &p2, SubPopModel const &spm, PopulationData const &popdata)
{
	return genmatch(p1, p2, MatchType(sibling_t), spm, popdata);
}

double
sibmatch(const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm)
{
	return genmatch(p1, p2, MatchType(sibling_t), spm, popdata);
}

double
relmatch(const MatchType &match_type, const Profile &p1, const Profile &p2, SubPopModel const &spm, PopulationData const &popdata)
{
	return genmatch(p1, p2, match_type, spm, popdata);
}

double
relmatch(const MatchType &match_type, const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm)
{
	return genmatch(p1, p2, match_type, spm, popdata);
}

HMatrix // Element (p, q) is B(pq|ij)
condProb(BackFreq const &hback, // 4.1 frequencies
     Allele const &i,
     Allele const &j,
     double fst)
{
    HMatrix B; // B(pq|ij)

    for (PMF< std::pair<Allele, Allele> >::const_iterator it = hback.m_pmf.begin();
         it != hback.m_pmf.end();
         ++it)    // for all i,j in hback
    {
        Allele p = it->first.first;
        Allele q = it->first.second;
        double fp = hback(p); // HW frequency of p
        double fq = hback(q); // HW frequency of q

        double x = probBNpq_giv_ij(p, q, i, j, fp, fq, fst);

        B.set(p, q, x);
    }

    return B;
}

HMatrix
FSTCorrection(
        HMatrix const &prof,
        BackFreq const &hback, // 4.1 frequencies
        double fst)
{
    HMatrix C;

    for (PMF< std::pair<Allele, Allele> >::const_iterator it = hback.m_pmf.begin();
         it != hback.m_pmf.end();
         ++it)    // for all i,j in hback
    {
        Allele i = it->first.first;
        Allele j = it->first.second;
//        double p = it->second;

        HMatrix Bpq_ij = condProb(hback, i, j, fst); // B(pq|ij)

        double x = match(prof, Bpq_ij, hback);

        C.set(i, j, (i==j) ? x : 2*x); // *2 because an HMatrix
    }

    return C;
}

TEST(BBSF0)
{
    // Check BBSF against known formulae for P(aa|aa) and P(ab|ab)

    double fst = 0;
    SubPopModel nrc410(SubPopModel::NRC4_10, fst);

    double hom = nrc410.prob(true, 0.1, 0.1);
    double het = nrc410.prob(false, 0.1, 0.2);

    double hom1 = BBSF(2, 2, fst, 0.1) * BBSF(3, 3, fst, 0.1);
    CHECK_EQUAL(hom, hom1);

    double het1 = 2 * BBSF(1, 2, fst, 0.1) * BBSF(1, 3, fst, 0.2); // 2 for heterozygote
    double het2 = 2 * BBSF(1, 2, fst, 0.2) * BBSF(1, 3, fst, 0.1);
    CHECK_EQUAL(het, het1);
    CHECK_EQUAL(het, het2);
}

TEST(BBSF)
{
    // Check BBSF against known formulae for P(aa|aa) and P(ab|ab)

    double fst = 0.02;
    SubPopModel nrc410(SubPopModel::NRC4_10, fst);

    double hom = nrc410.prob(true, 0.1, 0.1);
    double het = nrc410.prob(false, 0.1, 0.2);

    double hom1 = BBSF(2, 2, fst, 0.1) * BBSF(3, 3, fst, 0.1);
    CHECK_CLOSE(hom, hom1, 1e-6);

    double het1 = 2 * BBSF(1, 2, fst, 0.1) * BBSF(1, 3, fst, 0.2); // 2 for heterozygote
    double het2 = 2 * BBSF(1, 2, fst, 0.2) * BBSF(1, 3, fst, 0.1);
    CHECK_CLOSE(het, het1, 1e-6);
    CHECK_CLOSE(het, het2, 1e-6);

    double het0 = 2 * 0.1 * 0.2 * (1 - fst);     // P(AB)
    double hom0 = 0.1 * 0.1 + 0.1*(1 - 0.1)*fst; // P(AA)

    CHECK_EQUAL(het0, 2 * BBSF(0, 1, fst, 0.1) * BBSF(0, 0, fst, 0.2));
    CHECK_EQUAL(hom0,     BBSF(1, 1, fst, 0.1) * BBSF(0, 0, fst, 0.1));
}

struct ProfilesFixture
{
	ProfilesFixture()
	: id1(Identifiler, "P01")
	, id2(Identifiler, "P02")
	, p1(id1)
	, p2(id2)
	{
  	        // population database for FGA:
	        PMF<Allele> fga;
	        fga[Allele(18)]    = 0.02649;
	        fga[Allele(19)]    = 0.05298;
	        fga[Allele(20)]    = 0.12748;
	        fga[Allele(21)]    = 0.18543;
	        fga[Allele(21, 2)] = 0.00497;

	        // population database for D7S820:
	        PMF<Allele> d7s820;
	        d7s820[Allele(7)]    = 0.01821;
	        d7s820[Allele(8)]    = 0.15066;
	        d7s820[Allele(8, 1)] = 0.00166;
	        d7s820[Allele(9)]    = 0.17715;
	        d7s820[Allele(10)]   = 0.24338;
	        d7s820[Allele(11)]   = 0.20695;
	        d7s820[Allele(12)]   = 0.16556;
	        d7s820[Allele(13)]   = 0.03477;
	        d7s820[Allele(14)]   = 0.00166;

            // population database for D18S51:
            PMF<Allele> d18s51;
            d18s51[Allele(10)] = 0.00828;
            d18s51[Allele(11)] = 0.01656;
            d18s51[Allele(12)] = 0.12748;
            d18s51[Allele(13)] = 0.13245;
            d18s51[Allele(14)] = 0.13742;
            d18s51[Allele(14, 2)] = 0.00166;
            d18s51[Allele(15)] = 0.15894;;
            d18s51[Allele(16)] = 0.13907;
            d18s51[Allele(17)] = 0.12583;
            d18s51[Allele(18)] = 0.07616;
            d18s51[Allele(19)] = 0.03808;
            d18s51[Allele(20)] = 0.02152;
            d18s51[Allele(21)] = 0.00828;
            d18s51[Allele(22)] = 0.00828;

            // population database for AMEL:
            PMF<Allele> amel;
            amel[Allele(Allele::X)] = 0.5;
            amel[Allele(Allele::Y)] = 0.5;

            PopulationData testpop;
	        testpop.setLocus(FGA, fga);
            testpop.setLocus(D7S820, d7s820);
            testpop.setLocus(D18S51, d18s51);
	        testpop.setSampleSize(302);

	        popSet(testpop);
	}

	ProfileData id1;
	ProfileData id2;
	Profile p1;
	Profile p2;
};

TEST_FIXTURE(ProfilesFixture, match_1)
{
	// test the simplest match algorithm
	// Both alleles as all loci must match (but ignore Amelogenin)

	// empty profiles match
	CHECK_EQUAL(1, match(p1, p2));

	// This matches too. Only mismatched alleles present in both profiles cause failure
	p1[D7S820] = AlleleSet(Allele(26, 1), Allele(26, 1));
	CHECK(match(p1, p2) > 0);
	CHECK(match(p2, p1) > 0);

	// two identical alleles at each locus - BINGO
	p2[D18S51] = AlleleSet(Allele(14, 2), Allele(14, 2));
	CHECK(match(p1, p2) > 0);

	// check (X, Y) = (Y, X)
	double lr;
	p1[FGA] = AlleleSet(20, 21);
	p2[FGA] = AlleleSet(21, 20);
	CHECK((lr = match(p1, p2)) > 0);

	// AMEL does not affect the lr
	p1[AMEL] = AlleleSet(Allele::X, Allele::Y);
	p2[AMEL] = AlleleSet(Allele::Y, Allele::X);
	CHECK_EQUAL(lr, match(p1, p2));

	// AMEL does not affect the lr
	p1[AMEL] = AlleleSet(Allele::X, Allele::X);
	p2[AMEL] = AlleleSet(Allele::Y, Allele::X);
	CHECK_EQUAL(lr, match(p1, p2));

	// ! overwriting FGA
	p1[FGA] = AlleleSet(Allele(20, 1), 21);
	p2[FGA] = AlleleSet(21, Allele(20, 1)); // generates "unknown allele" warning (FGA)
	if (populationData().sampleSize() > 0)
	{
		CHECK_EQUAL(1, match(p1, p2)); // ignored: LR = 1;
	}

	// check (A, B) != (B, C)
	p1[D7S820] = AlleleSet(Allele(8), Allele(10));
	p2[D7S820] = AlleleSet(Allele(10), Allele(8, 1));
	CHECK_EQUAL(0, match(p1, p2));
}

TEST_FIXTURE(ProfilesFixture, match_nomatch_1)
{
	// check (A, A) != (A, B)
	p1[D7S820] = AlleleSet(Allele(10), Allele(10));
	p2[D7S820] = AlleleSet(Allele(10), Allele(8, 1));
	CHECK_EQUAL(0, match(p1, p2));
}

TEST_FIXTURE(ProfilesFixture, match_nomatch_2)
{
	// check (A, A) != (B, A) order doesn't matter
	p1[D7S820] = AlleleSet(Allele(8), Allele(10));
	p2[D7S820] = AlleleSet(Allele(8, 1), Allele(10));
	CHECK_EQUAL(0, match(p1, p2));
}

enum TestLocus
{
  LOCUS1 = 0
, LOCUS2
, LOCUS3
, LOCUS4
, LOCUS5
, LOCUS6
, LOCUS7
, LOCUS8
, LOCUS9
, LOCUS10
, LOCUS11
, LOCUS12
, LOCUS13
, AMEL_HERE_IGNORED
, LOCUS14
, LOCUS15
, LOCUS16
, test_size
};

enum TestAlleles
{
  F = 0 // unknown
, A = 10
, B
, C
, D
, Z
, U
, V
, W
};

static PMF<Allele>::POD LOCUS1_FREQ[] =
{
  { Allele(A), 0.1 } // Allele() optional
, { B, 0.2 }
, { C, 0.3 }
, { D, 0.4 }
};

struct SpreadsheetFixture
{
	SpreadsheetFixture()
	: id1(test, "P01")
	, id2(test, "P02")
    , id3(test, "P03")
	, p1(id1)
	, p2(id2)
    , p3(id3)
	, bg_locus1(LOCUS1_FREQ, sizeof(LOCUS1_FREQ)/sizeof(PMF<Allele>::POD))
	{
	    for (TestLocus i =  LOCUS1; i < test_size; i = TestLocus(i + 1))
	    {
	        popdata.setLocus(i, bg_locus1);
	    }
	}

	ProfileData id1;
	ProfileData id2;
    ProfileData id3;
	Profile p1;
	Profile p2;
    Profile p3;
	PMF<Allele> bg_locus1;
	PopulationData popdata;
};

static PMF<Allele>::POD LOCUS2_FREQ[] =
{
  { Allele(A), 0.2 } // Allele() optional
, { B, 0.16 }
, { C, 0.04 }
, { D, 0.6 }
};

struct SpreadsheetFixture2
{
	SpreadsheetFixture2()
	: id1(test, "P01")
	, id2(test, "P02")
	, p1(id1)
	, p2(id2)
	, bg_locus2(LOCUS2_FREQ, 4)
	{
		popdata.setLocus(LOCUS2, bg_locus2); //
	}

	ProfileData id1;
	ProfileData id2;
	Profile p1;
	Profile p2;
	PMF<Allele> bg_locus2;
	PopulationData popdata;
};

TEST_FIXTURE(SpreadsheetFixture, condProb)
{
    double fst = 0.02;
    SubPopModel nrc410(SubPopModel::NRC4_10, fst);

    double aaaa = nrc410.prob(true, 0.1, 0.1);
    double abab = nrc410.prob(false, 0.1, 0.2);

    BackFreq bf(bg_locus1, SubPopModel());
    CHECK_CLOSE(0.1, bf(A), 1e-6);
    CHECK_CLOSE(0.2, bf(B), 1e-6);
    CHECK_CLOSE(0.3, bf(C), 1e-6);
    CHECK_CLOSE(0.4, bf(D), 1e-6);

    CHECK_CLOSE(0.01, bf(make_pair(A,A)), 1e-6);
    CHECK_CLOSE(0.04, bf(make_pair(B,B)), 1e-6);
    CHECK_CLOSE(0.04, bf(make_pair(A,B)), 1e-6);

    HMatrix Bpq_ij = condProb(bf, A, A, fst); // Element (p, q) is B(pq|ij)
    CHECK_CLOSE(aaaa, Bpq_ij(A, A), 1e-6);       // P(AA|AA)
    double abaa = 2 * BBSF(2, 2, fst, 0.1) * BBSF(0, 3, fst, 0.2);
    CHECK_CLOSE(abaa, Bpq_ij(A, B), 1e-6);       // P(AB|AA)

    Bpq_ij = condProb(bf, A, B, fst);
    CHECK_CLOSE(abab, Bpq_ij(A, B), 1e-6);       // P(AB|AB)
    double aaab = BBSF(1, 2, fst, 0.1) * BBSF(2, 3, fst, 0.1);
    CHECK_CLOSE(aaab, Bpq_ij(A, A), 1e-6);       // P(AA|AB) = P(A|AB)*P(A|AAB)

    // background according to dirichlet (or NRC II 4.4)
#if 1
    double ab = 2 * 0.1 * 0.2 * (1 - fst);     // P(AB)
    double aa = 0.1 * 0.1 + 0.1*(1 - 0.1)*fst; // P(AA)
#else
    // HW
    double ab = 2 * 0.1 * 0.2;     // P(AB)
    double aa = 0.1 * 0.1;         // P(AA)
#endif

    // check we obey Bayes' theorem
    CHECK_EQUAL(abaa * aa, aaab * ab);
}

TEST_FIXTURE(SpreadsheetFixture, condProb2)
{

}

TEST_FIXTURE(SpreadsheetFixture, AliceAndBob)
{
    // This test is used as an example in the Algorithms Document
    p1[LOCUS1] = AlleleSet(A, B); // Alice
    p2[LOCUS1] = AlleleSet(A, F); // Bob

    double lr = 3; // calculated by hand
    CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), 1e-6);
    CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), 1e-6);
}

TEST_FIXTURE(SpreadsheetFixture2, pmatch_test_1)
{
	p1[LOCUS2] = AlleleSet(A, B);
	p2[LOCUS2] = AlleleSet(A, B);

	double lr = 5.5625; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, match_2)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	double lh = 0.1*0.1;
	CHECK_CLOSE(1/lh, match(p1, p2, popdata), 1e-3);
}

TEST_FIXTURE(SpreadsheetFixture, match_3alleles)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	p1[LOCUS2] = AlleleSet(A, B);
	p2[LOCUS2] = AlleleSet(A, B);

	p1[LOCUS3] = AlleleSet(C, C);
	p2[LOCUS3] = AlleleSet(C, C);

	p1[LOCUS4] = AlleleSet(C, D);
	p2[LOCUS4] = AlleleSet(C, D);

	double lh = 0.1*0.1 * 2*0.1*0.2 * 0.3*0.3 * 2*0.3*0.4;
	CHECK_CLOSE(1/lh, match(p1, p2, popdata), 1);
}

#define CHECK_NEAR(expected, actual, tolerance) \
	CHECK_CLOSE( (expected), (actual), (expected * tolerance) )

TEST_FIXTURE(SpreadsheetFixture, match_3alleles_nomatch)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	p1[LOCUS2] = AlleleSet(A, B);
	p2[LOCUS2] = AlleleSet(A, B);

	p1[LOCUS3] = AlleleSet(C, C); // homozygote
	p2[LOCUS3] = AlleleSet(C, D); // not a match at second position

	p1[LOCUS4] = AlleleSet(C, D);
	p2[LOCUS4] = AlleleSet(C, D);

	CHECK_EQUAL(0, match(p1, p2, popdata));

	//
	// test delta error term
	//

	// in the general algorithm, set delta non-zero and there
	// should be a non-zero likelihood of a match

	// likelihood for the three matching positions (with delta = 0)
	double lh = 0.1*0.1 * 2*0.1*0.2 * 2*0.3*0.4;

	double delta = 1e-3;
    float  near  = 1e-2; // the calculated values are approximate

	p1.setErrorRate(delta);
	p2.setErrorRate(0);

//    CHECK_NEAR((1/lh)*(delta1_exact(C, C, C, D, .3, .4, delta)), match(p2, p1, popdata), near);

    float lr_exact =
            delta2_exact(A, A, A, A, .1, .1, delta, 0) *
            delta2_exact(A, B, A, B, .1, .2, delta, 0) *
            delta2_exact(C, C, C, D, .3, .3, delta, 0) *
            delta2_exact(C, D, C, D, .3, .4, delta, 0);

    CHECK_NEAR(lr_exact, match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(C, C, C, D, .3, .3, delta, 0)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, C, C, D, .3, .3, delta, 0)), match(p1, p2, popdata), near);

    p1.setErrorRate(0);
    p2.setErrorRate(delta);

    CHECK_NEAR((1/lh)*(delta1_exact(C, D, C, C, .3, .3, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(C, C, C, D, .3, .3, 0, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, C, C, D, .3, .3, 0, delta)), match(p1, p2, popdata), near);

    p1.setErrorRate(delta);
    p2.setErrorRate(delta/2);

    CHECK_NEAR((1/lh)*(delta2_exact(C, C, C, D, .3, .3, delta, delta/2)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, C, C, D, .3, .3, delta, delta/2)), match(p1, p2, popdata), near);

//	CHECK_GREATER(1/(lh*delta), match(p1, p2, popdata));
//  CHECK_LESS(1/lh, match(p1, p2, popdata)); // match less likely than just seeing the three matching loci

    // NB if *both* are heterozygotes, we need a factor of .5 ???

    p1[LOCUS3] = AlleleSet(C, D); // heterozygote
    p2[LOCUS3] = AlleleSet(B, D); // not a match at first position
//    p2[LOCUS3] = AlleleSet(A, D); // not a match at first position
//    p2[LOCUS3] = AlleleSet(D, D); // not a match at first position

    p1.setErrorRate(delta);
    p2.setErrorRate(0);

    CHECK_NEAR((1/lh)*(delta1_exact(C, D, B, D, .2, .4, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(C, D, B, D, .3, .4, delta, 0)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, D, B, D, .2, .4, delta, 0)), match(p1, p2, popdata), near);

    p1.setErrorRate(0);
    p2.setErrorRate(delta);

    CHECK_NEAR((1/lh)*(delta1_exact(B, D, C, D, .3, .4, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(C, D, B, D, .3, .4, 0, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, D, B, D, .2, .4, 0, delta)), match(p1, p2, popdata), near);

    p1.setErrorRate(delta);
    p2.setErrorRate(delta/2);

    CHECK_NEAR((1/lh)*(delta2_exact(C, D, B, D, .3, .4, delta, delta/2)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(C, D, B, D, .2, .4, delta, delta/2)), match(p1, p2, popdata), near);

//    p1[LOCUS3] = AlleleSet(C, C);
    p1[LOCUS3] = AlleleSet(B, C);

    p2[LOCUS3] = AlleleSet(A, D); // not a match at both positions
//    p2[LOCUS3] = AlleleSet(D, D); // not a match at both positions

    p1.setErrorRate(delta);
    p2.setErrorRate(0);

    CHECK_NEAR((1/lh)*(delta1_exact(B, C, A, D, .1, .4, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(B, C, A, D, .2, .3, delta, 0)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(B, C, A, D, .1, .4, delta, 0)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta*delta), match(p1, p2, popdata), near);

    p1.setErrorRate(0);
    p2.setErrorRate(delta);

    CHECK_NEAR((1/lh)*(delta1_exact(A, D, B, C, .2, .3, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_exact(B, C, A, D, .2, .3, 0, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(B, C, A, D, .1, .4, 0, delta)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta*delta), match(p1, p2, popdata), near);

    p1.setErrorRate(delta);
    p2.setErrorRate(delta/2);

    CHECK_NEAR((1/lh)*(delta2_exact(B, C, A, D, .2, .3, delta, delta/2)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh)*(delta2_approx(B, C, A, D, .1, .4, delta, delta/2)), match(p1, p2, popdata), near);
    CHECK_NEAR((1/lh) * (delta + delta/2) * (delta + delta/2), match(p1, p2, popdata), near);

}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case1)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	double lr = 30.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case2_a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, B);

	double lr = 2.75; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case2_b)
{
	p1[LOCUS1] = AlleleSet(C, D);
	p2[LOCUS1] = AlleleSet(D, D);

	double lr = 0.875; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_3loci)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	p1[LOCUS2] = AlleleSet(A, A);
	p2[LOCUS2] = AlleleSet(A, B);

	p1[LOCUS3] = AlleleSet(C, D);
	p2[LOCUS3] = AlleleSet(D, D);

	double lr = 30.25 * 2.75 * 0.875; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case3_a)
{
	p1[LOCUS1] = AlleleSet(B, B);
	p2[LOCUS1] = AlleleSet(C, C);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case3_b)
{
	p1[LOCUS1] = AlleleSet(C, C);
	p2[LOCUS1] = AlleleSet(B, B);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case4_a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(B, C);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case4_b)
{
	p1[LOCUS1] = AlleleSet(C, D);
	p2[LOCUS1] = AlleleSet(B, B);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case5_a)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, B);

	double lr = 8.375; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case5_b)
{
	p1[LOCUS1] = AlleleSet(C, D);
	p2[LOCUS1] = AlleleSet(C, D);

	double lr = 2.020833; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case6_a)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, C);

	double lr = 1.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case6_b)
{
	p1[LOCUS1] = AlleleSet(A, D);
	p2[LOCUS1] = AlleleSet(B, D);

	double lr = 0.5625; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case6_c)
{
	p1[LOCUS1] = AlleleSet(A, D);
	p2[LOCUS1] = AlleleSet(B, D);

	double lr = 0.5625; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case7_a)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(C, D);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, sibmatch_case7_b)
{
	p1[LOCUS1] = AlleleSet(A, D);
	p2[LOCUS1] = AlleleSet(B, C);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

// test unknown allele
TEST_FIXTURE(SpreadsheetFixture, sibmatch_general_unknowns)
{
	p1[LOCUS1] = AlleleSet(A, A);

	/////
    p2[LOCUS1] = AlleleSet(A, A);
    CHECK_CLOSE(30.25, sibmatch(p1, p2, popdata), tol);
    p2[LOCUS1] = AlleleSet(A, B);
    CHECK_CLOSE(2.75, sibmatch(p1, p2, popdata), tol);
    p2[LOCUS1] = AlleleSet(A, C);
    CHECK_CLOSE(2.75, sibmatch(p1, p2, popdata), tol);
    p2[LOCUS1] = AlleleSet(A, D);
    CHECK_CLOSE(2.75, sibmatch(p1, p2, popdata), tol);
	/////

	p2[LOCUS1] = AlleleSet(A, F); // F = unknown

	// should give a weighted sum of (AA AA), (AA AB), (AA AC), (AA AD)

	double lr = 0.1 * 30.25 + 0.2 * 2.75 +  0.3 * 2.75 + 0.4 * 2.75; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

// test PMF of alleles
TEST_FIXTURE(SpreadsheetFixture, sibmatch_general_alternatives)
{
	p1[LOCUS1] = AlleleSet(A, A);

	PMF<Allele> Apmf, AorBpmf;
	Apmf.insert(make_pair(A, 1.0));
	AorBpmf.insert(make_pair(A, 0.9));
	AorBpmf.insert(make_pair(B, 0.1));

	vector< PMF<Allele> > av;
	av.push_back(Apmf);        // A
	av.push_back(AorBpmf);     // 0.9A + 0.1B

	p2[LOCUS1] = AlleleSet(av);

	// should give a weighted sum of (AA AA), (AA AB)

	double lr = 0.9 * 30.25 + 0.1 * 2.75; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mixtures_1)
{
	p1.setNumContributors(2);
	p1[LOCUS1] = AlleleSet(A, B, C, D);
	p2[LOCUS1] = AlleleSet(A, B);

	// match
	double lr = 4.166667; // from Sums spreadsheet
	CHECK_CLOSE(lr, match(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, match(p2, p1, popdata), tol);

	// sibling match
	lr = 2.229167; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mixtures_2)
{
	p1.setNumContributors(2);
	p1[LOCUS1] = AlleleSet(A, A, B, B);
	p2[LOCUS1] = AlleleSet(A, B);

	// match
	double lr = 16.666667; // from Sums spreadsheet
	CHECK_CLOSE(lr, match(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, match(p2, p1, popdata), tol);

	// sibling match
	lr = 6.291667; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mixtures_3)
{
	p1.setNumContributors(2);
	p1[LOCUS1] = AlleleSet(A, B, A, B); // order should not matter
	p2[LOCUS1] = AlleleSet(A, B);

	// match
	double lr = 16.666667; // from Sums spreadsheet
	CHECK_CLOSE(lr, match(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, match(p2, p1, popdata), tol);

	// sibling match
	lr = 6.291667; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mixtures_4)
{
	p1.setNumContributors(2);
	p1[LOCUS1] = AlleleSet(A, A, A, B);
	p2[LOCUS1] = AlleleSet(A, B);

	// match
	double lr = 12.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, match(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, match(p2, p1, popdata), tol);

	// sibling match
	lr = 5.5625; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mixtures_5)
{
	p1.setNumContributors(2);
	p1[LOCUS1] = AlleleSet(A, B, B, B);
	p2[LOCUS1] = AlleleSet(A, B);

	// match
	double lr = 12.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, match(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, match(p2, p1, popdata), tol);

	// sibling match
	lr = 4.9375; // from Sums spreadsheet
	CHECK_CLOSE(lr, sibmatch(p1, p2, popdata), tol);
	CHECK_CLOSE(lr, sibmatch(p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_degree) // degree functions
{
	HMatrix h1 = AlleleSet(A, A).getMatchHMatrix(popdata(LOCUS1));
	HMatrix h1a = degree1(h1, popdata(LOCUS1));
	HMatrix h1b = degree(1, h1, popdata(LOCUS1));

	CHECK_EQUAL(h1a, h1b);

	HMatrix h2 = AlleleSet(A, B).getMatchHMatrix(popdata(LOCUS1));
	HMatrix h2a = degree1(h2, popdata(LOCUS1));
	HMatrix h2b = degree(1, h2, popdata(LOCUS1));

	CHECK_EQUAL(h2a, h2b);
}

// degree (1,0) = parent or child
TEST_FIXTURE(SpreadsheetFixture, rel_10a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	MatchType mtype(degree_1_t);

	double lr = 10; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_10b)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, C);

	MatchType mtype(degree_1_t);

	double lr = 2.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_10c)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(D, D);

	MatchType mtype(degree_1_t);

	double lr = 0; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

// degree (2,0) = grandparent/grandchild/uncle/nephew
TEST_FIXTURE(SpreadsheetFixture, rel_20a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	MatchType mtype(degree_2_t);

	double lr = 5.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_20b)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, C);

	MatchType mtype(degree_2_t);

	double lr = 1.75; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_20c)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(D, D);

	MatchType mtype(degree_2_t);

	double lr = 0.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

// degree (1, 2) = simultaneous child and grandchild
TEST_FIXTURE(SpreadsheetFixture, rel_21a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	MatchType mtype(degree_pq_t, 1, 2);

	double lr = 55; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);

	// check Balding eqn 7.14 for LRs
	double Ri = match(p1, p2, popdata);
	double Rp = relmatch(MatchType(degree_1_t), p1, p2, popdata);
	CHECK_CLOSE(lr, (Ri + Rp)/2, tol);

// need to extend these functions to bilinear relationships
//    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::B11, 0)), tol);
//    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::B11, 0)), tol);
//
//    CHECK_CLOSE(lr, baldingNichols(mtype, p2, p3, popdata, 0), tol);
//    CHECK_CLOSE(lr, baldingNichols(mtype, p2, p3, popdata, 0), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_21b)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, C);

	MatchType mtype(degree_pq_t, 1, 2);

	double lr = 1.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);

	// check Balding eqn 7.14 for LRs
    double Ri = match(p1, p2, popdata);
    double Rp = relmatch(MatchType(degree_1_t), p1, p2, popdata);
    CHECK_CLOSE(lr, (Ri + Rp)/2, tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_21c)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(D, D);

	MatchType mtype(degree_pq_t, 1, 2);

	double lr = 0; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);

	// check Balding eqn 7.14 for LRs
    double Ri = match(p1, p2, popdata);
    double Rp = relmatch(MatchType(degree_1_t), p1, p2, popdata);
    CHECK_CLOSE(lr, (Ri + Rp)/2, tol);

}

TEST_FIXTURE(SpreadsheetFixture, rel_21d) // This is an asymmetric relationship
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, B);

	MatchType mtype(degree_pq_t, 1, 2);

	double lr = 2.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);

	MatchType mtype_gen(gen_t, 0.5, 0.5, 0.25, 0.25);
    CHECK_CLOSE(lr, relmatch(mtype_gen, p1, p2, popdata), tol);

    double lr_inv = 15; // inverse relationship (by hand)
    CHECK_CLOSE(lr_inv, relmatch(mtype, p2, p1, popdata), tol);
    CHECK_CLOSE(lr_inv, relmatch(mtype_gen, p2, p1, popdata), tol);

	MatchType mtype_inv(inv_t, 0.5, 0.5, 0.25, 0.25);
	CHECK_CLOSE(lr, relmatch(mtype_inv, p2, p1, popdata), tol);
    CHECK_CLOSE(lr_inv, relmatch(mtype_inv, p1, p2, popdata), tol);

    // check Balding eqn 7.14 for LRs (It does not work for bilinear relationships!)
    double Ri = match(p1, p2, popdata);
    double Rp = relmatch(MatchType(degree_1_t), p1, p2, popdata);
    CHECK_CLOSE(lr, (Ri + Rp)/2, tol);

    Ri = match(p2, p1, popdata);
    Rp = relmatch(MatchType(degree_1_t), p2, p1, popdata);
//    CHECK_CLOSE(15, (Ri + Rp)/2, tol); // inverse relationship is 15. This fails.
}

TEST_FIXTURE(SpreadsheetFixture, rel_21e)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(C, D);

	MatchType mtype(degree_pq_t, 1, 2);

	double lr = 0; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);

	// check Balding eqn 7.14 for LRs
    double Ri = match(p1, p2, popdata);
    double Rp = relmatch(MatchType(degree_1_t), p1, p2, popdata);
    CHECK_CLOSE(lr, (Ri + Rp)/2, tol);
}

// degree (2, 2) = grandchild by a son and daughter
TEST_FIXTURE(SpreadsheetFixture, rel_22a)
{
	p1[LOCUS1] = AlleleSet(A, A);
	p2[LOCUS1] = AlleleSet(A, A);

	MatchType mtype(degree_pq_t, 2, 2);

	double lr = 30.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_22b)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(A, C);

	MatchType mtype(degree_pq_t, 2, 2);

	double lr = 1.5; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_22c)
{
	p1[LOCUS1] = AlleleSet(A, B);
	p2[LOCUS1] = AlleleSet(D, D);

	MatchType mtype(degree_pq_t, 2, 2);

	double lr = 0.25; // from Sums spreadsheet
	CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata), tol);
	CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol);
}

// Check that MatchType(sibling_t) is equivalent to MatchType(gen_t, 0.5, 0, 0, 0.5)
void
testSibling(Profile const &p1, Profile const &p2, PopulationData const &popdata)
{
    MatchType mtype(sibling_t);
    double lr = relmatch(mtype, p1, p2, popdata);

    MatchType gtype(gen_t, 0.5, 0, 0, 0.5); // sibling ?

    CHECK_CLOSE(lr, relmatch(gtype, p1, p2, popdata), tol);
    CHECK_CLOSE(lr, relmatch(gtype, p2, p1, popdata), tol);
}

TEST_FIXTURE(SpreadsheetFixture, rel_sib)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, A);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, B);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(B, B);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(B, C);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, B);
    p2[LOCUS1] = AlleleSet(A, B);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, B);
    p2[LOCUS1] = AlleleSet(A, C);
    testSibling(p1, p2, popdata);

    p1[LOCUS1] = AlleleSet(A, B);
    p2[LOCUS1] = AlleleSet(C, D);
    testSibling(p1, p2, popdata);

}

TEST_FIXTURE(SpreadsheetFixture, rel_24a) // asymmetric ???
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, B);

    MatchType mtype(degree_pq_t, 2, 4);
    double lr = relmatch(mtype, p1, p2, popdata);
//    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol); // asymmetric

    MatchType mtype_gen(gen_t, 0.25, 0.25, 0.0625, 0.0625);
    CHECK_CLOSE(lr, relmatch(mtype_gen, p1, p2, popdata), tol);
//    CHECK_CLOSE(lr, relmatch(mtype_gen, p2, p1, popdata), tol); // asymmetric

    MatchType mtype_inv(inv_t, 0.25, 0.25, 0.0625, 0.0625);
    CHECK_CLOSE(lr, relmatch(mtype_inv, p2, p1, popdata), tol);
//    CHECK_CLOSE(lr, relmatch(mtype_inv, p1, p2, popdata), tol); // asymmetric
}

TEST_FIXTURE(SpreadsheetFixture, rel_pga) // one shared parent, one shared grandparent
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, B);

    MatchType mtype(gen_t, 0.125, 0, 0, 0.5); // NB NOT a D24
    double lr = relmatch(mtype, p1, p2, popdata);

    MatchType mtype24(degree_pq_t, 2, 4);
//    double lr_true = relmatch(mtype24, p2, p1, popdata); // AB/AC 1.6875
    double lr_true = relmatch(mtype24, p1, p2, popdata); // AB/AC 1.6875
//    CHECK_EQUAL(0, lr_true);

//    double lr_true = 2.9375; // no independent calculation!
    CHECK_CLOSE(lr_true, lr, tol);

    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata), tol); // symmetric

    MatchType mtype_inv(inv_t, 0.125, 0, 0, 0.5);
    CHECK_CLOSE(lr, relmatch(mtype_inv, p2, p1, popdata), tol); // symmetric
    CHECK_CLOSE(lr, relmatch(mtype_inv, p1, p2, popdata), tol); // symmetric
}

// Subpopulation models
//const double theta = 0.01;

TEST_FIXTURE(SpreadsheetFixture, spm_nrc44_ident)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, A);

    double theta = 0;
    MatchType mtype(ident_t);

    double lr = 1/(0.1*0.1); // HW case
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::NRC4_4, theta)), tol);
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::NRC4_4, theta)), tol);
}

// Direct calculation of LR by BN method for testing

void getIBD(MatchType const &mtype, double &k0, double &k1, double &k2)
{
    switch(mtype.m_rel_type)
    {
    case ident_t:
    {
        k0 = 0; k1 = 0; k2 = 1;
        break;
    }
    case sibling_t:
    {
        k0 = 0.25; k1 = 0.5; k2 = 0.25;
        break;
    }
    case degree_1_t:
    {
        k0 = 0; k1 = 1; k2 = 0;
        break;
    }
    case degree_2_t:
    {
        k0 = 0.5; k1 = 0.5; k2 = 0;
        break;
    }
    case degree_pq_t:
    {
        int n = mtype.m_path1steps, m = mtype.m_path2steps;
        double np = pow(2, 1-n), mp = pow(2, 1-m);
        k2 = np * mp;
        k0 = (1-np)*(1-mp);
        k1 = 1 - k0 - k2;
        break;
    }
    case gen_t:
    case inv_t:
    default:
        Assert2(false, "Unknown match type in getIBD()");
    }

    Assert(k0+k1+k2 == 1);
}

// LR for identity, Balding-Nichols model
double RiBN(
        MatchType const &mtype,
        Allele const &p,
        Allele const &q,
        Allele const &i,
        Allele const &j,
        PMF<Allele> const &back,
        double theta)
{
    if ( (p!=i) || (q!=j) ) return 0; // must be the same

    SubPopModel spm(SubPopModel::NRC4_10, theta);
    return 1/spm.prob(p==q, back.val(p), back.val(q));
}

//
// LR for paternity, Balding-Nichols model
// From Balding p121 (the other way up)
//

double RpBN_AA_AA(double pa, double theta)
{
    return (1 + 2*theta) / (3*theta + (1 - theta)*pa);
}

double RpBN_AA_AB(double pa, double theta)
{
    return (1 + 2*theta) / (2*(2*theta + (1 - theta)*pa));
}

double RpBN_AB_AB(double pa, double pb, double theta)
{
    return (1 + 2*theta) / (4*(theta + (1 - theta)*pa) *
                              (theta + (1 - theta)*pb) /
                              (2*theta + (1 - theta)*(pa+pb)) );
}

double RpBN_AB_AC(double pa, double theta)
{
    return (1 + 2*theta) / (4*(theta + (1 - theta)*pa) );
}

// The D1 relationship, Balding-Nichols method
double RpBN(
        MatchType const &mtype,
        Allele const &p,
        Allele const &q,
        Allele const &i,
        Allele const &j,
        PMF<Allele> const &back,
        double theta)
{
    if ( (p!=i) && (p!=j) && (q!=i) && (q!=j)) return 0; // 0 unless share one allele

    if (p==q) // AA
    {
        if (i==j)
        {
            return RpBN_AA_AA(back.val(p), theta);
        }
        else  // AA/AB
        {
            return RpBN_AA_AB(back.val(p), theta);
        }
    }
    else      // AB
    {
        if (i==j) // AA/AB
        {
            return RpBN_AA_AB(back.val(i), theta);
        }
        else
        {
            if (p==i && q==j) // AB/AB
            {
                return RpBN_AB_AB(back.val(p), back.val(q), theta);
            }
            else      // AB/AC
            {
                if (p==i || p==j)
                {
                    return RpBN_AB_AC(back.val(p), theta);
                }
                else
                {
                    return RpBN_AB_AC(back.val(q), theta);
                }
            }
        }
    }
}

// Direct BN formula for test purposes. Return LR.
double fBN(
        MatchType const &mtype,
        Allele const &p,
        Allele const &q,
        Allele const &i,
        Allele const &j,
        PMF<Allele> const &back,
        double theta)
{
    double k0=0, k1=0, k2=0; // ibd coefficients
    getIBD(mtype, k0, k1, k2);

    double Ri = RiBN(mtype, p, q, i, j, back, theta);
    double Rp = RpBN(mtype, p, q, i, j, back, theta);

    return k0 + k1*Rp + k2*Ri; // Balding p126 Eqn. 7.14 (Rp, Ri the other way up)
}

// Direct BN formulae for test purposes
// NB the summation over probability matrices is not exact in general (it gives
// a different answer from the B11 correction)
// - but it is exact in the non-probabilistic case
double baldingNichols(
        MatchType const &mtype,
        Profile const &p1,
        Profile const &p2,
        PopulationData const &popdata,
        double theta)
{
    double likelihood_ratio = 1;

    // If all loci **that are present in both profiles** match return LR, else 0

    for (LocusSet::const_iterator ip1 = p1.begin();
        ip1 != p1.end();
        ++ip1)
    {
        Locus locus = ip1->first;

        if (locus == AMEL) continue;

        LocusSet::const_iterator ip2 = p2.find(locus);

        // if second profile contains this locus
        if (ip2 != p2.end())
        {
            PMF<Allele> back = popdata(locus); // copy - may be modified
            HMatrix h1 = ip1->second.getMatchHMatrix(back, p1.data(), true /*sparse*/);  // cached
            HMatrix h2 = ip2->second.getMatchHMatrix(back, p2.data(), true /*sparse*/);

            // sum over both HMatrices
            double lr = 0;
            for (PMF< std::pair<Allele, Allele> >::const_iterator it1 = h1.m_pmf.begin();
                 it1 != h1.m_pmf.end();
                 ++it1)
            {
                Allele p = it1->first.first;
                Allele q = it1->first.second;
                double w1 = it1->second;

                for (PMF< std::pair<Allele, Allele> >::const_iterator it2 = h2.m_pmf.begin();
                     it2 != h2.m_pmf.end();
                     ++it2)
                {
                    Allele i = it2->first.first;
                    Allele j = it2->first.second;
                    double w2 = it2->second;

                    lr = w1 * w2 * fBN(mtype, p, q, i, j, back, theta);
                }
            }

            if (lr == 0)
            {
                return 0; // no match
            }
            else
            {
                likelihood_ratio *= lr;
            }
        }
    }

    return likelihood_ratio;
}

string
desc(int a, int b, int c, int d)
{
    string ret = "(";
    ret = ret + char('A' - 1 + a);
    ret = ret + char('A' - 1 + b);
    ret = ret + '/';
    ret = ret + char('A' - 1 + c);
    ret = ret + char('A' - 1 + d);
    ret = ret + ')';
    return ret;
}

TEST_FIXTURE(SpreadsheetFixture, BN)
{
    // Check B11 correction factor against Balding-Nichols method

    TestAlleles p1_1 = A;

    int n = 0;
    for (TestAlleles p1_2 = A; p1_2 <= D; p1_2 = (TestAlleles)(p1_2 + 1) ) {
    for (TestAlleles p2_1 = A; p2_1 <= D; p2_1 = (TestAlleles)(p2_1 + 1) ) {
    for (TestAlleles p2_2 = p2_1; p2_2 <= D; p2_2 = (TestAlleles)(p2_2 + 1) )
    {
        ++n;

        p1[LOCUS1] = AlleleSet(p1_1, p1_2);
        p2[LOCUS1] = AlleleSet(p2_1, p2_2);

//        string descr = desc(p1_1, p1_2, p2_1, p2_2);

        double theta_val[] = { 0, 0.01, 0.02, 0.03, 0.05 };
        const int ntheta = 5;

        RelType type[] = { ident_t, degree_1_t, sibling_t, degree_2_t, degree_pq_t };
        const int nrels = 5;

        for (int i=0; i<ntheta; ++i)
        {
            double theta = theta_val[i];

            for (int j=0; j<nrels; ++j)
            {
                MatchType mtype_num(type[j], 3, MatchType::INF); // last two args used in D3 case only
                double tol = 1e-5;

//                cout << descr << " : type = " << setw(5) << mtype_num << " theta = " << theta << endl;

                // Using direct Balding-Nichols calculation
                double lr_bn = baldingNichols(mtype_num, p1, p2, popdata, theta);
                double lr_bn2 = baldingNichols(mtype_num, p2, p1, popdata, theta);
                CHECK_CLOSE(lr_bn, lr_bn2, tol);

                // Using B11 correction factor
                double lr = relmatch(mtype_num, p1, p2, popdata, SubPopModel(SubPopModel::B11, theta));
                double lr2 = relmatch(mtype_num, p2, p1, popdata, SubPopModel(SubPopModel::B11, theta));
                CHECK_CLOSE(lr, lr2, tol);

                CHECK_CLOSE(lr_bn, lr, tol);
            }
        }
    }}}

    CHECK_EQUAL(40, n); // 40 test cases
}

TEST_FIXTURE(SpreadsheetFixture, B11_profile)
{
    // Check B11 correction factor against Balding-Nichols method for an entire profile

    p1[LOCUS1]  = AlleleSet(10, 10); p2[LOCUS1]  = AlleleSet(10, 10); p3[LOCUS1]  = AlleleSet(10, 10); // D8
    p1[LOCUS2]  = AlleleSet(10, 11); p2[LOCUS2]  = AlleleSet(10, 10); p3[LOCUS2]  = AlleleSet(10, 10); // D2
    p1[LOCUS3]  = AlleleSet(10, 12); p2[LOCUS3]  = AlleleSet(10, 12); p3[LOCUS3]  = AlleleSet(10, 12); // D7
    p1[LOCUS4]  = AlleleSet(10, 13); p2[LOCUS4]  = AlleleSet(11, 13); p3[LOCUS4]  = AlleleSet(11, 13); // CS
    p1[LOCUS5]  = AlleleSet(11, 10); p2[LOCUS5]  = AlleleSet(12, 13); p3[LOCUS5]  = AlleleSet(11, 11); // P3
    p1[LOCUS6]  = AlleleSet(11, 11); p2[LOCUS6]  = AlleleSet(11, 11); p3[LOCUS6]  = AlleleSet(11, 11); // TH
    p1[LOCUS7]  = AlleleSet(11, 12); p2[LOCUS7]  = AlleleSet(11, 11); p3[LOCUS7]  = AlleleSet(11, 11); // D13
    p1[LOCUS8]  = AlleleSet(11, 13); p2[LOCUS8]  = AlleleSet(11, 13); p3[LOCUS8]  = AlleleSet(11, 13); // D16
    p1[LOCUS9]  = AlleleSet(10, 10); p2[LOCUS9]  = AlleleSet(11, 11); p3[LOCUS9]  = AlleleSet(10, 10); // D25
    p1[LOCUS10] = AlleleSet(12, 11); p2[LOCUS10] = AlleleSet(10, 13); p3[LOCUS10] = AlleleSet(11, 12); // D19
    p1[LOCUS11] = AlleleSet(12, 12); p2[LOCUS11] = AlleleSet(13, 13); p3[LOCUS11] = AlleleSet(10, 12); // vWA
    p1[LOCUS12] = AlleleSet(12, 13); p2[LOCUS12] = AlleleSet(10, 10); p3[LOCUS12] = AlleleSet(12, 12); // TPOX
    p1[LOCUS13] = AlleleSet(13, 10); p2[LOCUS13] = AlleleSet(13, 10); p3[LOCUS13] = AlleleSet(13, 11); // D18
    p1[LOCUS14] = AlleleSet(13, 11); p2[LOCUS14] = AlleleSet(13, 10); p3[LOCUS14] = AlleleSet(13, 10); // D5
    p1[LOCUS15] = AlleleSet(13, 12); p2[LOCUS15] = AlleleSet(10, 11); p3[LOCUS15] = AlleleSet(13, 13); // P_D
    p1[LOCUS16] = AlleleSet(13, 13); p2[LOCUS16] = AlleleSet(10, 11); p3[LOCUS16] = AlleleSet(13, 13); // P_E

    double theta_val[] = { 0, 0.01, 0.02, 0.03, 0.05 };
    const int ntheta = 5;

    RelType type[] = { ident_t, degree_1_t, sibling_t, degree_2_t, degree_pq_t };
    const int nrels = 5;

    for (int i=0; i<ntheta; ++i)
    {
        double theta = theta_val[i];

        for (int j=0; j<nrels; ++j)
        {
            MatchType mtype_num(type[j], 3, MatchType::INF); // last two args used in D3 case only
            double tol = 1e-6;
            double lr=0, lr2=0, lr_bn=0, lr_bn2=0;

            // P1/P1
            // Using direct Balding-Nichols calculation
            lr_bn = baldingNichols(mtype_num, p1, p1, popdata, theta);

            // Using B11 correction factor
            lr = relmatch(mtype_num, p1, p1, popdata, SubPopModel(SubPopModel::B11, theta));

            CHECK_NEAR(lr_bn, lr, tol);

            bool print = false;
            if (print)
            {
                cout << "B11_profile p1/p1: type = " << setw(5) << mtype_num << " theta = " << theta << " lr = " << lr << endl;
            }

            // P1/P2
            lr_bn = baldingNichols(mtype_num, p1, p2, popdata, theta);
            lr_bn2 = baldingNichols(mtype_num, p2, p1, popdata, theta);
            CHECK_NEAR(lr_bn, lr_bn2, tol);

            lr = relmatch(mtype_num, p1, p2, popdata, SubPopModel(SubPopModel::B11, theta));
            lr2 = relmatch(mtype_num, p2, p1, popdata, SubPopModel(SubPopModel::B11, theta));
            CHECK_NEAR(lr, lr2, tol);

            CHECK_NEAR(lr_bn, lr, tol);
            if (print)
            {
                cout << "B11_profile p1/p2: type = " << setw(5) << mtype_num << " theta = " << theta << " lr = " << lr << endl;
            }

            // P1/P3
            lr_bn = baldingNichols(mtype_num, p1, p3, popdata, theta);
            lr_bn2 = baldingNichols(mtype_num, p3, p1, popdata, theta);
            CHECK_NEAR(lr_bn, lr_bn2, tol);

            lr = relmatch(mtype_num, p1, p3, popdata, SubPopModel(SubPopModel::B11, theta));
            lr2 = relmatch(mtype_num, p3, p1, popdata, SubPopModel(SubPopModel::B11, theta));
            CHECK_NEAR(lr, lr2, tol);

            CHECK_NEAR(lr_bn, lr, tol);
            if (print)
            {
                cout << "B11_profile p1/p3: type = " << setw(5) << mtype_num << " theta = " << theta << " lr = " << lr << endl;
            }

            // P2/P3
            lr_bn = baldingNichols(mtype_num, p2, p3, popdata, theta);
            lr_bn2 = baldingNichols(mtype_num, p3, p2, popdata, theta);
            CHECK_NEAR(lr_bn, lr_bn2, tol);

            lr = relmatch(mtype_num, p2, p3, popdata, SubPopModel(SubPopModel::B11, theta));
            lr2 = relmatch(mtype_num, p3, p2, popdata, SubPopModel(SubPopModel::B11, theta));
            CHECK_NEAR(lr, lr2, tol);

            CHECK_NEAR(lr_bn, lr, tol);
            if (print)
            {
                cout << "B11_profile p2/p3: type = " << setw(5) << mtype_num << " theta = " << theta << " lr = " << lr << endl;
            }
        }

        // Print results for ident, NRC4_4
        // P1/P1
#if 0
        MatchType mtype(ident_t);
        double lr = relmatch(mtype, p1, p1, popdata, SubPopModel(SubPopModel::NRC4_4, theta));
        cout << "B11_profile p1/p1: type = " << setw(5) << mtype << " theta = " << theta << " lr = " << lr << endl;
#endif
    }
}

TEST_FIXTURE(SpreadsheetFixture, spm_nrc44_sib_theta0)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, A);

    MatchType mtype(sibling_t);

    double lr = 30.25; // from Sums spreadsheet, HW case
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::NRC4_4, 0)), tol);
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::NRC4_4, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, spm_nrc44_d2_theta0)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(A, A);

    MatchType mtype(degree_2_t);

    double lr = 5.5; // from Sums spreadsheet, HW case
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::NRC4_4, 0)), tol);
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::NRC4_4, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, spm_nrc44_ff)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(F, F);

    double lr = 1;
    double tol = 1e-6;
    CHECK_CLOSE(lr, relmatch(MatchType(ident_t), p1, p2, popdata, SubPopModel(SubPopModel::NRC4_4, 0.02)), tol);
    CHECK_CLOSE(lr, relmatch(MatchType(ident_t), p2, p1, popdata, SubPopModel(SubPopModel::NRC4_4, 0.02)), tol);

    CHECK_CLOSE(lr, relmatch(MatchType(degree_2_t), p1, p2, popdata, SubPopModel(SubPopModel::NRC4_4, 0.02)), tol);
    CHECK_CLOSE(lr, relmatch(MatchType(degree_2_t), p2, p1, popdata, SubPopModel(SubPopModel::NRC4_4, 0.02)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, spm_B11_ff)
{
    p1[LOCUS1] = AlleleSet(A, A);
    p2[LOCUS1] = AlleleSet(F, F);

    double lr = 1;
    double tol = 1e-6;
    CHECK_CLOSE(lr, relmatch(MatchType(ident_t), p1, p2, popdata, SubPopModel(SubPopModel::B11, 0.02)), tol);
    CHECK_CLOSE(lr, relmatch(MatchType(ident_t), p2, p1, popdata, SubPopModel(SubPopModel::B11, 0.02)), tol);

    CHECK_CLOSE(lr, relmatch(MatchType(degree_2_t), p1, p2, popdata, SubPopModel(SubPopModel::B11, 0.02)), tol);
    CHECK_CLOSE(lr, relmatch(MatchType(degree_2_t), p2, p1, popdata, SubPopModel(SubPopModel::B11, 0.02)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mut_d1_HW)
{
    // Mutation is only important in the D1 case.
    // Note it is NOT symmetrical. It matters which profile belongs to the parent and
    // which to the child

    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(B, B);
    p2[LOCUS1] = AlleleSet(C, C);

    MatchType mtype(degree_1_t);

    double tol = 1e-6; // In the HW case we get the "exact" answer

    // p1 is the father
    double lr = 1e-3 / (2 * 0.3); // r / 2P(C)
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);

    // p2 is the father
    lr = 1e-3 / (2 * 0.2); // r / 2P(B)
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mut_sib_HW)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(B, B);
    p2[LOCUS1] = AlleleSet(C, C);

    MatchType mtype(sibling_t);

    double tol = 1e-2; // Mutation rate should make less than 1% difference in sibling case
                       // we don't really care about the exact answer

    double lr = 0.25; // hand calculation
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AF_AA_HW)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(A, F);
    p2[LOCUS1] = AlleleSet(A, A);

    MatchType mtype(ident_t);

    double tol = 1e-8;
    double Pa = LOCUS1_FREQ[0].m_value;

    double lr = 1/Pa;
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AF_AB_HW)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(A, F);
    p2[LOCUS1] = AlleleSet(A, B);

    MatchType mtype(ident_t);

    double tol = 1e-8;
    double Pa = LOCUS1_FREQ[0].m_value;

    double lr = 1/(2*Pa);
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AF_AF_HW)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(A, F);
    p2[LOCUS1] = AlleleSet(A, F);

    MatchType mtype(ident_t);

    double tol = 1e-8;
    double Pa = LOCUS1_FREQ[0].m_value;

    double lr = (1 + Pa)/(2*Pa);
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AF_BF_HW)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(A, F);
    p2[LOCUS1] = AlleleSet(B, F);

    MatchType mtype(ident_t);

    double tol = 1e-8;

    double lr = 0.5; // constant
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AA_point5_HW)
{
    p1[LOCUS1] = AlleleSet(A, A);

    PMF<Allele> Apmf, A_at_Point5pmf;
    Apmf.insert(make_pair(A, 1.0));
    A_at_Point5pmf.insert(make_pair(A, 0.5));

    vector< PMF<Allele> > av;
    av.push_back(Apmf);            // A
    av.push_back(A_at_Point5pmf);  // 0.5A (+ 0.5F)

    p2[LOCUS1] = AlleleSet(av);

    MatchType mtype(ident_t);

    double tol = 1e-8;

    double Pa = LOCUS1_FREQ[0].m_value;
    double lr = (1+Pa)/(2*Pa*Pa);

    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AA_normalization)
{
    p1[LOCUS1] = AlleleSet(A, A);

    PMF<Allele> Apmf, A_at_Point5pmf;
    Apmf.insert(make_pair(A, 1.0));
    A_at_Point5pmf.insert(make_pair(A, 1.2)); // over-normalized

    vector< PMF<Allele> > av;
    av.push_back(Apmf);            // A
    av.push_back(A_at_Point5pmf);  // 1.2A

    p2[LOCUS1] = AlleleSet(av);

    MatchType mtype(ident_t);

    double tol = 1e-8;

    double Pa = LOCUS1_FREQ[0].m_value;
    double lr = 1/(Pa*Pa); // same as AA v AA

    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, AB_point5_HW)
{
    p1[LOCUS1] = AlleleSet(A, C);

    PMF<Allele> Apmf, A_at_Point5pmf;
    Apmf.insert(make_pair(A, 1.0));
    A_at_Point5pmf.insert(make_pair(C, 0.5));

    vector< PMF<Allele> > av;
    av.push_back(Apmf);            // A
    av.push_back(A_at_Point5pmf);  // 0.5A (+ 0.5F)

    p2[LOCUS1] = AlleleSet(av);

    MatchType mtype(ident_t);

    double tol = 1e-8;

    double Pa = LOCUS1_FREQ[0].m_value;
    double Pb = LOCUS1_FREQ[2].m_value;
    double lr = (1+Pb)/(4*Pa*Pb);

    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::HW, 0)), tol);
}

// LR for BB is father of CC (See Balding Weight-of-evidence Table 7.4)
double
mutTest(double r,     // mutation rate B->C (*2)
        double theta, // Fst
        double fc)    // P(C)
{
    return r*(1+2*theta) / (2*(theta + (1-theta)*fc));
}

// LR for CD is father of AB (See Balding Weight-of-evidence Table 7.4)
// NB the only possible one-step mutation here is C->B
double
mutTest2(double r,     // mutation rate C->B (*2)
         double theta, // Fst
         double fb)    // P(B)
{
    return r*(1+2*theta) / (8*(1-theta)*fb);
}

TEST_FIXTURE(SpreadsheetFixture, mut_d1_B11)
{
    // Mutation is only important in the D1 case.
    // Note it is NOT symmetrical. It matters which profile belongs to the parent and
    // which to the child

    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(B, B);
    p2[LOCUS1] = AlleleSet(C, C);

    MatchType mtype(degree_1_t);

    double tol = 1e-6;

    double theta = 0.02;

    // p1 is the father
    double lr = mutTest(1e-3, theta, 0.3);
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::B11, theta)), tol);

    // p2 is the father
    lr = mutTest(1e-3, theta, 0.2);
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::B11, theta)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mut_sib_B11)
{
    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(B, B);
    p2[LOCUS1] = AlleleSet(C, C);

    MatchType mtype(sibling_t);

    double tol = 1e-2; // Mutation rate should make less than 1% difference in sibling case
                       // we don't really care about the exact answer

    double lr = 0.25; // hand calculation
    CHECK_NEAR(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::B11, 0.01)), tol);
    CHECK_NEAR(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::B11, 0.01)), tol);
}

TEST_FIXTURE(SpreadsheetFixture, mut_d1_B11_b)
{
    // Mutation is only important in the D1 case.
    // Note it is NOT symmetrical. It matters which profile belongs to the parent and
    // which to the child

    p1.setMutRate(1e-3);
    p2.setMutRate(1e-3);

    p1[LOCUS1] = AlleleSet(A, B);
    p2[LOCUS1] = AlleleSet(C, D);

    MatchType mtype(degree_1_t);

    double tol = 1e-6;

    double theta = 0.02;

    // p1 is the father
    double lr = mutTest2(1e-3, theta, 0.3);
    CHECK_CLOSE(lr, relmatch(mtype, p1, p2, popdata, SubPopModel(SubPopModel::B11, theta)), tol);

    // p2 is the father
    lr = mutTest2(1e-3, theta, 0.2);
    CHECK_CLOSE(lr, relmatch(mtype, p2, p1, popdata, SubPopModel(SubPopModel::B11, theta)), tol);
}

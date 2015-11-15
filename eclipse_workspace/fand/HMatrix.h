/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * HMatrix.h
 *
 *  Created on: Jan 7, 2010
 *      Author: gareth
 */

#ifndef HMATRIX_H_
#define HMATRIX_H_

#include <memory>
#include "PMF.h"
#include "Allele.h"
#include "SubPopModel.h"

 // forward declarations
class BackFreq;

struct MatchType; class Profile; class PopulationData;

typedef double termFunc(Allele const &, Allele const &, BackFreq const &, Allele const &, Allele const &);

// HMatrix (Half-Matrix) is the Upper Triangular Square Matrix representing
// the PMF of a diploid genotype at a particular locus.
class HMatrix
{
public:
	HMatrix();

	// construct as product of pair of PMF<Allele>
	HMatrix(const PMF<Allele> &pmf1, const PMF<Allele> &pmf2);

	virtual ~HMatrix() {}

	int size() const { return m_pmf.size(); }
	int empty() const { return m_pmf.empty(); }
	void clear() { m_pmf.clear(); ive_changed(); }
    double sum() const { return m_pmf.sum(); }

	// get/set a value
	double operator()(const Allele &a1, const Allele &a2) const;
	void set(const Allele &a1, const Allele &a2, double d);

	// return marginal probability of an allele
	double margP(const Allele &a) const;

	bool normalize();

	// Output as vector of floats. NB in sorted order.
	const std::vector<float> &vec() const;

	// raw float array
	const float *bin() const;

	HMatrix& operator=(const HMatrix &h2);
    HMatrix& operator+=(const HMatrix &h2);

	bool operator==(const HMatrix &h2) const { return m_pmf == h2.m_pmf; }

protected:

	// construct from one PMF<Allele> as a background distribution
    HMatrix(const PMF<Allele> &pmf, SubPopModel const &spm /* = SubPopModel() */);

    PMF< std::pair<Allele, Allele> >::const_iterator
    find(const Allele &a1, const Allele &a2) const;

    // all non-const functions should end with a call to ive_changed
    // so that derived classes can re-calculate stuff
    virtual void ive_changed() { }

    PMF< std::pair<Allele, Allele> > m_pmf;

private:
//    HMatrix(const HMatrix&); // use default copy ctor

	void makeBin() const;

	// element-wise add
	friend HMatrix operator+(const HMatrix &h1, const HMatrix &h2);

	// element-wise subtract
	friend HMatrix operator-(const HMatrix &h1, const HMatrix &h2);

	// multiply by constant
	friend HMatrix operator*(const HMatrix &h1, double d);
    friend HMatrix operator/(const HMatrix &h1, double d);
	friend HMatrix operator*(double d, const HMatrix &h1);

	// element-wise multiply
	friend HMatrix operator*(const HMatrix &h1, const HMatrix &h2);

	friend HMatrix & operator*=(HMatrix & h, double d);
    friend HMatrix & operator/=(HMatrix & h, double d);
    friend HMatrix & operator+=(HMatrix & h, double d);
    friend HMatrix & operator-=(HMatrix & h, double d);

	// element-wise divide
	friend HMatrix operator/(const HMatrix &h1, const HMatrix &h2);

	// element-wise multiply then sum
	friend double dot(const HMatrix &h1, const HMatrix &h2);
    friend HMatrix halveOffDiag(HMatrix const &h);
    friend double hdot(const HMatrix &h1, const HMatrix &h2);

	// TODO rethink this?
    friend HMatrix mutate_1s(HMatrix const &prof, const PMF<Allele> &background, double r);
	friend HMatrix sib(HMatrix const &prof, const PMF<Allele> &back);
	friend HMatrix relNHW(HMatrix const &prof, const BackFreq &back, termFunc rel_ij);
	friend HMatrix degree1(HMatrix const &prof, const PMF<Allele> &back);
	friend HMatrix degree(int n, HMatrix const &prof, const PMF<Allele> &back);
	friend void pcontrib(HMatrix const &prof, PMF<Allele> &contrib);
	friend HMatrix genRel(float a1, float b1, float a2, float b2, HMatrix const &prof, const PMF<Allele> &back);
    friend HMatrix invRel(float a1, float b1, float a2, float b2, HMatrix const &prof, const PMF<Allele> &back);
    friend HMatrix makeHMatrixNHW(const PMF<Allele> &pmf1, const PMF<Allele> &pmf2, const BackFreq &hback, float delta, bool sparse);

    friend HMatrix FSTCorrection(HMatrix const &prof, BackFreq const &hback, double fst);
    friend HMatrix condProb(BackFreq const &hback, Allele const &i, Allele const &j, double fst);

    friend std::ostream & operator<<(std::ostream &os, const HMatrix &a);

    friend BackFreq calcBackForAAAATEST(PMF<Allele> const &back, SubPopModel const &spm);
    friend double dirichletModBackground(HMatrix const &h1, Allele const &i, double fi, double fst);
    friend double baldingNichols(
            MatchType const &mtype,
            Profile const &p1,
            Profile const &p2,
            PopulationData const &popdata,
            double theta);
    friend HMatrix
    relBN(MatchType type,
          HMatrix const &prof,
          const PMF<Allele> &back,
          double fst);
    friend HMatrix
    relBN(MatchType type,
          Allele const &p,
          Allele const &q,
          const PMF<Allele> &back,
          double fst);
    friend double
    sibBN(Allele const &p,
          Allele const &q,
          Allele const &r,
          Allele const &s,
          double fp,
          double fq,
          double fst);
    friend double
    DnBN(int n,
         Allele const &p,
         Allele const &q,
         Allele const &r,
         Allele const &s,
         double fr,
         double fs,
         double fst);


	mutable std::vector<float> m_bin; // cached raw array of values (in sorted order)
};

// this is just to allow comparison of HMatrices in unit tests
inline std::ostream &
operator<<(std::ostream &os, const HMatrix &a)
{
    os << a.m_pmf;
	return os;
}

// class representing genotype background frequencies (as calculated from a subpopulation model)
// Backfreq is an HMatrix and is therefore symmetrical

class BackFreq : public HMatrix
{
public:

    explicit
    BackFreq(const PMF<Allele> &pmf, SubPopModel const &spm = SubPopModel());

    // ? should we have copy ctor and assign that call ive_changed ?

    virtual ~BackFreq() {}

    // probability of a (unconditionally, same as original pmf)
    double operator()(Allele const &a) const;

    // probability of g (allele pair, in either order)
    double operator()(std::pair<Allele, Allele> const &g) const;

    // probability of g (allele pair, in given order)
    double pOrdered(std::pair<Allele, Allele> const &g) const;

    // probability of p given q (conditional) // NB overrides HMatrix operator(p,q)
    double operator()(Allele const &p, Allele const &q) const;

    PMF<Allele> const & alleleFreqs() const;

protected:
    PMF<Allele> m_allele_pmf;

    virtual void ive_changed();
};

#endif /* HMATRIX_H_ */

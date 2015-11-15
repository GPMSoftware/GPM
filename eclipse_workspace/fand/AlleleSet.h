/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * AlleleSet.h
 *
 *  Created on: Dec 22, 2009
 *      Author: gareth
 */

#ifndef ALLELESET_H_
#define ALLELESET_H_

#include "PMF.h"
#include "Allele.h"
#include "HMatrix.h"
#include "ProfileData.h"

#include <iostream>

// The set of Alleles at a particular locus in a Profile.
//
// AlleleSet holds a PMF for each allele believed to be present in the simple.
// NB there may be more than two alleles in the case of mixtures or tri-allele patterns
//

struct MatchType;

class AlleleSet
{
private:

public:
	// needs default c'tor because used as a value type in a map that uses operator[] to create new elements
	AlleleSet();

	// simple case - construct from two alleles
	AlleleSet(const Allele a1, const Allele a2);

	// mixture - construct from 4 alleles
	AlleleSet(const Allele a1, const Allele a2, const Allele a3, const Allele a4);

	// general case - construct from a vector of PMFs
	// (must contain at least 2 PMFs)
	AlleleSet(std::vector< PMF<Allele> > const &pmf_list);

	// the number of alleles in the profile
	int size() const { return m_pmfs.size(); }

	bool checkBackground(PMF<Allele> const &background) const;

	// return an allele as a PMF
	const PMF<Allele>& getAllele(int n) const;

    // return AlleleSet as a vector of PMFs
    const std::vector< PMF<Allele> >& getAllAlleles() const;

    // set an allele PMF
    void setAllele(int n, PMF<Allele> const &pmf);

    // check if this is a 'simple' AlleleSet
	bool simple(Allele &a1, Allele &a2) const;

	// add an allele PMF
	void pushBack(const PMF<Allele> a) { m_pmfs.push_back(a); reset(); cacheSimple(); }

	// get the HMatrix for this profile
	// NB result will be cached until reset() is called

	// get the HMatrix for this profile, non-HW case
	const HMatrix &getMatchHMatrixNHW(
	        const BackFreq &hback,
	        const ProfileData &data = ProfileData(),
	        bool sparse = false) const;

	// get the HMatrix for this profile, HW case (cached)
    const HMatrix &getMatchHMatrix(
            const PMF<Allele> &background,
            const ProfileData &data = ProfileData(),
            bool sparse = false) const;

	// get the HMatrix for sibling of this profile, where background gives HW gene frequencies
	// NB result will be cached until reset() is called
	const HMatrix &getSibHMatrix(
	        HMatrix const &prof,
	        const PMF<Allele> &background,
	        const ProfileData &data = ProfileData(),
	        bool sparse = false) const;

	// get the HMatrix for relatives (HW)
	HMatrix
	getRelHMatrix(
	        MatchType const &match_type,
	        const PMF<Allele> &background,
	        const ProfileData &data,
	        bool sparse) const;

	// get the B11 subpop model correction matrix
    // NB result will be cached until theta changes or reset() is called
	const HMatrix &
	getSpmcHMatrix(
	        const ProfileData &data,
	        const PMF<Allele> &background,
	        SubPopModel const &spm) const;

	// reset cached HMatrices
    // this cacheing strategy depends upon the Profile calling reset whenever the ProfileData changes
	// NB thread safety for cacheing and recalculating in GPU solutions is provided by Profile
	void reset() const;

   // use these functions to set the parameters data, or cached matrices will not be recalculated
    void setTheta(double val) const;
    void setData(const ProfileData &data) const;

private:
    void makeMatchHMatrix(const PMF<Allele> &background, bool sparse) const;

	void makeSibHMatrix(HMatrix const &prof, const PMF<Allele> &background, bool sparse) const;

    HMatrix identNHW(const BackFreq &hback, bool sparse) const;

	void cacheSimple();

	// A PMF for each Allele believed to exist in the sample
	// (there may be more then two in the case of mixtures or tri-allele patterns)
	std::vector< PMF<Allele> > m_pmfs;

	// simple case optimization
	bool m_simple;
	Allele m_a1, m_a2;

	// Cached HMatrices. These are cleared when the data they depend on changes.
	// NB the match matrices are used in calculating other matrices, so it is worth cacheing them.
	// The subpopulation correction matrix is used many times and there is a big impact on speed
	// if it is not cached.
	// We cache the sibling matrix, and currently no other relationship matrices. This could be done -
	// it is a trade-off between speed and memory. They may not be used more then once.

	mutable HMatrix m_match;        // Ident matrix. Depends upon error rate
    mutable HMatrix m_match_sparse; // Sparse Ident matrix. Depends upon error rate
    mutable HMatrix m_match_nhw;    // Sparse non-HW Ident matrix. Depends upon error rate

	mutable HMatrix m_sib;          // Sib matrix (HW). Depends upon error rate and mutation rate
	mutable HMatrix m_sib_sparse;   // Sparse Sib matrix (HW). Depends upon error rate and mutation rate

	mutable HMatrix m_spmc;         // Subpopulation model correction matrix (sparse)
                                    // Depends upon error rate and theta
                                    // (Cacheing this matrix is very important!)
                                    // NB the GPU solution does not use m_spmc or m_theta

	// Data the cached matrices depend on
	mutable double m_error_rate;
	mutable double m_mutation_rate;
	mutable double m_theta;
	mutable int m_num_contributors;

    // called by setData
    void setErrorRate(double val) const;
    void setMutRate(double val) const;
    void setNumContribs(int val) const;
};

#endif /* ALLELESET_H_ */

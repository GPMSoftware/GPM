/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * match.h
 *
 *  Created on: Dec 9, 2009
 *      Author: gareth
 *
 *      Match profiles on the CPU
 *      (Allows parts of the calculation to be done on the GPUs for testing purposes only - see #defines at the top of match.cpp)
 *
 *      NB For full GPU solutions see cuda_accel library
 */

#ifndef MATCH_H_
#define MATCH_H_

#include "Profile.h"
#include "populationdata.h"
#include "loci.h"
#include "cuda_accel/cuda_accel.h"

struct MatchType : public CUDAMatchType
{
	enum  { INF = INT_MAX };

	explicit MatchType(RelType rel, int p1=INF, int p2=INF);  // Dnm with integer n, m
    explicit MatchType(RelType rel, float a1, float b1, float a2, float b2);
    explicit MatchType(RelType rel, double n, double m);  // Dnm with fractional n, m

	std::string string() const;

	int generations() const; // the number of generations (for mutation purposes)
};

std::ostream &
operator<<(std::ostream &os, const MatchType &a);

bool
operator<(const MatchType &m1, const MatchType &m2);

// form a degree(p,q) HMatrix for prof
HMatrix
degree(
	int p, int q,        // degree of relationship
	HMatrix const &prof,
	const PMF<Allele> &back);

// form a generalized R.C. HMatrix for prof
HMatrix
genRel(
    float a1,
    float b1,
    float a2,
    float b2,
    HMatrix const &prof,
    const PMF<Allele> &back);

HMatrix
invRel(
    float a1,
    float b1,
    float a2,
    float b2,
    HMatrix const &prof,
    const PMF<Allele> &back);

double // likelihood ratio
match(const Profile &p1, const Profile &p2, SubPopModel const &spm = SubPopModel(), PopulationData const &popdata = populationData());

double // likelihood ratio
match(const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm = SubPopModel());

double // likelihood ratio
sibmatch(const Profile &p1, const Profile &p2, SubPopModel const &spm = SubPopModel(), PopulationData const &popdata = populationData());

double // likelihood ratio
sibmatch(const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm = SubPopModel());

double
relmatch(const MatchType &match_type, const Profile &p1, const Profile &p2, SubPopModel const &spm = SubPopModel(), PopulationData const &popdata = populationData());

double
relmatch(const MatchType &match_type, const Profile &p1, const Profile &p2, PopulationData const &popdata, SubPopModel const &spm = SubPopModel());

void
getIBD(MatchType const &mtype, double &k0, double &k1, double &k2);

#endif /* MATCH_H_ */

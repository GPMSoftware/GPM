/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Profile.h
 *
 *  Created on: Nov 30, 2009
 *      Author: gareth
 */

#ifndef PROFILE_H_
#define PROFILE_H_

#include "ProfileData.h"
#include "AlleleSet.h"
#include "loci.h"
#include "Cached.h"

#include <map>

typedef std::map<Locus, AlleleSet> LocusSet;

class Profile // Identifiler profile
{
public:
	explicit Profile(ProfileData prof_data, const LocusSet &loci = LocusSet());
	Profile(Profile const &other);
	Profile & operator=(Profile const &other);

	// get binary representation of profile (cached)
//	CachedArray<float>::ConstTempRef binRef(PopulationData const &popdata) const;
    std::vector<float> const & bin(PopulationData const &popdata, SubPopModel const &spm = SubPopModel()) const;

	AlleleSet& operator[](int loc);

	const AlleleSet& operator[](int loc) const;

	bool empty() const { return m_loci.empty(); }
	int numLoci() const { return m_loci.size(); }
	LocusSet::const_iterator begin() const { return m_loci.begin(); }
	LocusSet::const_iterator end()   const { return m_loci.end(); }
	LocusSet::const_iterator find(Locus locus) const { return m_loci.find(locus); }
	void clear() { m_loci.clear(); }

	const ProfileData& data() const    { return m_data; }
	void setData(ProfileData const &d);
	void setErrorRate(double d);
	void setMutRate(double d);
    void setNumContributors(int n);

private:
	void makeBin(PopulationData const &popdata, SubPopModel const &spm) const;
    void reset() const;
	void resetH() const;

    ProfileData m_data;

	LocusSet  m_loci;

	static const AlleleSet m_empty;

	// cached state
//	mutable CachedArray<float> m_bin; // Cached binary data
    mutable std::vector<float> m_bin; // binary data
	mutable long m_hash;              // Population data hash (to see if cached data is still current)
};

#endif /* PROFILE_H_ */

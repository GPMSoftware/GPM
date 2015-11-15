/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * populationdata.h
 *
 *  Created on: Dec 7, 2009
 *      Author: gareth
 */

#ifndef POPULATIONDATA_H_
#define POPULATIONDATA_H_

#include "PMF.h"
#include "Allele.h"
//#include "ProfileData.h"
#include "loci.h"
#include "util.h"
#include "cuda_accel/cuda_accel.h"
#include <vector>
#include <iostream>
#include <fstream>

// how to treat unknown alleles
enum Unknowns
{
    read_from_env = 0,
	ignore,
	add_as_rare,
};

class PopulationData;

class LocusInfo: public CudaLocusInfo
{
public:
    LocusInfo() {}

private:
    void update(PopulationData const &p);

    friend class PopulationData;
};

std::ostream & operator<<(std::ostream &os, const CudaLocusInfo& loc_info);

class PopulationData
{
public:
	PopulationData(Unknowns policy = read_from_env);

	PopulationData(std::string const &path, Unknowns policy = read_from_env);

	bool read(std::string const &path);

    bool write(std::string const &path) const;

    int numLoci() const;

    // normalize all Loci (return true on success)
    bool normalize();

	bool hasLocus(int const &loc) const;

	const PMF<Allele>& operator()(int const &loc) const;

	void setLocus(int const &loc, PMF<Allele> const &pmf);

	// NB if allele is not in the population database it may be added (depending on policy of the database) so NOT const
	// If added_total is non-null the number of alleles added is added to *added_total
	double getFrequency(int const &loc, Allele const &a, int *added_total = 0);

	// get the locus info associated with this population database
	const LocusInfo& getLocusInfo() const;

	void clear();

//	bool operator<(PopulationData const &other) const;

	int sampleSize() const;

	void setSampleSize(int n);

    long int getHash() const;

    void setPolicy(Unknowns policy);

    Unknowns getPolicy() const;

    std::string path() const;

private:
	void getPolicyFromEnv();

	std::map< Locus, PMF<Allele> > m_data;
	int m_sample_size;
	Unknowns m_unknowns_policy; // policy for handling unknown alleles
	LocusInfo m_locus_info;
    long int m_hash;
    std::string m_path;

	void reHash();
};

// get default population data
PopulationData const & populationData();

// get set of all known Loci (with no allele data)
PopulationData const & knownLoci();

// get a test population data set
PopulationData & testPopulationData();

// set default population data
void popSet(const PopulationData &newpop);

// clear default population data
void popClear();

// print population database to file (debugging)
void popPrint(const PopulationData &pop);

#endif /* POPULATIONDATA_H_ */

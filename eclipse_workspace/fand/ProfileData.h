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

#ifndef PROFILEDATA_H_
#define PROFILEDATA_H_

#include <string>

enum ProfileType
{
	prof_type_unknown = 0,
	Identifiler,
	PowerPlex16,
	SGM_Plus,
	mixed,
	test,
	num_profile_types
};

enum EvidenceType
{
    ev_type_unknown = 0,
    crime,
    reference
};

//// known loci (all kits)
//enum Locus
//{
//	D8S1179 = 0,
//	D21S11,
//	D7S820,
//	CSF1PO,
//	D3S1358,
//	THO1,
//	D13S317,
//	D16S539,
//	D2S1338,
//	D19S433,
//	vWA,
//	TPOX,
//	D18S51,
//	AMEL,
//	D5S818,
//	FGA,
//	PENTA_D,
//	PENTA_E,
////	D1S, D12S, D2S, D10S, D22S, SE33, DYS391
//	D1S1656,
//	D12S391,
//	D10S1248,
//	D22S1045,
//	SE33,
//	DYS391,
//	num_loci,
//	locus_none = -1
//};

struct ProfileData
{
	ProfileData(ProfileType         type            = prof_type_unknown,
			    const std::string  &name            = "No name",
			    int                 ncontrib        = 1,
			    float               delta           = 0,
			    float               mut             = 0,
                EvidenceType        evidence_type   = ev_type_unknown,
                const std::string   dataset         = "",
                const std::string   sample_id       = "",
                const std::string   profile_id      = "")
	: m_kit_type(type)
	, m_evidence_type(evidence_type)
	, m_id(name)
	, m_dataset(dataset)
	, m_sample_id(sample_id)
	, m_profile_id(profile_id)
	, m_num_contributors(ncontrib)
	, m_error_rate(delta)
	, m_mutation_rate(mut)
	{
	}

	ProfileType  m_kit_type;
	EvidenceType m_evidence_type;
	std::string  m_id;
    std::string  m_dataset;
    std::string  m_sample_id;
    std::string  m_profile_id;
	int          m_num_contributors;
	float        m_error_rate;
	float        m_mutation_rate;
};

ProfileType
kitID(std::string kitname);

std::string
kitName(ProfileType prof_type);

EvidenceType
evidenceID(std::string evidence_type);

std::string
evidenceName(EvidenceType ev_type);

std::ostream &
operator<<(std::ostream &os, ProfileData const &p);

#endif /* PROFILE_DATA_H_ */

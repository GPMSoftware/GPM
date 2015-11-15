/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * cuda_match.h
 *
 *  Created on: Mar 17, 2010
 *      Author: gareth
 *
 *  CUDA match functions
 */

#ifndef CUDA_MATCH_H_
#define CUDA_MATCH_H_

#include "cuda_accel.h"
#include "GPUDevices.h"
#include "fand/Profile.h"
#include "fand/match.h"

#include <vector>

int cuda_match_n(
		ProfileRange & db,
		const Profile& p,
		struct NResults &results,
		double lr_threshold,
	    double delta,
		const MatchType &match_type,
	    SubPopModel const &spm);

int cuda_match_n2(
		ProfileRange &db,
		struct N2Results &results,
		double lr_threshold,
		double delta,
		MatchType const &match_type,
	    SubPopModel const &spm);

int cuda_match_nm(
		ProfileRange &db1,
		ProfileRange &db2,
		struct N2Results &results,
		double lr_threshold,
		double delta1,
		double delta2,
        MatchType const &match_type,
	    SubPopModel const &spm,
		ArrayPart part);

// streaming solution (n2, nm, ident, kinship)
void cuda_stream_match(
		ProfileRange &db1,
		ProfileRange &db2,
		MatchType const &match_type,    // Relative type
		DProfile back,                  // Background
		CudaLocusInfo const &loc_info,  // Alleles per locus etc
		struct N2Results &results,
		float lr_threshold,
		double delta1,
		double delta2,
		ArrayPart part,
		PopulationData const &mpopdata,
		SubPopModel spm);

#endif /* CUDA_MATCH_H_ */

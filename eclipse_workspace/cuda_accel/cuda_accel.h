/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * cuda_accel.h
 *
 *  Created on: Jan 19, 2010
 *      Author: gareth
 *
 *  Cuda kernels and associated functions
 *
 */

#ifndef CUDA_ACCEL_H_
#define CUDA_ACCEL_H_

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <math.h>
#include <memory>
#include <vector>

// representation of profile on device

// NB we allocate two arrays of this size PER THREAD (on the device) in make_sib and make_2path.
// Automatic array variables are stored in global memory (not registers), but need to be a fixed size,
// so we can afford to make these biggish.
// NBB we also need 2 arrays of size CUDA_MAX_ALLELES * (CUDA_MAX_ALLELES + 1) / 2 in make_invrc
static const int CUDA_MAX_ALLELES  = 80; // the biggest we can allocate in make_invrc

//static const int cuda_num_loci  = 18; // 16 for Identifiler
static const int cuda_num_loci  = 25; // now including NGM loci

static const int blockSize = 8;

// DProfile points to a concatenated array of half-matrices of size cuda_locus_size (one for each locus)
struct CUDAData
{
	float *data;
	int size;
};

struct ConstCUDAData
{
	const float *data;
	int size;
};

struct CudaLocusInfo
{
	int num_alleles[32];  // number of alleles at each locus (everything can be calculated from this - see LocusInfo::update)
	int locus_size[32];   // number of frequencies in half-matrix representation of each locus
	int locus_offset[32]; // offset to start of each locus in DProfile
	int back_offset[32];  // offset to start of each locus in background vector
	int profile_size;     // number of frequencies in DProfile
	int back_size;        // number of frequencies in background vector
};

struct CudaSubPopModel
{
	enum Type
	{
		HW      = 0,  // Hardy-Weinberg (no sub-population model)
		NRC4_4,       // NRC II equations 4.4 (TODO: currently this means use 4.4 for IDENT and B11 otherwise)
		NRC4_10,      // NRC II equations 4.10
		B11,          // Full treatment using correction factor from "theorem 1"
	};

	Type type;
	float theta_bar;
};

// public derivation seems to work despite being C++!
struct DProfile : public CUDAData
{
	DProfile(float *d=0, int n=0) { data=d; size=n; }
};

struct ConstDProfile : public ConstCUDAData { };

// DBackground points to a concatenated array of background vectors of size given by
// CudaLocusInfo::num_alleles (one for each locus)
struct DBackground : public CUDAData { };

// the part of an N*N array to calculate on
enum ArrayPart
{
	lower = -1,  // includes diagonal
	full  =  0,
	upper =  1   // excludes diagonal
};

// results of an N search. Room for up to nresult_max results.
// NB count > nresult_max indicates overflow has occurred: the first nresult_max are saved.
static const int nresult_max = 100000;
struct NResults
{
	int count;

	int index[nresult_max];
	float lr[nresult_max];
};

// results of an N2 search. Room for up to n2result_max results.
// NB count > n2result_max indicates overflow has occurred: the first n2result_max are saved.
static const int n2result_max = 100000;
struct N2Results
{
	int count;

	int index1[n2result_max];
	int index2[n2result_max];
	float lr[n2result_max];
};

enum RelType
{
	ident_t,       // identity (or identical twin)
	degree_1_t,    // parent, child
	sibling_t,     // full sibling
	degree_2_t,    // grandparent, grandchild, half-sib
	degree_pq_t,   // degree(p,q)
	gen_t,         // general (4-number) relationship
	inv_t,         // inverse of general (4-number) relationship
	none_t,        // unrelated
	unknown_t
};

struct CUDAMatchType
{
	enum  { INF = INT_MAX };

	RelType m_rel_type;
	int m_path1steps;
	int m_path2steps;
	float m_a1, m_b1, m_a2, m_b2; // generalized relationship coefficient
};

int pindex(int i, int j, int n);

// triangular number
int tri(int i);

// inverse triangle functions.
// lower_tri: return greatest i such that tri(i) <= k
// NB 0 < k < tri(N-1)

// implementation 1: binary chop
// NB if we don't know N we use k+1 as our first high value
int lower_tri1(int k, int n=0);

// implementation 2: quadratic (solve n(n+1) = 2k for n)
// (this looks simpler but is actually slower)
int lower_tri2(int k);

// map from k (the linear index in an upper-triangular DProfile) to i, j coordinates in an n*n matrix
// (the inverse of pindex)
//
//   j->
// i 0 1 2 3
// |   4 5 6
// v     7 8
//         9
#define lower_tri lower_tri1
inline void inv_p(int n, int k, int &i, int &j)
{
	// method: number in reverse order:
	// k + k_ = tri(n)-1
	// i + i_ = j + j_ = n-1

	int k_ = tri(n) - k - 1;
	int i_ = lower_tri(k_, n);
    int j_ = k_ - tri(i_);
    j = n - 1 - j_;
    i = n - 1 - i_;
}
#undef lower_tri

// construct sibling
void
make_sib(
	DProfile prof,            // profile dataset
    DBackground back,         // background vector
	DProfile psib,            // sibling dataset (to construct)
	CudaLocusInfo const *loc_info); // locus info

void
make_sib_BN(
	DProfile prof,            // profile dataset
    DBackground back,         // background vector
	DProfile psib,            // sibling dataset (to construct)
	CudaLocusInfo const *loc_info, // locus info
	float theta);

// construct two path relative (other than sibling)
void
make_2path(
	DProfile prof,            // profile dataset
	int n,                    // first degree of relationship
	int m,                    // second degree of relationship (MatchType::INF == F)
    DBackground back,         // background vector
	DProfile prel,            // relative dataset (to construct)
	CudaLocusInfo const *loc_info); // locus info

// construct relative from generalized relationship coefficient
void
make_genrc(
	DProfile prof,       // profile dataset
	float a1,            // proportion of 'a' inherited on side 1
	float b1,            // proportion of 'b' inherited on side 1
	float a2,            // proportion of 'a' inherited on side 2
	float b2,            // proportion of 'b' inherited on side 2
    DBackground back,    // background vector
	DProfile prel,       // relative dataset (to construct)
	CudaLocusInfo const *loc_info);

// construct the inverse of the above
void
make_invrc(
		DProfile prof,       // profile dataset
		float a1,            // proportion of 'a' inherited on side 1
		float b1,            // proportion of 'b' inherited on side 1
		float a2,            // proportion of 'a' inherited on side 2
		float b2,            // proportion of 'b' inherited on side 2
		DBackground back,    // background vector
		DProfile prel,       // relative dataset (to construct)
		CudaLocusInfo const *loc_info);

// n-match
void cuda_n_match(
		std::vector<DProfile> &prof1,
		DProfile prof2,
		DProfile back,
		CudaLocusInfo const &loc_info,
		struct NResults &results,
		float lr_threshold);

// n2-match
#define STREAMED
#ifndef STREAMED
void cuda_n2match(
		std::vector<DProfile> &prof_db, // Profile database
		DProfile back,                  // Background
		CudaLocusInfo const &loc_info,  // Alleles per locus etc
		struct N2Results &results,
		float lr_threshold);

void cuda_nm_match(
		std::vector<DProfile> &prof_db1, // Profile 1 database
		std::vector<DProfile> &prof_db2, // Profile 2 database
		DProfile back,                   // Background
		CudaLocusInfo const &loc_info,   // Alleles per locus etc
		struct N2Results &results,
		float lr_threshold,
		ArrayPart part);

void cuda_n2relmatch(
		std::vector<DProfile> &prof_db,    // Profile database
		CUDAMatchType   const &match_type, // Relative type
		DProfile               backh,      // Background as a half-matrix
		CudaLocusInfo const &loc_info,     // Alleles per locus etc
        DBackground            backv,      // Background as a vector (to allow siblings to be constructed)
		struct N2Results      &results,
		float                  lr_threshold);

// match REL(prof_db1) against prof_db2
void cuda_nm_relmatch(
	std::vector<DProfile> &prof_db1,   // Profile database
	std::vector<DProfile> &prof_db2,   // Profile database
	CUDAMatchType   const &match_type, // Relative type
	DProfile               backh,      // Background as a half-matrices
	CudaLocusInfo const &loc_info,     // Alleles per locus etc
    DBackground            backv,      // Background as vectors
	N2Results &results,
	float lr_threshold,
	ArrayPart part);
#endif

void setDevice(int device);

bool isContext();

int getDevice();

void copyOffsets(CudaLocusInfo const &locus_info);

void
runCudaMatch6(
		int 			nBlocks1,
		int 			nBlocks2,
		int 			blockSize,
		cudaStream_t	stream,
		DProfile 	   *profdb1_d,
		int 			db1_chunk_size,
		DProfile 	   *profdb2_d,
		int 			db2_chunk_size,
		DProfile 		back,
		ArrayPart 		part,
		N2Results	   *results_d,
		float 			lr_threshold,
		int 			i,
		int 			j,
		DProfile 	   *spmc1_dp = 0); // SPM correction matrices (for profdb1)

void
runCudaSPMC(
	int 		nBlocks,
	int 		blockSize,
	DProfile 	*prof_db,    // profile dataset
	DProfile 	*spmc_db,    // correction factor matrices (to construct)
	int 		n,           // size of dataset
	DProfile    back,        // background matrix (4.4)
    DBackground backv,       // background vector (HW)
    float       theta);      // Fst

void
runCudaSib(
		int 		nBlocks,
		int 		blockSize,
		DProfile 	*prof_db,    // profile dataset
		DProfile 	*sib_db,     // sibling dataset (to construct)
		int 		n,           // size of dataset
		DBackground back,        // background
		const CudaSubPopModel &spm); // subpopulation model

void
runCudaGenrc(
	int 		nBlocks,
	int 		blockSize,
	DProfile 	*prof_db,    // profile dataset
	DProfile 	*rel_db,     // relative dataset (to construct)
	int 		n,           // size of dataset
	float 		a1,          // proportion of 'a' inherited on side 1
	float 		b1,          // proportion of 'b' inherited on side 1
	float 		a2,          // proportion of 'a' inherited on side 2
	float 		b2,          // proportion of 'b' inherited on side 2
	DBackground back,        // background
	bool inverse = false);

void
runCudaRel(
	int 		nBlocks,
	int 		blockSize,
	DProfile 	*prof_db,    // profile dataset
	DProfile 	*rel_db,     // relative dataset (to construct)
	int 		n,           // size of dataset
	int 		d1,          // degree first path
	int 		d2,          // degree second path
	DBackground back,        // background
	const CudaSubPopModel &spm); // subpopulation model



#endif /* CUDA_ACCEL_H_ */

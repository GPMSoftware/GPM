/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * cuda_accel.cpp
 *
 *  Created on: Jan 19, 2010
 *      Author: gareth
 *
 *  Cuda kernels and associated functions
 *
 *  ! must be compiled with nvcc
 *  NB atomic...() functions require -arch=sm_11 flag to nvcc. This is not documented.
 */

#include <cuda.h>
#include "fand/Assert.h"
#include "fand/util.h"
#include "cuda_accel.h"
#include <iostream>
#include <math.h>

#include "fand/MessageStream.h"
INIT_MESSAGES("cuda_accel");
#include "fand/messages.h"

__constant__ CudaLocusInfo loc_info;

#define COPY_OFFSETS() Assert(cudaMemcpyToSymbol("loc_info", &locus_info, sizeof(CudaLocusInfo), 0, cudaMemcpyHostToDevice) == cudaSuccess)

// delta function
#define d(a, b) ((int)(a==b))

__host__
void copyOffsets(CudaLocusInfo const &locus_info)
{
	COPY_OFFSETS();
}

__host__
bool isContext() // is there a CUDA context current in ths thread?
{
//	std::cout << "isContext" << std::endl;
	CUresult cur;
	CUcontext ctx;
	cur = cuCtxPopCurrent(&ctx);
    return (cur == CUDA_SUCCESS);
}

__host__
void setDevice(int device)
{
//	std::cout << "setDevice" << std::endl;

	// check if there is a CUDA context already
    Assert2( ! isContext(), "cudaSetDevice: there is already a context");

    cudaError_t err;
    err = cudaSetDevice(device);
    Assert2(err == cudaSuccess, "cudaSetDevice failed");
}

__host__
int getDevice()
{
	int device = -1;
    cudaError_t err;
    err = cudaGetDevice(&device);
    Assert2(err == cudaSuccess, "cudaGetDevice failed");
    return device;
}

// triangular number
__host__ __device__
int tri(int i)
{
    return i * (i+1) / 2;
}

// inverse triangle function.
// lower_tri: return greatest i such that tri(i) <= k
// NB 0 < k < tri(N-1)

// implementation 1: binary chop
// NB if we don't know N we use k+1 as our first high value
__host__ __device__
int lower_tri1(int k, int n)
{
    int i_low = 0;
    int i_high = n? (n-1) : (k+1);
    int i;
    while ((i_high - i_low) > 1)
    {
        i = (i_low + i_high)/2;
        if (tri(i) > k)
        {
            i_high = i;
        } else {
            i_low = i;
        }
    }
    return i_low;
}

// implementation 2: quadratic (solve n(n+1) = 2k for n)
// (this looks simpler but is actually slower)
__host__ __device__
int lower_tri2(int k)
{
    float n = (sqrtf(1.0 + 8*k) - 1)/2.0;
    return (int)n;
}

// map from k (the linear index in an upper-triangular DProfile) to i, j coordinates in an n*n matrix
// NB this is NOT the inverse of pindex - but uses a different numbering!
// choose implementation
#define lower_tri lower_tri1
__device__
void getij(int n, int k, int &i, int &j)
{
    i = lower_tri(k, n) + 1;
    j = k - tri(i-1);
}

// map from i, j coordinates to the index in a DProfile
// ( j is along the top: in the upper triangle i <= j)
// Illustrated for n=5:
//  j-->
//     0  1  2  3  4
// i 0 0  1  2  3  4
// | 1    5  6  7  8
// \/2       9 10 11
//   3         12 13
//   5            14
//
// NB this is the same as the order in which elements in an HMatrix h are stored, where
// Allele i = h.m_pmf.first.first
// Allele j = h.m_pmf.first.second

__host__ __device__
int pindex(int i, int j, int n)
{
    return  j + n*i - tri(i);
}

// the folded (upper triangular) index
__host__ __device__
int pindex_ut(int i, int j, int n)
{
    if (i <= j)
    {
        return pindex(i, j, n);
    }
    else
    {
        return pindex(j, i, n);
    }
}

__device__
void storeNResult(NResults *results, int i, float lr)
{
	// first atomically increment the results index
	int n = atomicAdd(&(results->count), 1); // returns old value

	if (n < nresult_max)
	{
		// store results
		results->index[n] = i;
		results->lr[n] = lr;
	}
}

__device__
void storeN2Result(N2Results *results, int i1, int i2, float lr)
{
	// first atomically increment the results index
	int n = atomicAdd(&(results->count), 1);  // returns old value

	if (n < n2result_max)
	{
		// store results
		results->index1[n] = i1;
		results->index2[n] = i2;
		results->lr[n] = lr;
	}
}

// compute element-wise locus1 * locus2 / background and sum
__global__
void cuda_match1(float *loc1, float *loc2, float *back, int n, float *result)
{
// will work in emulator only!
//  printf("blockIdx.x = %d blockDim.x = %d threadIdx.x = %d\n",
//          blockIdx.x, blockDim.x, threadIdx.x);

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<n)
  {
	  // to sum this on the device we need to do a reduction (see SDK reduction example)
	  result[idx] = loc1[idx] * loc2[idx] / back[idx];
  }
}

// Calculate LR for two profiles
// Assume all loci are present in array.
// Values <1 signify locus is not present in the array and should be ignored
__device__
float pp_match(DProfile prof1, DProfile prof2, DProfile back, DProfile spmc1 = DProfile())
{
	float profile_lr = 1;

	for (int i=0; i< cuda_num_loci; ++i)
	{
		int locus_size = loc_info.locus_size[i];

		if (locus_size == 0)
		{
			continue; // locus not present in the population database. Skip it.
		}

		int offset = loc_info.locus_offset[i];
		int n_alleles = loc_info.num_alleles[i];

		if ( (prof1.data[offset] < 0) || (prof2.data[offset] < 0) || (back.data[offset] < 0))
		{
			// locus not present in one or both profiles, or the population database. Skip it.
		}
		else
		{

			float locus_lr = 0;

			for (int j=0; j<locus_size; ++j)
			{
				if (back.data[offset + j] > 0)
				{
					locus_lr += prof1.data[offset + j] * prof2.data[offset + j] / back.data[offset + j];
				}
			}

			// Subpopulation correction
			if (spmc1.data != 0)
			{

				float f = 0; // correction factor

				int  p = 0, q = 0; // indices into Dprofiles
				for (int j=0; j<locus_size; ++j)
				{
					// Because we have half-matrices we must divide the off-diagonal elements by 2
					f += spmc1.data[offset + j] * prof2.data[offset + j] / (2 - d(p,q));

					if (q == n_alleles-1)
					{
						++p; q = p;
					}
					else
					{
						++q;
					}
				}

				if (f>0) locus_lr /= f;
			}

			profile_lr *= locus_lr;
		}
	}

	return profile_lr;
}

// this kernel just does one profile/profile match
__global__ void cuda_match2(DProfile prof1, DProfile prof2, DProfile back, float *result)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx>0) return;

	*result = pp_match(prof1, prof2, back);
}

// This kernel compares prof1 with one of the profiles in prof_db
__global__
void cuda_match3(
	DProfile *prof_db,        // profile database
	int n,                    // size of profile database
	DProfile prof1,           // test profile
	DProfile back,            // background
	NResults *results,
	float lr_threshold)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx<n)
    {
    	// compare profile idx with prof1

    	DProfile prof2 = prof_db[idx];

		float profile_lr = pp_match(prof1, prof2, back);

		// NB fixed size results table
		// A better solution here is a hash table. We still need atomic writes.
		if (profile_lr > lr_threshold)
		{
			storeNResult(results, idx, profile_lr);
		}
    }
}

// This kernel compares prof1 with ALL of the profiles in prof_db
__global__ void cuda_match4(DProfile *prof_db,        // profile database
		                    int n,                    // size of profile database
		                    DProfile back,            // background
		                    N2Results *results,
		                    float lr_threshold)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx<n)
    {
    	// compare profile idx with all other profiles in database

    	DProfile prof2 = prof_db[idx]; // copy into local memory?

    	for (int j=idx+1; j<n; ++j)    // some threads do lots more work than others!
    	{
        	DProfile prof1 = prof_db[j];

        	float profile_lr = pp_match(prof1, prof2, back);

			// NB fixed size results table
			// A better solution here is a hash table. We still need atomic writes.
			if (profile_lr > lr_threshold) // or suitable threshold
			{
				storeN2Result(results, idx, j, profile_lr);
			}
    	}
    }
}

// This kernel can be called in full or half addressing (in full mode half the threads do nothing)
//#define CUDA_MATCH5_HALF

// This kernel compares ONE of the profiles in prof_db with ONE other
__global__ void cuda_match5(DProfile *prof_db,        // profile database
		                    int n,                    // size of profile database
		                    DProfile back,            // background
		                    N2Results *results,
		                    float lr_threshold)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;

	int i = idx;
	int j = idy;

#ifdef CUDA_MATCH5_HALF
	// In "full-size" mode the kernel is called for an array of size n*n.
	// However only n(n-1)/2 matches are needed, so over half the threads do nothing.

	// In "half-size" mode the kernel is called for an array of size (n/2) * (n-1).
	// Then it is necessary to translate the coordinates of this array into the upper-triangular part of the actual n*n array

	// in fact this is slower! - probably because entire warps do nothing and finish quickly, so waste less
	// time than calling getij

	int k = (n-1) * idx + idy;
	getij(n, k, i, j);
#endif

#ifdef __DEVICE_EMULATION__

	  printf("cuda_match5 kernel: nprofiles = %d\n", n);

	  printf("blockIdx.x = %d blockDim.x = %d threadIdx.x = %d\n",
	          blockIdx.x, blockDim.x, threadIdx.x);

	  printf("blockIdx.y = %d blockDim.y = %d threadIdx.y = %d\n",
	          blockIdx.y, blockDim.y, threadIdx.y);

	  printf("i = %d, j = %d\n", i, j);
#endif

    if (i<n && j<i) // in full-size mode half the threads not used!
    {
    	// re-order to "diagonal traversal" to reduce contention for global memory?
    	// In fact this slows us down!
//    	j = i - (j + 1);

    	// compare profile i with j

		DProfile prof1 = prof_db[j];
    	DProfile prof2 = prof_db[i];

		float profile_lr = pp_match(prof1, prof2, back);

		// NB fixed size results table
		// A better solution here is a hash table. We still need atomic writes.
		if (profile_lr > lr_threshold) // or suitable threshold
		{
			storeN2Result(results, i, j, profile_lr);
		}
    }
}

// This kernel compares ONE of the profiles in prof_db1 with ONE in prof_db2
// This may the same profile, or a derivative of it (such as the corresponding sibling profiles)
// in which case we should perform only the upper-triangular comparison (and n1 == n2)
// Or it may be a different data set, in which case we must do the whole matrix
__global__ void cuda_match6(DProfile *prof_db1,       // profile dataset 1
		                    int n1,                   // size of dataset 1
		                    DProfile *prof_db2,       // profile dataset 2
		                    int n2,                   // size of dataset 2
		                    DProfile back,            // background
		                    ArrayPart part,           //  0 = whole array
		                                              //  1 = upper (excluding diagonal)
		                                              // -1 = lower (including diagonal)
		                    N2Results *results,
		                    float lr_threshold,
		                    int i_start = 0,
		                    int j_start = 0,
		            		DProfile *spmc1_dp = 0)   // SPM correction matrices (for profdb1)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;

	int i = idx;
	int j = idy;

#ifdef __DEVICE_EMULATION__

	  printf("cuda_match6 kernel: n1 = %d n2 = %d\n", n1, n2);

	  printf("blockIdx.x = %d blockDim.x = %d threadIdx.x = %d\n",
	          blockIdx.x, blockDim.x, threadIdx.x);

	  printf("blockIdx.y = %d blockDim.y = %d threadIdx.y = %d\n",
	          blockIdx.y, blockDim.y, threadIdx.y);
#endif

	// NB in upper/lower mode half the threads not used
   	if ( !(i<n1 && j<n2)                 // not in array
         || (part == upper && !(i_start + i < j_start + j))    // not in upper triangle
         || (part == lower &&  (i_start + i < j_start + j)) )  // not in lower triangle
    {
    	return;
    }

    // compare profile i with j

	DProfile prof1 = prof_db1[i];
	DProfile prof2 = prof_db2[j];
	DProfile spmc1 = spmc1_dp ? spmc1_dp[i] : DProfile();
	float profile_lr = pp_match(prof1, prof2, back, spmc1);

	// NB fixed size results table
	// A better solution here is a hash table. We still need atomic writes.
	if (profile_lr > lr_threshold)
	{
		storeN2Result(results, i_start + i, j_start + j, profile_lr);
	}
}

__host__
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
		DProfile 	   *spmc1_dp)
{
	cuda_match6 <<< dim3(nBlocks1, nBlocks2), dim3(blockSize, blockSize), 0, stream >>>
			(profdb1_d, db1_chunk_size, profdb2_d, db2_chunk_size, back, part, results_d, lr_threshold, i, j, spmc1_dp);
}

// like cuda_match6 but tiled. One locus at a time to allow shared memory.
//
// 2D Grid 1D block
//
//#define SHARED
__global__
void cuda_match7(
	DProfile *prof_db1,       // profile dataset 1
	int n1,                   // size of dataset 1
	DProfile *prof_db2,       // profile dataset 2
	int n2,                   // size of dataset 2
	DProfile back,            // background
	bool upper,
	N2Results *results,
	float lr_threshold)
{
	int i = threadIdx.x;
	int idx = blockIdx.x*blockDim.x + i;

#ifdef __DEVICE_EMULATION__

	  printf("\ncuda_match7 kernel: n1 = %d n2 = %d\n", n1, n2);

	  printf("blockIdx.x = %d blockDim.x = %d threadIdx.x = %d\n",
	          blockIdx.x, blockDim.x, threadIdx.x);

	  printf("blockIdx.y = %d blockDim.y = %d threadIdx.y = %d\n",
	          blockIdx.y, blockDim.y, threadIdx.y);

	  printf("i = %d idx = %d\n", i, idx);
#endif


	// do one locus at a time

	// shared memory copies of the current locus
#ifdef SHARED
	__shared__ float locus_p1[cuda_locus_size * blockSize];
	__shared__ float locus_p2[cuda_locus_size * blockSize];
	__shared__ float locus_back[cuda_locus_size * blockSize];
#endif

	// results
	__shared__ float lr[blockSize][blockSize];
	for (int j=0; j<blockSize; ++j)
	{
		lr[i][j] = 1;
	}

	for (int loc=0; loc<cuda_num_loci; ++loc)
	{
		int offset = loc_info.locus_offset[loc];
		int locus_size = loc_info.locus_size[loc];

		if (locus_size == 0)
		{
			// TODO: locus not present in population database. What to do?
		}

#ifdef __DEVICE_EMULATION__
		  printf("locus = %d offset = %d\n", loc, offset);
#endif
		//
		// copy locus loc of each profile into shared memory
		//
#ifdef SHARED
		// This memcpy is very slow! is there a faster way to do it?
#if 0
		memcpy(locus_p1 + cuda_locus_size * i, prof_db1[blockIdx.x*blockDim.x + i].data + offset, cuda_locus_size * sizeof(float));
		memcpy(locus_p2 + cuda_locus_size * i, prof_db2[blockIdx.y*blockDim.x + i].data + offset, cuda_locus_size * sizeof(float));
		memcpy(locus_back + cuda_locus_size * i, back.data + offset, cuda_locus_size * sizeof(float));
#else
		float *d1 = locus_p1 + cuda_locus_size * i;
		float *d2 = locus_p2 + cuda_locus_size * i;
		float *d3 = locus_back + cuda_locus_size * i;
		float *s1 = prof_db1[blockIdx.x*blockDim.x + i].data + offset;
		float *s2 = prof_db2[blockIdx.y*blockDim.x + i].data + offset;
		float *s3 = back.data + offset;

		int c = cuda_locus_size;
		while (c--)
		{
			*d1++ = *s1++;
			*d2++ = *s2++;
			*d3++ = *s3++;
		}
#endif

#endif
		// when everyone has caught up all the data for this locus, in this grid element will have been copied
		__syncthreads();

		//
		// do a row of matches
		//
		for (int j=0; j<blockSize; ++j)
		{
			int idy = blockIdx.y*blockDim.x + j; // NB blockDim.x not blockDim.y

#ifdef __DEVICE_EMULATION__
			printf("j = %d idy = %d\n", j, idy);
#endif

			if (upper && ((n1!=n2) || !(idy<n1 && idx<idy))) // in upper mode half the threads not used
			{
				lr[i][j] = 0;
#ifdef __DEVICE_EMULATION__
			printf("nowt to do\n");
#endif
				continue;
			}

#ifndef SHARED
			//
			// Compare loci in global memory
			//
			DProfile prof1 = prof_db1[idx];
			DProfile prof2 = prof_db2[idy];

			if ( (prof1.data[offset] < 0) || (prof2.data[offset] < 0) )
			{
				// locus not present in one or both profiles
	//			continue;
			}
			else
			{

				float locus_lr = 0;

				for (int k=0; k<locus_size; ++k)
				{
					if (back.data[offset + k] > 0)
					{
						locus_lr += prof1.data[offset + k] * prof2.data[offset + k] / back.data[offset + k];
					}
				}

				lr[i][j] *= locus_lr;
			}
#else
			//
			// Compare loci in shared memory
			//
			float *loc1 = locus_p1 + cuda_locus_size * i;
			float *loc2 = locus_p2 + cuda_locus_size * j;


			if ( (*loc1 < 0) || (*loc2 < 0) )
			{
				// locus not present in one or both profiles
			}
			else
			{
				float locus_lr = 0;

				for (int k=0; k<cuda_locus_size; ++k)
				{
					if (locus_back[k] > 0)
					{
						locus_lr += loc1[k] * loc2[k] / locus_back[k];
					}
				}

				lr[i][j] *= locus_lr;
			}
#endif
		}
	}

	// row complete - report results
	for (int j=0; j<blockSize; ++j)
	{
		int idy = blockIdx.y*blockDim.x + j; // NB blockDim.x not blockDim.y

#ifdef __DEVICE_EMULATION__
	  printf("lr[%d][%d] = %f ", i, j, lr[i][j]);
#endif
		if (lr[i][j] > lr_threshold)
		{
			storeN2Result(results, idx, idy, lr[i][j]);
		}
	}
#ifdef __DEVICE_EMULATION__
	  printf("\n");
#endif
}

// construct the contribution of prof to an n-degree descendant
// (equivalently, the contribution from an n-degree ancestor)
// This is the PMF of the allele contributed by prof.
// NB the other allele will be 'F' (background) unless this
// is a two-path descendant. n == 0 represents infinite dilution,
// and returns 'F'
//
// This amounts to evaluating (1/2^n) A + (1/2^n) B + (1 - 2/(2^n)) F
//
__host__ __device__
void allele_pmf(
	int n,                      // degree of relationship
    int i,                      // allele index
    int j,                      // allele index
    int size,                   // size of PMF
    const float back[],         // background
    float pmf[])                // allele PMF (to construct)
{
	float wa = (n==CUDAMatchType::INF) ? 0 : 1.0/pow(2.0, n); // weight of alleles
	float wf = 1 - 2*wa;                                      // weight of background

	for (int k=0; k<size; ++k)
	{
		pmf[k] = back[k] * wf;
	}

	pmf[i] += wa;
	pmf[j] += wa;
}

// construct the vector a * A + b * B + (1-(a+b)) F
__host__ __device__
void gen_allele_pmf(
	float wi,                   // proportion of 'Ai'
	float wj,                   // proportion of 'Aj'
    int i,                      // allele index
    int j,                      // allele index
    int size,                   // size of PMF
    const float back[],         // background
    float pmf[])                // allele PMF (to construct)
{
	float wf = 1 - wi - wj;     // weight of background

	for (int k=0; k<size; ++k)
	{
		pmf[k] = back[k] * wf;
	}

	pmf[i] += wi;
	pmf[j] += wj;
}

__host__ __device__
void
make_sib(
	DProfile prof,       // profile dataset
    DBackground back,    // background vector
	DProfile psib,       // sibling dataset (to construct)
	CudaLocusInfo const *loc_info)
{
    //
    // construct sib_db from prof_db. Method is as in match.cpp: sib()
    //
	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue;
		}

		int loc_offset = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];

		// if the input profile is -1 (locus not present), then so is the sibling
		if (prof.data[loc_offset] < 0)
		{
			for (int j=0; j<locus_size; ++j)
			{
				psib.data[loc_offset + j] = -1;
			}
			continue;
		}

		int n_alleles = loc_info->num_alleles[loc];

		float Pq[CUDA_MAX_ALLELES]; // Sum_q Pqj
		float Pr[CUDA_MAX_ALLELES]; // Sum_r Pir

		for (int i=0; i<n_alleles; ++i)
		{
			Pq[i] = 0; Pr[i] = 0;
		}

		// Sum rows and columns and zero the result
		// (Here we need only consider non-zero elements of P,
		// i.e. in the upper triangle).
		int k = 0; // index into prof
		for (int i=0; i<n_alleles; ++i)
		{
			for (int j=i; j<n_alleles; ++j)
			{
				float p = prof.data[loc_offset + k];

				Pq[j] += p;
				Pr[i] += p;

				psib.data[loc_offset + k] = 0;
				++k;
			}
		}

		// This needs to be a sum over all (i, j).
		// But remember prof and sib are upper-triangular
		for (int j=0; j<n_alleles; ++j)
		{
			for (int i=0; i<n_alleles; ++i)
			{
                int k = pindex_ut(i, j, n_alleles); // the 'folded' index into prof/sib

				float Sij = 0;

				// b(i)b(j)
				Sij += back.data[back_offset + i] * back.data[back_offset + j];

				// b(i) * Sum_q Pqj
				Sij += back.data[back_offset + i] * Pq[j];

				// b(j) * Sum_r Pir
				Sij += back.data[back_offset + j] * Pr[i];

				// Pij. NB since prof is a half-matrix we must take the value as 0 below the diagonal,
				// or we will double-count
				if (i <= j)
				{
					Sij += prof.data[loc_offset + k];
				}

				Sij *= 0.25;

				psib.data[loc_offset + k] += Sij;
			}
		}
	}
}

// Beta-binomial sampling formula:
// Probability of Allele A given that y out of n alleles observed were A; HW freq P(A) = f
__host__ __device__
float bbsf(int y,
           int n,
           float fst,
           float f)
{
    float ret = (y * fst + (1-fst) * f) / (1 + (n-1) * fst);
    return ret;
}

// probability a new allele is p given i, j observed
__host__ __device__
float
prob_BNp_giv_ij(
		int p,
		int i,
		int j,
        float fp,
        float fst)
{
    // P(p|ij)

    int y = d(p,i) + d(p,j);
    return bbsf(y, 2, fst, fp);
}

// probability p and a new allele is i, j given p, q observed
__host__ __device__
float
prob_BNpX_is_ij(
	int p,
	int q,
	int i,
	int j,
    float fi,
    float fj,
    float fst)
{
    // P(pX == ij) | pq where X is drawn from the background

    if (p==i)
    {
        return prob_BNp_giv_ij(j,p,q,fj,fst); // P(j|pq)
    }
    else if (p==j)
    {
        return prob_BNp_giv_ij(i,p,q,fi,fst); // P(i|pq)
    }
    else
    {
        return 0;
    }
}

// probability two new alleles are p,q given i,j observed
// (NB this is a generalized version of NRC4_10)
__host__ __device__
float
prob_BNpq_giv_ij(
		int p,
		int q,
		int i,
		int j,
        float fp,
        float fq,
        float fst)
{
    // P(pq|ij) = P(p|ij) * P(q|pij) : *2 if p!=q

    int y1 = d(p,i) + d(p,j);
    int y2 = d(q,p) + d(q,i) + d(q,j);

    float x = bbsf(y1, 2, fst, fp) * bbsf(y2, 3, fst, fq);

    if (p != q) x *= 2;

    return x;
}

// P(Rel(pq) = ij)
__host__ __device__
float
rel_BN(
      int p,
      int q,
      int i,
      int j,
      float fi,
      float fj,
      float k0,    // no alleles Identical By Descent
      float k1,    // one allele IBD
      float k2d,   // two alleles (different) IBD
//    float k2s,   // two alleles (the same) IBD // TODO - needed for bilineal
      float fst)
{
    return    k2d   * d(p,i) * d(q,j)                        // P(pq == ij)
           + (k1/2) * (prob_BNpX_is_ij(p,q,i,j,fi,fj,fst) +  // P(p? == ij)
                       prob_BNpX_is_ij(q,p,i,j,fi,fj,fst) )  // P(q? == ij)
           +  k0    *  prob_BNpq_giv_ij(i,j,p,q,fi,fj,fst);  // P(?? == ij) == P(ij|pq)
}

__host__ __device__
void
make_rel_BN(
	DProfile prof,       // profile dataset
    DBackground back,    // background vector
	DProfile psib,       // sibling dataset (to construct)
	CudaLocusInfo const *loc_info,
	float k0,              // no alleles Identical By Descent
	float k1,              // one allele IBD
	float k2d,             // two alleles (different) IBD
//      float k2s,         // two alleles (the same) IBD // TODO - needed for bilineal
	float theta)
{
    //
    // construct sib_db from prof_db.
	//
    // For all (loc):
    //   For all (ij):
	//     psib(ij) = 0
    //
    //   For all (pq): // elements of prof
	//     w = prof(pq)
	//     if (w>0)
    // 	     For all (ij): // elements of psib
	//         fi = back(i)
	//         fj = back(j)
    //         psib(ij) += w * P(Rel(pq) = ij)
	//
	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue; // next locus
		}

		int loc_offset  = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];
		int n_alleles   = loc_info->num_alleles[loc];

		// if the input profile is -1 (locus not present), then so is the sibling
		if (prof.data[loc_offset] < 0)
		{
			for (int k=0; k<locus_size; ++k)
			{
				psib.data[loc_offset + k] = -1;
			}
			continue; // next locus
		}

		// zero psib
		for (int k=0; k<locus_size; ++k)
		{
			psib.data[loc_offset + k] = 0;
		}

		// loop over elements of prof (which is upper-triangular)
		// *in the order of the DProfile index*
		int n = 0; // index into prof
		for (int p=0; p<n_alleles; ++p)
		{
			for (int q=p; q<n_alleles; ++q)
			{
				float w = prof.data[loc_offset + n]; // weight of this element of prof

				if (w > 0)
				{
					// loop over elements of psib *in the order of the DProfile index*
					int k = 0; // index into psib
					for (int i=0; i<n_alleles; ++i)
					{
						for (int j=i; j<n_alleles; ++j)
						{
							float fi = back.data[back_offset + i];
							float fj = back.data[back_offset + j];

							// P( Rel(pq) = ij )
							float x = rel_BN(p, q, i, j, fi, fj, k0, k1, k2d, theta);
							psib.data[loc_offset + k] += w * x;
							++k;
						}
					}
				}
				++n;
			}
		}
	}
}

// construct the (n, m) two-path relative of prof
// Use m=MatchType::INF for a one-path relationship

__host__ __device__
void
make_2path(
	DProfile prof,       // profile dataset
	int n,               // first degree of relationship
	int m,               // second degree of relationship (MatchType::INF == F)
    DBackground back,    // background vector
	DProfile prel,       // relative dataset (to construct)
	CudaLocusInfo const *loc_info)
{
	// At each locus, for each element of prof, construct a degree n and a degree m and multiply them
	// TODO this is the 'long hand' method - on CUDA we should reduce the sum algebraically as we do for sibling and this will be faster
	// We can use this implementation for comparison.

	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int loc_offset = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue;
		}

		// if the input profile is -1 (locus not present), then so is the relative
		if (prof.data[loc_offset] < 0)
		{
			for (int j=0; j<locus_size; ++j)
			{
				prel.data[loc_offset + j] = -1;
			}
			continue;
		}

		int n_alleles = loc_info->num_alleles[loc];

		float Pn[CUDA_MAX_ALLELES]; // n-degree distribution
		float Pm[CUDA_MAX_ALLELES]; // m-degree distribution

		// loop over elements of prof (which is upper-triangular)
		for (int ip=0; ip<n_alleles; ++ip)
		{
			for (int jp=ip; jp<n_alleles; ++jp)
			{
				int p = pindex(ip, jp, n_alleles);
				float wp = prof.data[loc_offset + p]; // weight of this element of prof

				if (wp > 0)
				{
					// construct n-degree and m-degree allele distributions for this element
					allele_pmf(n, ip, jp, n_alleles, back.data + back_offset, Pn);
					allele_pmf(m, ip, jp, n_alleles, back.data + back_offset, Pm);

					// multiply and add in to the result
					for (int j=0; j<n_alleles; ++j)
					{
						for (int i=0; i<n_alleles; ++i)
						{
			                int k = pindex_ut(i, j, n_alleles); // the 'folded' index into prof/prel

							prel.data[loc_offset + k] += Pn[i] * Pm[j] * wp;
						}
					}
				}
			}
		}
	}
}

__host__ __device__
void
make_genrc_elem(
	float wp,            // element weight
	int i,               // element i
	int j,               // element j (>=i)
	int n_alleles,       // alleles at this locus
	float a1,            // proportion of 'a' inherited on side 1
	float b1,            // proportion of 'b' inherited on side 1
	float a2,            // proportion of 'a' inherited on side 2
	float b2,            // proportion of 'b' inherited on side 2
    float *pback,        // background element
    float *prel)         // relative element (to construct)
{
	float P1[CUDA_MAX_ALLELES];   // vector for side 1, as given
	float P2[CUDA_MAX_ALLELES];   // vector for side 2, as given
	float P1_s[CUDA_MAX_ALLELES]; // vector for side 1, a's and b's swapped
	float P2_s[CUDA_MAX_ALLELES]; // vector for side 2, a's and b's swapped

	if (wp > 0)
	{
		wp /= 2; // because we must average two components to form the matrix

		// construct vectors as given
		gen_allele_pmf(a1, b1, i, j, n_alleles, pback, P1);
		gen_allele_pmf(a2, b2, i, j, n_alleles, pback, P2);

		// construct vectors with a's and b's swapped
		gen_allele_pmf(b1, a1, i, j, n_alleles, pback, P1_s);
		gen_allele_pmf(b2, a2, i, j, n_alleles, pback, P2_s);

		// multiply and add in to the result
		for (int p=0; p<n_alleles; ++p)
		{
			for (int q=0; q<n_alleles; ++q)
			{
				int k = pindex_ut(p, q, n_alleles); // the 'folded' index into prof/prel

				prel[k] += (P1[p] * P2[q] + P1_s[p] * P2_s[q]) * wp;
			}
		}
	}
}

// construct the [(a1, b1), (a2, b2)] generalized relationship coefficient relative of prof
//
// NB since we do not know which parent an allele comes from, we must always
// 1) form a matrix with the numbers as given
// 2) form a matrix with a1/b1 swapped and a2/b2 swapped
// 3) average them
//
// NB the matrices are symmetrical (upper triangular) which takes care of the side-1/side-2 symmetry
//
__host__ __device__
void
make_genrc(
	DProfile prof,       // profile dataset
	float a1,            // proportion of 'a' inherited on side 1
	float b1,            // proportion of 'b' inherited on side 1
	float a2,            // proportion of 'a' inherited on side 2
	float b2,            // proportion of 'b' inherited on side 2
    DBackground back,    // background vector
	DProfile prel,       // relative dataset (to construct)
	CudaLocusInfo const *loc_info)
{
	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int loc_offset = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue;
		}

		// if the input profile is -1 (locus not present), then so is the relative
		if (prof.data[loc_offset] < 0)
		{
			for (int j=0; j<locus_size; ++j)
			{
				prel.data[loc_offset + j] = -1;
			}
			continue;
		}

		int n_alleles = loc_info->num_alleles[loc];

		// loop over elements of prof (which is upper-triangular)
		for (int ip=0; ip<n_alleles; ++ip)
		{
			for (int jp=ip; jp<n_alleles; ++jp)
			{
				int p = pindex(ip, jp, n_alleles);
				float wp = prof.data[loc_offset + p]; // weight of this element of prof

				make_genrc_elem(wp, ip, jp, n_alleles, a1, b1, a2, b2, back.data + back_offset, prel.data + loc_offset);
			}
		}
	}
}

//
// make the inverse of the given 4-number relationship
// (Bayes theorem - see maths paper)
//
__host__ __device__
void
make_invrc(
	DProfile prof,       // profile dataset
	float a1,            // proportion of 'a' inherited on side 1
	float b1,            // proportion of 'b' inherited on side 1
	float a2,            // proportion of 'a' inherited on side 2
	float b2,            // proportion of 'b' inherited on side 2
    DBackground back,    // background vector
	DProfile prel,       // relative dataset (to construct)
	CudaLocusInfo const *loc_info)
{
	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int loc_offset = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue;
		}

		// if the input profile is -1 (locus not present), then so is the relative
		if (prof.data[loc_offset] < 0)
		{
			for (int j=0; j<locus_size; ++j)
			{
				prel.data[loc_offset + j] = -1;
			}
			continue;
		}

		int n_alleles = loc_info->num_alleles[loc];

		for (int i=0; i<n_alleles; ++i)
		{
			for (int j=i; j<n_alleles; ++j)
			{
				const int size = CUDA_MAX_ALLELES * (CUDA_MAX_ALLELES + 1) / 2; // auto arrays must be const size
				float R[size] = { 0 }; // should initialize whole array to zero ...

				// ... but doesn't, so
				int actual_size = n_alleles * (n_alleles + 1) / 2;
				for (int k=0; k<actual_size; ++k)
				{
				    R[k] = 0;
				}

				int k_ij = pindex(i, j, n_alleles);

				float b_ij = (back.data + back_offset)[i] * (back.data + back_offset)[j]
														  * (i==j ? 1 : 2);

				const float AiAj = 1.0;
				make_genrc_elem(AiAj, i, j, n_alleles, a1, b1, a2, b2, back.data + back_offset, R);

				float prel_elem = 0;

				for (int p=0; p<n_alleles; ++p)
				{
					for (int q=p; q<n_alleles; ++q)
					{
						int k_pq = pindex(p, q, n_alleles);

						float b_pq = (back.data + back_offset)[p] * (back.data + back_offset)[q]
																  * (p==q ? 1 : 2);

						prel_elem += prof.data[loc_offset + k_pq] * R[k_pq] / b_pq;
					}
				}

				prel.data[loc_offset + k_ij] = prel_elem * b_ij;
			}
		}
	}
}

// Construct Generalized Relationship Coefficient relative profiles (non-sibling)
__global__ void cuda_genrc(DProfile *prof_db,      // profile dataset
						  DProfile *rel_db,        // relative dataset (to construct)
		                  int n,                   // size of dataset
 		             	  float a1,                // proportion of 'a' inherited on side 1
		             	  float b1,            	   // proportion of 'b' inherited on side 1
		             	  float a2,                // proportion of 'a' inherited on side 2
		             	  float b2,                // proportion of 'b' inherited on side 2
		                  DBackground back,        // background
		                  bool inverse = false)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i>=n)
    {
    	return;
    }

    // initialize rel_db[i] to zero
	for (int j=0; j<prof_db[i].size; ++j)
	{
		rel_db[i].data[j] = 0;
	}

    // construct rel_db[i] from prof_db[i]
	if (inverse)
	{
		make_invrc(prof_db[i], a1, b1, a2, b2, back, rel_db[i], &loc_info);
	}
	else
	{
		make_genrc(prof_db[i], a1, b1, a2, b2, back, rel_db[i], &loc_info);
	}
}

// Construct (d1, d2) relative profiles (non-sibling)
__global__ void cuda_rel(DProfile *prof_db,       // profile dataset
						 DProfile *rel_db,        // relative dataset (to construct)
		                 int n,                   // size of dataset
		                 int d1,                  // degree first path
		                 int d2,                  // degree second path (NOT HANDLED)
		                 DBackground back,        // background
		                 CudaSubPopModel spm)
//		                 float theta)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i>=n)
    {
    	return;
    }

    if (spm.type == CudaSubPopModel::HW)
    {
        // initialize rel_db[i] to zero
        // NB memset needs a constant size. We could use cuda_prof_size but watch out if this is no longer constant!
    //	memset(rel_db[i].data, 0, cuda_prof_size * sizeof(float));
    	for (int j=0; j<prof_db[i].size; ++j)
    	{
    		rel_db[i].data[j] = 0;
    	}

        // construct rel_db[i] from prof_db[i].
        make_2path(prof_db[i], d1, d2, back, rel_db[i], &loc_info);

    }
    else // Balding-Nichols (we do not support 4.4 for relatives)
    {
    	// IBD coefficients for 1-path:
    	float k1 = pow(2.0, 1-d1);
    	float k0 = 1 - k1;
    	float k2d = 0;

    	make_rel_BN(prof_db[i], back, rel_db[i], &loc_info, k0, k1, k2d, spm.theta_bar);
    }
}

__device__
void
make_spmc(
	DProfile prof,       // profile dataset
	DProfile back,       // background matrix (4.4)
    DBackground backv,   // background vector (HW)
	DProfile spmc,       // SPM correction matrices (to construct)
	float theta,         //  Fst
	CudaLocusInfo const *loc_info)
{
    // for each locus, construct a correction matrix
	for (int loc=0; loc< cuda_num_loci; ++loc)
	{
		int locus_size = loc_info->locus_size[loc];

		if (locus_size == 0) // locus not in population database
		{
			continue;
		}

		int loc_offset = loc_info->locus_offset[loc];
		int back_offset = loc_info->back_offset[loc];

		// if the input profile is -1 (locus not present), then so is the correction matrix
		if (prof.data[loc_offset] < 0)
		{
			for (int j=0; j<locus_size; ++j)
			{
				spmc.data[loc_offset + j] = -1;
			}
			continue;
		}

		int n_alleles = loc_info->num_alleles[loc];

		// loop over all genotypes
		int k = 0; // index into prof, spmc etc
		for (int i=0; i<n_alleles; ++i)
		{
			for (int j=i; j<n_alleles; ++j)
			{
				float x = 0;
				int m = 0;
				for (int p=0; p<n_alleles; ++p)
				{
					for (int q=p; q<n_alleles; ++q)
					{
				        float fp = backv.data[back_offset + p]; // HW frequency of p
				        float fq = backv.data[back_offset + q]; // HW frequency of q

						float Bpq_ij = prob_BNpq_giv_ij(p, q, i, j, fp, fq, theta);

						x += prof.data[loc_offset + m] * Bpq_ij / back.data[loc_offset + m];
						++m;
					}
				}

				spmc.data[loc_offset + k] = (i==j) ? x : 2*x; // *2 because a half-matrix
				++k;
			}
		}
	}
}

// Construct correction matrices
__global__ void cuda_spmc(DProfile *prof_db,       // profile dataset
						  DProfile *spmc_db,       // SPM correction matrix (to construct)
		                  int n,                   // size of dataset
		              	  DProfile back,           // background matrix (4.4)
		                  DBackground backv,       // background vector (HW)
		                  float theta)             // Fst
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i>=n)
    {
    	return;
    }

    make_spmc(prof_db[i], back, backv, spmc_db[i], theta, &loc_info);
}

// Construct sibling profiles
__global__ void cuda_sib(DProfile *prof_db,       // profile dataset
						 DProfile *sib_db,        // sibling dataset (to construct)
		                 int n,                   // size of dataset
		                 DBackground back,        // background
		                 CudaSubPopModel spm)
//		                 float theta)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i>=n)
    {
    	return;
    }

    // construct sib_db[i] from prof_db[i].
    if (spm.type == CudaSubPopModel::HW)
    {
    	make_sib(prof_db[i], back, sib_db[i], &loc_info);
    }
    else // Balding-Nichols (we do not support 4.4 for relatives)
    {
    	// IBD coefficients for sibling:
    	float k0 = 0.25;
    	float k1 = 0.5;
    	float k2d = 0.25;

    	make_rel_BN(prof_db[i], back, sib_db[i], &loc_info, k0, k1, k2d, spm.theta_bar);
    }

}

#ifndef STREAMED
void cuda_n2match(std::vector<DProfile> &prof_db, // Profile database
				  DProfile back,                  // Background
				  CudaLocusInfo const &locus_info,
				  N2Results &results,
				  float lr_threshold)
{
	int n_floats_per_profile = back.size;
	int n_profiles = prof_db.size();
	int n_floats = n_floats_per_profile * n_profiles;
	Assert(n_floats_per_profile == locus_info.profile_size);

	//
	// allocate DProfiles data and result on device
	//
	Timer t;

	// the result
	N2Results *results_d;
	Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
			"cuda_n2match: A: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: B: cudaThreadSynchronize failed");
	Assert(results_d);

	// the data in the profiles
	float *profdb_data = 0;
	Assert2(cudaMalloc((void **) &profdb_data, sizeof(float)*n_floats) != cudaErrorMemoryAllocation,
			"cuda_n2match: C: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: D: cudaThreadSynchronize failed");
	Assert(profdb_data);

	// the DProfiles themselves (containing pointers to the data)
	DProfile *profdb_d = 0;
	Assert2(cudaMalloc((void **) &profdb_d, sizeof(DProfile)*n_profiles) != cudaErrorMemoryAllocation,
			"cuda_n2match: E: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: F: cudaThreadSynchronize failed");
	Assert(profdb_d);

	// the data for the background
	float *back_data = 0;
	Assert2(cudaMalloc((void **) &back_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
			"cuda_n2match: G: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: H: cudaThreadSynchronize failed");
	Assert(back_data);

    t.stop();
    info << startl << "cuda_n2match(): Allocating memory on device took " << t << " seconds" << std::endl;

	//
	// copy to device once only
	//
    t.start();

    // copy offsets
    COPY_OFFSETS();

	// zero the result
	results.count = 0;
	cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// Data for each DProfile. After copying each data array, copy the device address into the DProfile.
	for (int i=0; i<n_profiles; ++i)
	{
		float *addr = profdb_data + (i * n_floats_per_profile);
		cudaMemcpy(addr, prof_db[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
		prof_db[i].data = addr;
	}

	// the DProfiles themselves (containing pointers to the data)
	cudaMemcpy(profdb_d, &(prof_db[0]), sizeof(DProfile)*n_profiles, cudaMemcpyHostToDevice);

	// the data for the background
	cudaMemcpy(back_data,   back.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
	back.data  = back_data;

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: I: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_n2match(): Copying data to device took " << t << " seconds" << std::endl;

    //
	// Call match kernel
    //
    t.start();

    info << startl << "launching kernel cuda_match6" << std::endl;
    int nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;

#if 1
    cuda_match6 <<< dim3(nBlocks, nBlocks), dim3(blockSize, blockSize) >>>
    		(profdb_d, n_profiles, profdb_d, n_profiles, back, upper, results_d, lr_threshold);
#else
    cuda_match7 <<< dim3(nBlocks, nBlocks), dim3(blockSize) >>>
    		(profdb_d, n_profiles, profdb_d, n_profiles, back, true, results_d, lr_threshold);
#endif

    cudaError_t cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    cudaError_t err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_n2match: J: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_n2match(): Kernel cuda_match6 took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: K: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb_data);
	cudaFree(profdb_d);
	cudaFree(back_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2match: L: cudaThreadSynchronize failed");
}

// match one db against another (not necessarily the same size)
void cuda_nm_match(std::vector<DProfile> &prof_db1, // Profile database
				   std::vector<DProfile> &prof_db2, // Profile database
				   DProfile back,                   // Background
 	 			   CudaLocusInfo const &locus_info,
				   N2Results &results,
				   float lr_threshold,
				   ArrayPart part)
{
	int n_floats_per_profile = back.size;
	int n1 = prof_db1.size();
	int n2 = prof_db2.size();
	int n_floats1 = n_floats_per_profile * n1;
	int n_floats2 = n_floats_per_profile * n2;
	Assert(n_floats_per_profile == locus_info.profile_size);

//	Assert(part == full || n1 == n2);

	//
	// allocate DProfiles data and result on device
	//
	Timer t;

	// the result
	N2Results *results_d;
	Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
			"cuda_nm_match: AA: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: A: cudaThreadSynchronize failed");
	Assert(results_d);

	// the data in the profiles
	float *profdb1_data = 0;
	Assert2(cudaMalloc((void **) &profdb1_data, sizeof(float)*n_floats1) != cudaErrorMemoryAllocation,
			"cuda_nm_match: B: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: C: cudaThreadSynchronize failed");
	Assert(profdb1_data);

	float *profdb2_data = 0;
	Assert2(cudaMalloc((void **) &profdb2_data, sizeof(float)*n_floats2) != cudaErrorMemoryAllocation,
			"cuda_nm_match: D: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: E: cudaThreadSynchronize failed");
	Assert(profdb2_data);

	// the DProfiles themselves (containing pointers to the data)
	DProfile *profdb1_d = 0;
	Assert2(cudaMalloc((void **) &profdb1_d, sizeof(DProfile)*n1) != cudaErrorMemoryAllocation,
			"cuda_nm_match: F: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: G: cudaThreadSynchronize failed");
	Assert(profdb1_d);

	DProfile *profdb2_d = 0;
	Assert2(cudaMalloc((void **) &profdb2_d, sizeof(DProfile)*n2) != cudaErrorMemoryAllocation,
			"cuda_nm_match: H: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: I: cudaThreadSynchronize failed");
	Assert(profdb2_d);

	// the data for the background
	float *back_data = 0;
	Assert2(cudaMalloc((void **) &back_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
			"cuda_nm_match: J: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: K: cudaThreadSynchronize failed");
	Assert(back_data);

    t.stop();
    info << startl << "cuda_nm_match(): Allocating memory on device took " << t << " seconds" << std::endl;

	//
	// copy to device
	//
    t.start();

    // copy offsets
    COPY_OFFSETS();

	// zero the result
	results.count = 0;
	cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// Data for each DProfile. After copying each data array, copy the device address into the DProfile.
	for (int i=0; i<n1; ++i)
	{
		float *addr = profdb1_data + (i * n_floats_per_profile);
		cudaMemcpy(addr, prof_db1[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
		prof_db1[i].data = addr;
	}

	for (int i=0; i<n2; ++i)
	{
		float *addr = profdb2_data + (i * n_floats_per_profile);
		cudaMemcpy(addr, prof_db2[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
		prof_db2[i].data = addr;
	}

	// the DProfiles themselves (containing pointers to the data)
	cudaMemcpy(profdb1_d, &(prof_db1[0]), sizeof(DProfile)*n1, cudaMemcpyHostToDevice);
	cudaMemcpy(profdb2_d, &(prof_db2[0]), sizeof(DProfile)*n2, cudaMemcpyHostToDevice);

	// the data for the background
	cudaMemcpy(back_data,   back.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
	back.data  = back_data;

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: L: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_nm_match(): Copying data to device took " << t << " seconds" << std::endl;

    //
	// Call solution kernel
    //
    t.start();

    info << startl << "launching kernel cuda_match6" << std::endl;
    int nBlocks1 = n1/blockSize + (n1%blockSize == 0?0:1);
    int nBlocks2 = n2/blockSize + (n2%blockSize == 0?0:1);
    info << alignl << "nBlocks1 = " << nBlocks1 << std::endl;
    info << alignl << "nBlocks2 = " << nBlocks2 << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n1 = " << n1 << std::endl;
    info << alignl << "n2 = " << n2 << std::endl;
    cuda_match6 <<< dim3(nBlocks1, nBlocks2), dim3(blockSize, blockSize) >>>
    		(profdb1_d, n1, profdb2_d, n2, back, part, results_d, lr_threshold);

    cudaError_t cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    cudaError_t err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_nm_match: M: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_nm_match(): Kernel cuda_match6 took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: N: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb1_data);
	cudaFree(profdb2_data);
	cudaFree(profdb1_d);
	cudaFree(profdb2_d);
	cudaFree(back_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_match: P: cudaThreadSynchronize failed");
}

// match REL(prof_db1) against prof_db2 (not necessarily the same size)
void cuda_nm_relmatch(
	std::vector<DProfile> &prof_db1,   // Profile database
	std::vector<DProfile> &prof_db2,   // Profile database
	CUDAMatchType   const &match_type, // Relative type
	DProfile               backh,      // Background as a half-matrices
    CudaLocusInfo   const &locus_info,
    DBackground            backv,      // Background as vectors
	N2Results &results,
	float lr_threshold,
	ArrayPart part)
{
	int n_floats_per_profile = backh.size;
	int n1 = prof_db1.size();
	int n2 = prof_db2.size();
	int n_floats1 = n_floats_per_profile * n1;
	int n_floats2 = n_floats_per_profile * n2;
	Assert(n_floats_per_profile == locus_info.profile_size);

//	Assert(part == full || n1 == n2);

	//
	// allocate DProfiles data and result on device
	//
    Timer t;

	// the result
	N2Results *results_d;
	Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: AA: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: A: cudaThreadSynchronize failed");
	Assert(results_d);

    t.stop();
    info << startl << "cuda_nm_relmatch(): cudaMalloc at 1680 took " << t << " seconds" << std::endl;

	// the data in the profiles
	// NB we will start by using db2 as temporary storage for db1, so create it at size n1
	float *profdb1_data = 0;
	Assert2(cudaMalloc((void **) &profdb1_data, sizeof(float)*n_floats1) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: B: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: C: cudaThreadSynchronize failed");
	Assert(profdb1_data);

	float *profdb2_data = 0;
	Assert2(cudaMalloc((void **) &profdb2_data, sizeof(float)*n_floats1) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: D: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: E: cudaThreadSynchronize failed");
	Assert(profdb2_data);

	// the DProfiles themselves (containing pointers to the data)
	DProfile *profdb1_d = 0;
	Assert2(cudaMalloc((void **) &profdb1_d, sizeof(DProfile)*n1) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: F: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: G: cudaThreadSynchronize failed");
	Assert(profdb1_d);

	DProfile *profdb2_d = 0;
	Assert2(cudaMalloc((void **) &profdb2_d, sizeof(DProfile)*n1) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: H: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: I: cudaThreadSynchronize failed");
	Assert(profdb2_d);

	// the data for the background
	float *backh_data = 0;
	Assert2(cudaMalloc((void **) &backh_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: K: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: L: cudaThreadSynchronize failed 6");
	Assert(backh_data);

	float *backv_data = 0;
	Assert2(cudaMalloc((void **) &backv_data, sizeof(float)*locus_info.back_size) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: M: cudaMalloc failed");

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: N: cudaThreadSynchronize failed 7");
	Assert(backv_data);

    t.stop();
    info << startl << "cuda_nm_relmatch(): Allocating memory on device took " << t << " seconds" << std::endl;

	//
	// copy to device
	//
    t.start();

    // copy offsets
    COPY_OFFSETS();

	// zero the result
	results.count = 0;
	cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// copy db1 data into db2 on the device
	for (int i=0; i<n1; ++i)
	{
		float *addr = profdb2_data + (i * n_floats_per_profile);

		cudaMemcpy(addr, prof_db1[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);

		// after copying, overwrite the data pointer with the device address
		prof_db1[i].data = addr;
	}

	// copy the DProfiles themselves (containing pointers to the data) for db2
	cudaMemcpy(profdb2_d, &(prof_db1[0]), sizeof(DProfile)*n1, cudaMemcpyHostToDevice);

	// create array of DProfiles for db1 on the device
	for (int i=0; i<n1; ++i)
	{
		float *addr = profdb1_data + (i * n_floats_per_profile);

		// copy no data into db1 - this is where the relatives profiles will go

		// re-use the db1 data pointers
		prof_db1[i].data = addr;
	}

	// copy the DProfiles themselves (containing pointers to the data) for db1
	cudaMemcpy(profdb1_d, &(prof_db1[0]), sizeof(DProfile)*n1, cudaMemcpyHostToDevice);

	// the data for the background DProfile
	cudaMemcpy(backh_data,   backh.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
	backh.data  = backh_data;

	// the data for the background DBackground
	cudaMemcpy(backv_data,   backv.data,  sizeof(float)*locus_info.back_size, cudaMemcpyHostToDevice);
	backv.data  = backv_data;

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: L: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_nm_relmatch(): Copying first dataset to device took " << t << " seconds" << std::endl;

	// call kernel to calculate relative profile of first db. (This now goes in the db1 position on the device)
    t.start();
    if (match_type.m_rel_type == sibling_t)
    {
    	info << startl << "launching kernel cuda_sib" << std::endl;
    }
    else
    {
    	info << startl << "launching kernel cuda_rel" << std::endl;
    }

    int nBlocks = n1/blockSize + (n1%blockSize == 0?0:1);
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n1 = " << n1 << std::endl;

    Assert2(match_type.m_rel_type != ident_t, "cuda_nm_relmatch: P2: called with m_rel_type == ident_t");

    if (match_type.m_rel_type == sibling_t)
    {
		cuda_sib <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb2_d, profdb1_d, n1, backv);
    }
    else if (match_type.m_rel_type == gen_t)
    {
		cuda_genrc <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb2_d, profdb1_d, n1, match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, backv);
    }
    else if (match_type.m_rel_type == inv_t)
    {
    	bool inverse = true;
		cuda_genrc <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb2_d, profdb1_d, n1, match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, backv, inverse);
    }
    else
    {
		cuda_rel <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb2_d, profdb1_d, n1, match_type.m_path1steps, match_type.m_path2steps, backv);
    }

    cudaError_t cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    cudaError_t err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_nm_relmatch: Q: cudaThreadSynchronize failed 9");

    t.stop();
    info << startl << "Kernel " << ((match_type.m_rel_type == sibling_t) ? "cuda_sib" : "cuda_rel") << " took " << t << " seconds" << std::endl;

    // deallocate db2 on the device and reallocate it the right size
    t.start();

	cudaFree(profdb2_data);
	cudaFree(profdb2_d);

    profdb2_data = 0;
	Assert2(cudaMalloc((void **) &profdb2_data, sizeof(float)*n_floats2) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: D: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: E2: cudaThreadSynchronize failed");
	Assert(profdb2_data);

	profdb2_d = 0;
	Assert2(cudaMalloc((void **) &profdb2_d, sizeof(DProfile)*n2) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: H: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: I2: cudaThreadSynchronize failed");
	Assert(profdb2_d);

    t.stop();
    info << startl << "cuda_nm_relmatch(): Allocating memory on device took " << t << " seconds" << std::endl;

    // copy second database to db2 on the device
	t.start();

	for (int i=0; i<n2; ++i)
	{
		float *device_addr = profdb2_data + (i * n_floats_per_profile);
		cudaMemcpy(device_addr, prof_db2[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
		prof_db2[i].data = device_addr;
	}

	// the DProfiles themselves (containing pointers to the data)
	cudaMemcpy(profdb2_d, &(prof_db2[0]), sizeof(DProfile)*n2, cudaMemcpyHostToDevice);

    t.stop();
    info << startl << "cuda_nm_relmatch(): Copying second dataset to device took " << t << " seconds" << std::endl;

    // call solution kernel
    t.start();

    info << startl << "launching kernel cuda_match6" << std::endl;
    int nBlocks1 = n1/blockSize + (n1%blockSize == 0?0:1);
    int nBlocks2 = n2/blockSize + (n2%blockSize == 0?0:1);
    info << alignl << "nBlocks1 = " << nBlocks1 << std::endl;
    info << alignl << "nBlocks2 = " << nBlocks2 << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n1 = " << n1 << std::endl;
    info << alignl << "n2 = " << n2 << std::endl;
    cuda_match6 <<< dim3(nBlocks1, nBlocks2), dim3(blockSize, blockSize) >>>
    		(profdb1_d, n1, profdb2_d, n2, backh, part, results_d, lr_threshold);

    cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_nm_relmatch: M: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_nm_relmatch(): Kernel cuda_match6 took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: N: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

    // clean up
	cudaFree(results_d);
	cudaFree(profdb1_data);
	cudaFree(profdb2_data);
	cudaFree(profdb1_d);
	cudaFree(profdb2_d);
	cudaFree(backh_data);
	cudaFree(backv_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: P: cudaThreadSynchronize failed");
}

// load data and allocate equal space for the relative data
// call kernel to fill relative data
// call match kernel, upper mode, point to the separate arrays
void cuda_n2relmatch(
	std::vector<DProfile> &prof_db,    // Profile database
	CUDAMatchType   const &match_type, // Relative type
	DProfile               backh,      // Background as a half-matrices
    CudaLocusInfo   const &locus_info,
    DBackground            backv,      // Background as vectors
	N2Results             &results,
	float                  lr_threshold)
{
	int n_floats_per_profile = backh.size;
	int n_profiles = prof_db.size();
	int n_floats = n_floats_per_profile * n_profiles;
	Assert(n_floats_per_profile == locus_info.profile_size);
	Assert(backv.size == locus_info.back_size);
	//
	// allocate DProfiles data and result on device
	//
    Timer t;

	// the result
	N2Results *results_d;
	Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: A: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: B: cudaThreadSynchronize failed 1");
	Assert(results_d);

	// the data in the profiles (two copies: originals and relatives)
	float *profdb1_data = 0;
	Assert2(cudaMalloc((void **) &profdb1_data, sizeof(float)*n_floats) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: C: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: D: cudaThreadSynchronize failed 2");
	Assert(profdb1_data);

	float *profdb2_data = 0;
	Assert2(cudaMalloc((void **) &profdb2_data, sizeof(float)*n_floats) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: E: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: F: cudaThreadSynchronize failed 3");
	Assert(profdb2_data);

	// the DProfiles themselves, containing pointers to the data (two, ditto)
	DProfile *profdb1_d = 0;
	Assert2(cudaMalloc((void **) &profdb1_d, sizeof(DProfile)*n_profiles) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: G: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: H: cudaThreadSynchronize failed 4");
	Assert(profdb1_d);

	DProfile *profdb2_d = 0;
	Assert2(cudaMalloc((void **) &profdb2_d, sizeof(DProfile)*n_profiles) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: I: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: J: cudaThreadSynchronize failed 5");
	Assert(profdb2_d);

	// the data for the background
	float *backh_data = 0;
	Assert2(cudaMalloc((void **) &backh_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: K: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: L: cudaThreadSynchronize failed 6");
	Assert(backh_data);

	float *backv_data = 0;
	Assert2(cudaMalloc((void **) &backv_data, sizeof(float)*locus_info.back_size) != cudaErrorMemoryAllocation,
			"cuda_n2relmatch: M: cudaMalloc failed");

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: N: cudaThreadSynchronize failed 7");
	Assert(backv_data);

    t.stop();
    info << startl << "cuda_n2relmatch(): Allocating memory on device took " << t << " seconds" << std::endl;

	//
	// copy to device (copy original, allocate space for relative)
	//
    t.start();

    // copy offsets
    COPY_OFFSETS();

	// zero the result
	results.count = 0;
	cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// Data for each DProfile. After copying each data array, copy the device address into the DProfile.
	for (int i=0; i<n_profiles; ++i)
	{
		float *addr = profdb1_data + (i * n_floats_per_profile);
		cudaMemcpy(addr, prof_db[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
		prof_db[i].data = addr;
	}

	// the DProfiles themselves (containing pointers to the data)
	cudaMemcpy(profdb1_d, &(prof_db[0]), sizeof(DProfile)*n_profiles, cudaMemcpyHostToDevice);

	// construct and copy over addresses for the sibling database
	for (int i=0; i<n_profiles; ++i)
	{
		float *addr = profdb2_data + (i * n_floats_per_profile);
		// no data to copy: profdb2 will be uninitialized on the device
		prof_db[i].data = addr;
	}
	cudaMemcpy(profdb2_d, &(prof_db[0]), sizeof(DProfile)*n_profiles, cudaMemcpyHostToDevice);

	// the data for the background DProfile
	cudaMemcpy(backh_data,   backh.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
	backh.data  = backh_data;

	// the data for the background DBackground
	cudaMemcpy(backv_data,   backv.data,  sizeof(float)*backv.size, cudaMemcpyHostToDevice);
	backv.data  = backv_data;

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: P: cudaThreadSynchronize failed 8");

    t.stop();
    info << startl << "cuda_n2relmatch(): Copying data to device took " << t << " seconds" << std::endl;

	// call kernel to calculate relative profile
    t.start();

    info << startl << "launching kernel cuda_sib" << std::endl;
    int nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;

    Assert2(match_type.m_rel_type != ident_t, "cuda_n2relmatch: P2: called with m_rel_type == ident_t");

    if (match_type.m_rel_type == sibling_t)
    {
		cuda_sib <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb1_d, profdb2_d, n_profiles, backv);
    }
    else if (match_type.m_rel_type == gen_t)
    {
		cuda_genrc <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb1_d, profdb2_d, n_profiles, match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, backv);
    }
    else if (match_type.m_rel_type == inv_t)
    {
    	bool inverse = true;
		cuda_genrc <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb1_d, profdb2_d, n_profiles, match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2, backv, inverse);

    }
    else
    {
		cuda_rel <<< dim3(nBlocks), dim3(blockSize) >>>
				(profdb1_d, profdb2_d, n_profiles, match_type.m_path1steps, match_type.m_path2steps, backv);
    }

    cudaError_t cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    cudaError_t err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_n2relmatch: Q: cudaThreadSynchronize failed 9");

    t.stop();
    info << startl << "Kernel " << ((match_type.m_rel_type == sibling_t) ? "cuda_sib" : "cuda_rel") << " took " << t << " seconds" << std::endl;

	// call match kernel. upper mode, point to original and sibling arrays
    t.start();
    info << startl << "launching kernel cuda_match6" << std::endl;
    nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;
    cuda_match6 <<< dim3(nBlocks, nBlocks), dim3(blockSize, blockSize) >>>
    		(profdb1_d, n_profiles, profdb2_d, n_profiles, backh, upper, results_d, lr_threshold);

    cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_n2relmatch: R: cudaThreadSynchronize failed 10");

    t.stop();
    info << startl << "cuda_n2relmatch(): Kernel cuda_match6 took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: S: cudaThreadSynchronize failed 11");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb1_data);
	cudaFree(profdb2_data);
	cudaFree(profdb1_d);
	cudaFree(profdb2_d);
	cudaFree(backh_data);
	cudaFree(backv_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_n2relmatch: T: cudaThreadSynchronize failed 12");
}

#endif

// this version calls a kernel to do an n_profile/n_profile match
void cuda_lr(std::vector<DProfile> &prof_db, // Profile database
	         DProfile back,                  // Background
		     CudaLocusInfo   const &locus_info,
	         N2Results &results,
	         float lr_threshold)
{
	int n_floats_per_profile = back.size;
	int n_profiles = prof_db.size();
	int n_floats = n_floats_per_profile * n_profiles;
	Assert(n_floats_per_profile == locus_info.profile_size);

	//
    // allocate DProfiles data and result on device
	//

	// the result
	N2Results *results_d;
    Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(results_d);

    // the data in the profiles
	float *profdb_data = 0;
    Assert2(cudaMalloc((void **) &profdb_data, sizeof(float)*n_floats) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(profdb_data);

    // the DProfiles themselves (containing pointers to the data)
    DProfile *profdb_d = 0;
    Assert2(cudaMalloc((void **) &profdb_d, sizeof(DProfile)*n_profiles) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(profdb_d);

    // the data for the background
	float *back_data = 0;
    Assert2(cudaMalloc((void **) &back_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(back_data);

    //
    // copy to device
    //

    // copy offsets
    COPY_OFFSETS();

    // zero the result
    results.count = 0;
    cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// Data for each DProfile. After copying each data array, copy the device address into the DProfile.
    for (int i=0; i<n_profiles; ++i)
    {
    	float *addr = profdb_data + (i * n_floats_per_profile);
    	cudaMemcpy(addr, prof_db[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
    	prof_db[i].data = addr;
    }

    // the DProfiles themselves (containing pointers to the data)
    cudaMemcpy(profdb_d, &(prof_db[0]), sizeof(DProfile)*n_profiles, cudaMemcpyHostToDevice);

    // the data for the background
    cudaMemcpy(back_data,   back.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
    back.data  = back_data;

    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: A: cudaThreadSynchronize failed");

	// compute LR
    Timer t;
#if 0
    // each thread does N matches
    int nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);
    info << startl << "launching kernel cuda_match4" << std::endl;
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;
    cuda_match4 <<< dim3(nBlocks), dim3(blockSize) >>> (profdb_d, n_profiles, back, results_d, lr_threshold);
#else
    // each thread does 1 match
#ifdef CUDA_MATCH5_HALF
    info << startl << "launching kernel cuda_match5 (Half addressing)" << std::endl;
    int nxBlocks = n_profiles/(2*blockSize) + ((n_profiles/2)%blockSize == 0?0:1);
    int nyBlocks = (n_profiles-1)/blockSize + ((n_profiles-1)%blockSize == 0?0:1);
    info << alignl << "nxBlocks = " << nxBlocks << std::endl;
    info << alignl << "nyBlocks = " << nyBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;
    cuda_match5 <<< dim3(nxBlocks, nyBlocks), dim3(blockSize, blockSize) >>> (profdb_d, n_profiles, back, results_d, lr_threshold);

#else
    info << startl << "launching kernel cuda_match5 (Full addressing)" << std::endl;
    int nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;
    cuda_match5 <<< dim3(nBlocks, nBlocks), dim3(blockSize, blockSize) >>> (profdb_d, n_profiles, back, results_d, lr_threshold);

#endif
#endif
    cudaError_t cts = cudaThreadSynchronize();
    info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
    cudaError_t err = cudaGetLastError();
    info << startl << "cudaGetLastError() == " << err << std::endl;

    Assert2(cts == cudaSuccess, "cuda_lr: B: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "Kernel took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: C: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb_data);
	cudaFree(profdb_d);
	cudaFree(back_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: D: cudaThreadSynchronize failed");
}

// calls a kernel to do an n_profile/profile match
void cuda_n_match(
	std::vector<DProfile> &prof_db, // Profile database
	DProfile prof,                  // Single profile to compare it with
	DProfile back,                  // Background
    CudaLocusInfo const &locus_info,
	NResults &results,
	float lr_threshold)
{
	int n_floats_per_profile = back.size;
	int n_profiles = prof_db.size();
	int n_floats = n_floats_per_profile * n_profiles;
	Assert(n_floats_per_profile == locus_info.profile_size);
	Assert(prof.size == locus_info.profile_size);

	//
    // allocate DProfiles data and result on device
	//
    Timer t;

	// the result
	NResults *results_d;
    Assert2(cudaMalloc((void **) &results_d, sizeof(NResults)*1) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(results_d);

    // the data in the profiles
	float *profdb_data = 0;
    Assert2(cudaMalloc((void **) &profdb_data, sizeof(float)*n_floats) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(profdb_data);

    // the DProfiles themselves (containing pointers to the data)
    DProfile *profdb_d = 0;
    Assert2(cudaMalloc((void **) &profdb_d, sizeof(DProfile)*n_profiles) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(profdb_d);

    // the data for the test profile
    float *prof_data = 0;
    Assert2(cudaMalloc((void **) &prof_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(prof_data);

    // the data for the background
	float *back_data = 0;
    Assert2(cudaMalloc((void **) &back_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(back_data);

    t.stop();
    info << startl << "cuda_n_match(): Allocating memory on device took " << t << " seconds" << std::endl;

    //
    // copy to device
    //
    t.start();

    // copy offsets
    COPY_OFFSETS();

	// zero the result
    results.count = 0;
    cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

	// Data for each DProfile. After copying each data array, copy the device address into the DProfile.
    for (int i=0; i<n_profiles; ++i)
    {
    	float *addr = profdb_data + (i * n_floats_per_profile);
    	cudaMemcpy(addr, prof_db[i].data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
    	prof_db[i].data = addr;
    }

    // the DProfiles themselves (containing pointers to the data)
    cudaMemcpy(profdb_d, &(prof_db[0]), sizeof(DProfile)*n_profiles, cudaMemcpyHostToDevice);

    // the data for the test profile
    cudaMemcpy(prof_data, prof.data, sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
    prof.data = prof_data;

    // the data for the background
    cudaMemcpy(back_data,   back.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
    back.data  = back_data;

    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: A: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_n_match(): Copying data to device took " << t << " seconds" << std::endl;

	// Call kernel to compute LRs
    t.start();

    int nBlocks = n_profiles/blockSize + (n_profiles%blockSize == 0?0:1);

    info << startl << "launching kernel cuda_match3" << std::endl;
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    info << alignl << "n_profiles = " << n_profiles << std::endl;
    cuda_match3 <<< dim3(nBlocks), dim3(blockSize) >>> (profdb_d, n_profiles, prof, back, results_d, lr_threshold);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: B: cudaThreadSynchronize failed");

    t.stop();
    info << startl << "cuda_n_match(): Kernel cuda_match3 took " << t << " seconds" << std::endl;

    // get result from device
    cudaMemcpy(&results, results_d, sizeof(NResults)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: C: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb_data);
	cudaFree(profdb_d);
	cudaFree(prof_data);
	cudaFree(back_data);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: D: cudaThreadSynchronize failed");
}

// this version calls a kernel to do a single profile/profile match
double cuda_lr(
		ConstDProfile prof1,
		ConstDProfile prof2,
		ConstDProfile back,
	    CudaLocusInfo const &locus_info)
{
	float result_h = 0;
	int n = back.size;
	Assert(n == locus_info.profile_size);
	Assert(prof1.size == locus_info.profile_size);
	Assert(prof2.size == locus_info.profile_size);

    // allocate DProfiles data and result on device
    float *result_d = 0;     // pointer to device memory
    Assert2(cudaMalloc((void **) &result_d, sizeof(float)*1) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(result_d);

	float *prof1_d = 0;
    Assert2(cudaMalloc((void **) &prof1_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(prof1_d);

	float *prof2_d = 0;
    Assert2(cudaMalloc((void **) &prof2_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(prof2_d);

	float *back_d = 0;
    Assert2(cudaMalloc((void **) &back_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(back_d);

    // copy offsets
    COPY_OFFSETS();

	// load onto device
    cudaMemcpy(result_d, &result_h, sizeof(float)*1, cudaMemcpyHostToDevice);
    cudaMemcpy(prof1_d, prof1.data, sizeof(float)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(prof2_d, prof2.data, sizeof(float)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(back_d,  back.data,  sizeof(float)*n, cudaMemcpyHostToDevice);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    // store the device addresses
    DProfile p1, p2, b;
    p1.data = prof1_d; p1.size = locus_info.profile_size; // cuda_prof_size;
    p2.data = prof2_d; p2.size = locus_info.profile_size; // cuda_prof_size;
    b.data  = back_d;  b.size  = back.size;

	// compute LR
//    int blockSize = 8;
//    int nBlocks = n/blockSize + (n%blockSize == 0?0:1);
    int blockSize = 1;
    int nBlocks = 1;
    info << startl << "launching kernel cuda_match2" << std::endl;
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    cuda_match2 <<< dim3(nBlocks), dim3(blockSize) >>> (p1, p2, b, result_d);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    // get result from device
    cudaMemcpy(&result_h, result_d, sizeof(float)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    info << startl << 1 << " results copied from device" << std::endl;

	// clean up
	cudaFree(result_d);
	cudaFree(prof1_d);
	cudaFree(prof2_d);
	cudaFree(back_d);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    return (double)result_h;
}

// Performs a single locus/locus match (calculates LR) on the CUDA hardware.
// NB it is not efficient to do just one at a time - this is a test
//
// loc1, loc2 and background must have corresponding entries
double cuda_lr(const float *loc1, const float *loc2, const float *background, int n)
{
	float sum = 0;

    // allocate arrays and result on device
    float *loc1_d=0, *loc2_d=0, *back_d=0, *result_d=0;     // pointers to device memory
    Assert2(cudaMalloc((void **) &loc1_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaMalloc((void **) &loc2_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaMalloc((void **) &back_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaMalloc((void **) &result_d, sizeof(float)*n) != cudaErrorMemoryAllocation,
    		"cuda_lr: cudaMalloc failed");
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");
    Assert(loc1_d && loc2_d && back_d && result_d);

    // allocate result on host
    float *result_h = (float*)malloc(sizeof(float)*n);

	// load onto device
    cudaMemcpy(loc1_d, loc1,      sizeof(float)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(loc2_d, loc2,      sizeof(float)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(back_d,  background, sizeof(float)*n, cudaMemcpyHostToDevice);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

	// compute LR
    int nBlocks = n/blockSize + (n%blockSize == 0?0:1);
    info << startl << "launching kernel cuda_match1" << std::endl;
    info << alignl << "nBlocks = " << nBlocks << std::endl;
    info << alignl << "blockSize = " << blockSize << std::endl;
    cuda_match1 <<< dim3(nBlocks), dim3(blockSize) >>> (loc1_d, loc2_d, back_d, n, result_d);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    // get result from device
    cudaMemcpy(result_h, result_d, sizeof(float)*n, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    info << startl << n << " results copied from device" << std::endl;

	// clean up
	cudaFree(loc1_d);
	cudaFree(loc2_d);
	cudaFree(back_d);
	cudaFree(result_d);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_lr: cudaThreadSynchronize failed");

    // sum on the CPU
    // to sum this on the device we need to do a reduction (see SDK reduction example)
	sum = 0;
	for (int i=0; i<n; ++i)
	{
		sum += result_h[i];
	}

	return sum;
}

__host__
void
runCudaSPMC(
	int 		nBlocks,
	int 		blockSize,
	DProfile 	*prof_db,    // profile dataset
	DProfile 	*spmc_db,    // correction factor matrices (to construct)
	int 		n,           // size of dataset
	DProfile    back,        // background matrix (4.4)
    DBackground backv,       // background vector (HW)
    float       theta)       // Fst
{
	cuda_spmc <<< dim3(nBlocks), dim3(blockSize) >>>
			(prof_db, spmc_db, n, back, backv, theta);
}

__host__
void
runCudaSib(
	int 		nBlocks,
	int 		blockSize,
	DProfile 	*prof_db,    // profile dataset
	DProfile 	*sib_db,     // sibling dataset (to construct)
	int 		n,           // size of dataset
	DBackground back,        // background
	const CudaSubPopModel &spm)
{
	cuda_sib <<< dim3(nBlocks), dim3(blockSize) >>>
			(prof_db, sib_db, n, back, spm);
}

__host__
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
	bool inverse)
{
	cuda_genrc <<< dim3(nBlocks), dim3(blockSize) >>>
			(prof_db, rel_db, n, a1, b1, a2, b2, back, inverse);
}

__host__
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
	const CudaSubPopModel &spm)
{
	cuda_rel <<< dim3(nBlocks), dim3(blockSize) >>>
			(prof_db, rel_db, n, d1, d2, back, spm);
}


/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * cuda_match.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: gareth
 *
 *  CUDA match functions
 */

#include "cuda_match.h"

#include <UnitTest++/UnitTest++.h>

#include "fand/MessageStream.h"
INIT_MESSAGES("cuda_match");
#include "fand/messages.h"

// A class to manage the memory of CUDA data structures on the host side.
// (DProfile is a C struct so can't have constructors and stuff).
// Use pointer-like syntax to access the underlying DProfile
// (The template parameter should derive from CUDAData)
template <class T>
class MemManager
{
public:

	// allocate memory and take responsibility for freeing it
	MemManager(int size)
	{
		pset(size, (float*)malloc(size*sizeof(float)), true);
		memset(m_cudadata.data, 0, size * sizeof(float));
	}

	// point to existing memory (don't take responsibility for freeing it)
	MemManager(int size, float *data)
	{
		pset(size, data, false);
	}

	// copy constructor (increments reference count)
	MemManager(const MemManager& other)
	{
		m_cudadata.size = other.m_cudadata.size;
		m_cudadata.data = other.m_cudadata.data;
		m_ref_count_ptr = other.m_ref_count_ptr;

		if (m_ref_count_ptr)
		{
			(*m_ref_count_ptr)++;
		}
	}

	//destructor
	~MemManager()
	{
		pcleanup();
	}

	// point to existing memory (optionally take responsibility for freeing it)
	void set(int size, float *data, bool freeme = false)
	{
	    pcleanup();
	    pset(size, data, freeme);
	}

	// assignment (increments reference count)
	MemManager& operator=(const MemManager& other)
	{
		// copy construct and swap
		MemManager tmp(other);
		pswap(tmp);
		return *this;
	}

	// access using pointer notation
	T operator*() { return m_cudadata; }
	T *operator->() { return &m_cudadata; }

private:
	T m_cudadata;
	int *m_ref_count_ptr;   // points to reference count, or 0 if not controlled

	void pcleanup()
	{
		if (m_ref_count_ptr)
		{
			(*m_ref_count_ptr)--;

			if (*m_ref_count_ptr == 0)
			{
				free(m_cudadata.data);
				m_cudadata.data = 0;
				m_cudadata.size = 0;

				delete m_ref_count_ptr;
				m_ref_count_ptr = 0;
			}
		}
	}

	void pset(int size, float *data, bool freeme)
	{
		m_cudadata.size = size;
		m_cudadata.data = data;

		if (freeme)
		{
			m_ref_count_ptr = new int;
			*m_ref_count_ptr = 1;
		}
		else
		{
			m_ref_count_ptr = 0;
		}
	}

	void pswap(MemManager& other)
	{
		std::swap(m_cudadata,    other.m_cudadata);
		std::swap(m_ref_count_ptr, other.m_ref_count_ptr);
	}
};

TEST(inv_p)
{
	int i, j;

	inv_p(4, 5, i, j);
	CHECK_EQUAL(1, i);
	CHECK_EQUAL(2, j);

	inv_p(4, 9, i, j);
	CHECK_EQUAL(3, i);
	CHECK_EQUAL(3, j);
}

void prepareDb(
		ProfileRange &db,
		PopulationData &mpopdata,
		std::vector<DProfile> &dp2vec,
		float *buf,
		double delta)
{
	Timer t1;
	MemMeter mem;

//	t1.stop();
//    info << startl << "*** merging took " << t1 << " seconds" << std::endl;
    info << alignl << "*** memory in use " << mem << std::endl;

    t1.start();

    float *p = buf;
	for(size_t i=0; i<db.size(); ++i)
	{
		db[i].setErrorRate(delta);

		DProfile dp2;

//		CachedArray<float>::ConstTempRef ref = db[i].binRef(mpopdata);
//		dp2.size = ref.numValues();
//		memcpy(p, ref.values(), dp2.size * sizeof(float));

		std::vector<float> const &ref = db[i].bin(mpopdata);
		dp2.size = ref.size();
		memcpy(p, &ref[0], dp2.size * sizeof(float));

        dp2.data = p;
        p += dp2.size;
		dp2vec.push_back(dp2);
	}

    t1.stop();
    info << startl << "*** constructing dp2vec took " << t1 << " seconds" << std::endl;
    info << alignl << "*** memory in use " << mem << std::endl;
}

Profile pBackground(PopulationData &mpopdata)
{
	Profile back = Profile(ProfileData(Identifiler, "Background"));

	// Here we need an entry for every locus that has population data.
	// miss out any that don't. Profile::makeBin() will sort it out.

	for (size_t j=0; j<cuda_num_loci; ++j)
	{
		if (mpopdata.hasLocus(j))
		{
			std::vector< PMF<Allele> > background;
			background.push_back(mpopdata(j));
			background.push_back(mpopdata(j));
			back[j] = AlleleSet(background);
		}
	}

	return back;
}

// construct background as vector (1D)
void
dBackground(PopulationData const &mpopdata, MemManager<DBackground> &back)
{
	// Here we need exactly cuda_num_loci entries, each of the right size
	// Pad with zero (which will be ignored).

	LocusInfo const &locus_info = mpopdata.getLocusInfo();
    Assert(back->size == locus_info.back_size);

	for (int j=0; j<cuda_num_loci; ++j)
	{
        int offset = locus_info.back_offset[j];
        int num_alleles = locus_info.num_alleles[j];

		if (mpopdata.hasLocus(j))
		{
			Assert(mpopdata(j).size() <= (int)CUDA_MAX_ALLELES);
			Assert(mpopdata(j).size() > 0);

			PMF<Allele>::const_iterator it = mpopdata(j).begin();
            for (int i=0; i<num_alleles; ++i)
			{
				if (it != mpopdata(j).end())
				{
					back->data[offset + i] = (it++)->second;
				}
				else
				{
					back->data[offset + i] = 0;
				}
			}
		}
		else
		{
			for (int i=0; i<num_alleles; ++i)
			{
				back->data[offset + i] = 0;
			}
		}
	}
}

// N2 match. Handles relative matches.
int cuda_match_n2(
			ProfileRange &db,
		    struct N2Results &results,
		    double lr_threshold,
		    double delta,
		    MatchType const &match_type,
		    SubPopModel const &spm)
{
	// time the preparation phase
    Timer t, t1, t2;
    MemMeter mem;

	// ensure consistent set of alleles over all the profiles
	PopulationData mpopdata = populationData(); // copy - may be modified
    LocusInfo const &locus_info = mpopdata.getLocusInfo();
    info2 << startl << locus_info << std::endl;

    // construct a profile with the (merged) background PMF (Hardy-Weinberg frequencies)
	Profile pb = pBackground(mpopdata);

	// construct DProfiles for background (HW)
	DProfile dback;
	std::vector<float> const &ref = pb.bin(mpopdata);
	dback.size = ref.size();
    dback.data = new float[dback.size];
    memcpy(dback.data, &ref[0], dback.size * sizeof(float));

    t.stop();
    info << startl << "*** preparation took " << t << " seconds" << std::endl;
    info << alignl << "*** memory in use " << mem << std::endl;

	// time the match phase
    t.start();

	// perform N^2 match on the CUDA device
	cuda_stream_match(db, db, match_type, dback, locus_info, results, lr_threshold, delta, delta, upper, mpopdata, spm);

	int nresults = results.count;
	if (results.count > n2result_max)
	{
		warn << startl << "cuda_lr: results overflow: " << nresults << " matches truncated to " << n2result_max << std::endl;
		std::cout << "WARNING: RESULTS OVERFLOW. " << nresults << " matches truncated to " << n2result_max << std::endl;
		std::cout << "(Increase LR threshold until this warning does not appear)" << std::endl;
		nresults = n2result_max;
	}

	t.stop();
	info << startl << "*** match on GPU took " << t << " seconds" << std::endl;
	info << alignl << "*** memory in use " << mem << std::endl;

	delete dback.data;

	return nresults;
}

// NM match.
int cuda_match_nm(
		ProfileRange &db1,
		ProfileRange &db2,
		struct N2Results &results,
		double lr_threshold,
		double delta1,
		double delta2,
		MatchType const &match_type,
	    SubPopModel const &spm,
		ArrayPart part)
{
	// time the preparation phase
    Timer t, t1, t2;
    MemMeter mem;

	// ensure consistent set of alleles over all the profiles
	PopulationData mpopdata = populationData(); // copy - may be modified

    LocusInfo const &locus_info = mpopdata.getLocusInfo();
    info2 << startl << locus_info << std::endl;

//	t1.stop();
//    info << startl << "*** merging took " << t1 << " seconds" << std::endl;
    info << alignl << "*** memory in use " << mem << std::endl;

    t1.start();
	std::vector<DProfile> dpvec1, dpvec2;

    // construct a profile with the (merged) background PMF (Hardy-Weinberg frequencies)
	Profile pb = pBackground(mpopdata);

	// construct DProfiles for background
	DProfile dback;
	std::vector<float> const &ref = pb.bin(mpopdata);
	dback.size = ref.size();
    dback.data = new float[dback.size];
    memcpy(dback.data, &ref[0], dback.size * sizeof(float));

    t.stop();
    info << startl << "*** preparation took " << t << " seconds" << std::endl;
    info << alignl << "*** memory in use " << mem << std::endl;

	// time the match phase
    t.start();

	// perform NM match on the CUDA device
	cuda_stream_match(db1, db2, match_type, dback, locus_info, results, lr_threshold, delta1, delta2, part, mpopdata, spm);
	
	int nresults = results.count;
	if (results.count > n2result_max)
	{
		warn << startl << "cuda_lr: results overflow: " << nresults << " matches truncated to " << n2result_max << std::endl;
		std::cout << "WARNING: RESULTS OVERFLOW. " << nresults << " matches truncated to " << n2result_max << std::endl;
        std::cout << "(Increase LR threshold until this warning does not appear)" << std::endl;
		nresults = n2result_max;
	}

	t.stop();
	info << startl << "*** match on GPU took " << t << " seconds" << std::endl;
	info << alignl << "*** memory in use " << mem << std::endl;

	delete dback.data;

	return nresults;
}

int
cuda_match_n(
    ProfileRange &db,
    const Profile& p,
    struct NResults &results,
    double lr_threshold,
    double delta,
    const MatchType &match_type,
    SubPopModel const &spm)
{
	// just call cuda_match_nm

    struct N2Results results2;
    results.count = 0;

    std::vector<Profile> pvec;
    pvec.push_back(p);
    ProfileRangeVec prvec(pvec);
    ProfileRange db1(prvec);

    int nresults = cuda_match_nm(db, db1, results2, lr_threshold, delta, 0, match_type, spm, full);

    for (int i=0; i<nresults; ++i)
    {
		results.index[i] = results2.index1[i];
		results.lr[i]    = results2.lr[i];
    }

	return nresults;
}

// match one db against another (not necessarily the same size)
//
// Streaming version (overlap writes and kernel executions)
//
// There are two reasons to use streaming:
// 1) we can deal with data sets too big to go on the CUDA devices all at once
//    (provided we have enough CPU RAM)
// 2) we can overlap writes and kernels, which should be faster
//
// The basic idea is to put a large CHUNK of DB1 on the device, and then 'stream'
// through smaller chunks of DB2. This is done using two CUDA streams (S0, S1), which will
// be scheduled to overlap copies and kernels.
//
// Pseudocode:
//
//	Allocate PINNED memory areas for each card (in main thread): (CHUNK of DB1) and (BIGCHUNK of DB2)
//
//	FOR EACH (CHUNK of DB1)
//
//		Copy CHUNK of DB1 to PINNED
//		ASYNC Copy CHUNK of DB1 to DEVICE S1
//
//		FOR EACH (BIGCHUNK of DB2)
//
//			FOR EACH (SMALLCHUNK of DB2 in BIGCHUNK)
//
//				Copy SMALLCHUNK of DB2 to PINNED
//				ASYNC Copy CHUNK of DB2 to DEVICE S0
//
//				Copy SMALLCHUNK of DB2 to PINNED
//				ASYNC Copy CHUNK of DB2 to DEVICE S1
//
//				ASYNC Kernel S0	CUDA C
//
//				ASYNC Kernel S1	CUDA C
//
//			END FOR
//
//			SYNCHRONIZE
//
//		END FOR
//
//	END FOR
//
//	Free (PINNED)
//
// It is important to note:
// * The Profile data obtained with Profile::getRef() is cached and must be
//   copied elsewhere to be used. We copy it to pinned memory.
// * The copies to the GPU and the kernels are launched asynchronously, i.e.
//   they are queued on the GPU and control returns to the CPU.
//   Only when we call cudaThreadSynchronize at the end of each BIGCHUNK loop does the CPU wait for all
//   copies and kernels to finish. We make use of overlap only within the SMALLCHUNK loops.
// * Asynchronous copy must be from pinned memory on the CPU.
// * On the GPU, (a) operations within a stream are executed in the order queued.
//               (b) copies are executed in the order queued
//               (c) kernels are executed in the order queued, but
//               (d) otherwise, copies and kernels may be done in either order **or at the same time**
// * This means on the CPU:
//     o We must have the whole of a BIGCHUNK of DB2 in pinned memory at once
//     o We interleave CPU copies of a SMALLCHUNK to pinned with launches of Asynchronous copies
//       to the GPU. That makes sure the data is in pinned memory before it is copied, and overlaps
//       allows CPU copies to overlap with GPU operations.
// * On the GPU:
//     o Since a cpu2gpu copy (in each stream) happens before the corresponding kernel, we need only
//       hold two SMALLCHUNKs of DB2 data on the GPU (one for each stream). However we need the whole of
//       the DB1 CHUNK
//

// Copy chunk of Profile data to device
// The DProfile array must be pinned and contiguous on the host, starting at dp_host[0]
// The data array must be pinned and contiguous on the host, starting at dp_host[0].data
// NB this copy is asynchronous.
void copyPinnedToDevice(
		ConstDProfile  *dp_host,     // host data (pinned, to copy)
		DProfile       *dp_device,   // device DProfile array
		float          *data_device, // device data
		int 			n,           // number of floats in a profile
		int 			size,        // size of chunk (Profiles)
		cudaStream_t   &stream)
{
	const float *data_host = dp_host[0].data; // the start of the data array

	// put the device addresses in the host DProfiles array
	// (Note that different chunks of DB2 are using different parts of this array, so it is OK to modify it)
	float *p = data_device;
	for (int k=0; k<size; ++k)
	{
		dp_host[k].data = p;
		int profile_size = dp_host[k].size;
		Assert(dp_host[k].size == n); // all profiles should be size n
		p += dp_host[k].size;
	}

	// copy the data
	cudaMemcpyAsync(data_device, data_host, sizeof(float) * size * n,
					cudaMemcpyHostToDevice, stream);

	// copy the DProfiles (ditto)
	cudaMemcpyAsync(dp_device, dp_host, sizeof(DProfile)*size,
					cudaMemcpyHostToDevice, stream);
}

// Copy the Profile range into pinned memory
void
copyProfilesToPinned(
		ProfileRange &prof_db,                // IN  - Profiles
		PopulationData mpopdata,              // IN  - population data (merged)
		float    *data_pinned,                // OUT - pinned copy of data
		ConstDProfile *dp_pinned,             // OUT - pinned copy of DProfiles
		double delta,
		SubPopModel const &spm)
{
	MemMeter m;

	float *data_p = data_pinned;
	int size = prof_db.size();
	for(size_t i=0; i<size; ++i)
	{
		ConstDProfile &dpx = (dp_pinned)[i];

//		info << alignl << "copyProfilesToPinned(): memory in use (C0): " << m << endl;

//		CachedArray<float>::ConstTempRef ref = prof_db[i].binRef(mpopdata);
		prof_db[i].setErrorRate(delta);

		std::vector<float> const &ref = prof_db[i].bin(mpopdata, spm);

//		info << alignl << "copyProfilesToPinned(): memory in use (C1): " << m << endl;

		dpx.size = ref.size();
		dpx.data = &ref[0];                                 // point to the cached/calculated data
		memcpy(data_p, dpx.data, dpx.size * sizeof(float)); // copy to pinned memory
		dpx.data = data_p;                                  // point to the copy
		data_p += dpx.size;
		Assert(dpx.size != 0);
		Assert(dpx.data != 0);
	}
}

// Allocate a contiguous data area and a contiguous array of ConstDProfile in pinned memory
void
allocPinned(
		int             db_size,               // IN  - number of profiles
		int             n_floats_per_profile,  // IN  - number of floats per profile
		float         **data_pinned,           // OUT - pinned copy of data
		ConstDProfile **dp_pinned)             // OUT - pinned copy of DProfiles
{
	unsigned long db_floats = n_floats_per_profile * db_size; // number of floats for db data
	cudaHostAlloc((void**)data_pinned, (unsigned long)db_floats * sizeof(float), cudaHostAllocDefault);
	cudaHostAlloc((void**)dp_pinned,   (unsigned long)db_size * sizeof(ConstDProfile), cudaHostAllocDefault);

	Assert(*data_pinned);
	Assert(*dp_pinned);
}

void cuda_stream_match(
		ProfileRange &prof_db1,
		ProfileRange &prof_db2,
		MatchType const &match_type, // Relative type
		DProfile back,               // Background (HW frequencies)
 	 	CudaLocusInfo const &locus_info,
		N2Results &results,
		float lr_threshold,
		double delta1,
		double delta2,
		ArrayPart part,
		PopulationData const &mpopdata,
		SubPopModel spm)
{
	if (spm.type == SubPopModel::B11 || spm.type == SubPopModel::NRC4_4)
	{
		if (match_type.m_rel_type != ident_t)
		{
			spm.type = SubPopModel::B11; // use NRC4_4 for ident only
		}

		if (match_type.m_rel_type == gen_t || match_type.m_rel_type == inv_t
				|| (match_type.m_rel_type == degree_pq_t && match_type.m_path2steps != MatchType::INF) )
		{
			Assert2(false, "Subpopulation model not currently supported for this relationship type");
		}
	}
	else if (spm.type != SubPopModel::HW)
	{
		Assert2(spm.type == SubPopModel::B11, "cuda_stream_match(): only HW, NRC4_4 and B11 spm supported");
	}

	// background as a vector
    MemManager<DBackground> vback(locus_info.back_size);
	dBackground(mpopdata, vback);
	DBackground backv = *vback; // a copy we can change

	// for SPMs we need the 4.4 background frequencies, not HW
	if (spm.type != SubPopModel::HW)
	{
		for (int loc=0; loc< cuda_num_loci; ++loc)
		{
			int loc_offset  = locus_info.locus_offset[loc];
			int locus_size  = locus_info.locus_size[loc];
			int back_offset = locus_info.back_offset[loc];
			int n_alleles   = locus_info.num_alleles[loc];

			if (locus_size == 0) // locus not in population database
			{
				continue;
			}

			int k = 0; // index into back
			for (int i=0; i<n_alleles; ++i)
			{
				for (int j=i; j<n_alleles; ++j)
				{
			        float fi = backv.data[back_offset + i]; // HW frequency of i
			        float fj = backv.data[back_offset + j]; // HW frequency of j

					back.data[loc_offset + k] = spm.prob((i==j), fi, fj);
					++k;
				}
			}
		}
	}

	TaskData const &my_task = gpuDevices().getTaskData();

	size_t n_floats_per_profile = back.size;
	size_t n1 = prof_db1.size();
	size_t n2 = prof_db2.size();
	size_t n_floats1 = n_floats_per_profile * my_task.db1_chunk;
	size_t n_floats2big = n_floats_per_profile * my_task.db2_bigchunk;
	size_t n_floats2small = n_floats_per_profile * my_task.db2_smallchunk;

	Assert(n_floats_per_profile == locus_info.profile_size);

	//
	// allocate DProfiles data and result on device
	//
	Timer t;
	MemMeter m;
	info << startl << "cuda_stream_match(): memory in use at start: " << m << endl;

	// the result
	N2Results *results_d;
	Assert2(cudaMalloc((void **) &results_d, sizeof(N2Results)*1) != cudaErrorMemoryAllocation,
			"cuda_stream_match: AA: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: A: cudaThreadSynchronize failed");
	Assert(results_d);

	// the data in the profiles
	float *profdb1_data = 0;
	Assert2(cudaMalloc((void **) &profdb1_data, sizeof(float)*n_floats1) != cudaErrorMemoryAllocation,
			"cuda_stream_match: B: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: C: cudaThreadSynchronize failed");
	Assert(profdb1_data);

	float *profdb2_data_s0 = 0;
	Assert2(cudaMalloc((void **) &profdb2_data_s0, sizeof(float)*n_floats2small) != cudaErrorMemoryAllocation,
			"cuda_stream_match: D0: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: E0: cudaThreadSynchronize failed");
	Assert(profdb2_data_s0);

	float *profdb2_data_s1 = 0;
	Assert2(cudaMalloc((void **) &profdb2_data_s1, sizeof(float)*n_floats2small) != cudaErrorMemoryAllocation,
			"cuda_stream_match: D1: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: E1: cudaThreadSynchronize failed");
	Assert(profdb2_data_s1);

	// the DProfiles themselves (containing pointers to the data)
	DProfile *profdb1_dp = 0;
	Assert2(cudaMalloc((void **) &profdb1_dp, sizeof(DProfile)*my_task.db1_chunk) != cudaErrorMemoryAllocation,
			"cuda_stream_match: F: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: G: cudaThreadSynchronize failed");
	Assert(profdb1_dp);

	DProfile *profdb2_dp_s0 = 0;
	Assert2(cudaMalloc((void **) &profdb2_dp_s0, sizeof(DProfile)*my_task.db2_smallchunk) != cudaErrorMemoryAllocation,
			"cuda_stream_match: H0: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: I0: cudaThreadSynchronize failed");
	Assert(profdb2_dp_s0);

	DProfile *profdb2_dp_s1 = 0;
	Assert2(cudaMalloc((void **) &profdb2_dp_s1, sizeof(DProfile)*my_task.db2_smallchunk) != cudaErrorMemoryAllocation,
			"cuda_stream_match: H1: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: I1: cudaThreadSynchronize failed");
	Assert(profdb2_dp_s1);

	// the data for the background
	float *back_data = 0;
	Assert2(cudaMalloc((void **) &back_data, sizeof(float)*n_floats_per_profile) != cudaErrorMemoryAllocation,
			"cuda_stream_match: J: cudaMalloc failed");
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: K: cudaThreadSynchronize failed");
	Assert(back_data);

	// background vector (used only for kinship matches)
	float *backv_data = 0;
	Assert2(cudaMalloc((void **) &backv_data, sizeof(float)*locus_info.back_size) != cudaErrorMemoryAllocation,
			"cuda_nm_relmatch: M: cudaMalloc failed");

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_nm_relmatch: N: cudaThreadSynchronize failed 7");
	Assert(backv_data);

	// Subpopulation correction factor matrices for DB1
	float *profdb1spmc_data = 0;
	DProfile *profdb1spmc_dp = 0;

	if (spm.type == SubPopModel::B11)
	{
		// space for the correction matrices
		Assert2(cudaMalloc((void **) &profdb1spmc_data, sizeof(float)*n_floats1) != cudaErrorMemoryAllocation,
				"cuda_stream_match: SA: cudaMalloc failed");
		Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: SB: cudaThreadSynchronize failed");
		Assert(profdb1spmc_data);

		Assert2(cudaMalloc((void **) &profdb1spmc_dp, sizeof(DProfile)*my_task.db1_chunk) != cudaErrorMemoryAllocation,
				"cuda_stream_match: SC: cudaMalloc failed");
		Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: SD: cudaThreadSynchronize failed");
		Assert(profdb1spmc_dp);

		// Make the DProfile point to the data on the device
		int n1 = my_task.db1_chunk;
		std::vector<DProfile> dp(n1);
		float *addr = profdb1spmc_data;
		for (int i=0; i<n1; ++i)
		{
			dp[i].size = n_floats_per_profile;
			dp[i].data = addr;
			addr += n_floats_per_profile;
		}

		cudaMemcpy(profdb1spmc_dp, &(dp[0]), sizeof(DProfile)*n1, cudaMemcpyHostToDevice);
		Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: SE: cudaThreadSynchronize failed");
	}

    // Allocate pinned memory areas
	float         *data1_pinned, *data2_pinned;
	ConstDProfile *dp1_pinned,   *dp2_pinned;
	allocPinned(my_task.db1_chunk,    n_floats_per_profile, &data1_pinned, &dp1_pinned);
	allocPinned(my_task.db2_bigchunk, n_floats_per_profile, &data2_pinned, &dp2_pinned);

    // copy locus info to device
	copyOffsets(locus_info);

	// copy data for the background DProfile
	cudaMemcpy(back_data,   back.data,  sizeof(float)*n_floats_per_profile, cudaMemcpyHostToDevice);
	back.data  = back_data;
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: L: cudaThreadSynchronize failed");

	// the data for the background DBackground (used for kinship matches and PM correction)
	cudaMemcpy(backv_data,   backv.data,  sizeof(float)*locus_info.back_size, cudaMemcpyHostToDevice);
	backv.data  = backv_data;
	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: L2: cudaThreadSynchronize failed");

	// zero the results
	results.count = 0;
	cudaMemcpy(results_d, &results, sizeof(NResults)*1, cudaMemcpyHostToDevice);

    t.stop();
    info << startl << "cuda_stream_match(): Allocating memory on device took " << t << " seconds (t = " << t.stop() << ")"<< std::endl;

    t.start();

	// Initialize streams
	cudaStream_t stream0, stream1;
	cudaError_t cts, err;
	cts = cudaStreamCreate(&stream0);
    Assert2(cts == cudaSuccess, "cuda_stream_match: M1: cudaStreamCreate failed");

	cts = cudaStreamCreate(&stream1);
    Assert2(cts == cudaSuccess, "cuda_stream_match: M2: cudaStreamCreate failed");

    // loop over chunks of DB1
	for (size_t i = 0; i < n1; i += my_task.db1_chunk)
	{
		info << startl << "i = " << i << ", cached memory = " << Cached::getBytesUsed() << " bytes" << endl;

		info << alignl << "cuda_stream_match(): memory in use (0): " << m << endl;

		size_t db1_chunk = min(my_task.db1_chunk, n1 - i);

		Timer t2;
	    t2.start();

		// Copy chunk of BD1 to pinned memory

	    ProfileRange prof_db1_chunk(prof_db1,
									prof_db1.begin() + i,
									prof_db1.begin() + i + db1_chunk); // chunk to copy

		info << alignl << "cuda_stream_match(): memory in use (1): " << m << endl;

		prof_db1_chunk.readMyData(); // read from database (if data is in the MySQL database)

		info << alignl << "cuda_stream_match(): memory in use (2): " << m << endl;

		copyProfilesToPinned(prof_db1_chunk, mpopdata, data1_pinned, dp1_pinned, delta1, spm);

		info << alignl << "cuda_stream_match(): memory in use (3): " << m << endl;

		t2.stop();
	    info << startl << "cuda_stream_match(): copy chunk of DB1 (" << db1_chunk << " profiles) to pinned took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;
	    t2.start();

		// Copy chunk of DB1 to device in stream1. (NB: so that it can't overlap with the stream 1 kernel
		// that is executed last in the j-loop below. Both streams kernels use DB1 data).
		copyPinnedToDevice(dp1_pinned, profdb1_dp, profdb1_data, n_floats_per_profile, db1_chunk, stream1);

		info << alignl << "cuda_stream_match(): memory in use (4): " << m << endl;

	    t2.stop();
	    info << startl << "cuda_stream_match(): copy chunk of DB1 to device took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;

		// calculate subpopulation correction matrices for chunk of DB1
	    // (NB before we transform the DB1 matrix into that of a relative)
		if (spm.type == SubPopModel::B11)
		{
		    int nBlocks = db1_chunk/blockSize + (db1_chunk%blockSize == 0?0:1);
	    	runCudaSPMC(nBlocks, blockSize, profdb1_dp, profdb1spmc_dp,
	    			    db1_chunk, back, backv, spm.theta_bar);

			cts = cudaThreadSynchronize(); // all streams
//				info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
			err = cudaGetLastError();
//				info << startl << "cudaGetLastError() == " << err << ": " << cudaGetErrorString(err) << std::endl;
			Assert2(cts == cudaSuccess, "cuda_stream_match: M5: cudaThreadSynchronize failed");

		}

		// Kinship match. Transform DB1 matrix into that of the required relative.
		if (match_type.m_rel_type != ident_t)
		{
		    t2.start();

			// Use profdb2_data_s0 array of size SMALLCHUNK on the device as a temporary storage for db1
			// data while we calculate the relatives profiles. We have the whole db1_chunk on the device.
			// Copy small chunks of it into profdb2_data_s0, and call a kernel to calculate the relative
			// profile, which goes back in the db1 array, where it came from.

			// TODO we are doing this asynchronously, but it could be streamed.

			// Make sure pinned copy has been done
			cts = cudaThreadSynchronize(); // all streams
			info << startl << "cudaThreadSynchronize() == " << cts << std::endl;

			// make profdb2_dp_s0 point to data in profdb2_data_s0
			float *addr = profdb2_data_s0;
			for (size_t j=0; j<my_task.db2_smallchunk; ++j)
			{
				dp2_pinned[j].data = addr;
				dp2_pinned[j].size = n_floats_per_profile;
				addr += n_floats_per_profile;
			}
			cudaMemcpy(profdb2_dp_s0, &(dp2_pinned[0]), sizeof(DProfile)*my_task.db2_smallchunk, cudaMemcpyHostToDevice);
			err = cudaGetLastError();
			info << startl << "cudaGetLastError() == " << err << ": " << cudaGetErrorString(err) << std::endl;
			Assert2(err == cudaSuccess, "cuda_stream_match: M3a: cudaMemcpy failed");

			// divide db1_chunk into my_task.db2_smallchunk sized bits
			for (size_t k = 0; k < db1_chunk; k += my_task.db2_smallchunk)
			{
				Timer t3;
			    t3.start();

				int small_chunk = min(my_task.db2_smallchunk, db1_chunk - k);
				float *db1_data_start = profdb1_data + k * n_floats_per_profile;
				DProfile *db1_dp_start = profdb1_dp + k;

				// copy small_chunk of DB1 array into profdb2_data_s0
				cudaMemcpy(profdb2_data_s0, db1_data_start, sizeof(float)*small_chunk*n_floats_per_profile,
						cudaMemcpyDeviceToDevice);
				err = cudaGetLastError();
//				info << startl << "cudaGetLastError() == " << err << ": " << cudaGetErrorString(err) << std::endl;
				Assert2(err == cudaSuccess, "cuda_stream_match: M3b: cudaMemcpy failed");

				// call kernel to calculate relative profiles and put back in DB1 array
			    int nBlocks = small_chunk/blockSize + (small_chunk%blockSize == 0?0:1);

			    if (match_type.m_rel_type == sibling_t)
			    {
			    	runCudaSib(nBlocks, blockSize, profdb2_dp_s0, db1_dp_start, small_chunk, backv, spm);
			    }
			    else if (match_type.m_rel_type == gen_t)
			    {
			    	bool inverse = false;
			    	runCudaGenrc(nBlocks, blockSize, profdb2_dp_s0, db1_dp_start, small_chunk,
			    			     match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2,
			    			     backv, inverse);
			    }
			    else if (match_type.m_rel_type == inv_t)
			    {
			    	bool inverse = true;
			    	runCudaGenrc(nBlocks, blockSize, profdb2_dp_s0, db1_dp_start, small_chunk,
			    			     match_type.m_a1, match_type.m_b1, match_type.m_a2, match_type.m_b2,
			    			     backv, inverse);
			    }
			    else
			    {
			    	runCudaRel(nBlocks, blockSize, profdb2_dp_s0, db1_dp_start, small_chunk,
			    			   match_type.m_path1steps, match_type.m_path2steps,
			    			   backv, spm);
			    }

				cts = cudaThreadSynchronize(); // all streams
//				info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
				err = cudaGetLastError();
//				info << startl << "cudaGetLastError() == " << err << ": " << cudaGetErrorString(err) << std::endl;
				Assert2(cts == cudaSuccess, "cuda_stream_match: M4: cudaThreadSynchronize failed");

			    t3.stop();
			    info << startl << "cuda_stream_match(): calculate relatives for small chunk " << k << " of DB1 took " << t3 << " seconds (t = " << t.stop() << ")"<< std::endl;
			}

		    t2.stop();
		    info << startl << "cuda_stream_match(): calculate relatives for chunk of DB1 took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;
		}


		// declared here so the debugger can see them (!)
		size_t j, j1, db2_chunks, db2_chunk_s0, db2_chunk_s1;

		info << alignl << "cuda_stream_match(): memory in use (5): " << m << endl;

	    // loop over big chunks of DB2
		for (size_t bigchunk_begin = 0; bigchunk_begin < n2; bigchunk_begin += my_task.db2_bigchunk)
		{
			size_t db2_bigchunk = min(my_task.db2_bigchunk, n2 - bigchunk_begin);
			size_t bigchunk_end = bigchunk_begin + db2_bigchunk;

		    float       *data_p = data2_pinned;  // ptr to small chunk data in pinned memory
		    ConstDProfile *dp_p = dp2_pinned;    // ptr to small chunk DProfiles in pinned memory

		    // loop over small chunks of DB2
		    for (/*int*/ j = bigchunk_begin; j < bigchunk_end; j += 2 * my_task.db2_smallchunk)
			{
				/*size_t*/ db2_chunks = min(2*my_task.db2_smallchunk, bigchunk_end - j); // both chunks
				/*size_t*/ db2_chunk_s0 = db2_chunks/2;              // small chunk for stream 0
				/*size_t*/ db2_chunk_s1 = db2_chunks - db2_chunk_s0; // small chunk for stream 1 (may be different for last chunk)

				// NB it does not matter in this case whether we do [copy0, copy1, kernel0, kernel1] or
				// [copy0, kernel0, copy1, kernel1]. But the former is generally safer.

				t2.start();

				// Copy chunk of BD2 to pinned memory for stream 0
				ProfileRange prof_db2_chunk_s0(prof_db2,
											   prof_db2.begin() + j,
											   prof_db2.begin() + j + db2_chunk_s0);

				prof_db2_chunk_s0.readMyData(); // read from database (if data is in the MySQL database)

				copyProfilesToPinned(prof_db2_chunk_s0, mpopdata, data_p, dp_p, delta2, spm);

			    t2.stop();
			    info << startl << "cuda_stream_match(): copy small chunk of DB2 to pinned (stream 0) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;
			    t2.start();

				// Copy chunk of DB2 to device in stream 0
				copyPinnedToDevice(dp_p, profdb2_dp_s0, profdb2_data_s0, n_floats_per_profile, db2_chunk_s0, stream0);

			    t2.stop();
			    info << startl << "cuda_stream_match(): copy small chunk of DB2 to device (stream 0) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;

				// next small chunk
				data_p += db2_chunk_s0 * n_floats_per_profile;
				dp_p   += db2_chunk_s0;

				t2.start();

				// Copy chunk of BD2 to pinned memory for stream 1
				/*size_t*/ j1 = j + db2_chunk_s0;

				ProfileRange prof_db2_chunk_s1(prof_db2,
											   prof_db2.begin() + j1,
											   prof_db2.begin() + j1 + db2_chunk_s1);

				prof_db2_chunk_s1.readMyData(); // read from database (if data is in the MySQL database)

				copyProfilesToPinned(prof_db2_chunk_s1, mpopdata, data_p, dp_p, delta2, spm);

			    t2.stop();
			    info << startl << "cuda_stream_match(): copy small chunk of DB2 to pinned (stream 1) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;
			    t2.start();

				// Copy chunk of DB2 to device stream 1
				copyPinnedToDevice(dp_p, profdb2_dp_s1, profdb2_data_s1, n_floats_per_profile, db2_chunk_s1, stream1);

			    t2.stop();
			    info << startl << "cuda_stream_match(): copy small chunk of DB2 to device (stream 1) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;

				// next small chunk
				data_p += db2_chunk_s1 * n_floats_per_profile;
				dp_p   += db2_chunk_s1;

				int nBlocks1 = db1_chunk/blockSize + (db1_chunk%blockSize == 0?0:1);
				info << alignl << "nBlocks1 = " << nBlocks1 << std::endl;
				info << alignl << "blockSize = " << blockSize << std::endl;
				info << alignl << "db1_chunk = " << db1_chunk << std::endl;

				// Run kernel for stream 0
				int nBlocks2 = db2_chunk_s0/blockSize + (db2_chunk_s0%blockSize == 0?0:1);
				info << alignl << "nBlocks2 = " << nBlocks2 << std::endl;
				info << alignl << "db2_chunk_s0 = " << db2_chunk_s0 << std::endl;

				info << startl << "launching kernel cuda_match6" << std::endl;
				t2.start();
				runCudaMatch6(nBlocks1, nBlocks2, blockSize, stream0, profdb1_dp, db1_chunk,
						      profdb2_dp_s0, db2_chunk_s0, back, part, results_d, lr_threshold,
						      i, j, profdb1spmc_dp);
			    t2.stop();
			    info << startl << "cuda_stream_match(): kernel cuda_match6 (stream 0) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;

				// Run kernel for stream 1
				nBlocks2 = db2_chunk_s1/blockSize + (db2_chunk_s1%blockSize == 0?0:1);
				info << alignl << "nBlocks2 = " << nBlocks2 << std::endl;
				info << alignl << "db2_chunk_s1 = " << db2_chunk_s1 << std::endl;

				info << startl << "launching kernel cuda_match6" << std::endl;

			    t2.start();
				runCudaMatch6(nBlocks1, nBlocks2, blockSize, stream1, profdb1_dp, db1_chunk,
						      profdb2_dp_s1, db2_chunk_s1, back, part, results_d, lr_threshold,
						      i, j1, profdb1spmc_dp);
			    t2.stop();
			    info << startl << "cuda_stream_match(): kernel cuda_match6 (stream 1) took " << t2 << " seconds (t = " << t.stop() << ")"<< std::endl;
			}

		    // We need all kernels to finish before copying a new big chunk of DB2
			cts = cudaThreadSynchronize(); // all streams
			info << startl << "cudaThreadSynchronize() == " << cts << std::endl;
			err = cudaGetLastError();
			info << startl << "cudaGetLastError() == " << err << std::endl;
			Assert2(cts == cudaSuccess, "cuda_stream_match: M4: cudaThreadSynchronize failed");
		}

	}

	info << alignl << "cuda_stream_match(): memory in use (6): " << m << endl;

    t.stop();
    info << startl << "cuda_stream_match(): took " << t << " seconds (t = " << t.stop() << ")"<< std::endl;

    // get results from device
    cudaMemcpy(&results, results_d, sizeof(N2Results)*1, cudaMemcpyDeviceToHost);
    Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: N: cudaThreadSynchronize failed");

    info << startl << min(results.count, n2result_max) << " results copied from device" << std::endl;

	// clean up
	cudaFree(results_d);
	cudaFree(profdb1_data);
	cudaFree(profdb2_data_s0);
	cudaFree(profdb2_data_s1);
	cudaFree(profdb1_dp);
	cudaFree(profdb2_dp_s0);
	cudaFree(profdb2_dp_s1);
	cudaFree(back_data);
	cudaFree(backv_data);
	cudaFree(profdb1spmc_data);
	cudaFree(profdb1spmc_dp);
	cudaFreeHost(data1_pinned);
	cudaFreeHost(dp1_pinned);
	cudaFreeHost(data2_pinned);
	cudaFreeHost(dp2_pinned);

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cuda_stream_match: P: cudaThreadSynchronize failed");
}

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * GPUDevices.h
 *
 *  Created on: Nov 25, 2010
 *      Author: gareth
 *
 * Singleton class representing the GPU computing resource
 */

#ifndef GPUDEVICES_H_
#define GPUDEVICES_H_

#include "fand/Profile.h"
#include "fand/ProfileRange.h"
#include "fand/populationdata.h"
#include "fand/match.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <vector>

static const int MAX_CPU_THREADS = 12;
static const int MAIN_THREAD_NO = MAX_CPU_THREADS - 1;
static const int DEFAULT_CUDA_DEV = 0;
static const int DB2_CHUNK_RATIO = 10;
static const float DB1_CHUNK_F = 0.8;       // fraction of GPU RAM to use for DB1_CHUNK
static const float DEV_F = 0.9;             // fraction of GPU RAM to use for this app

// Data for a calculation task performed on one GPU device in one CPU thread
struct TaskData
{
	int device_no;
	size_t db1_chunk;
	size_t db2_bigchunk;
	size_t db2_smallchunk;
};

extern TaskData taskData[]; // TODO: might be better as a static map in GPUDevices ???

struct Result
{
    Result(const MatchType &type, double lr) : match_type(type), likelihood_ratio(lr) {}

    MatchType match_type;
    double    likelihood_ratio;
};

struct NMmatchResults : public std::map< std::pair<int, int>, std::vector<Result> >
{
    typedef std::map< std::pair<int, int>, std::vector<Result> > base;

    NMmatchResults(int ncomp = 0) : n_comparisons(ncomp) {}
    void clear() { n_comparisons = 0; base::clear(); }

    int n_comparisons;        // number of comparisons made (used for DB size correction)
};

struct NmatchResults : public std::map< int, std::vector<Result> >
{
    typedef std::map< int, std::vector<Result> > base;

    NmatchResults(int ncomp) : n_comparisons(ncomp)  {}
    void clear() { n_comparisons = 0; base::clear(); }

    int n_comparisons;        // number of comparisons made (used for DB size correction)
};

enum Split
{
	split_none = 0, // perform match on one CUDA device (device 0)
	split_parallel  // split task over multiple cuda devices, in parallel
};

typedef std::vector<Profile>::const_iterator ProfileIterator;

struct ThreadData
{
    ThreadData(
		int idx = MAIN_THREAD_NO,
		int dev = DEFAULT_CUDA_DEV,
        size_t ram = 0)
    : task(idx), device(dev), cpu_ram(ram)
    {}

    int  task;    // task (thread) index
    int  device;  // CUDA device to use
    size_t cpu_ram; // application memory to use (B)
};

struct N2ThreadData : public ThreadData
{
	N2ThreadData() : match_type(unknown_t) {}

    MatchType match_type;

    NMmatchResults *matches;
    double lr_threshold;
    double delta;
    SubPopModel spm;

    ProfileRange dbA;
    ProfileRange dbB;
};

struct NMThreadData : public ThreadData
{
	NMThreadData() : match_type(unknown_t) {}

	MatchType match_type;

    NMmatchResults *matches;
    double lr_threshold;
	double delta1;
	double delta2;
    SubPopModel spm;

    ProfileRange dbA1;
	ProfileRange dbA2;
	ProfileRange dbB1;
	ProfileRange dbB2;
};

typedef void *(*FUNC)(void *);

// return the GPUDevices singleton
class GPUDevices;
GPUDevices &gpuDevices();

class GPUDevices
{
public:

    // return number of CUDA-capable GPUs present
	// NB device IDs run from 0 to GPUsPresent() - 1
	//
    int GPUsPresent() const;

    // return number of CUDA-capable GPUs enabled for processing
    // (optionally, return set of device Ids)
    int GPUsEnabled(std::set<int> *device_ids = 0) const;

    // enable n GPUs. Return the number enabled.
    int enableGPUs(int n);

    // enable (only) the GPUs in device_ids. Return the number enabled.
    int enableGPUs(std::set<int> const &device_ids);

    // enable All GPUs. Return the number enabled.
    int enableAllGPUs();

    // enable a single GPU. Return the number enabled.
    int enableDevice(int dev);

    // disable a single GPU. Return the number enabled.
    int disableDevice(int dev);

    // return if device is enabled
    bool isEnabled(int dev) const;

    // return if device has a display attached
    bool hasDisplay(int dev) const;

    // return device global memory (bytes)
    size_t totalGlobalMem(int dev) const;

    // return cudaDeviceProp structure for the device
    cudaDeviceProp deviceProp(int dev) const;

    // match Profile against Profile vector
    int n_match(
            NmatchResults& matches,
            ProfileRange & db,
            Profile& p,
            double lr_threshold,
            double delta,
            MatchType const &match_type,
            SubPopModel const & spm);

    // match Profile vector against itself
    int n2_match(
            NMmatchResults &matches,
            ProfileRange &db,
            double lr_threshold,
            double delta,
            MatchType const &match_type,
            Split split,
            SubPopModel const & spm);

    // match Profile vector against Profile vector
    int nm_match(
            NMmatchResults &matches,
            ProfileRange & dbA,
            ProfileRange & dbB,
            double lr_threshold,
            double delta1,
            double delta2,
            MatchType const &match_type,
            Split split,
            SubPopModel const & spm);

	// Return task data for this thread.
	// NB must be called *within* the thread *after* initThread()
	static TaskData & getTaskData();

	// return task number
	static int getTask();

	// return device used by this task/thread
	static int getDev();

	// Calculate the parameters for running a task on a GPU
	// (Defaults to running on the default CUDA device in the main thread)
    // NB This must be called *within* the thread to ensure thread, task and device are linked
	static bool setupTask(
			int db1_size,
			int db2_size,
			ThreadData const &thread_data =	ThreadData(MAIN_THREAD_NO, DEFAULT_CUDA_DEV),
			GPUDevices const &gpu_devs = gpuDevices(),
			PopulationData const &popdata = populationData());

	static void cleanupTask();

private:
    GPUDevices(); // singleton
    GPUDevices(size_t tot_glob_mem); // test constructor

    // Link this thread to a task.
	static void initThread(ThreadData const &data);

    int n2_match_1gpu(
            NMmatchResults &matches,
            ProfileRange & db,
            double lr_threshold,
            double delta,
            MatchType const &match_type,
            SubPopModel const &spm,
    		ThreadData const &data = ThreadData(MAIN_THREAD_NO, DEFAULT_CUDA_DEV) );

    int n2_match_2gpus(
            NMmatchResults &matches,
            ProfileRange &db,
            double lr_threshold,
            double delta,
            MatchType const &match_type,
            SubPopModel const &spm);

    int n2_match_4gpus(
            NMmatchResults &matches,
            ProfileRange &db,
            double lr_threshold,
            double delta,
            MatchType const &match_type,
            SubPopModel const &spm);

    int nm_match_1gpu(
            NMmatchResults &matches,
            ProfileRange & db1,
            ProfileRange & db2,
            double lr_threshold,
            double delta1,
            double delta2,
            ArrayPart part,
            MatchType const &match_type,
            SubPopModel const &spm,
			ThreadData const &data = ThreadData(MAIN_THREAD_NO, DEFAULT_CUDA_DEV) );

    int nm_match_2gpus(
            NMmatchResults &matches,
            ProfileRange & dbA,
            ProfileRange & dbB,
            double lr_threshold,
            double delta1,
            double delta2,
            ArrayPart part,
            MatchType const &match_type,
            SubPopModel const &spm);

    int nm_match_4gpus(
            NMmatchResults &matches,
            ProfileRange & dbA,
            ProfileRange & dbB,
            double lr_threshold,
            double delta1,
            double delta2,
            ArrayPart part,
            MatchType const &match_type,
            SubPopModel const &spm);

    void n2matchcpy(NMmatchResults &matches, NMmatchResults &tile, int istart, int jstart);

    std::vector<cudaDeviceProp> m_devices;
    std::set<int> m_devices_enabled;

    friend GPUDevices & gpuDevices();
    friend void threadN2Split_4tasks(N2ThreadData *data);
    friend void threadNMSplit_4tasks(NMThreadData *data);
    friend class TestsetupTask;                   // unit tests
    friend class Pop_fixtureProfileRangeDBHelper; //

    static __thread int m_task; // Per thread variable. This is the index used to keep track of
								// the data held in TaskData.

    // An object that will be constructed at static initialization time
    // to ensure the main thread is initialized
    static struct Initializer
    {
    	Initializer()
    	{
    		// ensure main thread is initialized
    		gpuDevices().initThread(ThreadData(MAIN_THREAD_NO, DEFAULT_CUDA_DEV));
    	}
    } m_initializer;
};

#endif /* GPUDEVICES_H_ */

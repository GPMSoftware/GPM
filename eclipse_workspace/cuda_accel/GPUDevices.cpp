/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * GPUDevices.cpp
 *
 *  Created on: Nov 25, 2010
 *      Author: gareth
 */

#include "GPUDevices.h"
#include "cuda_match.h"
#include "fand/HostManager.h"
#include <UnitTest++/UnitTest++.h>

#include <algorithm>

#include "fand/MessageStream.h"
INIT_MESSAGES("GPUDevices");
#include "fand/messages.h"
using namespace std;

void threadN2Split_4tasks(N2ThreadData *data);
void threadNMSplit_4tasks(NMThreadData *data);

__thread int GPUDevices::m_task = -1; // per thread variable

GPUDevices::Initializer GPUDevices::m_initializer; // object that does stuff at static initialization time

TaskData taskData[MAX_CPU_THREADS];

GPUDevices &
gpuDevices()
{
    static GPUDevices theGPUDevices;

    return theGPUDevices;
}

GPUDevices::GPUDevices()
{
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);

    int dev;

    if (deviceCount == 0)
    {
        cout << "There is no device supporting CUDA: "; // emulation ?
    }
    else
    {
    	cout << "There are " << deviceCount << " devices supporting CUDA: ";

		for (dev = 0; dev < deviceCount; ++dev)
		{
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);
			m_devices.push_back(deviceProp);

			// by default, enable all devices that are not attached to the display
			if (deviceProp.kernelExecTimeoutEnabled == 0)
			{
				m_devices_enabled.insert(dev);
			}
		}

		// but if all cuda devices are running displays, enable the first device
		if (m_devices_enabled.size() == 0)
		{
	        m_devices_enabled.insert(0);
		}
    }
    cout << m_devices_enabled.size() << " devices enabled" << endl;
}

// test constructor
GPUDevices::GPUDevices(size_t tot_glob_mem)
{
	// create a fiction GPU for unit tests
	cudaDeviceProp deviceProp;
	deviceProp.totalGlobalMem = tot_glob_mem;
	m_devices.push_back(deviceProp);
}

// return number of CUDA-capable GPUs present
// (optionally, return vector of device Ids)
int
GPUDevices::GPUsPresent() const
{
    return m_devices.size();
}

// return number of CUDA-capable GPUs enabled for processing
// (optionally, return vector of device Ids)
int
GPUDevices::GPUsEnabled(set<int> *device_ids) const
{
    if (device_ids)
    {
        *device_ids = m_devices_enabled;
    }
	return m_devices_enabled.size();
}

// enable n GPUs. Return the number enabled.
int
GPUDevices::enableGPUs(int n)
{
	m_devices_enabled.clear();

	// first enable non-display devices
	for (int dev = 0; dev < m_devices.size(); ++dev)
	{
		if (m_devices[dev].kernelExecTimeoutEnabled == 0)
		{
			m_devices_enabled.insert(dev);
		}
		if (n == m_devices_enabled.size()) return n;
	}

	// then enable display devices
	for (int dev = 0; dev < m_devices.size(); ++dev)
	{
		if (m_devices[dev].kernelExecTimeoutEnabled > 0)
		{
			m_devices_enabled.insert(dev);
		}
		if (n == m_devices_enabled.size()) return n ;
	}

	// enabled all devices but it is not as many as requested
	return m_devices_enabled.size();
}

// enable (only) the GPUs in device_ids. Return the number enabled.
int
GPUDevices::enableGPUs(set<int> const &device_ids)
{
	m_devices_enabled.clear();

	set<int>::const_iterator it;
	for (it = device_ids.begin(); it != device_ids.end(); ++it)
	{
		int dev = *it;
		if (dev < m_devices.size())
		{
			m_devices_enabled.insert(dev);
		}
	}
	return m_devices_enabled.size();
}

// enable All GPUs. Return the number enabled.
int
GPUDevices::enableAllGPUs()
{
	for (int dev = 0; dev < m_devices.size(); ++dev)
	{
		m_devices_enabled.insert(dev);
	}
	return m_devices_enabled.size();
}

// enable a single GPU. Return the number enabled.
int
GPUDevices::enableDevice(int dev)
{
	if (dev < m_devices.size())
	{
		m_devices_enabled.insert(dev);
	}
	return m_devices_enabled.size();
}

// disable a single GPU. Return the number enabled.
int
GPUDevices::disableDevice(int dev)
{
	if (dev < m_devices.size())
	{
		m_devices_enabled.erase(dev);
	}
	return m_devices_enabled.size();
}

bool
GPUDevices::isEnabled(int dev) const
{
	set<int>::const_iterator it;
	if (m_devices_enabled.find(dev) == m_devices_enabled.end())
	{
		return false;
	}

	return true;
}

bool
GPUDevices::hasDisplay(int dev) const
{
	if (dev < m_devices.size() &&
	    m_devices[dev].kernelExecTimeoutEnabled > 0)
	{
		return true;
	}

	return false;
}

size_t
GPUDevices::totalGlobalMem(int dev) const
{
	if (dev < m_devices.size())
	{
		return m_devices[dev].totalGlobalMem;
	}

	return 0;
}

cudaDeviceProp
GPUDevices::deviceProp(int dev) const
{
	cudaDeviceProp ret;

	if (dev < m_devices.size())
	{
		ret = m_devices[dev];
	}

	return ret;
}

// match Profile against Profile database (GPU)
int
GPUDevices::n_match(
        NmatchResults& matches,
    	ProfileRange &pr,
        Profile& p,
        double lr_threshold,
        double delta,
        MatchType const &match_type,
        SubPopModel const & spm)
{

    p.setErrorRate(delta);

    // run on default device in main thread
    if ( ! setupTask(pr.size(), 1) )
    {
    	return 0;
    }

    struct NResults results;
    results.count = 0;

    int nresults = cuda_match_n(pr, p, results, lr_threshold, delta, match_type, spm);

	cleanupTask();

    for (int j=0; j<nresults; ++j)
    {
        matches[results.index[j]].push_back(Result(match_type, results.lr[j]));
    }

    return nresults;
}

// match Profile database against itself (GPU)
int
GPUDevices::n2_match(
        NMmatchResults &matches,
        ProfileRange &db,
        double lr_threshold,
        double delta,
        MatchType const &match_type,
        Split split,
        SubPopModel const & spm)
{

    if (split != split_parallel /*|| db.size() < 2*/)
    {
        // run on one card only (whichever is the default)
        (void)n2_match_1gpu(matches, db, lr_threshold, delta, match_type, spm);
    }
    else
    {
    	int ngpu = m_devices_enabled.size();
    	Assert(ngpu >= 1);

    	if (db.size() < 2)
    	{
    		ngpu = 1;
    	}

    	switch (ngpu)
    	{
    	case 1:
    	{
    		// run on the first enabled card
            (void)n2_match_1gpu(matches, db, lr_threshold, delta, match_type, spm);
            break;
    	}
    	case 2:
    	case 3:
    		// treat 3 as 2 for now
            (void)n2_match_2gpus(matches, db, lr_threshold, delta, match_type, spm);
            break;
    	case 4:
    	default:
    		// if more than 4 just use 4 for now
            (void)n2_match_4gpus(matches, db, lr_threshold, delta, match_type, spm);
            break;
    	}
    }

    return matches.size();
}

int
GPUDevices::nm_match(
        NMmatchResults &matches,
        ProfileRange & dbA,
        ProfileRange & dbB,
        double lr_threshold,
        double delta1,
        double delta2,
        MatchType const &match_type,
        Split split,
        SubPopModel const & spm)
{
    if (split != split_parallel /* || dbA.size() < 2 || dbB.size() < 2 */ )
    {
        // run on one card only (whichever is the default) in this (main) CPU thread
        (void)nm_match_1gpu(matches, dbA, dbB, lr_threshold, delta1, delta2, full, match_type, spm);
    }
    else
    {
    	int ngpu = m_devices_enabled.size();
    	Assert(ngpu >= 1);

    	if (dbA.size() < 4 || dbB.size() < 4) // TODO use multiple GPUs for (few) vs (many)
    	{
    		ngpu = 1;
    	}

    	switch (ngpu)
    	{
    	case 1:
    	{
    		// run on the first enabled card in this (main) CPU thread
            Assert (!dbA.gotData());
            Assert (!dbB.gotData());
            (void)nm_match_1gpu(matches, dbA, dbB, lr_threshold, delta1, delta2, full, match_type, spm);
            break;
    	}
    	case 2:
    	case 3:
    		// treat 3 as 2 for now
            (void)nm_match_2gpus(matches, dbA,  dbB, lr_threshold, delta1, delta2, full, match_type, spm);
            break;

    	case 4:
    	default:
    		// if more than 4 just use 4 for now
            (void)nm_match_4gpus(matches, dbA,  dbB, lr_threshold, delta1, delta2, full, match_type, spm);
            break;
    	}

    }

    return matches.size();
}

int
GPUDevices::n2_match_2gpus(
        NMmatchResults &matches,
        ProfileRange &db,
        double lr_threshold,
        double delta,
        MatchType const &match_type,
        SubPopModel const &spm)
{
	// split into four tasks and launch two on each GPU

	int ntasks = 4;
	Assert(m_devices_enabled.size() >= 2);

    int size1 = db.size() / 2;

    ProfileRange dbA(db, db.begin(),         db.begin() + size1);
    ProfileRange dbB(db, db.begin() + size1, db.end());

	// calculate how much RAM is available for use by all devices
    MemoryLimits memlims = theHostManager().getLimits();
    size_t cpu_ram_for_devices = memlims.total_app * 1024;

    pthread_t threadID[ntasks];
    NMmatchResults empty_result = NMmatchResults(0);
    vector<NMmatchResults> tile(ntasks, empty_result);
    N2ThreadData data[ntasks];

	set<int>::const_iterator it = m_devices_enabled.begin();
	int gpu1 = *it;
	int gpu2 = *++it;

    for(int i = 0; i < ntasks; i++)
    {
        data[i].task         = i;
        data[i].device       = (i % 2) ? gpu1 : gpu2;
        data[i].dbA          = dbA;
        data[i].dbB          = dbB;
        data[i].matches      = &tile[i];
        data[i].lr_threshold = lr_threshold;
        data[i].delta        = delta;
        data[i].spm          = spm;
        data[i].match_type   = match_type;
        data[i].cpu_ram      = size_t(cpu_ram_for_devices * 0.5);

        pthread_create(&threadID[i],
                       NULL,
                       (FUNC)threadN2Split_4tasks,
                       &data[i]);

//        cout << "thread " << i << " launched" << endl;
    }

    for(int i = 0; i < ntasks; i++)
    {
        pthread_join(threadID[i], NULL);
//        cout << "thread " << i << " joined" << endl;
    }

    // combine results from each thread
    n2matchcpy(matches, tile[0], 0, 0);
    n2matchcpy(matches, tile[1], 0, size1);
    n2matchcpy(matches, tile[2], 0, size1);
    n2matchcpy(matches, tile[3], size1, size1);

    return matches.size();
}

int
GPUDevices::n2_match_4gpus(
        NMmatchResults &matches,
        ProfileRange &db,
        double lr_threshold,
        double delta,
        MatchType const &match_type,
        SubPopModel const &spm)
{
    //
	// Split over 4 cards
    // launch a separate host thread for each card
    //
	int ngpu = 4;
	Assert(m_devices_enabled.size() >= ngpu);

    int size1 = db.size() / 2;
	ProfileRange dbA(db, db.begin(),         db.begin() + size1);
	ProfileRange dbB(db, db.begin() + size1, db.end());

    MemoryLimits memlims = theHostManager().getLimits();
    size_t cpu_ram_for_devices = memlims.total_app * 1024;

    pthread_t threadID[ngpu];
    NMmatchResults empty_result = NMmatchResults(0);
    vector<NMmatchResults> tile(ngpu, empty_result);
    N2ThreadData data[ngpu];

	set<int>::const_iterator it = m_devices_enabled.begin();

    for(int i = 0; i < ngpu; i++)
    {
        data[i].task         = i;
        data[i].device       = *it++;
        data[i].dbA          = dbA;
        data[i].dbB          = dbB;
        data[i].matches      = &tile[i];
        data[i].lr_threshold = lr_threshold;
        data[i].delta        = delta;
        data[i].spm          = spm;
        data[i].match_type   = match_type;
        data[i].cpu_ram      = size_t(cpu_ram_for_devices * 0.25);

        pthread_create(&threadID[i],
                       NULL,
                       (FUNC)threadN2Split_4tasks,
                       &data[i]);
    }

    for(int i = 0; i < ngpu; i++)
    {
        pthread_join(threadID[i], NULL);
    }

    // combine results from each thread
    n2matchcpy(matches, tile[0], 0, 0);
    n2matchcpy(matches, tile[1], 0, size1);
    n2matchcpy(matches, tile[2], 0, size1);
    n2matchcpy(matches, tile[3], size1, size1);

    return matches.size();
}

int
GPUDevices::nm_match_2gpus(
        NMmatchResults &matches,
        ProfileRange & dbA,
        ProfileRange & dbB,
        double lr_threshold,
        double delta1,
        double delta2,
        ArrayPart part,
        MatchType const &match_type,
        SubPopModel const & spm)
{
	// split into 4 tasks and run on 2 GPUs
	int ntasks = 4;
	Assert(m_devices_enabled.size() >= 2);

	int size1 = dbA.size()/2;
	int size2 = dbB.size()/2;

	ProfileRange dbA1(dbA, dbA.begin(),         dbA.begin() + size1);
	ProfileRange dbA2(dbA, dbA.begin() + size1, dbA.end());
	ProfileRange dbB1(dbB, dbB.begin(),         dbB.begin() + size2);
	ProfileRange dbB2(dbB, dbB.begin() + size2, dbB.end());

    MemoryLimits memlims = theHostManager().getLimits();
    size_t cpu_ram_for_devices = memlims.total_app * 1024;

	pthread_t threadID[ntasks];
	NMmatchResults empty_result = NMmatchResults(0);
	vector<NMmatchResults> tile(ntasks, empty_result);
	NMThreadData data[ntasks];

	set<int>::const_iterator it = m_devices_enabled.begin();
	int gpu1 = *it;
	int gpu2 = *++it;

    for(int i = 0; i < ntasks; i++)
    {
		data[i].task         = i;
        data[i].device       = (i % 2) ? gpu1 : gpu2;
		data[i].dbA1         = dbA1;
		data[i].dbA2         = dbA2;
		data[i].dbB1         = dbB1;
		data[i].dbB2         = dbB2;
		data[i].matches      = &tile[i];
		data[i].lr_threshold = lr_threshold;
		data[i].delta1       = delta1;
		data[i].delta2       = delta2;
		data[i].spm          = spm;
		data[i].match_type   = match_type;
        data[i].cpu_ram      = size_t(cpu_ram_for_devices * 0.25); // 0.5?

        pthread_create(&threadID[i],
					   NULL,
					   (FUNC)threadNMSplit_4tasks,
					   &data[i]);
	}

	for(int i = 0; i < ntasks; i++)
	{
		pthread_join(threadID[i], NULL);
	}

	// combine
	n2matchcpy(matches, tile[0], 0, 0);
	n2matchcpy(matches, tile[1], 0, size2);
	n2matchcpy(matches, tile[2], size1, 0);
	n2matchcpy(matches, tile[3], size1, size2);

    return matches.size();
}

int
GPUDevices::nm_match_4gpus(
        NMmatchResults  &matches,
        ProfileRange & dbA,
        ProfileRange & dbB,
        double         lr_threshold,
        double         delta1,
        double         delta2,
        ArrayPart      part,
        MatchType      const &match_type,
        SubPopModel const & spm)
{
    //
	// Split over 4 cards
	// launch a separate host thread for each card
	//
	int ngpu = 4;
	Assert(m_devices_enabled.size() >= ngpu);

	int size1 = dbA.size()/2;
	int size2 = dbB.size()/2;

	ProfileRange dbA1(dbA, dbA.begin(),         dbA.begin() + size1);
	ProfileRange dbA2(dbA, dbA.begin() + size1, dbA.end());
	ProfileRange dbB1(dbB, dbB.begin(),         dbB.begin() + size2);
	ProfileRange dbB2(dbB, dbB.begin() + size2, dbB.end());

    MemoryLimits memlims = theHostManager().getLimits();
    size_t cpu_ram_for_devices = memlims.total_app * 1024;

	pthread_t threadID[ngpu];
	NMmatchResults empty_result = NMmatchResults(0);
	vector<NMmatchResults> tile(ngpu, empty_result);

	NMThreadData data[ngpu];

	set<int>::const_iterator it = m_devices_enabled.begin();

	for(int i = 0; i < ngpu; i++)
	{
		data[i].task         = i;
		data[i].device       = *it++;
		data[i].dbA1         = dbA1;
		data[i].dbA2         = dbA2;
		data[i].dbB1         = dbB1;
		data[i].dbB2         = dbB2;
		data[i].matches      = &tile[i];
		data[i].lr_threshold = lr_threshold;
		data[i].delta1       = delta1;
		data[i].delta2       = delta2;
		data[i].spm          = spm;
		data[i].match_type   = match_type;
        data[i].cpu_ram      = size_t(cpu_ram_for_devices * 0.25);

		pthread_create(&threadID[i],
					   NULL,
					   (FUNC)threadNMSplit_4tasks,
					   &data[i]);
	}

	for(int i = 0; i < ngpu; i++)
	{
		pthread_join(threadID[i], NULL);
	}

	// combine
	n2matchcpy(matches, tile[0], 0, 0);
	n2matchcpy(matches, tile[1], 0, size2);
	n2matchcpy(matches, tile[2], size1, 0);
	n2matchcpy(matches, tile[3], size1, size2);

    return matches.size();
}

// run n2 match on one GPU (GPU id is set with setDevice() before this function is called)
int
GPUDevices::n2_match_1gpu(
        NMmatchResults &matches,
        ProfileRange & db,
        double lr_threshold,
        double delta,
        MatchType const &match_type,
        SubPopModel const &spm,
		ThreadData const &data)
{
    struct N2Results results;
    results.count = 0;

    // run on one card only
    if ( ! setupTask(db.size(), db.size(), data) )
    {
    	return 0;
    }

    int nresults = cuda_match_n2(db, results, lr_threshold, delta, match_type, spm);
	cleanupTask();

    for (int j=0; j<nresults; ++j)
    {
        matches[ make_pair(results.index1[j], results.index2[j]) ].push_back(Result(match_type, results.lr[j]));
    }

    return matches.size();
}

// run NM match on one GPU (GPU id is set with setDevice() before this function is called)
int
GPUDevices::nm_match_1gpu(
        NMmatchResults &matches,
        ProfileRange & db1,
        ProfileRange & db2,
        double lr_threshold,
        double delta1,
        double delta2,
        ArrayPart part,
        MatchType const &match_type,
        SubPopModel const &spm,
		ThreadData const &data)
{
    struct N2Results results;
    results.count = 0;

    // do the match with the largest database as db1 // TODO WHY? try it the other way?
    bool swap_dbs = (db1.size() < db2.size());

    int nresults = 0;
    if (!swap_dbs)
    {
	    if ( ! setupTask(db1.size(), db2.size(), data) )
	    {
	    	return 0;
	    }

        Assert (!db1.gotData());
        Assert (!db2.gotData());
		nresults = cuda_match_nm(db1, db2, results, lr_threshold, delta1, delta2, match_type, spm, part); // in cuda_match.cpp
    }
    else
    {
	    if ( ! setupTask(db2.size(), db1.size(), data) )
	    {
	    	return 0;
	    }

		nresults = cuda_match_nm(db2, db1, results, lr_threshold, delta2, delta1, match_type, spm, part);
    }

	cleanupTask();

	for (int j=0; j<nresults; ++j)
    {
	    if (!swap_dbs)
	    {
	    	matches[ make_pair(results.index1[j], results.index2[j]) ].push_back(Result(match_type, results.lr[j]));
	    }
	    else
	    {
	    	matches[ make_pair(results.index2[j], results.index1[j]) ].push_back(Result(match_type, results.lr[j]));
	    }
    }

	return matches.size();
}

// split an N2 match over 4 devices
// TODO split over all enabled devices
void
threadN2Split_4tasks(N2ThreadData *data)
{
    info << startl << "threadN2Split_4tasks(): device = " << data->device << endl;

    if (data->task == 0)
    {
        // dbA/dbA upper match
        info << startl << "*** Split : AA" << endl;
        (void) gpuDevices().n2_match_1gpu(*(data->matches), data->dbA, data->lr_threshold, data->delta, data->match_type, data->spm, *data);
    }

    // TODO would be more efficient in tasks 1 and 2 to match half of A against all of B
    if (data->task == 1)
    {
        // dbA/dbB upper match
        info << startl << "*** Split : AB" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA, data->dbB, data->lr_threshold, data->delta, data->delta, upper, data->match_type, data->spm, *data);
    }

    if (data->task == 2)
    {
        // dbA/dbB lower match
        info << startl << "*** Split : BA" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA, data->dbB, data->lr_threshold, data->delta, data->delta, lower, data->match_type, data->spm, *data);
    }

    if (data->task == 3)
    {
        // dbB/dbB upper match
        info << startl << "*** Split : BB" << endl;
        (void) gpuDevices().n2_match_1gpu(*(data->matches), data->dbB, data->lr_threshold, data->delta, data->match_type, data->spm, *data);
    }
}

// split an NM match over 4 devices
void
threadNMSplit_4tasks(NMThreadData *data)
{
    info << startl << "threadNMSplit_4tasks(): device = " << data->device << endl;

    if (data->task == 0)
    {
        // dbA/dbA upper match
        info << startl << "*** Split : A1.B1" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA1, data->dbB1, data->lr_threshold, data->delta1, data->delta2, full, data->match_type, data->spm, *data);
    }

    if (data->task == 1)
    {
        // dbA/dbB upper match
        info << startl << "*** Split : A1.B2" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA1, data->dbB2, data->lr_threshold, data->delta1, data->delta2, full, data->match_type, data->spm, *data);
    }

    if (data->task == 2)
    {
        // dbA/dbB lower match
        info << startl << "*** Split : A2.B1" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA2, data->dbB1, data->lr_threshold, data->delta1, data->delta2, full, data->match_type, data->spm, *data);
    }

    if (data->task == 3)
    {
        // dbB/dbB upper match
        info << startl << "*** Split : A2.B2" << endl;
        (void) gpuDevices().nm_match_1gpu(*(data->matches), data->dbA2, data->dbB2, data->lr_threshold, data->delta1, data->delta2, full, data->match_type, data->spm, *data);
    }
}

// copy results from tile into matches, at the given offset
void
GPUDevices::n2matchcpy(NMmatchResults &matches, NMmatchResults &tile, int istart, int jstart)
{
    for (NMmatchResults::const_iterator it = tile.begin(); it != tile.end(); ++it)
    {
        Assert(it->second.size() == 1);
        matches[make_pair(it->first.first + istart, it->first.second + jstart)].push_back(it->second[0]);
    }
}

// Link this thread to a task.
// NB This must be called *within* the thread to ensure thread, task and device are linked
void
GPUDevices::initThread(ThreadData const &data)
{
//	cout << "initThread" << endl;

	Assert2(data.task >= 0 && data.task < MAX_CPU_THREADS,
			"initThread: Bad thread index");

	m_task = data.task;
	setDevice(data.device);

	taskData[m_task].device_no = data.device;
	// ...
}

void
GPUDevices::cleanupTask()
{
//	cout << "cleanupTask" << endl;

	Assert2(cudaThreadSynchronize() == cudaSuccess, "cleanupThread: cudaThreadSynchronize failed");
	Assert2(cudaThreadExit() == cudaSuccess, "cleanupThread: cudaThreadExit failed");
}

// Return task data for this thread.
// NB must be called *within* the thread *after* initThread()
TaskData &
GPUDevices::getTaskData()
{
	Assert2(m_task >= 0 && m_task < MAX_CPU_THREADS,
			"getTaskData: Bad thread index");
	return taskData[m_task];
}

int
GPUDevices::getTask()
{
	Assert2(m_task >= 0 && m_task < MAX_CPU_THREADS,
			"getIndex: Bad thread index");
	return m_task;
}

int
GPUDevices::getDev()
{
	Assert2(m_task >= 0 && m_task < MAX_CPU_THREADS,
			"getDev: Bad thread index");
	return taskData[m_task].device_no;
}

// Calculate the parameters for running a task on a GPU
void
fitToDevice(
		float f,
		size_t  device_profs,
		size_t  &db1_chunk,
		size_t  &db2_bigchunk,
		size_t  &db2_smallchunk)
{
	db1_chunk      = size_t(f * DB1_CHUNK_F * device_profs);
	db2_smallchunk = size_t(f * (1-DB1_CHUNK_F)/2.0 * device_profs);
	db2_bigchunk   = db2_smallchunk * DB2_CHUNK_RATIO;
}

bool
GPUDevices::setupTask(
		int 			      db1_size,
		int 			      db2_size,
		ThreadData     const &data,
		GPUDevices     const &gpu_devs,
		PopulationData const &popdata)
{
	int task = data.task;

	size_t &db1_chunk        = taskData[task].db1_chunk;
	size_t &db2_bigchunk     = taskData[task].db2_bigchunk;
	size_t &db2_smallchunk   = taskData[task].db2_smallchunk;

	size_t bytes_per_profile = popdata.getLocusInfo().profile_size * sizeof(float);

	long cpu_ram = data.cpu_ram; // bytes
	size_t device_ram = size_t(gpu_devs.totalGlobalMem(data.device) * DEV_F); // bytes
	size_t device_profs = device_ram / bytes_per_profile;

	if (cpu_ram == 0)
	{
		// memory calculation has not been done
		// do it now on the assumption this GPU can use all the pinned memory for this thread
		cpu_ram = (long)theHostManager().getLimits().total_app * 1024;
	}

	// for data read from file: Need BIN_PROF_K for each profile

	long pinned_ram = device_ram; // no point having more pinned RAM than device RAM

	int pinned_max = min((long)(MAX_PINNED_PROFILES * bytes_per_profile), (long)(PINNED_F * cpu_ram));
	if (pinned_ram > pinned_max) // but not too much
	{
		pinned_ram = pinned_max;
	}

	// make sure we have a minimal amount of pinned RAM available
	if (pinned_ram < (long)(MIN_PINNED_PROFILES * bytes_per_profile) )
	{
		error << startl << "Not enough RAM for these datasets: " << db1_size << ", " << db2_size << " Profiles" << endl;
		return false;
	}

	size_t pinned_profs = pinned_ram / bytes_per_profile;

	// first, try putting all the data on the device at once
	db1_chunk      = db1_size;
	db2_bigchunk   = db2_size;
	db2_smallchunk = db2_bigchunk;

	// see if we overflow device memory
	if (2 * db1_chunk + 2 * db2_smallchunk > device_profs)
	{
		// try using smaller chunks of DB2
		db1_chunk      = db1_size;
		db2_bigchunk   = db2_size;
		db2_smallchunk = db2_bigchunk / DB2_CHUNK_RATIO;
	}

	// see if we still overflow device memory
	if (2 * db1_chunk + 2 * db2_smallchunk > device_profs)
	{
		// ok, scale to the size of device memory
		fitToDevice(1, device_profs, db1_chunk, db2_bigchunk, db2_smallchunk);

		db1_chunk /= 2;

		db2_bigchunk = min(db2_bigchunk, (size_t)db2_size);
		db2_smallchunk = min(db2_smallchunk, (size_t)db2_size);

		Assert(2 * db1_chunk + 2 * db2_smallchunk <= device_profs);
	}

	// see if we overflow CPU pinned memory limit
	if (db1_chunk + db2_bigchunk > pinned_profs)
	{
		// scale the allocation down to fit
		float f = (0.99 * pinned_profs) / (device_profs * (DB1_CHUNK_F + 0.5*DB2_CHUNK_RATIO*(1 - DB1_CHUNK_F)));
		fitToDevice(f, device_profs, db1_chunk, db2_bigchunk, db2_smallchunk);

		db2_bigchunk = min(db2_bigchunk, (size_t)db2_size);
		db2_smallchunk = min(db2_smallchunk, (size_t)db2_size);

		Assert(2 * db1_chunk + 2 * db2_smallchunk <= device_profs);
		Assert(db1_chunk + db2_bigchunk <= pinned_profs);
	}

    info << startl << "setupTask()" << endl;
    info << alignl << "    task              = " << task << endl;
    info << alignl << "    device            = " << data.device << endl;
    info << alignl << "    cpu_ram (B)       = " << cpu_ram << endl;
    info << alignl << "    pinned_ram (B)    = " << pinned_ram << endl;
    info << alignl << "    device_ram (B)    = " << device_ram << endl;
    info << alignl << "    bytes_per_profile = " << bytes_per_profile << endl;
    info << alignl << "    db1_size          = " << db1_size << endl;
    info << alignl << "    db2_size          = " << db2_size << endl;
    info << alignl << "    db1_chunk         = " << db1_chunk << " (" << db1_chunk * bytes_per_profile << ")" << endl;
    info << alignl << "    db2_bigchunk      = " << db2_bigchunk << " (" << db2_bigchunk * bytes_per_profile << ")" << endl;
    info << alignl << "    db2_smallchunk    = " << db2_smallchunk << " (" << db2_smallchunk * bytes_per_profile << ")" << endl;

	initThread(data);

	return true;
}

TEST(CUDAProperties)
{
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);

	if (deviceCount != 0)
	{

		for (int dev = 0; dev < deviceCount; ++dev)
		{
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);

			string name(deviceProp.name);
			if (name.find("Emulation") == string::npos)
			{
				cout << "device " << dev << " name " << setw(16) << deviceProp.name
										 << " totalGlobalMem = " << setw(12) << deviceProp.totalGlobalMem
										 << " display? " << (deviceProp.kernelExecTimeoutEnabled ? "Y" : "N")
										 << " deviceOverlap? " << (deviceProp.deviceOverlap ? "Y" : "N")
										 << " canMapHostMemory? " << (deviceProp.canMapHostMemory ? "Y" : "N") << endl;
				CHECK(deviceProp.deviceOverlap);
	//			CHECK(deviceProp.canMapHostMemory); // we don't use this feature
			}
		}
	}
}

#if 0 // disabled while we tweak the algorithm

TEST(setupTask)
{
	PopulationData popdata = testPopulationData();

	int bytes_per_profile = popdata.getLocusInfo().profile_size * sizeof(float);
	CHECK_EQUAL(76, bytes_per_profile); // = (10+3+3+3)*4 = 76

	int task     = 0;
	int dev      = 0;
	int db1_size = 1000; // Profiles
	int db2_size = 1000; // Profiles


	// plenty of RAM
	{
		// CPU RAM available as pinned for this device (B)
		size_t cpu_ram  = (db1_size + db2_size) * bytes_per_profile;

	    // Total global memory on device (B)
		// NB only a fraction DEV_F of this is available for Profiles, so double the required amount
		size_t gpu_ram = 2 * (db1_size + db2_size) * bytes_per_profile;

		ThreadData data(task, dev, cpu_ram);
		GPUDevices myGPUDevices(gpu_ram);
		GPUDevices::setupTask(db1_size, db2_size, data, myGPUDevices, popdata);

		size_t expected_db1_chunk      = db1_size;
		size_t expected_db2_smallchunk = db2_size;
		size_t expected_db2_bigchunk   = db2_size;

		CHECK_EQUAL(expected_db1_chunk,      taskData[0].db1_chunk);
		CHECK_EQUAL(expected_db2_bigchunk,   taskData[0].db2_bigchunk);
		CHECK_EQUAL(expected_db2_smallchunk, taskData[0].db2_smallchunk);

		// check we fit in GPU RAM
		CHECK(expected_db1_chunk + 2*expected_db2_smallchunk <= gpu_ram * DEV_F);

		// check we fit in CPU RAM
		CHECK(expected_db1_chunk + expected_db2_bigchunk <= cpu_ram);

	}

	// not quite enough GPU RAM
	{
		// CPU RAM available as pinned for this device (B)
		size_t cpu_ram  = (db1_size + db2_size) * bytes_per_profile;

	    // Total global memory on device (B)
		// NB only a fraction DEV_F of this is available for Profiles, this is not quite enough to
		// put all the profiles on the device, hence a SMALL_CHUNK will be used
		size_t gpu_ram = (db1_size + db2_size) * bytes_per_profile;

		ThreadData data(task, dev, cpu_ram);
		GPUDevices myGPUDevices(gpu_ram);
		GPUDevices::setupTask(db1_size, db2_size, data, myGPUDevices, popdata);

		CHECK_EQUAL(db1_size, taskData[0].db1_chunk);
		CHECK_EQUAL(db2_size, taskData[0].db2_bigchunk);
		CHECK_EQUAL(db2_size / DB2_CHUNK_RATIO, taskData[0].db2_smallchunk);
	}

	// nowhere near enough GPU RAM
	{
		// CPU RAM available as pinned for this device (B)
		size_t cpu_ram  = (db1_size + db2_size) * bytes_per_profile;

	    // Total global memory on device (B)
		// NB only a fraction DEV_F of this is available for Profiles, this is not quite half
		// what is needed to put all the profiles on the device.
		size_t gpu_ram = (db1_size + db2_size) * bytes_per_profile / 2;

		ThreadData data(task, dev, cpu_ram);
		GPUDevices myGPUDevices(gpu_ram);
		GPUDevices::setupTask(db1_size, db2_size, data, myGPUDevices, popdata);

		float profiles_on_device     = size_t(gpu_ram * DEV_F / bytes_per_profile);
		size_t expected_db1_chunk      = size_t(DB1_CHUNK_F * profiles_on_device);
		size_t expected_db2_smallchunk = size_t(0.5 * (1 - DB1_CHUNK_F) * profiles_on_device);
		size_t expected_db2_bigchunk   = expected_db2_smallchunk * DB2_CHUNK_RATIO;

		CHECK_EQUAL(expected_db1_chunk,      taskData[0].db1_chunk);
		CHECK_EQUAL(expected_db2_bigchunk,   taskData[0].db2_bigchunk);
		CHECK_EQUAL(expected_db2_smallchunk, taskData[0].db2_smallchunk);

		// check we fit in GPU RAM
		CHECK(expected_db1_chunk + 2*expected_db2_smallchunk <= gpu_ram * DEV_F);

		// check we fit in CPU RAM
		CHECK(expected_db1_chunk + expected_db2_bigchunk <= cpu_ram);
	}

	// enough GPU RAM, not enough CPU RAM
	{
		// CPU RAM available as pinned for this device (B)
		// Only half what we need to put all the Profiles in pinned memory at once
		int cpu_ram  = (db1_size + db2_size) * bytes_per_profile / 2;

	    // Total global memory on device (B)
		// NB only a fraction DEV_F of this is available for Profiles, so double the required amount
		size_t gpu_ram = 2 * (db1_size + db2_size) * bytes_per_profile;

		ThreadData data(task, dev, cpu_ram);
		GPUDevices myGPUDevices(gpu_ram);
		GPUDevices::setupTask(db1_size, db2_size, data, myGPUDevices, popdata);

		// first calculate what would fit nicely on the GPU (as previous test)
		float profiles_on_device     = size_t(gpu_ram * DEV_F / bytes_per_profile);
		size_t expected_db1_chunk      = size_t(DB1_CHUNK_F * profiles_on_device);
		size_t expected_db2_smallchunk = size_t(0.5 * (1 - DB1_CHUNK_F) * profiles_on_device);
		size_t expected_db2_bigchunk   = expected_db2_smallchunk * DB2_CHUNK_RATIO;

		// we overshoot CPU RAM by this fraction
		float f = (cpu_ram/bytes_per_profile) / float(expected_db1_chunk + expected_db2_bigchunk);
		CHECK(f < 1);

		// so scale everything back by a factor f
		expected_db1_chunk       = size_t(expected_db1_chunk * f);
		expected_db2_bigchunk    = size_t(expected_db2_bigchunk * f);
		expected_db2_smallchunk  = size_t(expected_db2_smallchunk * f);

		// This is an overestimate, due to safety factors and rounding
		CHECK(expected_db1_chunk      >= taskData[0].db1_chunk);
		CHECK(expected_db2_bigchunk   >= taskData[0].db2_bigchunk);
		CHECK(expected_db2_smallchunk >= taskData[0].db2_smallchunk);

		// Should not be over by more than 5%
		CHECK(0.95 * expected_db1_chunk      < taskData[0].db1_chunk);
		CHECK(0.95 * expected_db2_bigchunk   < taskData[0].db2_bigchunk);
		CHECK(0.95 * expected_db2_smallchunk < taskData[0].db2_smallchunk);

		// check we fit in GPU RAM
		CHECK(expected_db1_chunk + 2*expected_db2_smallchunk <= gpu_ram * DEV_F);

		// check we fit in CPU RAM
		CHECK(expected_db1_chunk + expected_db2_bigchunk <= cpu_ram);

	}

	// not enough of either
	{
		// CPU RAM available as pinned for this device (B)
		// Only half what we need to put all the Profiles in pinned memory at once
		int cpu_ram  = (db1_size + db2_size) * bytes_per_profile / 2;

	    // Total global memory on device (B)
		// NB only a fraction DEV_F of this is available for Profiles
		size_t gpu_ram = (db1_size + db2_size) * bytes_per_profile / 4;

		ThreadData data(task, dev, cpu_ram);
		GPUDevices myGPUDevices(gpu_ram);
		GPUDevices::setupTask(db1_size, db2_size, data, myGPUDevices, popdata);

		// first calculate what would fit nicely on the GPU (as previous test)
		float profiles_on_device     = size_t(gpu_ram * DEV_F / bytes_per_profile);
		size_t expected_db1_chunk      = size_t(DB1_CHUNK_F * profiles_on_device);
		size_t expected_db2_smallchunk = size_t(0.5 * (1 - DB1_CHUNK_F) * profiles_on_device);
		size_t expected_db2_bigchunk   = expected_db2_smallchunk * DB2_CHUNK_RATIO;

		// we overshoot CPU RAM by this fraction
		float f = (cpu_ram/bytes_per_profile) / float(expected_db1_chunk + expected_db2_bigchunk);

		if (f < 1)
		{
			// scale everything back by a factor f
			expected_db1_chunk       = size_t(expected_db1_chunk * f);
			expected_db2_bigchunk    = size_t(expected_db2_bigchunk * f);
			expected_db2_smallchunk  = size_t(expected_db2_smallchunk * f);

			// This is an overestimate, due to safety factors and rounding
			CHECK(expected_db1_chunk      >= taskData[0].db1_chunk);
			CHECK(expected_db2_bigchunk   >= taskData[0].db2_bigchunk);
			CHECK(expected_db2_smallchunk >= taskData[0].db2_smallchunk);

			// Should not be over by more than 5%
			CHECK(0.95 * expected_db1_chunk      < taskData[0].db1_chunk);
			CHECK(0.95 * expected_db2_bigchunk   < taskData[0].db2_bigchunk);
			CHECK(0.95 * expected_db2_smallchunk < taskData[0].db2_smallchunk);
		}
		else
		{
			CHECK_EQUAL(expected_db1_chunk,      taskData[0].db1_chunk);
			CHECK_EQUAL(expected_db2_bigchunk,   taskData[0].db2_bigchunk);
			CHECK_EQUAL(expected_db2_smallchunk, taskData[0].db2_smallchunk);
		}

		// check we fit in GPU RAM
		CHECK(expected_db1_chunk + 2*expected_db2_smallchunk <= gpu_ram * DEV_F);

		// check we fit in CPU RAM
		CHECK(expected_db1_chunk + expected_db2_bigchunk <= cpu_ram);

	}
}
#endif

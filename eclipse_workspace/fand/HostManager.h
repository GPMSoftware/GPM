/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * HostManager.h
 *
 *  Created on: Apr 15, 2011
 *      Author: gareth
 */

#ifndef HOSTMANAGER_H_
#define HOSTMANAGER_H_

#include "Cached.h"

#include <fstream>

static const float APP_F = 0.9;             // fraction of available CPU RAM to use for this app
static const float PINNED_F = 0.1;          // fraction of available RAM to use for pinned memory
static const size_t MIN_PINNED_PROFILES = 1000;   // minimum pinned profiles we can work with
static const size_t MAX_PINNED_PROFILES = 20000;  // max pinned profiles
static const float RAW_PROF_K = 20;         // RAM needed for a raw profile (KiB)
static const float BIN_PROF_K = 74;         // RAM needed for a binary profile (KiB)

struct MemoryLimits
{
    MemoryLimits(int total_app = 0 /*, int cache = 0, int pinned = 0*/ )
    : total_app(total_app) //, cache(cache), pinned(pinned)
    {}

    size_t total_app; // total memory for use by this application (kB)
//    size_t cache;     // limit for cache memory (kB)
//    size_t pinned;    // limit for pinned memory (kB)
};

// Manage the memory available to us on the host.
// * At the start of the program, work out how much memory we
//   have and allocate it between "cache", "pinned" and "other"
// * Update and return info about the memory
class HostManager
{
private:
    HostManager() {}
    friend HostManager & theHostManager();
public:

    struct MemInfo
    {
        size_t memTotal_k, memFree_k, buffers_k, cached_k;
    };

    MemoryLimits m_limits;

    // calculate the MemoryLimits for use by this process
    void init();

    // Return the MemoryLimits calculated by init()
    // (Must call init() before getLimits())
    MemoryLimits getLimits();

    // return the MemInfo structure
    MemInfo meminfo() const;

    // return the total physical memory on the system (k)
    size_t total();

    // return 'free' memory (according to linux) but see 'available' (k)
    size_t free();

    // return 'free' memory (according to linux) (k)
    size_t buffers();

    // return 'free' memory (according to linux) (k)
    size_t cached();

    // return memory available for programs. (k)
    // NB This is in fact free + buffers + cached, because buffers and cached are released
    // by the system when needed
    size_t available();
};

// return the HostManager singleton
HostManager & theHostManager();

#endif /* HOSTMANAGER_H_ */

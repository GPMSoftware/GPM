/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * HostManager.cpp
 *
 *  Created on: Apr 15, 2011
 *      Author: gareth
 */

#include "HostManager.h"

// return HostManager singleton
HostManager &
theHostManager()
{
    static HostManager the_host_manager;
    return the_host_manager;
}

void
HostManager::init()
{
    m_limits.total_app = size_t(available() * APP_F);    // leave some room for other apps.
//    m_limits.cache     = size_t(m_limits.total_app * CACHED_F);
//    m_limits.pinned    = size_t(m_limits.total_app * PINNED_F);

    // set the cache
//    (void)Cached::setLimit(m_limits.cache * 1024);
    (void)Cached::setLimit(1000000000); // Cached not currently used
}

// Must call init() before getLimits()
MemoryLimits
HostManager::getLimits()
{
    Assert(m_limits.total_app > 0);
    return m_limits;
}

HostManager::MemInfo
HostManager::meminfo() const
{
    MemInfo ret;

    std::ifstream f("/proc/meminfo"); // the memory in meminfo is in kB
    std::string dummy;
    f >> dummy >> ret.memTotal_k >> dummy;
    f >> dummy >> ret.memFree_k >> dummy;
    f >> dummy >> ret.buffers_k >> dummy;
    f >> dummy >> ret.cached_k >> dummy;

    return ret;
}

// return the total physical memory on the system (k)
size_t
HostManager::total()
{
    MemInfo mi = meminfo();
    return mi.memTotal_k;
}

size_t
HostManager::free()
{
    MemInfo mi = meminfo();
    return mi.memFree_k;
}

size_t
HostManager::buffers()
{
    MemInfo mi = meminfo();
    return mi.buffers_k;
}

size_t
HostManager::cached()
{
    MemInfo mi = meminfo();
    return mi.cached_k;
}

size_t
HostManager::available()
{
    MemInfo mi = meminfo();
    return mi.memFree_k + mi.buffers_k + mi.cached_k;
}

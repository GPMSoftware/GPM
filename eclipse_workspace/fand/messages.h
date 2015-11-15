/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * messages.h
 *
 *  Created on: Mar 25, 2010
 *      Author: gareth
 */

#ifndef MESSAGES_H_
#define MESSAGES_H_

#define warn  warn_f()
#define error error_f()
#define info  info_f()
#define info2 info2_f()
#define ok    ok_f()

static MessageStream& error_f()
{
    static MessageStream theError(thisid, thisfile, "error", "ERROR:   ", std::cout, "error_log.txt", true);  // to cerr???, logged, enabled by default
    return theError;
}

static MessageStream& warn_f()
{
    static MessageStream theWarning (thisid, thisfile, "warn",  "WARNING: ", std::cout); // to cout, not logged, enabled by default
    return theWarning;
}

static MessageStream& info_f()
{
    static MessageStream theInfo(thisid, thisfile, "info",  "INFO:    ", std::cout, "", false); // to cout, not logged, not enabled by default
    return theInfo;
}

static MessageStream& info2_f()
{
    static MessageStream theInfo2(thisid, thisfile, "info2", "INFO2:   ", std::cout, "", false); // to cout, not logged, not enabled by default
    return theInfo2;
}

static MessageStream& ok_f()
{
    static MessageStream theOk(thisid,    thisfile, "ok",    "OK:      ", std::cout, "", false); // to cout, not logged, not enabled by default
    return theOk;
}

#endif /* MESSAGES_H_ */

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * gmatch.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef GMATCH_H_
#define GMATCH_H_

#include "wx/wx.h"
#include <vector>
#include <set>
#include <string>

#include "fand/Database.h"
#include "fand/dbmatch.h"
#include "fand/Profile.h"

extern boost::shared_ptr<Database> db;

extern wxColour readOnlyGrey;

enum
{
    ID_QUIT = 1,
	ID_IMPORT,
    ID_EXPORT,
    ID_SINGLE,
    ID_MULTI,
    ID_PCOMBO,
    ID_DSLIST,
    ID_DISPLAY,
    ID_MATCH,
    ID_MATCHLIST,
    ID_RESULTLIST,
    ID_SAVE_RESULT,
    ID_CLEAR_RESULT,
    ID_META,
    ID_NEW,
    ID_EDIT,
    ID_DELETE,
    ID_CLEAR,
    ID_SAVE,
    ID_KIT,
    ID_GRID,
    ID_THRESH,
    ID_STRING0,
    ID_STRING1,
    ID_STRING2,
    ID_RELLIST,
    ID_CR,
    ID_BROWSE,
    ID_BROWSE_FDB,
    ID_CDELTA,
    ID_RDELTA,
    ID_FCOMBO,
    ID_VCOMBO,
    ID_VTEXT,
    ID_NEWB,
    ID_DELB,
    ID_LOGIC,
    ID_SRCH,
    ID_FILE,
    ID_FORMAT,
    ID_GPULIST,
    ID_CPU,
    ID_GPU,
    ID_HW,
    ID_NRC4_4,
    ID_NRC4_10,
    ID_THETA,
    ID_HOST,
    ID_USER,
    ID_PASS,
    ID_CHANGE,
    ID_KM,
    ID_MGRS,
    ID_SEARCH,
    ID_LIST,
    ID_PREV,
    ID_NEXT,
    ID_ADD,
    ID_REMOVE,
    ID_SCOMBO,
    ID_SAVE_SEARCH,
    ID_DEL_SEARCH,
    ID_MENU_NEW,
    ID_MENU_EDIT,
    ID_MENU_DELETE,
    ID_MENU_CLEAR,
    ID_MENU_SAVE,
    ID_MENU_BULK,
    ID_APPLY,
    ID_CHOICE,
    ID_NEWVAL,
    ID_NOMUT,
    ID_CRIME,
    ID_REF
};

void
setItemContainerStrings(wxItemContainer* ic, std::vector<std::string> const &s);

void
setItemContainerStrings(wxItemContainer* ic, std::set<std::string> const &s);

void
do_match(
    NMmatchResults          &matches,
    ProfileRange            &crime_db,
    ProfileRange            &ref_db,
    bool                    ref_is_crime,
    MatchParams     const & match_params);

void
nm_output(MatchResults &results,
          bool output_profs);

#endif /* GMATCH_H_ */

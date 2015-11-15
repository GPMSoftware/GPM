/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
//============================================================================
// Name        : gmatch.cpp
// Author      : DGW Software Consultants LTD
// Description : GUI Match Profiles/Profile Databases read from CSV files
//============================================================================

#include "gmatch.h"
#include "MyFrame.h"

#include "fand/dbmatch.h"
#include "fand/ProfileFilter.h"
#include "fand/ProfileData.h"
//#include "fand/populationdata.h"
#include "fand/loci.h"
#include "fand/Database.h"
#include "fand/Version.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("gmatch");
#include "fand/messages.h"

#include <UnitTest++/UnitTest++.h>

#include "wx/wx.h"

using namespace std;

#define SPLIT split_parallel
//#define SPLIT split_none

char *progname = "gmatch";

static char *db_name = "fand";  /* database name */

boost::shared_ptr<Database> db( new Database(db_name) );

wxColour readOnlyGrey(235, 235, 235);

class MyApp: public wxApp
{
    virtual bool OnInit();
//    virtual int FilterEvent(wxEvent& event);
};

IMPLEMENT_APP(MyApp)

#if 0
int MyApp::FilterEvent(wxEvent& event)
{
    if ( event.GetEventType()==wxEVT_KEY_DOWN)
    {
        int code = ((wxKeyEvent&)event).GetKeyCode();
        cout << "MyApp::FilterEvent: " << code << endl;

        wxKeyEvent *e = dynamic_cast<wxKeyEvent*>(&event);

        if (code == WXK_MULTIPLY)
        {
            e->m_keyCode = '@'; // This does not change the character that appears
        }
        else
        {
            // ...
            e->m_keyCode = 'X';
        }

    }

    return -1; // skip
}
#endif

bool runUnitTests()
{
	cout << "Running Unit Tests ..." << endl;
	MessageStream::enable(false); // suppress error/warning/etc during unit tests

	bool ret = (UnitTest::RunAllTests() == 0);

	MessageStream::enable(true);

	return ret;
}

bool MyApp::OnInit()
{
	if (!runUnitTests())
	{
		cout << "Unit tests failed - quitting" << endl;
		exit(1);
	}
	cout << "Unit tests completed successfully." << endl << endl;

    Assert(isMainThread());

    MyFrame *frame = new MyFrame( _(PROGRAM_NAME) + (" V") + _(MAJOR_VERSION) + _(".") + _(MINOR_VERSION), wxPoint(50, 50),
                                  wxSize(450,340) );
    frame->Show(true);
    SetTopWindow(frame);

    MyFrame::refreshMessageWindow();
    MyFrame::m_inputArea->emptyMode();

    return true;
}

void
setItemContainerStrings(wxItemContainer* ic, vector<string> const &s)
{
    ic->Clear();
    for (size_t i=0; i<s.size(); ++i)
    {
        wxString str(s[i].c_str(), *wxConvCurrent);
        ic->Append(str);
    }
}

void
setItemContainerStrings(wxItemContainer* ic, set<string> const &s)
{
    ic->Clear();
    set<string>::const_iterator it;
    for (it = s.begin(); it != s.end(); ++it)
    {
        wxString str(it->c_str(), *wxConvCurrent);
        ic->Append(str);
    }
}

void
do_match(
    NMmatchResults          &matches,
    ProfileRange            &crime_db,
    ProfileRange            &ref_db,
    bool                    ref_is_crime,
    MatchParams     const & match_params)
{
	do_gmatch(
		matches,
		error,        // NB output goes to 'error' to make it appear in the GUI window
		SPLIT,
		crime_db,
		ref_db,
		ref_is_crime,
		match_params.crime_delta,
		match_params.ref_delta,
		match_params.match_types,
		match_params.spm,
		match_params.mut,
		match_params.lr_threshold,
		match_params.n_cuda_accel,
		match_params.n2_cuda_accel);
}

const string&
locusNames()
{
    static string ret;

    if (ret.size() == 0)
    {
        ostringstream ss;

        ss << setw(24) << "Allele";
        for (int locus = 0; locus < num_loci; ++locus)
        {
            ss << setw(9) << locus_name[locus] /*<< setw(10) << locus_name[locus]*/ ;
        }

        ret = ss.str();
    }

    return ret;
}

// temporary
// direct output to 'error' to make it appear in the GUI window
// use the indexes c_ids, r_ids to get the profiles efficiently
void
nm_output(MatchResults &results,
          bool output_profs)
{
	NMmatchResults const &matches = results.results;

    int n = 0; // match count
    for (NMmatchResults::const_iterator it = matches.begin(); it != matches.end(); ++it)
    {
        ++n;
        int m = 0; // relationship count
        int i = it->first.first;
        int j = it->first.second;

        std::vector<Profile> pvec;
        std::map<int, std::string>::const_iterator it2 = results.c_ids.find(i);
        Assert2(it2 != results.c_ids.end(), "nm_output(): can't read profile");
		dbReadProfile(*db, it2->second, pvec);

		it2 = results.r_ids.find(j);
        Assert2(it2 != results.r_ids.end(), "nm_output(): can't read profile");
		dbReadProfile(*db, it2->second, pvec);

		const Profile &p1 = pvec[0], &p2 = pvec[1];

        for (vector<Result>::const_iterator rit = it->second.begin(); rit != it->second.end(); ++rit)
        {
            double lr = rit->likelihood_ratio;

            error << "match " << n << "." << ++m << ": " << setw(9) << p1.data().m_id << ", " << setw(9) << p2.data().m_id;
            error << ": " << setw(30) << rit->match_type << ": Likelihood ratio = " << setprecision(6) << lr;
            error << endl;
        }
        if (output_profs)
        {
            error << locusNames() << endl << p1 << endl << p2 << endl << endl;
        }
    }
}

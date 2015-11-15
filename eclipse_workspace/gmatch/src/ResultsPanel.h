/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ResultsPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef RESULTSPANEL_H_
#define RESULTSPANEL_H_

#include "fand/dbmatch.h"
#include "fand/Profile.h"
#include "MyGrid.h"

#include "wx/wx.h"

#include <vector>

static const int SIZE = 1024;

struct ResultDisplayItem
{
    int c_index;    // index into C-Set
    int r_index;    // index into R-Set
    double max_lr;  // max LR of any result

    // representation as C-ID | R-ID | max_lr
    std::string c_string;
    std::string r_string;
    std::string lr_string;
};

class MatchGrid :public MyGrid
{
public:

	MatchGrid( wxWindow *parent,
        wxWindowID id,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxWANTS_CHARS,
        const wxString& name = wxGridNameStr )
    : MyGrid(parent, id, pos, size, style, name)
    {
		setup();
    }

	void setup();
    void OnChar(wxKeyEvent& event);
    int GetSelection();
    void setMatches(std::vector<ResultDisplayItem> const &s);
    void setMessage(std::string const &msg);
};

class ResultGrid :public MyGrid
{
public:

	ResultGrid( wxWindow *parent,
        wxWindowID id,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxWANTS_CHARS,
        const wxString& name = wxGridNameStr )
    : MyGrid(parent, id, pos, size, style, name)
    {
		setup();
    }

	void setup();
    void setResults(std::vector<std::pair<std::string, std::string> > const &s);
};

class ResultsPanel : public wxPanel
{
public:
    // Constructor
    ResultsPanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
    , m_matchList(0)
    , m_resultsList(0)
    {
        setup(title);
    }

    void OnSelectMatch(wxGridRangeSelectEvent& event);
    void OnSave(wxCommandEvent& event);
    void OnClear(wxCommandEvent& event);
    void OnFormat(wxCommandEvent& event);

    void displayMatches(
        ProfileRange const         &cset,
        ProfileRange const         &rset,
        std::vector<MatchType>     &match_types,
        float                       lr,
        NMmatchResults       const &matches);

    void writeMatches(
            std::ofstream              &ofs,
            NMmatchResults const &matches);

    void DisplayResults(bool sci);

    void SelectResult();

    void clear();

    wxRadioBox   *m_format;         // display format (sci/fixed)
    MatchGrid    *m_matchList;      // List of pairs of profiles
    ResultGrid   *m_resultsList;    // List of results for the selected pair of profiles
    wxButton     *m_save;
    wxButton     *m_clear;

    // state
    bool m_got_results;
    bool m_saved_results;

    std::vector<ResultDisplayItem> m_rdi_array;   // Results to display
    std::vector<MatchType>         m_match_types; // The match types that were tested
    float                          m_lr_thresh;   // The LR threshold of the results
    std::map<int, std::string>     m_c_ids;       // IDs of matching profiles in the C-set
    std::map<int, std::string>     m_r_ids;       // IDs of matching profiles in the R-set

    DECLARE_EVENT_TABLE()

private:
    void setup(wxString const &title);
    bool formatSci();

};

#endif /* RESULTSPANEL_H_ */

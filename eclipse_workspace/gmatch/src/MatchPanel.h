/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MatchPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef MATCHPANEL_H_
#define MATCHPANEL_H_

#include "wx/wx.h"

class MatchPanel : public wxPanel
{
public:
    // Constructor
	MatchPanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
    , m_popTxt(0)
    , m_lrthreshTxt(0)
    , m_cdeltaTxt(0)
    , m_rdeltaTxt(0)
	, m_string0Txt(0)
	, m_string1Txt(0)
	, m_string2Txt(0)
    , m_relList(0)
    {
        setup(title);
    }

    void OnBrowse(wxCommandEvent& event);
    double crimeDelta();
    double refDelta();
    double LRThreshold();
    bool relIsSelected(int n);
    std::string relString(int n);


    wxTextCtrl          *m_popTxt;
    wxTextCtrl          *m_lrthreshTxt;
    wxTextCtrl          *m_cdeltaTxt;
    wxTextCtrl          *m_rdeltaTxt;
    wxTextCtrl          *m_string0Txt;
    wxTextCtrl          *m_string1Txt;
    wxTextCtrl          *m_string2Txt;
    wxCheckListBox      *m_relList;

    DECLARE_EVENT_TABLE();

private:
    void setup(wxString const &title);
};

#endif /* MATCHPANEL_H_ */

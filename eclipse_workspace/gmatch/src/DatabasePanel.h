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

#ifndef DATABASEPANEL_H_
#define DATABASEPANEL_H_

#include "wx/wx.h"

class DatabasePanel : public wxPanel
{
public:
    // Constructor
	DatabasePanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
    , m_hostTxt(0)
    , m_userTxt(0)
    , m_passwordTxt(0)
    {
        setup(title);
    }

    void OnChange(wxCommandEvent& event);

	std::string hostname();
	std::string username();
	std::string password();

    wxTextCtrl          *m_hostTxt;
    wxTextCtrl          *m_userTxt;
    wxTextCtrl          *m_passwordTxt;
    wxStaticText        *m_statusTxt;

    DECLARE_EVENT_TABLE();

private:
    void connect();

    void setup(wxString const &title);
};

class DatabaseDialog : public wxDialog
{
public:
    // Constructor
	DatabaseDialog(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxString& caption = _("Change database connection"),
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxCAPTION | wxRESIZE_BORDER | wxSYSTEM_MENU )
    : wxDialog(parent, winid, caption, pos, size, style)
    , m_hostTxt(0)
    , m_userTxt(0)
    , m_passwordTxt(0)
    {
        setup();
    }

    wxTextCtrl          *m_hostTxt;
    wxTextCtrl          *m_userTxt;
    wxTextCtrl          *m_passwordTxt;

//    DECLARE_EVENT_TABLE();

private:
    void setup();
};

#endif /* DATABASEPANEL_H_ */

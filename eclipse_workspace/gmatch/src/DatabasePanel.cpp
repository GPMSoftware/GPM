/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ResultsPanel.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "DatabasePanel.h"
#include "gmatch.h"
#include "MyFrame.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("ResultsPanel");
#include "fand/messages.h"

#include "wx/statline.h"

using namespace std;

BEGIN_EVENT_TABLE(DatabasePanel, wxPanel)
    EVT_BUTTON(ID_CHANGE, DatabasePanel::OnChange)
END_EVENT_TABLE();

extern boost::shared_ptr<Database> db;

void
DatabasePanel::setup(wxString const &title)
{
    // Widget Hierarchy
    wxPanel *panel = this;

    wxStaticText* labelA = new wxStaticText(panel, wxID_ANY, _("Server"));
    wxString s1; s1 << db->getHost();
    m_hostTxt = new wxTextCtrl(panel, wxID_ANY, s1);
    m_hostTxt->SetEditable(false);
    m_hostTxt->SetBackgroundColour(readOnlyGrey);

    wxButton *changeBtn = new wxButton(panel, ID_CHANGE, _("Change..."));

    wxStaticText* labelB = new wxStaticText(panel, wxID_ANY, _("User name"));
    wxString s2; s2 << db->getUser();
    m_userTxt = new wxTextCtrl(panel, wxID_ANY, s2);
    m_userTxt->SetEditable(false);
    m_userTxt->SetBackgroundColour(readOnlyGrey);

    wxStaticText* labelC = new wxStaticText(panel, wxID_ANY, _("Password"));
    wxString s3; s3 << db->getPass();
    m_passwordTxt = new wxTextCtrl(panel, wxID_ANY, s3);
    m_passwordTxt->SetEditable(false);
    m_passwordTxt->SetBackgroundColour(readOnlyGrey);

    wxStaticText* labelD = new wxStaticText(panel, wxID_ANY, _("Status:"));
    wxString s4; s4 << "Not connected";
    m_statusTxt = new wxStaticText(panel, wxID_ANY, s4);

    // Sizer Hierarchy
    wxBoxSizer *dbSizer  = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *dbHSizerA = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *dbHSizerB = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *dbHSizerC = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *dbHSizerD = new wxBoxSizer( wxHORIZONTAL );

    dbSizer->Add(dbHSizerA, 0, wxEXPAND | wxALL, 2);
		dbHSizerA->Add(labelA, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerA->Add(m_hostTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerA->Add(changeBtn, 1, wxALIGN_CENTER | wxLEFT, 10);

	dbSizer->Add(dbHSizerB, 0, wxEXPAND | wxALL, 2);
		dbHSizerB->Add(labelB, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerB->Add(m_userTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerB->AddStretchSpacer(1);

	dbSizer->Add(dbHSizerC, 0, wxEXPAND | wxALL, 2);
		dbHSizerC->Add(labelC, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerC->Add(m_passwordTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerC->AddStretchSpacer(1);

	dbSizer->Add(dbHSizerD, 0, wxEXPAND | wxALL, 2);
		dbHSizerD->Add(labelD, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerD->Add(m_statusTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerD->AddStretchSpacer(1);

	panel->SetSizer(dbSizer);

	// establish connection to database
    connect();
}

string
DatabasePanel::hostname()
{
	std::string tmp;  // we can assign from ToUTF8() but not construct!
	tmp = m_hostTxt->GetValue().ToUTF8();
	return tmp;
}

string
DatabasePanel::username()
{
	std::string tmp;
	tmp = m_userTxt->GetValue().ToUTF8();
	return tmp;
}

string
DatabasePanel::password()
{
	std::string tmp;
	tmp = m_passwordTxt->GetValue().ToUTF8();
	return tmp;
}

void
DatabasePanel::connect()
{

	db->setHost(hostname());
	db->setUser(username());
	db->setPass(password());

    if (!db->connect())
	{
        m_statusTxt->SetLabel("Not connected");

        wxMessageBox( _("Database connection failed:.\nCheck hostname, username and password."),
                      _("Error"),
                      wxOK | wxICON_EXCLAMATION, this);
	}
    else
    {
        wxString s; s << "Connected to " << hostname();
    	m_statusTxt->SetLabel(s);
    }
}

void
DatabasePanel::OnChange(wxCommandEvent& event)
{
    cout << "DatabasePanel::OnChange(): Id = " << event.GetId() << endl;

    DatabaseDialog dialog(this);

    // populate dialog with current values
    dialog.m_hostTxt->SetValue(m_hostTxt->GetValue());
    dialog.m_userTxt->SetValue(m_userTxt->GetValue());
    dialog.m_passwordTxt->SetValue(m_passwordTxt->GetValue());

    if(dialog.ShowModal() == wxID_OK)
    {
    	cout << "DatabasePanel::OnChange(): OK" << endl;

    	// populate tab with new values
        m_hostTxt->SetValue(dialog.m_hostTxt->GetValue());
        m_userTxt->SetValue(dialog.m_userTxt->GetValue());
        m_passwordTxt->SetValue(dialog.m_passwordTxt->GetValue());

        // re-open the database connection
    	(void)db->disconnectAll();
        connect();

        // discard anything we have read from the database
        // refresh selectors, etc
        MyFrame::refresh();
    }
    else
    {
    	cout << "DatabasePanel::OnChange(): CANCEL" << endl;
    }
}

void
DatabaseDialog::setup()
{
    // Widget Hierarchy
    wxDialog *dialog = this;

    wxStaticText* labelA = new wxStaticText(dialog, wxID_ANY, _("Server"));
    m_hostTxt = new wxTextCtrl(dialog, ID_HOST);

    wxStaticText* labelB = new wxStaticText(dialog, wxID_ANY, _("User name"));
    m_userTxt = new wxTextCtrl(dialog, ID_USER);

    wxStaticText* labelC = new wxStaticText(dialog, wxID_ANY, _("Password"));
    m_passwordTxt = new wxTextCtrl(dialog, ID_PASS);

    wxStaticLine* line1 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    wxButton* okbtn = new wxButton ( this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0 );
    wxButton* cancel = new wxButton ( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );

    // Sizer Hierarchy
    wxBoxSizer *dbSizer  = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *dbHSizerA = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *dbHSizerB = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *dbHSizerC = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *controlSizer = new wxBoxSizer( wxHORIZONTAL );

    dbSizer->Add(dbHSizerA, 0, wxEXPAND | wxALL, 2);
		dbHSizerA->Add(labelA, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerA->Add(m_hostTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerA->AddStretchSpacer(1);

	dbSizer->Add(dbHSizerB, 0, wxEXPAND | wxALL, 2);
		dbHSizerB->Add(labelB, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerB->Add(m_userTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerB->AddStretchSpacer(1);

	dbSizer->Add(dbHSizerC, 0, wxEXPAND | wxALL, 2);
		dbHSizerC->Add(labelC, 1, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerC->Add(m_passwordTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		dbHSizerC->AddStretchSpacer(1);

	dbSizer->Add(line1, 0, wxGROW|wxALL, 5);

	dbSizer->Add(controlSizer, 0, wxEXPAND | wxALL/*, 10*/);
		controlSizer->Add(okbtn, 0, wxALIGN_LEFT|wxALL, 5);
		controlSizer->Add(100, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);
		controlSizer->Add(cancel, 0, wxALIGN_RIGHT|wxALL, 5);

    SetSizer(dbSizer);
    GetSizer()->Fit(this);
    GetSizer()->SetSizeHints(this);
}

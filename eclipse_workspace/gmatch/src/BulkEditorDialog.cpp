/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * BulkEditorDialog.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "BulkEditorDialog.h"
#include "MyFrame.h"

#include "gmatch.h"

#include "wx/statline.h"

#include "fand/dbmatch.h"
#include "fand/geotrans.h"
#include "fand/util.h"


#include "fand/MessageStream.h"
INIT_MESSAGES("BulkEditorDialog");
#include "fand/messages.h"

using namespace std;

BEGIN_EVENT_TABLE(BulkEditorDialog, wxDialog)
    EVT_BUTTON(ID_APPLY, BulkEditorDialog::OnApply)
    EVT_CHOICE(ID_CHOICE, BulkEditorDialog::OnChange)
    EVT_TEXT(ID_NEWVAL, BulkEditorDialog::OnChange)
END_EVENT_TABLE();

BulkEditorDialog * BulkEditorDialog::theBulkEditorDialog = 0;

void BulkEditorDialog::profilesFound()
{
	if (theBulkEditorDialog)
	{
		bool enable = ( theBulkEditorDialog->m_sp_panel->m_selected_profiles.size() > 0 );
		BulkEditorDialog::theBulkEditorDialog->enableEdits(enable);
	}
}

void
BulkEditorDialog::setup()
{
//	std::map<std::string, DBField> meta_fields;
//
//	extern boost::shared_ptr<Database> db;
//
//	if (!db->listMetaFields(meta_fields))
//    {
//        error << "BulkEditorDialog::setup: can't read metadata field names" << endl;
//    }

    wxFont font(
    		8,                  // size
    		wxFONTFAMILY_SWISS, // family
    		wxNORMAL,           // style
    		wxNORMAL,           // weight
    		false,              // underlined
    		_("Sans")           // face
    		);

    SetFont(font);

	wxStaticText *label1 = new wxStaticText(this, wxID_ANY, _("Select Profiles:"));

    m_sp_panel = new SelectProfilesPanel(this, _("Select"), wxID_ANY, wxDefaultPosition, wxSize(100,200), wxSIMPLE_BORDER);
    m_sp_panel->setFindNotify(profilesFound);

    wxStaticLine* line1 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );

    wxStaticText *label2 = new wxStaticText(this, wxID_ANY, _("Edit Profiles:"));

    m_action = new wxChoice(this, ID_CHOICE);
    vector<string> vecstr;
    vecstr.push_back("Delete Profiles");
    vecstr.push_back("Move profiles to dataset:");
    setItemContainerStrings(m_action, vecstr);
    m_action->SetSelection(0);

    m_value = new wxTextCtrl(this, ID_NEWVAL);
    m_apply = new wxButton ( this, ID_APPLY, wxT("&Apply"), wxDefaultPosition, wxDefaultSize, 0 );

    wxStaticLine* line2 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    wxButton* okbtn = new wxButton ( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );

    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *editSizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *controlSizer = new wxBoxSizer( wxHORIZONTAL );

    topSizer->Add(label1, 0, wxGROW|wxALL, 5);
    topSizer->Add(m_sp_panel, 1, wxGROW|wxALL, 5);
    topSizer->Add(line1, 0, wxGROW|wxALL, 5);
    topSizer->Add(editSizer, 0, wxGROW|wxALL, 5);
		editSizer->Add(label2, 0, wxCENTRE|wxALL, 5);
		editSizer->Add(m_action, 0, wxGROW|wxALL, 5);
		editSizer->Add(m_value, 1, wxGROW|wxALL, 5);
		editSizer->Add(m_apply, 0, wxGROW|wxALL, 5);

    topSizer->Add(line2, 0, wxGROW|wxALL, 5);
    topSizer->Add(controlSizer, 0, wxEXPAND | wxALL/*, 10*/);
		controlSizer->Add(1, 1, 1);
        controlSizer->Add(okbtn, 0, wxCENTRE|wxALL, 5);
		controlSizer->Add(1, 1, 1);

    SetSizer(topSizer);
    topSizer->SetSizeHints(this);

    enableEdits(false);
}

void
BulkEditorDialog::enableEdits(bool enabled)
{
	bool delete_flag = (m_action->GetSelection() == 0);
	string val = std::string(m_value->GetValue().ToUTF8());
	cleanup(val);

	m_action->Enable(enabled);
	m_value->Enable(enabled && !delete_flag);
	m_apply->Enable(enabled && (delete_flag || !val.empty()) );
}

void
BulkEditorDialog::popup()
{
	theBulkEditorDialog = new BulkEditorDialog(NULL, wxID_ANY, _("Bulk Profile Editor"), wxDefaultPosition, wxSize(200,200));

    int ret;
    if ((ret = theBulkEditorDialog->ShowModal()) == wxID_OK)
    {
    	theBulkEditorDialog = 0;
	}
}

bool
BulkEditorDialog::doDelete()
{
	string sql_query = m_sp_panel->sqlQuery();
	if (!db->deleteProfiles(sql_query))
	{
		error << startl << "BulkEditorDialog::doDelete(): " << "sql_query" << " : DELETE failed " << endl;
		MyFrame::refreshMessageWindow();
		return false;
	}
	return true;
}

bool
BulkEditorDialog::doMove()
{
	string dataset = string(m_value->GetValue().ToUTF8());
	string where_clause = m_sp_panel->sqlWhereClause();

	if (!db->moveProfiles(dataset, where_clause))
	{
		error << startl << "BulkEditorDialog::doMove(): '" << "sql_query" << "' to " << dataset << " : move failed " << endl;
		MyFrame::refreshMessageWindow();
		return false;
	}
	return true;
}

void
BulkEditorDialog::apply()
{
	if (m_action->GetSelection() == 0) // Delete Profiles
	{
		ostringstream oss;
		int n = m_sp_panel->m_selected_profiles.size();
		oss << "Delete " << n << " profiles?\n (Cannot be undone)";

        if (wxYES == wxMessageBox(oss.str(), _("Delete Profiles"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
    		if (doDelete())
    		{
    			// void the selection
    			m_sp_panel->hasChanged();
    			m_sp_panel->refresh();
    		}
        }
	}
	else // Set dataset
	{
		string val = std::string(m_value->GetValue().ToUTF8());
		cleanup(val);

		ostringstream oss;
		int n = m_sp_panel->m_selected_profiles.size();
		oss << "Move " << n << " profiles?\n (Cannot be undone)";

        if (wxYES == wxMessageBox(oss.str(), _("Move to Dataset"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
    		if (doMove())
    		{
    			// void the selection
    			m_sp_panel->hasChanged();
    			m_sp_panel->refresh();
    		}
        }
	}
}

void
BulkEditorDialog::OnApply(wxCommandEvent& event)
{
    cout << "OnApply(): Id = " << event.GetId() << endl;
    apply();
    MyFrame::refreshMessageWindow();
}

void
BulkEditorDialog::OnChange(wxCommandEvent& event)
{
    cout << "OnChange(): Id = " << event.GetId() << endl;
    enableEdits(true);
}

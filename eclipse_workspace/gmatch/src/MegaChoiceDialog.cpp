/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MegaChoiceDialog.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "MegaChoiceDialog.h"
#include "MyFrame.h"

#include "gmatch.h"

#include "wx/statline.h"

#include "fand/dbmatch.h"
#include "fand/geotrans.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("MetadataDialog");
#include "fand/messages.h"

using namespace std;

BEGIN_EVENT_TABLE(MegaChoiceDialog, wxDialog)
    EVT_TEXT(ID_SEARCH, MegaChoiceDialog::OnSearch)
    EVT_BUTTON(ID_PREV, MegaChoiceDialog::OnPrev)
    EVT_BUTTON(ID_NEXT, MegaChoiceDialog::OnNext)
    EVT_BUTTON(ID_ADD, MegaChoiceDialog::OnAdd)
    EVT_BUTTON(ID_REMOVE, MegaChoiceDialog::OnRemove)
END_EVENT_TABLE();

static const int NDISPLAY = 100; // number of options to display at once

void
MegaChoiceDialog::setup()
{
	std::map<std::string, DBField> meta_fields;

	extern boost::shared_ptr<Database> db;

	if (!db->listMetaFields(meta_fields))
    {
        error << "MultiSelector::Selector::setup: can't read metadata field names" << endl;
    }

    m_text = new wxTextCtrl(this, ID_SEARCH, _(""));
    m_list = new wxListBox(this, ID_LIST, wxDefaultPosition, wxDefaultSize, 0);

    wxStaticLine* line2 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    wxButton* okbtn = new wxButton ( this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0 );
    wxButton* cancel = new wxButton ( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );

    wxBoxSizer *controlSizer = new wxBoxSizer( wxHORIZONTAL );

	wxStaticText *label1 = new wxStaticText(this, wxID_ANY, _("Starts with:"));
	ostringstream oss;
	oss << "< Previous " << NDISPLAY;
    m_prevBtn = new wxButton ( this, ID_PREV, oss.str(), wxDefaultPosition, wxDefaultSize, 0 );
    oss.str("");
	oss << "Next " << NDISPLAY << " >";
    m_nextBtn = new wxButton ( this, ID_NEXT, oss.str(), wxDefaultPosition, wxDefaultSize, 0 );

    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *hSizer   = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *v1Sizer  = new wxBoxSizer( wxVERTICAL );

    topSizer->Add(hSizer, 1, wxEXPAND | wxALL, 0);
    hSizer->Add(v1Sizer, 2, wxEXPAND | wxALL, 0);

    v1Sizer->Add(label1, 0, wxEXPAND|wxALL, 1);
    v1Sizer->Add(m_text, 0, wxEXPAND|wxALL, 1);
    v1Sizer->Add(m_prevBtn, 0, wxEXPAND|wxALL, 1);
    v1Sizer->Add(m_list, 1, wxEXPAND|wxALL, 1);
    v1Sizer->Add(m_nextBtn, 0, wxEXPAND|wxALL, 1);

    topSizer->Add(line2, 0, wxGROW|wxALL, 5);
    topSizer->Add(controlSizer, 0, wxEXPAND | wxALL/*, 10*/);
        controlSizer->Add(okbtn, 0, wxALIGN_LEFT|wxALL, 5);
        controlSizer->Add(20, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);
        controlSizer->Add(cancel, 0, wxALIGN_RIGHT|wxALL, 5);

	list("", NDISPLAY);

	if (m_mode == MULTI)
	{
	    wxBoxSizer *v2Sizer  = new wxBoxSizer( wxVERTICAL );
	    wxBoxSizer *v3Sizer  = new wxBoxSizer( wxVERTICAL );

	    hSizer->Add(v2Sizer, 1, wxEXPAND | wxALL, 0);
	    v2Sizer->AddSpacer(80);
	    v2Sizer->AddStretchSpacer(1);
	    m_addBtn = new wxButton ( this, ID_ADD, _(">"), wxDefaultPosition, wxDefaultSize, 0 );
	    v2Sizer->Add(m_addBtn, 0, wxEXPAND|wxALL, 1);
	    m_remBtn = new wxButton ( this, ID_REMOVE, _("<"), wxDefaultPosition, wxDefaultSize, 0 );
	    v2Sizer->Add(m_remBtn, 0, wxEXPAND|wxALL, 1);
	    v2Sizer->AddStretchSpacer(1);
	    v2Sizer->AddSpacer(30);

	    m_list2 = new wxListBox(this, ID_LIST, wxDefaultPosition, wxDefaultSize, 0);

	    hSizer->Add(v3Sizer, 2, wxEXPAND | wxALL, 0);
	    v3Sizer->AddSpacer(80);
	    v3Sizer->Add(m_list2, 1, wxEXPAND|wxALL, 1);
	    v3Sizer->AddSpacer(30);
	}

    SetSizer(topSizer);
    GetSizer()->Fit(this);
    GetSizer()->SetSizeHints(this);
}

// list (up to) n items starting at the one that matches hint
void
MegaChoiceDialog::list(string hint, int n)
{
    m_begin = m_vals.lower_bound(hint);
    list(n);
}


// list (up to) n items starting at m_begin
// make m_end point to the last displayed item
// update the prev/next buttons
void
MegaChoiceDialog::list(int n)
{
    m_end = m_begin;

    set<string>::const_iterator v_end = m_vals.end();

    for (int i=0; i<n && m_end != v_end; ++i)
    {
        ++m_end;
    }
    set<string> ids(m_begin, m_end);

    setItemContainerStrings(m_list, ids);

    updatePrevNextButtons();
}

void
MegaChoiceDialog::OnSearch(wxCommandEvent& event)
{
    cout << "OnSearch(): Id = " << event.GetId() << endl;

    // get auto-complete hint
    string hint(m_text->GetValue().ToUTF8());

    list(hint, NDISPLAY);
}

void
MegaChoiceDialog::updatePrevNextButtons()
{
    m_prevBtn->Enable(m_begin != m_vals.begin());
    m_nextBtn->Enable(m_end   != m_vals.end());
}

void
MegaChoiceDialog::OnPrev(wxCommandEvent& event)
{
    cout << "OnPrev(): Event Id = " << event.GetId() << endl;

    set<string>::const_iterator v_begin = m_vals.begin();

    for (int i=0; i<NDISPLAY && m_begin != v_begin; ++i)
    	m_begin--;

    list(NDISPLAY);
}

void
MegaChoiceDialog::OnNext(wxCommandEvent& event)
{
    cout << "OnNext(): Event Id = " << event.GetId() << endl;

    set<string>::const_iterator v_end = m_vals.end();

    for (int i=0; i<NDISPLAY &&  m_begin != v_end; ++i)
    	m_begin++;

    list(NDISPLAY);
}

void
MegaChoiceDialog::OnAdd(wxCommandEvent& event)
{
    cout << "OnAdd(): Event Id = " << event.GetId() << endl;

	wxString str = m_list->GetStringSelection();

	// check not null and not already in the list
	if (!str.empty() && m_list2->FindString(str) == wxNOT_FOUND)
    {
		m_list2->Append(str);
    }
}

void
MegaChoiceDialog::OnRemove(wxCommandEvent& event)
{
    cout << "OnRemove(): Event Id = " << event.GetId() << endl;

    int sel;
	if ((sel = m_list2->GetSelection()) != wxNOT_FOUND)
    {
		m_list2->Delete(sel);
    }
}

wxArrayString
MegaChoiceDialog::GetStringSelections()
{
	Assert(m_mode == MULTI);

	wxArrayString ret;
	size_t s = m_list2->GetCount();
	for (size_t i = 0; i<s; ++i)
	{
		ret.push_back( m_list2->GetString(i) );
	}

	return ret;
}

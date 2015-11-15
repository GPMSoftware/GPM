/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MetadataDialog.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "MetadataDialog.h"
#include "MyFrame.h"

#include "gmatch.h"

#include "wx/statline.h"

#include "fand/dbmatch.h"
#include "fand/geotrans.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("MetadataDialog");
#include "fand/messages.h"

using namespace std;

BEGIN_EVENT_TABLE(MetadataDialog, wxDialog)
    EVT_BUTTON(wxID_APPLY, MetadataDialog::OnApply)
END_EVENT_TABLE();

MetadataDialog::MetaControl::MetaControl(MetadataDialog *parent, const std::string &n, const std::string &f)
: wxPanel(parent, wxID_ANY)
, text(0)
, name(n)
, field(f)
{
	setup();
}

void
MetadataDialog::MetaControl::setup()
{
    wxStaticText *label  = new wxStaticText(this, wxID_STATIC, wxString(name.c_str(), *wxConvCurrent), wxDefaultPosition, wxDefaultSize, 0 );
    text                 = new wxTextCtrl(this, wxID_ANY, _(""), wxDefaultPosition, wxDefaultSize, 0 );
    wxBoxSizer   *hSizer = new wxBoxSizer( wxHORIZONTAL );

    hSizer->Add(label, 1, wxEXPAND|wxALL, 5);
    hSizer->Add(text,  1, wxEXPAND|wxALL, 5);

    SetSizer(hSizer);
//    GetSizer()->Fit(this);
//    GetSizer()->SetSizeHints(this);

}

void
MetadataDialog::setEditable(bool edit)
{
	std::map<std::string, MetaControl*>::iterator it;
	for (it = m_controls.begin(); it != m_controls.end(); ++it)
	{
		// This is a frig to make the lat/long fields not editable.
		// Obviously there should be some provision for this in the metadata definitions file
		if (it->second->name == "SAMPLE_LAT" || it->second->name == "SAMPLE_LONG")
		{
			it->second->text->Enable(false);
		}
		else
		{
			it->second->text->Enable(edit);
		}
	}
}

void
MetadataDialog::setup()
{
	std::map<std::string, DBField> meta_fields;

	extern boost::shared_ptr<Database> db;

	if (!db->listMetaFields(meta_fields))
    {
        error << "MultiSelector::Selector::setup: can't read metadata field names" << endl;
    }

    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

	// loop over all metadata fields, creating a control for each
    std::map<std::string, DBField>::const_iterator it;
	for (it = meta_fields.begin(); it != meta_fields.end(); ++it)
	{
		MetaControl *meta = new MetaControl(this, it->first, it->second.field_name);
	    topSizer->Add(meta, 0, wxEXPAND | wxALL/*, 10*/);
		m_controls[it->second.field_name] = meta;
	}

    wxStaticLine* line2 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    wxButton* okbtn = new wxButton ( this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0 );
    wxButton* apply = new wxButton ( this, wxID_APPLY, wxT("&Apply"), wxDefaultPosition, wxDefaultSize, 0 );
    wxButton* cancel = new wxButton ( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );

    topSizer->Add(line2, 0, wxGROW|wxALL, 5);

    wxBoxSizer *controlSizer = new wxBoxSizer( wxHORIZONTAL );

    topSizer->Add(controlSizer, 0, wxEXPAND | wxALL/*, 10*/);
        controlSizer->Add(okbtn, 0, wxALIGN_LEFT|wxALL, 5);
        controlSizer->Add(20, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);
        controlSizer->Add(apply, 0, wxALIGN_RIGHT|wxALL, 5);
        controlSizer->Add(20, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);
        controlSizer->Add(cancel, 0, wxALIGN_RIGHT|wxALL, 5);

    SetSizer(topSizer);
    GetSizer()->Fit(this);
    GetSizer()->SetSizeHints(this);
}

void
MetadataDialog::populate(MetaData const &metad)
{
	MetaData::const_iterator it;
	for (it = metad.begin(); it != metad.end(); ++it)
	{
		std::map<std::string, MetaControl*>::iterator p;
		p = m_controls.find(it->first);

		if (p == m_controls.end())
		{
			error << startl << "unknown metadata field in MetadataDialog::populate()" << endl;
		}
		else
		{
			p->second->text->SetValue(wxString(it->second.c_str(), *wxConvCurrent));
		}
	}
}

void
MetadataDialog::apply()
{
	// make the lat/lon fields depend on the MGRS field
	std::map<std::string, MetaControl*>::const_iterator p;

	if ((p = m_controls.find("mgrs")) == m_controls.end())
	{
		return;
	}

	string mgrs(p->second->text->GetValue().ToUTF8());

	if (mgrs.empty())
	{
		return;
	}

	LatLonStrings latlon = latLonStrings(mgrs);

	if ((p = m_controls.find("s_lat")) != m_controls.end())
	{
		p->second->text->SetValue(wxString(latlon.lat.c_str(), *wxConvCurrent));
	}

	if ((p = m_controls.find("s_long")) != m_controls.end())
	{
		p->second->text->SetValue(wxString(latlon.lon.c_str(), *wxConvCurrent));
	}
}

void
MetadataDialog::readback(MetaData &metad)
{
	MetaData::iterator it;
	for (it = metad.begin(); it != metad.end(); ++it)
	{
		std::map<std::string, MetaControl*>::const_iterator p;
		p = m_controls.find(it->first);

		if (p == m_controls.end())
		{
			error << startl << "unknown metadata field in MetadataDialog::readback()" << endl;
		}
		else
		{
			it->second = p->second->text->GetValue().ToUTF8();
		}
	}
}

void
MetadataDialog::OnApply(wxCommandEvent& event)
{
    cout << "OnApply(): Id = " << event.GetId() << endl;

    apply();

    MyFrame::refreshMessageWindow();
}

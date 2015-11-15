/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * SelectProfilesPanel.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "SelectProfilesPanel.h"
#include "MegaChoiceDialog.h"
#include "gmatch.h"
#include "MyFrame.h"

#include "fand/dbmatch.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("SelectProfilesPanel");
#include "fand/messages.h"

#include "wx/statline.h"

#include <UnitTest++/UnitTest++.h>

using namespace std;

const unsigned int COMBO_MAX = 20; 			// maximum number of items to display in a combo box
const unsigned int SINGLE_CHOICE_MAX = 100; // maximum number of items to display in a wxSingleChoiceDialog

extern boost::shared_ptr<Database> db;

std::map<std::string, std::string> SelectProfilesPanel::m_searches; // Named searches
std::set<SelectProfilesPanel*> SelectProfilesPanel::m_instances; // pointers to all instances
int MultiSelector::Selector::next_id = 0; // used to generate unique IDs for selectors

BEGIN_EVENT_TABLE(MultiSelector::Selector, wxPanel)
    EVT_BUTTON(ID_NEWB, MultiSelector::Selector::OnNew)
    EVT_BUTTON(ID_DELB, MultiSelector::Selector::OnDel)
    EVT_BUTTON(ID_LOGIC, MultiSelector::Selector::OnLogic)
	EVT_COMBOBOX(ID_FCOMBO, MultiSelector::Selector::OnChangeField)
	EVT_COMBOBOX(ID_VCOMBO, MultiSelector::Selector::OnChangeValue)
	EVT_TEXT(ID_VTEXT, MultiSelector::Selector::OnChangeValue)
    EVT_BUTTON(ID_BROWSE, MultiSelector::Selector::OnBrowse)
END_EVENT_TABLE();

BEGIN_EVENT_TABLE(SelectProfilesPanel, wxPanel)
    EVT_BUTTON(wxID_FIND, SelectProfilesPanel::OnFind)
	EVT_TEXT(ID_KM, SelectProfilesPanel::OnChangeValue)
	EVT_TEXT(ID_MGRS, SelectProfilesPanel::OnChangeValue)
	EVT_COMBOBOX(ID_SCOMBO, SelectProfilesPanel::OnChangeSearch)
    EVT_BUTTON(ID_SAVE_SEARCH, SelectProfilesPanel::OnSave)
    EVT_BUTTON(ID_DEL_SEARCH, SelectProfilesPanel::OnDelete)
END_EVENT_TABLE();

// This is a 'ctype facet' that defines which characters are treated
// as whitespace. We want to use ',' as a separator, and NOT ' ' so that
// values can contain ' '.
struct field_reader: std::ctype<char>
{
    field_reader(): std::ctype<char>(get_table()) {}

    static std::ctype_base::mask const*
    get_table()
    {
        static std::vector<std::ctype_base::mask>
            rc(table_size, std::ctype_base::mask());

        rc['\n'] = std::ctype_base::space;
        rc[','] = std::ctype_base::space;
        return &rc[0];
    }
};

SelectProfilesPanel::SelectProfilesPanel(
        wxWindow *parent,
        wxString const &title,
        wxWindowID winid,
        const wxPoint& pos,
        const wxSize& size,
        long style,
        const wxString& name)
: wxPanel(parent, winid, pos, size, style, name)
, m_saved(0)
, m_saveBtn(0)
, m_delBtn(0)
, m_selectp(0)
, m_findBtn(0)
, m_pcountLbl(0)
, m_km(0)
, m_mgrs(0)
, findNotify(0)
{
    setup(title);

    // add self to set of instances
    m_instances.insert(this);
}


SelectProfilesPanel::~SelectProfilesPanel()
{
	// remove self from list of instances
	std::set<SelectProfilesPanel*>::iterator it = m_instances.find(this);
	Assert2(it != m_instances.end(), "SelectProfilesPanel::~SelectProfilesPanel: oops");
	m_instances.erase(it);
}

MultiSelector::Selector::Selector(MultiSelector *parent)
: wxPanel(parent->m_scrolled, wxID_ANY)
, multisel(parent)
, sizer(0)
, field(0)
, value(0)
, newBtn(0)
, delBtn(0)
, logicBtn(0)
, logic(true) // OR logic
, mode(COMBO_BOX)
, id(next_id++)
{
	setup();
}

void
MultiSelector::Selector::setup()
{

	sizer  = new wxBoxSizer(wxHORIZONTAL);
	field  = new wxComboBox(this, ID_FCOMBO, _(""), wxDefaultPosition, wxSize(200, 30), multisel->fstrings(), wxCB_DROPDOWN | wxCB_READONLY);

	// start with comboBox visible
	value  = new wxComboBox(this, ID_VCOMBO, _(""), wxDefaultPosition, wxDefaultSize, multisel->vstrings(), wxCB_DROPDOWN | wxCB_READONLY);

	// ... and the text/browse controls hidden
	text   = new wxTextCtrl(this, ID_VTEXT, _(""));
	text->Hide();
	browse = new wxButton(this, ID_BROWSE,   _("..."), wxDefaultPosition, wxSize(30, 30));
	browse->Hide();
	browse->SetToolTip(wxT("Browse..."));

	newBtn = new wxButton(this, ID_NEWB,    _("+"), wxDefaultPosition, wxSize(30, 30));
	newBtn->SetToolTip(wxT("new Selector"));

	delBtn = new wxButton(this, ID_DELB,    _("-"), wxDefaultPosition, wxSize(30, 30));
	delBtn->SetToolTip(wxT("delete Selector"));

	logicBtn = new wxButton(this, ID_LOGIC, _("|"), wxDefaultPosition, wxSize(30, 30));
	logicBtn->SetBackgroundColour(*wxLIGHT_GREY);
	logicBtn->SetToolTip(wxT("and / or"));

	sizer->Add(logicBtn, 0, wxEXPAND | wxALL, 1);
	sizer->Add(field,  0, wxEXPAND | wxALL, 1);
	sizer->Add(value,  1, wxEXPAND | wxALL, 1);
	sizer->Add(text,   1, wxEXPAND | wxALL, 1);
	sizer->Add(browse, 0, wxEXPAND | wxALL, 1);
	sizer->Add(newBtn, 0, wxEXPAND | wxALL, 1);
	sizer->Add(delBtn, 0, wxEXPAND | wxALL, 1);

	this->SetSizer(sizer);

	// after everything else is constructed
	field->SetSelection(0); // default to Profile Key
}

wxArrayString
MultiSelector::fstrings() // array of field names
{
	wxArrayString ret;

	for (size_t i=0; i<m_sel_fields.size(); ++i)
	{
		ret.push_back( wxString(m_sel_fields[i].name.c_str(), *wxConvCurrent) );
	}

	return ret;
}

wxArrayString
MultiSelector::vstrings() // array of default values
{
	wxArrayString ret;
//	ret.push_back(_("One"));
//	ret.push_back(_("Two"));
//	ret.push_back(_("Three"));
	return ret;
}

void
MultiSelector::setup(wxString const &title)
{
	std::map<std::string, DBField> meta_fields;

	if (!db->listMetaFields(meta_fields))
    {
        error << "MultiSelector::Selector::setup: can't read metadata field names" << endl;
    }

	m_sel_fields.push_back(SelectorField("Profile Key", "profile_key"));
	m_sel_fields.push_back(SelectorField("Dataset",    "dataset"));

	std::map<std::string, DBField>::const_iterator it;
	for (it = meta_fields.begin(); it != meta_fields.end(); ++it)
	{
		m_sel_fields.push_back(SelectorField(it->first, it->second.field_name, it->second.field_type));
	}

	m_selsizer = new wxBoxSizer(wxVERTICAL);
	m_scrolled = new wxScrolledWindow(this, wxID_ANY, wxPoint(0,0), wxSize(465,150), wxVSCROLL);

	int pixelsPerUnitX = 10;
	int pixelsPerUnitY = 10;
	int unitsX         = 10;
	int unitsY         = 1000;

//	m_scrolled->SetBackgroundColour( "GREEN" );
	m_scrolled->SetScrollbars(pixelsPerUnitX, pixelsPerUnitY, unitsX, unitsY);

	addSelectorPos(0);

	m_scrolled->SetSizer(m_selsizer);
}

int
MultiSelector::findSelector(int id)
{
	for (size_t i=0; i<m_selectors.size(); ++i)
	{
		if (m_selectors[i]->id == id)
		{
			return i;
		}
	}
	return -1;
}

void
MultiSelector::addSelectorPos(int pos)
{
	Selector *sel = new Selector(this);

	if (pos == 0)
	{
		sel->logicBtn->Enable(false);
		sel->delBtn->Enable(false);
	}

	m_selectors.insert(m_selectors.begin()+pos, sel);

	m_selsizer->Insert(pos, sel, 0, wxEXPAND | wxALL, 2);

	m_selsizer->Layout();
	m_scrolled->FitInside();
}

void
MultiSelector::addSelectorId(int id)
{
	int pos = findSelector(id);
	Assert2(pos>=0, "oops");
	addSelectorPos(pos+1); // insert after
}

void
MultiSelector::delSelectorPos(int pos)
{
	Selector *sel = m_selectors[pos];
	m_selectors.erase(m_selectors.begin()+pos);
	sel->Destroy();
	m_selsizer->Layout();
	m_scrolled->FitInside();
}

void
MultiSelector::delSelectorId(int id)
{
	int pos = findSelector(id);
	Assert2(pos>=0, "oops");
	if (m_selectors.size() > 1) // don't delete the only item
	{
		delSelectorPos(pos);
	}
}

// remove all selectors including the first
void
MultiSelector::clear()
{
	int n = m_selectors.size();

	for (int i=n-1; i>=0; --i)
	{
		delSelectorPos(i);
	}
}

const DBField*
MultiSelector::findDBFieldFromName(const string& name)
{
	DBField *ret = 0;

	for (size_t i=0; i<m_sel_fields.size(); ++i)
	{
		if (m_sel_fields[i].name == name)
		{
			ret = &m_sel_fields[i].db_field;
			break;
		}
	}

	return ret;
}

// return the query from the MultiSelector
std::string
MultiSelector::sqlQuery()
{
	string sql_query = "";

	for (size_t i=0; i<m_selectors.size(); ++i)
	{
		string name = m_selectors[i]->getName();
		const DBField *dbf = findDBFieldFromName(name);
		Assert2(dbf, "MultiSelector::sqlQuery(): field not found");

		string value_str = m_selectors[i]->getValue();

		if (dbf->field_name.empty() || value_str.empty()) continue;

		if (! sql_query.empty())
		{
			sql_query = "( " + sql_query + (m_selectors[i]->logic ? " ) OR " : " ) AND ");
		}

		if (dbf->field_type == "DOUBLE") // do a numeric comparison
		{
			// if value field does not start with a comparison op then assume =
			if ( (value_str.find("=")  != 0) &&
				 (value_str.find("!=") != 0) &&
				 (value_str.find("<")  != 0) &&
				 (value_str.find("<=") != 0) &&
				 (value_str.find(">")  != 0) &&
				 (value_str.find(">=") != 0) &&
				 (value_str.find("<>") != 0) )
			{
				value_str = "= " + value_str;
			}

			sql_query += dbf->field_name + " " + value_str; // no quotes needed for numeric value
		}
		else // is a string type
		{
			// does value contain SQL metacharacters?
			if (value_str.find_first_of("_%") != value_str.npos)
			{
				// use LIKE
				sql_query += dbf->field_name + " LIKE '" + value_str + "'";
			}
			else
			{
				// use a string '=' search
				sql_query += dbf->field_name + " = '" + value_str + "'";
			}
		}
	}

	return sql_query;
}

void
MultiSelector::update()
{
	for (size_t i=0; i<m_selectors.size(); ++i)
	{
		m_selectors[i]->update();
	}
}

void
MultiSelector::Selector::OnNew(wxCommandEvent& event)
{
    cout << "OnNew(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;
    multisel->addSelectorId(id);
    multisel->hasChanged();
    MyFrame::refreshMessageWindow();
}

void
MultiSelector::Selector::OnDel(wxCommandEvent& event)
{
    cout << "OnDel(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;
    multisel->delSelectorId(id);
    multisel->hasChanged();
    MyFrame::refreshMessageWindow();
}

void
MultiSelector::Selector::OnLogic(wxCommandEvent& event)
{
    cout << "OnLogic(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;
    logic = ! logic;
	setLogicBtn(logic);

    multisel->hasChanged();
    MyFrame::refreshMessageWindow();
}

static std::string geoBoxQuery(std::string mylat, std::string mylon, std::string dist)
{
	return "( geoBox(s_lat, s_long, " + mylat + ", " + mylon + ", " + dist + ") = TRUE ) ";
}

static std::string geoCircleQuery(std::string mylat, std::string mylon, std::string dist)
{
	return "( geoCircle(s_lat, s_long, " + mylat + ", " + mylon + ", " + dist + ") = TRUE ) ";
}

static std::string distQuery(std::string mylat, std::string mylon, std::string dist)
{
	if (mylat.empty() || mylon.empty() || dist.empty())
	{
		return "";
	}

//	return geoBoxQuery(mylat, mylon, dist); // surprisingly, this is slower than geoCircleQuery
	return geoCircleQuery(mylat, mylon, dist);
}


LatLonStrings
SelectProfilesPanel::latLonStrings()
{
	LatLonStrings ret;
	string mgrs_string(m_mgrs->GetValue().ToUTF8());

	if ( !mgrs_string.empty() )
	{
		ret = ::latLonStrings(mgrs_string);
	}

	return ret;
}

// return the distance query
std::string
SelectProfilesPanel::distQuery()
{
	LatLonStrings latlon = latLonStrings();
	return ::distQuery(latlon.lat, latlon.lon, std::string(m_km->GetValue().ToUTF8()));
}

static std::string
sqlWhereClause(std::string dist_query, std::string sel_query)
{
	// construct where clause like:
	// dist_query AND ( meta_query OR  ... meta_query )

	if ( dist_query.empty() && sel_query.empty() )
	{
		return "";
	}

	string where_clause = dist_query;

	if ( ! dist_query.empty() && ! sel_query.empty() )
	{
		sel_query = "AND (" + sel_query + ")";
	}

	where_clause += sel_query;

	return where_clause;
}

static std::string
sqlFullQuery(std::string dist_query, std::string sel_query)
{
	// construct query like:
	// SELECT * FROM Profiles WHERE where_clause

	string where_clause = sqlWhereClause(dist_query, sel_query);

	if ( where_clause.empty() )
	{
		return "";
	}

	string sql_query = "SELECT * FROM Profiles WHERE " + where_clause;
	return sql_query;
}

// return the query for the entire panel
std::string
SelectProfilesPanel::sqlQuery()
{
	return sqlFullQuery(distQuery(), m_selectp->sqlQuery());
}

// return the query for the entire panel
std::string
SelectProfilesPanel::sqlWhereClause()
{
	return ::sqlWhereClause(distQuery(), m_selectp->sqlQuery());
}

MultiSelector *
SelectProfilesPanel::createMultiSelector(SelectProfilesPanel *parent, wxString const &title)
{
	MultiSelector* panel = new MultiSelector(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

void
SelectProfilesPanel::OnFind(wxCommandEvent& event)
{
    cout << "OnFind(): Id = " << event.GetId() << endl;

    string sql_query = sqlQuery();

    if (sql_query.empty())
    {
        error << "SelectProfilesPanel::OnFind(): SQL query is empty" << endl;
    }
    else
    {
		// don't read - just check the size and return a ProfileRangeDB
		if (!db->isConnected())
		{
			cout << "SelectProfilesPanel::OnFind(): No database connection" << endl;
		}
		else
		{
			ProfileRangeDB prdb(db, sql_query);
			m_selected_profiles = ProfileRange(prdb);

			// report how many found
			size_t len = 32;
			wxChar buf[len];
			int num_profiles_found = m_selected_profiles.size();
			wxSnprintf(buf, len, _("%d Profiles Selected"), num_profiles_found);
			m_pcountLbl->SetLabel(buf);

			if (num_profiles_found == 0)
			{
				error << "SelectProfilesPanel::OnFind(): 0 Profiles found" << endl;
			}
			else
			{
				ok << "SelectProfilesPanel::OnFind(): " << num_profiles_found << " Profiles found" << endl;
			}

			// grey out the button until a control is changed
			m_findBtn->Enable(false);
		}
    }

    if (findNotify) (*findNotify)();
}

void
SelectProfilesPanel::OnDelete(wxCommandEvent& event)
{
    cout << "OnDelete(): Event Id = " << event.GetId() << endl;

	string sel(m_saved->GetValue().ToUTF8());

	if (sel.empty())
	{
		// pop up error dialog
        wxMessageBox( _("Select a search to delete and try again."),
                      _("Delete search"),
                      wxOK | wxICON_EXCLAMATION, this);
		return;
	}

    std::map<std::string, std::string>::iterator it;
    if ( (it = m_searches.find(sel)) == m_searches.end() )
    {
		// pop up error dialog
		wxMessageBox( "Error: no search named " + sel,
					  _("Delete search"),
					  wxOK | wxICON_EXCLAMATION, this);
		return;
    }
    else
    {
    	// delete search
        if (wxYES == wxMessageBox("Delete search " + sel + "?", _("Delete search"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
			cout << "Delete search: " << sel << endl;
			m_searches.erase(it);
			int d = m_saved->FindString(wxString(sel.c_str(), *wxConvCurrent));
			m_saved->Delete(d);
			m_saved->SetValue("");
        }
    }
    saveSearchesToFile();
    updateSearches();
    MyFrame::refreshMessageWindow();
}

void
SelectProfilesPanel::OnSave(wxCommandEvent& event)
{
    cout << "OnSave(): Event Id = " << event.GetId() << endl;

	string sel(m_saved->GetValue().ToUTF8());

	if (sel.empty())
	{
		// pop up error dialog
        wxMessageBox( _("Enter a search name and try again."),
                      _("Save search"),
                      wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	ostringstream oss;
	oss << *this;

    std::map<std::string, std::string>::const_iterator it;
    if ( (it = m_searches.find(sel)) == m_searches.end() )
    {
    	// new saved search
        if (wxYES == wxMessageBox("Save search as " + sel + "?", _("Save search"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
			cout << "Save search: new search " << sel << endl;

			m_searches[sel] = oss.str();
			m_saved->Append(wxString(sel.c_str(), *wxConvCurrent));
        }
    }
    else
    {
    	// pop up confirm dialog to overwrite
        if (wxYES == wxMessageBox(sel + " exists. Overwrite?", _("Save search"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
    		cout << "Save search: overwriting " << sel << endl;
    		m_searches[sel] = oss.str();
        }
    }
    saveSearchesToFile();
    updateSearches();
    MyFrame::refreshMessageWindow();
}

void
SelectProfilesPanel::refresh()
{
	// update all selection lists
	m_selectp->update();
}

// update the choices in the value combo
void
MultiSelector::Selector::updateValues(DBField const *dbf)
{
	vals.clear();

    if (dbf->field_type == "DOUBLE")
    {
    	// use the text widget alone
    	mode = TEXT;
    	value->Hide();
		text->Show();
		text->Clear();
		browse->Hide();
    }
    else
    {

		// NB getting all the values can take a few seconds if we are querying for
		// profile IDs and there are say 1 million profiles.
		// We optimize by searching only for the first SINGLE_CHOICE_MAX + 1 elements.
		// In the case of a MEGA_CHOICE we must do a full search when the browse button is pressed.
		if (!db->listValues(dbf->field_name, vals, SINGLE_CHOICE_MAX + 1))
		{
			error << "SelectProfilesPanel::refreshDatasetList(): " << dbf->field_name << " could not be read: " << endl;
		}

		// remove NULL if present
		set<string>::const_iterator it = vals.find("");
		if (it != vals.end())
		{
			vals.erase(it);
		}

		if (vals.size() > COMBO_MAX)
		{
			if (vals.size() > SINGLE_CHOICE_MAX)
			{
				// use the text/browse buttons with MegaChoice dialog
				mode = MEGA_CHOICE;
			}
			else
			{
				// use the text/browse buttons with wxSingleChoiceDialog
				mode = SINGLE_CHOICE;
			}

			value->Hide();
			text->Show();
			text->Clear();
			browse->Show();
		}
		else
		{
			// use Combo box
			mode = COMBO_BOX;
			value->Show();
			text->Hide();
			browse->Hide();

			if (vals.empty())
			{
				value->Clear();
			}
			else
			{
				setItemContainerStrings(value, vals);
				value->SetSelection(0);
			}
		}
    }

    multisel->m_selsizer->Layout();
	multisel->m_scrolled->FitInside();
}

std::string
MultiSelector::Selector::getLogic() const
{
	return logic ? "OR" : "AND";
}

std::string
MultiSelector::Selector::getName() const
{
	string name(field->GetValue().ToUTF8());
	return name;
}

void
MultiSelector::Selector::setLogic(string logic_str)
{
	logic = logic_str == "OR" ? true : false;
	setLogicBtn(logic);
}

void
MultiSelector::Selector::setLogicBtn(bool logic)
{
    logicBtn->SetLabel(logic ? _("|") : _("&&"));
}

void
MultiSelector::Selector::setName(string const &name)
{
	field->SetValue(name);
}

std::string
MultiSelector::Selector::getValue() const
{
	string value_str;
	switch (mode)
	{
	case Selector::COMBO_BOX:
		value_str = value->GetValue().ToUTF8();
		break;
	case Selector::SINGLE_CHOICE:
	case Selector::TEXT:
	case Selector::MEGA_CHOICE:
		value_str = text->GetValue().ToUTF8();
		break;
	default:
		Assert2(0, "MultiSelector::Selector::getValue(): oops");
	}

	cleanup(value_str); // remove leading/trailing space

	return value_str;
}

std::string
MultiSelector::Selector::getStreamedValue() const
{
	// for streaming we represent the empty string as ' '
	string ret = getValue();
	if (ret.empty()) ret = " ";
	return ret;
}

void
MultiSelector::Selector::setValue(string const &val_str)
{
	string value_str;
	switch (mode)
	{
	case Selector::COMBO_BOX:
		value->SetValue(val_str);
		break;
	case Selector::SINGLE_CHOICE:
	case Selector::TEXT:
	case Selector::MEGA_CHOICE:
		text->SetValue(val_str);
		break;
	default:
		Assert2(0, "MultiSelector::Selector::setValue(): oops");
	}
}

void
MultiSelector::Selector::setStreamedValue(string const &val_str)
{
	// for streaming we represent the empty string as ' '
	setValue( val_str == " " ? "" : val_str);
}

const DBField *
MultiSelector::Selector::getDbf()
{
	string name(field->GetValue().ToUTF8());
	const DBField *dbf = multisel->findDBFieldFromName(name);
	if (!dbf)
	{
		cout << "MultiSelector::Selector::dbf(): field name = " << name << endl;
		Assert2(false, "MultiSelector::Selector::dbf(): field not found");
	}
	return dbf;
}

void
MultiSelector::Selector::update()
{
	updateValues(getDbf());
}

void
MultiSelector::Selector::OnChangeField(wxCommandEvent& event)
{
    cout << "value(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;

	// update the choices in the value combo
    update();

    multisel->hasChanged();
    MyFrame::refreshMessageWindow();
}

void
MultiSelector::Selector::OnChangeValue(wxCommandEvent& event)
{
    cout << "OnChangeValue(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;

    multisel->hasChanged();
    MyFrame::refreshMessageWindow();
}

void
MultiSelector::Selector::OnBrowse(wxCommandEvent& event)
{
    cout << "OnBrowse(): Event Id = " << event.GetId() << " Selector ID = " << id << endl;

    switch(mode)
    {
	case SINGLE_CHOICE:
		singleChoice();
		break;

    case MEGA_CHOICE:
    	megaChoice();
		break;

    default:
		Assert2(0, "MultiSelector::Selector::OnBrowse: mode error");
    }
}

void
MultiSelector::Selector::megaChoice()
{
	// make sure we have the entire list of values
	const DBField *dbf = getDbf();
    if (!db->listValues(dbf->field_name, vals))
    {
        error << "MultiSelector::Selector::megaChoice(): " << dbf->field_name << " could not be read: " << endl;
    }


	// pop up MegaChoiceDialog (SINGLE mode) to select profile key
	MegaChoiceDialog dialog(NULL, wxID_ANY, _("Select Profile Key"), vals);

	int ret;
	if ((ret = dialog.ShowModal()) == wxID_OK)
	{
		text->SetValue(dialog.GetStringSelection());

		multisel->hasChanged();
		MyFrame::refreshMessageWindow();
	}
}

void
MultiSelector::Selector::singleChoice()
{

    // pop up a single selection dialog
	wxArrayString choices;

	set<string>::const_iterator it;
	for (it = vals.begin(); it != vals.end(); ++it)
	{
		choices.push_back( wxString(it->c_str(), *wxConvCurrent) );
	}

    wxSingleChoiceDialog dialog(this, field->GetValue(), _("Select"), choices);

    dialog.SetSelection(0);

    if (dialog.ShowModal() == wxID_OK)
    {
		text->SetValue(dialog.GetStringSelection());

		multisel->hasChanged();
	    MyFrame::refreshMessageWindow();
    }
}

void
SelectProfilesPanel::OnChangeValue(wxCommandEvent& event)
{
    cout << "OnChangeValue(): Event Id = " << event.GetId() << endl;

    hasChanged();
    MyFrame::refreshMessageWindow();
}

void
SelectProfilesPanel::hasChanged()
{
    // when one of the controls is changed, reset the number of profiles to 0
    // and enable the Find Profiles button
    m_selected_profiles.clear();
    m_pcountLbl->SetLabel(_("0 Profiles Selected"));
    m_findBtn->Enable(true);

    m_saved->SetValue("");

    if (findNotify) (*findNotify)();
}

void addDBFunctions(Database & db)
{
	// The magic number 111 is the distance subtended by one degree at the surface of the earth in km
	// R = 6378km
	// D = R * 2pi / 360 = 111

	string geoBoxDefn = "CREATE FUNCTION geoBox (lat DOUBLE, lon DOUBLE, mylat DOUBLE, mylon DOUBLE, dist DOUBLE)\n"
			            "RETURNS BOOL DETERMINISTIC\n"
			            "BEGIN\n"
						"set @dlon = dist/abs(cos(radians(mylat))*111);\n"
						"set @dlat = dist/111;\n"
			            "RETURN lon BETWEEN (mylon - @dlon) AND (mylon + @dlon) AND\n"
						"       lat BETWEEN (mylat - @dlat) AND (mylat + @dlat);\n"
	                    "END\n";

	string geoDistDefn = "CREATE FUNCTION geoDist (lat DOUBLE, lon DOUBLE, mylat DOUBLE, mylon DOUBLE)\n"
			            "RETURNS DOUBLE DETERMINISTIC\n"
			            "BEGIN\n"
						"RETURN 6378 * 2 * ASIN(SQRT( POWER(SIN(radians(mylat - lat)/ 2), 2) + "
						"COS(radians(mylat)) * COS(radians(lat)) * "
						"POWER(SIN(radians(mylon - lon) / 2), 2) ));\n"
			            "END\n";

	string geoCircleDefn = "CREATE FUNCTION geoCircle (lat DOUBLE, lon DOUBLE, mylat DOUBLE, mylon DOUBLE, dist DOUBLE)\n"
			            "RETURNS BOOL DETERMINISTIC\n"
			            "RETURN (geoDist(lat, lon, mylat, mylon) <= dist)\n";

	(void)db.do_sql("DROP FUNCTION geoBox");
	(void)db.do_sql("DROP FUNCTION geoDist");
	(void)db.do_sql("DROP FUNCTION geoCircle");

	if (!db.do_sql(geoBoxDefn) || !db.do_sql(geoDistDefn) || !db.do_sql(geoCircleDefn))
    {
        error << "SelectProfilesPanel::setup: can't define geoBox/geoDist functions" << endl;
    }
}

void
SelectProfilesPanel::setup(wxString const &title)
{
	addDBFunctions(*db);

    wxPanel *panel = this;

    m_findBtn = new wxButton(panel, wxID_FIND, _("Find Profiles"));
	wxStaticText *label0 = new wxStaticText(panel, wxID_ANY, _("Searches:"));
    m_saved = new wxComboBox(panel, ID_SCOMBO, wxEmptyString, wxDefaultPosition, wxSize(80, wxDefaultCoord));
    m_saveBtn = new wxButton(panel, ID_SAVE_SEARCH, _("Save..."), wxDefaultPosition, wxSize(50, wxDefaultCoord));
    m_delBtn = new wxButton(panel, ID_DEL_SEARCH, _("Delete..."), wxDefaultPosition, wxSize(50, wxDefaultCoord));
	m_pcountLbl = new wxStaticText(panel, wxID_ANY, _("0 Profiles Selected"));
    wxStaticLine* line1 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	wxStaticText *label1 = new wxStaticText(panel, wxID_ANY, _(" Metadata:"));
    m_selectp = createMultiSelector(this, _("title"));
    wxStaticLine* line2 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
//	wxStaticText *label2 = new wxStaticText(panel, wxID_ANY, _(" AND Distance:"));

	wxStaticText *l_dist = new wxStaticText(panel, wxID_ANY, _("&& Distance"));
	m_km                 = new wxTextCtrl(this, ID_KM, _(""));
	wxStaticText *l_km   = new wxStaticText(panel, wxID_ANY, _("(km) from point"));
	m_mgrs               = new wxTextCtrl(this, ID_MGRS, _(""));
	wxStaticText *l_mgrs = new wxStaticText(panel, wxID_ANY, _("(MGRS)"));

	wxBoxSizer *vSizer  = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *hSizer  = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *hSizer2 = new wxBoxSizer( wxHORIZONTAL );

    vSizer->Add(hSizer, 0, wxEXPAND | wxALL, 5);
		hSizer->Add(m_findBtn, 0, wxCENTRE | wxRIGHT, 2);
		hSizer->Add(m_pcountLbl, 0, wxCENTRE | wxRIGHT, 2);
	    hSizer->Add(1, 1,   1, wxEXPAND | wxALL, 2);
		hSizer->Add(label0, 0, wxCENTRE | wxLEFT, 2);
		hSizer->Add(m_saved, 0, wxCENTRE | wxLEFT, 2);
		hSizer->Add(m_saveBtn, 0, wxCENTRE | wxLEFT, 2);
		hSizer->Add(m_delBtn, 0, wxCENTRE | wxLEFT, 2);
    vSizer->Add(line1,     0, wxEXPAND | wxUP, 5);
    vSizer->Add(label1,    0, wxEXPAND | wxUP, 5);
    vSizer->Add(m_selectp, 1, wxEXPAND | wxUP, 5);
    vSizer->Add(line2,     0, wxEXPAND | wxUP, 5);
//    vSizer->Add(label2,    0, wxEXPAND | wxUP, 5);
    vSizer->Add(hSizer2,   0, wxEXPAND | wxALL, 5);
		hSizer2->Add(l_dist, 0, wxCENTRE | wxRIGHT, 5);
		hSizer2->Add(m_km,   1, wxEXPAND | wxRIGHT, 5);
		hSizer2->Add(l_km,   0, wxCENTRE | wxRIGHT, 5);
		hSizer2->Add(m_mgrs, 1, wxEXPAND | wxRIGHT, 5);
		hSizer2->Add(l_mgrs, 0, wxCENTRE | wxRIGHT, 5);

    panel->SetSizer(vSizer);
    vSizer->SetSizeHints(this);

    loadSearchesFromFile();
}

void
SelectProfilesPanel::loadSearchesFromFile()
{
	ifstream ifs(".fand_searches");

    m_searches.clear();
    string line;
	while (getline(ifs, line))
	{
		string search_name, search;
		int n = line.find(',');

		search_name = line.substr(0, n); // up to first ','
		search      = line.substr(n+1);  // beyond first ','

		m_searches[search_name] = search;
	}

	refreshSearches();
}

void
SelectProfilesPanel::refreshSearches()
{
	// refresh the list of searches displayed in the ComboBox
	m_saved->Clear();
    std::map<std::string, std::string>::const_iterator it;
	for (it = m_searches.begin(); it != m_searches.end(); ++it)
	{
		m_saved->Append(wxString(it->first.c_str(), *wxConvCurrent));
	}
}

void
SelectProfilesPanel::saveSearchesToFile()
{
	// Save as one line per search

	ofstream ofs(".fand_searches");

    std::map<std::string, std::string>::const_iterator it;
	for (it = m_searches.begin(); it != m_searches.end(); ++it)
	{
		ofs << it->first << "," << it->second << "\n";
	}
	ofs.close();
}

void
SelectProfilesPanel::updateSearches()
{
	std::set<SelectProfilesPanel*>::iterator it;

	for (it = m_instances.begin(); it != m_instances.end(); ++it)
	{
		(*it)->refreshSearches();
	}
}

std::ostream &
operator<<(std::ostream &os, const MultiSelector &a)
{
	for (int i=a.m_selectors.size()-1; i>=0; --i) // Reverse order
	{
		MultiSelector::Selector *s = a.m_selectors[i];
		os << s->getLogic() << "," << s->getName() << "," << s->getStreamedValue() << ",";
	}
	return os;
}

std::istream &
operator>>(std::istream &is, MultiSelector &a)
{
	a.clear();

    is.imbue(std::locale(std::locale(), new field_reader()));

	string logic, name, value;
	while (is >> logic >> name >> value)
	{
		a.addSelectorPos(0);
		MultiSelector::Selector *s = a.m_selectors[0];

		s->setLogic(logic);
		s->setName(name);
	    s->update();
		s->setStreamedValue(value);

		if (a.m_selectors.size() > 1)
		{
			s = a.m_selectors[1];
			s->logicBtn->Enable(true);
			s->delBtn->Enable(true);
		}
	}

	if (a.m_selectors.empty())
	{
		a.addSelectorPos(0);
	}

    MyFrame::refreshMessageWindow();
	return is;
}

std::string
getStreamed(wxTextCtrl *p)
{
	string ret = std::string(p->GetValue().ToUTF8());
	return ret == "" ? " " : ret;
}

void
setStreamed(wxTextCtrl *p, string const &val)
{
	p->SetValue(val == " " ? "" : val);
}

std::ostream &
operator<<(std::ostream &os, const SelectProfilesPanel &a)
{
	// distance
	os << getStreamed(a.m_km) << "," << getStreamed(a.m_mgrs) << ",";

	// MultiSelector
	if (a.m_selectp)
	{
		os << *(a.m_selectp);
	}
	return os;
}

std::istream &
operator>>(std::istream &is, SelectProfilesPanel &a)
{
    is.imbue(std::locale(std::locale(), new field_reader()));

	// distance
	string km_str, mgrs_str;
	is >> km_str >> mgrs_str;

	setStreamed(a.m_km, km_str);
	setStreamed(a.m_mgrs, mgrs_str);

	// MultiSelector
	if (a.m_selectp)
	{
		is >> *(a.m_selectp);
	}
	return is;
}

void
SelectProfilesPanel::OnChangeSearch(wxCommandEvent& event)
{
    cout << "OnChangeSearch(): Event Id = " << event.GetId() << endl;
    recallSearch();
    MyFrame::refreshMessageWindow();
}

void
SelectProfilesPanel::recallSearch()
{
	string sel(m_saved->GetStringSelection().ToUTF8());

    std::map<std::string, std::string>::const_iterator it;
    if ( (it = m_searches.find(sel)) != m_searches.end() )
    {
		istringstream iss(it->second);
		iss >> *this;
    }

    hasChanged();

    // make sure we are still displaying the selected search name
    m_saved->SetValue(sel);
}

//                A
//   --------------------------
//   |                   \  B |
//   |                    \   |
//   |                 C   \  |
// F |            X           |   D
//   |                        |
//   |                        |
//   |            E           |
//   --------------------------
//
// dist = 111km; d_lat = 1; d_long = 0.5;
// X(60, -3)
// A(62, -3)       : 222 km N
// B(60.72, -1.56) : 113 km NE
// C(60.63, -1.74) :  99 km NE
// D(60,    -0.94) : 120 km E
// E(59.1,  -3)    : 100 km S
// F(60,    -5.16) : 120 km W
//
static char meta_defns[] =
		"SAMPLE_LAT,   s_lat,  DOUBLE\n"
		"SAMPLE_LONG,  s_long, DOUBLE\n"
        "WIS,          wis,    VARCHAR(16)\n"
        "CEXE,         cexe,   VARCHAR(16)\n"
        "LAB,          lab,    VARCHAR(16)\n"
        "TAT,          tat,    VARCHAR(16)\n"
        "DEF,          def,    VARCHAR(16)\n";

struct Dist_test
{
	Dist_test()
    : db(new Database("test", meta_defns))
    , idA(Identifiler,
            "SAMPLE_A-FULL",        // ID
            1,                       // contributors
            0.001,                   // delta
            0,                       // mutation rate
            reference,               // evidence_type
            "BADDIES",               // dataset
            "SAMPLE_A",              // sample_id
            "FULL")                 // profile_id
    , idB(Identifiler, "SAMPLE_B-FULL", 1, 0.001, 0, reference, "BADDIES", "SAMPLE_B", "FULL")
    , idC(Identifiler, "SAMPLE_C-FULL", 1, 0.001, 0, reference, "BADDIES", "SAMPLE_C", "FULL")
    , idD(Identifiler, "SAMPLE_D-FULL", 1, 0.001, 0, reference, "BADDIES", "SAMPLE_D", "FULL")
    , idE(Identifiler, "SAMPLE_E-FULL", 1, 0.001, 0, reference, "BADDIES", "SAMPLE_E", "FULL")
    , idF(Identifiler, "SAMPLE_F-FULL", 1, 0.001, 0, reference, "BADDIES", "SAMPLE_F", "FULL")
    , A(idA), B(idB), C(idC), D(idD), E(idE), F(idF)
    {
        A.m_metadata["s_lat"]  = "62";
        A.m_metadata["s_long"] = "-3";

        B.m_metadata["s_lat"]  = "60.72";
        B.m_metadata["s_long"] = "-1.56";

        C.m_metadata["s_lat"]  = "60.63";
        C.m_metadata["s_long"] = "-1.74";

        D.m_metadata["s_lat"]  = "60";
        D.m_metadata["s_long"] = "-0.94";

        E.m_metadata["s_lat"]  = "59.1";
        E.m_metadata["s_long"] = "-3";

        F.m_metadata["s_lat"]  = "60";
        F.m_metadata["s_long"] = "-5.16";
    }

	boost::shared_ptr<Database> db;
    ProfileData idA, idB, idC, idD, idE, idF;
    DBProfile A, B, C, D, E, F;
};

// check distance functions
TEST_FIXTURE(Dist_test, dist_test)
{
    CHECK_EQUAL(true, db->connect());
    CHECK_EQUAL(true, db->clear());
    CHECK_EQUAL(0, db->size());
    CHECK_EQUAL(true, db->insert(A));
    CHECK_EQUAL(true, db->insert(B));
    CHECK_EQUAL(true, db->insert(C));
    CHECK_EQUAL(true, db->insert(D));
    CHECK_EQUAL(true, db->insert(E));
    CHECK_EQUAL(true, db->insert(F));
    CHECK_EQUAL(6, db->size());

    addDBFunctions(*db);

    //
	// test geoBox
    //
	string sql_query1 = sqlFullQuery(geoBoxQuery("60", "-3", "111"), "");
//	cout << sql_query1 << endl;

	ProfileRangeDB prdb_box(db, sql_query1);

	// should have found B, C, E
	CHECK_EQUAL(3, prdb_box.size());

	prdb_box.readMyData();

	CHECK_EQUAL("SAMPLE_B-FULL", prdb_box[0].data().m_id.c_str());
	CHECK_EQUAL("SAMPLE_C-FULL", prdb_box[1].data().m_id.c_str());
	CHECK_EQUAL("SAMPLE_E-FULL", prdb_box[2].data().m_id.c_str());

	//
	// test geoCircle
	//
	string sql_query2 = sqlFullQuery(geoCircleQuery("60", "-3", "111"), "");
	ProfileRangeDB prdb_dist(db, sql_query2);

	// should have found C, E
	CHECK_EQUAL(2, prdb_dist.size());

	prdb_dist.readMyData();

	CHECK_EQUAL("SAMPLE_C-FULL", prdb_dist[0].data().m_id.c_str());
	CHECK_EQUAL("SAMPLE_E-FULL", prdb_dist[1].data().m_id.c_str());
}

void
setID(ProfileData &id, char c, int i)
{
    ostringstream oss;
    oss << "SAMPLE_ " << c << i;
	id.m_id = oss.str();
	id.m_profile_id = "FULL";
	id.m_sample_id = id.m_id + "-" + id.m_profile_id;
}

#if 0
// results:
//
// control: count profiles = 0.02
// control: read profiles = 4.04
// geoBox: count profiles = 0.91
// geoBox: read profiles = 3.16
// geoCircle: count profiles = 0.55
// geoCircle: read profiles = 2.44
//
// ie geoCircle is actually faster than geoBox
//
TEST_FIXTURE(Dist_test, dist_timings)
{
    CHECK_EQUAL(true, db->connect());
    CHECK_EQUAL(true, db->clear());
    CHECK_EQUAL(0, db->size());

	int N = 10000;
    for (int i=0; i<N; ++i)
    {
		setID(A.m_data, 'A', i);
		CHECK_EQUAL(true, db->insert(A));

		setID(B.m_data, 'B', i);
		CHECK_EQUAL(true, db->insert(B));

		setID(C.m_data, 'C', i);
		CHECK_EQUAL(true, db->insert(C));

		setID(D.m_data, 'D', i);
		CHECK_EQUAL(true, db->insert(D));

		setID(E.m_data, 'E', i);
		CHECK_EQUAL(true, db->insert(E));

		setID(F.m_data, 'F', i);
		CHECK_EQUAL(true, db->insert(F));
    }

    CHECK_EQUAL(6*N, db->size());

    addDBFunctions(*db);

    //
	// control: no distance calculation
    //
	string sql_query0 = sqlFullQuery("", "dataset = 'BADDIES'");

	Timer t;
	ProfileRangeDB prdb_control(db, sql_query0);
	t.stop();
	cout << "control: count profiles = " << t << endl;

	// should have found everything
	CHECK_EQUAL(6*N, prdb_control.size());

	t.start();
	prdb_control.readMyData();
	t.stop();
	cout << "control: read profiles = " << t << endl;

    //
	// test geoBox
    //
	string sql_query1 = sqlFullQuery(geoBoxQuery("60", "-3", "111"), "");

	t.start();
	ProfileRangeDB prdb_box(db, sql_query1);
	t.stop();
	cout << "geoBox: count profiles = " << t << endl;

	// should have found B, C, E
	CHECK_EQUAL(3*N, prdb_box.size());

	t.start();
	prdb_box.readMyData();
	t.stop();
	cout << "geoBox: read profiles = " << t << endl;

	//
	// time geoCircle
	//
	string sql_query2 = sqlFullQuery(geoCircleQuery("60", "-3", "111"), "");

	t.start();
	ProfileRangeDB prdb_dist(db, sql_query2);
	t.stop();
	cout << "geoCircle: count profiles = " << t << endl;

	// should have found C, E
	CHECK_EQUAL(2*N, prdb_dist.size());

	t.start();
	prdb_dist.readMyData();
	t.stop();
	cout << "geoCircle: read profiles = " << t << endl;

}
#endif

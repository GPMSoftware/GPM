/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * InputPanel.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "InputPanel.h"
#include "gmatch.h"
#include "MetadataDialog.h"
#include "MyFrame.h"
#include "MegaChoiceDialog.h"

#include "fand/loci.h"
#include "fand/Parser.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("InputPanel");
#include "fand/messages.h"

#include "wx/statline.h"

#include <sstream>

using namespace std;

extern boost::shared_ptr<Database> db;

// example code
//BEGIN_EVENT_TABLE(MyWindow, wxWindow)
//  EVT_ERASE_BACKGROUND(MyWindow::OnErase)
//END_EVENT_TABLE()

BEGIN_EVENT_TABLE(InputGrid, wxGrid)
//    EVT_CHAR(InputGrid::OnChar)
  EVT_ERASE_BACKGROUND(InputGrid::OnErase) // this does not work
END_EVENT_TABLE();

BEGIN_EVENT_TABLE(InputPanel, wxPanel)
    EVT_BUTTON(ID_META,   InputPanel::OnMeta)
    EVT_BUTTON(ID_BROWSE, InputPanel::OnBrowse)
    EVT_BUTTON(ID_NEW,    InputPanel::OnNew)
    EVT_BUTTON(ID_EDIT,   InputPanel::OnEdit)
    EVT_BUTTON(ID_DELETE, InputPanel::OnDelete)
    EVT_BUTTON(ID_CLEAR,  InputPanel::OnClear)
    EVT_BUTTON(ID_SAVE,   InputPanel::OnSave)
    EVT_GRID_CELL_CHANGE(InputPanel::OnEditCell)
    EVT_CHOICE(ID_KIT,    InputPanel::OnSelectKit)
//    EVT_GRID_SELECT_CELL(InputPanel::OnSelectCell)
END_EVENT_TABLE();

static bool wayToSort(int i, int j) { return i > j; }

void
InputGrid::setup()
{
    int grid_cols  = 4 + num_loci;
    int grid_width = GetSize().GetWidth();
    int col_width = grid_width/grid_cols;

    CreateGrid(3, grid_cols);

    wxGridCellTextEditor *editor = new wxGridCellTextEditor() ;
    SetDefaultEditor(editor);

    wxFont font = GetDefaultCellFont();
    SetLabelFont(font); // *wxSMALL_FONT

    SetDefaultColSize(col_width);
    SetColSize( 0, 80 );
    SetColLabelValue(0, _("Sample\nID"));
    SetColSize( 1, 80 );
    SetColLabelValue(1, _("Profile\nID"));
    SetColLabelValue(2, _("N in\nmix"));
    SetColFormatNumber(2);
    SetColLabelValue(3, _("Allele"));

    SetColLabelSize(60);
    SetColLabelTextOrientation(wxVERTICAL);

    // Set up mapping for loci and write Locus labels
    // (Default kit is Identifiler)
    mapLoci(identifiler_loci, num_identifiler_loci);

    SetRowLabelSize(0); // no row labels

}

void
InputGrid::OnErase(wxEraseEvent& event) // this function does not get called
{
	std::cout << "InputGrid::OnErase()" << std::endl;

    wxClientDC* clientDC = NULL;
    if (!event.GetDC())
        clientDC = new wxClientDC(this);

    wxDC* dc = clientDC ? clientDC : event.GetDC() ;

//    wxSize sz = GetClientSize();
//    wxEffects effects;
//    effects.TileBitmap(wxRect(0, 0, sz.x, sz.y), *dc, m_bitmap);

	wxPen pen(*wxRED, 1); // red pen of width 1
	dc->SetPen(pen);
	dc->DrawRectangle(wxRect(0,0,10,10));

    if (clientDC)
        delete clientDC;
}

void
InputGrid::mapLoci(Locus *kit_loci, int kit_len)
{
	// save allele values (using old m_gridLoc)
	std::vector< std::vector<string> > values;

	if (!m_gridLoc.empty())
	{
	    for (int i = 0; i < GetNumberRows(); ++i)
	    {
	    	values.push_back(std::vector<string>());

	    	for (int loc=0; loc<num_loci; ++loc)
	    	{
	    		int g = m_gridLoc[loc];
	    		values[i].push_back( std::string(GetCellValue(i, g+4)) );
	    	}
	    }
	}

	std::vector<int> tmp; // tmp[g] = loc

	tmp.clear();

	// start with the desired kit
	for (int i=0; i<kit_len; ++i)
	{
		tmp.push_back(kit_loci[i]);
	}

	// add other loci at the end
	for (int i=0; i<num_loci; ++i)
	{
		if (std::find(tmp.begin(), tmp.end(), i) == tmp.end())
		{
			tmp.push_back(i);
		}
	}

	Assert(tmp.size() == num_loci);
	m_gridLoc.resize(num_loci);     // m_gridloc[loc] = g;

	for (int i=0; i<num_loci; ++i)
	{
		m_gridLoc[tmp[i]] = i;
	}

	// (re-)write locus labels
	for (int g=0; g<num_loci; ++g)
	{
		wxString str(locus_name[tmp[g]], *wxConvCurrent);

//		if (g >= kit_len) // Alleles not part of the kit
//		{
//			std::transform(str.begin(), str.end(), str.begin(), ::tolower);
//		}
		SetColLabelValue(4+g, str);
	}

#if 0
	// draw a bar at the top of the column labels
    wxWindow* colLabelWindow = GetGridColLabelWindow();

	wxClientDC dc(colLabelWindow);
	wxPen pen(*wxRED, 1); // red pen of width 1
	dc.SetPen(pen);
	dc.DrawPoint(wxPoint(1,1));
	dc.DrawRectangle(wxRect(0,0,10,10));
	dc.SetPen(wxNullPen);

	colLabelWindow->SetEventHandler();

#endif

//	colLabelWindow->Refresh();

//    colLabelWindow->SetBackgroundColour(wxColour(235, 255, 235));

    //    AutoSizeColOrRow(0, true, wxGRID_COLUMN);

	// restore allele values (using new m_gridLoc)
	if (!values.empty())
	{
	    for (int i = 0; i < GetNumberRows(); ++i)
	    {
	    	for (int loc=0; loc<num_loci; ++loc)
	    	{
	    		int g = m_gridLoc[loc];
	    		std::string val = values[i][loc];
	    	    SetCellValue(wxGridCellCoords(i, g+4), wxString(val.c_str()));
	    	}
	    }
	}
}

void
InputPanel::setup()
{
// Widget Hierarchy
//    wxSashLayoutWindow* sash =
//        new wxSashLayoutWindow(parent, ID_INPUT_AREA,
//                               wxDefaultPosition, wxSize(200, 30),
//                               wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
//
//    sash->SetMinSize(wxSize(1000, 200));
//    sash->SetDefaultSize(wxSize(1000, 200));
//    sash->SetOrientation(wxLAYOUT_HORIZONTAL);
//    sash->SetAlignment(wxLAYOUT_TOP);
//    sash->SetBackgroundColour(wxColour(255, 0, 0));
//    sash->SetSashVisible(wxSASH_BOTTOM, true);

    wxPanel *panel = this;

    wxStaticText *label = new wxStaticText(panel, wxID_ANY, _("Dataset"));
    m_dataset           = new wxComboBox(panel, wxID_ANY, _(""));
    m_dataset->Enable(false);

    m_meta  = new wxButton(panel, ID_META, _("Metadata ...")/*, wxDefaultPosition, wxSize(100,10)*/);

    m_profileset     = new wxChoice(panel, wxID_ANY);
    vector<string> psets;
    psets.push_back("All Profiles");
    psets.push_back("C-Set");
    psets.push_back("R-Set");
    setItemContainerStrings(m_profileset, psets);
    m_profileset->SetSelection(0);

    wxButton* browse_profiles = new wxButton(panel, ID_BROWSE, _("Select Profiles from..."));

    wxButton* new_s = new wxButton(panel, ID_NEW, _("New Sample"));
    m_edit  = new wxButton(panel, ID_EDIT, _("Edit"));
    m_del   = new wxButton(panel, ID_DELETE, _("Delete"));
    m_clear = new wxButton(panel, ID_CLEAR, _("Clear"));
    m_save = new wxButton(panel, ID_SAVE, _("Save"));

    wxStaticLine* line1 = new wxStaticLine(panel, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );

    int grid_width = 1350; // 1245;
    m_grid = new InputGrid(panel, wxID_ANY, wxDefaultPosition, wxSize(grid_width, wxDefaultCoord));

    wxStaticText *label2 = new wxStaticText(panel, wxID_ANY, _("STR Kit"));

    m_kit     = new wxChoice(panel, ID_KIT);
    vector<string> kits;
    kits.push_back("Identifiler");
    kits.push_back("PowerPlex 16");
    kits.push_back("PowerPlex ESI 17");
    kits.push_back("NGM Select");
    kits.push_back("GlobalFiler");
    setItemContainerStrings(m_kit, kits);
    m_kit->SetSelection(0);

    // Sizer Hierarchy
    wxBoxSizer *vSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *hSizer = new wxBoxSizer( wxHORIZONTAL );

    vSizer->Add(hSizer, 0, wxEXPAND | wxALL/*, 10*/);

		hSizer->Add(browse_profiles, 0, wxALIGN_CENTER | wxALL, 2);
		hSizer->Add(m_profileset, 0, wxALIGN_CENTER | wxALL, 2);

		hSizer->AddSpacer(10);

		hSizer->Add(new_s,  0, wxALIGN_CENTER | wxALL, 2);
		hSizer->Add(m_edit,   0, wxALIGN_CENTER | wxALL, 2);
		hSizer->Add(m_del,    0, wxALIGN_CENTER | wxALL, 2);
		hSizer->Add(m_clear, 0, wxALIGN_CENTER | wxALL, 2);
		hSizer->Add(m_save,  0, wxALIGN_CENTER | wxALL, 2);

		hSizer->Add(1, 1,   1, wxEXPAND | wxALL, 2);

		hSizer->Add(label2,  0, wxALIGN_CENTER | wxALL, 5);
		hSizer->Add(m_kit, 0, wxALIGN_CENTER | wxALL, 2);

		hSizer->Add(1, 1,   1, wxEXPAND | wxALL, 2);

		hSizer->Add(label,  0, wxALIGN_CENTER | wxALL, 5);
		hSizer->Add(m_dataset,  0, wxALIGN_CENTER | wxALL, 5);
		hSizer->Add(m_meta,   0, wxALIGN_CENTER | wxALL, 5);

	vSizer->Add(line1, 0, wxEXPAND | wxALL, 0);
	vSizer->Add(m_grid, 1, wxEXPAND | wxALL, 10);

    panel->SetSizer(vSizer);

    refresh(); // too early to have any effect?

    enableOps(false);
}

void
InputPanel::setReadOnly(int i, int j, bool ro)
{
    m_grid->SetReadOnly(i, j, ro);

    wxColour col = readOnlyGrey;
    if (ro == false)
    {
        col = m_grid->GetCellValue(i,j).empty() ? wxColour(235, 255, 235) /*green*/: *wxWHITE;
    }
    m_grid->SetCellBackgroundColour(i, j, col);
}

void
InputPanel::setReadOnly(bool ro)
{
    for (int i = 0; i < m_grid->GetNumberRows(); ++i)
    {
        for (int j = 0; j < m_grid->GetNumberCols(); ++j)
        {
            bool read_only = (j==3 ? true : ro); // allele count col is always read only
            setReadOnly(i, j, read_only);
        }
    }
    m_grid->ForceRefresh();
}

// highlight cell and make it the only editable cell
void
InputPanel::setCellError(int i, int j)
{
    setReadOnly(true);

    m_grid->SetReadOnly(i, j, false);
    m_grid->SetCellBackgroundColour(i, j, wxColour(255, 225, 225)); // pink
}

bool
InputPanel::displayProfile(DBProfile const &dbp)
{
    // display single profile and put input area into view mode
	vector<DBProfile> dbp_vec;
	dbp_vec.push_back(dbp);
	return displayProfiles(dbp_vec);
}

bool
InputPanel::displayProfiles(vector<DBProfile> const &dbp_vec)
{
	// display vector of profiles and put input area into view mode
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_VIEW:
        return viewMode(dbp_vec);
        break;
    case MODE_EDIT:
    case MODE_NEW:
        if (wxYES == wxMessageBox(_("Discard current profile?"), _("Display Profile"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
            return viewMode(dbp_vec);
        }
        break;
    default:
        Assert2(0, "Mode not handled");
    }

    return false;
}

bool
InputPanel::viewMode(vector<DBProfile> const &dbp_vec)
{
	// for multiple profiles, insert them in the grid in reverse order
	// (because each new one goes at the top and pushes the rest down)
	bool is_first = true;
	int size = dbp_vec.size();
	for (int i = size-1; i >= 0; --i) // yes really must be an int!
	{
		viewMode(dbp_vec[i], is_first);
		is_first = false;
	}

	enableOps(size == 1);

	return true;
}

void
InputPanel::enableOps(bool enable)
{
    m_meta->Enable(enable);
    m_edit->Enable(enable);
    m_clear->Enable(enable);
    m_del->Enable(enable);
    m_save->Enable(false);
}

void
InputPanel::enableSave(bool enable)
{
    m_save->Enable(enable);
}

void
InputPanel::setAlleleTextByGrid(int allele, int g, std::string val)
{
    m_grid->SetCellValue(wxGridCellCoords(allele, g+4), wxString(val.c_str()));
}

void
InputPanel::setAlleleTextByLocus(int allele, int loc, std::string val)
{
	int g = m_grid->gPos(loc);
	setAlleleTextByGrid(allele, g, val);
}

std::string
InputPanel::getAlleleTextByGrid(int allele, int g)
{
    return std::string(m_grid->GetCellValue(allele, g+4));
}

std::string
InputPanel::getAlleleTextByLocus(int allele, int loc)
{
	int g = m_grid->gPos(loc);
    return getAlleleTextByGrid(allele, g);
}

bool
InputPanel::viewMode(DBProfile const &dbp, bool is_first)
{
    m_num_profiles = 1;
    m_num_contribs.clear();
    m_num_contribs.push_back(dbp.m_data.m_num_contributors);

    // set number of rows
    int nrows = 2*m_num_contribs[0]+1; // rows needed by this profile
    if (is_first)
    {
		m_grid->ClearGrid();
		setNumberRows(nrows);

		// record metadata
		m_metadata = dbp.m_metadata;
    }
    else
    {
    	// insert rows at the beginning for this profile, pushing previous profiles down
        m_grid->InsertRows(0, nrows);
    }

    // ID etc
    int col = 0;
    m_grid->SetCellValue(wxGridCellCoords(0,col++), wxString(dbp.m_data.m_sample_id.c_str(), *wxConvCurrent));
    m_grid->SetCellValue(wxGridCellCoords(0,col++), wxString(dbp.m_data.m_profile_id.c_str(), *wxConvCurrent));
    m_grid->SetCellValue(wxGridCellCoords(0,col++), wxString::Format(_T("%d"), m_num_contribs[0]));

    // Allele (count) column
    for (int i=0; i<2*m_num_contribs[0]; ++i)
    {
        m_grid->SetCellValue(i, col, wxString::Format(_T("%d"), i+1));
    }
    col++;

    // Data

    for (int loc=0; loc<num_loci; ++loc)
    {
        vector<string> allele_strings;
        DBProfile::LocusMap::const_iterator it = dbp.m_b11text.find((Locus)loc);

        if (it != dbp.m_b11text.end())
        {
            allele_strings = it->second;
        }

        for (int j=0; j<2*m_num_contribs[0]; ++j)
        {
           string data = (size_t)j<allele_strings.size() ? allele_strings[j] : "F";
		   setAlleleTextByLocus(j, loc, data);
        }
    }

    setReadOnly(true);

    m_dataset->SetValue(wxString(dbp.m_data.m_dataset.c_str(), *wxConvCurrent));
    m_dataset->Enable(false);

    m_mode = MODE_VIEW;

    return true;
}

void
InputPanel::refresh()
{
    set<string> datasets;
    if (!db->listDatasets(datasets))
    {
        error << startl << "InputPanel::refresh(): Datasets could not be read: " << endl;
        return;
    }

    setItemContainerStrings(m_dataset, datasets);
}

string
InputPanel::sampleID()
{
    return string(m_grid->GetCellValue(0,0).ToUTF8());
}

string
InputPanel::profileID(int prof)
{
    int row = profFirstRow(prof);
    if (row < 0)
    {
        return "";
    }
    else
    {
        return string(m_grid->GetCellValue(row,1).ToUTF8());
    }
}

string
InputPanel::dataset()
{
    wxString str = m_dataset->GetValue();

    if (str.size() == 0)
    {
        return "unknown";
    }
    else
    {
        return string(str.ToUTF8());
    }
}

int
InputPanel::numContribs(int prof)
{
    int row = profFirstRow(prof);
    if (row < 0)
    {
        return -1;
    }
    else
    {
        ulong ncontribs;
        m_grid->GetCellValue(row, 2).ToULong(&ncontribs);
        return (int)ncontribs;
    }
}

void
InputPanel::OnMeta(wxCommandEvent& event)
{
    cout << "OnMeta(): Id = " << event.GetId() << endl;

    MetadataDialog dialog(NULL, wxID_ANY);
    dialog.populate(m_metadata);
    dialog.setEditable(m_mode==MODE_EDIT || m_mode==MODE_NEW);

    int ret;
    if ((ret = dialog.ShowModal()) == wxID_OK)
    {
    	dialog.apply();
        dialog.readback(m_metadata);
    }
}

void
InputPanel::OnNew(wxCommandEvent& event)
{
    cout << "OnNew(): Id = " << event.GetId() << endl;
    newSample();
}

void
InputPanel::newSample()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_VIEW:
        newSampleMode();
        break;
    case MODE_EDIT:
    case MODE_NEW:
        if (wxYES == wxMessageBox(_("Discard current profile?"), _("New Sample"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
            newSampleMode();
        }
        break;
    default:
        Assert2(0, "Mode not handled");
    }

    MyFrame::refreshMessageWindow();
}

void
InputPanel::OnBrowse(wxCommandEvent& event)
{
    cout << "OnBrowse(): Id = " << event.GetId() << endl;

    selectProfiles();
}

void
InputPanel::selectProfiles()
{

    string query;
	set<string> vals;

	switch (m_profileset->GetSelection())
    {
    case 0: // all profiles

		// empty query
		break;

    case 1: // C-set
    {
		if (MyFrame::cSet().size() == 0)
        {
            wxMessageBox( _("C-Set: No profiles selected."),
                          _("Find Profile"),
                          wxOK | wxICON_EXCLAMATION, this);
            return;
        }
        else
        {
        	query = MyFrame::cQuery();
        }
        break;
    }
    case 2: // R-set
    {
		if (MyFrame::rSet().size() == 0)
        {
            wxMessageBox( _("R-Set: No profiles selected."),
                          _("Find Profile"),
                          wxOK | wxICON_EXCLAMATION, this);
            return;
		}
		else
		{
			query = MyFrame::rQuery();
		}
        break;
    }
    default:
    	Assert2(false, "InputPanel::OnBrowse: logic error");
    }

	if (!db->listProfileIDs(vals, query))
	{
		error << startl << "InputPanel::OnBrowse(): " << "profile_key" << " could not be read: " << endl;
		MyFrame::refreshMessageWindow();
		return;
	}

    // pop up a MegaChoiceDialog (MULTI mode) to choose a profile key
	MegaChoiceDialog dialog(NULL, wxID_ANY, _("Select Profiles"), vals, MegaChoiceDialog::MULTI);

	int ret;
	if ((ret = dialog.ShowModal()) == wxID_OK)
	{
		wxArrayString selections = dialog.GetStringSelections();

		if (!selections.empty())
		{
			vector<string> p_keys;
			for (size_t i=0; i<selections.size(); ++i)
			{
				p_keys.push_back(string(selections[i].ToUTF8()));
			}

			viewProfiles(p_keys);
		}
	}

    MyFrame::refreshMessageWindow();
}

#if 1
void
InputPanel::viewProfiles(vector<string> const &p_keys)
{
	vector<DBProfile> dbpvec;

	for (size_t i=0; i<p_keys.size(); ++i)
	{
	    ProfileData data;
	    DBProfile dbp(data);
		if (db->readProfile(p_keys[i], dbp))
		{
			dbpvec.push_back(dbp);
		}
		else
		{
			error << startl << "Reading profile: " << p_keys[i] << endl;
		}
	}

	displayProfiles(dbpvec);
}

#else
// This version does all the queries in one go but the results are
// returned in random order (NB COULD USE SORT BY)
//
void
InputPanel::viewProfiles(vector<string> const &p_keys)
{
	string sql_query;

	for (size_t i=0; i<p_keys.size(); ++i)
	{
		if (sql_query.empty())
		{
			sql_query = "select * from Profiles WHERE ";
		}
		else
		{
			sql_query += " or ";
		}

		sql_query += string("profile_key") + " = '" + p_keys[i] + "'";
	}

	vector<DBProfile> dbpvec;
	if (db.readQuery(sql_query, dbpvec))
	{
		viewMode(dbpvec);
	}
	else
	{
		error << startl << "Reading query: " << sql_query << endl;
	}
}
#endif

int
InputPanel::profFirstRow(int nprof)
{
    if (nprof < 0 || nprof >= m_num_profiles) return -1;

    int n = 0;
    for (int i=0; i<nprof; ++i)
    {
        n += 2 * m_num_contribs[i];
    }

    return n;
}

int
InputPanel::profLastRow(int nprof)
{
    if (nprof < 0 || nprof >= m_num_profiles) return -1;

    return profFirstRow(nprof) + 2 * m_num_contribs[nprof] - 1;
}

pair<int, int>
InputPanel::profAndAlleleAtRow(int row)
{
    int prof = -1, allele = -1;

    int n = 0;
    for (int i=0; i<m_num_profiles; ++i)
    {
        n += 2 * m_num_contribs[i]; // one past profRowLast(i)

        if (row < n)
        {
            prof = i;
            allele = n - row;
            break;
        }
    }

    return make_pair(prof, allele);
}

// user has changed the number of contributors for a profile
void
InputPanel::changeNumContribs(int nprof, int ncontribs)
{
    if (ncontribs > m_num_contribs[nprof]) // increased
    {
        int nadd = 2 * (ncontribs - m_num_contribs[nprof]);

        if (nprof == m_num_profiles - 1)
        {
            m_grid->AppendRows(nadd);
        }
        else
        {
            m_grid->InsertRows(profLastRow(nprof) + 1, nadd);
        }

    }
    else if (ncontribs < m_num_contribs[nprof]) // decreased
    {
        int ndel = 2 * (m_num_contribs[nprof] - ncontribs);
        m_grid->DeleteRows(profLastRow(nprof) - ndel + 1, ndel);
    }

    m_num_contribs[nprof] = ncontribs;

    enableCellsForEdit(m_mode == MODE_NEW);
}

// user has typed in a new profile ID in the last available position
void
InputPanel::newProfile()
{
    m_num_profiles++;
    m_num_contribs.push_back(1);

    // add rows to grid
    m_grid->AppendRows(2);

    enableCellsForEdit(true);
}

// Set cells (non-)editable as required. Also update the allele count fields
// The new_sample flags indiates we are entering a new sample, as opposed to
// editing an existing profile. In this mode multiple profiles are allowed.
void
InputPanel::enableCellsForEdit(bool new_sample)
{
    setReadOnly(true);      // entire grid disabled unless enabled below

    if (new_sample)
    {
    	// enable entry of Sample ID only
        setReadOnly(0,0,false); // Sample ID
    }
    else
    {
    	// enable editing of Sample and Profile IDs
        setReadOnly(0,0,false); // Sample ID
        setReadOnly(0,1,false); // Profile ID
    }

    if (!sampleID().empty())
    {
        for (int prof=0; prof<m_num_profiles; ++prof)
        {
            int row = profFirstRow(prof);
            if (new_sample)
            {
                setReadOnly(row,1,false); // Profile ID
            }

            // ncontribs
            m_grid->SetCellValue(wxGridCellCoords(row,2), wxString::Format(_T("%d"), m_num_contribs[prof]));
            setReadOnly(row,2,false);

            for (int k = 0; k<(2 * m_num_contribs[prof]); ++k)
            {
                // set the allele counter
                m_grid->SetCellValue(wxGridCellCoords(row+k,3), wxString::Format(_T("%d"), k+1));

                for (int g = 0; g<num_loci; ++g)
                {
                    setReadOnly(row+k,g+4,false); // Data

                    // if cell is empty, enter a default value of 'F'
					if (getAlleleTextByGrid(row+k, g).empty())
                    {
						setAlleleTextByGrid(row+k, g, "F");
                    }
                }
            }
        }

        // the 'new profile' cell
        if (new_sample)
        {
            int row = 0;
            if (m_num_profiles)
            {
                row = profLastRow(m_num_profiles - 1) + 1;
            }
            Assert2(row < m_grid->GetNumberRows(), "wrong number of rows");
            setReadOnly(row,1,false); // Profile ID
        }
    }

    m_grid->ForceRefresh();

    m_dataset->Enable(true);
}

void
InputPanel::setNumberRows(int rows)
{
	m_grid->setNumberRows(rows);
}

void
InputPanel::newSampleMode()
{
    // set up grid for new sample
    // Initially only the Sample ID box is editable

    m_grid->ClearGrid();
    m_num_profiles = 0;
    m_num_contribs.clear();
    m_num_contribs.push_back(1);

    // set number of rows to 3
    setNumberRows(3);

    bool new_sample = true;
    enableCellsForEdit(new_sample);

    m_dataset->Clear();
    m_grid->SetGridCursor(0, 0);
    m_grid->EnableCellEditControl();
    m_grid->ForceRefresh();
    m_grid->SetFocus();

	enableOps(true);
	enableSave(true);
    refresh();

    m_mode = MODE_NEW;
}

void
InputPanel::OnEdit(wxCommandEvent& event)
{
    cout << "OnEdit(): Id = " << event.GetId() << endl;
    editProfile();
}


void
InputPanel::editProfile()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_EDIT:
    case MODE_NEW:
        // null operation
        break;
    case MODE_VIEW:
        editMode();
        break;
    default:
        Assert2(0, "Mode not handled");
    }

    MyFrame::refreshMessageWindow();
}

void
InputPanel::editMode()
{
    bool new_sample;
    enableCellsForEdit(new_sample = false);
    refresh();
    enableSave(true);
    m_mode = MODE_EDIT;
}

void
InputPanel::OnDelete(wxCommandEvent& event)
{
    cout << "OnDelete(): Id = " << event.GetId() << endl;
    deleteProfile();
}

void
InputPanel::deleteProfile()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_NEW:
        // null operation
        break;
    case MODE_EDIT:
    case MODE_VIEW:
        doDeleteProfile();
        break;
    default:
        Assert2(0, "Mode not handled");
    }

    MyFrame::refreshMessageWindow();
}

void
InputPanel::doDeleteProfile()
{
    // data
    ProfileData data;
    DBProfile dbp_tmp(data);

    // Construct Key from IDs
    string prof_key = sampleID() + "-" + profileID(0);

    // Check if Profile already exists in database
    bool exists = db->readProfile(prof_key, dbp_tmp);

    // confirm dialog
    if (exists)
    {
        wxString message = wxString::Format(wxT("Delete record %s?"), wxString(prof_key.c_str(), *wxConvCurrent).c_str());
        if (wxNO == wxMessageBox( message,
                      _("Delete"),
                      wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
            return;
        }
    }

    // delete
    if (!db->deleteProfile(prof_key))
    {
        cerr << "Delete Failed for " << prof_key << " (continuing)" << endl;
    }
    else
    {
        emptyMode();
        MyFrame::refresh();
    }
}

void
InputPanel::OnClear(wxCommandEvent& event)
{
    cout << "OnClear(): Id = " << event.GetId() << endl;
    clearProfile();
}

void
InputPanel::clearProfile()
{
    (void)tryClear();
    MyFrame::refreshMessageWindow();
}

void
InputPanel::clearView()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_VIEW:
        (void) emptyMode();
        break;
    case MODE_EDIT:
    case MODE_NEW:
        break;
    default:
        Assert2(0, "Mode not handled");
    }
}

bool
InputPanel::tryClear()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_VIEW:
        return emptyMode();
        break;
    case MODE_EDIT:
    case MODE_NEW:
        if (wxYES == wxMessageBox(_("Discard current profile?"), _("Clear Profile"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
        {
            return emptyMode();
        }
        break;
    default:
        Assert2(0, "Mode not handled");
    }

    return false;
}

bool
InputPanel::emptyMode()
{
    setNumberRows(3);
    m_grid->ClearGrid();
    setReadOnly(true);
    m_dataset->Clear();
    m_dataset->Enable(false);
	enableOps(false);
    m_mode = MODE_EMPTY;
    return true;
}

void
InputPanel::OnSave(wxCommandEvent& event)
{
    cout << "OnSave(): Id = " << event.GetId() << endl;

    saveProfile();
}

void
InputPanel::saveProfile()
{
    switch (m_mode)
    {
    case MODE_EMPTY:
    case MODE_VIEW:
        break;       // nothing to save
    case MODE_EDIT:
    case MODE_NEW:
    {
        if (m_error)
        {
            wxMessageBox( _("Data has syntax errors.\nFix before saving."),
                          _("Save"),
                          wxOK | wxICON_EXCLAMATION, this);
            return;
        }

        // count empty data cells
        int nempty = 0;
        for (int i = 0; i <= profLastRow(m_num_profiles - 1); ++i)
        {
            for (int g = 0; g < num_loci; ++g)
            {
				if (getAlleleTextByGrid(i, g).empty()) nempty++;
            }
        }

        if (nempty) // empty data cells
        {
            wxMessageBox( _("Empty data cells.\nComplete before saving."),
                          _("Save"),
                          wxOK | wxICON_EXCLAMATION, this);
            return;
        }

        saveSample();
        break;
    }
    default:
        Assert2(0, "Mode not handled");
    }

    MyFrame::refreshMessageWindow();
}

// parse allele data (in NEW mode this replaces 'D' statements with {S}@B)
void
InputPanel::transformAlleleData(map<pair<int, int>, string > &allele_data)
{
    for (int loc=0; loc<num_loci; ++loc)
    {
        // Build the set of all alleles in the sample, at this locus (ignore the p values)
        PMF<Allele> sample_all;
        for (int i=0; i <= profLastRow(m_num_profiles - 1); ++i)
        {
            PMF<Allele> pmf;

		    string allele_string = getAlleleTextByLocus(i,loc);
    	    allele_data[make_pair(i,loc)] = allele_string;
            if (Parser::theParser().parse(allele_string, pmf) == Parser::yacc_ok)
            {
                sample_all += pmf;
            }
        }

        // Build a string representing 'D'
        ostringstream dstream;
        int dstr_count = 0;
        dstream << "(";

        PMF<Allele>::iterator it, current;
        for(it = sample_all.begin(); it != sample_all.end();)
        {
            current = it++;

            if (current->first != Allele::unknown) // not 'D'
            {
                dstream << current->first;
                ++dstr_count;

                if (it != sample_all.end())
                {
                    dstream << "/";
                }
            }
        }

        dstream << ")@B";

        // for each row, substitute 'D' and store
        for (int i=0; i <= profLastRow(m_num_profiles - 1); ++i)
        {
            // search for a 'D' (should be first character if present)
            size_t pos = allele_data[make_pair(i,loc)].find('D');
            if (pos != string::npos)
            {
                Assert(pos == 0);
                if (dstr_count > 0)
                {
                    allele_data[make_pair(i,loc)].replace(0, 1, dstream.str());
                }
                else
                {
                    warn << startl << "'D' used with no alleles in the sample (assuming 'F')" << endl;
                    allele_data[make_pair(i,loc)] = "F";
                }
            }
        }
    }
}

// saves the currently edited sample
// (used in both NEW and EDIT modes)
void
InputPanel::saveSample()
{
    map<pair<int, int>, string > allele_data;
    transformAlleleData(allele_data);

    for (int prof=0; prof<m_num_profiles; ++prof)
    {
        ProfileData data;
        DBProfile dbp(data);
        DBProfile dbp_tmp(data);

        // data
        dbp.m_data.m_sample_id           = sampleID();
        dbp.m_data.m_profile_id          = profileID(prof);

        // Construct Key from IDs
        string prof_key                  = dbp.m_data.m_sample_id.substr(0, Database::SAMPLE_ID_LEN)
        		                         + "-"
        		                         + dbp.m_data.m_profile_id.substr(0, Database::PROFILE_ID_LEN);

        dbp.m_data.m_id                  = prof_key;
        dbp.m_data.m_num_contributors = numContribs(prof);

        // TODO
//        dbp.m_data.m_kit_type;
//        dbp.m_data.m_evidence_type;
        dbp.m_data.m_dataset = dataset();
//        dbp.m_data.m_error_rate;

        // save metadata
        dbp.m_metadata = m_metadata;

        // add allele data
        for (int loc=0; loc<num_loci; ++loc)
        {
            for (int i=profFirstRow(prof); i <= profLastRow(prof); ++i)
            {
            	dbp.m_b11text[(Locus)loc].push_back(allele_data[make_pair(i, loc)]);
            }
        }

        // Check if Profile already exists in database
        bool exists = db->readProfile(prof_key, dbp_tmp);

        // confirm dialog
        if (exists)
        {
            wxString message = wxString::Format(wxT("Overwrite record %s?"), wxString(prof_key.c_str(), *wxConvCurrent).c_str());
            if (wxNO == wxMessageBox( message,
                          _("Save"),
                          wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
            {
                break;
            }
        }

        // write to database (overwrite mode)
        if (!db->overwrite(dbp))
        {
            cerr << "Write Failed for " << prof_key << " (continuing)" << endl;
        }
        else
        {
        	ok << "Profile " << prof_key << " saved to database" << endl;
            MyFrame::refresh();
        }
    }
}

bool
InputPanel::checkData(string allele_string)
{
    // don't check empty strings (we have another check for that)
    if (allele_string.empty()) return true;

    // in MODE_EDIT don't allow 'D' constructs
    if (m_mode == MODE_EDIT && allele_string.find("D") != string::npos)
    {
        return false;
    }

    // parse the data
    PMF<Allele> pmf;
    return (Parser::theParser().parse(allele_string, pmf) == Parser::yacc_ok);
}

void
InputPanel::OnEditCell(wxGridEvent& event)
{
    cout << "OnEditCell(): Id = " << event.GetId() << endl;

    // let the Grid process its own events before we try to override its behaviour
//    cout << "calling wxSafeYield ... " << endl;
//    wxSafeYield();
//    cout << "... wxSafeYield finished." << endl;

    // make sure we do not get into an infinite loop
    // by processing edits we have generated from this callback
    // (?? doesn't WxWidgets do this for you?)
    static bool processing = false;

    if (! processing)
    {
        processing = true;

        int i = event.GetRow();
        int j = event.GetCol();

        if (i == 0 && j == 0) // the sample ID field
        {
        	bool new_sample = (m_mode == MODE_NEW);
            enableCellsForEdit(new_sample);

            // DOESNT WORK
            // set focus to profile ID field
//            m_accept_select = true;
//            m_grid->SetGridCursor(0, 1);
//            m_grid->GoToCell(0, 1);
//            m_grid->EnableCellEditControl();
//            m_grid->ForceRefresh();

        }
        else if ( (i == profLastRow(m_num_profiles - 1) + 1) && (j == 1) ) // the last (empty) profile ID field
        {
            if (! m_grid->GetCellValue(i,j).empty()) // now non-empty
            {
                newProfile();
            }
        }
        else if (j == 2) // a ncontribs field
        {
            pair<int, int> prof_allele = profAndAlleleAtRow(i);

            unsigned long ncontribs;
            if(m_grid->GetCellValue(i,j).ToULong(&ncontribs) && ncontribs > 0)
            {
                changeNumContribs(prof_allele.first, (int)ncontribs);
            }
            else
            {
                // restore previous value
                m_grid->SetCellValue(wxGridCellCoords(i,j), wxString::Format(_T("%d"), m_num_contribs[prof_allele.first]));
                wxBell();
            }
        }
        else if (j > 3) // a data field
        {
        	int g = j - 4;
        	string allele_str = getAlleleTextByGrid(i,g);

            // do any key translations
            bool changed = false;
            for (size_t k=0; k<allele_str.size(); ++k)
            {
                if (allele_str[k] == '*')
                {
                    allele_str[k] = '@';
                    changed = true;
                }
                else if (allele_str[k] == 'f' || allele_str[k] == '+')
                {
                    allele_str[k] = 'F';
                    changed = true;
                }
                else if (allele_str[k] == 'd' || allele_str[k] == '-')
                {
                    allele_str[k] = 'D';
                    changed = true;
                }
            }

            if (changed)
            {
            	setAlleleTextByGrid(i,g, allele_str);
            }

            if (checkData(allele_str))
            {
                // clear error state
                // NB refresh cells even if we are not currently in error state - colours may have changed.
            	bool new_sample = (m_mode == MODE_NEW);
                enableCellsForEdit(new_sample);
                m_error = false;
            }
            else
            {
                // user has entered an error: set error state
                setCellError(i, j);
                wxBell();
                m_error = true;
            }
        }
    }
    processing = false;

    MyFrame::refreshMessageWindow();
}

void
InputPanel::OnSelectKit(wxCommandEvent& event)
{
    cout << "OnSelectKit(): Id = " << event.GetId() << endl;

	switch (m_kit->GetSelection())
    {
    case 0: // Identifiler
		m_grid->mapLoci(identifiler_loci, num_identifiler_loci);
		break;

    case 1: // PowerPlex 16
		m_grid->mapLoci(powerplex16_loci, num_powerplex16_loci);
		break;

    case 2: // PowerPlex ESI 17
		m_grid->mapLoci(powerplexesi17_loci, num_powerplexesi17_loci);
		break;

    case 3: // NGM Select
		m_grid->mapLoci(ngm_loci, num_ngm_loci);
		break;

    case 4: // GlobalFiler
		m_grid->mapLoci(globalfiler_loci, num_globalfiler_loci);
		break;

    default:
    	;
    }
}

#if 0
void
InputPanel::OnSelectCell(wxGridEvent& event)
{
    int row = m_grid->GetGridCursorRow();
    int col = m_grid->GetGridCursorCol();

    int i = event.GetRow();
    int j = event.GetCol();
    wxPoint p = event.GetPosition();

    cout << "OnSelectCell(): Id = " << event.GetId() << " in (" << row << ", " << col << ")" << endl;
    cout << "Selecting ( " << i << ", " << j << ") at (" << p.x << ", " << p.y << ")" << endl;

    if (! m_grid->IsReadOnly(i, j))
    {
        m_grid->ShowCellEditControl();
    }

    event.Skip();
}
#endif

void
InputGrid::OnChar(wxKeyEvent& event)
{
    cout << "OnChar(): Id = " << event.GetId() << endl;
}


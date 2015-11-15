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

#include "ResultsPanel.h"
#include "gmatch.h"
#include "MyFrame.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("ResultsPanel");
#include "fand/messages.h"

const int DEFAULT_ROWS = 11;

using namespace std;

BEGIN_EVENT_TABLE(ResultsPanel, wxPanel)
	EVT_RADIOBOX(ID_FORMAT, ResultsPanel::OnFormat)
	EVT_GRID_RANGE_SELECT(ResultsPanel::OnSelectMatch)
    EVT_BUTTON(ID_SAVE_RESULT, ResultsPanel::OnSave)
    EVT_BUTTON(ID_CLEAR_RESULT, ResultsPanel::OnClear)
END_EVENT_TABLE();

void MatchGrid::setup()
{
	CreateGrid(DEFAULT_ROWS, 3, wxGrid::wxGridSelectRows);
	EnableEditing(false);

	SetRowLabelSize(60);
	SetDefaultColSize(100);

	SetColLabelValue(0, _("C-Profile"));
	SetColLabelValue(1, _("R-Profile"));
	SetColLabelValue(2, _("Max LR"));
}

int
MatchGrid::GetSelection()
{
    wxArrayInt rows = GetSelectedRows();

    cout << "GetSelection: rows.size() = " << rows.size() << endl;
    if (rows.size() != 1)
    {
    	return -1;
    }
    else
    {
    	return rows[0];
    }
}

void
ResultGrid::setup()
{
	CreateGrid(DEFAULT_ROWS, 2, wxGrid::wxGridSelectRows);
	EnableEditing(false);

//	HideRowLabels();
    SetRowLabelSize(0); // no row labels

	SetDefaultColSize(180);
//	SetColSize(0, 70);
//	SetColSize(1, 290); // space for 40 characters in the chosen font

	SetColLabelValue(0, _("Match"));
	SetColLabelValue(1, _("Likelihood ratio"));

	SetDefaultCellAlignment(wxALIGN_RIGHT, wxALIGN_CENTRE);
}

void
ResultsPanel::setup(wxString const &title)
{
    // Widget Hierarchy

    wxPanel *panel = this;

    wxStaticText *label = new wxStaticText(panel, wxID_ANY, title);

	wxArrayString bstrings;
	bstrings.Add(wxT("Scientific"));
	bstrings.Add(wxT("Fixed"));

	m_format = new wxRadioBox(panel, ID_FORMAT, _("Numeric format"), wxDefaultPosition, wxDefaultSize, bstrings, wxRA_SPECIFY_COLS, wxRA_HORIZONTAL);
	m_format->SetSelection(0);

	// set the size explicitly because resizing fails horribly
    int grid_width = 375;
    int grid_height = 275;
    m_matchList         = new MatchGrid(panel, ID_MATCHLIST, wxDefaultPosition, wxSize(grid_width, grid_height) );

    m_resultsList       = new ResultGrid(panel, ID_RESULTLIST, wxDefaultPosition, wxDefaultSize);

    m_save              = new wxButton(panel, ID_SAVE_RESULT, _("Save"));
    m_save->Enable(false);
    m_clear             = new wxButton(panel, ID_CLEAR_RESULT, _("Clear"));
    m_clear->Enable(false);

    // set a fixed width font
    wxFont font(
    		9,                  // size
    		wxFONTFAMILY_SWISS, // family
    		wxNORMAL,           // style
    		wxNORMAL,           // weight
    		false,              // underlined
    		_("Monospace")      // face
    		);

    m_matchList->SetFont(font);
    m_resultsList->SetFont(font);

    // Sizer Hierarchy
    wxBoxSizer *vSizer  = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *hSizer0  = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *hSizer1  = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *hSizer2  = new wxBoxSizer(wxHORIZONTAL);

    vSizer->Add(hSizer0,           0, wxEXPAND | wxALL, 0);
    hSizer0->AddStretchSpacer(14);
    hSizer0->Add(label,      0, wxALIGN_CENTER | wxALL, 5);
    hSizer0->AddStretchSpacer(8);
    hSizer0->Add(m_format,   0, wxALIGN_CENTER | wxALL, 5);
    hSizer0->AddStretchSpacer(1);

    vSizer->Add(hSizer1,           0, wxEXPAND | wxALL, 0);
    hSizer1->Add(m_matchList,      0, wxEXPAND | wxALL, 5);
    hSizer1->Add(m_resultsList,    0, wxEXPAND | wxALL, 5);

    vSizer->Add(hSizer2,          0, wxALIGN_RIGHT | wxALL, 5);
    hSizer2->Add(m_save,          1, wxEXPAND | wxLEFT, 5);
    hSizer2->Add(m_clear,         1, wxEXPAND | wxLEFT, 5);

    panel->SetSizer(vSizer);

    m_got_results = false;
    m_saved_results = false;
}

bool operator<(ResultDisplayItem const &rdi1, ResultDisplayItem const &rdi2)
{
    return rdi1.max_lr < rdi2.max_lr;
}

bool rdiComp(ResultDisplayItem const &rdi1, ResultDisplayItem const &rdi2)
{
    return rdi2 < rdi1; // sort in descending order
}

void
ResultsPanel::displayMatches(
        ProfileRange const    &cset,
        ProfileRange const    &rset,
        vector<MatchType>     &match_types,
        float                  lr,
        NMmatchResults const  &matches)
{
    m_rdi_array.clear();
    m_c_ids.clear();
    m_r_ids.clear();
    m_match_types = match_types;
    m_lr_thresh = lr;

    // first create the indices with "" for the ids
    NMmatchResults::const_iterator it;
    for(it = matches.begin(); it != matches.end(); ++it)
    {
        ResultDisplayItem rdi;

    	rdi.c_index = it->first.first;
    	rdi.r_index = it->first.second;

        m_c_ids[rdi.c_index] = "";
        m_r_ids[rdi.r_index] = "";

        const vector<Result> &res = it->second;

        rdi.max_lr = 0;
        for(size_t i=0; i<res.size(); ++i)
        {
            if (res[i].likelihood_ratio > rdi.max_lr)
            {
                rdi.max_lr = res[i].likelihood_ratio;
            }
        }

        m_rdi_array.push_back(rdi);

    }

    // Now read the profiles from the database in chunks and fill in the profile IDs.
    static int NPROF = 1000;
    getProfileNames(m_c_ids, cset, NPROF);
    getProfileNames(m_r_ids, rset, NPROF);

    // now step through the RDIs and fill in the details
    std::vector<ResultDisplayItem>::iterator it2;
    for(it2 = m_rdi_array.begin(); it2 != m_rdi_array.end(); ++it2)
    {
        string p1_id = m_c_ids[it2->c_index];
        string p2_id = m_r_ids[it2->r_index];

        it2->c_string = p1_id.c_str();
        it2->r_string = p2_id.c_str();

		std::ostringstream ss;
		ss << setprecision(3) << /*fixed <<*/ it2->max_lr;
		it2->lr_string = ss.str();
    }

    // sort on max_lr
    sort(m_rdi_array.begin(), m_rdi_array.end(), rdiComp);

    m_matchList->setMatches(m_rdi_array);

    m_resultsList->Clear();

    m_save->Enable(true);
    m_clear->Enable(true);

    m_got_results = true;
    m_saved_results = false;
}

void
ResultsPanel::writeMatches(
        ofstream              &ofs,
        NMmatchResults  const &matches)
{
	MatchResults results(m_match_types, matches, m_c_ids, m_r_ids, m_lr_thresh);
	outputMatches(ofs, results);
}

bool compResult(Result const &r1, Result const &r2)
{
    return r1.likelihood_ratio > r2.likelihood_ratio;
}

void
ResultsPanel::OnSelectMatch(wxGridRangeSelectEvent& event)
{
	int id = event.GetId();

    cout << "OnSelectMatch(): Id = " << event.GetId() << endl;

    if (id == ID_MATCHLIST)
    {
    	if (m_got_results) DisplayResults(formatSci());
    }
    else if (id == ID_RESULTLIST)
    {
    	if (m_got_results) SelectResult();
    }
    else
    {
    	Assert2(false, "unknown grid event");
    }
}

void
ResultsPanel::DisplayResults(bool sci)
{
    // find which match has been selected
    int sel = m_matchList->GetSelection();
    cout << "sel = " << sel << endl;

    if (sel<0 || sel > ((int)m_rdi_array.size())-1) return;

    // get the result
    int c = m_rdi_array[sel].c_index;
    int r = m_rdi_array[sel].r_index;

    vector<Result> res = MyFrame::m_match_results[make_pair(c, r)];

    //
    // Display the results in the result window
    //

    // sort the results on lr
//    sort(res.begin(), res.end(), compResult); // sorting disabled ! see below

    if (sci)
    {
    	m_resultsList->SetColSize(0, 180);
    	m_resultsList->SetColSize(1, 180);
    }
    else
    {
    	m_resultsList->SetColSize(0, 70);
    	m_resultsList->SetColSize(1, 290); // space for 40 characters in the chosen font
    }

    // display
    // NB we assume the match_types array is in the same order as the results array

    vector<pair<string, string> > result_strings;
    int j = 0;
    for (size_t i=0; i<m_match_types.size(); ++i)
    {
        char buf1[SIZE], buf2[SIZE];
        if (res[j].match_type.m_rel_type == m_match_types[i].m_rel_type)
        {
            snprintf(buf1, SIZE, "%s", m_match_types[i].string().c_str());
            if (sci)
            {
            	snprintf(buf2, SIZE, "%9.3g", res[j++].likelihood_ratio); // scientific
            }
            else
            {
            	snprintf(buf2, SIZE, "%0.0f", res[j++].likelihood_ratio); // fixed
            }
        }
        else // no result reported for i
        {
            snprintf(buf1, SIZE, "%s", m_match_types[i].string().c_str());
            if (m_lr_thresh > 0)
            {
                snprintf(buf2, SIZE, "<%6.2g", m_lr_thresh);
            }
            else
            {
                snprintf(buf2, SIZE, "%8.2g", 0.0);
            }
        }

        result_strings.push_back(make_pair(string(buf1), string(buf2)));
    }

    m_resultsList->setResults(result_strings);

    //
    // Display the profiles in the Input window
    //
    cout << m_c_ids[c] << " " << m_r_ids[r] << endl;
	vector<string> p_keys;
	p_keys.push_back(m_c_ids[c]);
	p_keys.push_back(m_r_ids[r]);
    MyFrame::m_inputArea->viewProfiles(p_keys);
}

void
ResultsPanel::SelectResult()
{
    cout << "SelectResult()" << endl;
}

void
ResultsPanel::OnSave(wxCommandEvent& event)
{
    if (!m_got_results) return;

    // pop up save dialog
    wxString caption = _("Save results");
    wxString widlcard = _("Results (.res.csv)|*.res.csv");
    wxString defaultDir = wxString(homeDir(), *wxConvCurrent);
    wxString defaultFilename = wxEmptyString;

    wxFileDialog dialog(this, caption, defaultDir, defaultFilename, widlcard, wxFD_SAVE);
    if(dialog.ShowModal() == wxID_OK)
    {
        wxString path = dialog.GetPath();
        string filename(path.ToUTF8());

        // open file
         std::ofstream ofs(filename.c_str());
         if (! ofs.good()){
             error << startl << "error opening file: " << filename << endl;
             return;
         }

        cout << "Saving results to " << filename << "..." << endl;

        writeMatches(ofs, MyFrame::m_match_results);

        m_saved_results = true;
    }
}

void
ResultsPanel::OnClear(wxCommandEvent& event)
{
    cout << "ResultsPanel::OnClear(): Id = " << event.GetId() << endl;
    clear();
}

bool
ResultsPanel::formatSci()
{
    return (m_format->GetSelection() == 0);
}

void
ResultsPanel::OnFormat(wxCommandEvent& event)
{
    cout << "ResultsPanel::OnFormat(): Id = " << event.GetId() << endl;

    DisplayResults(formatSci());
}

void
ResultsPanel::clear()
{
    m_matchList->setNumberRows(DEFAULT_ROWS);
    m_matchList->Clear();

    m_resultsList->setNumberRows(DEFAULT_ROWS);
    m_resultsList->Clear();

    m_save->Enable(false);
    m_clear->Enable(false);

    m_got_results = false;
    m_saved_results = false;
}

void
MatchGrid::setMessage(string const &msg)
{
	Clear();
	setNumberRows(1);
	SetCellValue(0, 0, wxString(msg.c_str(), *wxConvCurrent));
}

void
MatchGrid::setMatches(vector<ResultDisplayItem> const &s)
{
	Clear();
	setNumberRows(s.size());

    vector<ResultDisplayItem>::const_iterator it;
    int row = 0;
    for (it = s.begin(); it != s.end(); ++it, ++row)
    {
        SetCellValue(row, 0, wxString(it->c_string.c_str(), *wxConvCurrent));
        SetCellValue(row, 1, wxString(it->r_string.c_str(), *wxConvCurrent));
        SetCellValue(row, 2, wxString(it->lr_string.c_str(), *wxConvCurrent));
    }
}

void
ResultGrid::setResults(vector<pair<string, string> > const &s)
{
	Clear();
	setNumberRows(s.size());

    vector<pair<string, string> >::const_iterator it;
    int row = 0;
    for (it = s.begin(); it != s.end(); ++it, ++row)
    {
        SetCellValue(row, 0, wxString(it->first.c_str(), *wxConvCurrent));
        SetCellValue(row, 1, wxString(it->second.c_str(), *wxConvCurrent));
    }
}

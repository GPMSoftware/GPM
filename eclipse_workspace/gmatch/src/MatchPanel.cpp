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

#include "MatchPanel.h"
#include "gmatch.h"
#include "MyFrame.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("ResultsPanel");
#include "fand/messages.h"

using namespace std;

BEGIN_EVENT_TABLE(MatchPanel, wxPanel)
    EVT_BUTTON(ID_BROWSE_FDB, MatchPanel::OnBrowse)
END_EVENT_TABLE();

void
MatchPanel::setup(wxString const &title)
{
    // Widget Hierarchy

    wxPanel *panel = this;

    wxStaticText* labelA = new wxStaticText(panel, wxID_ANY, _("Crime delta"));
    wxString s1; s1 << MyFrame::m_match_params->crime_delta;
    m_cdeltaTxt          = new wxTextCtrl(panel, ID_CDELTA, s1);

    wxStaticText* labelB = new wxStaticText(panel, wxID_ANY, _("Reference delta"));
    wxString s2; s2 << MyFrame::m_match_params->ref_delta;
    m_rdeltaTxt          = new wxTextCtrl(panel, ID_RDELTA, s2);

//    wxStaticText* labelS0 = new wxStaticText(panel, wxID_ANY, _("Extra rels"));
//    m_string0Txt         = new wxTextCtrl(panel, ID_STRING0, _(""));
//    m_string1Txt         = new wxTextCtrl(panel, ID_STRING0, _(""));
//    m_string2Txt         = new wxTextCtrl(panel, ID_STRING0, _(""));

    wxStaticText* label1 = new wxStaticText(panel, wxID_ANY, _("LR threshold"));
    wxString s3; s3 << MyFrame::m_match_params->lr_threshold;
    m_lrthreshTxt        = new wxTextCtrl(panel, ID_THRESH, s3);

    wxStaticText* label2 = new wxStaticText(panel, wxID_ANY, _("Relationships"));
    wxArrayString cstrings;
    cstrings.Add(wxT("IDENT"));
    cstrings.Add(wxT("D1"));
    cstrings.Add(wxT("SIB"));
    cstrings.Add(wxT("D2"));
    cstrings.Add(wxT("D3"));
//    cstrings.Add(wxT("D14"));
//    cstrings.Add(wxT("I14"));
//    cstrings.Add(wxT("D23"));
//    cstrings.Add(wxT("I23"));
//    cstrings.Add(wxT("SIBCSN"));
//    cstrings.Add(wxT("DFC"));
    m_relList = new wxCheckListBox(panel, ID_RELLIST, wxDefaultPosition,
//                wxSize(wxDefaultCoord, 50),
                wxSize(30, 30),
                cstrings, wxLB_MULTIPLE);

    // Default relationship selection
    std::vector<MatchType>::iterator it;
    for (it  = MyFrame::m_match_params->match_types.begin();
         it != MyFrame::m_match_params->match_types.end();
         ++it)
    {
    	switch (it->m_rel_type)
		{
    		case ident_t:
    		    m_relList->Check(0, true);
    		    break;
    		case degree_1_t:
    		    m_relList->Check(1, true);
    		    break;
    		case sibling_t:
    		    m_relList->Check(2, true);
    		    break;
    		case degree_2_t:
    		    m_relList->Check(3, true);
    		    break;
    		case degree_pq_t:
    		case gen_t:
    		case inv_t:
    		case none_t:
    		case unknown_t:
    		default:
				break;
		}
    }

    wxStaticText *label3 = new wxStaticText(panel, wxID_ANY, _("Frequency DB"));
    m_popTxt             = new wxTextCtrl(panel, wxID_ANY, _("(default)"));

	(void)populationData(); // check POPDATA is set and can be read
    string popdata = getStringEnv("POPDATA");

    m_popTxt->SetValue(wxString(popdata.c_str(), *wxConvCurrent));
    m_popTxt->SetInsertionPointEnd();
    m_popTxt->SetEditable(false);
    m_popTxt->SetBackgroundColour(readOnlyGrey);

    wxButton     *browseFDB = new wxButton(panel, ID_BROWSE_FDB, _("Browse..."));

    // Sizer Hierarchy
    wxBoxSizer *matchSizer  = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *matchHSizerA = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *matchHSizerB = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *matchHSizerZ = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *matchHSizerC = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *matchHSizerD = new wxBoxSizer( wxHORIZONTAL );

    matchSizer->Add(matchHSizerA, 0, wxEXPAND | wxALL, 2);
		matchHSizerA->Add(label3, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerA->Add(m_popTxt, 3, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerA->Add(browseFDB, 1, wxALIGN_CENTER | wxLEFT, 10);

    matchSizer->Add(matchHSizerB, 0, wxEXPAND | wxALL, 2);
		matchHSizerB->Add(labelA, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerB->Add(m_cdeltaTxt, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerB->AddStretchSpacer(3);

	matchSizer->Add(matchHSizerZ, 0, wxEXPAND | wxALL, 2);
		matchHSizerZ->Add(labelB, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerZ->Add(m_rdeltaTxt, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerZ->AddStretchSpacer(3);

    matchSizer->Add(matchHSizerC, 1, wxEXPAND | wxALL, 2);
		matchHSizerC->Add(label2, 1, wxALIGN_TOP | wxLEFT, 10);
		matchHSizerC->Add(m_relList, 3, wxEXPAND | wxLEFT, 10);
		matchHSizerC->AddStretchSpacer(1);

// Extra Rels
//	wxBoxSizer *matchHSizerS0 = new wxBoxSizer( wxHORIZONTAL );
//	matchSizer->Add(matchHSizerS0, 0, wxSTRETCH_NOT | wxALL, 2);
//		matchHSizerS0->Add(labelS0, 1, wxALIGN_CENTER | wxLEFT, 10);
//		matchHSizerS0->Add(m_string0Txt, 1, wxALIGN_CENTER | wxLEFT, 13);
//		matchHSizerS0->Add(m_string1Txt, 1, wxALIGN_CENTER | wxLEFT, 0);
//		matchHSizerS0->Add(m_string2Txt, 1, wxALIGN_CENTER | wxLEFT, 0);

	matchSizer->Add(matchHSizerD, 0, wxEXPAND | wxALL, 2);
		matchHSizerD->Add(label1, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerD->Add(m_lrthreshTxt, 1, wxALIGN_CENTER | wxLEFT, 10);
		matchHSizerD->AddStretchSpacer(3);

	panel->SetSizer(matchSizer);
}

void
MatchPanel::OnBrowse(wxCommandEvent& event)
{
    cout << "MatchPanel::OnBrowse(): Id = " << event.GetId() << endl;

    wxString caption = _("Select Frequency Database");
    wxString defaultDir = wxString(homeDir(), *wxConvCurrent);

    wxDirDialog dialog(this, caption, defaultDir);

    if(dialog.ShowModal() == wxID_OK)
    {
        wxString path = dialog.GetPath();
        string pop_dir(path.ToUTF8());

        PopulationData pop;

        if (pop.read(pop_dir))
        {
        	popSet(pop);
        	m_popTxt->SetValue(path);
            m_popTxt->SetInsertionPointEnd();

        	// discard anything we have read from the database
            MyFrame::refresh();
        }
        else
        {
            wxMessageBox( _("Unable to read Frequency Database\n(reverting to previous Database)"),
                          _("Select Frequency Database"),
                          wxOK | wxICON_EXCLAMATION, this);
            return;
        }
    }
}

double
MatchPanel::crimeDelta()
{
	double ret = 0;

    if (! m_cdeltaTxt->GetValue().ToDouble(&ret))
    {
        warn << startl << "Invalid value for Crime Delta (using 0)" << endl;
        ret = 0;
    }

    return ret;
}

double
MatchPanel::refDelta()
{
	double ret = 0;

    if (! m_rdeltaTxt->GetValue().ToDouble(&ret))
    {
        warn << startl << "Invalid value for Reference Delta (using 0)" << endl;
        ret = 0;
    }

    return ret;
}

double
MatchPanel::LRThreshold()
{
	double ret = 0;

    if (! m_lrthreshTxt->GetValue().ToDouble(&ret))
    {
        warn << startl << "Invalid value for LR Threshold (using 0)" << endl;
        ret = 0;
    }

    return ret;
}

bool
MatchPanel::relIsSelected(int n)
{
	return m_relList->IsChecked(n);
}

string
MatchPanel::relString(int n)
{
	std::string tmp;  // we can assign from ToUTF8() but not construct!
	if      (n==0 && m_string0Txt) tmp = m_string0Txt->GetValue().ToUTF8();
	else if (n==1 && m_string1Txt) tmp = m_string1Txt->GetValue().ToUTF8();
	else if (m_string2Txt)         tmp = m_string2Txt->GetValue().ToUTF8();
	return tmp;
}

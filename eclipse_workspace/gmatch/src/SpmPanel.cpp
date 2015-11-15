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

#include "SpmPanel.h"
#include "gmatch.h"
#include "MyFrame.h"
#include <set>

#include "fand/MessageStream.h"
INIT_MESSAGES("ResultsPanel");
#include "fand/messages.h"

using namespace std;

BEGIN_EVENT_TABLE(SpmPanel, wxPanel)
    EVT_RADIOBUTTON(ID_HW, SpmPanel::OnChange)
    EVT_RADIOBUTTON(ID_NRC4_4, SpmPanel::OnChange)
    EVT_RADIOBUTTON(ID_NRC4_10, SpmPanel::OnChange)
    EVT_RADIOBUTTON(ID_NOMUT, SpmPanel::OnChange)
    EVT_RADIOBUTTON(ID_CRIME, SpmPanel::OnChange)
    EVT_RADIOBUTTON(ID_REF, SpmPanel::OnChange)
END_EVENT_TABLE();

void
SpmPanel::setup(wxString const &title)
{
    // Widget Hierarchy

    wxPanel *panel = this;

    // Subpopulation model
    wxStaticText *spm_label = new wxStaticText(panel, wxID_ANY, _("Subpopulation Model"));
    m_hwB     = new wxRadioButton(panel, ID_HW, _("Hardy-Weinberg (no subpopulation correction)"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    m_nrc44B  = new wxRadioButton(panel, ID_NRC4_4, _("NRC II equations 4.4 for Identity; Balding-Nichols for Familial"), wxDefaultPosition, wxDefaultSize);
    m_nrc410B = new wxRadioButton(panel, ID_NRC4_10, _("NRC II equations 4.10 for Identity; Balding-Nichols for Familial"), wxDefaultPosition, wxDefaultSize);
    m_hwB->SetValue(true);

    m_thetaLbl = new wxStaticText(panel, wxID_ANY, _("Fst or theta"));
    wxString s1; s1 << MyFrame::m_match_params->spm.theta_bar;
    m_thetaTxt = new wxTextCtrl(panel, ID_THETA, s1);

    // subpopulation Default selection
    switch (MyFrame::m_match_params->spm.type)
    {
    case SubPopModel::HW:
    	m_hwB->SetValue(true);
    	m_thetaTxt->SetValue(_("0"));
    	break;

    case SubPopModel::NRC4_4:
    	m_nrc44B->SetValue(true);
    	break;

    case SubPopModel::NRC4_10:
    case SubPopModel::B11:
    	m_nrc410B->SetValue(true);
    	break;

    default:
    	Assert2(false, "Unknown subpopulation model");
    	break;
    }

    // Mutation model
    wxStaticText *mut_label = new wxStaticText(panel, wxID_ANY, _("Mutation Model (one-to-one matches only)"));

    m_nomutB = new wxRadioButton(panel, ID_NOMUT, _("No mutation model"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    m_crimeB = new wxRadioButton(panel, ID_CRIME, _("Stepwise mutation model: crime profile is ancestor"), wxDefaultPosition, wxDefaultSize);
    m_refB   = new wxRadioButton(panel, ID_REF, _("Stepwise mutation model: reference profile is ancestor"), wxDefaultPosition, wxDefaultSize);

    m_mutLbl = new wxStaticText(panel, wxID_ANY, _("Mutation rate\n(per allele, per generation)"));
    wxString s2; s1 << MyFrame::m_match_params->mut.mutation_rate;
    m_mutTxt = new wxTextCtrl(panel, ID_THETA, s1);

    // mutation Default selection
    switch (MyFrame::m_match_params->mut.type)
    {
    case MutModel::NO_MUTATIONS:
    	m_nomutB->SetValue(true);
    	m_mutTxt->SetValue(_("0"));
    	break;

    case MutModel::CRIME_IS_ANCESTOR:
    	m_crimeB->SetValue(true);
    	break;

    case MutModel::REF_IS_ANCESTOR:
    	m_refB->SetValue(true);
    	break;

    default:
    	Assert2(false, "Unknown mutation model");
    	break;
    }

    doLogic();

    // Sizer Hierarchy
    wxBoxSizer  *topSizer = new wxBoxSizer( wxVERTICAL );

    topSizer->AddSpacer(10);
    topSizer->Add(spm_label, 0, wxEXPAND | wxALL, 5);
//    topSizer->AddSpacer(10);

    topSizer->Add(m_hwB, 0, wxEXPAND | wxRIGHT | wxLEFT, 5);
    topSizer->Add(m_nrc44B, 0, wxEXPAND | wxRIGHT | wxLEFT, 5);
    topSizer->Add(m_nrc410B, 0, wxEXPAND | wxRIGHT | wxLEFT, 5);
    topSizer->AddSpacer(2);
    wxBoxSizer *hSizer = new wxBoxSizer( wxHORIZONTAL );
    topSizer->Add(hSizer, 0, wxEXPAND | wxALL, 5);
		hSizer->Add(m_thetaLbl, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL  | wxALL, 5);
		hSizer->Add(m_thetaTxt, 0, wxLEFT | wxALL, 5);
	topSizer->AddSpacer(10);

	int space = 5;
    topSizer->Add(mut_label, 0, wxEXPAND | wxALL, space);
//	topSizer->AddSpacer(10);
    topSizer->Add(m_nomutB, 0, wxEXPAND | wxRIGHT | wxLEFT, space);
    topSizer->Add(m_crimeB, 0, wxEXPAND | wxRIGHT | wxLEFT, space);
    topSizer->Add(m_refB, 0, wxEXPAND | wxRIGHT | wxLEFT, space);
    topSizer->AddSpacer(2);
    wxBoxSizer *hSizer2 = new wxBoxSizer( wxHORIZONTAL );
    topSizer->Add(hSizer2, 0, wxEXPAND | wxALL, space);
		hSizer2->Add(m_mutLbl, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL  | wxALL, space);
		hSizer2->Add(m_mutTxt, 0, wxLEFT | wxALL, space);

	panel->SetSizer(topSizer);
}

void
SpmPanel::OnChange(wxCommandEvent& event)
{
    cout << "OnChange(): Id = " << event.GetId() << endl;

    doLogic();
}


void
SpmPanel::doLogic()
{
	bool hw = m_hwB->GetValue();

    // enable/disable theta
	m_thetaLbl->Enable(!hw);
	m_thetaTxt->Enable(!hw);

	bool mut = ! m_nomutB->GetValue();

    // enable/disable mutation rate
	m_mutLbl->Enable(mut);
	m_mutTxt->Enable(mut);
}

SubPopModel
SpmPanel::subPopModel() const
{
	SubPopModel::Type type = SubPopModel::HW;

	if (m_nrc44B->GetValue())
	{
		type = SubPopModel::NRC4_4;
	}

	if (m_nrc410B->GetValue())
	{
		type = SubPopModel::B11;
	}

	double theta = 0;

    if (! m_thetaTxt->GetValue().ToDouble(&theta) ||
    		theta < 0 || theta > 1)
    {
        warn << startl << "Invalid value for theta (using 0)" << endl;
        theta = 0;
    }

	return SubPopModel(type, theta);
}

MutModel
SpmPanel::mutModel() const
{
	MutModel::Type type = MutModel::NO_MUTATIONS;

	if (m_crimeB->GetValue())
	{
		type = MutModel::CRIME_IS_ANCESTOR;
	}

	if (m_refB->GetValue())
	{
		type = MutModel::REF_IS_ANCESTOR;
	}

	double mut_rate = 0;

    if (! m_mutTxt->GetValue().ToDouble(&mut_rate) ||
    		mut_rate < 0 || mut_rate >= 1)
    {
        warn << startl << "Invalid value for mutation rate (using 0)" << endl;
        mut_rate = 0;
    }

	return MutModel(type, mut_rate);
}

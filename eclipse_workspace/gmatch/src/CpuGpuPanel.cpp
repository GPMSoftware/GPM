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

#include "CpuGpuPanel.h"
#include "gmatch.h"
#include "MyFrame.h"
#include <set>

#include "fand/MessageStream.h"
INIT_MESSAGES("ResultsPanel");
#include "fand/messages.h"

using namespace std;

//BEGIN_EVENT_TABLE(CpuGpuPanel, wxPanel)
//    EVT_RADIOBUTTON(ID_CPU, CpuGpuPanel::OnCpu)
//    EVT_RADIOBUTTON(ID_GPU, CpuGpuPanel::OnGpu)
//END_EVENT_TABLE();

void
CpuGpuPanel::setup(wxString const &title)
{
    // Widget Hierarchy

    wxPanel *panel = this;

    wxStaticText* s_label = new wxStaticText(panel, wxID_ANY, _("One to one:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
    m_s_cpuB  = new wxRadioButton(panel, ID_CPU, _("CPU"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    m_s_gpuB  = new wxRadioButton(panel, ID_GPU, _("GPU"), wxDefaultPosition, wxDefaultSize);
    m_s_gpuB->Enable(false);

    wxStaticText* n_label = new wxStaticText(panel, wxID_ANY, _("One to many:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
    m_n_cpuB  = new wxRadioButton(panel, ID_CPU, _("CPU"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    m_n_gpuB  = new wxRadioButton(panel, ID_GPU, _("GPU"), wxDefaultPosition, wxDefaultSize);

    wxStaticText* n2_label = new wxStaticText(panel, wxID_ANY, _("Many to many:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
    m_n2_cpuB = new wxRadioButton(panel, ID_CPU, _("CPU"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    m_n2_gpuB = new wxRadioButton(panel, ID_GPU, _("GPU"), wxDefaultPosition, wxDefaultSize);

    wxArrayString cstrings;

    int n_devices = gpuDevices().GPUsPresent();

    for (int i=0; i< n_devices; ++i)
    {
    	cudaDeviceProp props = gpuDevices().deviceProp(i);
    	cstrings.Add(props.name);
    	if (gpuDevices().hasDisplay(i))
    	{
    		cstrings.back() += _(" (display)");
    	}
    }

    m_cudaProcsList = new wxCheckListBox(panel, ID_GPULIST, wxDefaultPosition,
                wxSize(wxDefaultCoord, 150), cstrings, wxLB_MULTIPLE);

    // Default selection
    std::set<int> gpus;
    (void)gpuDevices().GPUsEnabled(&gpus);
    std::set<int>::iterator it;
    for (it = gpus.begin(); it != gpus.end(); ++it)
    {
        m_cudaProcsList->Check(*it, true);
    }

    if (MyFrame::m_match_params->n_cuda_accel)
    {
    	m_n_gpuB->SetValue(true);
    }

    if (MyFrame::m_match_params->n2_cuda_accel)
    {
    	m_n2_gpuB->SetValue(true);
    }

    // Sizer Hierarchy
    wxBoxSizer  *topSizer = new wxBoxSizer( wxVERTICAL );
    wxGridSizer *gSizer   = new wxGridSizer( 3, 3, 0, 0 );
    wxBoxSizer  *hSizer   = new wxBoxSizer( wxHORIZONTAL );
    m_gpuBox = new wxStaticBoxSizer( wxHORIZONTAL, panel, _("Enable CUDA GPUs") );

    topSizer->Add(gSizer, 0, wxALL, 5);
		gSizer->Add(s_label,  0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 2);
		gSizer->Add(m_s_cpuB, 0, wxEXPAND | wxALL, 2);
		gSizer->Add(m_s_gpuB, 0, wxEXPAND | wxALL, 2);

		gSizer->Add(n_label,  0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 2);
		gSizer->Add(m_n_cpuB, 0, wxEXPAND | wxALL, 2);
		gSizer->Add(m_n_gpuB, 0, wxEXPAND | wxALL, 2);

		gSizer->Add(n2_label,  0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 2);
		gSizer->Add(m_n2_cpuB, 0, wxEXPAND | wxALL, 2);
		gSizer->Add(m_n2_gpuB, 0, wxEXPAND | wxALL, 2);

    topSizer->Add(hSizer, 0, wxEXPAND | wxALL, 2);
    hSizer->AddSpacer(20);
    hSizer->Add(m_gpuBox, 0, wxEXPAND | wxALL, 2);
    m_gpuBox->Add(m_cudaProcsList, 0, wxEXPAND | wxALL, 2);

	panel->SetSizer(topSizer);
}

bool
CpuGpuPanel::nCudaAccel() const
{
	return m_n_gpuB->GetValue();
}

bool
CpuGpuPanel::n2CudaAccel() const
{
	return m_n2_gpuB->GetValue();
}

int
CpuGpuPanel::numGPUsListed() const
{
	return m_cudaProcsList->GetCount();
}

bool
CpuGpuPanel::gpuIsSelected(int n) const
{
	return m_cudaProcsList->IsChecked(n);
}


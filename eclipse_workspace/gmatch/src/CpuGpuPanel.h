/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * CpuGpuPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef CPUGPUPANEL_H_
#define CPUGPUPANEL_H_

#include "wx/wx.h"

class CpuGpuPanel : public wxPanel
{
public:
    // Constructor
	CpuGpuPanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
    , m_s_cpuB(0)
    , m_s_gpuB(0)
    , m_n_cpuB(0)
    , m_n_gpuB(0)
    , m_n2_cpuB(0)
    , m_n2_gpuB(0)
	, m_gpuBox(0)
    , m_cudaProcsList(0)
    {
        setup(title);
    }

	wxRadioButton         *m_s_cpuB;
	wxRadioButton         *m_s_gpuB;
	wxRadioButton         *m_n_cpuB;
	wxRadioButton         *m_n_gpuB;
	wxRadioButton         *m_n2_cpuB;
	wxRadioButton         *m_n2_gpuB;
	wxStaticBoxSizer      *m_gpuBox;
    wxCheckListBox        *m_cudaProcsList;

    bool nCudaAccel() const;
    bool n2CudaAccel() const;
    int  numGPUsListed() const;
    bool gpuIsSelected(int n) const;

//    DECLARE_EVENT_TABLE();

private:
    void setup(wxString const &title);

};

#endif /* CPUGPUPANEL_H_ */

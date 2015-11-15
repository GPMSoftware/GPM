/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * SpmPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef SPMPANEL_H_
#define SPMPANEL_H_

#include "fand/SubPopModel.h"
#include "fand/MutModel.h"

#include "wx/wx.h"

class SpmPanel : public wxPanel
{
public:
    // Constructor
	SpmPanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
	, m_hwB(0)
	, m_nrc44B(0)
	, m_nrc410B(0)
	, m_nomutB(0)
	, m_crimeB(0)
	, m_refB(0)
	, m_thetaLbl(0)
	, m_thetaTxt(0)
	, m_mutLbl(0)
	, m_mutTxt(0)
    {
        setup(title);
    }

	wxRadioButton         *m_hwB;
	wxRadioButton         *m_nrc44B;
	wxRadioButton         *m_nrc410B;

	wxRadioButton		  *m_nomutB;
	wxRadioButton         *m_crimeB;
	wxRadioButton         *m_refB;

	wxStaticText          *m_thetaLbl;
	wxTextCtrl            *m_thetaTxt;
	wxStaticText          *m_mutLbl;
	wxTextCtrl            *m_mutTxt;


	SubPopModel subPopModel() const;
	MutModel mutModel() const;

    DECLARE_EVENT_TABLE();

private:
    void setup(wxString const &title);
	void OnChange(wxCommandEvent& event);
	void doLogic();
};

#endif /*SPMPANEL_H_ */

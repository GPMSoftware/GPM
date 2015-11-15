/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MegaChoiceDialog.h
 *
 *  Created on: Oct 18, 2011
 *      Author: gareth
 *
 *      A dialog for choosing between a large number of strings.
 *      I.e. more than it is sensible to put in a SingleChoiceDialog
 */

#ifndef MEGACHOICEDIALOG_H_
#define MEGACHOICEDIALOG_H_

#include "fand/Database.h"

#include "wx/wx.h"

class MegaChoiceDialog : public wxDialog
{
public:
	enum Mode
	{
		SINGLE,
		MULTI,
	};

    // Constructor
	MegaChoiceDialog(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxString& caption = _("Metadata"),
            std::set<std::string> const &vals = std::set<std::string>(),
            const Mode &mode = SINGLE,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxCAPTION | wxRESIZE_BORDER | wxSYSTEM_MENU )
    : wxDialog(parent, winid, caption, pos, size, style)
	, m_vals(vals)
	, m_mode(mode)
    , m_text(0)
    , m_list(0)
    , m_prevBtn(0)
    , m_nextBtn(0)
    {
        setup();
    }

    wxString GetStringSelection() const { return m_list->GetStringSelection(); }
    wxArrayString GetStringSelections();


    DECLARE_EVENT_TABLE()

private:
    void setup();
    void list(std::string s, int n);
    void list(int n);
    void OnSearch(wxCommandEvent& event);
    void OnPrev(wxCommandEvent& event);
    void OnNext(wxCommandEvent& event);
    void OnAdd(wxCommandEvent& event);
    void OnRemove(wxCommandEvent& event);
    void updatePrevNextButtons();

    std::set<std::string> const &m_vals;                  // values to choose from
    const Mode &m_mode;                                   //

    std::set<std::string>::const_iterator m_begin, m_end; // displayed values range

    wxTextCtrl *m_text;
    wxListBox  *m_list;
    wxButton   *m_prevBtn;
    wxButton   *m_nextBtn;
    wxButton   *m_addBtn;
    wxButton   *m_remBtn;
    wxListBox  *m_list2;

};

#endif /* MEGACHOICEDIALOG_H_ */

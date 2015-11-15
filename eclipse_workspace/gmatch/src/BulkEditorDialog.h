/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * BulkEditorDialog.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef BULKEDITORDIALOG_H_
#define BULKEDITORDIALOG_H_

#include "fand/Database.h"
#include "SelectProfilesPanel.h"

#include "wx/wx.h"

class BulkEditorDialog : public wxDialog
{
public:
    // Constructor
    BulkEditorDialog(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxString& caption = _("BulkEditorDialog"),
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxCAPTION | wxRESIZE_BORDER | wxSYSTEM_MENU )
    : wxDialog(parent, winid, caption, pos, size, style)
    {
        setup();
    }

    static void popup();

private:
    SelectProfilesPanel* m_sp_panel;
    wxChoice   *m_action;
    wxTextCtrl *m_value;
    wxButton   *m_apply;

    DECLARE_EVENT_TABLE()
//    DECLARE_CLASS()

private:
    void OnSize(wxSizeEvent& event);

    void setup();
    void OnApply(wxCommandEvent& event);
    void OnChange(wxCommandEvent& event);
    void apply();
    bool doDelete();
    bool doMove();
    void enableEdits(bool enabled);

    static BulkEditorDialog *theBulkEditorDialog;

    static void profilesFound();

};

#endif /* BULKEDITORDIALOG_H_ */

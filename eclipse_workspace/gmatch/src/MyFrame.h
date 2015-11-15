/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MyFrame.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef MYFRAME_H_
#define MYFRAME_H_

#include "InputPanel.h"
#include "SelectProfilesPanel.h"
#include "MatchPanel.h"
#include "SpmPanel.h"
#include "CpuGpuPanel.h"
#include "DatabasePanel.h"
#include "ResultsPanel.h"
#include "MetadataDialog.h"

#include "fand/dbmatch.h"

#include "wx/wx.h"

class MyFrame: public wxFrame
{
public:

    MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
    void OnQuit(wxCommandEvent& event);
    void OnAction(wxCommandEvent& event);
    void OnMatch(wxCommandEvent& event);
    void OnSelfOrRef(wxCommandEvent& event);
//    void OnSashDrag(wxSashEvent& event);

    InputPanel *createInputArea(wxWindow *parent);
    SelectProfilesPanel *createSelectArea(wxWindow *parent, wxString const &title);
    MatchPanel   *createMatchArea(wxWindow *parent, wxString const &title);
    SpmPanel     *createSpmArea(wxWindow *parent, wxString const &title);
    CpuGpuPanel  *createCpuGpuArea(wxWindow *parent, wxString const &title);
    DatabasePanel  *createDatabaseArea(wxWindow *parent, wxString const &title);
    wxTextCtrl   *createTextArea(wxWindow *parent);
    ResultsPanel *createResultsArea(wxWindow *parent, wxString const &title);

    static SelectProfilesPanel *m_cSetSelect;
    static SelectProfilesPanel *m_rSetSelect;
    static MatchPanel          *m_matchp;
    static SpmPanel            *m_subpopmodel;
    static CpuGpuPanel         *m_processors;
    static DatabasePanel       *m_database;
    static ResultsPanel        *m_resultsArea;

    static InputPanel     *m_inputArea;
    static wxTextCtrl     *m_messages;

    static wxStaticText   *m_Clabel;
    static wxStaticText   *m_Rlabel;
    static wxRadioBox     *m_radioB;
    static wxButton       *m_match;

    static NMmatchResults m_match_results;
    static MatchParams   *m_match_params;

    static void refreshMessageWindow();
    static void refreshMatchButton();
    static bool displayProfile(DBProfile const &p);

    static void refresh()
    {
        m_cSetSelect->refresh();
        m_rSetSelect->refresh();
        m_inputArea->refresh();
        refreshMessageWindow();
    }

    static void importFile(wxWindow *parent);
    static void exportFile(wxWindow *parent);

    static ProfileRange & cSet();
    static ProfileRange & rSet();

    static std::string cQuery();
    static std::string rQuery();

    DECLARE_EVENT_TABLE();

private:
    static bool m_constructed;

    static bool refIsCrime();
};

class ImportDialog : public wxDialog
{
public:
    // Constructor
    ImportDialog(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxString& caption = _("Import profiles"),
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxCAPTION | wxRESIZE_BORDER | wxSYSTEM_MENU )
    : wxDialog(parent, winid, caption, pos, size, style)
    , m_file(0)
    , m_dataset(0)
    , m_key1(0)
    , m_key2(0)
    {
        setup();
    }

    wxTextCtrl *m_file;
    wxComboBox *m_dataset;
    wxTextCtrl *m_key1;
    wxTextCtrl *m_key2;

    void OnBrowse(wxCommandEvent& event);
    void OnChange(wxCommandEvent& event);
    int keyCol1();
    int keyCol2();

    DECLARE_EVENT_TABLE();

private:
    void setup();
};

#endif /* MYFRAME_H_ */

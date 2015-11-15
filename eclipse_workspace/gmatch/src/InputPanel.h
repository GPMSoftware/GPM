/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * InputPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef INPUTPANEL_H_
#define INPUTPANEL_H_

#include "fand/Database.h"
#include "MyGrid.h"

#include "wx/wx.h"
#include "wx/grid.h"

#include <vector>
#include <map>
#include <string>

class InputGrid : public MyGrid
{
public:

	InputGrid( wxWindow *parent,
        wxWindowID id,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxWANTS_CHARS,
        const wxString& name = wxGridNameStr )
    : MyGrid(parent, id, pos, size, style, name)
    {
		setup();
    }

	void setup();
    void OnChar(wxKeyEvent& event);

    // convert locus to grid position
    int gPos(int i) { return m_gridLoc[i]; }

    // map or re-map the loci onto the grid
    void mapLoci(Locus *kit_loci, int kit_len);

    void OnErase(wxEraseEvent& event);

private:
    std::vector<int> m_gridLoc; // m_gridLoc[loc] = grid position of loc

    DECLARE_EVENT_TABLE()

};

class InputPanel : public wxPanel
{
public:
    // Constructor
    InputPanel(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(parent, winid, pos, size, style, name)
    , m_error(false)
    , m_unsaved_edits(false)
    , m_accept_select(true)
    , m_grid(0)
    , m_dataset(0)
    , m_meta(0)
    , m_edit(0)
    , m_del(0)
    , m_save(0)
    , m_mode(MODE_EMPTY)
    , m_num_profiles(0)
    {
        setup();
    }

    bool emptyMode();
    void clearView();

    bool displayProfile(DBProfile const &p);
    bool displayProfiles(std::vector<DBProfile> const &dbp_vec);
    void refresh();
    void viewProfiles(std::vector<std::string> const &p_keys);

    void selectProfiles();
    void newSample();
    void editProfile();
    void deleteProfile();
    void clearProfile();
    void saveProfile();

private:
    void setup();

    void OnMeta(wxCommandEvent& event);
    void OnBrowse(wxCommandEvent& event);
    void OnNew(wxCommandEvent& event);
    void OnEdit(wxCommandEvent& event);
    void OnDelete(wxCommandEvent& event);
    void OnClear(wxCommandEvent& event);
    void OnSave(wxCommandEvent& event);
    void OnEditCell(wxGridEvent& event);
    void OnSelectKit(wxCommandEvent& event);
    void OnSelectCell(wxGridEvent& event);

    void setReadOnly(bool ro);
    void setReadOnly(int i, int j, bool ro);

    void newSampleMode();
    void editMode();
    bool viewMode(std::vector<DBProfile> const &dbp_vec);
    bool viewMode(DBProfile const &dbp, bool is_first=true);
    bool tryClear();

    std::string sampleID();
    std::string profileID(int prof);
    std::string dataset();
    int numContribs(int prof);

    int profFirstRow(int nprof);
    int profLastRow(int nprof);
    std::pair<int, int> profAndAlleleAtRow(int row);
    void changeNumContribs(int nprof, int ncontribs);
    void newProfile();
    void enableCellsForEdit(bool new_sample);
    void setNumberRows(int rows);
    void setCellError(int i, int j);

    void saveSample();
    bool checkData(std::string allele_string);
    void transformAlleleData(std::map<std::pair<int, int>, std::string > &allele_data);
    void doDeleteProfile();
    void enableOps(bool enable);
    void enableSave(bool enable);

    void setAlleleTextByGrid(int allele, int g, std::string val);
    void setAlleleTextByLocus(int allele, int loc, std::string val);
    std::string getAlleleTextByGrid(int allele, int g);
    std::string getAlleleTextByLocus(int allele, int loc);


    bool m_error;
    bool m_unsaved_edits;   // unused
    bool m_accept_select;

    InputGrid  *m_grid;    // Input/edit grid
    wxComboBox *m_dataset;
    wxChoice   *m_profileset;
    wxButton   *m_meta;
    wxButton   *m_edit;
    wxButton   *m_clear;
    wxButton   *m_del;
    wxButton   *m_save;
    wxButton   *m_browse;
    wxChoice   *m_kit;

    enum Mode
    {
        MODE_EMPTY = 0,
        MODE_VIEW,
        MODE_EDIT,
        MODE_NEW,
    };

    Mode m_mode;
    int m_num_profiles;             // >1 only in MODE_NEW
    std::vector<int> m_num_contribs;     // for each profile

    MetaData m_metadata; // the metadata of the first (or only) profile we are viewing

    DECLARE_EVENT_TABLE()

};


#endif /* INPUTPANEL_H_ */

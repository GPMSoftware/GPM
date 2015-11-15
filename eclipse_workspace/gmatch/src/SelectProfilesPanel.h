/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * SelectProfilesPanel.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef SELECTPROFILESPANEL_H_
#define SELECTPROFILESPANEL_H_

#include "fand/Profile.h"
#include "fand/Database.h"
#include "fand/geotrans.h"

#include "wx/wx.h"

#include <vector>
#include <set>
#include <string>

class MultiSelector;

class SelectProfilesPanel : public wxPanel
{
public:
    // Constructor
    SelectProfilesPanel(
            wxWindow *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr);

    ~SelectProfilesPanel();

    void OnFind(wxCommandEvent& event);
    void OnSave(wxCommandEvent& event);
    void OnDelete(wxCommandEvent& event);

	void refresh();
	int selected() { return m_selected_profiles.size(); }
    void hasChanged();

	std::string sqlQuery();       // SQL query to get the selected profiles
	std::string sqlWhereClause(); // just the WHERE clause

	void setFindNotify(void (*f)(void) ) { findNotify = f; } // set find callback function

    MultiSelector *createMultiSelector(SelectProfilesPanel *parent, wxString const &title);

    wxComboBox *m_saved;
    wxButton *m_saveBtn, *m_delBtn;

    MultiSelector *m_selectp;

    wxButton         *m_findBtn;   // find button
    wxStaticText     *m_pcountLbl; // label showing count of profiles selected

    wxTextCtrl *m_km;
    wxTextCtrl *m_mgrs;

    ProfileRange m_selected_profiles;
    std::set<std::string> m_profile_ids;              // all the IDs in the database

    enum { NPROFS = 10 }; // number of profiles to display in combo box

    static std::map<std::string, std::string> m_searches; // Named searches

    static std::set<SelectProfilesPanel*> m_instances;

    DECLARE_EVENT_TABLE()

private:
    void setup(wxString const &title);

    void loadSearchesFromFile();
    static void saveSearchesToFile();
    void refreshSearches();
    static void updateSearches();
    void OnChangeValue(wxCommandEvent& event);
    void OnChangeSearch(wxCommandEvent& event);
    void recallSearch();
	std::string distQuery();
    LatLonStrings latLonStrings();

    void (*findNotify)(void);

	friend std::istream & operator>>(std::istream &is, SelectProfilesPanel &a);
};

// Many selectors
class MultiSelector : public wxPanel
{
public:
	MultiSelector(
            SelectProfilesPanel *parent,
            wxString const &title,
            wxWindowID winid = wxID_ANY,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxTAB_TRAVERSAL | wxNO_BORDER,
            const wxString& name = wxPanelNameStr)
    : wxPanel(dynamic_cast<wxWindow*>(parent), winid, pos, size, style, name)
	, m_profiles_panel(parent)
	, m_scrolled(0)
	, m_selsizer(0)
    {
        setup(title);
    }

	std::string sqlQuery();
	void update();

private:

    wxArrayString fstrings();
    wxArrayString vstrings();

    // A single line of controls in the SelectProfilesPanel, allowing
    // a single criterion for the selection of profiles to be specified
	struct Selector : public wxPanel
	{
		Selector(MultiSelector *parent);
		void setup();
		MultiSelector *multisel; // parent
		wxBoxSizer    *sizer;
		wxComboBox    *field;
		wxComboBox    *value;    // COMBO_BOX mode
		wxTextCtrl    *text;     // SINGLE_CHOICE mode
		wxButton      *browse;   // SINGLE_CHOICE mode
		wxButton      *newBtn;
		wxButton      *delBtn;
		wxButton      *logicBtn;

		bool logic;

		enum
		{
			COMBO_BOX,
			SINGLE_CHOICE,
			TEXT,
			MEGA_CHOICE,
		} mode;

		int id;
		static int next_id;
		std::set<std::string> vals; // possible values of the field

		std::string getLogic() const;
		std::string getName() const;
		std::string getValue() const;
		std::string getStreamedValue() const;

		void setLogic(std::string logic_str);
		void setName(std::string const &name);
		void setValue(std::string const &value);
		void setStreamedValue(std::string const &val_str);
		void setLogicBtn(bool logic);

		void OnNew(wxCommandEvent& event);
		void OnDel(wxCommandEvent& event);
		void OnLogic(wxCommandEvent& event);
		void OnChangeField(wxCommandEvent& event);
		void OnChangeValue(wxCommandEvent& event);
		void OnBrowse(wxCommandEvent& event);

		void update();
		void updateValues(DBField const *dbf);

		DECLARE_EVENT_TABLE();

	private:
		const DBField * getDbf();
		void megaChoice();
		void singleChoice();
	};

    void setup(wxString const &title);
    int  findSelector(int id);
	void addSelectorPos(int pos);
	void addSelectorId(int id);
	void delSelectorPos(int pos);
	void delSelectorId(int id);
	void hasChanged() { m_profiles_panel->hasChanged();	}
	void clear();

	SelectProfilesPanel   *m_profiles_panel;
	wxScrolledWindow      *m_scrolled;
	wxBoxSizer            *m_selsizer;

	std::vector<Selector*> m_selectors;

	struct SelectorField // a field on which you can search
	{
		SelectorField(std::string const &n, std::string const &f, std::string const &t = "")
		: name(n), db_field(f, t) {}
		std::string name;
		DBField     db_field;
	};
	std::vector<SelectorField> m_sel_fields;

	const DBField* findDBFieldFromName(const std::string& name);

	friend std::ostream & operator<<(std::ostream &os, const MultiSelector &a);
	friend std::istream & operator>>(std::istream &is, MultiSelector &a);

	friend class SelectProfilesPanel;
    friend class TestStream_MultiSelector_SelectorField; // unit tests
};

std::ostream &
operator<<(std::ostream &os, const MultiSelector &a);

std::istream &
operator>>(std::istream &is, MultiSelector &a);

std::ostream &
operator<<(std::ostream &os, const SelectProfilesPanel &a);

std::istream &
operator>>(std::istream &is, SelectProfilesPanel &a);

#endif /* SELECTPROFILESPANEL_H_ */

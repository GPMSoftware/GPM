/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MyFrame.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#include "MyFrame.h"
#include "gmatch.h"
#include "ResultsPanel.h"
#include "BulkEditorDialog.h"

#include "fand/Version.h"
#include "fand/HostManager.h"

#include "fand/MessageStream.h"
INIT_MESSAGES("MyFrame");
#include "fand/messages.h"

#include "wx/statline.h"
#include "wx/notebook.h"
#include "wx/imaglist.h"

#include <fenv.h>

using namespace std;

wxTextCtrl          * MyFrame::m_messages = 0;
InputPanel          * MyFrame::m_inputArea = 0;
SelectProfilesPanel * MyFrame::m_cSetSelect = 0;
SelectProfilesPanel * MyFrame::m_rSetSelect = 0;
MatchPanel          * MyFrame::m_matchp = 0;
SpmPanel            * MyFrame::m_subpopmodel = 0;
CpuGpuPanel         * MyFrame::m_processors = 0;
DatabasePanel       * MyFrame::m_database = 0;
ResultsPanel        * MyFrame::m_resultsArea = 0;
wxRadioBox          * MyFrame::m_radioB = 0;
wxButton            * MyFrame::m_match = 0;
wxStaticText        * MyFrame::m_Clabel = 0;
wxStaticText        * MyFrame::m_Rlabel = 0;

bool                  MyFrame::m_constructed(false);
NMmatchResults        MyFrame::m_match_results(0);
MatchParams         * MyFrame::m_match_params = 0;

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(ID_QUIT, MyFrame::OnQuit)
    EVT_MENU(ID_IMPORT, MyFrame::OnAction)
    EVT_MENU(ID_EXPORT, MyFrame::OnAction)
    EVT_MENU(ID_MENU_NEW, MyFrame::OnAction)
    EVT_MENU(ID_MENU_EDIT, MyFrame::OnAction)
    EVT_MENU(ID_MENU_DELETE, MyFrame::OnAction)
    EVT_MENU(ID_MENU_CLEAR, MyFrame::OnAction)
    EVT_MENU(ID_MENU_SAVE, MyFrame::OnAction)
    EVT_MENU(ID_MENU_BULK, MyFrame::OnAction)
    EVT_RADIOBOX(ID_CR, MyFrame::OnSelfOrRef)
    EVT_BUTTON(ID_MATCH, MyFrame::OnMatch)

//    EVT_SIZE(MyFrame::OnSize)
//    EVT_SASH_DRAGGED(ID_INPUT_AREA, MyFrame::OnSashDrag)
END_EVENT_TABLE();

static void profilesFound()
{
	MyFrame::refreshMessageWindow();
	MyFrame::refreshMatchButton();
}

MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
: wxFrame ( NULL, -1, title, pos, size )
{
    // enable floating point exceptions
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // print program name, version and date
    cout << PROGRAM_NAME << " (gmatch): Version " << MAJOR_VERSION << "." << MINOR_VERSION;
    time_t t = time(0);
    cout << " - Starting at " << ctime(&t) << endl;

    // Read environment variables
    Env env;
    env.read();
    cout << env << endl;

    // create default match parameters from env
    m_match_params = new MatchParams(env);

    wxMenu *fileMenu = new wxMenu;
    fileMenu->Append( ID_IMPORT, _("&Import profiles...") );
    fileMenu->Append( ID_EXPORT, _("&Export database...") );
    fileMenu->AppendSeparator();
    fileMenu->Append( ID_QUIT, _("E&xit") );

    wxMenu *editMenu = new wxMenu;
    // TODO select profiles from SUBMENU all profiles, cset, rset
    editMenu->Append( ID_MENU_NEW,    _("&New sample...") );
    editMenu->Append( ID_MENU_EDIT,   _("&Edit profile") );
    editMenu->Append( ID_MENU_DELETE, _("&Delete profile") );
    editMenu->Append( ID_MENU_CLEAR,  _("&Clear profile") );
    editMenu->Append( ID_MENU_SAVE,   _("&Save profile") );
    editMenu->Append( ID_MENU_BULK,   _("&Bulk editor") );

    wxMenuBar *menuBar = new wxMenuBar;
    menuBar->Append( fileMenu, _("&File") );
    menuBar->Append( editMenu, _("&Edit") );

    SetMenuBar( menuBar );

    CreateStatusBar();
    SetStatusText( _("Ready") );

    // Widget Hierarchy
    wxPanel* panel = new wxPanel(this, wxID_ANY, wxDefaultPosition, wxSize(100,100));
    wxFont font(
    		8,                  // size
    		wxFONTFAMILY_SWISS, // family
    		wxNORMAL,           // style
    		wxNORMAL,           // weight
    		false,              // underlined
    		_("Sans")           // face
    		);

    panel->SetFont(font);

    panel->Show(true);

    wxNotebook *tabs1 = new wxNotebook(panel, wxID_ANY);

    m_database    = createDatabaseArea(tabs1, _("Database")); // create first because opens db connection
    m_cSetSelect  = createSelectArea(tabs1, _("C-Set"));
    m_cSetSelect->setFindNotify(profilesFound);
    m_rSetSelect  = createSelectArea(tabs1, _("R-Set"));
    m_rSetSelect->setFindNotify(profilesFound);
    m_matchp      = createMatchArea(tabs1, _("Match"));
    m_subpopmodel = createSpmArea(tabs1, _("SPMs"));
    m_processors  = createCpuGpuArea(tabs1, _("Processors"));

    m_inputArea = createInputArea(panel);
    m_inputArea->Show(true);

    tabs1->AddPage(m_cSetSelect, wxT("C-Set"), true);
    tabs1->AddPage(m_rSetSelect, wxT("R-Set"), false);
    tabs1->AddPage(m_matchp, wxT("Match parameters"), false);
    tabs1->AddPage(m_subpopmodel, wxT("Subpopulation/Mutation"), false);
    tabs1->AddPage(m_processors, wxT("Processors"), false);
    tabs1->AddPage(m_database, wxT("Database"), false);

    m_resultsArea = createResultsArea(panel, _("Results"));
    m_messages = createTextArea(panel);

    m_Clabel = new wxStaticText(panel, wxID_ANY, _("0 Selected")); // space for big numbers
    m_Rlabel = new wxStaticText(panel, wxID_ANY, _("0 Selected"));

    wxArrayString bstrings;
    bstrings.Add(wxT("Itself"));
    bstrings.Add(wxT("R-Set"));

    m_radioB = new wxRadioBox(panel, ID_CR, _("Match C-Set against"), wxDefaultPosition, wxDefaultSize, bstrings, wxRA_SPECIFY_COLS, wxRA_HORIZONTAL);
    m_radioB->SetSelection(1);

    m_match = new wxButton(panel, ID_MATCH, _("Match"));
    m_match->Enable(false);

    // Sizer Hierarchy
    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *matchSizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *matchVSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *matchHSizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *cbox        = new wxStaticBoxSizer( wxHORIZONTAL, panel, _("C Profiles") );
    wxBoxSizer *rbox        = new wxStaticBoxSizer( wxHORIZONTAL, panel, _("R Profiles") );

    topSizer->Add(m_inputArea, 1, wxEXPAND | wxALL/*, 10*/);
    topSizer->Add(matchSizer, 0, wxEXPAND | wxALL/*, 10*/);
        matchSizer->Add(matchVSizer, 0, wxEXPAND | wxALL/*, 10*/);
			matchVSizer->Add(tabs1, 0, wxEXPAND | wxALL/*, 10*/);
			matchVSizer->Add(matchHSizer, 0, wxEXPAND | wxALL, 10);
				matchHSizer->Add(cbox, 1, wxALIGN_CENTER | wxRIGHT, 10);
					cbox->Add(m_Clabel, 0, wxALIGN_CENTER | wxALL, 4);
				matchHSizer->Add(rbox, 1, wxALIGN_CENTER | wxRIGHT, 10);
					rbox->Add(m_Rlabel, 0, wxALIGN_CENTER | wxALL, 4);
				matchHSizer->Add(m_radioB, 1, wxALIGN_CENTER);
				matchHSizer->Add(m_match,  1, wxALIGN_CENTER | wxALL, 4);
        matchSizer->Add(m_resultsArea, 2, wxALIGN_CENTER | wxALL/*, 10*/);
    topSizer->Add(m_messages, 0, wxEXPAND | wxALL/*, 10*/);

    topSizer->Fit(this);           // this!
    topSizer->SetSizeHints(this);  // this!

    panel->SetSizer(topSizer);     // panel!
    panel->Layout();

    m_constructed = true;

    // make sure messages generated at static initialization time are displayed
    refreshMessageWindow();

    // We have now grabbed the memory we need to run the GUI.
    // calculate memory limits for solutions
    theHostManager().init();
}

InputPanel *
MyFrame::createInputArea(wxWindow *parent)
{
    InputPanel* panel = new InputPanel(parent, ID_GRID, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}


SelectProfilesPanel *
MyFrame::createSelectArea(wxWindow *parent, wxString const &title)
{
    SelectProfilesPanel* panel = new SelectProfilesPanel(parent, title, wxID_ANY, wxDefaultPosition, wxSize(100,200), wxSIMPLE_BORDER);
    return panel;
}

MatchPanel *
MyFrame::createMatchArea(wxWindow *parent, wxString const &title)
{
    MatchPanel* panel = new MatchPanel(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

SpmPanel *
MyFrame::createSpmArea(wxWindow *parent, wxString const &title)
{
	SpmPanel* panel = new SpmPanel(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

CpuGpuPanel *
MyFrame::createCpuGpuArea(wxWindow *parent, wxString const &title)
{
	CpuGpuPanel* panel = new CpuGpuPanel(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

DatabasePanel *
MyFrame::createDatabaseArea(wxWindow *parent, wxString const &title)
{
	DatabasePanel* panel = new DatabasePanel(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

ResultsPanel *
MyFrame::createResultsArea(wxWindow *parent, wxString const &title)
{
    // Widget Hierarchy
    ResultsPanel* panel = new ResultsPanel(parent, title, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER);
    return panel;
}

wxTextCtrl *
MyFrame::createTextArea(wxWindow *parent)
{
    wxTextCtrl *text = new wxTextCtrl(parent, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize(100,100),
                                      wxTE_MULTILINE | wxTE_READONLY);

    wxFont font(
    		8,                  // size
    		wxFONTFAMILY_SWISS, // family
    		wxNORMAL,           // style
    		wxNORMAL,           // weight
    		false,              // underlined
    		_("Monospace")      // face
    		);

    text->SetFont(font);

    return text;
}

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	MemMeter m;
    cout << "OnQuit: this process is now using " << m << endl;

    Close(TRUE);
}

void MyFrame::OnAction(wxCommandEvent& event)
{
    cout << "OnAction(): id = " << event.GetId() << endl;

    const wxChar *message = 0;
    switch (event.GetId())
    {
    case ID_IMPORT:
        message = _( "Import profiles" );
        importFile(this);
        break;
    case ID_EXPORT:
        message = _( "Export database" );
        exportFile(this);
        break;

	// TODO select profiles

    case ID_MENU_NEW:
    	m_inputArea->newSample();
    	break;
    case ID_MENU_EDIT:
    	m_inputArea->clearProfile();
    	break;
    case ID_MENU_DELETE:
    	m_inputArea->deleteProfile();
    	break;
    case ID_MENU_CLEAR:
    	m_inputArea->clearProfile();
    	break;
    case ID_MENU_SAVE:
    	m_inputArea->saveProfile();
    	break;
    case ID_MENU_BULK:
    	BulkEditorDialog::popup();
    	break;
    default:
        error << startl << "unknown state in OnAction()" << endl;
        return;
    }

    refreshMessageWindow();
}

void
MyFrame::refreshMessageWindow()
{
	if (m_messages)
	{
		wxString str(MessageStream::read().c_str(), *wxConvCurrent);
		m_messages->AppendText(str);
	}
}

bool
MyFrame::refIsCrime()
{
    return (m_radioB->GetSelection() == 0);
}

void
MyFrame::refreshMatchButton()
{
	if (!m_constructed) return;

    if (m_cSetSelect->m_selected_profiles.size() &&
            ( refIsCrime() ||
              m_rSetSelect->m_selected_profiles.size()))
    {
        m_match->Enable(true);
    }
    else
    {
        m_match->Enable(false);
    }

    size_t len = 32;
    wxChar buf[len];

    int c_profiles_found = m_cSetSelect->selected();
    wxSnprintf(buf, len, _("%d Selected"), c_profiles_found);
    m_Clabel->SetLabel(buf);

    int r_profiles_found = m_rSetSelect->selected();
    wxSnprintf(buf, len, _("%d Selected"), r_profiles_found);
    m_Rlabel->SetLabel(buf);
}

bool
MyFrame::displayProfile(DBProfile const &p)
{
    // display profile and put input area into view mode
    return m_inputArea->displayProfile(p);
}

int
ImportDialog::keyCol1()
{
	long ret = 0;

    if (! m_key1->GetValue().ToLong(&ret))
    {
        warn << startl << "Invalid value for keyCol1 (using 1)" << endl;
        ret = 1;
    }

    return ret - 1;
}

int
ImportDialog::keyCol2()
{
	long ret = 0;

    if (!m_key2->IsEmpty() && !m_key2->GetValue().ToLong(&ret))
    {
        warn << startl << "Invalid value for keyCol2 (using none)" << endl;
        ret = 0;
    }

    return ret - 1;
}

void
MyFrame::importFile(wxWindow *parent)
{
    ImportDialog dialog(parent);

    if(dialog.ShowModal() == wxID_OK)
    {
        wxString path = dialog.m_file->GetValue();
        string filename(path.ToUTF8());
        ProfileFilter::FileType ftype = fileType(filename);

        string dataset = baseName(filename); // default if not given

        wxString ds = dialog.m_dataset->GetValue();
        if (! ds.empty())
        {
            dataset = ds.ToUTF8();
        }

        switch(ftype)
        {
            case ProfileFilter::B11:
            case ProfileFilter::Identifiler:
            case ProfileFilter::PowerPlex16:
            case ProfileFilter::other_csv:

                break; // OK

            case ProfileFilter::SGM_Plus:
            case ProfileFilter::type_unknown:
            default:
                error << startl << "unrecognised file type: " << filename << endl;
                error << alignl << "(should be *.b11.csv *.idf.csv *.p16.csv or *.csv)" << endl;
                return;
        }

        // open file
        std::ifstream ifs(filename.c_str());
        if (! ifs.good()){
            error << startl << "error opening file: " << filename << endl;
            return;
        }

        // read key columns
        int col1 = dialog.keyCol1();
        int col2 = dialog.keyCol2();

        // TODO: get these values from somewhere
        bool overwrite         = true;
        string kit             = "unknown";
        string evtype          = "unknown";
        float  delta           = 0;

        ProfileType prof_type = kitID(kit);
        EvidenceType evidence_type = evidenceID(evtype);

        ok << "Importing profiles from " << filename << "..." << endl;

        // import file into database
        bool op_ok = false;
        try
        {
        	op_ok = db->import(
                    overwrite,
                    ifs,
                    ftype,
                    prof_type,
                    dataset,
                    evidence_type,
                    delta,
                    col1,
                    col2);
        }
        catch(...)
        {
        	op_ok = false;
        }

        if (op_ok)
        {
            ok << "Done." << endl;
        }
        else
        {
            error << "Errors importing profiles (see log file)." << endl;
        }

        // refresh dialogs
        refresh();

        m_inputArea->clearView(); //  displayed profile may have changed
    }
}

void MyFrame::exportFile(wxWindow *parent)
{
    wxString caption = _("Export profiles");
    wxString widlcard = _("B11 (.b11.csv)|*.b11.csv");
    wxString defaultDir = wxString(homeDir(), *wxConvCurrent);
    wxString defaultFilename = wxEmptyString;

    wxFileDialog dialog(parent, caption, defaultDir, defaultFilename, widlcard, wxFD_SAVE);
    if(dialog.ShowModal() == wxID_OK)
    {
        wxString path = dialog.GetPath();
        string filename(path.ToUTF8());

        // open file
         std::ofstream ofs(filename.c_str());
         if (! ofs.good()){
             error << startl << "error opening file: " << filename << endl;
             return;
         }

        ok << "Exporting database to " << filename << "..." << endl;

        // export database to file
        if (db->exportDB(ofs))
        {
            ok << "Done." << endl;
        }
        else
        {
            error << "Errors exporting database (see log file)." << endl;
        }
    }
}

ProfileRange &
MyFrame::cSet()
{
    return m_cSetSelect->m_selected_profiles;
}

ProfileRange &
MyFrame::rSet()
{
    if (refIsCrime())
    {
        return m_cSetSelect->m_selected_profiles;
    }
    else
    {
        return m_rSetSelect->m_selected_profiles;
    }
}

string
MyFrame::cQuery()
{
    return m_cSetSelect->sqlQuery();
}

string
MyFrame::rQuery()
{
    if (refIsCrime())
    {
    	return cQuery();
    }
    else
    {
    	return m_rSetSelect->sqlQuery();
    }
}

void
MyFrame::OnSelfOrRef(wxCommandEvent& event)
{
    cout << "OnSelfOrRef(): Id = " << event.GetId() << endl;

    if (refIsCrime())
    {
        // self - disable the reference set
    	m_rSetSelect->Enable(false);
    }
    else
    {
        // reference - enable the reference set
    	m_rSetSelect->Enable(true);
    }

    refreshMatchButton();
}
void
MyFrame::OnMatch(wxCommandEvent& event)
{
    cout << "OnMatch(): Id = " << event.GetId() << endl;

    // prompt to save results
    if (m_resultsArea->m_got_results &&
        ! m_resultsArea->m_saved_results)
    {
        if (! (wxYES == wxMessageBox(_("Discard current results?"), _("Match"),
                wxNO_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this)) )
        {
            return;
        }
    }

    // Get datasets. Make copies because we do not want modifications
    // (such as reading from the database) to be remembered in subsequent calculations
    ProfileRange c_set = cSet();
    ProfileRange r_set = rSet();

    if (c_set.size() == 0)
    {
        error << startl << "C-Set: no records selected (nothing to do)" << endl;
        return;
    }

    if (r_set.size() == 0)
    {
        error << startl << "R-Set: no records selected (nothing to do)" << endl;
        return;
    }

    // get parameters from Match tab
    m_match_params->crime_delta  = m_matchp->crimeDelta();
    m_match_params->ref_delta    = m_matchp->refDelta();
    m_match_params->lr_threshold = m_matchp->LRThreshold();

    vector<string> relatives;
    if (m_matchp->relIsSelected(0)) relatives.push_back("IDENT");
    if (m_matchp->relIsSelected(1)) relatives.push_back("D1");
    if (m_matchp->relIsSelected(2)) relatives.push_back("SIB");
    if (m_matchp->relIsSelected(3)) relatives.push_back("D2");
    if (m_matchp->relIsSelected(4)) relatives.push_back("D30");

    if (m_matchp->relIsSelected(5)) relatives.push_back("D14"); // D14
//    if (m_matchp->relIsSelected(5)) relatives.push_back("GEN[1/2,1/2,1/16,1/16]"); // D14
    if (m_matchp->relIsSelected(6)) relatives.push_back("INV[1/2,1/2,1/16,1/16]"); // I14

    if (m_matchp->relIsSelected(7)) relatives.push_back("GEN[1/4,1/4,1/8,1/8]");   // D23
    if (m_matchp->relIsSelected(8)) relatives.push_back("INV[1/4,1/4,1/8,1/8]");   // I23

    if (m_matchp->relIsSelected(9))  relatives.push_back("GEN[1/4,0,0,1/2]");      // sib-cousins
    if (m_matchp->relIsSelected(10)) relatives.push_back("GEN[1/16,0,0,1/16]");    // double first cousins

    for (int i=0; i<3; ++i)  // Extra relationships specified as strings
    {
    	string s = m_matchp->relString(i);
    	if (s != "") relatives.push_back(s);
    }

    m_match_params->match_types.clear();
    matchTypesFromString(m_match_params->match_types, relatives);

    // get parameters from subpopulation tab
    m_match_params->spm  = m_subpopmodel->subPopModel();
    m_match_params->mut  = m_subpopmodel->mutModel();

    // get parameters from Processors tab
    m_match_params->n_cuda_accel  = m_processors->nCudaAccel();
    m_match_params->n2_cuda_accel = m_processors->n2CudaAccel();

    // check we can do the requested stuff, and issue warnings if not

    wxString warnings;

    if (m_match_params->mut.type != MutModel::NO_MUTATIONS)
    {
    	// can do mutation model on one-to-one only
    	if ( (c_set.size() > 1) || (r_set.size() > 1) )
    	{
    		m_match_params->mut = MutModel(MutModel::NO_MUTATIONS, 0);
    		string s("Mutation model is supported for one-to-one matches only: *Mutation model disabled*.");
    		warnings << s << "\n\n";
        	warn << startl << s << endl;
    	}
    }

    if (m_match_params->spm.type != SubPopModel::HW)
    {
    	// can do a restricted set of familial matches with SPMs
    	std::vector<MatchType>::const_iterator it;
    	for (it = m_match_params->match_types.begin(); it != m_match_params->match_types.end(); ++it)
    	{
    		RelType type = it->m_rel_type;
			if (!(type == ident_t    ||
				  type == sibling_t  ||
				  type == degree_1_t ||
				  type == degree_2_t ||
				  type == degree_pq_t && it->m_path2steps == MatchType::INF)) // Dn
			{
	    		m_match_params->spm = SubPopModel(SubPopModel::HW, 0);
	    		ostringstream oss;
	    		oss << "Subpopulation model not currently supported for: " <<  it->string() << " using Hardy-Weinberg for *everything*.";
	    		warnings << oss.str() << "\n\n";
	        	warn << startl << oss.str() << endl;
			}
    	}
    }

    std::set<int> device_ids;

    int n_listed = m_processors->numGPUsListed();
    for (int i=0; i<n_listed; ++i)
    {
    	if (m_processors->gpuIsSelected(i))
    	{
    		device_ids.insert(i);
    	}
    }

    int n_enabled = gpuDevices().enableGPUs(device_ids);

    Assert2(n_enabled == (int)device_ids.size(), "failed to enable GPUs");

    if (n_enabled == 0 && (m_match_params->n_cuda_accel || m_match_params->n2_cuda_accel))
    {
		string s("no CUDA GPUs selected: using CPU.");
		warnings << s << "\n\n";
    	warn << startl << s << endl;
    }

    // if necessary, display a warning before proceeding
    if (warnings.size() == 0 ||
    	wxYES == wxMessageBox( warnings + " Continue?", _("Warning"), wxYES_DEFAULT | wxYES_NO | wxICON_EXCLAMATION, this))
    {

        // clear results
        m_match_results.clear();

        // display "working" message in results area
        m_resultsArea->clear();
        m_resultsArea->m_matchList->setMessage("Working...");
        wxSafeYield(); // yield control to pending messages so the text is displayed

		// do match
		cout << "c_set.size() = " << c_set.size() << " r_set.size() = " << r_set.size() << endl;

		do_match(
			m_match_results,
			c_set,
			r_set,
			refIsCrime(),
			*m_match_params);

        // Print matches to stdout
		cout << "c_set.size() = " << c_set.size() << " r_set.size() = " << r_set.size() << endl;

		// Display results
		m_resultsArea->displayMatches(c_set, r_set, m_match_params->match_types, m_match_params->lr_threshold, m_match_results);

		MatchResults results(m_match_params->match_types, m_match_results, m_resultsArea->m_c_ids, m_resultsArea->m_r_ids, m_match_params->lr_threshold);

        nm_output(results, true);
    }

    refreshMessageWindow();
}

//void MyFrame::OnSashDrag(wxSashEvent& event)
//{
//    cout << "OnSashDrag(): id = " << event.GetId() << endl;
//
//    if (event.GetDragStatus() == wxSASH_STATUS_OUT_OF_RANGE)
//        return;
//
//    cout << "OnSashDrag(): OK so far" << endl;
//
//    switch (event.GetId())
//    {
//        case ID_INPUT_AREA:
//        {
//            cout << "OnSashDrag(): setting size: " << event.GetDragRect().height << endl;
//            dynamic_cast<wxSashLayoutWindow*>(m_inputArea)->SetDefaultSize(wxSize(1000, event.GetDragRect().height));
////            m_inputArea->SetSize(wxSize(1000, event.GetDragRect().height));
//            break;
//        }
//    }
//
//    this->Refresh();
//
//}

BEGIN_EVENT_TABLE(ImportDialog, wxDialog)
    EVT_BUTTON(ID_BROWSE, ImportDialog::OnBrowse)
    EVT_TEXT(ID_FILE, ImportDialog::OnChange)
END_EVENT_TABLE();

void
ImportDialog::setup()
{
    // Widget hierarchy

    wxStaticText* label1 = new wxStaticText ( this, wxID_STATIC, wxT("File"), wxDefaultPosition, wxDefaultSize, 0 );
    wxStaticText* label2 = new wxStaticText ( this, wxID_STATIC, wxT("Dataset (default)"), wxDefaultPosition, wxDefaultSize, 0 );

    m_file           = new wxTextCtrl(this, ID_FILE, wxT(""), wxDefaultPosition, wxSize(150, wxDefaultCoord), 0 );
    wxButton* browse = new wxButton  (this, ID_BROWSE, wxT("Browse..."), wxDefaultPosition, wxDefaultSize, 0 );
    m_dataset        = new wxComboBox(this, wxID_ANY, _(""));

    wxStaticText* key1lbl = new wxStaticText ( this, wxID_STATIC, wxT("Profile key [1] in column"), wxDefaultPosition, wxDefaultSize, 0 );
    wxStaticText* key2lbl = new wxStaticText ( this, wxID_STATIC, wxT("Profile key [2] in column"), wxDefaultPosition, wxDefaultSize, 0 );

    m_key1 = new wxTextCtrl(this, wxID_ANY, wxT("1"), wxDefaultPosition, wxSize(50, wxDefaultCoord), 0 );
    m_key2 = new wxTextCtrl(this, wxID_ANY, wxT(""), wxDefaultPosition, wxSize(50, wxDefaultCoord), 0 );
    m_key1->Enable(false);
    m_key2->Enable(false);

    set<string> datasets;
    if (!db->listDatasets(datasets))
    {
        error << "ImportDialog::setup(): Datasets could not be read: " << endl;
        return;
    }

    setItemContainerStrings(m_dataset, datasets);

    wxStaticLine* line1 = new wxStaticLine ( this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    wxButton* okbtn = new wxButton ( this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0 );
    wxButton* cancel = new wxButton ( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );

    // Sizer hierarchy

    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
    wxBoxSizer *h1Sizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *h2Sizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *h3Sizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *h4Sizer = new wxBoxSizer( wxHORIZONTAL );
    wxBoxSizer *controlSizer = new wxBoxSizer( wxHORIZONTAL );

    topSizer->Add(h1Sizer, 0, wxEXPAND | wxALL/*, 10*/);
        h1Sizer->Add(label1, 1, wxEXPAND|wxALL, 5);
        h1Sizer->Add(m_file, 1, wxEXPAND|wxALL, 5);
        h1Sizer->Add(browse, 1, wxEXPAND|wxALL, 5);

    topSizer->Add(h2Sizer, 0, wxEXPAND | wxALL/*, 10*/);
        h2Sizer->Add(label2, 1, wxEXPAND|wxALL, 5);
        h2Sizer->Add(m_dataset, 2, wxEXPAND|wxALL, 5);

	topSizer->Add(h3Sizer, 0, wxEXPAND | wxALL/*, 10*/);
		h3Sizer->Add(key1lbl, 1, wxEXPAND|wxALL, 5);
		h3Sizer->Add(m_key1, 2, wxEXPAND|wxALL, 5);

	topSizer->Add(h4Sizer, 0, wxEXPAND | wxALL/*, 10*/);
		h4Sizer->Add(key2lbl, 1, wxEXPAND|wxALL, 5);
		h4Sizer->Add(m_key2, 2, wxEXPAND|wxALL, 5);

    topSizer->Add(line1, 0, wxGROW|wxALL, 5);

    topSizer->Add(controlSizer, 0, wxEXPAND | wxALL/*, 10*/);
        controlSizer->Add(okbtn, 0, wxALIGN_LEFT|wxALL, 5);
        controlSizer->Add(100, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);
        controlSizer->Add(cancel, 0, wxALIGN_RIGHT|wxALL, 5);

    SetSizer(topSizer);
    GetSizer()->Fit(this);
    GetSizer()->SetSizeHints(this);
}

void
ImportDialog::OnBrowse(wxCommandEvent& event)
{
    cout << "OnBrowse(): Id = " << event.GetId() << endl;

    wxString caption = _("Import profiles");
//    wxString wildcard = _("Profiles (*.idf.csv;*.p16.csv;*.b11.csv)|*.idf.csv;*.p16.csv;*.b11.csv");
    wxString wildcard = _("Profiles (*.csv)|*.csv");
    wxString defaultDir = wxString(homeDir(), *wxConvCurrent);
    wxString defaultFilename = wxEmptyString;

    wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN | wxFD_FILE_MUST_EXIST);

    if(dialog.ShowModal() == wxID_OK)
    {
        wxString path = dialog.GetPath();
        string filename(path.ToUTF8());

        m_file->SetValue(path);

        if (m_dataset->GetValue().empty())
        {
            string dataset = baseName(filename);
            m_dataset->SetValue(wxString(dataset.c_str(), *wxConvCurrent));
        }
    }
}

void
ImportDialog::OnChange(wxCommandEvent& event)
{
    cout << "OnChange(): Id = " << event.GetId() << endl;

    wxString path = m_file->GetValue();
    string filename(path.ToUTF8());
    ProfileFilter::FileType ftype = fileType(filename);

    switch(ftype)
    {
        case ProfileFilter::B11:         // 1, 2
        case ProfileFilter::Identifiler: // 1, 2
        	m_key1->SetValue(_("1")); m_key1->Enable(false);
        	m_key2->SetValue(_("2")); m_key2->Enable(false);
        	break;

        case ProfileFilter::PowerPlex16: // 1
        	m_key1->SetValue(_("1")); m_key1->Enable(false);
        	m_key2->SetValue(_(""));  m_key2->Enable(false);
        	break;

        case ProfileFilter::other_csv:   // user defined
        	m_key1->SetValue(_("1")); m_key1->Enable(true);
        	m_key2->SetValue(_(""));  m_key2->Enable(true);
        	break;

        case ProfileFilter::SGM_Plus:
        case ProfileFilter::type_unknown:
        default:
            error << startl << "unrecognised file type: " << filename << endl;
            error << alignl << "(should be *.b11.csv *.idf.csv *.p16.csv or *.csv)" << endl;
            return;
    }

}

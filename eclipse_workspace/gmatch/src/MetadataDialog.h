/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MetadataDialog.h
 *
 *  Created on: Jul 29, 2010
 *      Author: gareth
 */

#ifndef METADATADIALOG_H_
#define METADATADIALOG_H_

#include "fand/Database.h"

#include "wx/wx.h"

class MetadataDialog : public wxDialog
{
public:
    // Constructor
    MetadataDialog(
            wxWindow *parent,
            wxWindowID winid = wxID_ANY,
            const wxString& caption = _("Metadata"),
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
            long style = wxCAPTION | wxRESIZE_BORDER | wxSYSTEM_MENU )
    : wxDialog(parent, winid, caption, pos, size, style)
    {
        setup();
    }

    void setEditable(bool edit);
    void populate(MetaData const &metad);
    void apply();

    void readback(MetaData &metad);

    DECLARE_EVENT_TABLE()
//    DECLARE_CLASS()

    void OnApply(wxCommandEvent& event);

private:
    void setup();

	struct MetaControl : public wxPanel
	{
		MetaControl(MetadataDialog *parent, const std::string &name, const std::string &field);
		void setup();
		wxTextCtrl    *text;

		std::string name;
		std::string field;
	};

	std::map<std::string, MetaControl*> m_controls;
};

#endif /* METADATADIALOG_H_ */

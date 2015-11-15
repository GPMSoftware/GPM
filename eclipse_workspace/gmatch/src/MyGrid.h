/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MyGrid.h
 *
 *  Created on: Nov 18, 2010
 *      Author: gareth
 */

#ifndef MYGRID_H_
#define MYGRID_H_

#include "wx/wx.h"
#include "wx/grid.h"

class MyGridCellEditor :public wxGridCellTextEditor
{
public:

    MyGridCellEditor() {}

//    virtual void Create(wxWindow* parent,
//                        wxWindowID id,
//                        wxEvtHandler* evtHandler);

//    void OnChar(wxKeyEvent& event);

//    DECLARE_EVENT_TABLE()

};

class MyGrid : public wxGrid
{
public:

    MyGrid( wxWindow *parent,
        wxWindowID id,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxWANTS_CHARS,
        const wxString& name = wxGridNameStr )
    : wxGrid(parent, id, pos, size, style, name)
      {}

    void Clear() { ClearGrid(); }
    void setNumberRows(int rows);
};

#endif /* MYGRID_H_ */

/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MyGrid.cpp
 *
 *  Created on: Nov 18, 2010
 *      Author: gareth
 */

#include "MyGrid.h"

void
MyGrid::setNumberRows(int rows)
{
    int n = GetNumberRows();
    if (n > rows)
    {
        DeleteRows(0, n-rows);
    }
    else if (n < rows)
    {
        InsertRows(0, rows-n);
    }
}

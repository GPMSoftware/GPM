/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Assert.cpp
 *
 *  Created on: Dec 15, 2009
 *      Author: gareth
 */

#include "Assert.h"
#include <UnitTest++/AssertException.h>

#include "MessageStream.h"
INIT_MESSAGES("assert")
#include "messages.h"

using namespace std;

// the purpose of this function is that you can put a breakpoint in it

void
throwAssertException(char const* description, char const* filename, int lineNumber)
{
	error << "AssertException thrown at " << filename << ": " << lineNumber << ": " << description << endl;
    cout << "AssertException thrown at " << filename << ": " << lineNumber << ": " << description << endl;
	throw UnitTest::AssertException(description, filename, lineNumber);
}

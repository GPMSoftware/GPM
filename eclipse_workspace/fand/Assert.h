/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Assert.h
 *
 *  Created on: Nov 30, 2009
 *      Author: gareth
 */

#ifndef ASSERT_H_
#define ASSERT_H_

// an assert macro suitable for use with UnitTest++
// Assert must throw an AssertException
// An Assert failure will then cause a test fail, unless it occurs inside the CHECK_ASSERT check macro
//

void
throwAssertException(char const* description, char const* filename, int lineNumber);

#define Assert( exp ) \
	if (!(exp)) throwAssertException("Assert: " #exp, __FILE__, __LINE__);

#define Assert2( exp, desc ) \
	if (!(exp)) throwAssertException("Assert: " #exp ": " desc, __FILE__, __LINE__);

#endif /* ASSERT_H_ */

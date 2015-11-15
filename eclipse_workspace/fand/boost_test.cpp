/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * boost_test.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: gareth
 */

#include <UnitTest++/UnitTest++.h>
#include <boost/regex.hpp>

#include <iostream>
#include <string>

//TEST(boost0)
//{
//
//    boost::regex pat( "xyz" );
//    boost::smatch matches;
//
//    std::string line = "hello.xyz\n";
//    CHECK ( boost::regex_match(line, matches, pat) ); // CRASH!
//
////    CHECK (! boost::regex_match(std::string("hello.xyz"), matches, pat));
//}

TEST(boost1)
{
    boost::shared_ptr<int> p(new int(42));
    CHECK_EQUAL(42, *p);
    *p = 1;
    CHECK_EQUAL(1, *p);
}

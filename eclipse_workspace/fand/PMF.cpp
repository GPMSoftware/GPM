/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * PMF.cpp
 *
 *  Created on: Dec 4, 2009
 *      Author: gareth
 */

#include "PMF.h"
#include <UnitTest++/UnitTest++.h>

PMF<std::string>::POD data[] =
{
  { "one",   0.2 }
, { "two",   0.4 }
, { "three", 0.6 }
, { "four",  0.8 }
};

//TEST (PMF)
//{
//	CHECK_EQUAL(0, PMF<std::string>::m_created);
//	CHECK_EQUAL(0, PMF<std::string>::m_destroyed);
//	CHECK_EQUAL(0, PMF<std::string>::count());
//}

TEST (PMF0)
{
	PMF<std::string> empty;
	CHECK(empty.begin() == empty.end());
	CHECK_EQUAL(false, empty.normalize());

	std::string s;
	CHECK_THROW(s = empty.mode(), std::out_of_range);
	CHECK_THROW(s = empty.pick(), std::out_of_range);
}

TEST (PMF1)
{
	PMF<std::string> singleton;
	singleton["only"] = 0.2;

	CHECK_EQUAL(1, singleton.size());
	CHECK_EQUAL(0.2, singleton["only"]);
	CHECK_EQUAL("only", singleton.mode());
	CHECK_EQUAL("only", singleton.pick());
}

TEST(PMF2)
{
	PMF<std::string> pmf(data, sizeof(data)/sizeof(PMF<std::string>::POD));

	pmf.normalize();

	CHECK_EQUAL(0.1, pmf["one"]);
	CHECK_EQUAL(0.2, pmf["two"]);
	CHECK_EQUAL(0.3, pmf["three"]);
	CHECK_EQUAL(0.4, pmf["four"]);

	PMF<std::string> results;
	const int N = 100;
	for(int i=0; i<N; ++i)
	{
		std::string s = pmf.pick();
		results[s]++;
	}

	CHECK_EQUAL(100, results["one"] + results["two"] + results["three"] + results["four"]);

	results.normalize();
	CHECK_CLOSE(1, results["one"] + results["two"] + results["three"] + results["four"], 1e-15);

	CHECK(results["one"]   < results["two"]);
	CHECK(results["two"]   < results["three"]);
	CHECK(results["three"] < results["four"]);

	results.normalize();

	CHECK(pmf.size() == 4);
	CHECK(pmf["not there"] == 0);
	CHECK(pmf.size() == 5);
}

TEST(PMF3)
{
	PMF<std::string> pmf1(data, sizeof(data)/sizeof(PMF<std::string>::POD));
	PMF<std::string> pmf2;

	pmf2["one"] = 0.6;
	pmf2["three"]  = 0.8;

	pmf1 += pmf2;

	CHECK_EQUAL(0.2 + 0.6, pmf1["one"]);
	CHECK_EQUAL(0.4,       pmf1["two"]);
	CHECK_EQUAL(0.6 + 0.8, pmf1["three"]);
}


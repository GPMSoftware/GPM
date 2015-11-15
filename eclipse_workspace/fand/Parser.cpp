/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Parser.cpp
 *
 *  Created on: Jun 14, 2010
 *      Author: gareth
 */

#include "populationdata.h"
#include "Parser.h"
//#include "y.tab.h" // nothing we need in here

#include <UnitTest++/UnitTest++.h>

// functions defined in the YACC grammar
void start_scan_bytes(char const *buf, size_t size);
void end_scan_bytes();
int yyparse(void);

// static definitions
PMF<Allele> *Parser::m_pmf = 0;
PopulationData *Parser::m_pop = 0;
Locus Parser::m_loc;
Parser::PMFPool Parser::m_pmf_pool;

pthread_mutex_t parser_mutex = PTHREAD_MUTEX_INITIALIZER;

Parser::Parser()
{
}

Parser::~Parser()
{
}

Parser &
Parser::theParser()
{
	static Parser p;
	return p;
}

int
Parser::parse(
	std::string const &s,
	PMF<Allele> &pmf_out,
	PopulationData *popdata,
	Locus loc)
{
    pthread_mutex_lock( &parser_mutex );

	pmf_out.clear();
	clearPMFPool();

	m_pmf  = &pmf_out;
	m_pop = popdata;
	m_loc = loc;

	// set up lex input to point to string
	start_scan_bytes(s.data(), s.size());

	// call parser (which writes the result into m_pmf)
	int ret = yyparse();

	//clean up
	end_scan_bytes();

    pthread_mutex_unlock( &parser_mutex );

	return ret;
}

double
Parser::background(Allele const &a)
{
	if (m_pop && m_loc != locus_none)
	{
		return m_pop->getFrequency(m_loc, a);
	}

	return 1; // we are not interested in the probability
}

PMF<Allele>*
Parser::newPMF()
{
	// Look for an unused PMF
	PMFPool::iterator it;
	for (it = m_pmf_pool.begin(); it != m_pmf_pool.end(); ++it)
	{
		if (it->second == false)
		{
		    // Mark as in use and return
			it->second = true;
			return it->first;
		}
	}

	// None free: add a new one to the pool and return it
	PMF<Allele> *p = new PMF<Allele>;
	m_pmf_pool.insert(std::make_pair(p, true));
	return p;
}

void
Parser::clearPMFPool()
{
	PMFPool::iterator it;
	for (it = m_pmf_pool.begin(); it != m_pmf_pool.end(); ++it)
	{
		it->first->clear();
		it->second = false;
	}
}

void
Parser::deletePMF(PMF<Allele> *pmf)
{
	// find the PMF, clear it, and mark as unused
	PMFPool::iterator it = m_pmf_pool.find(pmf);
	Assert2(it != m_pmf_pool.end(), "Parser::deletePMF: deleting a PMF that was not allocated");
	Assert2(it->second, "Parser::deletePMF: deleting a PMF that was not in use");

	it->first->clear();
	it->second = false;
}


TEST(Parser1)
{
//	CHECK_EQUAL(0, PMF<Allele>::m_created);
//	CHECK_EQUAL(0, PMF<Allele>::m_destroyed);
//	CHECK_EQUAL(0, PMF<Allele>::count());


	Parser &p = Parser::theParser();
	PMF<Allele> pmf;
	CHECK(p.parse("", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("F@.1", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("B@.1", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("a", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11&.1", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@0", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@1", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@1.0", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@.1/", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11/F", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("(11@.1", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11@.1)", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("D/11@.1)", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("11/D)", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("D@B)", pmf, 0) == Parser::yacc_error);
	CHECK(p.parse("10/(11/(12/(13@.1)@.2)@.3))", pmf, 0) == Parser::yacc_error);
}

TEST(Parser2)
{
	Parser &p = Parser::theParser();
	PMF<Allele> pmf;
	PMF<Allele> locus_pmf;
	locus_pmf[11] = 0.1;
	locus_pmf[12] = 0.2;
	locus_pmf[13] = 0.3;
	PopulationData back;
	Locus loc = (Locus)1;
	back.setLocus(loc, locus_pmf);

	CHECK(p.parse("F", pmf, 0) == Parser::yacc_ok);
	CHECK(pmf.empty());
	pmf.clear();

	CHECK(p.parse("11", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(11), pmf.begin()->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	pmf.clear();

	CHECK(p.parse("11.2", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(11, 2), pmf.begin()->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	pmf.clear();

	CHECK(p.parse("11.2@.1", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(11, 2), pmf.begin()->first);
	CHECK_CLOSE(0.1, pmf.begin()->second, 1e-6);
	pmf.clear();

	CHECK(p.parse("11.2@0.1", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(11, 2), pmf.begin()->first);
	CHECK_CLOSE(0.1, pmf.begin()->second, 1e-6);
	pmf.clear();

	// No longer legal
//	CHECK(p.parse("11.2@0.1@.2", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(1, pmf.size());
//	CHECK_EQUAL(Allele(11, 2), pmf.begin()->first);
//	CHECK_CLOSE(0.02, pmf.begin()->second, 1e-6);
//	pmf.clear();

	CHECK(p.parse("11.2@.123", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(11, 2), pmf.begin()->first);
	CHECK_CLOSE(0.123, pmf.begin()->second, 1e-6);
	pmf.clear();

	CHECK(p.parse("D", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(Allele::unknown), pmf.begin()->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	pmf.clear();

	CHECK(p.parse("D@0.9", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(1, pmf.size());
	CHECK_EQUAL(Allele(Allele::unknown), pmf.begin()->first);
	CHECK_CLOSE(0.9, pmf.begin()->second, 1e-6);
	pmf.clear();

	CHECK(p.parse("11.2@.2/13@.25", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(2, pmf.size());
	PMF<Allele>::iterator it = pmf.begin();
	CHECK_EQUAL(Allele(11, 2), it->first);
	CHECK_CLOSE(0.2, it->second, 1e-6);
	CHECK_EQUAL(Allele(13), (++it)->first);
	CHECK_CLOSE(0.25, it->second, 1e-6);
	pmf.clear();

	CHECK(p.parse("11/12/13", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(3, pmf.size());
	it = pmf.begin();
	CHECK_EQUAL(Allele(11), it->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	CHECK_EQUAL(Allele(12), (++it)->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	CHECK_EQUAL(Allele(13), (++it)->first);
	CHECK_EQUAL(1, pmf.begin()->second);
	pmf.clear();

    // No longer legal
//	CHECK(p.parse("(11/12/13)", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(3, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(11), it->first);
//	CHECK_EQUAL(1, pmf.begin()->second);
//	CHECK_EQUAL(Allele(12), (++it)->first);
//	CHECK_EQUAL(1, pmf.begin()->second);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_EQUAL(1, pmf.begin()->second);
//	pmf.clear();

    // No longer legal
//	CHECK(p.parse("(11/12/13)@0.5", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(3, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(11), it->first);
//	CHECK_CLOSE(0.5, it->second, 1e-6);
//	CHECK_EQUAL(Allele(12), (++it)->first);
//	CHECK_CLOSE(0.5, it->second, 1e-6);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_CLOSE(0.5, it->second, 1e-6);
//	pmf.clear();

	CHECK(p.parse("(11/12/13)@B", pmf, &back, loc) == Parser::yacc_ok);
	CHECK_EQUAL(3, pmf.size());
	it = pmf.begin();
	CHECK_EQUAL(Allele(11), it->first);
	CHECK_CLOSE(1.0/6, it->second, 1e-6);
	CHECK_EQUAL(Allele(12), (++it)->first);
	CHECK_CLOSE(2.0/6, it->second, 1e-6);
	CHECK_EQUAL(Allele(13), (++it)->first);
	CHECK_CLOSE(3.0/6, it->second, 1e-6);
	pmf.clear();

    // No longer legal
//	CHECK(p.parse("(11/12/13)@B/(14/15)@.2", pmf, &back, loc) == Parser::yacc_ok);
//	CHECK_EQUAL(5, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(11), it->first);
//	CHECK_CLOSE(0.1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(12), (++it)->first);
//	CHECK_CLOSE(0.2, it->second, 1e-6);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_CLOSE(0.3, it->second, 1e-6);
//	CHECK_EQUAL(Allele(14), (++it)->first);
//	CHECK_CLOSE(0.2, it->second, 1e-6);
//	CHECK_EQUAL(Allele(15), (++it)->first);
//	CHECK_CLOSE(0.2, it->second, 1e-6);
//	pmf.clear();
//
//	CHECK(p.parse("(11/12/13) @B / ( 14 / 15 ) @ .2", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(5, pmf.size());
//	pmf.clear();

	CHECK(p.parse("11.2@.2/10@.2/13@.25", pmf, 0) == Parser::yacc_ok);
	CHECK_EQUAL(3, pmf.size());
	it = pmf.begin();
	CHECK_EQUAL(Allele(10), it->first);
	CHECK_CLOSE(0.2, it->second, 1e-6);
	CHECK_EQUAL(Allele(11, 2), (++it)->first);
	CHECK_CLOSE(0.2, it->second, 1e-6);
	CHECK_EQUAL(Allele(13), (++it)->first);
	CHECK_CLOSE(0.25, it->second, 1e-6);
	pmf.clear();

    // No longer legal
//	CHECK(p.parse("(11.2@.2/10@0.2/13@.25)@.5", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(3, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(10), it->first);
//	CHECK_CLOSE(0.1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(11, 2), (++it)->first);
//	CHECK_CLOSE(0.1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_CLOSE(0.125, it->second, 1e-6);
//	pmf.clear();

    // No longer legal
//	CHECK(p.parse("(13@.25/(11.2/(10))@0.2)@.5", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(3, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(10), it->first);
//	CHECK_CLOSE(0.1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(11, 2), (++it)->first);
//	CHECK_CLOSE(0.1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_CLOSE(0.125, it->second, 1e-6);
//	pmf.clear();
//
    // No longer legal
//	CHECK(p.parse("10/(11/(12/(13@.1)@.2)@.3)", pmf, 0) == Parser::yacc_ok);
//	CHECK_EQUAL(4, pmf.size());
//	it = pmf.begin();
//	CHECK_EQUAL(Allele(10), it->first);
//	CHECK_CLOSE(1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(11), (++it)->first);
//	CHECK_CLOSE(1, it->second, 1e-6);
//	CHECK_EQUAL(Allele(12), (++it)->first);
//	CHECK_CLOSE(0.3, it->second, 1e-6);
//	CHECK_EQUAL(Allele(13), (++it)->first);
//	CHECK_CLOSE(0.006, it->second, 1e-6);
//	pmf.clear();

}

TEST(Parser3)
{
    Parser &p = Parser::theParser();
    PMF<Allele> pmf;
    PMF<Allele> locus_pmf;
    locus_pmf[Allele(14,3)] = 0.1;
    locus_pmf[Allele(15,3)] = 0.2;
    locus_pmf[Allele(16,3)] = 0.3;
    locus_pmf[Allele(17,3)] = 0.4;
    PopulationData back;
    Locus loc = (Locus)1;
    back.setLocus(loc, locus_pmf);

    CHECK(p.parse("14.3/15.3/16.3", pmf, 0) == Parser::yacc_ok);
//    CHECK(!pmf.empty());
    CHECK_EQUAL(3, pmf.size());
    CHECK_EQUAL(3, pmf.sum()); // pmf is unnormalized
    pmf.clear();

    // background is not added at this stage?
    CHECK(p.parse("14.3/15.3/16.3@0.6", pmf, 0) == Parser::yacc_ok);
    CHECK_EQUAL(3, pmf.size());
    CHECK_CLOSE(2.6, pmf.sum(), 1e-6); // pmf is unnormalized
    pmf.clear();

}


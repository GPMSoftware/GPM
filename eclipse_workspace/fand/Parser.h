/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Parser.h
 *
 *  Created on: Jun 14, 2010
 *      Author: gareth
 *
 *      Parser provides an interface between the YACC grammar that parses the B11 input language
 *      and the rest of the application.
 *
 *      NB Parser is s singleton class with mostly static functions and data. This is because the
 *      YACC parser is essentially a singleton - it uses global variables to maintain state.
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <string>
#include "PMF.h"
#include "Allele.h"

class PopulationData;

#include "loci.h"

class Parser {
public:
	enum { yacc_error = 1, yacc_ok = 0 };

	// call yacc grammar to parse a string s containing a single allele specification, such as
	// "(11/12/13)@.2/D@.8"
	// The result is returned in pmf_out.
	// NB any 'D' is represented by Allele::unknown (with the correct frequency).
	// It is up to the caller to test for and interpret 'D'.
	int parse(std::string const &s, PMF<Allele> &pmf_out, PopulationData *popdata = 0, Locus loc = locus_none);

	// get background frequency of allele (called by YACC grammar)
	static double background(Allele const &a);

	// return a new Allele PMF (called by YACC grammar)
	static PMF<Allele>* newPMF();

	// dispose of Allele PMF (called by YACC grammar)
	static void deletePMF(PMF<Allele> *pmf);

	// return the singleton parser
	static Parser &theParser();

	static PMF<Allele> *m_pmf;     // output PMF (set by YACC grammar)
	static PopulationData *m_pop;  // background frequency distribution
	static Locus m_loc;            // locus

private:
	Parser();
	virtual ~Parser();

	// A small pool of PMFs to allocate to the yacc grammar
	// The bool indicates 'in use'
	// NB these are never deleted, but kept for future use.
	// The number in the pool should stay small
	typedef std::map< PMF<Allele>*, bool > PMFPool;
	static PMFPool m_pmf_pool;

	// Clear all PMFs in the pool and mark as unused.
	// (Use of this function means the parser does not have to free up its
	// resources in the case of an error - and keeps the grammar simple).
	static void clearPMFPool();
};


#endif /* PARSER_H_ */

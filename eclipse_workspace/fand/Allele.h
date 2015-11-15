/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Allele.h
 *
 *  Created on: Dec 22, 2009
 *      Author: gareth
 */

#ifndef ALLELE_H_
#define ALLELE_H_

#include <string>
#include <iostream>

// NB FGA has max 51 repeats with 80 observed alleles
//    D21S11 has max 41 repeats with 89 observed alleles.
//
struct Allele
{
	enum { unknown = 0, X = 106, Y = 112 }; // Amelogenin: X = 106bp, Y=112bp, not an STR!

	// NB this c'tor provides conversion from int.
	Allele(unsigned int repeats = unknown, unsigned int variant = 0);

	std::string string() const; // output as string, e.g. "26.2"

	bool operator<(const Allele &other) const;
	bool operator==(const Allele &other) const;
	bool operator!=(const Allele &other) const;

	unsigned char m_repeats;
	unsigned char m_variant;
};

std::ostream &
operator<<(std::ostream &os, const Allele &a);

std::istream &
operator>>(std::istream &os, Allele &a);

#endif /* ALLELE_H_ */

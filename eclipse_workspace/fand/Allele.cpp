/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Allele.cpp
 *
 *  Created on: Dec 22, 2009
 *      Author: gareth
 */

#include "Allele.h"
#include "PMF.h"
//#include "CheckMacros.h" // hacked version of UnitTest++ header
#include <UnitTest++/UnitTest++.h>
#include "Assert.h"
#include <iostream>
#include <sstream>
#include <algorithm>

Allele::Allele(unsigned int repeats, unsigned int variant)
: m_repeats(repeats)
, m_variant(variant)
{
	Assert2(repeats < 256, "Allele::Allele(unsigned int repeats, unsigned int variant): repeat number too large");
	Assert2(variant < 256, "Allele::Allele(unsigned int repeats, unsigned int variant): variant number too large");
}

std::ostream &
operator<<(std::ostream &os, const Allele &a)
{
	os << a.string();
	return os;
}

std::istream &
operator>>(std::istream &is, Allele &a)
{
    std::string s;
    is >> s;

    if (s == "X")
    {
        a = Allele::X;
    }
    else if (s == "Y")
    {
        a = Allele::Y;
    }
    else
    {
        // Allele names look like repeat.variant e.g. 16 or 16.2
        // Read as a double and extract the integer values
        double allele_name = atof(s.c_str());
        int repeats = int(allele_name);
        int variant = int( (allele_name - repeats + 0.01) * 10); // nasty!

        a = Allele(repeats, variant);
    }

    return is;
}

TEST(Allele_err)
{
//  does two-byte allele result in more memory use?
//	std::cout << "sizeof(Allele)      = " << sizeof(Allele) << std::endl;          // 2
//	std::cout << "sizeof(Allele[10])  = " << sizeof(Allele[10]) << std::endl;      // 20
//	std::cout << "sizeof(PMF<Allele>) = " << sizeof(PMF<Allele>) << std::endl;     // 48 (same as with 1-byte allele!)

    CHECK_ASSERT( Allele a4(26, 256) );
    CHECK_ASSERT( Allele a5(256) );
}

std::string
Allele::string() const
{
	switch(m_repeats)
	{
	case unknown:
		return "F";
	case X:
		return "X";
	case Y:
		return "Y";
	default:
		{
			std::ostringstream ss;
			ss << (int)m_repeats;
			if (m_variant)
			{
				ss << "." << (int)m_variant;
			}
			return ss.str();
		}
	}
}

TEST(Allele_string)
{
	Allele a0;
    CHECK_EQUAL("F", a0.string());

    Allele a1(26);
    CHECK_EQUAL("26", a1.string());

	Allele a2(26, 0);
    CHECK_EQUAL("26", a2.string());

    Allele a3(26, 2);
    CHECK_EQUAL("26.2", a3.string());

    Allele a4(Allele::X);
    CHECK_EQUAL("X", a4.string());

    Allele a5(Allele::Y, 2);
    CHECK_EQUAL("Y", a5.string());

    Allele a6(Allele::unknown);
    CHECK_EQUAL("F", a6.string());
}

bool
Allele::operator<(const Allele &other) const
{
	if (m_repeats == other.m_repeats)
	{
		return m_variant < other.m_variant;
	}
	else
	{
		return m_repeats < other.m_repeats;
	}
}

bool
Allele::operator==(const Allele &other) const
{
	return m_repeats == other.m_repeats &&
		   m_variant == other.m_variant;
}

bool
Allele::operator!=(const Allele &other) const
{
	return ! (*this == other);
}

TEST(Allele_less)
{
    Allele a0(25);
    Allele a1(26);
	Allele a2(26, 1);
	Allele a3(26, 2);
    CHECK(a0 < a1);
    CHECK(!(a1 < a0));
    CHECK(a0 < a2);
    CHECK(!(a2 < a0));
    CHECK(a1 < a2);
    CHECK(!(a2 < a1));
    CHECK(a2 < a3);
    CHECK(!(a3 < a2));
}




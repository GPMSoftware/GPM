/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * SubPopModel.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: gareth
 */

#include "SubPopModel.h"
#include "Assert.h"

// returns P(p1,p2) according to the Sub-population Model (default HW)
double
SubPopModel::prob(bool homozygote, double p1, double p2) const
{
    Assert (p1==p2 || !homozygote);

    // default: HW
    double val = homozygote ? p1*p2 : 2*p1*p2;

    switch (type)
    {
    case SubPopModel::HW:
        break;

    case SubPopModel::NRC4_4:
    case SubPopModel::B11: // Uses NRC4_4 frequencies for background

        if (homozygote)
        {
            val = p1*p1 + p1*(1-p1)*theta_bar;
        }
        else
        {
            // NRC 4.1 says to use HW in heterozygote case
            // but we use 4.4 (properly normalized form)
            val = 2*p1*p2 * (1-theta_bar);
        }
        break;

    case SubPopModel::NRC4_10:
    {
        // optimization: calculate constants once
        double F  = theta_bar;
        double F_ = 1-F;
        double denom = (1+F)*(1+2*F);

        if (homozygote)
        {
            val = ((2*F + F_*p1) * (3*F + F_*p1)) / denom; // NRC 4.2a
        }
        else
        {
            val = 2*((F + F_*p1) * (F + F_*p2)) / denom;   // NRC 4.2b
        }
        break;
    }

    default:
        Assert2(false, "SubPopModel::operator(): Unknown Sub-population model");
    }

    return val;
}

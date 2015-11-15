/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * geotrans.h
 *
 *  Created on: Sep 15, 2011
 *      Author: gareth
 */

#ifndef GEOTRANS_H_
#define GEOTRANS_H_

#include "Precision.h"
#include <string>

void convertMGRSToGeodetic(
        const char* mgrsString,
        double& lat,
        double& lon,
        double& height);

void convertGeodeticToMGRS(
        double lat,
        double lon,
        double height,
        char* &mgrsString,
        MSP::CCS::Precision::Enum& precision);

struct LatLonStrings
{
    std::string lat;
    std::string lon;
};

LatLonStrings
latLonStrings(std::string const &mgrs);

#endif /* GEOTRANS_H_ */


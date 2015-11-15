/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * geotrans.cpp
 *
 *  Created on: Sep 15, 2011
 *      Author: gareth
 *
 *      Functions using the Geotrans CCS library: http://earth-info.nga.mil/GandG/geotrans/
 */

#include "geotrans.h"

#include "CoordinateConversionService.h"
#include "CoordinateSystemParameters.h"
#include "GeodeticParameters.h"
#include "CoordinateTuple.h"
#include "GeodeticCoordinates.h"
#include "CartesianCoordinates.h"
#include "Accuracy.h"
#include "MGRSorUSNGCoordinates.h"
#include "UTMParameters.h"
#include "UTMCoordinates.h"
#include "CoordinateType.h"
#include "HeightType.h"
#include "CoordinateConversionException.h"

// where those headers really are:
//
//#include "CoordinateConversion/CoordinateConversionService.h"
//#include "dtcc/CoordinateSystemParameters/CoordinateSystemParameters.h"
//#include "dtcc/CoordinateSystemParameters/GeodeticParameters.h"
//#include "dtcc/CoordinateTuples/CoordinateTuple.h"
//#include "dtcc/CoordinateTuples/GeodeticCoordinates.h"
//#include "dtcc/CoordinateTuples/CartesianCoordinates.h"
//#include "dtcc/CoordinateTuples/Accuracy.h"
//#include "dtcc/CoordinateTuples/MGRSorUSNGCoordinates.h"
//#include "dtcc/CoordinateSystemParameters/UTMParameters.h"
//#include "dtcc/CoordinateTuples/UTMCoordinates.h"
//#include "dtcc/Enumerations/CoordinateType.h"
//#include "dtcc/Enumerations/HeightType.h"
//#include "dtcc/Exception/CoordinateConversionException.h"

#include "util.h"

#include <iostream>
#include <UnitTest++/UnitTest++.h>

#include "fand/MessageStream.h"
INIT_MESSAGES("geotrans");
#include "fand/messages.h"

using namespace std;
using MSP::CCS::Precision;
using MSP::CCS::CoordinateConversionException;

// This function is copied from CCS/SampleCode/testCoordinateConversionSample.cpp
void convertGeocentricToMGRS(
        const double x,
        const double y,
        const double z,
        char*& mgrsString,
        Precision::Enum& precision)
{
    try
    {
        MSP::CCS::CoordinateSystemParameters geocentricParameters(MSP::CCS::CoordinateType::geocentric);
        MSP::CCS::CoordinateSystemParameters mgrsParameters(MSP::CCS::CoordinateType::militaryGridReferenceSystem);
        MSP::CCS::CoordinateConversionService ccs( "WGE", &geocentricParameters, "WGE", &mgrsParameters );
        MSP::CCS::Accuracy sourceAccuracy;
        MSP::CCS::Accuracy targetAccuracy;
        MSP::CCS::CartesianCoordinates sourceCoordinates(MSP::CCS::CoordinateType::geocentric, x, y, z);
        MSP::CCS::MGRSorUSNGCoordinates targetCoordinates;
        ccs.convertSourceToTarget( &sourceCoordinates, &sourceAccuracy, targetCoordinates, targetAccuracy );
        mgrsString = targetCoordinates.MGRSString();
        precision = targetCoordinates.precision();
    }
    catch(MSP::CCS::CoordinateConversionException e)
    {
        warn << startl << "CoordinateConversionException: " << e.getMessage() << endl;
        throw e;
    }
}

TEST(geotrans_convertGeocentricToMGRS)
{
    char *mgrs;
    Precision::Enum precision;

    try
    {
        convertGeocentricToMGRS(6378137, 0, 0, mgrs, precision); // the point where the Greenwich meridian crosses the equator
//        cout << "convertGeocentricToMGRS: MGRS = " << mgrs << " Precision = " << precision << endl;

        string expected("31NAA6602100000");
        string actual(mgrs);
        CHECK_EQUAL(expected, actual);
    }
    catch (...)
    {
        CHECK(false);
    }
}

// This function is copied from CCS/SampleCode/testCoordinateConversionSample.cpp
void convertGeocentricToGeodetic(
        const double x,
        const double y,
        const double z,
        double& lat,
        double& lon,
        double& height)
{
    try
    {
        MSP::CCS::CoordinateSystemParameters geocentricParameters(MSP::CCS::CoordinateType::geocentric);
        MSP::CCS::GeodeticParameters geodeticParameters(MSP::CCS::CoordinateType::geodetic, MSP::CCS::HeightType::ellipsoidHeight);
        MSP::CCS::CoordinateConversionService ccs( "WGE", &geocentricParameters, "WGE", &geodeticParameters );
        MSP::CCS::Accuracy sourceAccuracy;
        MSP::CCS::Accuracy targetAccuracy;
        MSP::CCS::CartesianCoordinates sourceCoordinates(MSP::CCS::CoordinateType::geocentric, x, y, z);
        MSP::CCS::GeodeticCoordinates targetCoordinates;
        ccs.convertSourceToTarget( &sourceCoordinates, &sourceAccuracy, targetCoordinates, targetAccuracy );
        lat = targetCoordinates.latitude();
        lon = targetCoordinates.longitude();
        height = targetCoordinates.height();
    }
    catch(MSP::CCS::CoordinateConversionException e)
    {
        warn << startl  << "CoordinateConversionException: " << e.getMessage() << endl;
        throw e;
    }
}

// This function is written by DGW
void convertMGRSToGeodetic(
        const char* mgrsString,
        double& lat,
        double& lon,
        double& height)
{
    try
    {
        MSP::CCS::CoordinateSystemParameters mgrsParameters(MSP::CCS::CoordinateType::militaryGridReferenceSystem);
        MSP::CCS::GeodeticParameters geodeticParameters(MSP::CCS::CoordinateType::geodetic, MSP::CCS::HeightType::ellipsoidHeight);
        MSP::CCS::CoordinateConversionService ccs("WGE", &mgrsParameters, "WGE", &geodeticParameters);
        MSP::CCS::Accuracy sourceAccuracy;
        MSP::CCS::Accuracy targetAccuracy;

        MSP::CCS::MGRSorUSNGCoordinates sourceCoordinates(MSP::CCS::CoordinateType::militaryGridReferenceSystem, mgrsString);
        MSP::CCS::GeodeticCoordinates targetCoordinates;

        ccs.convertSourceToTarget( &sourceCoordinates, &sourceAccuracy, targetCoordinates, targetAccuracy );
        lat = targetCoordinates.latitude();
        lon = targetCoordinates.longitude();
        height = targetCoordinates.height();
    }
    catch(MSP::CCS::CoordinateConversionException e)
    {
        warn << startl  << "CoordinateConversionException: " << e.getMessage() << endl;
        throw e;
    }
}

LatLonStrings
latLonStrings(string const &mgrs)
{
    LatLonStrings ret;

    try
    {
        double lat, lon, height;

        convertMGRSToGeodetic(mgrs.c_str(), lat, lon, height);

        ostringstream oss;
        oss << degrees(lat);
        ret.lat = oss.str();
        oss.str("");
        oss << degrees(lon);
        ret.lon = oss.str();
    }
    catch( ... )
    {
        warn << startl  << "latLonStrings: conversion failed: mgrs = " << mgrs << endl;
    }

    return ret;
}

TEST(geotrans_convertMGRSToGeodetic)
{
    char *mgrs = "31NAA6602100000"; // the point where the Greenwich meridian crosses the equator
//    char *mgrs = "31NAA6602100001"; // one metre north

    double lat = -1e42, lon = -1e42, height = -1e42;

    try
    {
        convertMGRSToGeodetic(mgrs, lat, lon, height);
//        cout << "convertMGRSToGeodetic: Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        double prec = 1e-7; // 1e-8 not only fails but causes a crash at ProfileRangeVec unit test !!!
        CHECK_CLOSE(0, lat, prec);
        CHECK_CLOSE(0, lon, prec);
        CHECK_CLOSE(0, height, prec);
    }
    catch (...)
    {
        CHECK(false);
    }
}

// This function is written by DGW
void convertGeodeticToMGRS(
        double lat,
        double lon,
        double height,
        char* &mgrsString,
        Precision::Enum& precision)
{
    try
    {
        MSP::CCS::CoordinateSystemParameters mgrsParameters(MSP::CCS::CoordinateType::militaryGridReferenceSystem);
        MSP::CCS::GeodeticParameters geodeticParameters(MSP::CCS::CoordinateType::geodetic, MSP::CCS::HeightType::ellipsoidHeight);
        MSP::CCS::CoordinateConversionService ccs("WGE", &geodeticParameters, "WGE", &mgrsParameters);
        MSP::CCS::Accuracy sourceAccuracy;
        MSP::CCS::Accuracy targetAccuracy;

        MSP::CCS::GeodeticCoordinates sourceCoordinates(MSP::CCS::CoordinateType::geodetic, lon, lat, height);
        MSP::CCS::MGRSorUSNGCoordinates targetCoordinates;

        ccs.convertSourceToTarget( &sourceCoordinates, &sourceAccuracy, targetCoordinates, targetAccuracy );
        mgrsString = targetCoordinates.MGRSString();
        precision = targetCoordinates.precision();
    }
    catch(MSP::CCS::CoordinateConversionException e)
    {
        warn << startl  << "CoordinateConversionException: " << e.getMessage() << endl;
        throw e;
    }
}

TEST(geotrans_convertGeodeticToMGRS)
{
    char *mgrs;
    Precision::Enum precision;

    try
    {
        convertGeodeticToMGRS(0, 0, 0, mgrs, precision);
//        cout << "convertGeodeticToMGRS: MGRS = " << mgrs << " Precision = " << precision << endl;

        string expected("31NAA6602100000");
        string actual(mgrs);
        CHECK_EQUAL(expected, actual);
    }
    catch (...)
    {
        CHECK(false);
    }
}

TEST(geotrans_convertGeodeticToMGRSAndBackToGeodetic)
{
    char *mgrs;
    Precision::Enum precision;

    try
    {
        convertGeodeticToMGRS(0, 0, 0, mgrs, precision);
//        cout << "convertGeodeticToMGRSAndBackToGeodetic: MGRS = " << mgrs << " Precision = " << precision << endl;

        string expected("31NAA6602100000");
        string actual(mgrs);
        CHECK_EQUAL(expected, actual);

        double lat = -1e42, lon = -1e42, height = -1e42;

        convertMGRSToGeodetic(mgrs, lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "convertGeodeticToMGRSAndBackToGeodetic: Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        double prec = 1e-5; // ~ one metre in degrees
        CHECK_CLOSE(0, lat, prec);
        CHECK_CLOSE(0, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        // Ross-on-wye
        convertGeodeticToMGRS(radians(51.9145), radians(-2.5824), 0, mgrs, precision);
//        cout << "convertGeodeticToMGRSAndBackToGeodetic: MGRS = " << mgrs << " Precision = " << precision << endl;

        string expected2("30UWC2872351611");
        string actual2(mgrs);
        CHECK_EQUAL(expected2, actual2);

        lat = -1e42; lon = -1e42; height = -1e42;

        convertMGRSToGeodetic(mgrs, lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "convertGeodeticToMGRSAndBackToGeodetic: Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.9145, lat, prec);
        CHECK_CLOSE(-2.5824, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);
    }
    catch (...)
    {
        CHECK(false);
    }
}

TEST(geotrans_truncations)
{
    // Results confirmed by GEOTRANS app
    //
    //        geotrans_truncations: 30UWC2872351611 : Lat = 51.9145, Lon = -2.5824,  Height = 0
    //        geotrans_truncations: 30UWC28725161   : Lat = 51.9145, Lon = -2.58244, Height = 0
    //        geotrans_truncations: 30UWC287516     : Lat = 51.9144, Lon = -2.58274, Height = 0
    //        geotrans_truncations: 30UWC2851       : Lat = 51.909,  Lon = -2.59296, Height = 0
    //        geotrans_truncations: 30UWC25         : Lat = 51.9004, Lon = -2.70931, Height = 0
    //        geotrans_truncations: 30UWC           : Lat = 51.4512, Lon = -3,       Height = 0

    try
    {
        double lat = -1e42, lon = -1e42, height = -1e42;
        double prec = 1e-4; // ~ 10 meters in degrees

        // Ross-on-wye
        string Ross_on_Wye("30UWC2872351611"); // MGRS

        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.9145, lat, prec);
        CHECK_CLOSE(-2.5824, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        Ross_on_Wye = string("30UWC") + string("2872") + string("5161");
        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.9145, lat, prec);
        CHECK_CLOSE(-2.5824, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        Ross_on_Wye = string("30UWC") + string("287") + string("516");
        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.9144, lat, prec);
        CHECK_CLOSE(-2.5827, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        Ross_on_Wye = string("30UWC") + string("28") + string("51");
        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.909, lat, prec);
        CHECK_CLOSE(-2.59296, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        Ross_on_Wye = string("30UWC") + string("2") + string("5");
        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.9004, lat, prec);
        CHECK_CLOSE(-2.70931, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

        Ross_on_Wye = "30UWC";
        convertMGRSToGeodetic(Ross_on_Wye.c_str(), lat, lon, height);
        lat = degrees(lat);
        lon = degrees(lon);
//        cout << "geotrans_truncations: " << Ross_on_Wye << " : Lat = " << lat << ", Lon = " << lon << ", Height = " << height << endl;
        CHECK_CLOSE(51.4512, lat, prec);
        CHECK_CLOSE(-3, lon, prec);
        CHECK_CLOSE(0, height, 1e-2);

    }
    catch (...)
    {
        CHECK(false);
    }
}

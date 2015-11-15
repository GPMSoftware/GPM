/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * Loci.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: gareth
 */

#include "loci.h"

// Names of all loci we know about
char locus_name[][name_size] =
{
  "D8S1179" //  D8S1179 = 0,
, "D21S11"  //  D21S11  = 1,
, "D7S820"  //  D7S820  = 2,
, "CSF1PO"  //  CSF1PO  = 3,
, "D3S1358" //  D3S1358 = 4,
, "TH01"    //  TH01    = 5,
, "D13S317" //  D13S317 = 6,
, "D16S539" //  D16S539 = 7,
, "D2S1338" //  D2S1338 = 8,
, "D19S433" //  D19S433 = 9,
, "vWA"     //  vWA     = 10,
, "TPOX"    //  TPOX    = 11,
, "D18S51"  //  D18S51  = 12,
, "AMEL"    //  AMEL    = 13
, "D5S818"  //  D5S818  = 14,
, "FGA"     //  FGA     = 15,
, "Penta D" //          = 16,
, "Penta E" //          = 17,
// new Jan 2015
, "D1S1656" //          = 18,
, "D12S391" //          = 19,
, "D10S1248"//          = 20,
, "D22S1045"//          = 21,
, "D2S441"  //          = 22
, "SE33"    //          = 23,
, "DYS391"  //          = 24,  // a new XY locus, This is in the GlobalFiler kit

, "AMXXXXXXEL" //       = 13   // alternative name
, "THO1"    //          = 5    // alternative name
};

const int locus_name_size = ( sizeof(locus_name) / name_size );

// Index of locus_name array (corresponds with symbols in enum Locus)
// NB this is also the order they go in the database!
int locus_index[] =
{
  0
, 1
, 2
, 3
, 4
, 5
, 6
, 7
, 8
, 9
, 10
, 11
, 12
, 13
, 14
, 15
, 16
, 17
, 18
, 19
, 20
, 21
, 22
, 23
, 24
// alternative names
, 13
, 5
};

// names of files in population database (could we not just uppercase and add .txt?)
char popdata_file[][filename_size] =
{
  "D8S1179.txt" //  D8S1179 = 0,
, "D21S11.txt"  //  D21S11,
, "D7S820.txt"  //  D7S820,
, "CSF1PO.txt"  //  CSF1PO,
, "D3S1358.txt" //  D3S1358,
, "TH01.txt"    //  THO1,     // NB filename correct
, "D13S317.txt" //  D13S317,
, "D16S539.txt" //  D16S539,
, "D2S1338.txt" //  D2S1338,
, "D19S433.txt" //  D19S433,
, "VWA.txt"     //  vWA,      // TODO filename misspelled ?
, "TPOX.txt"    //  TPOX,
, "D18S51.txt"  //  D18S51,
, "AMEL.txt"    //  AMEL,
, "D5S818.txt"  //  D5S818,
, "FGA.txt"     //  FGA,
, "Penta_D.txt" //  Penta D
, "Penta_E.txt" //  Penta E
, "D1S1656.txt"
, "D12S391.txt"
, "D10S1248.txt"
, "D2S441.txt"
, "D22S1045.txt"
, "SE33.txt"
, "DYS391.txt"
};

Locus identifiler_loci[] =
{
  D8S1179
, D21S11
, D7S820
, CSF1PO
, D3S1358
, TH01
, D13S317
, D16S539
, D2S1338
, D19S433
, VWA
, TPOX
, D18S51
, AMEL
, D5S818
, FGA
};

const int num_identifiler_loci = sizeof(identifiler_loci) / sizeof(Locus);

Locus powerplex16_loci[] =
{
  D3S1358 //
, TH01    //
, D21S11  //
, D18S51  //
, PENTA_E
, D5S818  //
, D13S317 //
, D7S820  //
, D16S539 //
, CSF1PO  //
, PENTA_D
, AMEL
, VWA     //
, D8S1179 //
, TPOX    //
, FGA     //
};

const int num_powerplex16_loci = sizeof(powerplex16_loci) / sizeof(Locus);

Locus powerplexesi17_loci[] =
{
  AMEL
, D3S1358
, D19S433
, D2S1338
, D22S1045
, D16S539
, D18S51
, D1S1656
, D10S1248
, D2S441
, TH01
, VWA
, D21S11
, D12S391
, D8S1179
, FGA
, SE33
};

const int num_powerplexesi17_loci = sizeof(powerplexesi17_loci) / sizeof(Locus);

// NGM Select
Locus ngm_loci[] =
{
    D10S1248,
    VWA,
    D16S539,
    D2S1338,
    AMEL,
    D8S1179,
    D21S11,
    D18S51,
    D22S1045,
    D19S433,
    TH01,
    FGA,
    D2S441,
    D3S1358,
    D1S1656,
    D12S391,
    SE33,
};

const int num_ngm_loci = sizeof(ngm_loci) / sizeof(Locus);

Locus globalfiler_loci[] =
{
    D3S1358
,   VWA
,   D16S539
,   CSF1PO
,   TPOX
//,   Yindel // = 1, 2
,   AMEL
,   D8S1179
,   D21S11
,   D18S51
,   DYS391
,   D2S441
,   D19S433
,   TH01
,   FGA
,   D22S1045
,   D5S818
,   D13S317
,   D7S820
,   SE33
,   D10S1248
,   D1S1656
,   D12S391
,   D2S1338
};

const int num_globalfiler_loci = sizeof(globalfiler_loci) / sizeof(Locus);

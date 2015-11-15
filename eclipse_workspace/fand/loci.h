/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * loci.h
 *
 *  Created on: Jan 6, 2015
 *      Author: gareth
 */

#ifndef LOCI_H_
#define LOCI_H_

// known loci (all kits)
enum Locus
{
    D8S1179 = 0,
    D21S11,
    D7S820,
    CSF1PO,
    D3S1358,
    TH01,
    D13S317,
    D16S539,
    D2S1338,
    D19S433,
    VWA,    // linked with D12
    TPOX,
    D18S51,
    AMEL,   // XY
    D5S818,
    FGA,
    PENTA_D,
    PENTA_E,
    // new Jan 2015
    D1S1656,
    D12S391,    // linked with D12
    D10S1248,
    D22S1045,
    D2S441,
    SE33,
    DYS391,     // Y chromosome
    num_loci,
    locus_none = -1
};

static const int filename_size = 20;
extern char popdata_file[][filename_size];

const int name_size = 15;            // length of locus name
extern char locus_name[][name_size]; // All loci we know about
extern const int locus_name_size;    // = ( sizeof(locus_name) / name_size )

extern int locus_index[];

extern Locus identifiler_loci[];
extern const int num_identifiler_loci; // = sizeof(identifiler_loci);

extern Locus powerplex16_loci[];
extern const int num_powerplex16_loci; // = sizeof(powerplex16_loci);

extern Locus powerplexesi17_loci[];
extern const int num_powerplexesi17_loci; // = sizeof(powerplexesi17_loci);

extern Locus ngm_loci[];
extern const int num_ngm_loci; // = sizeof(powerplex16_loci);

extern Locus globalfiler_loci[];
extern const int num_globalfiler_loci; // = sizeof(powerplex16_loci);

#endif /* LOCI_H_ */

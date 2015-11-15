/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
//============================================================================
// Name        : FAND1.cpp
// Author      : Gareth Williams
// Description : Run the unit tests and acceptance tests in the fand library
//============================================================================

// There are three sets of acceptance tests that may be run from eclipse as
// separate configurations:
// 1) Integration tests
// 2) BULK tests
// 3) BULK2 tests

#include "fand/dbmatch.h"
#include "fand/Profile.h"
#include "fand/MessageStream.h"
#include "fand/HostManager.h"
#include "fand/util.h"

#include <UnitTest++/UnitTest++.h>
#include <UnitTest++/TestReporterStdout.h>


using namespace std;

//#define SPLIT split_parallel // Run GPU solutions on as many GPUs as are available
#define SPLIT split_none // Run GPU solutions on one GPU

static vector<string> failures;
static int npassed=0, ntested=0;

static char *db_name = "fand";  /* database name (default=none) */

// run a regression/acceptance test
//
// matches a crime file against a reference file (or against itself),
// write the results to a results file and compares this with an expected
// results file
//
int regtest(
		string const &test_name,
		string const &crime_source, // file or dataset
		string const &ref_source,    // file or dataset
		bool          ref_is_crime,
		string const &results_file,
		string const &expected_file,
		double        lr_threshold,
		Split         cuda_split,
		string const &relatives,
        bool          n_cuda_accel,
        bool          n2_cuda_accel,
        string        prec = "1e-5",
        bool          from_db = false,
        SubPopModel const &spm = SubPopModel())
{
	failures.push_back(test_name); // remove if test passes
	ntested++;

	ProfileRange crime_db, ref_db_x;

    // read crime
	if (! readProfiles(crime_db, crime_source, from_db ? db_name : "") )
    {
    	cout << test_name << ": Failed to read " << crime_source << endl;
    	return -1;
    }

    if (!ref_is_crime)
    {
    	// read reference
        if (! readProfiles(ref_db_x, ref_source, from_db ? db_name : "") )
        {
        	cout << test_name << ": Failed to read " << ref_source << endl;
        	return -1;
        }
    }

    ProfileRange &ref_db = ref_is_crime ? crime_db : ref_db_x;

	if (crime_db.size() == 0 || ref_db.size() == 0)
	{
		cout << test_name << ": nothing to match. Goodbye." << endl;
		return -1;
	}

	//
	// do match
	//
	NMmatchResults matches(0);
	vector<MatchType> match_types;
	matchTypesFromString(match_types, split(relatives,   ":"));

	Timer t;

	cout << "Running " << test_name << " ..." << endl;

	do_gmatch(
			matches,
			std::cout,
			cuda_split,
			crime_db,
			ref_db,
			ref_is_crime,
			0,                 // crime_delta
			0,                 // ref_delta
			match_types,
			spm,               // subpopulation
			MutModel(),        // mutation
			lr_threshold,
	        n_cuda_accel,
	        n2_cuda_accel);

//	t.stop();

	cout << " test completed. Time elapsed = " << t.stop() << " copying results ... " << endl;

    std::map<int, std::string> c_ids, r_ids;

    // first create the indices with "" for the ids
    NMmatchResults::const_iterator it;
    for(it = matches.begin(); it != matches.end(); ++it)
    {
        int c_index = it->first.first;
        int r_index = it->first.second;

        c_ids[c_index] = "";
        r_ids[r_index] = "";
    }

    // Now read the profiles in chunks and fill in the correct values.
    // For efficiency, read 1000 profiles at a time
	static int NPROF = 1000;
	getProfileNames(c_ids, crime_db, NPROF);
	getProfileNames(r_ids, ref_db, NPROF);

	cout << " results copied. Time elapsed = " << t.stop() << " writing results to file ... " << endl;

	// write results
    std::ofstream ofs(results_file.c_str());
    if (! ofs.good()){
        cout << "error opening file: " << results_file << endl;
        return -1;
    }

	MatchResults results(match_types, matches, c_ids, r_ids, lr_threshold);
	outputMatches(ofs, results);
	ofs.close();

	cout << " results written. Time elapsed = " << t.stop() << " checking results ... " << endl;

	// compare results
	string command = "./numdiff.pl " + results_file + " " + expected_file + " " + prec;

	int ret = system(command.c_str());
	if (ret == 0)
	{
		cout << test_name << " passed. Time elapsed = " << t.stop() << endl;
		npassed++;
		failures.pop_back();
	}
	else
	{
		cout << test_name << " failed with " << ret << endl;
	}

	return ret;
}

char *progname = 0;

extern int SIMPLE_OPTIMIZATION;

int main(int argc, char *argv[])
{
    progname = argv[0];

    SIMPLE_OPTIMIZATION = 0; // turn off simple_optimization for testing (try it both ways)

    // see how much memory we have got to play with
    MemMeter m;
    cout << "this process is now using " << m << endl;
    theHostManager().init();
    cout << "available    = " << theHostManager().available()           << " kB" << endl;
    cout << "total_app    = " << theHostManager().getLimits().total_app << " kB" << endl;
    cout << "Cached limit = " << Cached::getLimit() << endl;

//    if (progname == 0) // Dont't run unit tests, but still compile.
//                          (this test is False but the compiler does not know it!)
//    {///////////////////////

	cout << "Running Unit Tests ..." << endl;
	MessageStream::enable(false); // suppress error/warning/etc during unit tests

#if 1
	// Run all tests
	if (int test_result = UnitTest::RunAllTests())
#else
	// Run a named test suite
	UnitTest::TestReporterStdout reporter;
	UnitTest::TestRunner runner(reporter);
	if ( int test_result = runner.RunTestsIf(
							UnitTest::Test::GetTestList(),
							"KitsSuite",
							UnitTest::True(), 0))
#endif
	{
		cout << "Unit tests failed - quitting" << endl;
		return test_result;
	}



	MessageStream::enable(true);
	cout << "Unit tests completed successfully." << endl << endl;

	if (getBoolEnv("UNIT_TESTS_ONLY", false))
	{
		return 0;
	}
//    }///////////////////

	int ret = 0;

	int NSOAK = getIntEnv("NSOAK", 1);

	for (int i = 1; i <= NSOAK; ++i) // SOAK TEST
	{

	if (NSOAK > 1)
	{
		cout << "SOAK TEST: " << i << "/" << NSOAK << endl;
	}

	bool MYSQL_TESTS = getBoolEnv("MYSQL_TESTS", false);

	// set rare allele policy to "ignore"
	PopulationData pop = populationData();
    pop.setPolicy(ignore);
    popSet(pop);

	if (! getBoolEnv("BULK2_TESTS", false))
	{

	cout << "Running Regression Tests ..." << endl;



    if ( getBoolEnv("BULK_TESTS", false))
    {

	// In this series of tests, populations of relatives generated by
	// the simstudy program are tested.
	//
	// GPU results are compared against CPU results

	// fam16 (96 profiles) against itself (read from file, CPU, HW)
	ret |= regtest("fam16a",
				   "FAND_test/fam16.idf.csv",
				   "FAND_test/fam16.idf.csv",
				   false,
				   "fam16a.res.csv",
				   "FAND_test/fam16_cpu.res.csv",
				   1e2,
				   SPLIT,
				   "IDENT:D1:SIB:D2:D03",
				   false,
				   false);

	// fam16 (96 profiles) against itself (read from file, GPU, HW)
	ret |= regtest("fam16b",
				   "FAND_test/fam16.idf.csv",
				   "FAND_test/fam16.idf.csv",
				   false,
				   "fam16b.res.csv",
				   "FAND_test/fam16_gpu.res.csv",
				   1e2,
				   SPLIT,
				   "IDENT:D1:SIB:D2:D03",
				   true,
				   true);

	// fam16 (96 profiles) against itself (read from file, CPU, B11-theta = 0.02)
	ret |= regtest("fam16c",
				   "FAND_test/fam16.idf.csv",
				   "FAND_test/fam16.idf.csv",
				   false,
				   "fam16c.res.csv",
				   "FAND_test/fam16_cpu_theta002.res.csv",
				   1e2,
				   SPLIT,
				   "IDENT:D1:SIB:D2:D03",
				   false,
				   false,
				   "1e-5",                                        // precision
				   false,                                         // read from MySQL database
				   SubPopModel(SubPopModel::B11, 0.02));          // Subpopulation model

	// fam16 (96 profiles) against itself (read from file, GPU, B11-theta = 0.02)
	ret |= regtest("fam16d",
				   "FAND_test/fam16.idf.csv",
				   "FAND_test/fam16.idf.csv",
				   false,
				   "fam16d.res.csv",
				   "FAND_test/fam16_gpu_theta002.res.csv",
				   1e2,
				   SPLIT,
				   "IDENT:D1:SIB:D2:D03",
				   true,
				   true,
				   "1e-5",                                    // precision
				   false,                                     // read from MySQL database
				   SubPopModel(SubPopModel::B11, 0.02));      // Subpopulation model


	// fam166 (996 profiles) against itself (read from file, GPU)
	ret |= regtest("fam166b",
				   "FAND_test/fam166.idf.csv",
				   "",
				   true,
				   "fam166b.res.csv",
				   "FAND_test/fam166_gpu.res.csv",
				   1e4,
				   SPLIT,
				   "IDENT:D1:SIB:D2:D03",
				   true,
				   true);

	// In this series of tests crime and reference datasets generated by nmatch1 are compared.
	// These datasets contain copies of a single profile at fixed positions. There are no actual relatives.
	// Large dataset matches - GPU only.
	// The data is read from file, and from the MySQL database
	//

	// Test of generated data, 1000 x 1000 profiles (GPU, read from file)
	// (files in identifiler format)
	ret |= regtest("idf1",
			       "FAND_test/crime_test_1000.idf.csv",
			       "FAND_test/ref_test_1000.idf.csv",
			       false,
			       "idf1.res.csv",
			       "FAND_test/crime_ref_1000_idf_gpu_ignore.res.csv",
			       1e4,
			       SPLIT,
			       "IDENT:D1:SIB:D2:D03",
			       true,
			       true);

	// crime 1000 against itself (GPU, read from file)
	ret |= regtest("idf2",
			       "FAND_test/crime_test_1000.idf.csv",
			       "",
			       true,
			       "idf2.res.csv",
			       "FAND_test/crime_crime_1000_idf_gpu_ignore.res.csv",
			       1e4,
			       SPLIT,
			       "IDENT:D1:SIB:D2:D03",
			       true,
			       true);

	// same thing (read from database) This can take along time if the database is large
	if (MYSQL_TESTS)
	{
	ret |= regtest("idf2db",
			       "crime_test_1000",
			       "",
			       true,
			       "idf2.res.csv",
			       "FAND_test/crime_crime_1000_idf_gpu_ignore.res.csv",
			       1e4,
			       SPLIT,
			       "IDENT:D1:SIB:D2:D03",
			       true,
			       true,
			       "1e-3",                                                   // precision
			       true);                                                    // read from MySQL database
	}
#if 0
	// some even larger tests that can be enabled when required

  //	 repeat with 10,000 x 10,000
	ret |= regtest("idf3",
			       "FAND_test/crime_test_10000.idf.csv",
			       "FAND_test/ref_test_10000.idf.csv",
			       false,
			       "idf3.res.csv",
			       "FAND_test/crime_ref_10000_idf_gpu_ignore.res.csv",
			       1e6,
			       SPLIT,
			       "IDENT:D1:SIB:D2:D03",
			       true,
			       true);


	// 10,000 against itself
	ret |= regtest("idf4",
			       "FAND_test/crime_test_10000.idf.csv",
			       "",
			       true,
			       "idf4.res.csv",
			       "FAND_test/crime_crime_10000_idf_gpu_ignore.res.csv",
			       1e6,
			       SPLIT,
			       "IDENT:D1:SIB:D2:D03",
			       true,
			       true);

// #else


    // repeat with 100,000 x 1,000
    ret |= regtest("idf5",
                   "FAND_test/ref_test_1e5.idf.csv",
                   "FAND_test/crime_test_1000.idf.csv",
                   false,
                   "idf5.res.csv",
                   "FAND_test/ref1e5_crime1000_ident_idf_gpu_ignore.res.csv",
                   1e6,
                   SPLIT,
                   "IDENT",
                   true,
                   true);

    // repeat with 200,000 x 1,000
    ret |= regtest("idf6",
                   "FAND_test/ref_test_2e5.idf.csv",
                   "FAND_test/crime_test_1000.idf.csv",
                   false,
                   "idf6.res.csv",
                   "FAND_test/ref2e5_crime1000_ident_idf_gpu_ignore.res.csv",
                   1e6,
                   SPLIT,
                   "IDENT",
                   true,
                   true);

    // repeat with 500,000 x 1,000
    ret |= regtest("idf7",
                   "FAND_test/ref_test_5e5.idf.csv",
                   "FAND_test/crime_test_1000.idf.csv",
                   false,
                   "idf7.res.csv",
                   "FAND_test/ref5e5_crime1000_ident_idf_gpu_ignore.res.csv",
                   1e6,
                   SPLIT,
                   "IDENT",
                   true,
                   true);

//    // repeat with 1,000,000 x 1,000
//    ret |= regtest("idf8",
//                   "FAND_test/ref_test_1e6.idf.csv",
//                   "FAND_test/crime_test_1000.idf.csv",
//                   false,
//                   "idf8.res.csv",
//                   "FAND_test/ref1e6_crime1000_ident_idf_gpu_ignore.res.csv",
//                   1e6,
//                   SPLIT,
//                   "IDENT",
//                   true,
//                   true);
//
//    // repeat with 2,000,000 x 1,000
//    ret |= regtest("idf9",
//                   "FAND_test/ref_test_2e6.idf.csv",
//                   "FAND_test/crime_test_1000.idf.csv",
//                   false,
//                   "idf8.res.csv",
//                   "FAND_test/ref2e6_crime1000_ident_idf_gpu_ignore.res.csv",
//                   1e6,
//                   SPLIT,
//                   "IDENT",
//                   true,
//                   true);

#endif

//	// "other CSV" format. NB background is not merged!
//	ret |= regtest("oth1", "FAND_test/crime_test_10.csv", "FAND_test/ref_test_10.csv",  false,
//			       "oth1.res.csv", "FAND_test/crime_ref_10_idf_cpu_ignore.res.csv",
//			       100, SPLIT, "IDENT:D1:SIB:D2", false, false);


} else { // BULK_TESTS == false

	// Main acceptance tests: check identity and familial LRs against independently calculated results

	// Full profile tests, HWE
	//
	// This series of tests uses the generated files sibtest1.b11.csv and sibtest2.b11.csv
	// sibtest1.b11.csv contains four profiles MUM-FP, DAD-FP, SIB1-FP, SIB2-FP.
	// MUM and DAD and randomly generated, and SIB1 and SIB2 are randomly generated children of MUM and DAD
	// sibtest2.b11.csv is a copy of sibtest1.b11.csv with the profiles renamed
	// MUM-1-FP, DAD-1-FP, SIB1-1-FP, SIB2-1-FP.
	//
	// The results are calculated using formulae in sibtest_formula.ods
	//
	// NB for the following tests environment variable POPDATA should be set to "Example"
	//

	// NM match for IDENT, D1, SIB, D2, on CPU
    ret |= regtest("reg1",                                                 // test name
    		       "FAND_test/sibtest1.b11.csv",                            // crime file
    		       "FAND_test/sibtest2.b11.csv",                            // reference file (optional)
    		       false,                                                  // use crime as reference
			       "reg1.res.csv",                                         // results file
			       "FAND_test/sibtest1_sibtest2_cpu_ignore.res.csv",        // expected results file
			       100,                                                    // LR threshold
			       SPLIT,                                                  // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                      // matches
			       false,                                                  // CUDA accelerate N matches
			       false);                                                 // CUDA accelerate N2 matches

	// As reg1 using the GEN specifications of the relationships
    ret |= regtest("reg1g",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg1g.res.csv",                                          // results file
			       "FAND_test/sibtest1_sibtest2_cpu_ignore_g.res.csv",        // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "GEN[1,0,0,1]:GEN[0.5,0.5,0,0]:GEN[0.5,0,0,0.5]:GEN[0.25,0.25,0,0]", // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// As reg1 using the INV specifications (these are symmetric relationships)
    ret |= regtest("reg1i",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg1i.res.csv",                                          // results file
			       "FAND_test/sibtest1_sibtest2_cpu_ignore_i.res.csv",        // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "INV[1,0,0,1]:INV[0.5,0.5,0,0]:INV[0.5,0,0,0.5]:INV[0.25,0.25,0,0]", // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

    // As reg1, on GPUs
    ret |= regtest("reg2",                                                   // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2.res.csv",                                           // results file
			       "FAND_test/sibtest1_sibtest2_gpu_ignore.res.csv",          // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

    // As reg2 using the GEN specifications of the relationships
    ret |= regtest("reg2g",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2g.res.csv",                                          // results file
			       "FAND_test/sibtest1_sibtest2_gpu_ignore_g.res.csv",        // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "GEN[1,0,0,1]:GEN[0.5,0.5,0,0]:GEN[0.5,0,0,0.5]:GEN[0.25,0.25,0,0]", // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

    // As reg2 using the INV specifications (these are symmetric relationships)
    ret |= regtest("reg2i",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2i.res.csv",                                          // results file
			       "FAND_test/sibtest1_sibtest2_gpu_ignore_i.res.csv",        // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "INV[1,0,0,1]:INV[0.5,0.5,0,0]:INV[0.5,0,0,0.5]:INV[0.25,0.25,0,0]", // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

    // Testing the D12 asymmetric relationship: D12, the GEN specification of D12, and the inverse of D12, on CPU
    ret |= regtest("reg1as",                                                 // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg1as.res.csv",                                         // results file
			       "FAND_test/sibtest1_sibtest2_cpu_ignore_d12.res.csv",      // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "D12:GEN[1/2,1/2,1/4,1/4]:INV[1/2,1/2,1/4,1/4]",          // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

    // As reg1as, on GPU
    ret |= regtest("reg2as",                                                 // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2as.res.csv",                                         // results file
			       "FAND_test/sibtest1_sibtest2_gpu_ignore_d12.res.csv",      // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "D12:GEN[1/2,1/2,1/4,1/4]:INV[1/2,1/2,1/4,1/4]",          // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

    // N2 match, sibtest1.b11.csv against itself, CPU
    ret |= regtest("reg1a",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "reg1a.res.csv",                                          // results file
			       "FAND_test/sibtest1_itself_cpu_ignore.res.csv",            // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// As reg1a, GPUs
    ret |= regtest("reg2a",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "reg2a.res.csv",                                          // results file
			       "FAND_test/sibtest1_itself_gpu_ignore.res.csv",            // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

    // N match, sib1-fp.b11.csv (single profile) against sibtest2.b11.csv, CPU
    ret |= regtest("reg1b",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg1b.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sibtest2_cpu_ignore.res.csv",           // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// As reg1b, GPUs
    ret |= regtest("reg2b",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2b.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sibtest2_gpu_ignore.res.csv",           // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// profile/profile match (CPU only)
    ret |= regtest("reg1c",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sib2-a-fp.b11.csv",                             // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg1c.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sib2-a-fp_cpu_ignore.res.csv",          // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// As reg1c, GPUs specified, but will run on the CPU anyway
    ret |= regtest("reg2c",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sib2-a-fp.b11.csv",                             // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg2c.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sib2-a-fp_cpu_ignore.res.csv",          // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// set rare allele policy to "add as rare", and repeat the tests for "IDENT:D1:SIB:D2"
    PopulationData pop = populationData();
    pop.setPolicy(add_as_rare);
    popSet(pop);

	// NM, CPU
    ret |= regtest("reg3",                                                   // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg3.res.csv",                                           // results file
			       "FAND_test/sibtest1_sibtest2_cpu_add.res.csv",             // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

    // NM, GPU
    ret |= regtest("reg4",                                                   // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg4.res.csv",                                           // results file
			       "FAND_test/sibtest1_sibtest2_gpu_add.res.csv",             // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// N2, CPU
	ret |= regtest("reg3a",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "reg3a.res.csv",                                          // results file
			       "FAND_test/sibtest1_itself_cpu_add.res.csv",               // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// N2, GPU
	ret |= regtest("reg4a",                                                  // test name
    		       "FAND_test/sibtest1.b11.csv",                              // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "reg4a.res.csv",                                          // results file
			       "FAND_test/sibtest1_itself_gpu_add.res.csv",               // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// N, CPU
	ret |= regtest("reg3b",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg3b.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sibtest2_cpu_add.res.csv",              // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// N, GPU
	ret |= regtest("reg4b",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sibtest2.b11.csv",                              // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg4b.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sibtest2_gpu_add.res.csv",              // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// profile/profile match (CPU only)
	ret |= regtest("reg3c",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sib2-a-fp.b11.csv",                             // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg3c.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sib2-a-fp_cpu_add.res.csv",             // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       false,                                                    // CUDA accelerate N matches
			       false);                                                   // CUDA accelerate N2 matches

	// GPUs specified, but will run on the CPU anyway
	ret |= regtest("reg4c",                                                  // test name
    		       "FAND_test/sib1-fp.b11.csv",                               // crime file
    		       "FAND_test/sib2-a-fp.b11.csv",                             // reference file (optional)
    		       false,                                                    // use crime as reference
			       "reg4c.res.csv",                                          // results file
			       "FAND_test/sib1-fp_sib2-a-fp_cpu_add.res.csv",             // expected results file
			       100,                                                      // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2",                                        // matches
			       true,                                                     // CUDA accelerate N matches
			       true);                                                    // CUDA accelerate N2 matches

	// Full profile tests, including subpopulation models, compared with
	// results calculated using the Balding-Nichols formulae
	// (The results are calculated in match.cpp in the unit test B11_profile
	// and may be output by setting the variable print to true).

	// use test population database
	string poppath = getenv("POPDATA");
	string testpop = path(poppath) + "Test";

    pop = populationData();
    Assert2(pop.read(testpop), "Population database 'Test' can't be read");
    popSet(pop);

//	cout << "numLoci = " << populationData().numLoci() << endl;
//	cout << populationData().getLocusInfo() << endl;

	// CPU, theta = 0
	ret |= regtest("prof_a",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_a.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_000.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0));                        // Subpopulation model

	// CPU, theta = 0.01
	ret |= regtest("prof_b",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_b.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_001.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.01));                     // Subpopulation model

	// CPU, theta = 0.05
	ret |= regtest("prof_c",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_c.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_005.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.05));                     // Subpopulation model

	// GPU, theta = 0
	ret |= regtest("prof_d",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_d.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_000.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0));                        // Subpopulation model

	// GPU, theta = 0.01
	ret |= regtest("prof_e",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_e.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_001.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.01));                     // Subpopulation model

	// GPU, theta = 0.05
	ret |= regtest("prof_f",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_f.res.csv",                                         // results file
			       "FAND_test/testset_cpu_B11_005.res.csv",                   // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.05));                     // Subpopulation model

	// GPU, theta = 0.05, NRC4_4 for IDENT, B11 for relatives
	ret |= regtest("prof_g",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_g.res.csv",                                         // results file
			       "FAND_test/testset_cpu_NRC44_005.res.csv",                 // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::NRC4_4, 0.05));                  // 4_4 for IDENT only, B11 otherwise

	// CPU, theta = 0.05
	ret |= regtest("prof_h",                                                 // test name
    		       "FAND_test/b11_testset.b11.csv",                           // crime file
    		       "",                                                       // reference file (optional)
    		       true,                                                     // use crime as reference
			       "prof_h.res.csv",                                         // results file
			       "FAND_test/testset_cpu_NRC44_005.res.csv",                 // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::NRC4_4, 0.05));                  // 4_4 for IDENT only, B11 otherwise

	// CPU, theta = 0.05, 1 vs n profiles
	ret |= regtest("prof_i",                                                 // test name
				   "FAND_test/b11_testset_c.b11.csv",                         // crime file
				   "FAND_test/b11_testset_r.b11.csv",                         // reference file (optional)
    		       false,                                                    // use crime as reference
			       "prof_i.res.csv",                                         // results file
			       "FAND_test/testset_cpu_1n_B11_005.res.csv",                // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.05));                     // Subpopulation model

	// CPU, theta = 0.05, 1 vs n profiles
	ret |= regtest("prof_j",                                                 // test name
				   "FAND_test/b11_testset_c.b11.csv",                         // crime file
				   "FAND_test/b11_testset_r.b11.csv",                         // reference file (optional)
    		       false,                                                    // use crime as reference
			       "prof_j.res.csv",                                         // results file
			       "FAND_test/testset_cpu_1n_NRC44_005.res.csv",              // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       false,                                                    // CUDA accelerate N matches
			       false,                                                    // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::NRC4_4, 0.05));                  // 4_4 for IDENT only, B11 otherwise

	// GPU, theta = 0.05, B11 subpopulation model
	ret |= regtest("prof_k",                                                 // test name
    		       "FAND_test/b11_testset_c.b11.csv",                         // crime file
    		       "FAND_test/b11_testset_r.b11.csv",                         // reference file (optional)
    		       false,                                                    // use crime as reference
			       "prof_k.res.csv",                                         // results file
			       "FAND_test/testset_cpu_1n_B11_005.res.csv",                // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::B11, 0.05));                     // Subpopulation model

	// GPU, theta = 0.05, B11 subpopulation model
	ret |= regtest("prof_l",                                                 // test name
    		       "FAND_test/b11_testset_c.b11.csv",                         // crime file
    		       "FAND_test/b11_testset_r.b11.csv",                         // reference file (optional)
    		       false,                                                    // use crime as reference
			       "prof_l.res.csv",                                         // results file
			       "FAND_test/testset_cpu_1n_NRC44_005.res.csv",              // expected results file
			       0,                                                        // LR threshold
			       SPLIT,                                                    // use multiple CUDA devices
			       "IDENT:D1:SIB:D2:D03",                                    // matches
			       true,                                                     // CUDA accelerate N matches
			       true,                                                     // CUDA accelerate N2 matches
                   "1e-5",                                                   // Precision for results comparison
			       false,                                                    // From database
                   SubPopModel(SubPopModel::NRC4_4, 0.05));                  // Subpopulation model

        // Test of probabilistic profiles with 'F'
	    // CPU
	    ret |= regtest("probFs",                                                 // test name
	                   "FAND_test/b11_testsetFs.b11.csv",                         // crime file
	                   "",                                                       // reference file (optional)
	                   true,                                                     // use crime as reference
	                   "probFs.res.csv",                                         // results file
	                   "FAND_test/testsetFs_cpu_n2.res.csv",                      // expected results file
	                   0,                                                        // LR threshold
	                   SPLIT,                                                    // use multiple CUDA devices
	                   "IDENT",                                                  // matches
	                   false,                                                    // CUDA accelerate N matches
	                   false);                                                   // CUDA accelerate N2 matches

        // Test of probabilistic profiles with 'F'
        // GPU
        ret |= regtest("probFsb",                                                // test name
                       "FAND_test/b11_testsetFs.b11.csv",                         // crime file
                       "",                                                       // reference file (optional)
                       true,                                                     // use crime as reference
                       "probFsb.res.csv",                                        // results file
                       "FAND_test/testsetFs_cpu_n2.res.csv",                      // expected results file
                       0,                                                        // LR threshold
                       SPLIT,                                                    // use multiple CUDA devices
                       "IDENT",                                                  // matches
                       true,                                                     // CUDA accelerate N matches
                       true);                                                    // CUDA accelerate N2 matches

        // Test of probabilistic profiles with "Allele@0.5"
        // CPU
        ret |= regtest("prob_point5s",                                           // test name
                       "FAND_test/b11_testset_point5s.b11.csv",                   // crime file
                       "",                                                       // reference file (optional)
                       true,                                                     // use crime as reference
                       "prob_point5s.res.csv",                                   // results file
                       "FAND_test/testsetPoint5s_cpu_n2.res.csv",                 // expected results file
                       0,                                                        // LR threshold
                       SPLIT,                                                    // use multiple CUDA devices
                       "IDENT",                                                  // matches
                       false,                                                    // CUDA accelerate N matches
                       false);                                                   // CUDA accelerate N2 matches

        // Test of probabilistic profiles with "Allele@0.5"
        // GPU
        ret |= regtest("prob_point5sb",                                          // test name
                       "FAND_test/b11_testset_point5s.b11.csv",                   // crime file
                       "",                                                       // reference file (optional)
                       true,                                                     // use crime as reference
                       "prob_point5sb.res.csv",                                  // results file
                       "FAND_test/testsetPoint5s_cpu_n2.res.csv",                 // expected results file
                       0,                                                        // LR threshold
                       SPLIT,                                                    // use multiple CUDA devices
                       "IDENT",                                                  // matches
                       true,                                                     // CUDA accelerate N matches
                       true);                                                    // CUDA accelerate N2 matches

} // BULK_TESTS == false

	} // BULK2_TESTS
	else
	{
		// More bulk tests, files generated by
		cout << "Running BULK2 Tests ..." << endl;

		// 500 x 1000
		ret |= regtest("B2_1e3_file",                                            // test name
	    		       "FAND_test/BULK2_TEST/CT_500.idf.csv",                     // crime file
	    		       "FAND_test/BULK2_TEST/RT_1e3.idf.csv",                     // reference file (optional)
	    		       false,                                                    // use crime as reference
				       "CT500_vs_RT1e3_T1000_v015_gmatch_s.res.csv",             // results file
				       "FAND_test/BULK2_TEST/CT500_vs_RT1e3_T1000_v015_gmatch_s.res.csv", // expected results file
				       1000,                                                     // LR threshold
				       SPLIT,                                                    // use multiple CUDA devices
				       "IDENT:D1:SIB:D2:D03",                                    // matches
				       true,                                                     // CUDA accelerate N matches
				       true);                                                    // CUDA accelerate N2 matches

		// 500 x 10,000
		ret |= regtest("B2_1e4_file",                                            // test name
	    		       "FAND_test/BULK2_TEST/CT_500.idf.csv",                     // crime file
	    		       "FAND_test/BULK2_TEST/RT_1e4.idf.csv",                     // reference file (optional)
	    		       false,                                                    // use crime as reference
				       "CT500_vs_RT1e4_T1000_v015_gmatch_s.res.csv",             // results file
				       "FAND_test/BULK2_TEST/CT500_vs_RT1e4_T1000_v015_gmatch_s.res.csv", // expected results file
				       1000,                                                     // LR threshold
				       SPLIT,                                                    // use multiple CUDA devices
				       "IDENT:D1:SIB:D2:D03",                                    // matches
				       true,                                                     // CUDA accelerate N matches
				       true);                                                    // CUDA accelerate N2 matches

		// Repeat, reading from MySQL database
		if (MYSQL_TESTS)
		{
		ret |= regtest("B2_1e3_db",                                              // test name
	    		       "CT_500",                                                 // crime database
	    		       "RT_1e3",                                                 // reference database (optional)
	    		       false,                                                    // use crime as reference
				       "CT500_vs_RT1e3_T1000_v015_gmatch_s_db.res.csv",          // results file
				       "FAND_test/BULK2_TEST/CT500_vs_RT1e3_T1000_v015_gmatch_s.res.csv", // expected results file
				       1000,                                                     // LR threshold
				       SPLIT,                                                    // use multiple CUDA devices
				       "IDENT:D1:SIB:D2:D03",                                    // matches
				       true,                                                     // CUDA accelerate N matches
				       true,                                                     // CUDA accelerate N2 matches
				       "1e-3",                                                   // precision
				       true);                                                    // read from MySQL database

		ret |= regtest("B2_1e4_db",                                              // test name
					   "CT_500",                                                 // crime database
					   "RT_1e4",                                                 // reference database (optional)
	    		       false,                                                    // use crime as reference
				       "CT500_vs_RT1e4_T1000_v015_gmatch_s_db.res.csv",          // results file
				       "FAND_test/BULK2_TEST/CT500_vs_RT1e4_T1000_v015_gmatch_s.res.csv", // expected results file
				       1000,                                                     // LR threshold
				       SPLIT,                                                    // use multiple CUDA devices
				       "IDENT:D1:SIB:D2:D03",                                    // matches
				       true,                                                     // CUDA accelerate N matches
				       true,                                                     // CUDA accelerate N2 matches
				       "1e-3",                                                   // precision
				       true);                                                    // read from MySQL database
		}
	} // end BULK2_TESTS

	cout << npassed << "/" << ntested << " regression tests passed" << endl;

	if (ret)
	{
		cout << "regression tests failed: " << failures << endl;
	}

	} // SOAK TEST

	return ret;
}

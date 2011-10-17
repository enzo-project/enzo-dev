/*
 * A simple testing framework for C++.
 * http://simplectest.sourceforge.net
 * Release 0.31
 *
 * See "readme.txt" for usage instructions.
 *
 * ---------------------
 * Configuration options: (To enable, #define these anywhere.)
 *
 *   TESTS_BREAK_ON_FAILURE
 *     If a failure occurs, break on the first failure and do not continue
 *     executing the tests. This skips the rest of the test suite (if there
 *     is any) but does not exit the program. Defaults to off.
 *
 *   TESTS_EXIT_ON_FAILURE
 *     If a failure occurs, break on the first failure and do not continue
 *     executing the tests. This skips the rest of ALL tests and suites.
 *     Defaults to off.
 *
 *   TESTS_IGNORE_EPSILON
 *     If this is defined, SET_EPSILON() calls will be silently ignored.
 *
 *   TEST_INDIVIDUAL
 *     If this is defined, the suite contained in this test file will
 *     become a self-runnable executable test. This requires exactly one
 *     suite (START_SUITE) being present in the test file.
 *
 * ---------------------
 * Simple C++ Testing Framework
 * Copyright (C) 2004 Jevon Wright
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#ifndef _SIMPLECTEST_TESTS_H
#define _SIMPLECTEST_TESTS_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <map>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <hdf5.h>

#ifdef TRANSFER
#include "preincludes.h"
#endif
#include "svn_version.def"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "communication.h"
#ifdef TRANSFER
#include "PhotonCommunication.h"
#include "gFLDProblem.h"
#endif
#undef DEFINE_STORAGE

// Define some defines that will exist when this project is in testing mode
#ifndef SIMPLECTEST
#define SIMPLECTEST
#endif

#ifndef TESTS
#define TESTS
#endif

#ifndef USE_MPI
#define MyProcessorNumber 0
#define ROOT_PROCESSOR 0
#endif

#define TESTABS(x) ((x) < 0 ? -(x) : (x))

void my_exit(int status);

int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);

// --- Main test framework ---
// The int main(void) declaration.
#define START_TESTS() \
	Eint32 main(Eint32 argc, char *argv[]) { \
		Eint32 fails = 0; \
		Eint32 tests = 0; \
		Eint32 testfails = 0; \
		Eint32 passes = 0; \
		char *name; \
		char *reason = NULL; \
		char info[80]; \
		double epsilon; \
		Eint32 continueTests = 1; \
		SET_EPSILON_DEFAULT(); \
		sprintf(info, " ");

#define END_TESTS() \
		name = name; \
		reason = reason; \
		info[0] = info[0]; \
                if (MyProcessorNumber == ROOT_PROCESSOR) { \
		    printf("\n--- Results ---\n"); \
                    printf("Tests run:%5d\n", tests);		   \
                    printf("Passes:   %5d\n", passes); \
                    printf("Failures: %5d\n", fails); \
                } \
		return fails; \
	} \
void my_exit(int status) \
{\
  CommunicationFinalize(); \
  exit(status);\
}

// Wrappers for creating independent executables for each suite
// (Use -DTEST_INDIVIDUAL to make "suites" runnable - see note above)
#ifdef TEST_INDIVIDUAL
// Create the int main(void) combo to make the suite runnable,
// and run the suite
#define START_SUITE_WRAPPER(name) \
		START_TESTS(); \
		SUITE(name); \
		END_TESTS();

#else
// Or, the wrapper definition does nothing
#define START_SUITE_WRAPPER(name)

#endif // test_individual

// --- Test suites ---
// Start a test suite (and possibly make it runnable)
// We open a brace, and close it with END_SUITE.
// Returns 1 if the test ended successfully.
#define START_SUITE(x) \
	START_SUITE_WRAPPER(x); \
	Eint32 test_suite_##x (Eint32 &fails, Eint32 &tests, Eint32 &testfails, Eint32 &passes) { \
                if (MyProcessorNumber == ROOT_PROCESSOR) { \
                   printf("[%s]\n", __FUNCTION__); \
                } \
		char *name; \
		char *reason = NULL; \
		char info[80]; \
		double epsilon; \
		SET_EPSILON_DEFAULT(); \

// End a test suite.
// Also close the brace we opened up.
// Returns 1 to indicate the test ended successfully
// (We also "use" reason/info/name so that compilation does not
// complain about unused variables)
#define END_SUITE() \
		name = name; \
		reason = reason; \
		info[0] = info[0]; \
		return 1; \
	}

// Defines the suite - use this definition in header files
#define DEFINE_SUITE(x) Eint32 test_suite_##x(Eint32 &fails, Eint32 &tests, Eint32 &testfails, Eint32 &passes);

// The actual suite execution
// We only want to run the suite of tests, if we haven't been told to cancel yet ('continueTests')
#define SUITE(x) \
		if (continueTests) { \
			Eint32 fails_##x = 0; \
			Eint32 tests_##x = 0; \
			Eint32 testfails_##x = 0; \
			Eint32 passes_##x = 0; \
			continueTests = test_suite_##x(fails_##x, tests_##x, testfails_##x, passes_##x); \
			fails += fails_##x; \
			tests += tests_##x; \
			testfails += testfails_##x; \
			passes += passes_##x; \
		}

// --- Individual tests ---
// Start a test
// We open up the brace so we have locality of reference for variables
#define START_TEST(x) \
	testfails = 0; \
        if (MyProcessorNumber == ROOT_PROCESSOR) { \
           printf("> %s...\n", x); \
        } \
	if (1) { \
		tests++; \
		name = x;

// End a test
// Close the brace we opened up (and also reset the epsilon)
#define END_TEST() \
		SET_EPSILON_DEFAULT(); \
	}

// Test description override
#define TEST(str) reason = (str);

// --- Assertion commands ---
// Asserts that a == b
#define ASSERT_EQUALS(a, b) \
	TEST(#a " is supposed to equal " #b); \
	ASSERT((a) == (b));

// Asserts that a != b
#define ASSERT_NOT_EQUALS(a, b) \
	TEST(#a " is not supposed to equal " #b); \
	ASSERT((a) != (b));

// Asserts that a >= b
#define ASSERT_GREATER_EQUAL(a, b) \
	TEST(#a " >= " #b); \
	sprintf(info, "\n\t(%f is not >= %f)", (float) (a), (float) (b)); \
	ASSERT((float) (a) >= (float) (b) || FLOAT_EQUALS(a, b));

// Asserts that a <= b
#define ASSERT_LESSTHAN_EQUAL(a, b) \
	TEST(#a " <= " #b); \
	sprintf(info, "\n\t(%f is not <= %f)", (float) (a), (float) (b)); \
	ASSERT((float) (a) <= (float) (b) || FLOAT_EQUALS(a, b));

// Asserts that a > b
#define ASSERT_GREATER(a, b) \
	TEST(#a " > " #b); \
	sprintf(info, "\n\t(%f is not > %f)", (float) (a), (float) (b)); \
	ASSERT((float) (a) > (float) (b));

// Asserts that a < b
#define ASSERT_LESSTHAN(a, b) \
	TEST(#a " < " #b); \
	sprintf(info, "\n\t(%f is not < %f)", (float) (a), (float) (b)); \
	ASSERT((float) (a) < (float) (b));

// Syntactic sugar for the macros above
#define ASSERT_EQ(a, b) ASSERT_EQUALS(a, b);
#define ASSERT_LT(a, b) ASSERT_LESSTHAN(a, b);
#define ASSERT_GT(a, b) ASSERT_GREATER(a, b);
#define ASSERT_LE(a, b) ASSERT_LESSTHAN_EQUAL(a, b);
#define ASSERT_GE(a, b) ASSERT_GREATER_EQUAL(a, b);
#define ASSERT_NE(a, b) ASSERT_NOT_EQUALS(a, b);

// Note about floating point comparison:
// A floating point value 'a' is equal to a floating point value 'b'
// if the two values are not too far apart (defined by an epsilon value, eps).

// We check that either the value difference is equal, or for large values,
// the difference is below the eps percentage (to allow for some error).
#define FLOAT_EQUALS(a, b) (TESTABS((float)(a) - (float)(b)) < epsilon \
		|| TESTABS((float) (a) - (float) (b)) <= TESTABS((float) (b) * epsilon))

// Asserts that a == b, given note above
#define ASSERT_EQUALS_FLOAT(a, b) \
	TEST(#a " is supposed to equal (float) " #b); \
	sprintf(info, "\n\t(%"FSYM" != %"FSYM")", (float) (a), (float) (b)); \
 	ASSERT(FLOAT_EQUALS(a, b));

// Asserts that a == b
#define ASSERT_EQUALS_INT(a, b) \
	TEST(#a " is supposed to equal (int) " #b); \
	sprintf(info, "\n\t(%"ISYM" != %"ISYM")", (int) (a), (int) (b)); \
	ASSERT_EQUALS(a, b);

// Asserts that a != b, given note above
#define ASSERT_NOT_EQUALS_FLOAT(a, b) \
	TEST(#a " is not supposed to equal (float) " #b); \
	sprintf(info, "\n\t(%f == %f)", (float) (a), (float) (b)); \
	ASSERT(!(FLOAT_EQUALS(a, b)));

// The main assertion
// If the assertion fails, we print out the test information
#define ASSERT(test) \
	if (reason == NULL) { \
		TEST(#test " fails"); \
	} \
	if(!(test)) {  \
	  printf("Proc %"ISYM": [FAIL] %s:%d : (%s) : %s %s\n", MyProcessorNumber, \
		 __FILE__, __LINE__, name, reason, info);		\
		TESTFAIL(); \
	} else { \
		passes++; \
	} \
	reason = NULL; \
	sprintf(info, " ");

// Automatic assertion failure, with a reason
// Takes test execution configuration into consideration
#ifdef TESTS_EXIT_ON_FAILURE
#define TESTFAIL() \
	fails++; \
	testfails++; \
	TESTS_HALT();
#else
#ifdef TESTS_BREAK_ON_FAILURE
#define TESTFAIL() \
	fails++; \
	testfails++; \
	TESTS_BREAK();
#else
#define TESTFAIL() \
	fails++; \
	testfails++;
#endif // tests_break_on_failure
#endif // tests_exit_on_failure

// --- Execution control ---
// Break execution of a test suite (only use this in a test suite!)
#define TESTS_BREAK() \
        if (MyProcessorNumber == ROOT_PROCESSOR) { \
          printf("[break] Test failure; aborting current test suites.\n"); \
        } \
	return 2;

// Break execution of ALL test suites (only use this in a test suite!)
#define TESTS_HALT() \
        if (MyProcessorNumber == ROOT_PROCESSOR) { \
           printf("[break] Test failure; aborting all test suites.\n"); \
        } \
	return 0;

// --- Miscellaneous ---
// Set the epsilon value
#ifdef TESTS_IGNORE_EPSILON
#define SET_EPSILON(f) epsilon = (f); SET_EPSILON_DEFAULT();
#else
#define SET_EPSILON(f) epsilon = (f);
#endif	// tests_ignore_epsilon

// Set the epsilon value to the default
#define SET_EPSILON_DEFAULT() epsilon = GET_EPSILON_DEFAULT();

// Get the epsilon value
#define GET_EPSILON() (epsilon)

// Get the default epsilon value
#define GET_EPSILON_DEFAULT() (0.0001)			// 0.01%

#ifdef SHUT_UP
#define stdout nowhere
#define stderr nowhere
#endif

#endif // _SIMPLECTEST_TESTS_H


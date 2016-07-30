/*
* Copyright (C) 2011 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


/*! \file
\brief Testing adaptive histograms, including arithmetic.

Run the executable and then use the shell script testing_hists.sh
to run checks out output.

 */

#include "testing_tools.hpp"
#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data
#include "subpaving_exception.hpp"


#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <algorithm> // equal
#include <stdexcept>
#include <cassert>

using namespace cxsc;
using namespace std;
using namespace subpavings;

void testHistArithmetic();
RVecData& getData1(RVecData& data);
RVecData& getData2(RVecData& data);
RVecData& getData3(RVecData& data);


int main()
{

	try {
		cout << "\nDefault constructor" << endl;
		
		AdaptiveHistogram adh;
		
		assert( adh.getDimensions() == 0 );
		cout << "Passed assert that getDimensions() == 0" << endl;
		
		string s = "defaultConstructor.txt";
		outputADH(s, adh);
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
		assert(!adh.getHoldAllStats());
		cout << "Passed assert holdAllStats is false " << endl;
		
	}
	catch (std::exception& ee) {
		cout << "\nFailed to do default constructor:\n" << ee.what() << endl;
		throw;
	}

	try {
			
		cout << "\nTry to split with default box (this should fail)" << endl;
		
		AdaptiveHistogram temp;
		
		temp.splitToShape("2,2,1");
		
		throw std::logic_error ("Should not be able to do this");

	}
	
	catch (std::exception const& ee) {
		cout << "\nException:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nMean with default constructor (this should fail)" << endl;
		
		AdaptiveHistogram adh;
		
		rvector mean = adh.getRootPavingMean();
		
		throw std::logic_error ("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do mean with default constructor:\n" << msg << endl;
	}
	try {
		cout << "\nVariance-covariance with default constructor (this should fail)" << endl;
		
		AdaptiveHistogram adh;
		
		RealVec varcov = adh.getRootPavingVarCovar();
		
		throw std::logic_error ("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do variance-covariance with default constructor:\n" << msg << endl;
	}
	try {
		cout << "\nCopy constructor with default constructor" << endl;
		
		AdaptiveHistogram adh1;
		AdaptiveHistogram adh2(adh1);
		
		string s = "copyOfDefaultConstructor.txt";
		outputADH(s, adh2);
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy constructor with default constructor:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nCopy assignment with default constructor" << endl;
		
		AdaptiveHistogram adh1;
		AdaptiveHistogram adh2 = adh1;
		
		string s = "copyAssignmentOfDefaultConstructor.txt";
		outputADH(s, adh2);
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy assignment with default constructor:\n" << msg << endl;
		throw;
	}

	try {
			
		cout << "\nConstructor with default (malformed) box (this should fail)" << endl;
		
		ivector pavingBox; // Ub < Ub
		AdaptiveHistogram temp(pavingBox);
		
		throw std::logic_error ("Should not be able to do this");

	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor with default (malformed) box:\n" << ee.what() << endl;
	}
	
	try {
			
		cout << "\nConstructor with thin-interval (malformed) box (this should fail)" << endl;
		
		int d = 2;
		ivector pavingBoxT1(d);
		interval pavingInterval1(-2,2);
		interval pavingInterval2(2,2); // thin
		pavingBoxT1[1] = pavingInterval1;
		pavingBoxT1[2] = pavingInterval2;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		throw std::logic_error ("Should not be able to do this");	
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor with thin-interval (malformed) box:\n" << ee.what() << endl;
	}
	
	try {
			
		cout << "\nConstructor with 1-d box" << endl;
		
		int d = 1;
		ivector pavingBoxT1(d);
		interval pavingIntervalT1(-2,2);
		for(int k=1; k <= d; k++) pavingBoxT1[k] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		assert( temp.getDimensions() == d );
		cout << "Passed assert that getDimensions() == " << d << endl;
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor 1-d box:\n" << ee.what() << endl;
		throw;
	}
	
	try {
			
		cout << "\nConstructor with 1-d box and split to shape" << endl;
		
		int d = 1;
		ivector pavingBoxT1(d);
		interval pavingIntervalT1(-2,2);
		for(int k=1; k <= d; k++) pavingBoxT1[k] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		assert( temp.getDimensions() == d );
		
		temp.splitToShape("2,2,1");
		
		size_t shouldBe = 3;
		assert(temp.getRootLeaves() == shouldBe);
		cout << "Passed assert that getRootLeaves() == " << shouldBe << endl;
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor 1-d box and splitToShape:\n" << ee.what() << endl;
		throw;
	}

	try {
			
		cout << "\nConstructor with 2-d box" << endl;
		
		int d = 2;
		ivector pavingBoxT1(d);
		interval pavingIntervalT1(-2,2);
		for(int k=1; k <= d; k++) pavingBoxT1[k] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		assert( temp.getDimensions() == d );
		cout << "Passed assert that getDimensions() == " << d << endl;
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor 2-d box:\n" << ee.what() << endl;
		throw;
	}
	
	try {
			
		cout << "\nConstructor with 2-d box and split to shape" << endl;
		
		int d = 2;
		ivector pavingBoxT1(d);
		interval pavingIntervalT1(-2,2);
		for(int k=1; k <= d; k++) pavingBoxT1[k] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		assert( temp.getDimensions() == d );
		
		temp.splitToShape("2,2,1");
		size_t shouldBe = 3;
		assert(temp.getRootLeaves() == shouldBe);
		cout << "Passed assert that getRootLeaves() == " << shouldBe << endl;
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor 2-d box and splitToShape:\n" << ee.what() << endl;
		throw;
	}
	
	try {
			
		cout << "\nMalformed split to shape instruction" << endl;
		cout << "\nWe should get an unchanged box" << endl;
		
		int d = 1;
		ivector pavingBoxT1(d);
		interval pavingIntervalT1(-2,2);
		for(int k=1; k <= d; k++) pavingBoxT1[k] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		assert( temp.getDimensions() == d );
		
		temp.splitToShape("");
		
		size_t shouldBe = 1;
		assert(temp.getRootLeaves() == shouldBe);
		cout << "Passed assert that getRootLeaves() == " << shouldBe << endl;
		
	}
	
	catch (std::exception& ee) {
		cout << "\nFailed to do Constructor 1-d box and splitToShape:\n" << ee.what() << endl;
		throw;
	}


	// ------- prepare to generate some data for more tests -----------
	
    // set up a random number generator
    const gsl_rng_type * T;
    gsl_rng * r;

    double sigma_x=1;   // distribution parameter for Biv Gaussian
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    const int n=10;    // number to generate
	const int n4=7;		// for hist4
    //create a generator chosen by the environment variable GSL_RNG_TYPE

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    string samplesFileName; // for samples
    string outputFileName;// for output file
    ofstream oss;         // ofstream object
    oss << scientific;  // set formatting for input to oss
    oss.precision(5);

    int d1 = 2; // dimension of the box to sample data from
    ivector pavingBox1(d1);
    interval pavingInterval1(-5,5);
    for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;

    RVecData theData1;   // a container for all the points generated

    // make a simulated data set to sample from
    for (int i = 0; i < n; i++) {

        rvector thisrv(d1);
        double x = 0;
        double y = 0;

        gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                rho, &x, &y);
        thisrv[1] = x;
        thisrv[2] = y;

        // put points generated into container
        theData1.push_back(thisrv);

    } // data  should be in theData

    
	RVecData theData2;   // a container for all the points generated

	// make a simulated data set allData to sample from
	for (int i = 0; i < n; i++) {

		rvector thisrv(d1);
		double x = 0;
		double y = 0;

		gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
								rho, &x, &y);
		thisrv[1] = x;
		thisrv[2] = y;

		// put points generated into container
		theData2.push_back(thisrv);

		}  // data  should be in theData

	
	RVecData theData3;   // a container for all the points generated

	// make a simulated data set allData to sample from
	for (int i = 0; i < n; i++) {

		rvector thisrv(d1);
		double x = 0;
		double y = 0;

		gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
								rho, &x, &y);
		thisrv[1] = x;
		thisrv[2] = y;

		// put points generated into container
		theData3.push_back(thisrv);

	}  // data  should be in theData

	RVecData theData4;   // a container for all the points generated

	// make a simulated data set allData to sample from
	for (int i = 0; i < n4; i++) { // note different size

		rvector thisrv(d1);
		double x = 0;
		double y = 0;

		gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
								rho, &x, &y);
		thisrv[1] = x;
		thisrv[2] = y;

		// put points generated into container
		theData4.push_back(thisrv);

	}  // data  should be in theData


    // free the random number generator
    gsl_rng_free (r);
	
	try {
		cout << "\nInsertion of data into histogram with no box" << endl;
		
		AdaptiveHistogram noBox;
		bool insertion = noBox.insertFromRVec(theData1);
		assert(insertion);
		cout << "Passed assert that insertion was successful" << endl;
		string s = "insertionNoBox.txt";
		outputADH(s, noBox);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do insertion of data into histogram with no box:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nInsertion of empty container of data into histogram" << endl;
		
		AdaptiveHistogram noBox;
		RVecData data;
		bool insertion = noBox.insertFromRVec(data);
		assert(!insertion);
		cout << "Passed assert that insertion was not successful" << endl;
		string s = "insertionEmptyData.txt";
		outputADH(s, noBox);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do insertion of empty container of data into histogram:\n" << msg << endl;
		throw;
	}

	try {
			
		cout << "\nConstruct a histogram that does not fit any of the data" << endl;
		
		ivector pavingBoxT1(d1);
	
		for(int k=1; k < d1; k++) pavingBoxT1[k] = pavingInterval1;
		
		interval pavingIntervalT1(100,101);
		pavingBoxT1[d1] = pavingIntervalT1;
		
		AdaptiveHistogram temp(pavingBoxT1);
		
		temp.splitToShape("2,2,1");
		
		bool insertion = temp.insertFromRVec(theData1);
		
		assert(!insertion);
		cout << "Passed assert that insertion was not successful" << endl;
		assert(!temp.getRootCounter());
		cout << "Passed assert that root counter == 0" << endl;
		
	}
	
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct a histogram that does not fit any of the data:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nConstruct a histogram that only fits some of the data" << endl;
		
		AdaptiveHistogram temp(pavingBox1);
		
		temp.splitToShape("2,2,1");
		
		RVecData data;
		data = getData3(data);
		
		bool insertion = temp.insertFromRVec(data);
		cout << "rootCounter is " << temp.getRootCounter() << endl;
		cout << "size of dataCollection is " << temp.getDataCollection().size() << endl;
		assert(insertion);
		cout << "Passed assert that insertion was successful" << endl;
		assert(temp.getRootCounter() == temp.getDataCollection().size() );
		cout << "Passed assert that root counter == size of data collection" << endl;
		
	}
	
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct a histogram that only fits some of the data:\n" << msg << endl;
		throw;
	}
	
	

	try {
		cout << "\nInsertion of data into histograms for further tests" << endl;
		cout << "\nhistograms heapHist1 and heapHeap2 are on the heap and have holdAllStats = false and are only used for an initial test" << endl;
		cout << "\nhistogram 1 has holdAllStats = true, Histograms 2, 3, 4 have default holdAllStats = false" << endl;

		AdaptiveHistogram * heapHistPtr1 = new AdaptiveHistogram(pavingBox1); 
		bool successfulInsertionHeapHist1 = false;
		
		try {
			heapHistPtr1->splitToShape("3,4,4,2,2,4,4,3");

			// put in the data in a 'pulse' with no further splitting
			successfulInsertionHeapHist1 = heapHistPtr1->insertFromRVec(theData1);
			if (!successfulInsertionHeapHist1) cout << "unsuccessful insertion heap hist1" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split heapHist1 to shape and insert data:\n" << msg << endl;
			
			throw;
		}
		
		AdaptiveHistogram * heapHistPtr2 = new AdaptiveHistogram(pavingBox1); 
		bool successfulInsertionHeapHist2 = false;
		
		try {
			heapHistPtr2->splitToShape("3,4,4,2,2,4,4,3");

			// put in the data in a 'pulse' with no further splitting
			successfulInsertionHeapHist2 = heapHistPtr2->insertFromRVec(theData1);
			if (!successfulInsertionHeapHist2) cout << "unsuccessful insertion heap hist2" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split heapHist1 to shape and insert data:\n" << msg << endl;
			
			throw;
		}
		
		AdaptiveHistogram myHistFirst(pavingBox1, true); // keep all stats for this one
		bool successfulInsertionHistFirst = false;
		
		try {
			cout << "\nHist1: Split to shape and then data insertion\n" << endl;
			
			myHistFirst.splitToShape("3,4,4,2,2,4,4,3");

			// put in the data in a 'pulse' with no further splitting
			successfulInsertionHistFirst = myHistFirst.insertFromRVec(theData1);
			if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split Hist1 to shape and insert data:\n" << msg << endl;
			
			throw;
		}
		
		AdaptiveHistogram myHistSecond(pavingBox1);
		bool successfulInsertionHistSecond = false;
		try {

			// applying SEB heuristics for k to satisfy k/n -> 0 as n -> +oo
			int k_int = n/2;
			cout << "\nHist2: Split with insertion, splitOnK, with k = " << k_int << endl;

			bool successfulInsertion = false;
			bool successfulPQSplit = false;

			// make the function object to get max k_int data members in each box
			SplitOnK splitK(k_int);

			// insert data into the histogram, splitting as we go, no logging
			successfulInsertionHistSecond = myHistSecond.insertFromRVec(theData2, splitK, NOLOG);
			if (!successfulInsertionHistSecond) cout << "unsuccessful insertion 2" << endl;

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split Hist2 using splitOnK:\n" << msg << endl;
			
			throw;
		}
		
		AdaptiveHistogram myHistThird(pavingBox1);
		bool successfulInsertionHistThird = false;
		
		try {
			successfulInsertionHistThird = myHistThird.insertFromRVec(theData3);
			if (successfulInsertionHistThird) {

				// set up priority queue parameters
				// split box with most counts
				CompCount nodeCompCount;
				// stop when number of leaves is critLeavesGTE
				size_t maxLeaves = 5;
				CritLeaves_GTE critLeavesGTE(maxLeaves);
				
				cout << "\nHist3: : Priority queue insertion, splitting nodes with most in first, until number of nodes is >= " << maxLeaves << endl;
				
				// now split with priority queue
			   // no minPoints or minVolB limitations on splittable nodes
			   double minvol = 0.0;
				successfulInsertionHistThird = 
					myHistThird.prioritySplit(nodeCompCount,
										critLeavesGTE, NOLOG, minvol); //no logs
				if (!successfulInsertionHistThird) cout << "Priority queue failed" << endl;
				
			}
			else  cout << "unsuccessful insertion 3" << endl;
		}

		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split Hist3 using priority queue:\n" << msg << endl;
			
			throw;
		}
		
		AdaptiveHistogram myHistFourth(pavingBox1); // 
		bool successfulInsertionHistFourth = false;
		
		try {
			myHistFourth.splitToShape("3,4,4,2,2,4,4,3"); // same splits as for 1
			successfulInsertionHistFourth = myHistFourth.insertFromRVec(theData4);
		}

		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to split Hist4 to shape and insert data:\n" << msg << endl;
			
			throw;
		}
		
		if (!successfulInsertionHeapHist1 || !successfulInsertionHeapHist2) {
			throw std::runtime_error("Could not insert data into heap hists");
		}

		try {
			cout << "\nTry copy constructing from heapHist1, deleting the original, and checking that we can access the data" << endl;
			
			assert(!heapHistPtr1->getHoldAllStats());
			
			string s1("doubleCheckHeapHist1.txt");
			doubleCheckOutput(s1, *heapHistPtr1);
			
			
			AdaptiveHistogram temp(*heapHistPtr1);
			string s2("doubleCheckCopyHeapHist1BeforeResettingHoldAllStats.txt");
			doubleCheckOutput(s2, temp);
			
			assert(!heapHistPtr1->getHoldAllStats());
			
			delete heapHistPtr1;
			heapHistPtr1 = NULL;
			
			temp.setHoldAllStats(true);
			string s3("doubleCheckCopyHeapHist1AfterResettingHoldAllStats.txt");
			doubleCheckOutput(s3, temp);
			
			cout << "Sucessful" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to copy constructing from heapHist1, deleting the original, etc:\n" << msg << endl;
			
			throw;
		}
		
		try {
			cout << "\nTry assignment copy constructing from heapHist2, deleting the original, and checking that we can access the data" << endl;
			
			assert(!heapHistPtr2->getHoldAllStats());
			
			string s1("doubleCheckHeapHist2.txt");
			doubleCheckOutput(s1, *heapHistPtr2);
			
			
			AdaptiveHistogram temp(*heapHistPtr2);
			string s2("doubleCheckCopyHeapHist2BeforeResettingHoldAllStats.txt");
			doubleCheckOutput(s2, temp);
			
			assert(!heapHistPtr2->getHoldAllStats());
			
			delete heapHistPtr2;
			heapHistPtr2 = NULL;
			
			temp.setHoldAllStats(true);
			string s3("doubleCheckCopyHeapHist2AfterResettingHoldAllStats.txt");
			doubleCheckOutput(s3, temp);
			
			cout << "Sucessful" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to copy constructing from heapHist2, deleting the original, etc:\n" << msg << endl;
			
			throw;
		}
		
		if (!(successfulInsertionHistFirst && successfulInsertionHistSecond 
						&& successfulInsertionHistThird && successfulInsertionHistFourth)) {
			throw std::runtime_error("Could not insert data into all hists");
		}
		cout << "\nChecking other histograms:" << endl;
							
		assert(myHistFirst.getHoldAllStats());
		cout << "Passed assert hist1 holdAllStats is true " << n << endl;
		assert(!myHistSecond.getHoldAllStats());
		assert(!myHistThird.getHoldAllStats());
		assert(!myHistFourth.getHoldAllStats());
		cout << "Passed assert hist2, hist3, hist4 holdAllStats are false " << n << endl;
						
		assert(myHistFirst.getRootCounter() == n);
		assert(myHistSecond.getRootCounter() == n);
		assert(myHistThird.getRootCounter() == n);
		cout << "Passed asserts hist1, hist2, hist3 root counters are " << n << endl;
		assert(myHistFourth.getRootCounter() == n4);
		cout << "Passed asserts hist4 root counter is " << n4 << endl;
		
		{
			rvector mean = myHistFourth.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
		}
		{
			rvector mean = myHistSecond.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
		}
		{
			rvector mean = myHistThird.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
		}
		cout << "Passed asserts that root paving means for hists 2, 3 and 4 all have length " << d1 << " and are all NaN" << endl;
		
		{
			RealVec vcov = myHistFourth.getRootPavingVarCovar();
			assert(vcov.size() == d1*d1);
			assert(checkAllNaN(vcov));
		}
		{
			RealVec vcov = myHistSecond.getRootPavingVarCovar();
			assert(vcov.size() == d1*d1);
			assert(checkAllNaN(vcov));
		}
		{
			RealVec vcov = myHistThird.getRootPavingVarCovar();
			assert(vcov.size() == d1*d1);
			assert(checkAllNaN(vcov));
		}
		cout << "Passed asserts that root paving variance-covariances for hists 2, 3 and 4 all have length " << d1*d1 << " and are all NaN" << endl;
		
		rvector mean1 = myHistFirst.getRootPavingMean();
		
		assert(VecLen(mean1) == d1);
		assert(!checkAllNaN(mean1));
		rvector calcMean1 = checkMean(theData1);
		
		assert(checkSame(mean1, calcMean1, n));
		cout << "Passed asserts that mean for histogram 1 is length " << d1 << " and calculation of mean matches check calc" << endl;
		
		rvector mean1dc = myHistFirst.getDataCollectionMean();
		
		assert(VecLen(mean1dc) == d1);
		assert(!checkAllNaN(mean1dc));
		
		assert(checkSame(mean1dc, calcMean1, n));
		cout << "Passed asserts that mean for histogram 1 from data collection is length " << d1 << " and calculation of mean matches check calc" << endl;
		
		
		RealVec vcov1 = myHistFirst.getRootPavingVarCovar();
		
		assert(vcov1.size() == d1*d1);
		assert(!checkAllNaN(vcov1));
		RealVec calcVCov1 = checkVarCov(theData1);
		assert(calcVCov1.size() == d1*d1);
		assert(checkSame(vcov1, calcVCov1, n) );
		cout << "Passed asserts that variance-covariance for histogram 1 is length " << d1*d1 << " calculation of variance-covariance matches check calc" << endl;

		RealVec vcov1dc = myHistFirst.getDataCollectionVarCovar();
		
		assert(vcov1dc.size() == d1*d1);
		assert(!checkAllNaN(vcov1dc));
		assert(checkSame(vcov1dc, calcVCov1, n) );
		cout << "Passed asserts that variance-covariance for histogram 1 from data collection is length " << d1*d1 << " calculation of variance-covariance matches check calc" << endl;

		rvector mean4dc = myHistFourth.getDataCollectionMean();
		
		assert(VecLen(mean4dc) == d1);
		assert(!checkAllNaN(mean4dc));
		rvector calcMean4dc = checkMean(theData4);
		
		assert(checkSame(mean4dc, calcMean4dc, n4));
		cout << "Passed asserts that mean for histogram 4 from data collection is length " 
			<< d1 << " and calculation of mean matches check calc" << endl;
	
		RealVec vcov4dc = myHistFourth.getDataCollectionVarCovar();
		
		assert(vcov4dc.size() == d1*d1);
		assert(!checkAllNaN(vcov4dc));
		RealVec calcVCov4dc = checkVarCov(theData4);
		assert(calcVCov4dc.size() == d1*d1);
		assert(checkSame(vcov4dc, calcVCov4dc, n) );
		cout << "Passed asserts that variance-covariance for histogram 4 from data collection is length " 
			<< d1*d1 << " calculation of variance-covariance matches check calc" << endl;

		string outputFileName1 = "Hist1basic.txt";
		outputADH(outputFileName1, myHistFirst);
		
		string outputFileName2 = "Hist2basic.txt";
		outputADH(outputFileName2, myHistSecond);
		
		string outputFileName3 = "Hist3basic.txt";
		outputADH(outputFileName3, myHistThird);
		
		string outputFileName4 = "Hist4basic.txt";
		outputADH(outputFileName4, myHistFourth);
		
		try {
			
			cout << "\ncopy construct copy of Hist1 (holdAllStats = true)" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
			assert( temp.getRootCounter() == myHistFirst.getRootCounter() );
			cout << "Passed assert that root counters of original and copy are the same" << endl;
			assert( temp.getHoldAllStats() == myHistFirst.getHoldAllStats() );
			cout << "Passed assert holdAllStats is same " << n << endl;
	
			rvector mn = temp.getRootPavingMean();
		
			assert(VecLen(mn) == d1);
			assert(!checkAllNaN(mn));
			assert(checkSame(mean1, mn, n));
			cout << "Passed asserts that mean for copy of histogram 1 is length " << d1 << " and is same as for original" << endl;
			
			
			RealVec vc = temp.getRootPavingVarCovar();
			
			assert(vc.size() == d1*d1);
			assert(!checkAllNaN(vc));
			assert(checkSame(vcov1, vc, n) );
			cout << "Passed asserts that variance-covariance for copy of histogram 1 is length " << d1*d1 << " and is same as for original" << endl;

			string outputFileName = "copyConstructHist1basic.txt";
			outputADH(outputFileName, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction for hist1:\n" << msg << endl;
			throw;
		}
		try {
			
			cout << "\ncopy construct copy of Hist4 (hold all stats = false)" << endl;
			
			AdaptiveHistogram temp(myHistFourth);
			
			assert( temp.getRootCounter() == myHistFourth.getRootCounter() );
			cout << "Passed assert that root counters of original and copy are the same" << endl;
			
			assert( temp.getHoldAllStats() == myHistFourth.getHoldAllStats() );
			cout << "Passed assert that holdAllStats of original and copy are the same" << endl;
			rvector mn = temp.getRootPavingMean();
			assert(VecLen(mn) == d1);
			assert(checkAllNaN(mn));
			RealVec vc = temp.getRootPavingVarCovar();
			assert(vc.size() == d1*d1);
			assert(checkAllNaN(vc));
			cout << "Passed asserts that root counters and mean and variance-covariance of original and copy are the same" << endl;

			string outputFileName = "copyConstructHist4basic.txt";
			outputADH(outputFileName, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction for hist4:\n" << msg << endl;
			throw;
		}
		
		
		try {
			
			cout << "\ncopy construct copy of Hist4 (hold all stats = false) and then reset holdAllStats to true" << endl;
			
			string s1("doubleCheckHist4.txt");
			doubleCheckOutput(s1, myHistFourth);
			
			
			AdaptiveHistogram temp(myHistFourth);
			
			assert( temp.getRootCounter() == myHistFourth.getRootCounter() );
			cout << "Passed assert that root counters of original and copy are the same" << endl;
			
			assert( temp.getHoldAllStats() == myHistFourth.getHoldAllStats() );
			cout << "Passed assert that holdAllStats of original and copy are the same" << endl;
			
			string s2("doubleCheckCopyHist4Before.txt");
			doubleCheckOutput(s2, temp);
			
			temp.setHoldAllStats(true);
			
			string s3("doubleCheckCopyHist4After.txt");
			doubleCheckOutput(s3, temp);
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is now true" << endl;
			
			rvector mn = temp.getRootPavingMean();
			assert(VecLen(mn) == d1);
			rvector mean4 = checkMean(theData4);
			assert( checkSame( mn, mean4, temp.getRootCounter() ) );
			RealVec vc = temp.getRootPavingVarCovar();
			RealVec vc4 = checkVarCov(theData4);
			assert(vc.size() == d1*d1);
			assert( checkSame( vc, vc4, temp.getRootCounter() ) );
			cout << "Passed asserts that root counter is same as original and mean and variance-covariance match check calcs" << endl;

			string outputFileName = "copyConstructHist4ResetHoldAllStatsTrue.txt";
			outputADH(outputFileName, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction for hist4 and reset holdAllStats to true:\n" << msg << endl;
			throw;
		}
		 
		AdaptiveHistogram copyHist1;
		AdaptiveHistogram copyHist2;
		
		try {
			
			cout << "\ncopy assignment copies of Hist1 and Hist 2" << endl;
			
			copyHist1 = myHistFirst;
			copyHist2 = myHistSecond;

			{
				assert( copyHist1.getRootCounter() == myHistFirst.getRootCounter() );
				assert( copyHist1.getHoldAllStats() == myHistFirst.getHoldAllStats() );
				rvector mn = copyHist1.getRootPavingMean();
				assert(VecLen(mn) == d1);
				assert(!checkAllNaN(mn));
				assert(checkSame(mean1, mn, n));
				RealVec vc = copyHist1.getRootPavingVarCovar();
				assert(vc.size() == d1*d1);
				assert(!checkAllNaN(vc));
				assert(checkSame(vcov1, vc, n) );
				cout << "Passed asserts that root counters, mean, and variance-covariance for copy of histogram 1 are same as for original" << endl;
			}
			{
				assert( copyHist2.getRootCounter() == myHistSecond.getRootCounter() );
				assert( copyHist2.getHoldAllStats() == myHistSecond.getHoldAllStats() );
				rvector mn = copyHist2.getRootPavingMean();
				assert(VecLen(mn) == d1);
				assert(checkAllNaN(mn));
				RealVec vc = copyHist2.getRootPavingVarCovar();
				assert(vc.size() == d1*d1);
				assert(checkAllNaN(vc));
				cout << "Passed asserts that root counters, mean, and variance-covariance for copy of histogram 1 are same as for original" << endl;
			}
			
			string outputFileName1 = "copyAssignHist1basic.txt";
			outputADH(outputFileName1, copyHist1);
			string outputFileName2 = "copyAssignHist2basic.txt";
			outputADH(outputFileName2, copyHist2);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy assignment for hist1 and hist2:\n" << msg << endl;
			throw;
		}
		
		try {
			
			cout << "\nCopy assignment copy of Hist4 (hold all stats = false) and then reset holdAllStats to true" << endl;
			
			AdaptiveHistogram temp = myHistFourth;
			
			assert( temp.getRootCounter() == myHistFourth.getRootCounter() );
			cout << "Passed assert that root counters of original and copy are the same" << endl;
			
			assert( temp.getHoldAllStats() == myHistFourth.getHoldAllStats() );
			cout << "Passed assert that holdAllStats of original and copy are the same" << endl;
			
			string s2("doubleCheckAssignmentCopyHist4Before.txt");
			doubleCheckOutput(s2, temp);
			
			temp.setHoldAllStats(true);
			
			string s3("doubleCheckAssignmentCopyHist4After.txt");
			doubleCheckOutput(s3, temp);
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is now true" << endl;
			
			rvector mn = temp.getRootPavingMean();
			assert(VecLen(mn) == d1);
			rvector mean4 = checkMean(theData4);
			assert( checkSame( mn, mean4, temp.getRootCounter() ) );
			RealVec vc = temp.getRootPavingVarCovar();
			RealVec vc4 = checkVarCov(theData4);
			assert(vc.size() == d1*d1);
			assert( checkSame( mn, mean4, temp.getRootCounter() ) );
			cout << "Passed asserts that root counter is same as original and mean and variance-covariance match check calcs" << endl;

			string outputFileName = "copyAssignmentHist4ResetHoldAllStatsTrue.txt";
			outputADH(outputFileName, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction for hist4 and reset holdAllStats to true:\n" << msg << endl;
			throw;
		}
	
		
		// clear all data
	
		AdaptiveHistogram clearedHist1;
		try {
			cout << "\nClear a copy of hist1 and check counts etc" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
			
			temp.clearAllHistData();
			
			assert(temp.getRootCounter() == 0);
			cout << "Passed assert that root counter is 0" << endl;
			BigDataCollection tmp = temp.getDataCollection();
			assert(tmp.empty());
			cout << "Passed assert that dataCollection is empty" << endl;
			assert( temp.getHoldAllStats() == myHistFirst.getHoldAllStats() );
			cout << "Passed assert that holdAllStats is unchanged" << endl;	
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			clearedHist1 = temp;
			
			string s = "clearedCopyHist1.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do clear a copy of hist1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nReinsert same data into cleared copy of hist1 and check counts etc" << endl;
			
			AdaptiveHistogram temp(clearedHist1);
			
			temp.insertFromRVec(theData1);
			
			assert(temp.getRootCounter() == myHistFirst.getRootCounter());
			cout << "Passed assert that root counter is same as for hist1" << endl;
			assert(temp.getDataCollection().size() == myHistFirst.getDataCollection().size());
			cout << "Passed assert that dataCollection is same as for hist1" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			rvector mean1 = myHistFirst.getRootPavingMean();
			assert( checkSame( mean, mean1, temp.getRootCounter() ) );
			cout << "Passed asserts that mean is same as for hist1" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			RealVec varcov1 = myHistFirst.getRootPavingVarCovar();
			assert( checkSame( varcov, varcov1, temp.getRootCounter() ) );
			
			cout << "Passed asserts that variance-covariance is same as for hist1" << endl;
			
			string s = "clearedCopyHist1WithDataReinserted.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reinsert same data into cleared copy of hist1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nReinsert hist4 data into cleared copy of hist1 and check counts etc" << endl;
			
			AdaptiveHistogram temp(clearedHist1);
			
			temp.insertFromRVec(theData4);
			
			assert(temp.getRootCounter() == n4);
			cout << "Passed assert that root counter is n4 = " << n4 << endl;
			assert(temp.getDataCollection().size() == n4);
			cout << "Passed assert that dataCollection is n4 = " << n4 << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			rvector mean4 = checkMean(theData4);
			assert( checkSame( mean, mean4, temp.getRootCounter() ) );
			cout << "Passed asserts that mean is as calculated by check calcs" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame( varcov, varcov4, temp.getRootCounter() ) );
			
			cout << "Passed asserts that variance-covariance is is as calculated by check calcs" << endl;
			
			string s = "clearedCopyHist1WithDataForHist4Reinserted.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reinsert hist4 data into cleared copy of hist1:\n" << msg << endl;
			throw;
		}
		
		AdaptiveHistogram clearedHist4;
		try {
			cout << "\nClear a copy of hist4 and check counts etc" << endl;
			
			AdaptiveHistogram temp(myHistFourth);
			
			temp.clearAllHistData();
			
			assert(temp.getRootCounter() == 0);
			cout << "Passed assert that root counter is 0" << endl;
			BigDataCollection tmp = temp.getDataCollection();
			assert(tmp.empty());
			cout << "Passed assert that dataCollection is empty" << endl;
			assert( temp.getHoldAllStats() == myHistFourth.getHoldAllStats() );
			cout << "Passed assert that holdAllStats is unchanged" << endl;	
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			clearedHist4 = temp;
			
			string s = "clearedCopyHist4.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do clear a copy of hist4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nReinsert same data into cleared copy of hist4 and check counts etc" << endl;
			
			AdaptiveHistogram temp(clearedHist4);
			
			temp.insertFromRVec(theData4);
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter());
			cout << "Passed assert that root counter is same as for hist4" << endl;
			assert(temp.getDataCollection().size() == myHistFourth.getDataCollection().size());
			cout << "Passed assert that dataCollection is same as for hist4" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is not held" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			
			cout << "Passed asserts that variance-covariance is not held" << endl;
			
			string s = "clearedCopyHist4WithDataReinserted.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reinsert same data into cleared copy of hist1:\n" << msg << endl;
			throw;
		}
		AdaptiveHistogram clearedHist4WithData1;
		try {
			cout << "\nReinsert hist1 data into cleared copy of hist4 and check counts etc" << endl;
			
			AdaptiveHistogram temp(clearedHist4);
			
			temp.insertFromRVec(theData1);
			
			assert(temp.getRootCounter() == n);
			cout << "Passed assert that root counter is n = " << n << endl;
			assert(temp.getDataCollection().size() == n);
			cout << "Passed assert that dataCollection is n = " << n << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is not held" << endl;
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is not held" << endl;
			
			clearedHist4WithData1 = temp;
			
			string s = "clearedCopyHist4WithDataForHist1Reinserted.txt";
			outputADH(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reinsert hist1 data into cleared copy of hist4:\n" << msg << endl;
			throw;
		}
		
		// change holdAllStats
		
		try {
			cout << "\nreset holdAllStats for copy of hist1" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
				
			{
				cout << "\nChange holdAllStats for copy of hist1 to false" << endl;
			
				bool newHoldAllStats = false;
				temp.setHoldAllStats(newHoldAllStats);
				assert(temp.getHoldAllStats() == newHoldAllStats);
				cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;
				assert(temp.getRootCounter() == myHistFirst.getRootCounter());
				cout << "Passed assert that root counter is unchanged" << endl;
				rvector mean = temp.getRootPavingMean();
				assert(VecLen(mean) == d1);
				assert(checkAllNaN(mean));
				cout << "Passed asserts that mean is not held" << endl;
				RealVec varcov = temp.getRootPavingVarCovar();
				assert(varcov.size() == d1*d1);
				assert(checkAllNaN(varcov));
				cout << "Passed asserts that variance-covariance is not held" << endl;
				
				string s = "holdAllStatsResetToFalseForCopyHist1.txt";
				outputADH(s, temp);
			}
			{
				cout << "\nChange holdAllStats for copy of hist1 to back to true" << endl;
			
				bool newHoldAllStats = true;
				temp.setHoldAllStats(newHoldAllStats);
				assert(temp.getHoldAllStats() == newHoldAllStats);
				cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;
				assert(temp.getRootCounter() == myHistFirst.getRootCounter());
				cout << "Passed assert that root counter is unchanged" << endl;
				rvector mean = temp.getRootPavingMean();
				assert(VecLen(mean) == d1);
				assert(!checkAllNaN(mean));
				assert( checkSame (mean, myHistFirst.getRootPavingMean(), temp.getRootCounter() ) );
				cout << "Passed asserts that mean is same as for hist 1" << endl;
				RealVec varcov = temp.getRootPavingVarCovar();
				RealVec varcov1 = myHistFirst.getRootPavingVarCovar();
				assert(varcov.size() == d1*d1);
				assert(!checkAllNaN(varcov));
				assert( checkSame (varcov, varcov1, temp.getRootCounter() ) );
				cout << "Passed asserts that variance-covariance is same as for hist1" << endl;
				
				string s = "holdAllStatsResetToTrueForCopyHist1.txt";
				outputADH(s, temp);
			}
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reset holdAllStats for copy of hist1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nreset holdAllStats for copy of hist4" << endl;
			
			AdaptiveHistogram temp(myHistFourth);
				
			{
				cout << "\nChange holdAllStats for copy of hist4 to true" << endl;
			
				bool newHoldAllStats = true;
				temp.setHoldAllStats(newHoldAllStats);
				assert(temp.getHoldAllStats() == newHoldAllStats);
				cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;
				assert(temp.getRootCounter() == myHistFourth.getRootCounter());
				cout << "Passed assert that root counter is unchanged" << endl;
				rvector mean = temp.getRootPavingMean();
				assert(VecLen(mean) == d1);
				assert(!checkAllNaN(mean));
				rvector mean4 = checkMean(theData4);
				assert( checkSame( mean, mean4, temp.getRootCounter() ) );
				cout << "Passed asserts that mean is as calculated by check calcs" << endl;
				RealVec varcov = temp.getRootPavingVarCovar();
				assert(varcov.size() == d1*d1);
				assert(!checkAllNaN(varcov));
				RealVec varcov4 = checkVarCov(theData4);
				assert( checkSame( varcov, varcov4, temp.getRootCounter() ) );
				
				cout << "Passed asserts that variance-covariance is as calculated by check calcs" << endl;
				
				string s = "holdAllStatsResetToTrueForCopyHist4.txt";
				outputADH(s, temp);
			}
			{
				cout << "\nChange holdAllStats for copy of hist4 to back to false" << endl;
			
				bool newHoldAllStats = false;
				temp.setHoldAllStats(newHoldAllStats);
				assert(temp.getHoldAllStats() == newHoldAllStats);
				cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;
				assert(temp.getRootCounter() == myHistFourth.getRootCounter());
				cout << "Passed assert that root counter is unchanged" << endl;
				rvector mean = temp.getRootPavingMean();
				assert(VecLen(mean) == d1);
				assert(checkAllNaN(mean));
				cout << "Passed asserts that mean is not held" << endl;
				RealVec varcov = temp.getRootPavingVarCovar();
				assert(varcov.size() == d1*d1);
				assert(checkAllNaN(varcov));
				cout << "Passed asserts that variance-covariance is not held" << endl;
				
				string s = "holdAllStatsResetToFalseForCopyHist4.txt";
				outputADH(s, temp);
			}
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reset holdAllStats for copy of hist4:\n" << msg << endl;
			throw;
		}
		
		//swap
						
		try {
			
			cout << "\nswap copies of hist1 and hist2" << endl;
			
			AdaptiveHistogram temp1(myHistFirst);	
			AdaptiveHistogram temp2(myHistSecond);	
			string swapcheck11("copyOfHist1detailsSwap1.txt");
			swapCheckOutput(swapcheck11, temp1);
			string swapcheck12("copyOfHist2detailsSwap1.txt");
			swapCheckOutput(swapcheck12, temp2);
			cout << "Extra details for copy of hist1 before swap are in " << swapcheck11 << endl;
			cout << "Extra details for copy of hist2 before swap are in " << swapcheck12 << endl;
			cout << endl;
		
			std::swap(temp1, temp2);
			
			assert( temp1.getRootCounter() == myHistSecond.getRootCounter() );
			assert( temp2.getRootCounter() == myHistFirst.getRootCounter() );
			
			//temp 1 mean and varcov should all be NaNs
			{
				rvector mean = temp1.getRootPavingMean();
				RealVec varcov = temp1.getRootPavingVarCovar();
				assert( checkAllNaN(mean) );
				assert (checkAllNaN( varcov ) );
			}
			//temp 2 mean and varcov should all be same as for hist 1
			{
				rvector mean = temp2.getRootPavingMean();
				RealVec varcov = temp2.getRootPavingVarCovar();
				rvector mean1 = myHistFirst.getRootPavingMean();
				RealVec varcov1 = myHistFirst.getRootPavingVarCovar();
				assert( checkSame(mean, mean1, temp1.getRootCounter() ) );
				assert( checkSame(varcov, varcov1, temp1.getRootCounter() ) );
			}
			
			cout << "Passed asserts that root counters and mean and variance-covariance of originals correspond to swaps" << endl;

			string s1 = "shouldNowBeCopyOfHist1Swap1.txt";
			outputADH(s1, temp2);
			
			string s2 = "shouldNowBeCopyOfHist2Swap1.txt";
			outputADH(s2, temp1);
			
			string swapcheck21 = "shouldNowBeCopyOfHist1detailsSwap1.txt";
			swapCheckOutput(swapcheck21, temp2);
			string swapcheck22 = "shouldNowBeCopyOfHist2detailsSwap1.txt";
			swapCheckOutput(swapcheck22, temp1);
			cout << "Extra details for should now be copy of hist1 after swap are in " << swapcheck21 << endl;
			cout << "Extra details for should now be copy of hist2 after swap are in " << swapcheck22 << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do swap of copies of hist1 and hist2 pavings:\n" << msg << endl;
			throw;
		}
		try {
			
			cout << "\nswap copies of hist2 and hist4" << endl;
			
			AdaptiveHistogram temp1(myHistSecond);	
			AdaptiveHistogram temp2(myHistFourth);	
			string swapcheck11("copyOfHist2detailsSwap2.txt");
			swapCheckOutput(swapcheck11, temp1);
			string swapcheck12("copyOfHist4detailsSwap2.txt");
			swapCheckOutput(swapcheck12, temp2);
			cout << "Extra details for copy of hist2 before swap are in " << swapcheck11 << endl;
			cout << "Extra details for copy of hist4 before swap are in " << swapcheck12 << endl;
			cout << endl;
		
			std::swap(temp1, temp2);
			
			assert( temp1.getRootCounter() == myHistFourth.getRootCounter() );
			assert( temp2.getRootCounter() == myHistSecond.getRootCounter() );
			
			//temp 1 mean and varcov should all be NaNs
			{
				rvector mean = temp1.getRootPavingMean();
				RealVec varcov = temp1.getRootPavingVarCovar();
				assert( checkAllNaN(mean) );
				assert (checkAllNaN( varcov ) );
			}
			//temp 2 mean and varcov should all be NaNs
			{
				rvector mean = temp2.getRootPavingMean();
				RealVec varcov = temp2.getRootPavingVarCovar();
				assert( checkAllNaN(mean) );
				assert (checkAllNaN( varcov ) );
			}
			
			cout << "Passed asserts that root counters and mean and variance-covariance of originals correspond to swaps" << endl;

			string s1 = "shouldNowBeCopyOfHist2Swap2.txt";
			outputADH(s1, temp2);
			
			string s2 = "shouldNowBeCopyOfHist4Swap2.txt";
			outputADH(s2, temp1);
			
			string swapcheck21 = "shouldNowBeCopyOfHist2detailsSwap2.txt";
			swapCheckOutput(swapcheck21, temp2);
			string swapcheck22 = "shouldNowBeCopyOfHist4detailsSwap2.txt";
			swapCheckOutput(swapcheck22, temp1);
			cout << "Extra details for should now be copy of hist2 after swap are in " << swapcheck21 << endl;
			cout << "Extra details for should now be copy of hist4 after swap are in " << swapcheck22 << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do swap of copies of hist2 and hist4 pavings:\n" << msg << endl;
			throw;
		}
		try {
			
			cout << "\nL1 distance with hist with no paving (should fail)" << endl;
			
			const AdaptiveHistogram temp1;
			const AdaptiveHistogram temp2;
			
			cxsc::real dis = temp1.getL1Distance(temp2);
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance with hist with no paving:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance from copy of hist 1 to hist with no paving (should fail)" << endl;
			
			const AdaptiveHistogram temp1(myHistFirst);
			const AdaptiveHistogram temp2;
			
			cxsc::real dis = temp1.getL1Distance(temp2);
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance from copy of hist 1 to hist with no paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance from copy of hist 1 to hist with incorrect paving box (this should fail)" << endl;
			
			const AdaptiveHistogram temp1(myHistFirst);
			
			ivector pavingBoxT1(d1);
		
			for(int k=1; k < d1; k++) pavingBoxT1[k] = pavingInterval1;
			
			interval pavingIntervalT1(-5,4);
			pavingBoxT1[d1] = pavingIntervalT1;
			
			const AdaptiveHistogram temp2(pavingBoxT1);
			
			cxsc::real dis = temp1.getL1Distance(temp2);
			throw std::logic_error("Should not be able to do this");
			
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance from copy of hist 1 to hist with incorrect paving box:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance from copy of hist 1 to hist with correct paving box but no data" << endl;
			
			const AdaptiveHistogram temp1(myHistFirst);
			
			const AdaptiveHistogram temp2(pavingBox1);
			
			cxsc::real dis = temp1.getL1Distance(temp2);
			
			cxsc::real shouldBe(1.0);
			assert(dis == shouldBe);
			cout << "Passed assert that distance from copy of hist 1 to hist with correct paving box but no data is " << shouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance from copy of hist 1 to hist with correct paving box but no data:\n" << msg << endl;
			throw;
		}
		
		try {
			
			cout << "\nL1 distance from copy of hist 1 to itself (should be 0)" << endl;
			
			const AdaptiveHistogram temp1(myHistFirst);
			//const AdaptiveHistogram temp2;
			
			cxsc::real dis = temp1.getL1Distance(temp1);
			
			cxsc::real shouldBe(0.0);
			assert(dis == shouldBe);
			cout << "Passed assert that distance between hist 1 and itself is " << shouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance from copy of hist 1 to itself:\n" << msg << endl;
			throw;
		}
		
		cxsc::real dis2_3(0.0);
		try {
			
			cout << "\nL1 distance between copies of hist2 and hist3" << endl;
			
			const AdaptiveHistogram temp1(myHistSecond);
			const AdaptiveHistogram temp2(myHistThird);
			
			dis2_3 = temp1.getL1Distance(temp2);
			
			cxsc::real shouldBe(10.0/10.0);
			assert(dis2_3 == shouldBe);
			cout << "Passed assert that distance between hist 2 and 3 is " << shouldBe << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance between copies of hist2 and hist3:\n" << msg << endl;
			throw;
		}
		
		try {
			
			cout << "\nL1 distance between copies of hist3 and hist2" << endl;
			
			const AdaptiveHistogram temp1(myHistThird);
			const AdaptiveHistogram temp2(myHistSecond);
			
			cxsc::real dis = temp1.getL1Distance(temp2);
			
			assert(dis == dis2_3);
			cout << "Passed assert that distance between hist 3 and 2 is same as distance between hist 2 and 3" << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance between copies of hist2 and hist3:\n" << msg << endl;
			throw;
		}
		
		
		try {
		
			cout << "\nCoverage with hist with no paving (should fail)" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -6;
			pt[2] = 0;
			
			AdaptiveHistogram temp;
			
			double cov = temp.findCoverage(pt);
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage with hist with no paving:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nDensity with hist with no paving" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -6;
			pt[2] = 0;
			
			AdaptiveHistogram temp;
			
			double ed = temp.findEmpiricalDensity(pt);
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density with hist with no paving:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage with hist1 with 3-d point (should fail)" << endl;
			
			cxsc::rvector pt(3);
			pt[1] = -6;
			pt[2] = 0;
			pt[2] = 1;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			
			double cov = myHistFirst.findCoverage(pt);
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage with hist with 3-d point:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nDensity with hist1 with 3-d point (should fail)" << endl;
			
			cxsc::rvector pt(3);
			pt[1] = -6;
			pt[2] = 0;
			pt[2] = 1;
			
			AdaptiveHistogram temp;
			
			double ed = myHistFirst.findEmpiricalDensity(pt);
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density with hist with 3-d point:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nCoverage and density with histogram with no data (should be 0)" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = 0;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			
			AdaptiveHistogram temp(pavingBox1);
			
			double cov = temp.findCoverage(pt);
			double ed = temp.findEmpiricalDensity(pt);
			
			double covShouldBe = 0.0;
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			
			double edShouldBe = 0.0;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density with hist with no data:\n" << msg << endl;
			throw;
		}
		
		try {
		
			cout << "\nCoverage and density for point not in hist 1 (should be 0)" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -6;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = 0.0;
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			
			double edShouldBe = 0.0;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point not in hist1:\n" << msg << endl;
			throw;
		}
		
		try {
		
			cout << "\nCoverage and density for point in hist 1" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -2.5;
			pt[2] = -2.5; // in XLLRR
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = (10.0/10);
			double edShouldBe = (4.0/(10*6.25));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point in hist1:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for point in hist 1" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = 0;
			pt[2] = 0; // in XRRLL
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = ((10.0-4.0)/10);
			double edShouldBe = (3.0/(10*6.25));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point in hist1:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for point in hist 1" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -5;
			pt[2] = 0; // in XLR
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = ((10.0-4.0-3.0)/10);
			double edShouldBe = (2.0/(10*25.0));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point in hist1:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for point in hist 1" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = 5;
			pt[2] = -5; // in XRL
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = ((10.0-4.0-3.0-2.0)/10);
			double edShouldBe = (1.0/(10*25.0));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point in hist1:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for point in hist 1" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -5;
			pt[2] = -5; // in XLLL
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = myHistFirst.findCoverage(pt);
			double ed = myHistFirst.findEmpiricalDensity(pt);
			
			double covShouldBe = ((10.0-4.0-3.0-2.0-1.0)/10);
			double edShouldBe = (0.0/(10*12.5));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for point in hist1:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for hist with 3 leaves all with same density" << endl;
			
			AdaptiveHistogram temp(pavingBox1);
			temp.splitToShape("2,2,1");

			// put in the data in a 'pulse' with no further splitting
			RVecData data;
			data = getData1(data);
			bool successfulInsertion = temp.insertFromRVec(data);
			if (successfulInsertion) {

			
				cxsc::rvector pt(d1);
				pt[1] = -5;
				pt[2] = -5; // in XLL
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
			
				double cov = temp.findCoverage(pt);
				double ed = temp.findEmpiricalDensity(pt);
				
				/* although pt is in XLL, which in the 
				 * ordering of the leaves by height happens
				 * to come last, coverage should use first
				 * histogram element of same height as XLL
				 * so in this case coverage is 1, ie any 
				 * chance ordering amongst leaves of the same height
				 * should not make difference to the answer*/
				
				double covShouldBe = ((4.0)/4);
				double edShouldBe = (1.0/(4*25.0));
				
				assert(cov == covShouldBe);
				cout << "Passed assert that coverage is " << covShouldBe << endl;
				assert(ed == edShouldBe);
				cout << "Passed assert that empirical density is " << edShouldBe << endl;
			}
			else cout << "data insertion unsuccessful" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for hist with 3 leaves all with same density:\n" << msg << endl;
			throw;
		}
		try {
		
			cout << "\nCoverage and density for hist with 3 leaves, smallest two with same density" << endl;
			
			AdaptiveHistogram temp(pavingBox1);
			temp.splitToShape("2,2,1");

			// put in the data in a 'pulse' with no further splitting
			RVecData data;
			data = getData2(data);
			bool successfulInsertion = temp.insertFromRVec(data);
			if (successfulInsertion) {

			
				cxsc::rvector pt(d1);
				pt[1] = -5;
				pt[2] = -5; // in XLL
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
			
				double cov = temp.findCoverage(pt);
				double ed = temp.findEmpiricalDensity(pt);
				
				/* although pt is in XLL, which in the 
				 * ordering of the leaves by height happens
				 * to come last, coverage should use first
				 * histogram element of same height as XLL,
				 * which is XR in the ordering used, 
				 * so in this case coverage only skips
				 * over the tallest node XLR (2 data points), ie any 
				 * chance ordering amongst leaves of the same height
				 * should not make difference to the answer*/
				
				
				double covShouldBe = ((5.0-2.0)/5);
				double edShouldBe = (1.0/(5*25.0));
				
				assert(cov == covShouldBe);
				cout << "Passed assert that coverage is " << covShouldBe << endl;
				assert(ed == edShouldBe);
				cout << "Passed assert that empirical density is " << edShouldBe << endl;
			}
			else cout << "data insertion unsuccessful" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage and density for hist with 3 leaves, smallest two with same density:\n" << msg << endl;
			throw;
		}
		
		cout << "\nEnd of test\n" << endl;

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed somewhere with histograms with data:\n" << msg << endl;
	}
	

	cout << endl;
	cout << "Test hist arithmetic " << endl;
	testHistArithmetic();
	

    return 0;
	

} // end of test program





RVecData& getData1(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getData2(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -1;
        thisrv[2] = 1;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getData3(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -6;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -6;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	
	return data;
}

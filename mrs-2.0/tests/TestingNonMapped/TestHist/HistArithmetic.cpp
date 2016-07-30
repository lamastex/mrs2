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
\brief Testing StatsSubPavings (aka SPSnodes) with histogram arithmetic
 */

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <cassert>
#include <stdexcept>

#include "testing_tools.hpp"
#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data

using namespace cxsc;
using namespace std;
using namespace subpavings;

void testHistArithmetic()
{

    // ------- prepare to generate some data for the tests -----------

    // set up a random number generator
    const gsl_rng_type * T;
    gsl_rng * r;

    double sigma_x=1;   // distribution parameter for Biv Gaussian
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    const int n=10;    // number to generate for hists 1, 2, 3
	const int n4 = 7; // for histogram 4
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
	for (int i = 0; i < n4; i++) { // note n4 not n

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

    AdaptiveHistogram myHistFirst(pavingBox1, true); // hold all stats
    myHistFirst.splitToShape("3,4,4,2,2,4,4,3");

    AdaptiveHistogram myHistSecond(pavingBox1);
    //myHistSecond.splitToShape("1,2,3,3");
    myHistSecond.splitToShape("2,3,4,4,2,3,3");

    AdaptiveHistogram myHistThird(pavingBox1);
    myHistThird.splitToShape("1,2,3,3");
	
	AdaptiveHistogram myHistFourth(pavingBox1); 
    myHistFourth.splitToShape("3,4,4,2,2,4,4,3");

    // put in the data in a 'pulse' with no further splitting
    bool successfulInsertionHistFirst = myHistFirst.insertFromRVec(theData1);
    if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;

	bool successfulInsertionHistSecond = myHistSecond.insertFromRVec(theData2);
    if (!successfulInsertionHistSecond) cout << "unsuccessful insertion 2" << endl;

    bool successfulInsertionHistThird = myHistThird.insertFromRVec(theData3);
    if (!successfulInsertionHistThird) cout << "unsuccessful insertion 3" << endl;
	
	bool successfulInsertionHistFourth = myHistFourth.insertFromRVec(theData4); 
    if (!successfulInsertionHistFourth) cout << "unsuccessful insertion 4" << endl;

    if (successfulInsertionHistFirst && successfulInsertionHistSecond 
					&& successfulInsertionHistThird && successfulInsertionHistFourth) {
						
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
		
		
		RealVec vcov1 = myHistFirst.getRootPavingVarCovar();
		
		assert(vcov1.size() == d1*d1);
		assert(!checkAllNaN(vcov1));
		RealVec calcVCov1 = checkVarCov(theData1);
		assert(calcVCov1.size() == d1*d1);
		assert(checkSame(vcov1, calcVCov1, n) );
		cout << "Passed asserts that variance-covariance for histogram 1 is length " << d1*d1 << " calculation of variance-covariance matches check calc" << endl;

			
		cout << endl << endl;
		
		string outputFileName1 = "Hist1.txt";
        myHistFirst.outputToTxtTabs(outputFileName1);
		cout << "Histogram 1 is in " << outputFileName1 << endl;
		
		string outputFileName2 = "Hist2.txt";
        myHistSecond.outputToTxtTabs(outputFileName2);
		cout << "Histogram 2 is in " << outputFileName2 << endl;
		
		string outputFileName3 = "Hist3.txt";
        myHistThird.outputToTxtTabs(outputFileName3);
		cout << "Histogram 3 is in " << outputFileName3 << endl;
		
		string outputFileName4 = "Hist4.txt";
        myHistFourth.outputToTxtTabs(outputFileName4);
		cout << "Histogram 4 is in " << outputFileName4 << endl;

        string outputFileName1c = "Hist1Check.txt";
        swapCheckOutput(outputFileName1c, *(myHistFirst.getSubPaving()));
		cout << "Histogram 1 paving details are in " << outputFileName1c << endl;
		
		string outputFileName2c = "Hist2Check.txt";
        swapCheckOutput(outputFileName2c, *(myHistSecond.getSubPaving()));
		cout << "Histogram 2 paving details are in " << outputFileName2c << endl;
		
		string outputFileName3c = "Hist3Check.txt";
       swapCheckOutput(outputFileName3c, *(myHistThird.getSubPaving()));
		cout << "Histogram 3 paving details are in " << outputFileName3c << endl;
		
		string outputFileName4c = "Hist4Check.txt";
        swapCheckOutput(outputFileName4c, *(myHistFourth.getSubPaving()));
		cout << "Histogram 4 paving details are in " << outputFileName4c << endl;
		
		try {
			cout << "\nTry to add histograms with incompatible dimensions (this should fail)" << endl;
		
			// make a box with the wrong dimensions
			int dT1 = 3; // dimension of the box to sample data from
			ivector pavingBoxT1(dT1);
			interval pavingIntervalT1(-5,5);
			for(int k=1; k <= dT1; k++) pavingBoxT1[k] = pavingIntervalT1;
		
			AdaptiveHistogram HistT1(pavingBoxT1);
			AdaptiveHistogram newHist = myHistFirst + HistT1;
			
			throw std::logic_error("Should not be able to do constructor this");
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add histograms with incompatible dimensions:\n" << msg << endl;
		}
		
		try {
			cout << "\nTry to add histograms with incompatible box side lengths (this should fail)" << endl;
		
			//make a box of wrong size
			int dT2 = 2; // dimension of the box to sample data from
			ivector pavingBoxT2(dT2);
			interval pavingIntervalT2(-5,6);
			for(int k=1; k <= dT2; k++) pavingBoxT2[k] = pavingIntervalT2;
		
			AdaptiveHistogram HistT2(pavingBoxT2);
			AdaptiveHistogram newHist = HistT2 + myHistFirst;
			
			throw std::logic_error("Should not be able to do constructor with default histogram");
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add histograms with incompatible box lengths:\n" << msg << endl;
		}

		try {
			cout << "\nTry to add two histograms with no subpavings" << endl;
			
			AdaptiveHistogram myHistTp1;
			AdaptiveHistogram myHistTp2;
			AdaptiveHistogram temp = myHistTp1 + myHistTp2;
			string s = "AdditionBothNoPaving.txt";
			temp.outputToTxtTabs(s);
			cout << "\nResults of addition of two histograms with no subpavings are in " 
					<< s << endl;
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
				
			
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add two histograms with no subpavings:\n" << msg << endl;
			throw;
		}
		
        try {
			cout << "\nTry to add hist with no subpaving and hist 1" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp = myHistTp + myHistFirst;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT1 = "AdditionNoPavingAndHist1.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of addition of no subpaving and Histogram 1 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add no subpaving and Histogram 1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist1 and hist with no subpaving" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp = myHistFirst + myHistTp;
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1,temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT2 = "AdditionHist1AndNoPaving.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of addition of Histogram 1 and no subpaving are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 1 and no subpaving:\n" << msg << endl;
			throw;
		}
		
		
		try {
			cout << "\nTry to add hist with no subpaving and hist 4" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp = myHistTp + myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN (temp.getRootPavingMean() ) );
			assert( checkAllNaN (temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats not held" << endl;

			string outputFileNameT1 = "AdditionNoPavingAndHist4.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of addition of no subpaving and Histogram 4 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add no subpaving and Histogram 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist4 and hist with no subpaving" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp = myHistFourth + myHistTp;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN (temp.getRootPavingMean() ) );
			assert( checkAllNaN (temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats not held" << endl;

			string outputFileNameT2 = "AdditionHist4AndNoPaving.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of addition of Histogram 4 and no subpaving are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 4 and no subpaving:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist with no subpaving but holdAllStats = 1, and hist 4" << endl;
			
			AdaptiveHistogram myHistTp(true);  // holdAllStats = true;
			AdaptiveHistogram temp = myHistTp + myHistFourth;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT1 = "AdditionNoPavingHoldAllStatsTrueAndHist4.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of addition of no subpaving with no subpaving but holdAllStats = 1, and Histogram 4 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add no subpaving but holdAllStats = 1, and Histogram 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist4 and hist with no subpaving but holdAllStats = 1" << endl;
			
			AdaptiveHistogram myHistTp(true);
			AdaptiveHistogram temp = myHistFourth + myHistTp;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT2 = "AdditionHist4AndNoPavingHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of addition of Histogram 4 and no subpaving but holdAllStats = 1are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 4 and no subpaving but holdAllStats = 1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist with no subpaving and label 1, and hist 1" << endl;
			
			AdaptiveHistogram myHistTp (false, 1); // hold allstats false, label 1
			AdaptiveHistogram temp = myHistTp + myHistFirst;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert( temp.getLabel() == 0);
			cout << "Passed assert that label of new his is 0" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add hist with no subpaving and label 1, and hist 1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist with no subpaving and label 1, and copy of hist 1 with label reset to 1" << endl;
			
			int label = 1;
			AdaptiveHistogram myHistTp1 (false, label); // hold allstats false, label 1
			AdaptiveHistogram myHistTp2 (myHistFirst); 
			myHistTp2.setLabel(label);
			AdaptiveHistogram temp = myHistTp1 + myHistTp2;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert( temp.getLabel() == label);
			cout << "Passed assert that label of new his is " << label << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add hist with no subpaving and label 1, and copy of hist 1 with label reset to 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry to add hist1 and hist with empty box" << endl;
			
			AdaptiveHistogram myHistTp4(pavingBox1);
			AdaptiveHistogram temp = myHistFirst + myHistTp4;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT4 = "AdditionHist1AndEmptyBox.txt";
			temp.outputToTxtTabs(outputFileNameT4);

			cout << "\nResults of addition of Histogram 1 and hist with empty box are in " 
					<< outputFileNameT4 << endl;
			
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 1 and hist with empty box:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist with empty box and hist 1" << endl;
			
			AdaptiveHistogram myHistTp4(pavingBox1);
			AdaptiveHistogram temp = myHistTp4 + myHistFirst;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT4 = "AdditionEmptyBoxAndHist1.txt";
			temp.outputToTxtTabs(outputFileNameT4);

			cout << "\nResults of addition of hist with empty box and Histogram 1 are in " 
					<< outputFileNameT4 << endl;
			
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 1 and hist with empty box:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist4 and hist with empty box" << endl;
			
			AdaptiveHistogram myHistTp4(pavingBox1);
			AdaptiveHistogram temp = myHistFourth + myHistTp4;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert(checkAllNaN( temp.getRootPavingMean() ));
			assert ( checkAllNaN( temp.getRootPavingVarCovar()) ); 
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats not held" << endl;

			string outputFileNameT4 = "AdditionHist4AndEmptyBox.txt";
			temp.outputToTxtTabs(outputFileNameT4);

			cout << "\nResults of addition of Histogram 4 and hist with empty box are in " 
					<< outputFileNameT4 << endl;
			
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 4 and hist with empty box:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist with empty box and hist 4" << endl;
			
			AdaptiveHistogram myHistTp4(pavingBox1);
			AdaptiveHistogram temp = myHistTp4 + myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert(checkAllNaN( temp.getRootPavingMean() ));
			assert ( checkAllNaN( temp.getRootPavingVarCovar()) ); 
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats not held" << endl;

			string outputFileNameT4 = "AdditionEmptyBoxAndHist4.txt";
			temp.outputToTxtTabs(outputFileNameT4);

			cout << "\nResults of addition of hist with empty box and Histogram 4 are in " 
					<< outputFileNameT4 << endl;
			
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 4 and hist with empty box:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry to add hist with empty box but holdAllStats = 1, and hist 4" << endl;
			
			AdaptiveHistogram myHistTp(pavingBox1, true);  // holdAllStats = true;
			AdaptiveHistogram temp = myHistTp + myHistFourth;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT1 = "AdditionEmptyBoxAllStatsTrueAndHist4.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of addition of hist with empty box with no subpaving but holdAllStats = 1, and Histogram 4 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add hist with with empty box but holdAllStats = 1, and Histogram 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to add hist4 and hist with empty box but holdAllStats = 1" << endl;
			
			AdaptiveHistogram myHistTp(pavingBox1, true);
			AdaptiveHistogram temp = myHistFourth + myHistTp;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT2 = "AdditionHist4AndNoPavingHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of addition of Histogram 4 with empty box but holdAllStats = 1are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histogram 4 with emptyBox but holdAllStats = 1:\n" << msg << endl;
			throw;
		}
		
		// try +=
		try {
			cout << "\nTry to hist with no subpaving += hist 1 ( result should not hold stats)" << endl;
			
			AdaptiveHistogram temp;
			temp += myHistFirst;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFirst but stats not held" << endl;

			string outputFileNameT1 = "NoPavingPlusEqualsHist1.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist with no subpaving += hist 1 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with no subpaving += hist 1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to do hist1 += hist with no subpaving (should hold all stats" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp(myHistFirst);
			temp += myHistTp;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1,temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT2 = "Hist1PlusEqualsHistNoPaving.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of hist1 += hist with no subpaving are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist1 += hist with no subpaving:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry to hist with no subpaving and label 1 += hist 1 ( result should have label 1)" << endl;
			
			int label = 1;
			AdaptiveHistogram temp (false, label);
			temp += myHistFirst;
			
			assert( temp.getLabel() == label);
			cout << "Passed assert that label is " << label << endl;	
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFirst but stats not held" << endl;
			
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with no subpaving and label 1 += hist 1:\n" << msg << endl;
			throw;
		}
				
		try {
			cout << "\nTry to do hist with no subpaving += hist 4" << endl;
			
			AdaptiveHistogram temp;
			temp += myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFourth but stats not held" << endl;

			string outputFileNameT1 = "NoPavingPlusEqualsHist4.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist with no subpaving += hist 4 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with no subpaving += hist 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to do hist4 += hist with no subpaving" << endl;
			
			AdaptiveHistogram myHistTp;
			AdaptiveHistogram temp(myHistFourth);
			temp += myHistTp;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFourth but stats not held" << endl;

			string outputFileNameT1 = "Hist4PlusEqualsNoPaving.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist 4 += hist with no paving are in " 
					<< outputFileNameT1 << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist4 += hist with no subpaving:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry do hist4 += hist with no subpaving but holdAllStats = 1" << endl;
			
			AdaptiveHistogram myHistTp(true);
			AdaptiveHistogram temp(myHistFourth);
			temp += myHistTp;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert(checkAllNaN( temp.getRootPavingMean() ) );
			assert(checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats are not held" << endl;

			string outputFileNameT2 = "Hist4PlusEqualsNoPavingHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of Histogram 4 += no subpaving but holdAllStats = 1 are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 4 += no subpaving but holdAllStats = 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry do hist with no subpaving but holdAllStats = 1, += hist4" << endl;
			
			AdaptiveHistogram temp(true);
			temp += myHistFourth;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT2 = "NoPavingHoldAllStatsTruePlusEqualsHist4.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of hist with no subpaving but holdAllStats = 1, += hist4 are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with no subpaving but holdAllStats = 1, += hist4:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry to hist with empty box += hist 1 ( result should not hold stats)" << endl;
			
			AdaptiveHistogram temp(pavingBox1);
			temp += myHistFirst;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFirst but stats not held" << endl;

			string outputFileNameT1 = "EmptyBoxPlusEqualsHist1.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist with empty box += hist 1 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with empty box += hist 1:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to do hist1 += hist with empty box (should hold all stats" << endl;
			
			AdaptiveHistogram myHistTp(pavingBox1);
			AdaptiveHistogram temp(myHistFirst);
			temp += myHistTp;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() );
			assert( checkSame( temp.getRootPavingMean(), mean1, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, vcov1,temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist counter, mean and var covar is same as for myHistFirst" << endl;

			string outputFileNameT2 = "Hist1PlusEqualsHistEmptyBox.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of hist1 += hist with empty box are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist1 += hist with empty box:\n" << msg << endl;
			throw;
		}
		
		
		try {
			cout << "\nTry to do hist with empty box += hist 4" << endl;
			
			AdaptiveHistogram temp(pavingBox1);
			temp += myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFourth but stats not held" << endl;

			string outputFileNameT1 = "EmptyBoxPlusEqualsHist4.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist with empty box += hist 4 are in " 
					<< outputFileNameT1 << endl;
		}
		
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with empty box += hist 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry to do hist4 += hist with empty box" << endl;
			
			AdaptiveHistogram myHistTp(pavingBox1);
			AdaptiveHistogram temp(myHistFourth);
			temp += myHistTp;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist counter same as for myHistFourth but stats not held" << endl;

			string outputFileNameT1 = "Hist4PlusEqualsEmptyBox.txt";
			temp.outputToTxtTabs(outputFileNameT1);

			cout << "\nResults of hist 4 += hist with empty box are in " 
					<< outputFileNameT1 << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist4 += hist with empty box:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry do hist4 += hist with empty box but holdAllStats = 1" << endl;
			
			AdaptiveHistogram myHistTp(pavingBox1, true);
			AdaptiveHistogram temp(myHistFourth);
			temp += myHistTp;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			assert(checkAllNaN( temp.getRootPavingMean() ) );
			assert(checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats are not held" << endl;

			string outputFileNameT2 = "Hist4PlusEqualsEmptyBoxHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of Histogram 4 += empty box but holdAllStats = 1 are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 4 += empty box but holdAllStats = 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry do hist with empty box but holdAllStats = 1, += hist4" << endl;
			
			AdaptiveHistogram temp(pavingBox1, true);
			temp += myHistFourth;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() );
			rvector mean = temp.getRootPavingMean();
			RealVec varcov = temp.getRootPavingVarCovar();
			rvector mean4 = checkMean(theData4);
			RealVec varcov4 = checkVarCov(theData4);
			assert( checkSame (mean, mean4, temp.getRootCounter() ) );
			assert( checkSame (varcov, varcov4, temp.getRootCounter() ) );
			cout << "Passed asserts that new hist counter is same as for myHistFourth and stats match check calculations" << endl;

			string outputFileNameT2 = "EmptyBoxHoldAllStatsTruePlusEqualsHist4.txt";
			temp.outputToTxtTabs(outputFileNameT2);

			cout << "\nResults of hist with empty box but holdAllStats = 1, += hist4 are in " 
					<< outputFileNameT2 << endl;
        }
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do hist with empty box but holdAllStats = 1, += hist4:\n" << msg << endl;
			throw;
		}
			
		try {
			cout << "\nTry copyOfHist1 += Hist2 (should have stats since Hist1 has stats)" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
				
			temp += myHistSecond;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist2" << endl;
			RVecData combData = combineData(theData1, theData2);
			cxsc::rvector checkMn = checkMean(combData);
			assert( checkSame( temp.getRootPavingMean(), checkMn, temp.getRootCounter() ) );
			RealVec checkVc = checkVarCov(combData);
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, checkVc, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist mean and var covar are the same as for check calculations" << endl;

			string outputFileName = "CopyHist1PlusEqualHist2.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copyOfHist1 += Hist2 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 1 += Histogram 2:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry copyOfHist2 += Hist1 (this should not keep stats)" << endl;
			
			AdaptiveHistogram temp(myHistSecond);
				
			temp += myHistFirst;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist2" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist mean and var covar are not held" << endl;

			string outputFileName = "CopyHist2PlusEqualHist1.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copyOfHist2 += Hist1 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 2 +=  Histogram1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copyOfHist2 += copy of Hist1 and then set holdAllStats to true" << endl;
				
			AdaptiveHistogram temp(myHistSecond);
			AdaptiveHistogram temp1(myHistFirst);
			
			temp += temp1;
			
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFirst.getRootCounter());
			
			bool newHoldAllStats = true;
			temp.setHoldAllStats(newHoldAllStats);
			assert(temp.getHoldAllStats() == newHoldAllStats);
			cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;

			string s1 = "doubleCheckHist2PlusEqualsHist1AfterResetHoldAllStatsToTrue.txt";
			doubleCheckOutput(s1, temp);
			
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFirst.getRootCounter());
			cout << "Passed assert that root counter is unchanged" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data1and2 = combineData(theData1, theData2);
			rvector mean1and2 = checkMean(data1and2);
			RealVec varcov1and2 = checkVarCov(data1and2);
			assert( checkSame (mean, mean1and2, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov1and2, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			string outputFileName = "Hist2PlusEqualsHist1ResetHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of Histograms 2 plus equals hist 1 and then reset holdAllStatsare true " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histograms 2 plus equals hist 1 and reset holdAllStats to true:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copyOfHist1 += Hist4 (this should keep stats)" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
				
			temp += myHistFourth;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistFourth.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist4" << endl;
			RVecData combData = combineData(theData1, theData4);
			cxsc::rvector checkMn = checkMean(combData);
			assert( checkSame( temp.getRootPavingMean(), checkMn, temp.getRootCounter() ) );
			RealVec checkVc = checkVarCov(combData);
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, checkVc, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist mean and var covar are the same as for check calculations" << endl;

			
			string outputFileName = "CopyHist1PlusEqualHist4.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copyOfHist1 += Hist4 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 1 += Histogram 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry copyOfHist4 += Hist1 (this should not keep stats)" << endl;
			
			AdaptiveHistogram temp(myHistFourth);
				
			temp += myHistFirst;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistFourth.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist4" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkAllNaN( varcov ) ); 
			cout << "Passed asserts that new hist mean and var covar are not held" << endl;

			string outputFileName = "CopyHist4PlusEqualHist1.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copyOfHist4 += Hist1 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Histogram 4 += Histogram 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copy of Hist1 with label 1 += Hist4 (this should have label 1)" << endl;
			
			AdaptiveHistogram temp(myHistFirst);
			int label = 1;
			temp.setLabel(label);	
			temp += myHistFourth;
			
			assert( temp.getLabel() == label);
			cout << "Passed assert that label is " << label << endl;	
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to copy of Hist1 with label 1 += Hist4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry copyOfHist4 with label 2 += copy of Hist1 with label 1 (this should have label 2)" << endl;
			
			AdaptiveHistogram temp4(myHistFourth);
			AdaptiveHistogram temp1(myHistFirst);
			
			int label4 = 2;
			int label1 = 1;
			temp4.setLabel(label4);	
			temp1.setLabel(label1);	
			temp4 += temp1;
			
			assert( temp4.getLabel() == label4);
			cout << "Passed assert that label is " << label4 << endl;	
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to copyOfHist4 with label 2 += copy of Hist1 with label 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry default constructor histogram += (copyOfHist1 += Hist2) (should not keep stats)" << endl;
			
			AdaptiveHistogram temp;
			AdaptiveHistogram newhist(myHistFirst);
			newhist += myHistSecond;
			temp+=newhist;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter() );
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is sum of counter for hists 1 and 2 and mean and var covar are not held" << endl;

			
			string outputFileName = "DefaultPlusEqualsCopyHist1PlusEqualHist2.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of default += (copyOfHist1 += Hist2) are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do default constructor histogram += (copyOfHist1 += Hist2):\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry addition hist1 + hist2 " << endl;
				
			AdaptiveHistogram temp = myHistFirst + myHistSecond;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist2" << endl;
			RVecData combData = combineData(theData1, theData2);
			cxsc::rvector checkMn = checkMean(combData);
			assert( checkSame( temp.getRootPavingMean(), checkMn, temp.getRootCounter() ) );
			RealVec checkVc = checkVarCov(combData);
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, checkVc, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist mean and var covar are the same as for check calculations" << endl;

			string outputFileName = "AdditionHist1AndHist2.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 1 and 2 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 1 and 2:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry addition hist2 + hist1 " << endl;
				
			AdaptiveHistogram temp = myHistSecond + myHistFirst;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter());
			cout << "Passed assert that new hist counter is counter Hist1 + counter Hist2" << endl;
			RVecData combData = combineData(theData1, theData2);
			cxsc::rvector checkMn = checkMean(combData);
			assert( checkSame( temp.getRootPavingMean(), checkMn, temp.getRootCounter() ) );
			RealVec checkVc = checkVarCov(combData);
			RealVec varcov = temp.getRootPavingVarCovar();
			assert ( checkSame( varcov, checkVc, temp.getRootCounter() ) ); 
			cout << "Passed asserts that new hist mean and var covar are the same as for check calculations" << endl;

			string outputFileName = "AdditionHist2AndHist1.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 2 and 1 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 2 and 1:\n" << msg << endl;
			throw;
		}
		
		
		try {
			cout << "\nTry addition hist2 + hist4 " << endl;
				
			AdaptiveHistogram temp = myHistSecond + myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFourth.getRootCounter());
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is sum of counter for hists 2 and 4 and stats are not held" << endl;
			
			string outputFileName = "AdditionHist2AndHist4.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 2 and 4 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 2 and 4:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry addition hist4 + hist2 " << endl;
				
			AdaptiveHistogram temp = myHistFourth + myHistSecond;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFourth.getRootCounter());
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that new hist counter is sum of counter for hists 2 and 4 and stats are not held" << endl;
			
			
			string outputFileName = "AdditionHist4AndHist2.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 4 and 2 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 4 and 2:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry addition hist2 + hist4 and then set holdAllStats to true" << endl;
				
			AdaptiveHistogram temp = myHistSecond + myHistFourth;
			
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFourth.getRootCounter());
			
			bool newHoldAllStats = true;
			temp.setHoldAllStats(newHoldAllStats);
			assert(temp.getHoldAllStats() == newHoldAllStats);
			cout << "Passed assert that getHoldAllStats() = " << newHoldAllStats << endl;

			string s1 = "doubleCheckHist1PlusHist2AfterResetHoldAllStatsToTrue.txt";
			doubleCheckOutput(s1, temp);
			
			assert(temp.getRootCounter() == myHistSecond.getRootCounter() + myHistFourth.getRootCounter());
			cout << "Passed assert that root counter is unchanged" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data2and4 = combineData(theData2, theData4);
			rvector mean2and4 = checkMean(data2and4);
			RealVec varcov2and4 = checkVarCov(data2and4);
			assert( checkSame (mean, mean2and4, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov2and4, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			string outputFileName = "AdditionHist2AndHist4ResetHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 2 and 4 and then reset holdAllStatsare true " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 2 and 4 and reset holdAllStats to true:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry addition hist3 + hist2 + hist1 " << endl;
				
			AdaptiveHistogram temp = myHistThird + myHistSecond + myHistFirst;
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data1and2and3 = combineData(theData1, theData2, theData3);
			rvector mean1and2and3 = checkMean(data1and2and3);
			RealVec varcov1and2and3 = checkVarCov(data1and2and3);
			assert( checkSame (mean, mean1and2and3, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov1and2and3, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			string outputFileName = "AdditionHist3AndHist2AndHist1.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 3 and 2 and 1 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 3 and 2 and 1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry addition (copy of hist1 with holdAllStats reset to false) + hist2 + hist3 " << endl;
			
			AdaptiveHistogram temp1(myHistFirst);
			temp1.setHoldAllStats(false);
			
			AdaptiveHistogram temp = temp1 + myHistSecond + myHistThird;
				
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that mean and variance-covariance are not held" << endl;
	

			string outputFileName = "AdditionHist1HoldAllStatsFalseAndHist2AndHist3.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 1 and 2 and 3 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add (copy of hist1 with holdAllStats reset to false) and hist2 and hist3:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry addition hist3 + hist2 + hist4 " << endl;
				
			AdaptiveHistogram temp = myHistThird + myHistSecond + myHistFourth;
			
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			
			assert(temp.getRootCounter() == myHistThird.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistFourth.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed asserts that mean and variance-covariance are not held" << endl;
	
			string outputFileName = "AdditionHist3AndHist2AndHist4.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 3 and 2 and 4 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 3 and 2 and 4:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry addition hist3 + hist2 + hist4 and then set hold all stats to true" << endl;
				
			AdaptiveHistogram temp = myHistThird + myHistSecond + myHistFourth;
			temp.setHoldAllStats(true);
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistThird.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistFourth.getRootCounter());
			cout << "Passed assert that root counter is unchanged" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data3and2and4 = combineData(theData3, theData2, theData4);
			rvector mean3and2and4 = checkMean(data3and2and4);
			RealVec varcov3and2and4 = checkVarCov(data3and2and4);
			assert( checkSame (mean, mean3and2and4, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov3and2and4, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			string outputFileName = "AdditionHist3AndHist2AndHist4ResetHoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of addition of Histograms 3 and 2 and 4 and reset holdAllStats to true are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to add Histograms 3 and 2 and 4 and reset holdAllStats to true:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copy of hist1 += copy of hist2 += copy of hist3 " << endl;
			
			AdaptiveHistogram temp(myHistFirst);
			AdaptiveHistogram temp2(myHistSecond);
			AdaptiveHistogram temp3(myHistThird);
				
			temp += temp2 += temp3;

			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
			
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data1and2and3 = combineData(theData1, theData2, theData3);
			rvector mean1and2and3 = checkMean(data1and2and3);
			RealVec varcov1and2and3 = checkVarCov(data1and2and3);
			assert( checkSame (mean, mean1and2and3, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov1and2and3, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			string outputFileName = "Hist1PlusEqualHist2PlusEqualHist3.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copy of hist1 += copy of hist2 += copy of hist3 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy of hist1 += copy of hist2 += copy of hist3:\n" << msg << endl;
			throw;
		}
		try {
			cout << "\nTry copy of hist3 += copy of hist2 += copy of hist1 " << endl;
			
			AdaptiveHistogram temp1(myHistFirst);
			AdaptiveHistogram temp2(myHistSecond);
			AdaptiveHistogram temp(myHistThird);
				
			temp += temp2 += temp1;

			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			assert(temp.getRootCounter() == myHistFirst.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed assert that stats are not held" << endl;
			
			string outputFileName = "Hist3PlusEqualHist2PlusEqualHist1.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copy of hist3 += copy of hist2 += copy of hist1 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy of hist2 += copy of hist2 += copy of hist1:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copy of hist4 += copy of hist2 += copy of hist3 " << endl;
			
			AdaptiveHistogram temp(myHistFourth);
			AdaptiveHistogram temp2(myHistSecond);
			AdaptiveHistogram temp3(myHistThird);
			
			temp += temp2 += temp3;
				
			assert( !temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is false" << endl;	
			
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			cout << "Passed assert that root counter is sum for three hists" << endl;
			assert( checkAllNaN( temp.getRootPavingMean() ) );
			assert( checkAllNaN( temp.getRootPavingVarCovar() ) );
			cout << "Passed assert that stats are not held" << endl;
			
			string outputFileName = "Hist4PlusEqualHist2PlusEqualHist3.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copy of hist4 += copy of hist2 += copy of hist3 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy of hist4 += copy of hist2 += copy of hist3:\n" << msg << endl;
			throw;
		}
		
		try {
			cout << "\nTry copy of hist4 += copy of hist2 += copy of hist3 and then reset holdAllStats to true" << endl;
			
			AdaptiveHistogram temp(myHistFourth);
			AdaptiveHistogram temp2(myHistSecond);
			AdaptiveHistogram temp3(myHistThird);
			
			temp += temp2 += temp3;
			assert(temp.getRootCounter() == myHistFourth.getRootCounter() + myHistSecond.getRootCounter()
												+ myHistThird.getRootCounter());
			
			temp.setHoldAllStats(true);
			
			assert( temp.getHoldAllStats());
			cout << "Passed assert that holdAllStats is true" << endl;	
				
			rvector mean = temp.getRootPavingMean();
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			RVecData data3and2and4 = combineData(theData3, theData2, theData4);
			rvector mean3and2and4 = checkMean(data3and2and4);
			RealVec varcov3and2and4 = checkVarCov(data3and2and4);
			assert( checkSame (mean, mean3and2and4, temp.getRootCounter() ) );
			RealVec varcov = temp.getRootPavingVarCovar();
			assert(varcov.size() == d1*d1);
			assert(!checkAllNaN(varcov));
			assert( checkSame (varcov, varcov3and2and4, temp.getRootCounter() ) );
			cout << "Passed asserts that mean and variance-covariance are as for check calcs" << endl;
	
			
			string outputFileName = "Hist4PlusEqualHist2PlusEqualHist3HoldAllStatsTrue.txt";
			temp.outputToTxtTabs(outputFileName);

			cout << "\nResults of copy of hist4 += copy of hist2 += copy of hist3 are in " 
						<< outputFileName << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy of hist4 += copy of hist2 += copy of hist3 and set holdAllStatsTrue:\n" << msg << endl;
			throw;
		}
		
		cout << endl;
		cout << "End of arithmetic tests\n\n" << endl;

    }
    else cout << "unsuccessful insertion" << endl;

} // end of test program

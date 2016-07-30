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
\brief Testing StatsSubPavings (aka SPSnodes)

Run the executable and then use the shell script testing_sps.sh
to run checks out output.

 */

#include "testing_tools.hpp"
#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data
#include "subpaving_exception.hpp"
#include "cxsc.hpp"

#include <fstream>  // input and output streams
#include <algorithm> // equal
#include <cassert>
#include <iomanip>


using namespace cxsc;
using namespace std;
using namespace subpavings;

std::ostream& swapCheck(std::ostream& os, const SPSnode& spn, int level = 0);
RVecData& getDataSinglePoint(RVecData& data);
RVecData& getDataTwoFiveDPoints(RVecData& data);
RVecData& getData2d(RVecData& data);
void checkNode(const string& s, const SPSnode * const spn);
		
int main()
{
	try {
		cout << "\nDefault constructed node" << endl;
		
		SPSnode temp;
		
		string s = "defaultNode.txt";
		outputSPS(s, temp);
		assert(temp.isEmpty());
		cout << "Passed assert that node is empty" << endl;;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default constructed SPSnode:\n" << msg << endl;
	}
	
	try {
		cout << "\nNode constructed with default box" << endl;
		
		cxsc::ivector uselessBox;
		SPSnode temp(uselessBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& ee) {
		cout << "\nFailed to construct node SPSnode with default box:\n" << ee.what() << endl;
	}

	try {
		
		cout << "\nConstruct node with box, split, and do test on nodes in tree" << endl;
		int d1 = 2; // dimension of the box to sample data from
		ivector pavingBox1(d1);
		interval pavingInterval1(-5,5);
		for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;

		AdaptiveHistogram temp(pavingBox1);
		temp.splitToShape("2,3,3,1");

		// put in the data in a 'pulse' with no further splitting
		RVecData data;
		data = getData2d(data);
		bool successfulInsertion = temp.insertFromRVec(data);
		if (successfulInsertion) {

			SPSnode* spn = temp.getSubPaving();
			
			string s = "testSPStree.txt";
			outputSPS(s, *spn);
			
			string scheck("testSPSnodeByNodeOutput.txt");
			checkNode(scheck, spn);
			cout << "node by node check output sent to file " << scheck << endl;
		}
		else cout << "data insertion unsuccessful" << endl;
		
	}
	catch (std::exception& ee) {
		cout << "\nFailed to do Construct node with box, split, and do test on nodes in tree:\n" << ee.what() << endl;
	}

	
	
	// ------- prepare to generate some data for the tests -----------

    // set up a random number generator
    const gsl_rng_type * T;
    gsl_rng * r;

    double sigma_x=1;   // distribution parameter for Biv Gaussian
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    const int n=10;    // number to generate 10000
    int n4 = 7; // for histogram 4
       
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

	string theData1Filename = "theData1.txt";
	outputFileStart(theData1Filename);
	outputFile(theData1Filename, theData1Filename, theData1);
    
    int d2 = 2; // dimension of the box to sample data from
    ivector pavingBox2(d2);
    interval pavingInterval2(-5,5);
    for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval2;

        RVecData theData2;   // a container for all the points generated

        // make a simulated data set allData to sample from
        for (int i = 0; i < n; i++) {

            rvector thisrv(d2);
            double x = 0;
            double y = 0;

            gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                    rho, &x, &y);
            thisrv[1] = x;
            thisrv[2] = y;

            // put points generated into container
            theData2.push_back(thisrv);

            }  // data  should be in theData

    int d3 = 2; // dimension of the box to sample data from
    ivector pavingBox3(d3);
    interval pavingInterval3(-5,5);
    for(int k=1; k <= d3; k++) pavingBox3[k] = pavingInterval3;

        RVecData theData3;   // a container for all the points generated

        // make a simulated data set allData to sample from
        for (int i = 0; i < n; i++) {

            rvector thisrv(d3);
            double x = 0;
            double y = 0;

            gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                    rho, &x, &y);
            thisrv[1] = x;
            thisrv[2] = y;

            // put points generated into container
            theData3.push_back(thisrv);

		}  // data  should be in theData

	int d4 = 2; // dimension of the box to sample data from
    ivector pavingBox4(d4);
    interval pavingInterval4(-5,5);
    for(int k=1; k <= d4; k++) pavingBox4[k] = pavingInterval4;

        RVecData theData4;   // a container for all the points generated

		 // make a simulated data set allData to sample from
        for (int i = 0; i < n4; i++) {

            rvector thisrv(d4);
            double x = 0;
            double y = 0;

            gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                    rho, &x, &y);
            thisrv[1] = x;
            thisrv[2] = y;

            // put points generated into container
            theData4.push_back(thisrv);

		}  // data  should be in theData4



    // free the random number generator
    gsl_rng_free (r);

    AdaptiveHistogram myHistFirst(pavingBox1, true); // keep all stats for hist1
    myHistFirst.splitToShape("3,4,4,2,2,4,4,3");

    AdaptiveHistogram myHistSecond(pavingBox2);
    //myHistSecond.splitToShape("1,2,3,3");
    myHistSecond.splitToShape("2,3,4,4,2,3,3");

    AdaptiveHistogram myHistThird(pavingBox3);
    myHistThird.splitToShape("1,2,3,3");
	
	AdaptiveHistogram myHistFourth(pavingBox4);
    //myHistSecond.splitToShape("1,2,3,3");
    myHistFourth.splitToShape("2,3,4,4,2,3,3"); // same splits as for 2


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
		assert(myHistFirst.getRootCounter() == n);
		assert(myHistSecond.getRootCounter() == n);
		assert(myHistThird.getRootCounter() == n);
		assert(myHistFourth.getRootCounter() == n4);
			
		cout << endl << endl;
		cout << "\nStart Test\n:" << endl;
		
		string outputFileName1 = "Hist1SPS.txt";
        myHistFirst.outputToTxtTabs(outputFileName1);
		cout << "Histogram 1 is in " << outputFileName1 << endl;
		
		string outputFileName2 = "Hist2SPS.txt";
        myHistSecond.outputToTxtTabs(outputFileName2);
		cout << "Histogram 2 is in " << outputFileName2 << endl;
		
		string outputFileName3 = "Hist3SPS.txt";
        myHistThird.outputToTxtTabs(outputFileName3);
		cout << "Histogram 3 is in " << outputFileName3 << endl;
		
		string outputFileName4 = "Hist4SPS.txt";
        myHistFourth.outputToTxtTabs(outputFileName4);
		cout << "Histogram 4 is in " << outputFileName4 << endl;

        // 
		
		// try looking directly at the subpavings
		cout << "\nTest subpavings directly\n" << endl;
		SPSnode spsHist1;
		SPSnode spsHist2;
		SPSnode spsHist3;
		SPSnode spsHist4;
		
		/* Note that the dataItrs in these nodes are iterators to the
		 * BigDataCollections in the original hists:  any operations
		 * we do here that actually use the dataItrs (like resetting
		 * countsOnly to false) work only because the hists 
		 * themselves remain in scope
		 * */
		
		try {
			cout << "\nCopy constructor of default constructed node" << endl;
			
			SPSnode cpy(spsHist3);
			
			string s = "copyConstructorDefaultNode.txt";
			outputSPS(s, cpy);
			assert(cpy.isEmpty());
			cout << "Passed assert that node is empty" << endl;;
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction of default constructed SPSnode:\n" << msg << endl;
		}
		
		try {
			cout << "\nCopy constructor of hist1 paving" << endl;
			
			SPSnode cpy1(*(myHistFirst.getSubPaving()));
			assert(!cpy1.isEmpty());
			cout << "Passed assert that node is not empty" << endl;

			string s = "copyConstructorSPSnodeHist1.txt";
			outputSPS(s, cpy1);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction of SPSnode for hist1:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nCopy assignment" << endl;
			SPSnode cpy1(*(myHistFirst.getSubPaving()));
			spsHist1 = cpy1;
			
			string s = "copyAssignmentSPSnodeHist1.txt";
			outputSPS(s, spsHist1);

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy assignment of SPSnode for hist1:\n" << msg << endl;
		}
		
		try {
			
			spsHist2 = SPSnode(*(myHistSecond.getSubPaving()));
			
			string s = "copyAssignmentSPSnodeHist2.txt";
			outputSPS(s, spsHist2);
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy assignment of SPSnode for hist2:\n" << msg << endl;
		}
		
		try {
			
			spsHist4 = SPSnode(*(myHistFourth.getSubPaving()));
			
			string s = "copyAssignmentSPSnodeHist4.txt";
			outputSPS(s, spsHist4);
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy assignment of SPSnode for hist4:\n" << msg << endl;
		}
		
		SPSnode spsHist1c;
		SPSnode spsHist2c;
		
		try {
			
			cout << "\ncopy assignment copies of hist1 and hist2 pavings" << endl;
			
			spsHist1c = spsHist1;
			spsHist2c = spsHist2;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy assignment for hist1 and hist2 pavings:\n" << msg << endl;
		}

	// try to split to shape from vector
	
		try {
			
			cout << "\ntry splitRootAtLeastToShapeSPS, empty reqDepths" << endl;
			
			std::vector < size_t > reqDepths;
			
			SPSnode spn(pavingBox1);
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths);
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to do splitRootAtLeastToShapeSPS, empty reqDepths:\n" << msg << endl;
		}
		
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS, minpoints not doable" << endl;
			
			size_t depths[] = {1, 1};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(pavingBox1);
			
			size_t minPoints = 1;
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths, minPoints);
			
			if (success) throw std::logic_error("Should not be able to do that");
			
		}
		
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS with less split than exists already" << endl;
			
			size_t depths[] = {1,2,3,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(spsHist1);
			
			size_t minPoints = 1;
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths, minPoints);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
		}
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS with exactly same split as exists already" << endl;
			
			size_t depths[] = {3,4,4,2,2,4,4,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(spsHist1);
			
			size_t minPoints = 1;
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths, minPoints);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
		}
		
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS with more splits in some places" << endl;
			
			size_t depths[] = {3,4,4,2,2,4,5,5,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(spsHist1);
			
			//size_t minPoints = 1;
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
		}
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS with more splits in some places but minPoints will not allow the full required split" << endl;
			
			size_t depths[] = {3,5,5,4,2,2,5,5,4,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(spsHist1);
			
			size_t minPoints = 1;
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths, minPoints);
			
			if (success) throw std::logic_error("Should not have been able to do that");
			
			cout << "Success = " << success << endl;
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
			
		}
		{
			
			cout << "\nsplitRootAtLeastToShapeSPS with more splits everywhere" << endl;
			
			size_t depths[] = {4,4,5,5,5,6,6,3,3,3,4,4,5,5,6,6,5,4,5,5};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPSnode spn(spsHist1);
			
			bool success = spn.splitRootAtLeastToShapeSPS(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			cout << "Success = " << success << endl;
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
			
		}
		
		
		
	// tree likelihood
	
		try {
			
			cout << "\nTree likelihood for empty paving (should fail)" << endl;
			
			cxsc::real logl = spsHist3.getUnscaledTreeLogLik();
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (subpavings::NoBox_Error& nbe) {
			std::string msg(nbe.what());
			cout << "\nFailed to do tree likelihood for empty paving:\n" << msg << endl;
		}

		try {
			
			cout << "\nTree likelihood for spsHist1 paving" << endl;
			
			cxsc::real logl = spsHist1.getUnscaledTreeLogLik();
			
			cout << "\nTree likelihood = " << logl << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do tree likelihood for spsHist1:\n" << msg << endl;
			throw;
		}
		
		try {
			
			cout << "\nTree likelihood for spsHist2 paving" << endl;
			
			cxsc::real logl = spsHist2.getUnscaledTreeLogLik();
			
			cout << "\nTree likelihood = " << logl << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do tree likelihood for spsHist2:\n" << msg << endl;
			throw;
		}
		
		try {
			
			cout << "\nTree likelihood for spsHist4 paving" << endl;
			
			cxsc::real logl = spsHist4.getUnscaledTreeLogLik();
			
			cout << "\nTree likelihood = " << logl << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do tree likelihood for spsHist4:\n" << msg << endl;
			throw;
		}
	
	// means and covariances
	
		try {
			cout << "\nMean of default constructed node" << endl;
			
			SPSnode temp;
			
			rvector mean = temp.getMean();
			
			assert(VecLen(mean) == 0);
			cout << "Passed assert that mean is empty" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of default constructed SPSnode:\n" << msg << endl;
		}
		try {
			cout << "\nVarCov of default constructed node" << endl;
			
			SPSnode temp;
			
			RealVec varcov = temp.getVarCovar();
			
			assert(varcov.empty());
			cout << "Passed assert that variance-covariance is empty" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of default constructed SPSnode:\n" << msg << endl;
		}
		
		try {
			cout << "\nMean of node with box but no data" << endl;
			
			SPSnode temp(pavingBox1);
			
			rvector mean = temp.getMean();
			
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of node with box but no data:\n" << msg << endl;
		}
		try {
			cout << "\nVarCov of node with box but no data" << endl;
			
			SPSnode temp(pavingBox1);
			
			RealVec varcov = temp.getVarCovar();
			
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of node with box but no data:\n" << msg << endl;
		}
		try {
			cout << "\nMean of node with box and single data point, keeping all stats" << endl;
			
			AdaptiveHistogram adh(pavingBox1, true ); // keep stats
			RVecData data;
			data = getDataSinglePoint(data);
			int nn = data.size();
			bool inserted = adh.insertFromRVec(data);
			if (inserted) {
				
				SPSnode temp(*adh.getSubPaving());
			
				rvector mean = temp.getMean();
				
				assert(VecLen(mean) == d1);
				assert(!checkAllNaN(mean));
				rvector calcMean = checkMean(data);
				assert(checkSame(mean, calcMean, nn));
				cout << "Passed asserts that mean is length " << d1 << " and calculation of mean matches check calc" << endl;
			}
			else cout << "failed to insert single data point" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of node with box and single data point, keeping all stats:\n" << msg << endl;
		}
		
		try {
			cout << "\nVarCov of node with box and single data point, keeping all stats" << endl;
			
			AdaptiveHistogram adh(pavingBox1, true ); // keep stats
			RVecData data;
			data = getDataSinglePoint(data);
			bool inserted = adh.insertFromRVec(data);
			if (inserted) {
				
				SPSnode temp(*adh.getSubPaving());
			
				RealVec vcov = temp.getVarCovar();
				
				assert(vcov.size() == d1*d1);
				assert(checkAllNaN(vcov));
				cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
				
			}
			else cout << "failed to insert single data point" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of node with box and single data point, keeping all stats:\n" << msg << endl;
		}

		try {
			cout << "\nMean of node with box and two 5-d data points, keeping all stats" << endl;
			
			int d = 5; // dimension of the box to sample data from
			ivector pavingBox5d(d);
			interval pavingInterval5d(-5,5);
			for(int k=1; k <= d; k++) pavingBox5d[k] = pavingInterval5d;

			AdaptiveHistogram adh(pavingBox5d, true ); // keep stats
			RVecData data;
			data = getDataTwoFiveDPoints(data);
			int nn = data.size();
			bool inserted = adh.insertFromRVec(data);
			if (inserted) {
				
				SPSnode temp(*adh.getSubPaving());
			
				rvector mean = temp.getMean();
				
				assert(VecLen(mean) == d);
				assert(!checkAllNaN(mean));
				rvector calcMean = checkMean(data);
				assert(checkSame(mean, calcMean, nn));
				cout << "Passed asserts that mean is length " << d << " and calculation of mean matches check calc" << endl;
			}
			else cout << "failed to insert two 5-d points" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of node with box and two 5-d data points, keeping all stats:\n" << msg << endl;
		}
		
		try {
			cout << "\nVarCov of node with box and two 5-d data points, keeping all stats" << endl;
			
			int d = 5; // dimension of the box to sample data from
			ivector pavingBox5d(d);
			interval pavingInterval5d(-5,5);
			for(int k=1; k <= d; k++) pavingBox5d[k] = pavingInterval5d;

			AdaptiveHistogram adh(pavingBox5d, true ); // keep stats
			RVecData data;
			data = getDataTwoFiveDPoints(data);
			int nn = data.size();
			bool inserted = adh.insertFromRVec(data);
			if (inserted) {
				
				SPSnode temp(*adh.getSubPaving());
			
				RealVec vcov = temp.getVarCovar();
				
				cout << "variance-covariance is " << vcov << endl;
			
				assert(vcov.size() == d*d);
				assert(!checkAllNaN(vcov));
				RealVec calcVCov = checkVarCov(data);
				assert(calcVCov.size() == d*d);
				cout << "variance-covariance from checkVarCov is " << calcVCov << endl;
				assert(checkSame(vcov, calcVCov, nn) );
				cout << "Passed asserts that variance-covariance is length " << d*d << " and calculation of variance-covariance matches check calc" << endl;
				
			}
			else cout << "failed to insert two 5-d points" << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of node with box and two 5-d data points, keeping all stats:\n" << msg << endl;
		}


		try {
			cout << "\nMean of hist 4 node (keep all stats was false)" << endl;
			
			rvector mean = spsHist4.getMean();
				
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of hist 4 node (keep all stats was false):\n" << msg << endl;
		}
		
		try {
			cout << "\nVarCov of hist 4 node (keep all stats was false)" << endl;
			
			RealVec vcov = spsHist4.getVarCovar();
			
			assert(vcov.size() == d1*d1);
			assert(checkAllNaN(vcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
					
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of hist 4 node (keep all stats was false):\n" << msg << endl;
		}
		
		try {
			cout << "\nMean of hist 1 node (keep all stats was true)" << endl;
			
			rvector mean = spsHist1.getMean();
			cout << "mean from spsHist1 is ";
			prettyPrint(cout, mean);
			cout << endl;
				
			assert(VecLen(mean) == d1);
			assert(!checkAllNaN(mean));
			rvector calcMean = checkMean(theData1);
			
			cout << "mean from checkMean is ";
			prettyPrint(cout, calcMean);
			cout << endl;
			
			assert(checkSame(mean, calcMean, n));
			cout << "Passed asserts that mean is length " << d1 << " and calculation of mean matches check calc" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do mean of hist 1 node (keep all stats was true):\n" << msg << endl;
		}
		
		try {
			cout << "\nVarCov of hist 1 node (keep all stats was true)" << endl;
			
			RealVec vcov = spsHist1.getVarCovar();
			
			cout << "variance-covariance from spsHist1 is " << vcov << endl;
			
			assert(vcov.size() == d1*d1);
			assert(!checkAllNaN(vcov));
			RealVec calcVCov = checkVarCov(theData1);
			assert(calcVCov.size() == d1*d1);
			cout << "variance-covariance from checkVarCov is " << calcVCov << endl;
			assert(checkSame(vcov, calcVCov, n) );
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " calculation of variance-covariance matches check calc" << endl;

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do varcov of hist 1 node (keep all stats was true):\n" << msg << endl;
		}
		
		// clear all data and resetting counts only
		
		try {
			cout << "\nClear a copy of spsHist1 and check counts etc" << endl;
			
			SPSnode temp(spsHist1);
			
			temp.clearAllDataHeld();
			
			assert(temp.getCountsOnly() == spsHist1.getCountsOnly());
			cout << "Passed assert that countsOnly is as before" << endl;
			assert(temp.getCounter() == 0);
			cout << "Passed assert that counter is 0" << endl;
			rvector mean = temp.getMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
			RealVec varcov = temp.getVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			string s = "clearedCopySPSnodeHist1.txt";
			outputSPS(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do clear a copy of spsHist1:\n" << msg << endl;
		}
		try {
			cout << "\nClear a copy of spsHist4 and check counts etc" << endl;
			
			SPSnode temp(spsHist4);
			
			temp.clearAllDataHeld();
			
			assert(temp.getCountsOnly() == spsHist4.getCountsOnly());
			cout << "Passed assert that countsOnly is as before" << endl;
			assert(temp.getCounter() == 0);
			cout << "Passed assert that counter is 0" << endl;
			rvector mean = temp.getMean();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
			RealVec varcov = temp.getVarCovar();
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			string s = "clearedCopySPSnodeHist4.txt";
			outputSPS(s, temp);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do clear a copy of spsHist4:\n" << msg << endl;
		}

		// when testing with histograms, reinsert same data and check means etc match, then 
		// reinsert different data and check means etc match a different hist
		
		try {
			cout << "\nReset counts only on a copy of spsHist1 and check counts etc" << endl;
			
			SPSnode temp(spsHist1);
			
			cout << "Reset counts only to true" << endl;
				{bool newCountsOnly = true;
				temp.setCountsOnly(newCountsOnly);
				
				assert(temp.getCountsOnly() == newCountsOnly);
				cout << "Passed assert that countsOnly is " << newCountsOnly << endl;
				assert(temp.getCounter() == spsHist1.getCounter());
				cout << "Passed assert that counter is as before" << endl;
				rvector mean = temp.getMean();
				assert(VecLen(mean) == d1);
				assert(checkAllNaN(mean));
				cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
				RealVec varcov = temp.getVarCovar();
				assert(varcov.size() == d1*d1);
				assert(checkAllNaN(varcov));
				cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
				
				string s = "countsOnlyResetToTrueCopySPSnodeHist1.txt";
				outputSPS(s, temp);
			}
			
			cout << "\nReset counts only to false" << endl;
			{
				bool newCountsOnly = false;
				temp.setCountsOnly(newCountsOnly);
				
				assert(temp.getCountsOnly() == newCountsOnly);
				cout << "Passed assert that countsOnly is " << newCountsOnly << endl;
				assert(temp.getCounter() == spsHist1.getCounter());
				cout << "Passed assert that counter is as before" << endl;
				rvector mean = temp.getMean();
				cout << " mean is ";
				prettyPrint(cout, mean);
				cout << endl;
				rvector mean1 = spsHist1.getMean();
				cout << " mean for spsHist1 is ";
				prettyPrint(cout, mean1);
				cout << endl;
				assert(VecLen(mean) == d1);
				assert( checkSame( mean, spsHist1.getMean(), temp.getCounter() ) );
				cout << "Passed asserts that mean is length " << d1 << " and all elements are as in spsHist1" << endl;
				RealVec varcov = temp.getVarCovar();
				cout << "variance-covariance is " << varcov << endl;
				RealVec varcov1 = spsHist1.getVarCovar();
				cout << "variance-covariance for spsHist1 is " << varcov1 << endl;
				assert(varcov.size() == d1*d1);
				assert( checkSame( varcov, varcov1, temp.getCounter() ) );
				cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are as in spsHist1" << endl;
				
				string s = "countsOnlyResetToFalseCopySPSnodeHist1.txt";
				outputSPS(s, temp);
			}
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reset counts only on a copy of spsHist1:\n" << msg << endl;
		}
		
		try {
			cout << "\nReset counts only on a copy of spsHist4 and check counts etc" << endl;
			
			SPSnode temp(spsHist4);
			
			cout << "Reset counts only to false" << endl;
				{bool newCountsOnly = false;
				temp.setCountsOnly(newCountsOnly);
				
				assert(temp.getCountsOnly() == newCountsOnly);
				cout << "Passed assert that countsOnly is " << newCountsOnly << endl;
				assert(temp.getCounter() == spsHist4.getCounter());
				cout << "Passed assert that counter is as before" << endl;
				cxsc::rvector mean = temp.getMean();
				cout << " mean is ";
				prettyPrint(cout, mean);
				cout << endl;
				assert(VecLen(mean) == d1);
				assert(!checkAllNaN(mean));
				rvector calcMean = checkMean(theData4);
				cout << " mean from checkMean is ";
				prettyPrint(cout, calcMean);
				cout << endl;
				assert( checkSame( mean, calcMean, temp.getCounter() ) );
			
				cout << "Passed asserts that mean is length " << d1 << " and is same as check calculation" << endl;
				RealVec varcov = temp.getVarCovar();
				assert(varcov.size() == d1*d1);
				assert(!checkAllNaN(varcov));
				RealVec calcVCov = checkVarCov(theData4);
				cout << "variance-covariance is " << varcov << endl;
				cout << "variance-covariance for from checkVarCov is " << calcVCov << endl;
				assert(checkSame(varcov, calcVCov, temp.getCounter()) );
				cout << "Passed asserts that variance-covariance is length " << d1*d1 << " calculation of variance-covariance matches check calc" << endl;
				
				string s = "countsOnlyResetToFalseCopySPSnodeHist4.txt";
				outputSPS(s, temp);
			}
			
			cout << "\nReset counts only to true" << endl;
			{
				bool newCountsOnly = true;
				temp.setCountsOnly(newCountsOnly);
				
				assert(temp.getCountsOnly() == newCountsOnly);
				cout << "Passed assert that countsOnly is " << newCountsOnly << endl;
				assert(temp.getCounter() == spsHist4.getCounter());
				cout << "Passed assert that counter is as before" << endl;
				rvector mean = temp.getMean();
				assert(VecLen(mean) == d1);
				assert( checkAllNaN( mean ) );
				cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
				RealVec varcov = temp.getVarCovar();
				RealVec varcov1 = spsHist1.getVarCovar();
				assert(varcov.size() == d1*d1);
				assert( checkAllNaN( varcov ) );
				cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
				
				string s = "countsOnlyResetToTrueCopySPSnodeHist4.txt";
				outputSPS(s, temp);
			}
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do reset counts only on a copy of spsHist4:\n" << msg << endl;
		}
		
		
		// swap
		
		try {
			
			cout << "\nswap copies of copies of hist1 and hist2 pavings" << endl;
			
			string sc11 = "swapCheckForHist1PavingSwap1.txt";
			string sc21 = "swapCheckForHist2PavingSwap1.txt";
			
			swapCheckOutput(sc11, spsHist1c);
			swapCheckOutput(sc21, spsHist2c);
			
			std::swap(spsHist1c, spsHist2c);
			
			string s1 = "shouldNowBeCopyOfHist1PavingSwap1.txt";
			outputSPS(s1, spsHist2c);
			
			string s2 = "shouldNowBeCopyOfHist2PavingSwap1.txt";
			outputSPS(s2, spsHist1c);
			
			string sc12 = "swapCheckForShouldNowBeCopyOfHist1PavingSwap1.txt";
			string sc22 = "swapCheckForShouldNowBeCopyOfHist2PavingSwap1.txt";
			
			swapCheckOutput(sc12, spsHist2c);
			swapCheckOutput(sc22, spsHist1c);
			
			cout << "swapCheckOutput for copy of hist1 paving before swap is in " << sc11 << " and should-be-copy-of-hist-1-paving after swap in " << sc12 << endl;
			cout << "swapCheckOutput for copy of hist2 paving before swap is in " << sc21 << " and should-be-copy-of-hist-2-paving after swap in " << sc22 << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do swap of copies of copies of hist1 and hist2 pavings:\n" << msg << endl;
		}
		try {
			
			cout << "\nswap more copies of hist4 and hist2 pavings" << endl;
			
			string sc11 = "swapCheckForHist4PavingSwap2.txt";
			string sc21 = "swapCheckForHist2PavingSwap2.txt";
			SPSnode temp1(spsHist4);
			SPSnode temp2(spsHist2);
			swapCheckOutput(sc11, temp1);
			swapCheckOutput(sc21, temp2);
			
			std::swap(temp1, temp2);
			
			string s1 = "shouldNowBeCopyOfHist4PavingSwap2.txt";
			outputSPS(s1, temp2);
			
			string s2 = "shouldNowBeCopyOfHist2PavingSwap2.txt";
			outputSPS(s2, temp1);
			
			string sc12 = "swapCheckForShouldNowBeCopyOfHist4PavingSwap2.txt";
			string sc22 = "swapCheckForShouldNowBeCopyOfHist2PavingSwap2.txt";
			
			swapCheckOutput(sc12, temp2);
			swapCheckOutput(sc22, temp1);
			
			cout << "swapCheckOutput for copy of hist4 paving before swap is in " << sc11 << " and should-be-copy-of-hist-4-paving after swap in " << sc12 << endl;
			cout << "swapCheckOutput for copy of hist2 paving before swap is in " << sc21 << " and should-be-copy-of-hist-2-paving after swap in " << sc22 << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do swap of more copies of hist4 and hist2 pavings:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure with empty paving spsHist3 (this should fail)" << endl;
			
			spsHist3.unionTreeStructure(&spsHist1);
			
			throw std::logic_error("Should not be able to do this");

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure with empty paving spsHist3:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure copy of hist1 paving with empty paving" << endl;
			
			SPSnode temp1(spsHist1);
			SPSnode temp2;
			
			temp1.unionTreeStructure(&temp2);
			
			assert(temp1.getRootCounter() == 0 );
			cout << "Passed assert that new node data counter == 0" << endl;
			
			cxsc::rvector mean = temp1.getMean();
			RealVec varcov = temp1.getVarCovar();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
		
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			string s = "unionOfHist1PavingAndEmptyPaving.txt";
			outputSPS(s, temp1);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure copy of hist1 paving with empty paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure cpy of empty hist3 paving with empty paving (this should fail)" << endl;
			
			SPSnode temp1(spsHist3);
			SPSnode temp2;
			
			temp1.unionTreeStructure(&temp2);
			
			throw std::logic_error("Should not be able to do this");

			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure empty hist3 paving with empty paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure copy of hist2 paving with null paving ptr" << endl;
			
			SPSnode* tempPtr = NULL;
			SPSnode temp1(spsHist2);
			
			temp1.unionTreeStructure(tempPtr);

			assert(temp1.getRootCounter() == 0 );
			cout << "Passed assert that new node data counter == 0" << endl;
			cxsc::rvector mean = temp1.getMean();
			RealVec varcov = temp1.getVarCovar();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
		
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			string s = "unionOfHist2PavingAndNullPavingPtr.txt";
			outputSPS(s, temp1);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure copy of hist2 paving with null paving ptr:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure copy of hist 1 and hist 2 pavings" << endl;
			
			SPSnode temp1(spsHist1);
			const SPSnode temp2(spsHist2);
			
			temp1.unionTreeStructure(&temp2);
			assert(temp1.getRootCounter() == 0 );
			cout << "Passed assert that new node data counter == 0" << endl;
			cxsc::rvector mean = temp1.getMean();
			RealVec varcov = temp1.getVarCovar();
			assert(VecLen(mean) == d1);
			assert(checkAllNaN(mean));
			cout << "Passed asserts that mean is length " << d1 << " and all elements are NaN" << endl;
		
			assert(varcov.size() == d1*d1);
			assert(checkAllNaN(varcov));
			cout << "Passed asserts that variance-covariance is length " << d1*d1 << " and all elements are NaN" << endl;
			
			string s = "unionOfHist1AndHist2Pavings.txt";
			outputSPS(s, temp1);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure copy of hist1 and hist2 pavings:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nunion tree structure hist 1 copy swapped with hist 2 copy pavings" << endl;
			
			SPSnode temp1(spsHist1c);
			const SPSnode temp2(spsHist2c);
			
			temp1.unionTreeStructure(&temp2);
			
			string s = "unionOfHist1CopySwappedAndHist2CopySwappedPavings.txt";
			outputSPS(s, temp1);
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do union tree structure hist1 and hist2 pavings:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance empty paving with hist1 paving (this should fail)" << endl;
			
			SPSnode temp;
			
			cxsc::real dis = temp.getL1Distance(&spsHist1);
			
			throw std::logic_error("Should not be able to do this");


			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance empty paving with hist1 paving:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance hist1 paving with null pointer paving (this should fail)" << endl;
			
			const SPSnode temp(spsHist1);
			SPSnode* tempPtr = NULL;
			
			cxsc::real dis = temp.getL1Distance(tempPtr);
			
			throw std::logic_error("Should not be able to do this");

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with null pointer paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance hist1 paving with empty paving (this should fail)" << endl;
			
			const SPSnode temp1(spsHist1);
			SPSnode temp2;
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with empty paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance empty paving with empty paving (this should fail)" << endl;
			
			const SPSnode temp1;
			SPSnode temp2;
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
		throw std::logic_error("Should not be able to do this");

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance empty paving with empty paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance hist1 paving with paving with wrong box dimensions (this should fail)" << endl;
			
			const SPSnode temp1(spsHist1);
			
			int dT1 = 3; // dimension of the box to sample data from
			ivector pavingBoxT1(dT1);
			//interval pavingIntervalT1(-5,5);
			for(int k=1; k <= dT1; k++) pavingBoxT1[k] = pavingInterval1;
		
			SPSnode temp2(pavingBoxT1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);

			throw std::logic_error("Should not be able to do this");

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with paving with wrong dimensions:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance hist1 paving with paving with wrong box lengths (this should fail)" << endl;
			
			const SPSnode temp1(spsHist1);
			
			ivector pavingBoxT1(d1);
			
			for(int k=1; k < d1; k++) pavingBoxT1[k] = pavingInterval1;
			
			interval pavingIntervalT1(-5,4);
			pavingBoxT1[d1] = pavingIntervalT1;
			
			SPSnode temp2(pavingBoxT1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with paving with wrong box lengths:\n" << msg << endl;
		}
		
		
		
		std::string L1Filename("unionOfEmptyHist3PavingAndEmptyPaving.txt");
		try {
			
			cout << "\nL1 distance hist1 paving with itself (should be 0)" << endl;
			
			const SPSnode temp1(spsHist1);
			//const SPSnode temp2(pavingBox1);
			
			cxsc::real dis = temp1.getL1Distance(&temp1);
			
			assert(dis == 0.0);
			cout << "Passed assert that this distance is 0" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with itself:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance hist1 paving with copy of itself (should be 0)" << endl;
			
			const SPSnode temp1(spsHist1);
			const SPSnode temp2(temp1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			assert(dis == 0.0);
			cout << "Passed assert that this distance is 0" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with copy of itself:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance paving with box but no data with another paving with right box, no data (should be 0)" << endl;
			
			const SPSnode temp1(pavingBox1);
			const SPSnode temp2(pavingBox1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			assert(dis == 0.0);
			cout << "Passed assert that this distance is 0" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance paving with box but no data with another paving with right box, no data:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance hist1 paving with paving with right box, no data (should be 1)" << endl;
			
			const SPSnode temp1(spsHist1);
			const SPSnode temp2(pavingBox1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			assert(dis == 1.0);
			cout << "Passed assert that this distance is 1" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance hist1 paving with paving with right box, no data:\n" << msg << endl;
		}
		try {
			
			cout << "\nL1 distance paving with right box, no data, with hist1 paving (should be 1) " << endl;
			
			const SPSnode temp1(spsHist1);
			const SPSnode temp2(pavingBox1);
			
			cxsc::real dis = temp2.getL1Distance(&temp1);
			
			assert(dis == 1.0);
			cout << "Passed assert that this distance is 1" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance paving with right box, no data, with hist1 paving:\n" << msg << endl;
		}
		
		cxsc::real dis1_2 = 0.0;
		
		try {
			
			cout << "\nL1 distance copy hist1 paving with copy hist2 paving " << endl;
			
			const SPSnode temp1(spsHist1);
			const SPSnode temp2(spsHist2);
			
			dis1_2 = temp1.getL1Distance(&temp2);
			
			cxsc::real shouldBe(1.5);
			bool isNotEqual = cxsc::abs(dis1_2 - shouldBe)
				> DBL_EPSILON * cxsc::max(dis1_2, shouldBe);
			assert(!isNotEqual);
			
			cout << "Passed assert that this distance is " << shouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance copy hist1 paving with copy hist2 paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance copy hist2 paving with copy hist1 paving " << endl;
			
			const SPSnode temp1(spsHist2);
			const SPSnode temp2(spsHist1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			bool isNotEqual = cxsc::abs(dis - dis1_2) 
				> DBL_EPSILON * cxsc::max(dis, dis1_2);
			assert(!isNotEqual);
			cout << "Passed assert that this distance is the same as between 1 and 2 " << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance copy hist2 paving with copy hist1 paving:\n" << msg << endl;
		}
		
		cxsc::real dis1_4 = 0.0;
		try {
			
			cout << "\nL1 distance copy hist1 paving with copy hist4 paving " << endl;
			
			const SPSnode temp1(spsHist1);
			const SPSnode temp2(spsHist4);
			
			dis1_4 = temp1.getL1Distance(&temp2);
			
			cxsc::real shouldBe(101.0/70.0);
			bool isNotEqual = cxsc::abs(dis1_4 - shouldBe) 
				> DBL_EPSILON * cxsc::max(dis1_4, shouldBe);
			assert(!isNotEqual);
			
			cout << "Passed assert that this distance is " << shouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance copy hist1 paving with copy hist4 paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance copy hist4 paving with copy hist1 paving " << endl;
			
			const SPSnode temp1(spsHist4);
			const SPSnode temp2(spsHist1);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			bool isNotEqual = cxsc::abs(dis - dis1_4) 
				> DBL_EPSILON * cxsc::max(dis, dis1_4);
			assert(!isNotEqual);
			
			cout << "Passed assert that this distance is the same as between 1 and 4 " << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance copy hist4 paving with copy hist1 paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nL1 distance copy hist2 paving with copy hist4 paving " << endl;
			
			const SPSnode temp1(spsHist2);
			const SPSnode temp2(spsHist4);
			
			cxsc::real dis = temp1.getL1Distance(&temp2);
			
			cxsc::real shouldBe(38.0/70.0);
			
			bool isNotEqual = cxsc::abs(dis - shouldBe) 
				> DBL_EPSILON * cxsc::max(dis, shouldBe);
			assert(!isNotEqual);
			
			cout << "Passed assert that this distance is " << shouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do L1 distance copy hist2 paving with copy hist4 paving:\n" << msg << endl;
		}
		try {
			
			cout << "\nFind containing node with empty paving (should fail)" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = 0;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const SPSnode temp;
			
			const SPSnode* node = temp.findContainingNode(pt);
			
			throw std::logic_error("Should not be able to do this");

		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do find containing node with empty paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nFind containing node with hist1 paving (should not find containing node)" << endl;
			
			cxsc::rvector pt(d1);
			pt[1] = -6;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const SPSnode* node = spsHist1.findContainingNode(pt);
			
			assert(node == NULL);
			cout << "Passed assert that no containing node found" << endl;
			
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do find containing node with hist1 paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nFind containing node with hist1 paving (should find all)" << endl;
			
			{
				cxsc::rvector pt(d1);
				pt[1] = -5;
				pt[2] = -5;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XLLL");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = -5;
				pt[2] = 5;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XLR");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = 5;
				pt[2] = -5;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XRL");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = 5;
				pt[2] = 5;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XRRR");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do find containing node with hist1 paving:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nMore find containing node with hist1 paving (should find)" << endl;
			
			{
				cxsc::rvector pt(d1);
				pt[1] = 0;
				pt[2] = 0;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XRRLL");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = 2.5;
				pt[2] = 0;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XRRR");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = -5;
				pt[2] = 0;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XLR");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			{
				cxsc::rvector pt(d1);
				pt[1] = -2.5;
				pt[2] = -2.5;
				
				cout << "Point ";
				prettyPrint(cout, pt);
				cout << endl;
				
				const SPSnode* node = spsHist1.findContainingNode(pt);
				assert(node != NULL);
				//cout << "Containing node is " << node->getNodeName() << endl;
				std::string expected("XLLRR");
				assert(node->getNodeName() == expected);
				cout << "Passed assert that containing node is " << expected << endl;
			}
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do find containing node with hist1 paving:\n" << msg << endl;
		}
		
		
		
		
	}
    else cout << "unsuccessful insertion" << endl;

	cout << "\nEnd test\n" << endl;

    return 0;

} // end of test program




std::ostream& swapCheck(std::ostream& os, const SPSnode& spn, int level)
{
	// do me
	for (int i = 0; i < level; ++i) {
		os << "\t";
	}
	
	os << spn.nodeStringSummary() << endl;
	/*
		 
	os << "I am " << spn.getNodeName() << ", my parent is ";
	if (spn.getParent() != NULL) os << spn.getParent()->getNodeName();
	else os << "NULL"; 
	
	os << ", my left child is ";
	if (spn.getLeftChild() != NULL) os << spn.getLeftChild()->getNodeName();
	else os << "NULL"; 
	
	os << ", my right child is ";
	if (spn.getRightChild() != NULL) os << spn.getRightChild()->getNodeName();
	else os << "NULL"; 
	
	os << endl;
	*/
	// recurse
	
	if ( spn.getLeftChild() ) swapCheck(os, *spn.getLeftChild(), level+1);
	if ( spn.getRightChild() ) swapCheck(os, *spn.getRightChild(), level+1);
	
	return os;

}

RVecData& getDataSinglePoint(RVecData& data) 
{
    int d = 2;
	{ 
		rvector thisrv(d);
		thisrv[1] = 1.0;
        thisrv[2] = 1.0;
		data.push_back(thisrv);
	}
	return data;
}

RVecData& getDataTwoFiveDPoints(RVecData& data) 
{
    int d = 5;
	{ 
		rvector thisrv(d);
		thisrv[1] = 1.0;
        thisrv[2] = 1.1;
		thisrv[3] = 1.2;
		thisrv[4] = 1.3;
		thisrv[5] = 1.4;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = 0.9;
        thisrv[2] = -0.5;
		thisrv[3] = 1.0;
		thisrv[4] = 1.1;
		thisrv[5] = 0.5;
		data.push_back(thisrv);
	}
	return data;
}

RVecData& getData2d(RVecData& data) 
{
    int d = 2;
	{ 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -1;
        thisrv[2] = 1;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 3;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -4;
        thisrv[2] = 4;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	
	return data;
}


void checkNode(const string& s, const SPSnode * const spn)
{
	ofstream os;
	if (spn->getParent() == NULL) os.open(s.c_str()); // don't append
	else os.open(s.c_str(), ios_base::app); // append
	if (os.is_open()) {
	
		os <<  endl;
		cxsc::real vol = spn->nodeRealVolume();
		size_t count = spn->getCounter();
		size_t rootCounter = spn->getRootCounter();
		int nodedepth = spn->getNodeDepth();
		cxsc::real sumLeafCountOverVol = spn->getSumLeafCountOverVol();
		size_t getSmallestLeafCount = spn->getSmallestLeafCount();
		size_t getLargestLeafCount = spn->getLargestLeafCount();
		//double smallestVol = spn->getSmallestLeafVol();
		//double largestVol = spn->getLargestLeafVol();
		cxsc::real getLogLik = spn->getLogLik(rootCounter);
		cxsc::dotprecision getSplitChangeLogLik = spn->getSplitChangeLogLik();
		cxsc::dotprecision getEMPSumCOPERR = spn->getEMPSumCOPERR(rootCounter);
		cxsc::dotprecision getEMPSumAIC = spn->getEMPSumAIC(rootCounter);
		cxsc::real getEMPContributionCOPERR = spn->getEMPContributionCOPERR(rootCounter);
		cxsc::real getEMPContributionAIC = spn->getEMPContributionAIC(rootCounter);
		
		cxsc::dotprecision getSplitChangeEMPCOPERR = spn->getSplitChangeEMPCOPERR(rootCounter);
		cxsc::dotprecision getSplitChangeEMPAIC = spn->getSplitChangeEMPAIC();
		cxsc::dotprecision getBestSplitChangeEMPCOPERR = spn->getBestSplitChangeEMPCOPERR(rootCounter);
		cxsc::dotprecision getBestSplitChangeEMPAIC = spn->getBestSplitChangeEMPAIC();
		
		
		os << spn->getNodeName() << " node depth is " << nodedepth << " and isLeaf " << (spn->isLeaf()? "Yes" : "No") << endl;
		os <<  "node count is " << count << " and node volume " << vol << endl;
		os <<  "sumLeafCountOverVol is " << sumLeafCountOverVol << endl;
		os <<  "getLargestLeafCount is " << getLargestLeafCount << " and getSmallestLeafCount is " << getSmallestLeafCount << endl;
		os <<  "getLogLik is " << getLogLik << endl;
		os <<  "getSplitChangeLogLik is " << cxsc::rnd(getSplitChangeLogLik) << endl;
		try {
			cxsc::dotprecision getMergeChangeLogLik = spn->getMergeChangeLogLik();
			os <<  "getMergeChangeLogLik is " << cxsc::rnd(getMergeChangeLogLik) << endl;
			
		}
		catch (subpavings::UnfulfillableRequest_Error const& ure) {
			os <<  "cannot do getMergeChangeLogLik since node is a leaf" << endl;
		}
		os <<  "getEMPSumCOPERR is " << cxsc::rnd(getEMPSumCOPERR) << " and getEMPSumAIC is " << cxsc::rnd(getEMPSumAIC) << endl;
		os <<  "getEMPContributionCOPERR is " << getEMPContributionCOPERR << " and getEMPContributionAIC is " << getEMPContributionAIC << endl;
		
		os <<  "getSplitChangeEMPCOPERR is " << cxsc::rnd(getSplitChangeEMPCOPERR) << " and getSplitChangeEMPAIC is " << cxsc::rnd(getSplitChangeEMPAIC) << endl;
		try {
			cxsc::dotprecision getMergeChangeEMPCOPERR = spn->getMergeChangeEMPCOPERR(rootCounter);
			os <<  "getMergeChangeEMPCOPERR is " << cxsc::rnd(getMergeChangeEMPCOPERR);
			
		}
		catch (subpavings::UnfulfillableRequest_Error const& ure) {
			os <<  "cannot do getMergeChangeEMPCOPERR since node is a leaf";
		}
		try {
			cxsc::dotprecision getMergeChangeEMPAIC = spn->getMergeChangeEMPAIC();
			os <<  " and getMergeChangeEMPAIC is " << cxsc::rnd(getMergeChangeEMPAIC) << endl;
			
		}
		catch (subpavings::UnfulfillableRequest_Error const& ure) {
			os <<  " and cannot do getMergeChangeEMPAIC since node is a leaf" << endl;
		}
		
		os <<  "getBestSplitChangeEMPCOPERR is " << cxsc::rnd(getBestSplitChangeEMPCOPERR) << " and getBestSplitChangeEMPAIC is " << cxsc::rnd(getBestSplitChangeEMPAIC) << endl;
		try {
			cxsc::dotprecision getBestMergeChangeEMPCOPERR = spn->getBestMergeChangeEMPCOPERR(rootCounter);
			os <<  "getBestMergeChangeEMPCOPERR is " << cxsc::rnd(getBestMergeChangeEMPCOPERR);
			
		}
		catch (subpavings::UnfulfillableRequest_Error const& ure) {
			os <<  "cannot do getBestMergeChangeEMPCOPERR since node is a leaf";
		}
		try {
			cxsc::dotprecision getBestMergeChangeEMPAIC = spn->getBestMergeChangeEMPAIC();
		os <<  " and getBestMergeChangeEMPAIC is " << cxsc::rnd(getBestMergeChangeEMPAIC) << endl;
		
		}
		catch (subpavings::UnfulfillableRequest_Error const& ure) {
			os <<  " and cannot do getBestMergeChangeEMPAIC since node is a leaf" << endl;
		}
		
		os.close();
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
		
	// recurse
	if (spn->hasLCwithBox() ) checkNode(s, spn->getLeftChild());
	if (spn->hasRCwithBox() ) checkNode(s, spn->getRightChild());
}

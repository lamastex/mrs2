/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
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
* MERCHANTABILITY or FsITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*! \file
\brief Get the Hellinger distance between actual Gaussian data and Gaussian data 
generated from a function estimate
 */

#include "histall.hpp"  // headers for the histograms
#include "intervalmappedspnode_measurers.hpp" // ordering for pq split
#include "functionestimator_interval.hpp"
#include "piecewise_constant_function.hpp"  

#include "GaussianFobj.hpp" // fobj
#include "toolz.hpp"

#include <vector>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iostream>

#include <limits> // to use negative infinity

// to use assert
#include "assert.h"

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main(int argc, char* argv[])
{
	// sort out user-defined parameters------------------//
	if ( argc < 7 ) {
		cerr << "Syntax: " << argv[0] << 
		" d maxLeavesEst n dataSeed critLeaves maxCheck" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "d n states, symmetryIndicator");
	}
	
	int d = atoi(argv[1]);  // ds
	size_t maxLeavesEst = atoi(argv[2]);  // number of leaves in estimator
	const int n = atoi(argv[3]);  // number of points to generate
	int dataSeed = atoi(argv[4]); // seed for data generation
    size_t critLeaves = atoi(argv[5]); //maximum number of leaves for PQ to stop splitting 
	int maxCheck = atoi(argv[6]); //maximum number of checks based on delta
	size_t startLeaves = atoi(argv[7]);
    
	cout << argv[0] << " : process id is " << getpid() << std::endl;

	// for output purposes
	// string formatting
	ofstream oss;         // ofstream object
   oss << scientific;  // set formatting for input to oss
   oss.precision(10);
   ostringstream stm;

	//=======generate actual data and get the root box==============//
	// set up a random number generator and use mt19937 for generator
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
	//long unsigned int seed = 1234;
	gsl_rng_set (r, dataSeed); // change the seed
	cout << "Data seed is " << dataSeed << endl;

	//cout << "\n Generate " << n << " data from actual density." << endl;

//	RVecData* actualDataPtr = new RVecData;

/*	for (size_t i = 0; i < n; i++) {
		rvector thisrv(d);
		for (size_t j = 1; j <= d; j++) {
			//double z = gsl_rng_uniform(r);
			double z = gsl_ran_gaussian(r, 1.0); // generate a normal r.v.
			thisrv[j] = (z);
		}
		// put points generated into container
		actualDataPtr->push_back(thisrv);
	}

	AdaptiveHistogram* actualHist = new AdaptiveHistogram(true); 
	actualHist->insertFromRVec(*actualDataPtr);
	ivector pavingBoxEst = actualHist->getRootBox();
    
   	delete actualHist;
	delete actualDataPtr;
	*/
	//=========end of getting a root box from the actual data=============
	
	//============== make estimate ============//
	cout << "\nMake the function estimator to " << maxLeavesEst << " leaves" << endl;
	
	// specify function object (from /examples/MappedTargets
	GaussianFobj fobj;
	
	/* function estimate is going to use same box as the histograms */
	// Use fobj and pavingBox to get an estimator
	
	//data generating partition
   ivector pavingBoxEst(d);
   interval pavingInterval(-5,5);
   for(int i=1; i <= d; i++) { pavingBoxEst[i] = pavingInterval; }
	
	FunctionEstimatorInterval estimator(pavingBoxEst, fobj);
	//cout << estimator << endl;
	
	LOGGING_LEVEL logEst = NOLOG; // logging for making estimator
	
	#if(1)
	size_t maxLeavesEstDown = static_cast<size_t>(1.2*maxLeavesEst); 
	// go down to 1.2 x max
	#endif
	#if(0)
		size_t maxLeavesEstDown = maxLeavesEst;
	#endif
	
	cout << "pq down to max leaves " << maxLeavesEstDown << endl;
	
	clock_t startEst = clock();
	
	// priority split driven by splitting leaf with max reimann diff
	ReimannDiffMeasurer measurer;
	estimator.prioritySplit(measurer, maxLeavesEstDown, logEst);
			
	// stop recording time here
	clock_t endEst = clock();
	cout << "Number of leaves in estimate: " << estimator.getRootLeaves() << " s."<< endl;	
	cout << "After split, getTotalAreaOfIntervalBand() = "
		<< estimator.getTotalAreaOfIntervalBand() << endl;
	double timingEst1 = ((static_cast<double>(endEst - startEst)) / CLOCKS_PER_SEC);
	cout << "Computing time for pq split in estimate: " << timingEst1 << " s."<< endl;
	
	startEst = clock();
	#if(1) 
		cout << "Hull propagation" << endl;
		estimator.hullPropagation();
		
		cout << "Priority merge to " << maxLeavesEst << " leaves" << endl;
		#if(0)
		// priority merge driven by minimising increase the reimann diff
		estimator.priorityMergeOnLoss(maxLeavesEst, logEst);
		#endif
		#if(1)
		// priority merge driven by merging cherry with minimum reimann diff
		estimator.priorityMerge(maxLeavesEst, logEst);
		#endif
					
		// stop recording time here
		endEst = clock();	
		double timingEst2 = ((static_cast<double>(endEst - startEst)) / CLOCKS_PER_SEC);
		cout << "Computing time for hull propagate and merge up in estimate: " << timingEst2 << " s."<< endl;
		
		cout << "After propagation and priority merge, getTotalAreaOfIntervalBand() = " 
					<< estimator.getTotalAreaOfIntervalBand() << endl;
		cout << "number of leaves is = " << estimator.getRootLeaves() << endl;
	#endif
	
	cout << "Making estimate and normalising" << endl;
	// Make PiecewiseConstantFunction estimate from estimator
	PiecewiseConstantFunction estimate = estimator.makePiecewiseConstantFunction();
	cout << "estimate has integral " << estimate.getTotalIntegral() << " before normalizing" << endl;
	real before = estimate.getTotalIntegral();
	estimate.normalise();
	cout << "estimate has integral " << estimate.getTotalIntegral() << endl;
	
	//optional
	/*estimate.outputToTxtTabs("PCF.txt");
	string Integral = "Integral.txt";
	oss.open(Integral.c_str());
	oss << before << "\t" << estimate.getTotalIntegral() << endl;
	oss << flush;
	oss.close();
*/
	//===========end of estimating function using PCF=========================//

	//===========generate data==============================================//
	// Use PiecewiseConstantFunction to generate data, supplying our own rng
	cout << "\nGenerating data for simulation" << endl;

	RVecData* theDataPtr = new RVecData;   // a container for all the points generated

	clock_t startData = clock();

	// Gaussian data
	int N = 1.5*n;
	estimate.simulateData(*theDataPtr, N, r);

	// stop recording time here
	clock_t endData = clock();	
	double timingData = ((static_cast<double>(endData - startData)) / CLOCKS_PER_SEC);
	cout << "Computing time for simulating data: " << timingData << " s."<< endl;

	cout << (*theDataPtr).size() << " points generated" << endl;
	
	try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow
	
	/* //optional 
	string dataFileName = "ActualData";
	dataFileName += stm.str(); 
	dataFileName += ".txt"; 
	oss.open(dataFileName.c_str());
	for (size_t i = 0; i < n; i++) { 
		for (size_t j = 1; j <= d; j++) {
				oss << (*actualDataPtr)[i][j] << "\t";
		}
		oss << "\n";
		//cout << "\n";
	}
	oss << flush;
	oss.close();

	cout << "Actual data written to  " << dataFileName << endl;
	
	//optional 
	dataFileName = "MappedData";
	dataFileName += stm.str(); 
	dataFileName += ".txt"; 
	oss.open(dataFileName.c_str());
	for (size_t i = 0; i < n; i++) { 
		for (size_t j = 1; j <= d; j++) {
				oss << (*theDataPtr)[i][j] << "\t";
		}
		oss << "\n";
		//cout << "\n";
	}
	oss << flush;
	oss.close();
	cout << "Estimated data written to  " << dataFileName << endl;
    */
	//===========end of generating data=================================//

    //=====Start MDE============================================//
   	cout << "========================================================" << endl;
	cout << "Run hold out estimation..." << endl;
    cout << "\nStart example: n = " << n << " and d = " << d << endl;

	//=========insert data into an AdaptiveHistogramValidation object=========//
	//containers for output needed
	vector<real> IAEV;
	vector<int> NumLeafNodesV;
	vector<double> timings;
	double timing = 0;
	
	// stopping criteria 
	bool stopCrit = false; // if true, allow for flag checking 

	//maximum number of leaf nodes allowed
	int trainCount = n; 
	int holdOutCount = N - n;
	cout << trainCount << " training data and " 
			<< holdOutCount << " validation data inserted." << endl; 
	//size_t maxLeafNodes = int(trainCount/log(trainCount*(2*d+1))); // temporarily
	size_t maxLeafNodes = 1000000;
	cout << "max leaf nodes: " << maxLeafNodes << endl;

	// indicators
   bool successfulInsertion1 = false;   
   bool successfulPQSplit1 = false;   
   
   // comparison objects
   CompCountVal compCount;
   CompVolVal compVol;
   SplitNever sn;   

	//what is this? 
	//int finalK = 1;
	//CritLargestCount_LTEV critCount(finalK);
	size_t minChildPoints = 0;
	CritLeaves_GTEV critToStop(critLeaves);
	cout << "split until there are " << critLeaves << " leaf nodes" << endl;
	
	//container for scheffe tournament candidates
	vector<AdaptiveHistogramValidation> optHist;
	
	cout << "PQ1" << endl;
	// Put the data from the container into the histogram  
	//AdaptiveHistogramValidation myHistVal1(pavingBoxEst);
	AdaptiveHistogramValidation optHist1(pavingBoxEst);
	int minTheta = 0;
	int m = 1;
	
	successfulInsertion1 = FALSE;
	//myHistVal1.insertFromRVecForHoldOut
		//				(*theDataPtr, sn, holdOutCount, NOLOG);
/*
	if (successfulInsertion1) {
		clock_t start, end;
		start = clock();
				
		successfulPQSplit1 = myHistVal1.prioritySplitAndEstimate
							(compCount, critToStop, NOLOG, 
							minChildPoints, 0.0, stopCrit, estimate, 
							dataSeed, 
							maxLeafNodes, maxCheck, minTheta);
	
		
		end = clock();
		timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
		cout << "Computing time : " << timing << " s."<< endl;
	} //temporary
	*/
		//if (successfulPQSplit1) { 
			//cout << "Getting histogram with " << minTheta << " splits." << endl;
			//temporary: set up a for loop to get the IAE of each histogram
			vector<real> vecIAE; 
			optHist1.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
							
			for (size_t i = startLeaves; i <= critLeaves; i++){
				cout << i << endl;				
				size_t critL = i;
				CritLeaves_GTEV critCountOpt(critL);
				clock_t start, end;
				start = clock();
				
				optHist1.prioritySplit(compCount, critCountOpt, NOLOG, minChildPoints, 0.0, maxLeafNodes); 
		
				PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(optHist1);
				real IAE = estimate.getIAE(*tempPCF);
				delete tempPCF; 
				
				vecIAE.push_back(IAE);
				
				end = clock(); 
				timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
				cout << "Computing time : " << timing << " s."<< endl;
		
				
				
				//output IAE to .txt file
/*			   ostringstream stm;
			   stm << i;	
			   string outputFileName;// for output file
			   outputFileName = "MappedIAE";
			   outputFileName += stm.str();
			   outputFileName += ".txt";
			   oss.open(outputFileName.c_str());
			   oss << IAE << endl;
			   oss << flush;
			   oss.close();
			   std::cout << "IAE output to " << outputFileName << endl;*/
			}
									
           //output IAE to .txt file
	//		ostringstream stm;
			stm << dataSeed;
			string outputFileName;// for output file
			outputFileName = "MappedIAE";
			outputFileName += stm.str();
			outputFileName += ".txt";
			oss.open(outputFileName.c_str());
			for (size_t i = 0; i < vecIAE.size(); i++){
				oss << vecIAE[i] << endl;
			}
			oss << flush;
			oss.close();
			std::cout << "IAE output to " << outputFileName << endl;
			
		//} temporary
	// } temporary
    //==============end of MDE=======================================//

	delete theDataPtr;

	return 0;
} // end of MCMC test program

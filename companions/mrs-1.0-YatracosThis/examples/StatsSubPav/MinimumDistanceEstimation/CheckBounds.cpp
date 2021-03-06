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
\brief Get histogram estimates for mapped Gaussian densities using minimum distance estimation (MDE) with hold-out (Devroye and Lugosi, 2001).
*/


#include "histall.hpp"  // headers for the histograms
#include "intervalmappedspnode_measurers.hpp" // ordering for pq split
#include "functionestimator_interval.hpp"
#include "piecewise_constant_function.hpp"  

#include "GaussianFobj.hpp" //function estimator object for Gaussian densitiesi
#include "UniformFobj.hpp"
#include "toolz.hpp"

#include <vector>
#include <algorithm>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iostream>

#include <limits> // to use negative infinity

#include "testDenCommon.hpp" // to use density testing tools
#include "testDenTools.hpp"
#include "mdeTools.hpp"

// to use assert
#include "assert.h"

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main(int argc, char* argv[])
{
	// User-defined parameters------------------//
	if ( argc < 7 ) {
		cerr << "Syntax: " << argv[0] << 
		" data seed, d, n, holdOutPercent, maxLeavesEst, critLeaves, num_checks " << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "data seed, d, n, holdOutPercent, maxLeavesEst, critLeaves, num_checks");
	}

	int dataSeed = atoi(argv[1]); // seed for data generation
	int d = atoi(argv[2]);  // dimension
	const int n = atoi(argv[3]);  // number of points to generate
	double holdOutPercent = atof(argv[4]);
	size_t maxLeavesEst = atoi(argv[5]);  // number of leaves in estimator
	size_t critLeaves = atoi(argv[6]); //maximum number of leaves for PQ to stop splitting 
	int num_checks = atoi(argv[7]); // check k histograms
	//size_t num_iters = atoi(argv[8]); // ...to zoom in
	
	cout << argv[0] << " : process id is " << getpid() << std::endl;
	// End of user-defined parameters--------//

	// string formatting for output purposes
	ofstream oss;       // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(10);
	ostringstream stm_seed, stm_d, stm_n;
	stm_seed << dataSeed; // index the txt file produced by stm
	stm_d << d;
	stm_n << n;

	// Set up a random number generator and use mt19937 for generator
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
	gsl_rng_set (r, dataSeed); // change the seed
	cout << "Data seed is " << dataSeed << endl;

	// Make the function estimator--------//
	cout << "\nMake the function estimator to " << maxLeavesEst << " leaves" << endl;

	//data generating partition
	ivector pavingBoxEst(d);
	interval pavingInterval(0,1);
	for(int i=1; i <= d; i++) { pavingBoxEst[i] = pavingInterval; }

	// specify function object (from /examples/MappedTargetsTrunk)
	//GaussianFobj fobj;
	UniformFobj fobj(pavingBoxEst);

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

	// start clock to record time for pq split in estimate
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

	// start clock to record time for hull propogation and merging up
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

	// optional - remove comments to output function estimator 
	//estimate.outputToTxtTabs("PCF.txt");	
	
	// End of making function estimator--------//

	// Use PiecewiseConstantFunction to generate data, supplying our own rng---//
	cout << "\nGenerating data for simulation" << endl;

	RVecData* theDataPtr = new RVecData;   // a container for all the points generated

	// start clock to record time taken to simulate data
	clock_t startData = clock();

	// Gaussian data
	estimate.simulateData(*theDataPtr, n, r);

	// stop recording time here
	clock_t endData = clock();	
	double timingData = ((static_cast<double>(endData - startData)) / CLOCKS_PER_SEC);
	cout << "Computing time for simulating data: " << timingData << " s."<< endl;

	cout << (*theDataPtr).size() << " points generated" << endl;

	// optional - remove comments to output simulated data 	
	/*
	string dataFileName = "simulated_data_for_bounds";
	dataFileName += stm.str(); 
	dataFileName += ".txt"; 
	oss.open(dataFileName.c_str());
	for (size_t i = 0; i < n; i++) { 
		for (size_t j = 1; j <= d; j++) {
				oss << (*theDataPtr)[i][j] << "\t";
		}
		oss << "\n";
	}
	oss << flush;
	oss.close();
	cout << "Simulated data written to  " << dataFileName << endl;
	*/
	// End of generating data--------//
	
	// Minimum distance estimation with hold-out--------//
	cout << "\nRunning minimum distance estimation with hold-out..." << endl;
	
	int holdOutCount = round(n*holdOutPercent);
	cout << holdOutCount << " points held out." << endl; 

	// parameters for function insertRVectorForHoldOut()
	SplitNever sn; 

	// parameters for getMDETheorem
	CompCountVal compCount; 
	CritLeaves_GTEV he(critLeaves); //the PQ will stop after critLeaves are reached
	size_t minChildPoints = 0;
	size_t maxLeafNodes = 1000000; 
	
	vector<int> sequence; //to store all the thetas
	size_t startLeaves = 0; 
	sequence.push_back(startLeaves + 1);
	sequence.push_back(critLeaves);
	//sequence to be used
	int increment = (critLeaves-startLeaves)/(num_checks);
	cout << "Increment by : " << increment << endl;
	int temp = startLeaves;
	getSequence(sequence, temp, critLeaves, increment);
	//for ( vector<int>::iterator it = sequence.begin(); it != sequence.end(); it++)
	
	vector<double>* vecMaxDeltaTheta = new vector<double>;	
	vector<real>* vecMaxDelta = new vector<real>;	
	vector<real>* vecIAEHoldOut = new vector<real>;		
	vector<real>* vecIAEAllPoints = new vector<real>;	
	
	// to get the required values for Theorem 2	
	// start the clock here
	double timing = 0;
	clock_t start, end;
	start = clock();

	// insert simulated data into an AdaptiveHistogramValidation object
	AdaptiveHistogramValidation myHistVal(pavingBoxEst);
	myHistVal.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
			
 	//run MDE
 	myHistVal.getMDETheoremValues
	 					(compCount, he, NOLOG, 
	 					minChildPoints, 0.0, estimate, 
	 					maxLeafNodes, sequence,	
	 					*vecMaxDeltaTheta, *vecMaxDelta,
	 					*vecIAEHoldOut); 
		
	end = clock();
	timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
	cout << "Computing time for MDE: " << timing << " s."<< endl;
	
	// to get the IAEs for Theorem 3
	cout << "\nGetting the IAEs needed for Theorem 3..." << endl;
	// start the clock here
	start = clock();

	// insert simulated data into an AdaptiveHistogram object
	AdaptiveHistogram myHist(pavingBoxEst);
	myHist.insertFromRVec(*theDataPtr, NOLOG);
	
	CompCount compCount1; 
	CritLeaves_GTE he1(critLeaves); //the PQ will stop after critLeaves are reached

 	myHist.prioritySplitMappedIAE
	 					(compCount1, he1,
             NOLOG, maxLeafNodes,
  					 estimate, *vecIAEAllPoints, sequence); 
	
	end = clock();
	timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
	cout << "Computing time for prioritySplit: " << timing << " s."<< endl;
			
	//find the minimum Delta theta
	double minDeltaTheta = *min_element((*vecMaxDeltaTheta).begin(), (*vecMaxDeltaTheta).end());	

	//find the position of the minimum delta
	size_t minPos = min_element((*vecMaxDeltaTheta).begin(), (*vecMaxDeltaTheta).end()) - (*vecMaxDeltaTheta).begin();
	int numLeavesDelta = sequence[minPos];
				
	//get the IAE using vecIAEHoldOut
	real IAEforMinDelta = (*vecIAEHoldOut)[minPos];
		
	//get minimum IAEHoldOut
	real minIAE = *min_element((*vecIAEHoldOut).begin(), (*vecIAEHoldOut).end());
		
	//find the position of the minimum IAEHoldOut	
	minPos = min_element((*vecIAEHoldOut).begin(), (*vecIAEHoldOut).end()) - (*vecIAEHoldOut).begin();
	int numLeavesIAE = sequence[minPos];
	
	//find the maximum Delta
	real maxDelta = *max_element((*vecMaxDelta).begin(), (*vecMaxDelta).end());	
	
	//get minimum IAEAllPoints
	real minIAEAllPoints = *min_element((*vecIAEAllPoints).begin(), (*vecIAEAllPoints).end());
		
	//find the position of the minimum IAEHoldOut	
	minPos = min_element((*vecIAEAllPoints).begin(), (*vecIAEAllPoints).end()) - (*vecIAEAllPoints).begin();
	int numLeavesIAEAllPoints = sequence[minPos];

	cout << "The minimum Delta theta is " << minDeltaTheta << " at " << numLeavesDelta << " leaf nodes with IAE" << IAEforMinDelta <<endl;
	cout << "The minimum IAEHoldOut is"  << minIAE << " at " << numLeavesIAE << " leaf nodes." << endl;
	cout << "The minimum IAE for all points is"  << minIAEAllPoints << " at " << numLeavesIAEAllPoints << " leaf nodes." << endl;
	cout << "The maximum Delta value is " << maxDelta << endl;
	
	// check theorem 2
	bool check = IAEforMinDelta <= 3*minIAE + 4*maxDelta;
	cout << "RHS: " << IAEforMinDelta << "\t" << "LHS: " << 3*minIAE + 4*maxDelta << endl;
	cout << boolalpha;
	cout << "Inequality for Theorem 2 is " << check << "." << endl;
	
	// optional - remove comments to txt file
	string outputName;
	outputName = "uniform_";
	outputName += stm_d.str();
	outputName += "d_";
	outputName += stm_n.str();
	outputName += "n_theorem2_check";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	oss << IAEforMinDelta << "\t" << (3*minIAE + 4*maxDelta) << "\t" << check << endl;
	oss << flush;
	oss.close();
	cout << "Checks for Theorem 2 output to " << outputName << endl;
	
	// optional - remove comments to txt file
	outputName = "uniform_";
	outputName += stm_d.str();
	outputName += "d_";
	outputName += stm_n.str();
	outputName += "n_theorem3_iaes";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	oss << IAEforMinDelta << "\t" << minIAEAllPoints << endl;
	oss << flush;
	oss.close();
	cout << "Required values for Theorem 3 output to " << outputName << "\n" << endl;

		// optional - remove comments to output the deltas and iaes to txt
	// outputs are: delta_theta   iaeholdout    iaeallpoints 
	/*
	outputName = "delta_theta_and_iaes";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	for (size_t i = 0; i < (*vecMaxDeltaTheta).size(); i++){
			oss << (*vecMaxDeltaTheta)[i] << "\t" << (*vecIAEHoldOut)[i] << "\t" << (*vecIAEAllPoints)[i] << endl;
	}		
	oss << flush;
	oss.close();
	cout << "Delta theta values and IAEs output to " << outputName << endl;
	*/
	
	// optional - remove comments to output the deltas to txt
	/*
	outputName = "delta";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	for (size_t i = 0; i < (*vecMaxDelta).size(); i++){
			oss << (*vecMaxDelta)[i] << endl;
	}		
	oss << flush;
	oss.close();
	cout << "Delta values output to " << outputName << endl;
	*/
		
	try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow
			
	//delete pointers;
	delete vecIAEHoldOut, vecIAEAllPoints;
	delete vecMaxDelta, vecMaxDeltaTheta;	
	delete theDataPtr;
	
		
	return 0;

} // end of program

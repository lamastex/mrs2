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

#include "UniformFobj.hpp" //function estimator object for Gaussian densities

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
		" dataSeed d n maxLeavesEst critLeaves maxCheck" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "data seed, d, n, maxLeavesEst, critLeaves, num_checks, num_iters");
	}

	int dataSeed = atoi(argv[1]); // seed for data generation
	int d = atoi(argv[2]);  // dimension
	const int n = atoi(argv[3]);  // number of points to generate
	double holdOutPercent = atof(argv[4]);
	size_t maxLeavesEst = atoi(argv[5]);  // number of leaves in estimator
	size_t critLeaves = atoi(argv[6]); //maximum number of leaves for PQ to stop splitting 
	int num_checks = atoi(argv[7]); // check k histograms
	size_t num_iters = atoi(argv[8]); // ...to zoom in
	
	cout << argv[0] << " : process id is " << getpid() << std::endl;
	// End of user-defined parameters--------//

	//--future?
	/*MIXSHAPE="1,1"
	MIXSHAPE="2,2,1"
	MIXSHAPE="3,3,2,1"
	MIXSHAPE="3,4,4, 2, 2, 3, 3"
	*///----

	// string formatting for output purposes
	ofstream oss;         // ofstream object
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
  interval pavingInterval(0,1); //standard uniform distribution
	for(int i=1; i <= d; i++) { pavingBoxEst[i] = pavingInterval; }

	// specify function object (from /examples/MappedTargetsTrunk)
	UniformFobj fobj(pavingBoxEst);
	FunctionEstimatorInterval estimator(pavingBoxEst, fobj);

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
	/*estimate.outputToTxtTabs("PCF.txt");
	string Integral = "Integral.txt";
	oss.open(Integral.c_str());
	oss << before << "\t" << estimate.getTotalIntegral() << endl;
	oss << flush;
	oss.close();
	*/
	// End of making function estimator--------//

	// Use PiecewiseConstantFunction to generate data, supplying our own rng---//
	cout << "\nGenerating data for simulation" << endl;

	RVecData* theDataPtr = new RVecData;   // a container for all the points generated

	// start clock to record time taken to simulate data
	clock_t startData = clock();

	estimate.simulateData(*theDataPtr, n, r);

	// stop recording time here
	clock_t endData = clock();	
	double timingData = ((static_cast<double>(endData - startData)) / CLOCKS_PER_SEC);
	cout << (*theDataPtr).size() << " points generated" << endl;
	cout << "Computing time for simulating data: " << timingData << " s."<< endl;

	//optional - remove comments to output simulated data 
	/*
	string dataFileName = "simulated_uniform_data";
	dataFileName += stm_seed.str(); 
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

	// parameters for prioritySplitAndEstimate
	CompCountVal compCount; 
	CritLeaves_GTEV he(critLeaves); //the PQ will stop after critLeaves are reached
	size_t minChildPoints = 0;
	size_t maxLeafNodes = 1000000; 
	bool computeIAE = FALSE; // do not compute the IAE first
	
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
	//	cout << *it << endl;
		
	cout << "Perform " << num_iters << " iterations" << endl; 
	vector<double>* vecMaxDelta = new vector<double>;	
	vector<real>* vecIAE = new vector<real>;		
	size_t k = 3;
	size_t iters = 0;

	// start the clock here
	double timing = 0;
	clock_t start, end;
	start = clock();


	while ( (increment) > 1 && iters < num_iters && (critLeaves - startLeaves) > num_checks) {				
		cout << "\nIteration " << iters << "......" << endl;

		// insert simulated data into an AdaptiveHistogramValidation object
		AdaptiveHistogramValidation myHistVal(pavingBoxEst);
		myHistVal.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
			
	 	//run MDE
	 	myHistVal.prioritySplitAndEstimate
	 					(compCount, he, NOLOG, 
	 					minChildPoints, 0.0, estimate, 
	 					maxLeafNodes, computeIAE, sequence,	
	 					*vecMaxDelta, *vecIAE); 
	
		//get the best 3 delta max values			
		vector<int> indtop;
		topk(*vecMaxDelta, indtop, 3);
		(*vecMaxDelta).clear();
		(*vecIAE).clear();
		//cout << "Best three indices: " << endl;
		//for ( vector<int>::iterator it = indtop.begin(); it != indtop.end(); it++)
			//*it = position //*sequence[*it] = leaves
		//	{ cout << *it << "\t" << sequence[*it] << endl;}
		
		//update final_sequence
		startLeaves = sequence[indtop[0]];
		critLeaves = sequence[indtop[2]];
		if ( (critLeaves - startLeaves) < num_checks ) 
			{ num_checks = critLeaves - startLeaves; }
		increment = (critLeaves-startLeaves)/(num_checks);
		//cout << " Increment by: " << increment << endl;
		
		temp = startLeaves;
		getSequence(sequence, temp, critLeaves, increment);		
		//cout << "updated sequence: " << endl;
		//for ( vector<int>::iterator it = sequence.begin(); it != sequence.end(); it++)
		//	cout << *it << endl;	
				
		//increment iters
		iters++;
	 } //end of while loop


	//Run MDE with the final sequence after breaking out of the loop	
	cout << "\nRun MDE with the final sequence..." << endl;
	computeIAE = TRUE;
	AdaptiveHistogramValidation finalHist(pavingBoxEst);
	finalHist.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
	finalHist.prioritySplitAndEstimate
	 				(compCount, he, NOLOG, 
	 				minChildPoints, 0.0, estimate, 
	 				maxLeafNodes, computeIAE, sequence,	
	 				*vecMaxDelta, *vecIAE);
						
	end = clock();
	timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
	cout << "Computing time for MDE: " << timing << " s."<< endl;
		
	//find the minimum delta
	double minDelta = *min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end());	
	
	//find the position of the minimum delta
	size_t minPos = min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end()) - (*vecMaxDelta).begin();
	int numLeavesDelta = sequence[minPos];
				
	//get the IAE using vecIAE
	real IAEforMinDelta = (*vecIAE)[minPos];
		
	// get minimum IAE
	real minIAE = *min_element((*vecIAE).begin(), (*vecIAE).end());
		
	//find the position of the minimum IAE	
	int numLeavesIAE = sequence[min_element((*vecIAE).begin(), (*vecIAE).end()) - (*vecIAE).begin()];
	//cout << (*vecIAE).size() << "\t" << (*vecMaxDelta).size() << endl;
	
	// difference of minIAE and IAEforMinDelta
	real diffIAE = abs(IAEforMinDelta - minIAE);
	
	cout << "The minimum max delta is " << minDelta << " at " << numLeavesDelta << " leaf nodes with IAE" << IAEforMinDelta << endl;
	cout << "The minimum IAE is"  << minIAE << " at " << numLeavesIAE << " leaf nodes." << endl;
	
		// optional - remove comments to output IAE to txt file
	string outputName;
	outputName = "uniform_";
	outputName += stm_d.str();
	outputName += "d_";
	outputName += stm_n.str();
	outputName += "n_iaes_and_diff";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	oss << IAEforMinDelta << "\t" << minIAE << "\t" << diffIAE << endl;
	oss << flush;
	oss.close();
	cout << "Main results output to " << outputName << "\n" << endl;
		
	// optional - remove comments to output the sequence of leaf nodes
	/*outputName = "sequence";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	for (size_t i = 0; i < (sequence).size(); i++){
		oss << (sequence)[i] << endl;
	}			 
	oss << flush;
	oss.close();
	cout << "Sequence of histograms used output to " << outputName << endl;

	// optional - remove comments to output the deltas to txt
	outputName = "deltas";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	for (size_t i = 0; i < (*vecMaxDelta).size(); i++){
			oss << (*vecMaxDelta)[i] << endl;
	}		
	oss << flush;
	oss.close();
	cout << "Delta values output to " << outputName << endl;
	
	// optional - remove comments to output the IAEs to txt
	outputName = "iaes";
	outputName += stm_seed.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	for (size_t i = 0; i < (*vecIAE).size(); i++){
			oss << (*vecIAE)[i] << endl;
	}		
	oss << flush;
	oss.close();
	cout << "IAEs output to " << outputName << endl;
	*/

	try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow
		
	//delete pointers;
	delete vecIAE;
	delete vecMaxDelta;	
	delete theDataPtr;
	
		
	return 0;

} // end of program

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

#include "GaussianFobj.hpp" //function estimator object for Gaussian densities
#include "toolz.hpp"

#include <vector>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iostream>

#include <limits> // to use negative infinity

#include "testDenCommon.hpp" // to use density testing tools
#include "testDenTools.hpp"

// to use assert
#include "assert.h"

using namespace cxsc;
using namespace std;
using namespace subpavings;

// return a vector of the top k indices of a
// ?.....more desc on what this is for
void topk(vector<double> a, vector<int> & indtop, size_t k){
	multimap<double, size_t> m; // mapping from value to its index
	vector<double>::iterator it;

	for (it = a.begin(); it != a.end(); ++it)
		m.insert(make_pair(*it, it - a.begin()));

	multimap<double, size_t>::iterator itm; // mapping from value to its index
	size_t indx=0;
	double val =0;
	for (itm = m.begin(); itm != m.end(); ++itm){
		//cout << itm->first <<" , "<< itm->second << endl;
		//if (itm->first != val) { 
			indtop.push_back(itm->second);
			indx++;
		//}
		//val = itm->first;
		if ( indx == k) break; 
	}	
	sort(indtop.begin(), indtop.end());
} // end of topk

int main(int argc, char* argv[])
{
	// User-defined parameters------------------//
	if ( argc < 7 ) {
		cerr << "Syntax: " << argv[0] << 
		" dataSeed d n maxLeavesEst critLeaves maxCheck" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "d n states, symmetryIndicator");
	}

	int dataSeed = atoi(argv[1]); // seed for data generation
	int d = atoi(argv[2]);  // dimension
	const int n = atoi(argv[3]);  // number of points to generate
	size_t maxLeavesEst = atoi(argv[4]);  // number of leaves in estimator
	size_t critLeaves = atoi(argv[5]); //maximum number of leaves for PQ to stop splitting 
	int maxCheck = atoi(argv[6]); //? maximum number of checks based on delta
	  
	//   bool flexi = atoi(argv[7]); //?
	//   size_t startLeaves = atoi(argv[8]); //?
	//   int padding = atoi(argv[9]); //?
	// 	size_t num_iters = atoi(argv[10]); //?
	  
	//    vector<int> sequence; //? input sequence
	//    if (argc > 10) {
	// 	for (size_t i=11; i<argc; i++) {
	// 		//cout << argv[i] << endl;
	// 		sequence.push_back(atoi(argv[i]));
	// 	}
	// }
	 //cout << sequence.size() << endl;

	cout << argv[0] << " : process id is " << getpid() << std::endl;
	// End of user-defined parameters--------//

	// string formatting for output purposes
	ofstream oss;       // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(10);
	ostringstream stm;
		stm << dataSeed; // index the txt file produced by stm

	// Set up a random number generator and use mt19937 for generator
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
	gsl_rng_set (r, dataSeed); // change the seed
	cout << "Data seed is " << dataSeed << endl;

	// Make the function estimator--------//
	cout << "\nMake the function estimator to " << maxLeavesEst << " leaves" << endl;

	// specify function object (from /examples/MappedTargetsTrunk)
	GaussianFobj fobj;

	//data generating partition
	ivector pavingBoxEst(d);
	interval pavingInterval(-5,5);
	for(int i=1; i <= d; i++) { pavingBoxEst[i] = pavingInterval; }

	FunctionEstimatorInterval estimator(pavingBoxEst, fobj);
	cout << estimator << endl;

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
	/*
	estimate.outputToTxtTabs("PCF.txt");
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

	// Gaussian data
	int N = 1.5*n;
	estimate.simulateData(*theDataPtr, N, r);

	// stop recording time here
	clock_t endData = clock();	
	double timingData = ((static_cast<double>(endData - startData)) / CLOCKS_PER_SEC);
	cout << "Computing time for simulating data: " << timingData << " s."<< endl;

	cout << (*theDataPtr).size() << " points generated" << endl;

	// optional - remove comments to output simulated data 
	/*
	string dataFileName = "SimulatedData";
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
	cout << "Simulated data written to  " << dataFileName << endl;
	*/
	// End of generating data--------//

	int trainCount = n; 
	int holdOutCount = N-n;
	cout << n << " training data and " << holdOutCount << " validation data inserted." << endl; 
	
	// Minimum distance estimation with hold-out--------//
	cout << "\nRunning minimum distance estimation with hold-out..." << endl;

	
	// parameters for XX function
	CritLeaves_GTEV critToStop(critLeaves); // ? 
	cout << "split until there are " << critLeaves << " leaf nodes" << endl;

  CompCountVal compCount;  // splut according to ...
  SplitNever sn; //??

	size_t maxLeafNodes = 1000000; //what is this for vs critLeaves
	bool stopCrit = flexi; // stopping criteria for ?
	//temp: flexi scheme or every k-th - what is flexi scheme?
	// if true, allow for flag checking 
	int finalK = 1; // ?
	size_t minChildPoints = 0; //?
	int minTheta = 0; //?
	int m = 1; //?

	vector<int> final_sequence; //to store all the thetas
	final_sequence.push_back(startLeaves);
	final_sequence.push_back(critLeaves);
	size_t k = 3; //take the best three

	// //sequence to be used
	// int increment = (critLeaves-startLeaves)/(maxCheck);
	// //cout << "Increment by : increment" << endl;
	// int temp = startLeaves;
	// while ( temp < critLeaves) {
	// 	temp += increment;
	// 	final_sequence.push_back(temp); 
	// }

	// sort(final_sequence.begin(), final_sequence.end());
	// final_sequence.erase( unique( final_sequence.begin(), final_sequence.end() ), final_sequence.end() );
	// //for ( vector<int>::iterator it = final_sequence.begin(); it != final_sequence.end(); it++)
	// //	cout << *it << endl;

	// //adaptive sequencing starts here - what is adaptive sequencing?!
	double timing = 0;
	clock_t start, end;
	start = clock();

	// cout << "Perform " << num_iters << " iterations" << endl;
	// vector<double>* vecMaxDelta = new vector<double>;	
	// vector<real>* vecIAE = new vector<real>;		
	// size_t i = 0;

	while ( (increment) > 1 && i < num_iters && (critLeaves - startLeaves) > maxCheck) {			
	
		cout << "Iteration: " << i << endl;

		// insert simulated data into an AdaptiveHistogramValidation object
 		AdaptiveHistogramValidation myHistVal1(pavingBoxEst);
 		myHistVal1.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
		
		// optional - remove comments to output ...
		/*
		dataSeed = i*1;
		ostringstream stm;
		stm << dataSeed;
		string sequenceName;
		sequenceName = "Sequence";
		sequenceName += stm.str();
		sequenceName += ".txt";
		oss.open(sequenceName.c_str());
		for ( vector<int>::iterator it = final_sequence.begin(); it != final_sequence.end(); it++) {
			oss << *it << endl;
		}			 
		oss << flush;
		oss.close();
		*/
		
	// 	//run MDE
	// 	myHistVal1.prioritySplitAndEstimate
	// 					(compCount, critToStop, NOLOG, 
	// 					minChildPoints, 0.0, stopCrit, estimate, 
	// 					dataSeed, 
	// 					maxLeafNodes, maxCheck, minTheta, final_sequence,	
	// 					*vecMaxDelta, *vecIAE);

	// 	//get the best 3 delta max values			
	// 	vector<int> indtop;
	// 	topk(*vecMaxDelta, indtop, 3);
	// 	(*vecMaxDelta).clear();
	// 	(*vecIAE).clear();
	// 	//cout << "Best three indices: " << endl;
	// 	//for ( vector<int>::iterator it = indtop.begin(); it != indtop.end(); it++)
	// 	//	{ cout << *it << "\t" << final_sequence[*it] << endl;}

		
	// 	//adjust indtop with padding 
	// 	//if ( final_sequence[indtop[0]] == 1) { final_sequence[indtop[0]]; } 
	// 	//else { final_sequence[indtop[0]] = final_sequence[indtop[0]] - padding; }
	// 	//final_sequence[indtop[2]] = final_sequence[indtop[2]] + padding;
	// 	//cout << "updated sequence with padding: " << endl;			
		
	// 	//update final_sequence
	// 	startLeaves = final_sequence[indtop[0]];
	// 	critLeaves = final_sequence[indtop[2]];
	// 	if ( (critLeaves - startLeaves) < maxCheck ) 
	// 	{ maxCheck = critLeaves - startLeaves; }
		
	// 	increment = (critLeaves-startLeaves)/(maxCheck);
	// 	cout << startLeaves << " Increment by: " << increment << " " << critLeaves << endl;
		
	// 	temp = startLeaves;
	// 	while ( temp < critLeaves) {
	// 		temp += increment;
	// 		//cout << "temp: "<< temp << endl;
	// 		final_sequence.push_back(temp); 
	// 	}
	// 	sort(final_sequence.begin(), final_sequence.end());
	// 	final_sequence.erase( unique( final_sequence.begin(), final_sequence.end() ), final_sequence.end() );
		
	// 	//cout << "updated sequence: " << endl;
	// 	//for ( vector<int>::iterator it = final_sequence.begin(); it != final_sequence.end(); it++)
	// 	//	cout << *it << endl;	
			
	// 	//increment i
	// 	i++;
	// } //end of while loop

	// //Run MDE with the final sequence after breaking out of the loop	
	// AdaptiveHistogramValidation myHistVal1(pavingBoxEst);
	// myHistVal1.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);

	// //run MDE
	// //dataSeed = i;
	// myHistVal1.prioritySplitAndEstimate
	// 				(compCount, critToStop, NOLOG, 
	// 				minChildPoints, 0.0, stopCrit, estimate, 
	// 				dataSeed, 
	// 				maxLeafNodes, maxCheck, minTheta, final_sequence,	
	// 				*vecMaxDelta, *vecIAE);
						
	// //temporary misusing maxcheck so that will only compute the IAEs once		 
	end = clock();
	timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
	cout << "Computing time for MDE: " << timing << " s."<< endl;
		
	// //find the minimum delta
	// //double minDelta = *min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end());	

	// //find the position of the minimum delta
	// size_t minPos = min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end()) - (*vecMaxDelta).begin();
	// int numLeavesDelta = final_sequence[minPos];
		
	// //get the IAE using vecIAE
	// real IAEforMinDelta = (*vecIAE)[numLeavesDelta - 1];
		
	// AdaptiveHistogramValidation* optHist = new AdaptiveHistogramValidation(pavingBoxEst);
	// (*optHist).insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
	// CritLeaves_GTEV critCountOpt(numLeavesDelta);
	// (*optHist).prioritySplit(compCount, critCountOpt, NOLOG, minChildPoints, 0.0, maxLeafNodes); 	
	// PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(*optHist);
	// IAEforMinDelta = estimate.getIAE(*tempPCF);
	// delete optHist; 

	// get minimum IAE
	real minIAE = *min_element((*vecIAE).begin(), (*vecIAE).end());
		
	//find the position of the minimum IAE	
	int numLeavesIAE = min_element((*vecIAE).begin(), (*vecIAE).end()) - (*vecIAE).begin() + 1;

	//delete pointers;
	delete vecIAE;
	delete vecMaxDelta;	
	delete theDataPtr;
	delete tempPCF;

	try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow

	// output results to txt file
	cout << IAEforMinDelta << "\t" << numLeavesDelta << "\t" << minIAE << "\t" << numLeavesIAE << endl;
	string outputName;
	outputName = "results";
	outputName += stm.str();
	outputName += ".txt";
	oss.open(outputName.c_str());
	oss << IAEforMinDelta << "\t" << numLeavesDelta << "\t" << minIAE << "\t" << numLeavesIAE << "\t" << timing << endl;
	oss << flush;
	oss.close();
  cout << "Error computations output to " << outputName << endl;
	
	// optional - remove comments to output ...
	/*
	string sequenceName;
	sequenceName = "Sequence";
	sequenceName += stm.str();
	sequenceName += ".txt";
	oss.open(sequenceName.c_str());
	for (size_t i = 0; i < (final_sequence).size(); i++){
		oss << (final_sequence)[i] << endl;
	}			 
	oss << flush;
	oss.close();
	*/	

  // End of MDE---------//

		
	return 0;

} // end of program
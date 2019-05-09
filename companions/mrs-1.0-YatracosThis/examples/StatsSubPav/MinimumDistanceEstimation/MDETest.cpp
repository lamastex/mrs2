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
\brief Get histogram estimates for a dataset from a .txt file 
* - the data set should be without header
* - may be delimited by tabs or spaces
*/

#include "histall.hpp"  // headers for the histograms
#include "intervalmappedspnode_measurers.hpp" // ordering for pq split
#include "functionestimator_interval.hpp"
#include "piecewise_constant_function.hpp"  

#include "toolz.hpp"

#include <vector>
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
	if ( argc < 5 ) {
		cerr << "Syntax: " << argv[0] << 
		" dataSeed d n maxLeavesEst critLeaves maxCheck" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + 
		"inputFileName, holdOutPercent, critLeaves, num_checks, num_iters, buildRootBox");
	}

	string inputFileName = argv[1];
	string outputFileName = argv[2];
	double holdOutPercent = atof(argv[3]);
	size_t critLeaves = atoi(argv[4]); //maximum number of leaves for PQ to stop splitting 
	int num_checks = atoi(argv[5]); // check k histograms
	size_t num_iters = atoi(argv[6]); // to zoom in
	

	cout << "Processing file " << inputFileName << endl;
	// End of user-defined parameters--------//

	// string formatting for output purposes
	ofstream oss;       // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(10);
	ostringstream stm_out;
	stm_out << outputFileName;
	
	// uncomment this if want to define own root box	
	/*int d = 1; //change the value depending on dimension required
	ivector pavingBox(d); 
	interval pavingInterval(-5, 5);
	for(int i=1; i <= d; i++) { pavingBox[i] = pavingInterval; }
	AdaptiveHistogramValidation finalHist(pavingBox);
	*/

	// comment this line if want to define own root box	
	AdaptiveHistogramValidation finalHist; 	// the root box will be made based on the data set given	
		
	//insert data from text into a RVecData container and store as pointer
	vector<size_t> numNodes; 
	RVecData* theDataPtr = new RVecData; 
  	bool successfulInsertion = false;
	successfulInsertion = finalHist.insertRvectorsFromTxtForHoldOut(
	inputFileName, *theDataPtr, numNodes, holdOutPercent, 0, NOLOG);
	ivector pavingBox = finalHist.getRootBox();
	
	if (successfulInsertion) {	
		// Minimum distance estimation with hold-out--------//
		cout << "\nRunning minimum distance estimation with hold-out..." << endl;
		
		int n = (*theDataPtr).size(); 
		int holdOutCount = round(n*holdOutPercent);
		cout << holdOutCount << " points held out." << endl; 
	
		// parameters for function insertRVectorForHoldOut()
		SplitNever sn; 
	
		// parameters for prioritySplitAndEstimate
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
			//cout << *it << endl;
			
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
			AdaptiveHistogramValidation myHistVal(pavingBox); 
			myHistVal.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
				
		 	//run MDE
		 	myHistVal.prioritySplitAndEstimatePlain
		 					(compCount, he, NOLOG, 
		 					minChildPoints, 0.0, 
		 					maxLeafNodes, sequence,	
		 					*vecMaxDelta); 
		
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
		finalHist.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);
		finalHist.prioritySplitAndEstimatePlain
		 				(compCount, he, NOLOG, 
		 				minChildPoints, 0.0,
		 				maxLeafNodes, sequence,	
		 				*vecMaxDelta);
							
		end = clock();
		timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
		cout << "Computing time for MDE: " << timing << " s."<< endl;
			
		//find the minimum delta
		double minDelta = *min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end());	
	
		//find the position of the minimum delta
		size_t minPos = min_element((*vecMaxDelta).begin(), (*vecMaxDelta).end()) - (*vecMaxDelta).begin();
		int numLeavesDelta = sequence[minPos];
	
		cout << "The minimum max delta is " << minDelta << " at " << numLeavesDelta << " leaf nodes." << endl;
		
		//build an adaptivehistogram based on numLeavesDelta
		AdaptiveHistogram mdeHist(pavingBox);
		mdeHist.insertFromRVec(*theDataPtr, NOLOG);
		CompCount compCountMDE;
		CritLeaves_GTE heMDE(numLeavesDelta);
		mdeHist.prioritySplit(compCountMDE, heMDE, NOLOG, maxLeafNodes);	

		// make mdeHist into a PCF object
		PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(mdeHist);
		
		//output the histogram built using minimum distance estimate
		string outputName = stm_out.str();
		outputName += "_mdehist.txt";
		mdeHist.outputToTxtTabs(outputName);
		cout << "MDE histogram output to " << outputName << endl;
		
		// output the tail probabilties
		outputName = stm_out.str();
		outputName += "_coverage.txt";
		oss.open(outputName.c_str());
		for (size_t i = 0; i < (*theDataPtr).size(); i++){
				oss << (*tempPCF).findCoverage((*theDataPtr)[i]) << endl;
		}		
		oss << flush;
		oss.close();
		cout << "Coverage values output to " << outputName << endl;
		delete tempPCF;

		// optional - remove comments to output the sequence of leaf nodes
		outputName = stm_out.str();
		outputName += "_theta_and_delta_theta.txt";
		oss.open(outputName.c_str());
		for (size_t i = 0; i < (sequence).size(); i++){
			oss << (sequence)[i] <<  (*vecMaxDelta)[i] << endl;
		}			 
		oss << flush;
		oss.close();
		cout << "The delta_theta values and corresponding number of leaf nodes (theta) output to " << outputName << endl;
		
		//delete pointers;
		delete vecIAE;
		delete vecMaxDelta;	
		delete theDataPtr;
	
		
	} // end of successfulInsertion = TRUE
	
	
	return 0;

} // end of program

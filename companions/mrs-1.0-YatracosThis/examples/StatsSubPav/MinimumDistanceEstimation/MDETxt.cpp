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

// return a vector of the top k indices of ...
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
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "data seed, d, n, maxLeavesEst, critLeaves, num_checks, num_iters");
	}

	int dataSeed = atoi(argv[1]); // seed for data generation
	int d = atoi(argv[2]);  // dimension
	const int n = atoi(argv[3]);  // number of points to generate
	size_t maxLeavesEst = atoi(argv[4]);  // number of leaves in estimator
	size_t critLeaves = atoi(argv[5]); //maximum number of leaves for PQ to stop splitting 
	int num_checks = atoi(argv[6]); // check k histograms
	size_t num_iters = atoi(argv[7]); // ...to zoom in
	
	cout << argv[0] << " : process id is " << getpid() << std::endl;
	// End of user-defined parameters--------//

	// string formatting for output purposes
	ofstream oss;       // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(10);
	ostringstream stm;
	stm << dataSeed; // index the txt file produced by stm

	// code for txt input goes here...

	// Set up a random number generator and use mt19937 for generator
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
	gsl_rng_set (r, dataSeed); // change the seed
	cout << "Data seed is " << dataSeed << endl;

	// Minimum distance estimation with hold-out--------//
	cout << "\nRunning minimum distance estimation with hold-out..." << endl;
	
	int trainCount = n; 
	int holdOutCount = N-n;
	cout << n << " training data and " << holdOutCount << " validation data inserted." << endl; 

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
	while ( temp < critLeaves) {
		temp += increment;
	 	sequence.push_back(temp); 
	 }
	sort(sequence.begin(), sequence.end());
	sequence.erase( unique( sequence.begin(), sequence.end() ), sequence.end() );
	//for ( vector<int>::iterator it = sequence.begin(); it != sequence.end(); it++)
		//cout << *it << endl;
		
	cout << "Perform " << num_iters << " iterations" << endl; 
	vector<double>* vecMaxDelta = new vector<double>;	//which delta vector is this?
	vector<real>* vecIAE = new vector<real>;		//which IAE is this?
	
	size_t k = 3; //take the best three  delta - what is best?

	// start the clock here for ...?
	double timing = 0;
	clock_t start, end;
	start = clock();

	size_t iters = 0;

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
	 					*vecMaxDelta, *vecIAE); //don't compute vecIAE here - is it possible?
	
		//get the best 3 delta max values			
		vector<int> indtop;
		topk(*vecMaxDelta, indtop, 3);
		(*vecMaxDelta).clear();
		(*vecIAE).clear();
		//cout << "Best three indices: " << endl;
		//for ( vector<int>::iterator it = indtop.begin(); it != indtop.end(); it++)
			// //*it = position //*sequence[*it] = leaves
		 	//{ cout << *it << "\t" << sequence[*it] << endl;}
		
		//update final_sequence
		startLeaves = sequence[indtop[0]];
		critLeaves = sequence[indtop[2]];
		if ( (critLeaves - startLeaves) < num_checks ) 
		{ num_checks = critLeaves - startLeaves; }
		
		increment = (critLeaves-startLeaves)/(num_checks);
		//cout << " Increment by: " << increment << endl;
		
		temp = startLeaves;
		while ( temp < critLeaves) {
			temp += increment;
			//cout << "temp: "<< temp << endl;
			sequence.push_back(temp); 
		}
		sort(sequence.begin(), sequence.end());
		sequence.erase( unique( sequence.begin(), sequence.end() ), sequence.end() );
		
		//cout << "updated sequence: " << endl;
		//for ( vector<int>::iterator it = sequence.begin(); it != sequence.end(); it++)
			//cout << *it << endl;	
				
		//increment iters
		iters++;
	 } //end of while loop

	//Run MDE with the final sequence after breaking out of the loop	
	cout << "\nRun MDE with the final sequence..." << endl;
	AdaptiveHistogramValidation myHistVal(pavingBoxEst);
	myHistVal.insertFromRVecForHoldOut(*theDataPtr, sn, holdOutCount, NOLOG);

	computeIAE = TRUE;
	myHistVal.prioritySplitAndEstimate
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
	real IAEforMinDelta = (*vecIAE)[numLeavesDelta - 1];
		
	// get minimum IAE
	real minIAE = *min_element((*vecIAE).begin(), (*vecIAE).end());
		
	//find the position of the minimum IAE	
	int numLeavesIAE = min_element((*vecIAE).begin(), (*vecIAE).end()) - (*vecIAE).begin() + 1;

	//delete pointers;
	delete vecIAE;
	delete vecMaxDelta;	
	delete theDataPtr;

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
	oss << IAEforMinDelta << "\t" << numLeavesDelta << "\t" << minIAE << "\t" << numLeavesIAE << endl;
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
		
	return 0;

} // end of program

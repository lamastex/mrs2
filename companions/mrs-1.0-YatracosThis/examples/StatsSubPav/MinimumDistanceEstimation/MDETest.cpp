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

#include "GaussianFobj.hpp" //function estimator object for Gaussian densities
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
	if ( argc < 4 ) {
		cerr << "Syntax: " << argv[0] << 
		" dataSeed d n maxLeavesEst critLeaves maxCheck" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "inputFileName, critLeaves, num_checks, num_iters");
	}

	string inputFileName = argv[1];
	size_t critLeaves = atoi(argv[2]); //maximum number of leaves for PQ to stop splitting 
	int num_checks = atoi(argv[3]); // check k histograms
	size_t num_iters = atoi(argv[4]); // to zoom in
	
	cout << "Processing file " << inputFileName << endl;
	// End of user-defined parameters--------//

	// string formatting for output purposes
	ofstream oss;       // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(10);
	ostringstream stm;
		
	AdaptiveHistogramValidation myHistVal; //the root box will be made based on the data set given
	vector<size_t> numNodes;
  
  bool successfulInsertion = false;
 	//int trainCount = n; 
	//int holdOutCount = N-n;
	//cout << n << " training data and " << holdOutCount << " validation data inserted." << endl; 

	RVecData* theDataPtr = new RVecData; 
	
	//successfulInsertion = myHistVal.insertRvectorsFromTxt(inputFileName, 
		//																										numNodes, 
			//																									NOLOG);
	                                                   	
	if (successfulInsertion) {	
		// /* optional output histogram
		string histFileName = "hist.txt";
		myHistVal.outputToTxtTabs(histFileName);
		//*/
	}
	
	delete theDataPtr;
	
	
	return 0;

} // end of program

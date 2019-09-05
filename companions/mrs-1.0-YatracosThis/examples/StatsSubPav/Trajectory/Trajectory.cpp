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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*! \file Trajectory.cpp
\brief SRPs for enclosures of position data
*/

// to use std::vector
#include <vector>
// to use iterators
#include <iterator>

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams

#include "toolz.hpp"    // toolz headers
#include "histall.hpp"  // headers for the histograms
#include <gsl/gsl_rng.h>

using namespace cxsc;
using namespace std;

//int main()
int main(int argc, char* argv[])
{
	// string formatting
	ofstream oss;         // ofstream object
	oss << scientific;  // set formatting for input to oss
	oss.precision(5);
	bool successfulInsertion = false;
	  
	// input parameters
	string DataFiles = argv[1]; // this is a row vector of filenames
	double Vol = atof(argv[2]); //An approximate volume of the object
	     
	// set up to read in data files
	//create a vector object to store the filenames of simulated data
	vector<string> DataFilesVec;
	string fileName;
	cout << "Reading in file names for simulated data: " << endl;    
	ifstream file; // create file for input
	file.open(DataFiles.c_str());
	// check if this file exists or not
	if ( !file ) { // exit if file doesn't exists'
		cerr << "Could not open " << DataFiles << ". It does not exist." 
		     << endl;
		exit(1);
	}
  // else read in the filenames
	
	// store the filenames in the vector simDataFilesVec
	while ( !file.eof() ) { // read until end of file or error
		file >> fileName;
		cout << fileName << endl; 
		DataFilesVec.push_back(fileName);
	}
	// Somehow an extra line will be read in. Need to delete that extra line.
	DataFilesVec.pop_back();
	
	//container to keep data to make root box
	RVecData* dataPtr;
	dataPtr = new RVecData;

	// put all simulated data into container allData
	cout << "\nPut all data in a container to get rootbox: " << endl;
	for (size_t i = 0;  i < DataFilesVec.size(); i++) {
		cout << DataFilesVec[i] << endl;
		//read into allData
		bool retvalue = readRvectorsFromTxt((*dataPtr), DataFilesVec[i], 0);
		if (retvalue == false) {
			cerr << "Could not open " << DataFilesVec[i] << ". It does not exist." 
				<< endl;
			exit(1);
		}
	}
	
	//Make root box from all the data
	cout << "\n" << endl;
	AdaptiveHistogram* histRoot;
	histRoot = new AdaptiveHistogram;
	histRoot->insertFromRVec((*dataPtr));
	ivector pavingBox = histRoot->getSubPaving()->getBox();
	//find the data dimensions from the first datapoint
	size_t dataDim = Ub(*(*dataPtr).begin()) - Lb(*(*dataPtr).begin()) + 1;
	cout << "Data has " << dataDim << " dimensions." << endl;
	
	delete dataPtr; //we do not need this in memory
	delete histRoot; //we do not need this in memory
	
	// get minimum volume using the same scale as given position data
	AdaptiveHistogramValidation hist(pavingBox);
	cout << "getRootBoxVol" << endl;
	double rootBoxVol = hist.getSubPaving()->nodeVolume();
	double approxDepth = floor(log(rootBoxVol/Vol)/log(2));
	double approxMinVol = rootBoxVol/pow(2,approxDepth);
	cout << "Vol: " << Vol << "\tapproxMinVol: " << approxMinVol << endl; 
	
	// set up for getting histograms
	// split on count and volume to get tightest possible enclosure
	SplitOnKandVol splitVolCount(approxMinVol);

	//collator object for adding trajectories
	AdaptiveHistogramVCollator coll;	
	 
	// get individual trajectories and add into collator
	for (size_t j = 1; j <= DataFilesVec.size(); j++) {
		cout << "================" << j << "======================" << endl;
		ostringstream stm1;
			stm1 << j;
		// tell user which data is being used.
		string inputFileName = DataFilesVec[j];
		cout << "Processing file " << inputFileName << endl;
		// get a count of lines in the txt file
		int dataCount = countLinesInTxt(inputFileName);
		// tell user how many lines there are in the file
		cout << "The file " << inputFileName << " has " << dataCount
				  << " lines in it" << endl << endl;

		// create histograms
		clock_t start, end;
		double timeTaken;
		start=clock();
		cout << "Getting enclosure for this trajectory: " << endl; 
		AdaptiveHistogramValidation myHist(pavingBox);
		vector<size_t> numNodes;
		successfulInsertion = myHist.insertRvectorsFromTxt(inputFileName, 
																				numNodes, 
																			 splitVolCount, NOLOG);
		end=clock();
		timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
		cout << "Computing time : " <<timeTaken<< " s." << endl;
		
		if (successfulInsertion) {
			
			// add into collator
			size_t aggBox = 0;
			cout << "\n Add into collator" << endl;
			coll.addToCollationWithVal(myHist, 2, aggBox);
			
			// /* optional output histogram
			string histFileName = "Hist";
			histFileName += stm1.str();
			histFileName += ".txt";
			myHist.outputToTxtTabs(histFileName);
			// */

		} // end of successful insertion
	} // end of going through each file
		

	// /* optional output histogram
	coll.outputAccumulationToTxtTabs("coll.txt");
	// */
				
	return 0;
} 

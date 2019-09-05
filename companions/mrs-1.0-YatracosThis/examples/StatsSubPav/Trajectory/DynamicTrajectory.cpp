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

/*! \file DynamicTrajectory.cpp
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

int main(int argc, char* argv[])

{
	// string formatting
	ofstream oss;         // ofstream object
  oss << scientific;  // set formatting for input to oss
  oss.precision(5);
  bool successfulInsertion = false;

	// input parameters
  double craftVol = atof(argv[1]); //note this can be put into the for loop if we know each individual craft size
	int starttime = atoi(argv[2]);
	int totalTimeBlock = atoi(argv[3]);
	int totalFlight = atoi(argv[4]);
	
	// this is hardcoded for now  
	// make an Adaptive Histogram object with a specified box
	ivector pavingBox(2);
	/*
	interval pavingInterval1(550, 1350);
	interval pavingInterval2(810, 1230);
  */
 	interval pavingInterval1(-10, 10);
	interval pavingInterval2(-10, 10);
	pavingBox[1] = pavingInterval1;
	pavingBox[2] = pavingInterval2;	  
	cout << "Box is: " << pavingBox << endl;
	 
   //split on k and volume to get tightest possible enclosure
   AdaptiveHistogramValidation hist(pavingBox);
   cout << "getRootBoxVol" << endl;
   double rootBoxVol = hist.getSubPaving()->nodeVolume();
   double approxDepth = floor(log(rootBoxVol/craftVol)/log(2));
   double approxMinVol = rootBoxVol/pow(2,approxDepth);
   cout << "craftVol: " << craftVol << "\tapproxMinVol: " << approxMinVol << endl; 
   SplitOnKandVol splitVolCount(approxMinVol);
	
	//create a vector of AdaptiveHistogramValidation objects
	vector<AdaptiveHistogramValidation> histVec(totalFlight, AdaptiveHistogramValidation(pavingBox)); 

	//create collator objects
	AdaptiveHistogramVCollator thisColl;
	AdaptiveHistogramVCollator currColl;
	AdaptiveHistogramVCollator spaceColl;
		
	for (size_t t=starttime; t < totalTimeBlock; t++) {
		std::ostringstream stm2;
		stm2 << t;

		int checkHist = 1;
		for (size_t j = 1; j <= totalFlight; j++) {
			//name of file to be read in
			std::ostringstream stm1;
			stm1 << j;
			string inputFileName = "Time"; 
			inputFileName += stm2.str();
			inputFileName += "Flight";
			inputFileName += stm1.str();
			inputFileName += ".txt";

			// tell user which data is being used.
			cout << "Processing file " << inputFileName << "\n" << endl;
			
			// get a count of lines in the txt file
			int dataCount = countLinesInTxt(inputFileName);
			
			// insert data into histogram
			AdaptiveHistogramValidation myHist(pavingBox);
			myHist = histVec[j-1];
			
			vector<size_t> numNodes;
			successfulInsertion = myHist.insertRvectorsFromTxt(inputFileName, 
										numNodes, splitVolCount, NOLOG);

			if (successfulInsertion) {
				
				histVec[j-1] = myHist;				
				
				size_t aggBox = 0;
				thisColl.addToCollationWithVal(myHist, 2, aggBox);
				
			} // end of successful insertion
			else { 
				checkHist++; 
				AdaptiveHistogramValidation newHist(pavingBox);
				histVec[j-1] = newHist;	
			}
		} // end of flights

		if (t==starttime) { 
			if ( (checkHist < totalFlight) ){
				thisColl.makeMinimal();
				currColl = thisColl;
				
				// /* optional
				string spaceCollFileName;
				spaceCollFileName = "spaceColl";
				spaceCollFileName += stm2.str();
				spaceCollFileName += ".txt";
				thisColl.outputAccumulationToTxtTabs(spaceCollFileName);	
				// */
			}
			// if there are no flights at starttime
			else { starttime++; }
		}
		
		else { // ( t > starttime) 
			spaceColl = thisColl - currColl;
			spaceColl.makeMinimal(); 
			
			// /* optional
			string spaceCollFileName;
			spaceCollFileName = "spaceColl";
			spaceCollFileName += stm2.str();
			spaceCollFileName += ".txt";
			spaceColl.outputAccumulationToTxtTabs(spaceCollFileName);
			// */ 
						
			currColl = spaceColl;
			thisColl = spaceColl;
			
		}
	} // end of time

	return 0;
} // end of dynamic air traffic example program

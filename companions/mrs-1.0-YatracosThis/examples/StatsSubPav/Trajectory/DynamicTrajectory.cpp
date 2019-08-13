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
	//------------------------string formatting-------------------------------- 
  ofstream oss;         // ofstream object
  oss << scientific;  // set formatting for input to oss
  oss.precision(5);
  bool successfulInsertion = false;

	// input parameters
	string DataFiles = argv[1]; // this is a row vector of filenames
	double Vol = atof(argv[2]); //An approximate volume of the object

	
	
   //get individual trajectories and add into collator 
	size_t starttime = atoi(argv[2]);
	size_t totalTimeBlock = atoi(argv[3]);
	size_t totalFlight = atoi(argv[4]);
	
	  
	 //-------------- make an Adaptive Histogram object with a specified box----  
	 
	int d = atoi(argv[5]); // dimension of the sample data
  ivector pavingBox(d);
  
	interval pavingInterval1(550, 1350);
	interval pavingInterval2(810, 1230);
  pavingBox[1] = pavingInterval1;
  pavingBox[2] = pavingInterval2;
   
	cout << "Box is: " << pavingBox << endl;
	 
   // get minimum volume--------------------------------------
   AdaptiveHistogramValidation hist(pavingBox);
   //0.1
   double craftVol = atof(argv[1]); //note this can be put into the for loop if we know each individual craft size
   cout << "getRootBoxVol" << endl;
   double rootBoxVol = hist.getSubPaving()->nodeVolume();
   double approxDepth = floor(log(rootBoxVol/craftVol)/log(2));
   double approxMinVol = rootBoxVol/pow(2,approxDepth);
   cout << "craftVol: " << craftVol << "\tapproxMinVol: " << approxMinVol << endl; 
 
   //split on k and volume to get tightest possible enclosure
   SplitOnKandVol splitVolCount(approxMinVol);
	
	//vector to store total number of nodes at each spaceColl
	vector<size_t> numAgg;
	vector<double> timings;
	//create totalTimeBlock collator objects
	AdaptiveHistogramVCollator currColl;
	AdaptiveHistogramVCollator updColl;
		
	//create AdaptiveHistogramValidation objects
   vector<AdaptiveHistogramValidation> histVec; 
   for (int i = 0; i < totalFlight; i++) {
		cout << i << endl;
		AdaptiveHistogramValidation myHist(pavingBox);
		histVec.push_back(myHist);
	}
	
	for (size_t t=starttime; t < totalTimeBlock; t++) {
		clock_t start, end;
		
		AdaptiveHistogramVCollator spaceColl;
		std::ostringstream stm2;
		stm2 << t;

		start = clock();
		int checkHist = 0;
		for (size_t j = 1; j <= totalFlight; j++) {
			//name of file to be read in
			std::ostringstream stm1;
			stm1 << j;
			string inputFileName = "Time"; 
			inputFileName += stm2.str();
			inputFileName += "Flight";
			inputFileName += stm1.str();
			inputFileName += ".txt";
			
		/*	if (d==2) {	inputFileName += "xy.txt"; }
		   else if (d==3) { inputFileName += "xyAlt.txt"; }
	      */

			// tell user which data is being used.
			cout << "--------------------------------------------" << endl;
			cout << "j = " << j << endl;
			cout << "Processing file " << inputFileName << endl;
			// get a count of lines in the txt file
			int dataCount = countLinesInTxt(inputFileName);
			// tell user how many lines there are in the file
			cout << "The file " << inputFileName << " has " << dataCount
					<< " lines in it" << endl << endl;
			AdaptiveHistogramValidation myHist(pavingBox);
			myHist = histVec[j-1];
			
			cout << "inserting and constructing histogram:" << endl;
			vector<size_t> numNodes;
			successfulInsertion = myHist.insertRvectorsFromTxt(inputFileName, 
										numNodes, splitVolCount, NOLOG);

			/* optional
			string histFileName;
			histFileName = "HistTime";
			histFileName += stm2.str();
			//histFileName += "Flight";
			//histFileName += stm1.str();
			histFileName += ".txt";
			myHist.outputToTxtTabs(histFileName);
		   */

			if (successfulInsertion) {
				histVec[j-1] = myHist;

				size_t aggBox = 0;
				cout << "adding myHist into collator" << endl;
				updColl.addToCollationWithVal(myHist, 2, aggBox);
			} // end of successful insertion
			else { 
				checkHist++; 
				AdaptiveHistogramValidation newHist(pavingBox);
				histVec[j-1] = newHist;
			}
		} // end of flights

		/* optional
		string collFileName = "updColl";
		collFileName += stm2.str();
		collFileName += ".txt";
		updColl.outputAccumulationToTxtTabs(collFileName);
		*/

		if (t==starttime) { 
			if ( (checkHist < totalFlight) ){
				cout << "get space coll at time " << t << endl;
				currColl = updColl; 
				currColl.makeMinimal();
				//string currCollFileName;
				//currCollFileName = "spaceColl";
				//currCollFileName += stm2.str();
				//currCollFileName += ".txt";
				//currColl.outputAccumulationToTxtTabs(currCollFileName);
				//string spaceCollFileName;
				//spaceCollFileName = "spaceColl";
				//spaceCollFileName += stm2.str();
				//spaceCollFileName += ".txt";
				//currColl.outputAccumulationToTxtTabs(spaceCollFileName);
				numAgg.push_back(currColl.getTotalNodes());
			}
			else { starttime++; }
		}
		
		else { // ( t > starttime) 
			cout << "get space coll at time " << t << endl;
			cout << "getDifference" << endl;
			spaceColl = updColl - currColl;
			//string diffCollFileName;
			//diffCollFileName = "diffColl";
			//diffCollFileName += stm2.str();
			//diffCollFileName += ".txt";
			//spaceColl.outputToTxtTabs("diffCollInd.txt");
			//spaceColl.outputAccumulationToTxtTabs(diffCollFileName);
						
			cout << "Make Minimal:" << endl;
			spaceColl.makeMinimal(); // this collator is the structure we want
			numAgg.push_back(spaceColl.getTotalNodes());
			//output this collator
			//string spaceCollFileName;
			//spaceCollFileName = "spaceColl";
			//spaceCollFileName += stm2.str();
			//spaceCollFileName += ".txt";
			//spaceColl.outputAccumulationToTxtTabs(spaceCollFileName);
						
			//only want the summary of the last column
			AdaptiveHistogramVCollator minimalColl(spaceColl, 2);
			currColl = minimalColl;
			//string currCollFileName;
			//currCollFileName = "currColl";
			//currCollFileName += stm2.str();
			//currCollFileName += ".txt";
			//currColl.outputAccumulationToTxtTabs(currCollFileName);
			updColl = minimalColl;
		}
		end = clock();
		double timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
		cout << "Computing time : " << timing << " s."<< endl;
		timings.push_back(timing);	
		cout << "**************time " << t << " done****************" << endl;
	} // end of time


   	vector<size_t>::iterator vecIt;
      string fileNameCount = "NumAgg.txt";
		ofstream os1;
      os1.open(fileNameCount.c_str());
      for (vecIt = numAgg.begin(); vecIt < numAgg.end(); vecIt++) {
         os1 << *vecIt << "\n";
      }
      os1 << flush;
      os1.close();
	 	
		vector<double>::iterator It;
      fileNameCount = "TimesDynamic.txt";
      os1.open(fileNameCount.c_str());
	os1 << starttime << "\n";
      for (It = timings.begin(); It < timings.end(); It++) {
         os1 << *It << "\n";
      }
      os1 << flush;
      os1.close();

	return 0;
} // end of dynamic air traffic example program

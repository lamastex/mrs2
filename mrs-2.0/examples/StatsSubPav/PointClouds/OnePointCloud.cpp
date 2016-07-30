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


/*! \file OnePointCloud.cpp
\brief Testing StatsSubPavings (aka SPSnodes) with Bivariate Gaussian / read-in data
 */

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams

#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data
/* what's in dataprep.hpp
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
*/
//for gsl permutations
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include <limits>

/* headers in Gaussian Carver Example
#include "adaptivehistogram.hpp" 
#include "histevalobj.hpp"
#include "piecewise_constant_function.hpp"
*/
#include "carver_seb.hpp"

#include <time.h>   // clock and time classes
#include <iostream>  // input and output streams
#include <fstream>  // file streams
#include<unistd.h>

// to be able to manipulate strings as streams
#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions
#include <iterator>

#include "testDenCommon.hpp" // to use density testing tools
#include "testDenTools.hpp"

#include "mixture_mvn.hpp" // to use MixtureMVN (Jenny's thesis)

using namespace cxsc;
using namespace std;
using namespace subpavings;
using namespace subpavings::kde;

int main(int argc, char ** argv) 
{
    /* setting up dimension and n */
	int reps = atoi(argv[1]); //replications	
    size_t dim = atoi(argv[2]); //dimensions
	size_t n = atoi(argv[3]); //sample size
	int K = atoi(argv[4]); //K-fold CV
	int MaxTempIterations= atoi(argv[5]); //how many temperatures
	double t_lo=atof(argv[6]);//lowest temperature in search
	double t_hi=atof(argv[7]); //highest temperature in search
	
	unsigned int thinout = 100;
	unsigned int samplesNeeded = 100;
	
	cout << "pid = " << getpid() << endl;
	
	cout << "\n" << dim << 
	"-d mixture example to get histogram errors, n = " << n 
	<< " , rep = " << reps << " with " << K << "-fold CV over " 
	<< MaxTempIterations << " temperatures" << endl;
    
    /* setting up filenames etc */
	std::string baseOutputDir("ResultsForBA/");
	
	std::string thisDir;
	{
		ostringstream oss;
		oss << "HistDen" << dim << "D";
		thisDir = oss.str();
	}
	
	/* make the output dir and get back path */
	std::string path = makeDir(baseOutputDir, thisDir);
	
	/* setting up filenames to use for output of log files
     * the program will add dimensions, number of samples, etc, 
     * to these base strings */
    //RPQ with optimal temperature
	string tempFilenameBase = path + "Temp"; // base file name for rep file
    string logFilenameBase = path + "LogMixTwo"; // base file name for log file
    
    //MCMC using LogTemperaturePrior logPrior(t_opt)
    string MCMCFilenameBase = path + "MCMC"; // base file name for rep file
    string logMCMCFilenameBase = path + "LogMixTwoMCMC"; // base file name for log file
    
    //RPQ using LogCatalanPrior
    string logCatFilenameBase = path + "LogMixTwoCat"; // base file name for log file
    
    /* points to use to try to assess approximation accuracy - the errors */
    //GT doesn't really get this part
	size_t intNbase = 1000000;
	size_t intN = (dim > 2 ? dim : 1) * intNbase;

    /* set up a temperature range */
   	std::vector<double> Temperatures;
	double t_Delta=(t_hi-t_lo)/double(MaxTempIterations-1);//change in temperature
	for(int i=0; i<MaxTempIterations; i++) Temperatures.push_back(t_lo+(double(i))*t_Delta);
	cout << "Selection over temperature range: " << Temperatures << endl; 
	
	/* prepare to generate some data from MVN for the tests */
    long unsigned int seed = 1234;
	MixtureMVN* mixMVNptr = makeMixture(dim , seed); 
	// seed gets changed in PriorSelect

	/* If true, will compute the L1-error and KL-information for each
	 * cross validation */
	bool ComputeCVError =false;
	
	/* A K-fold cross-validation will be performed over MaxTempIterations
	 * temperatures to obtain the optimal histogram. 
	 * 3 histograms will be built:
	 * 1. An RPQ histogram based on the optimal temperature.
	 * 2. An RQP histogram based on the LogCatalanPrior (for comparison 
	 *    purposes.
	 * 3. An MCMC histogram based on the optimal temperature.
	 */
	 doPriorSelect(
			dim,
			seed,
			reps,
			n,
			intN,
			K, 
			Temperatures,
			MaxTempIterations,
			tempFilenameBase,	
			MCMCFilenameBase,	
			logFilenameBase,	
			logMCMCFilenameBase,
			logCatFilenameBase,	
			mixMVNptr,
			ComputeCVError,
			thinout,
			samplesNeeded);
	
	bool RPQCat = false;
	if(RPQCat) {
	
	doRPQCat(
		dim,
		seed,
		reps,
		n,
		intN,
		logCatFilenameBase,	
		mixMVNptr);
	}
	
	//delete the data generator 
	delete mixMVNptr;	
	
    return 0;

} // end of bivariate gaussian test program

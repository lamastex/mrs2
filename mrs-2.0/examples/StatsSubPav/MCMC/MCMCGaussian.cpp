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


/*! \file
\brief Testing StatsSubPavings (aka SPSnodes) with MCMC new style (May 2012)
and multivariate (independent dimensions) gaussian data
 */


#include "histall.hpp"  // headers for the histograms
#include "piecewise_constant_function.hpp"  
#include "dataprep.hpp" // headers for getting data

#include <vector>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <unistd.h>

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main(int argc, char* argv[])
{
	// sort out user-defined parameters------------------//
	if ( argc < 5 ) {
		cerr << "Syntax: " << argv[0] << "d n states symmetryIndicator [log_full]" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "d n states, symmetryIndicator");
	}
	
	int d = atoi(argv[1]);  // dimensions
	const int n = atoi(argv[2]);  // number of points to generate
	unsigned int loops = atoi(argv[3]); // number of states to loop through
	int symmetryIndicator = atoi(argv[4]); // 1 for symmetric, 0 otherwise
	
	/* for logging use LOGSTATETRACE to get a log file just of the
	 * number of leaves at each state and the log posteriors for each state
	 * or LOGMCMCTRACE for the traces for each state and the 
	 * sample average, or TXT to get a log file
	 * showing the details of each change of state ie how decision
	 * to split or merge, or not, was taken and also a file showing
	 * the histogram in full at each state */
	LOGGING_LEVEL logging = LOGSTATETRACE;
	
	if ( argc > 5 && atoi( argv[5] ) ) logging = TXT;

		
	cout << argv[0] << " : process id is " << getpid() << std::endl;
    // ------- prepare to generate some data for the tests -----------

	interval pavingIntervalSym(-5,5);
	interval pavingIntervalNonSym(-5.5,5);
	
	interval pavingInterval = pavingIntervalSym; // if we are doing symmetric case
	// but if we've asked for non-symmetric, change this
	if (!symmetryIndicator) pavingInterval = pavingIntervalNonSym;
		
	// set up a random number generator and use mt19937 for generator
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
	long unsigned int seed = 1234;
	gsl_rng_set (r, seed); // change the seed
    
	string samplesFileName; // for samples
    string outputFileName;// for output file
    ofstream oss;         // ofstream object
    oss << scientific;  // set formatting for input to oss
    oss.precision(5);

    double sigma_x=1;   // distribution parameter for Biv Gaussian
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    RVecData* theDataPtr = new RVecData;   // a container for all the points generated

    // make a simulated data set allData to sample from
    for (int i = 0; i < n; i++) {

        rvector thisrv(d);
		{
			double x = 0;
			double y = 0;
		
			gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
								rho, &x, &y);
			thisrv[1] = x;
			thisrv[2] = y;
		}
		for(int j=3; j <= d; j++) {
			double x = 0;
			double y = 0;

			gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
								rho, &x, &y);
			thisrv[j] = y;
        }

		// put points generated into container
        theDataPtr->push_back(thisrv);

    }  // data  should be in theData

    // free the random number generator
    try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow

	/* these parameters are set up here just to get an idea of what
	 * is happening in each state - I don't care about the actual 
	 * sampling. variable loops came from user-supplied parameters.*/
	unsigned int burnin = 10; // don't care about this
	unsigned int thinout = loops/10; // don't care here either.
	
	size_t minPoints = 1;

	// set up proposal distribution object
	UniformProposal proposal;
	// set up prior distribution object
	LogCatalanPrior logPrior;
		
	cout << "\nStart example: n = " << n << " and d = " << d << endl;
	cout << "Paving interval is " << pavingInterval << endl;
		
    {
		// make the paving box
		ivector pavingBox(d);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
		// make an Adaptive Histogram with the given pavingBox and, by default,
		// holdAllStats = false so that the underlying rootPaving managed by the
		// myHistFirst will not maintain all available stats, only counts
		AdaptiveHistogram myHistFirst(pavingBox);

		// put in the data in a 'pulse' with no splitting, ie into root box
		    bool successfulInsertion = myHistFirst.insertFromRVec(*theDataPtr);

		if (successfulInsertion) {

			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			#if(1) // samples only
				std::vector < PiecewiseConstantFunction > samples;
				samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging);
			#endif
			
			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

			#if(0)
				cout << "Finished MCMC sample averaging" << endl;
				string mcmcAvFilename = "pcfAverage.txt";
				
				mcmcAv.outputToTxtTabs(mcmcAvFilename);
			#endif	
			
			
		}
		else cout << "Failed to insert data" << endl;
	}
	
	{
		// make the paving box
		ivector pavingBox(d);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
		// make an Adaptive Histogram with the given pavingBox and, by default,
		// holdAllStats = false so that the underlying rootPaving managed by the
		// myHistFirst will not maintain all available stats, only counts
		AdaptiveHistogram myHistFirst(pavingBox);

		// put in the data in a 'pulse' with no splitting, ie into root box
		bool successfulInsertion = myHistFirst.insertFromRVec(*theDataPtr);

		if (successfulInsertion) {
			
			
			// pq down to 1024 leaves
			LOGGING_LEVEL logPQ = NOLOG;
			CompCount compCount;
			double maxLeafPQ = std::pow(2.0, d);
			CritLeaves_GTE critStop(static_cast<size_t>(maxLeafPQ));
			cout << "PQ to "  << (static_cast<size_t>(maxLeafPQ)) << " leaves" << endl;
			
			bool successfulHist = myHistFirst.prioritySplit(compCount,
					 critStop, logPQ, minPoints);
					

			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			#if(1) // samples only
				std::vector < PiecewiseConstantFunction > samples;
				samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging);
			#endif
			
			
			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

			#if(0)
				cout << "Finished MCMC sample averaging" << endl;
				string mcmcAvFilename = "pcfAverage.txt";
				
				mcmcAv.outputToTxtTabs(mcmcAvFilename);
			#endif	
			
			
		}
		else cout << "Failed to insert data" << endl;
	}

	delete theDataPtr;

    return 0;

} // end of MCMC test program



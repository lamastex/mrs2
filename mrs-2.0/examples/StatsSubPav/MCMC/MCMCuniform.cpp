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
and uniform data
 */


#include "histall.hpp"  // headers for the histograms
#include "piecewise_constant_function.hpp"  
#include "dataprep.hpp" // headers for getting data

#include <vector>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main()
{
    // ------- prepare to generate some data for the tests -----------

// set up a random number generator
    const gsl_rng_type * T;
    gsl_rng * r;

    const int n=100000;    // number to generate
    //create a generator chosen by the environment variable GSL_RNG_TYPE

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    string samplesFileName; // for samples
    string outputFileName;// for output file
    ofstream oss;         // ofstream object
    oss << scientific;  // set formatting for input to oss
    oss.precision(5);

    int d = 5; // dimension of the box to sample data from
    ivector pavingBox(d);
    interval pavingInterval(0,1);
    for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

    RVecData theData;   // a container for all the points generated

    // make a simulated data set allData to sample from
    for (int i = 0; i < n; i++) {

        rvector thisrv(d);
        for(int j=1; j <= d; j++) {
			thisrv[j]  = gsl_rng_uniform(r);
        }

        // put points generated into container
        theData.push_back(thisrv);

    }  // data  should be in theData

    // free the random number generator
    gsl_rng_free (r);

/*
    // output data to a file
    string dataFile = "dataFile.txt";
    ofstream osd(dataFile.c_str());         // replace data
    if (osd.is_open()) {
        RVecDataItr dit;
        for (dit = theData.begin(); dit < theData.end(); dit++) {
            //osd << *dit << "\n";
            osd << (*dit)[1] << "\t" << (*dit)[2] << "\n";
        }
        osd.close();
    }
    else {
        std::cout << "Error: could not open file named "
            << dataFile << std::endl << std::endl;
    }
*/
	unsigned int loops = 1000;
	unsigned int burnin = 10;
	unsigned int thinout = 1;
	LOGGING_LEVEL logging = NOLOG;
	int minPoints = 1;

	// set up proposal distribution object
	UniformProposal proposal;
	// set up prior distribution object
	LogCatalanPrior logPrior;
		
	cout << endl << endl;
    cout << "\n\nStart example: n = " << n << " and d = " << d << endl;

    {
		// make an Adaptive Histogram with the given pavingBox and, by default,
		// holdAllStats = false so that the underlying rootPaving managed by the
		// myHistFirst will not maintain all available stats, only counts
		AdaptiveHistogram myHistFirst(pavingBox);

		// put in the data in a 'pulse' with no splitting, ie into root box
		    bool successfulInsertion = myHistFirst.insertFromRVec(theData);

		if (successfulInsertion) {

			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			// MCMC with 2000 states, burn in 1500, thinout every 100 etc
			// create a log file only for samples (no dot graphs)
			
			
				std::vector < PiecewiseConstantFunction > samples;
				samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging);
			#if(0)
			PiecewiseConstantFunction mcmcAv = myHistFirst.MCMC(
													loops, burnin, thinout,
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
		// make an Adaptive Histogram with the given pavingBox and, by default,
		// holdAllStats = false so that the underlying rootPaving managed by the
		// myHistFirst will not maintain all available stats, only counts
		AdaptiveHistogram myHistFirst(pavingBox);

		// put in the data in a 'pulse' with no splitting, ie into root box
		    bool successfulInsertion = myHistFirst.insertFromRVec(theData);

		if (successfulInsertion) {

			cout << "Starting MCMC AdaptiveHistogram samples "  << endl;
			clock_t start, end;
			start = clock();

			// MCMC with 2000 states, burn in 1500, thinout every 100 etc
			// create a log file only for samples (no dot graphs)
			
			std::vector < PiecewiseConstantFunction > samples;
			samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
												proposal, logPrior,
												minPoints, logging);
			
			#if(0)
				AdaptiveHistogramCollator tempColl;
				tempColl.addToCollation(samples);
				AdaptiveHistogramCollator mcmcAv = tempColl.makeAverage();
			#endif

			
			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

			#if(0)
				cout << "Finished MCMC sample averaging" << endl;
				string mcmcAvFilename = "adhAverage.txt";
				
				mcmcAv.outputToTxtTabs(mcmcAvFilename);
			#endif	
			
			
		}
		else cout << "Failed to insert data" << endl;
	}
	
    return 0;

} // end of MCMC test program



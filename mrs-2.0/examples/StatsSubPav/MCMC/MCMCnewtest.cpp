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
#include <cassert>  // input and output streams

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main()
{
    // ------- prepare to generate some data for the tests -----------

// set up a random number generator
    const gsl_rng_type * T;
    gsl_rng * r;

    const int n=1000;    // number to generate
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

	int prec = 10;
	bool confirm = true;
	
	unsigned int loops = 100;
	unsigned int burnin = 10;
	unsigned int thinout = 10;
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
			
			/* make a generator for this mcmc run */
			int seed = 1234;
			gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937);
				
			gsl_rng_set(rgsl, seed);

			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			// MCMC with 2000 states, burn in 1500, thinout every 100 etc
			// create a log file only for samples (no dot graphs)
			
			
			std::vector < PiecewiseConstantFunction > samples;
			samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
												proposal, logPrior,
												minPoints, logging,
												rgsl);

			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
			
			assert(!samples.empty());
			std::vector < PiecewiseConstantFunction >::iterator it = samples.begin();
			it++;
			for (; it < samples.end(); ++it) {
				samples.front() += (*it);
			}
			PiecewiseConstantFunction result = samples.front()/cxsc::real(1.0*samples.size());
			
			cout << "Finished MCMC averaging over samples" << endl;
			string mcmcAvFilename = "pcfAverage.txt";
			
			result.outputToTxtTabs(mcmcAvFilename, prec, confirm);
			
			
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
			
			/* make a generator for this mcmc run */
			int seed = 1234;
			gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937);
				
			gsl_rng_set(rgsl, seed);


			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			PiecewiseConstantFunction mcmcAv = myHistFirst.MCMC(
													loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging,
													rgsl);
			
			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
			
			cout << "Finished MCMC get average" << endl;
			string mcmcAvFilename = "pcfDirectAverage.txt";
			
			mcmcAv.outputToTxtTabs(mcmcAvFilename, prec, confirm);
					
			
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
			
			
			size_t nSamples = samples.size();
			assert (nSamples > 1);
			
			PiecewiseConstantFunction mcmcAv = samples[0];
			for (size_t i = 1; i < nSamples; ++i) {
				mcmcAv += samples[i];
			}
			mcmcAv /= static_cast<int>(nSamples);
			
			end = clock();

			cout << "Computing time : "
			   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

			cout << "Finished MCMC sample averaging" << endl;
			string mcmcAvFilename = "pcfSampleAverage.txt";
				
			mcmcAv.outputToTxtTabs(mcmcAvFilename, prec, confirm);
						
			
		}
		else cout << "Failed to insert data" << endl;
	}
	
    return 0;

} // end of MCMC test program



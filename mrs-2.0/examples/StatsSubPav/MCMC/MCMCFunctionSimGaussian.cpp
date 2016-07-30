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
and Gaussian data generated from a function estimate
 */


#include "histall.hpp"  // headers for the histograms
#include "intervalmappedspnode_measurers.hpp" // ordering for pq split
#include "functionestimator_interval.hpp"
#include "piecewise_constant_function.hpp"  

#include "GaussianFobj.hpp" // fobj

#include <vector>
#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <unistd.h> // required for getpid() call

using namespace cxsc;
using namespace std;
using namespace subpavings;

int main(int argc, char* argv[])
{
	// sort out user-defined parameters------------------//
	if ( argc < 6 ) {
		cerr << "Syntax: " << argv[0] << "d maxLeavesEst n states symmetryIndicator [log_full]" << endl;
		throw std::runtime_error("Syntax: " + std::string(argv[0]) + "d n states, symmetryIndicator");
	}
	
	int d = atoi(argv[1]);  // dimensions
	size_t maxLeavesEst = atoi(argv[2]);  // number of leaves in estimator
	const int n = atoi(argv[3]);  // number of points to generate
	unsigned int loops = atoi(argv[4]); // number of states to loop through
	int symmetryIndicator = atoi(argv[5]); // 1 for symmetric, 0 otherwise
	
	/* for logging use LOGSTATETRACE to get a log file just of the
	 * number of leaves at each state and the log posteriors for each state
	 * or LOGMCMCTRACE for the traces for each state and the 
	 * sample average, or TXT to get a log file
	 * showing the details of each change of state ie how decision
	 * to split or merge, or not, was taken and also a file showing
	 * the histogram in full at each state */
	LOGGING_LEVEL logging = LOGSTATETRACE;
	
	if ( argc > 6 && atoi( argv[6] ) ) logging = TXT;

		
	cout << argv[0] << " : process id is " << getpid() << std::endl;
	
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
    
	interval pavingIntervalSym(-5,5);
	interval pavingIntervalNonSym(-5.5,5);
	
	interval pavingInterval = pavingIntervalSym; // if we are doing symmetric case
	// but if we've asked for non-symmetric, change this
	if (!symmetryIndicator) pavingInterval = pavingIntervalNonSym;
	
	//-------------- make estimate ------------//
		
	cout << "\nMake the function estimator to " << maxLeavesEst << " leaves" << endl;
	
	// specify function object (from /examples/MappedTargets
	GaussianFobj fobj;
	
	/* function estimate is going to use same box as the histograms */
	ivector pavingBoxEst(d);
	for(int k=1; k <= d; k++) pavingBoxEst[k] = pavingInterval;
	
	// Use fobj and pavingBox to get an estimator
	FunctionEstimatorInterval estimator(pavingBoxEst, fobj);
	
	LOGGING_LEVEL logEst = NOLOG; // logging for making estimator
	
	#if(1)
	size_t maxLeavesEstDown = static_cast<size_t>(1.2*maxLeavesEst); // go down to 1.2 x max
	#endif
	#if(0)
		size_t maxLeavesEstDown = maxLeavesEst;
	#endif
	
	cout << "pq down to max leaves " << maxLeavesEstDown << endl;
	
	clock_t startEst = clock();
	
	// priority split driven by splitting leaf with max reimann diff
	ReimannDiffMeasurer measurer;
	estimator.prioritySplit(measurer, maxLeavesEstDown, logEst);
			
	// stop recording time here
	clock_t endEst = clock();
	cout << "Number of leaves in estimate: " << estimator.getRootLeaves() << " s."<< endl;	
	cout << "After split, getTotalAreaOfIntervalBand() = "
		<< estimator.getTotalAreaOfIntervalBand() << endl;
	double timingEst1 = ((static_cast<double>(endEst - startEst)) / CLOCKS_PER_SEC);
	cout << "Computing time for pq split in estimate: " << timingEst1 << " s."<< endl;
	startEst = clock();
	
	#if(1) 
		cout << "Hull propagation" << endl;
		estimator.hullPropagation();
		
		cout << "Priority merge to " << maxLeavesEst << " leaves" << endl;
		#if(0)
		// priority merge driven by minimising increase the reimann diff
		estimator.priorityMergeOnLoss(maxLeavesEst, logEst);
		#endif
		#if(1)
		// priority merge driven by merging cherry with minimum reimann diff
		estimator.priorityMerge(maxLeavesEst, logEst);
		#endif
					
		// stop recording time here
		endEst = clock();	
		double timingEst2 = ((static_cast<double>(endEst - startEst)) / CLOCKS_PER_SEC);
		cout << "Computing time for hull propagate and merge up in estimate: " << timingEst2 << " s."<< endl;
		
		cout << "After propagation and priority merge, getTotalAreaOfIntervalBand() = " 
					<< estimator.getTotalAreaOfIntervalBand() << endl;
		cout << "number of leaves is = " << estimator.getRootLeaves() << endl;
	#endif
	
	cout << "Making estimate and normalising" << endl;
	// Make PiecewiseConstantFunction estimate from estimator
	
	PiecewiseConstantFunction estimate = estimator.makePiecewiseConstantFunction();
	estimate.normalise();
	
	cout << "estimate has integral " << estimate.getTotalIntegral() << endl;
	
	// Use PiecewiseConstantFunction to generate data, supplying our own rng
		
	cout << "\nGenerating data for simulation" << endl;
	
	// set up a random number generator
	gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937);
	unsigned long int seed = 1234;
	gsl_rng_set(r, seed); // seed it
	
	RVecData* theDataPtr = new RVecData;   // a container for all the points generated

	clock_t startData = clock();
	
	estimate.simulateData(*theDataPtr, n, r);
	
	// stop recording time here
	clock_t endData = clock();	
	double timingData = ((static_cast<double>(endData - startData)) / CLOCKS_PER_SEC);
	cout << "Computing time for simulating data: " << timingData << " s."<< endl;

	cout << (*theDataPtr).size() << " points generated" << endl;
	
	try {
		gsl_rng_free (r);
		r = NULL;
	}
	catch(...) {}// catch and swallow

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

			#if(0) // samples only
				std::vector < PiecewiseConstantFunction > samples;
				samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging);
			#endif
			#if(1) // get an average so we can do IAE against estimate
				burnin = loops/2; // need to be slightly more careful	
				thinout = (loops - burnin)/10; // and here - should get 10 samples unless we get thinout = 0 which will cause an exception.
				PiecewiseConstantFunction mcmcAv = myHistFirst.MCMC(
														loops, burnin, thinout,
														proposal, logPrior,
														minPoints, logging);
				// now get the IAE of mcmcAv against the estimate
				real thisIAE = mcmcAv.getIAE(estimate);
				
				cout << "IAE against the estimate is " << thisIAE << endl;
			
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
			
			
			// make the histogram have the same shape as the estimate
			// except if minPoints will not allow this
			myHistFirst.reshapeToUnion(estimate, minPoints);
			cout << "After reshaping to shape of function, histogram has " << myHistFirst.getRootLeaves() << " leaves" << endl;
			
			cout << "Starting MCMC PiecewiseConstantFunction samples "  << endl;
			clock_t start, end;
			start = clock();

			
			#if(0) // samples only
				std::vector < PiecewiseConstantFunction > samples;
				samples = myHistFirst.MCMCsamples(samples, loops, burnin, thinout,
													proposal, logPrior,
													minPoints, logging);
			#endif
			#if(1) // get an average so we can do IAE against estimate
				burnin = loops/2; // need to be slightly more careful	
				thinout = (loops - burnin)/10; // and here - should get 10 samples unless we get thinout = 0 which will cause an exception.
				PiecewiseConstantFunction mcmcAv = myHistFirst.MCMC(
														loops, burnin, thinout,
														proposal, logPrior,
														minPoints, logging);
				// now get the IAE of mcmcAv against the estimate
				real thisIAE = mcmcAv.getIAE(estimate);
				
				cout << "IAE against the estimate is " << thisIAE << endl;
			
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



/*
* Copyright (C) 2014 Jennifer Harlow
*
*/


/*! \file
\brief Common routines for KDE tests

 */

#include "testDenCommon.hpp"
#include "testDenTools.hpp"

#include "mixture_mvn.hpp"

/*
#include "MCMCGRAutoNew.hpp"
#include "MCMCGRDiagnosticPSRFLeaves.hpp"
#include "MCMCGRDiagnosticPSRFLogpost.hpp"
#include "MCMCGRDiagnosticPSRFCherries.hpp"
#include "MCMCGRDiagnosticPSRFAvDepth.hpp"
#include "MCMCGRDiagnosticIntervalLeaves.hpp"
#include "MCMCGRDiagnosticIntervalCherries.hpp"
#include "MCMCGRDiagnosticIntervalAvDepth.hpp"

#include "real_kde_mid_estimate.hpp"
#include "functionestimator_kde.hpp"
*/

#include "adaptivehistogram.hpp" 
#include "histevalobj.hpp"
//#include "carver_seb.hpp"


//for gsl permutations
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>


#include <limits> // to use inifinity
#include <vector>
#include <string>
#include <ctime>   // clock and time classes
#include <iostream>  // input and output streams
#include <fstream>  // file streams
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator>  // output iterator
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions


using namespace cxsc;
using namespace subpavings;
using namespace std;

/*
PiecewiseConstantFunction _doHistMCMC(
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					LogMCMCPrior& logPrior);

					
std::vector< subpavings::AdaptiveHistogram* >& _getHistMCMCStarts(
					std::vector< subpavings::AdaptiveHistogram* >& hists,
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, 
					UniformSSMProposal proposal,
					size_t minPoints,	
					const string& histFilenameBase,
					long unsigned int seed,
					LogMCMCPrior& logPrior);
*/					
std::string makeFilename(const std::string& filenameBase,
						const std::string& insert,
						const std::string& suffix,
						size_t n,
						int rep);
						
std::string makeFilename(const std::string& filenameBase,
						const std::string& insert,
						const std::string& suffix,
						size_t n,
						int rep)
{
	/* Search for the last '/' in the scalarsFilenName, break
	 * it there, join again with scalarType and runIDstr in the middle*/

	std::string retValue;
	
	std::string file = filenameBase; 
	std::string path("");
	size_t found = file.find_last_of("/");
	if (found!=string::npos) {
		path = file.substr(0,found+1);
		file = file.substr(found+1);
	}
	// else file stays as full base
		
	ostringstream oss;
	oss << path << insert << file << "_n" << n << "_r" << (rep+1) << suffix;
		
	retValue  = oss.str();
		
	return retValue;
} 


/* overloaded version:
 * - makes its own root box
 * - use default thinout and samples */
/*
PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					LogMCMCPrior& logPrior)
		
{

	AdaptiveHistogram adhMCMC; // let it make its own root box
	bool successfulInsertion = 
				adhMCMC.insertRvectorsFromVectorOfVecDbls(simdata);
	if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
	
	unsigned int thinout = 100;
	unsigned int samplesNeeded = 100;
	

	return _doHistMCMC(adhMCMC,
				rep, n, minPoints,	
				histFilenameBase,
				timingMake,
				seed,
				thinout,
				samplesNeeded,
				logPrior);

}
*/

/* overloaded version:
 * - makes its own root box 
 * - use given thinout and samples*/
/*
PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					LogMCMCPrior& logPrior)		
{

	AdaptiveHistogram adhMCMC; // let it make its own root box
	bool successfulInsertion = 
				adhMCMC.insertRvectorsFromVectorOfVecDbls(simdata);
	if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
	
	
	return _doHistMCMC(adhMCMC,
				rep, n, minPoints,	
				histFilenameBase,
				timingMake,
				seed,
				thinout,
				samplesNeeded,
				logPrior);

}
*/

/* overloaded version:
 * - uses given root box
 * - use default thinout and samples */
/*
PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					LogMCMCPrior& logPrior)
		
{

	AdaptiveHistogram adhMCMC(box);
	bool successfulInsertion = 
				adhMCMC.insertRvectorsFromVectorOfVecDbls(simdata);
	if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
	
	unsigned int thinout = 100;
	unsigned int samplesNeeded = 100;
	

	return _doHistMCMC(adhMCMC,
				rep, n, minPoints,	
				histFilenameBase,
				timingMake,
				seed,
				thinout,
				samplesNeeded,
				logPrior);

}
*/

/* overloaded version:
 * - uses given root box
 * - use given thinout and samples */
/*
PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					LogMCMCPrior& logPrior)
		
{

	AdaptiveHistogram adhMCMC(box);
	bool successfulInsertion = 
				adhMCMC.insertRvectorsFromVectorOfVecDbls(simdata);
	if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
	
	

	return _doHistMCMC(adhMCMC,
				rep, n, minPoints,	
				histFilenameBase,
				timingMake,
				seed,
				thinout,
				samplesNeeded,
				logPrior);

}
*/
					
/* the actual code that does the MCMC histogram */
/*
PiecewiseConstantFunction _doHistMCMC(
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					LogMCMCPrior& logPrior)
{
	clock_t repStarttime = clock();
		
	// set up proposal distribution object - the stay split merge chain
	double probStay = 1.0E-6;
	UniformSSMProposal proposal(probStay);
	
	
	// a container for our histograms
	std::vector< subpavings::AdaptiveHistogram* > hists;
	
	double timingStarts = 0.0;
	{
		clock_t starttime = clock();
		

		// get the starts
		_getHistMCMCStarts(
						hists,
						adhMCMC,
						rep, n, 
						proposal,
						minPoints,	
						histFilenameBase,
						seed,
						logPrior);

		clock_t endtime = clock();	
		timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
		cout << "time to get starts = " << timingStarts << endl;
	}
			
	// -------------------- do the MCMC ---------------------

	/* MCMC */
	/*
	unsigned int maxLoops = 5000000;
	size_t d = adhMCMC.getDimensions();
	if (d > 3) maxLoops *= 2;
	
	real tolerance = 0.1;
	size_t samplingInterval = 100;
	double percent = 0.80;
	
	cout << "\nDoing autoMCMC: maxLoops = " << maxLoops << ", thinout = " << thinout
				<< ", samplesNeeded = " << samplesNeeded << endl;
	cout << " Testing convergence using all three (leaves, cherries, average leaf depth) (Interval)\n" << endl;
	
	std::vector< MCMCGRAuto::Diagnostic * > diagPtrs;
	//diagnostics have got to be in scope for the whole time!
	MCMCGRDiagnosticIntervalLeaves lDiag(tolerance, samplingInterval, percent);
	diagPtrs.push_back(&lDiag);


	MCMCGRDiagnosticIntervalCherries cDiag(tolerance, samplingInterval, percent);
	diagPtrs.push_back(&cDiag);


	MCMCGRDiagnosticIntervalAvDepth aDiag(tolerance, samplingInterval, percent);
	diagPtrs.push_back(&aDiag);
	
	unsigned long int seedMCMC = seed+rep+2;
		
	MCMCGRAutoNew autoMCMCNew(diagPtrs, seedMCMC);
	
	std::string scalarsFileName("Log.txt");
	clock_t startMCMC = clock();
	
	size_t chains = hists.size();
	cout << "chains: " << chains << endl; 
	
	//pass logPrior into this function
	//figure out how to switch on the LOGGING_LEVEL
	PiecewiseConstantFunction samplesAvg = 
	autoMCMCNew.doMCMCGRAuto(hists, maxLoops, 
						samplesNeeded, thinout, 
						minPoints,
						diagPtrs.size(), //flagthreshold
						rep+1,
						scalarsFileName,
						logPrior);
	
	clock_t endMCMC = clock();	
	
	double timingMCMC = (static_cast<double>(endMCMC-startMCMC)/CLOCKS_PER_SEC);	
	cout << "time to do MCMC = " << timingMCMC << endl;
	
	cout << "\nTotal time making histogram = " << (timingMCMC + timingStarts) << endl;
	timingMake.push_back(timingMCMC + timingStarts);
	
	
	for (size_t i = 0; i < hists.size(); ++i) {
		if (NULL != hists[i]) delete hists[i];
		hists[i] = NULL;
	}	

	
	cout << "\nEnd of making histogram with auto mcmc\n" << endl;
	
	return samplesAvg;
}
*/

/*
std::vector< subpavings::AdaptiveHistogram* >& _getHistMCMCStarts(
					std::vector< subpavings::AdaptiveHistogram* >& hists,
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, 
					UniformSSMProposal proposal,
					size_t minPoints,	
					const string& histFilenameBase,
					long unsigned int seed,
					LogMCMCPrior& logPrior)
{
	
	std::string postFileName("");
	std::string checkPostFileNameBase("");
		
	int precPQ = 5;
	
	/* some guesses for max points in a node to stop posterior queue */
	/*
	size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));
		
	/* some guesses for maximum leaves we'll let SEB queue go to */
	/*
	size_t maxLeavesSEB = (adhMCMC.getDimensions()) * n / critSEB; // integer division
	if (maxLeavesSEB > n/2) maxLeavesSEB = n/2;
	size_t maxLeavesCarving = (maxLeavesSEB*0.9)/8; // allow for increase 2^3
	if (n >= 5000000) maxLeavesCarving*= 2;
	if (n >= 10000000) maxLeavesCarving*= 4;
	
	SPSNodeMeasureVolMassMinus compCarving(adhMCMC.getRootCounter());

	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorCarving( compCarving, maxLeavesCarving);
	
	SPSNodeMeasureCount compSEB;
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorSEB( compSEB, critSEB, maxLeavesSEB);


	// behaviour of this has changed from first runs
	bool stopOnMaxPosterior = true;
	unsigned long int seedStarts = seed+rep+1;
	
	{				
		cout << "\n\nStarts chosen using\n\tmaxLeavesCarving\t" << maxLeavesCarving;
		cout << "\n\tcritSEB\t" << critSEB;
		cout << "\n\tmaxLeavesSEB\t" << maxLeavesSEB << endl;
		
	}

	int keep = 3;
	int carvingStarts = 10;
		
		CarverSEB::findStartingPointsOverdispersed(
			adhMCMC,
			hists,
			evaluatorCarving, 
			evaluatorSEB, 
			logPrior,
			minPoints,
			carvingStarts,
			keep,
			stopOnMaxPosterior, 
			postFileName,
			checkPostFileNameBase,	
			precPQ,
			seedStarts);


	assert(hists.size()==keep);
	
	return hists;
}
*/

std::vector < real >& getPCFDensities(const PiecewiseConstantFunction& pcf,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities)
{
	std::vector < double > timing;
	return getPCFDensities(pcf,
						intPts,
						pcfDensities,
						timing);
}

std::vector < real >& getPCFDensities(const PiecewiseConstantFunction& pcf,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities,
						std::vector < double >& timing)
{
		
	size_t N = intPts.size();
	
	std::vector < real > tmp(N);

	clock_t start = clock();
	
	if(N) {
		
		size_t dim = intPts.front().size(); // assume all the same
				
		for (size_t i = 0; i < N; ++i) { 
		
			rvector rv(dim);
			for (size_t j = 0; j < dim; ++j) rv[j+1] = intPts[i][j];
					
			tmp[i] = pcf.pointwiseExtension(rv);
			
		}
	}
	clock_t end = clock();
	timing.push_back(static_cast<double>(end-start)/CLOCKS_PER_SEC);
		
	pcfDensities.swap(tmp);
	
	return pcfDensities;
	
}


// censor the intPts data according to whether it has positive density in the pcf
std::vector < real >& getPCFDensitiesCensor(const PiecewiseConstantFunction& pcf,
						std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities)
{
	std::vector < double > timing;
	
	return getPCFDensitiesCensor(pcf,
								intPts,
								pcfDensities,
								timing);
}



// censor the intPts data according to whether it has positive density in the pcf
std::vector < real >& getPCFDensitiesCensor(const PiecewiseConstantFunction& pcf,
						std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities,
						std::vector < double >& timing)
{
	
	size_t N = intPts.size();
	
	std::vector < real > tmp;
	tmp.reserve(N);
	std::vector < std::vector < real > > tmpIntPts;
	tmpIntPts.reserve(N);

	clock_t start = clock();
	
	if(N) {
		
		size_t dim = intPts.front().size(); // assume all the same
				
		for (size_t i = 0; i < N; ++i) { 
		
			rvector rv(dim);
			for (size_t j = 0; j < dim; ++j){ rv[j+1] = intPts[i][j]; }
					
			real r = pcf.pointwiseExtension(rv);
			if (r > cxsc::MinReal) {
				tmpIntPts.push_back(intPts[i]);
				tmp.push_back(r);
			}
			
		}
	}
	clock_t end = clock();
	timing.push_back(static_cast<double>(end-start)/CLOCKS_PER_SEC);
		
	pcfDensities.swap(tmp);
	intPts.swap(tmpIntPts);
	
	cout << "\nNumber of importance sample points censored from " 
					<< N << " to " << intPts.size() << "\n" << endl;
					
	return pcfDensities;
	
}		


std::vector < real >& getPCFDensitiesCensor(const PiecewiseConstantFunction& pcf,
						subpavings::RVecData& intPts,
						std::vector < real >& pcfDensities, int dim)
{
		
	size_t N = intPts.size();
	
	std::vector < real > tmp;
	tmp.reserve(N);
	RVecData tmpIntPts;
	//std::vector < std::vector < real > > tmpIntPts;
	tmpIntPts.reserve(N);
	
	if(N) {
		
		//size_t dim = intPts.front().size(); // assume all the same
				
		for (size_t i = 0; i < N; ++i) { 
		
			rvector rv(dim);
			for (size_t j = 1; j <= dim; ++j){ rv[j] = intPts[i][j]; }
					
			real r = pcf.pointwiseExtension(rv);
			if (r > cxsc::MinReal) {
				tmpIntPts.push_back(intPts[i]);
				tmp.push_back(r);
			}
			
		}
	}
	
	pcfDensities.swap(tmp);
	intPts.swap(tmpIntPts);
	
	cout << "\nNumber of importance sample points censored from " 
					<< N << " to " << intPts.size() << "\n" << endl;
					
	return pcfDensities;
	
}	

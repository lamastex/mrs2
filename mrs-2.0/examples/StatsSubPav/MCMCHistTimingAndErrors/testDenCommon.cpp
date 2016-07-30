/*
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


#include "adaptivehistogram.hpp" 
#include "histevalobj.hpp"
#include "carver_seb.hpp"



#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

     


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
using namespace subpavings::kde;
using namespace std;


PiecewiseConstantFunction _doHistMCMC(
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded);

					
std::vector< subpavings::AdaptiveHistogram* >& _getHistMCMCStarts(
					std::vector< subpavings::AdaptiveHistogram* >& hists,
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, 
					UniformSSMProposal proposal,
					LogCatalanPrior logPrior,
					size_t minPoints,	
					const string& histFilenameBase,
					long unsigned int seed);
					
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
PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed)
		
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
				samplesNeeded);

}

/* overloaded version:
 * - makes its own root box 
 * - use given thinout and samples*/
PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded)		
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
				samplesNeeded);

}

/* overloaded version:
 * - uses given root box
 * - use default thinout and samples */
PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed)
		
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
				samplesNeeded);

}

/* overloaded version:
 * - uses given root box
 * - use given thinout and samples */
PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded)
		
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
				samplesNeeded);

}
						
/* the actual code that does the MCMC histogram */
PiecewiseConstantFunction _doHistMCMC(
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded)
{
	clock_t repStarttime = clock();
		
	// set up proposal distribution object - the stay split merge chain
	double probStay = 1.0E-6;
	UniformSSMProposal proposal(probStay);
	// set up prior distribution object
	LogCatalanPrior logPrior;

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
						logPrior,
						minPoints,	
						histFilenameBase,
						seed);

		clock_t endtime = clock();	
		timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
		cout << "time to get starts = " << timingStarts << endl;
	}
			
	// -------------------- do the MCMC ---------------------

	/* MCMC */
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
	
	std::string scalarsFileName("");
	clock_t startMCMC = clock();
		
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



std::vector< subpavings::AdaptiveHistogram* >& _getHistMCMCStarts(
					std::vector< subpavings::AdaptiveHistogram* >& hists,
					AdaptiveHistogram& adhMCMC,
					int rep, size_t n, 
					UniformSSMProposal proposal,
					LogCatalanPrior logPrior,
					size_t minPoints,	
					const string& histFilenameBase,
					long unsigned int seed)
{
	
	std::string postFileName("");
	std::string checkPostFileNameBase("");
		
	int precPQ = 5;
	
	/* some guesses for max points in a node to stop posterior queue */
	size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));
		
	/* some guesses for maximum leaves we'll let SEB queue go to */
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
			for (size_t j = 0; j < dim; ++j) rv[j+1] = intPts[i][j];
					
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


/* overloaded version to supply empty slice dims and slice pts */
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		MixtureMVN* mixMVNptr,
		unsigned int thinout,
		unsigned int samplesNeeded)
{
	std::vector < std::vector < int > > vecSliceDims;
	std::vector < std::vector < double > > vecSlicePts;
	
	doDenEst(
			dim,
			seed,
			reps,
			minPoints,
			n,
			intN,
			histFilenameBase,
			logFilenameBase,	
			mixMVNptr,
			vecSliceDims,
			vecSlicePts,
			thinout,
			samplesNeeded);
		
		
}

/* overloaded version to supply samples and thinout */
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		MixtureMVN* mixMVNptr,
		const std::vector < std::vector < int > >& vecSliceDims,
		const std::vector < std::vector < double > >& vecSlicePts)
{
	
	unsigned int thinout = 100;
	unsigned int samplesNeeded = 100;
	
	doDenEst(
			dim,
			seed,
			reps,
			minPoints,
			n,
			intN,
			histFilenameBase,
			logFilenameBase,	
			mixMVNptr,
			vecSliceDims,
			vecSlicePts,
			thinout,
			samplesNeeded);
	
}

/* code to do density estimation */
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		MixtureMVN* mixMVNptr,
		const std::vector < std::vector < int > >& vecSliceDims,
		const std::vector < std::vector < double > >& vecSlicePts,
		unsigned int thinout,
		unsigned int samplesNeeded)
{
	
	//only used if dim > 1
	vector < int > margdims1(1,1);
	vector < int > margdims2(1,2);
	
	/* coverage regions 
     * (only used if we are slicing) */
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	
			
	/* containers for stuff we will be storing */
	std::vector < std::vector < double > > timingMake(reps);
	std::vector < std::vector < double > > timingDensities(reps);
	std::vector < std::vector < double > > timingIntDensities(reps);
	std::vector < std::vector < real > > avLogDens(reps);
	std::vector < std::vector < real > > avLogDenRatios(reps);
	std::vector < std::vector < real > > estL1ErrorsQR(reps);
	std::vector < std::vector < size_t > > leaves(reps);
	
	size_t intNcensored = intN;
					
			
	for (int rep = 0; rep < reps; ++rep) { 
		
		cout << "\nrep number " << (rep+1) << endl;
		
		std::string endFilename;
		{
			ostringstream oss;
			oss << "_n" << n << "_r" << (rep+1) << ".txt";
			endFilename = oss.str();
		}
	
        /* a file for each repetition */
		string repLogFilename = histFilenameBase + "RepLog" + endFilename;
	
		long unsigned int dataseed = seed+rep;
		mixMVNptr->resetPRNG(dataseed);
	
		std::vector < std::vector < double > > simdata;
		
		cout << "\nGenerate " << n << " random values:" << endl;
		mixMVNptr->prn(simdata, n);
		
		cout << simdata[0] << endl; getchar();
		
		cout << "\nGet pcf from mcmc" << endl;
		
		PiecewiseConstantFunction mcmcPCF = doHistMCMC(rep, n, minPoints,	
					histFilenameBase,
					simdata,
					timingMake[rep],
					seed,
					thinout,
					samplesNeeded);
		
		ivector box = mcmcPCF.getRootBox();
		real boxVol = realVolume(box);
		
		cout << "\npcf box volume is " << boxVol << endl;
		
		
		
		{
			size_t l = mcmcPCF.getRootLeaves();
			cout << "\nMCMC average with " << l << " leaves"<< endl;
			
			leaves[rep].push_back(l);
					
			{
				ostringstream oss;
				oss << histFilenameBase << "_n" << n << "_l" <<l << "_r" << (rep+1) << ".txt";
				string filename = oss.str();
				mcmcPCF.outputToTxtTabs(filename);
				cout << "Average output to " << filename << endl;
			}
            
            /* if dim > 1, do marginals on first two dimensions and slices */
			if (dim > 1) {
				{
					PiecewiseConstantFunction tmpm = mcmcPCF.makeMarginal(margdims1);
					ostringstream oss;
					oss << histFilenameBase << "_marg_d" << margdims1[0] << "_n" << n << "_l" << l << "_r" << (rep+1) << ".txt";
					string filename = oss.str();
					tmpm.outputToTxtTabs(filename);
					cout << "Marginal 1 output to " << filename << endl;
				}
				{
					PiecewiseConstantFunction tmpm = mcmcPCF.makeMarginal(margdims2);
					ostringstream oss;
					oss << histFilenameBase << "_marg_d" << margdims2[0] << "_n" << n << "_l" << l << "_r" << (rep+1) << ".txt";
					string filename = oss.str();
					tmpm.outputToTxtTabs(filename);
					cout << "Marginal 2 output to " << filename << endl;
				}
				
				assert ( vecSliceDims.size() == vecSlicePts.size());
				
				for (size_t s = 0; s < vecSliceDims.size(); ++s) {
					
					std::vector < int > sliceDims = vecSliceDims[s];
					std::vector < double > slicePts = vecSlicePts.at(s);
				
					PiecewiseConstantFunction tmps = mcmcPCF.makeSlice(sliceDims, slicePts);
					//normalise it!
					tmps.normalise();
					ostringstream oss;
					oss << histFilenameBase << "_slice_d"; 
					for (size_t k = 0; k < sliceDims.size(); ++k) oss << "_" << sliceDims[k]; 
					oss << "_p";
					for (size_t k = 0; k < sliceDims.size(); ++k) oss << "_" << slicePts[k]; 
					oss << "_n" << n << "_l" << l << "_r" << (rep+1) << ".txt";
					string filename = oss.str();
					tmps.outputToTxtTabs(filename);
					cout << "Slice output to " << filename << endl;
				}

			}
			if (!vecSliceDims.empty()) {
				/* only do coverage regions if we are slicing */
				for (size_t ci = 0; ci < covs.size(); ++ci) {
					
					real cov = covs[ci];
					
					std::string covFileName;
					ostringstream oss;
					oss << histFilenameBase << "_coverage_c" << _double(cov) 
						<< "_n" << n << "_l" << l << "_r" << (rep+1) << ".txt";
					string filename = oss.str();
					
					mcmcPCF.outputCoverageRegion(filename, cov);
				
				}
			}
			
		}
		
		
			
		/*get quasi random points in the box */
		std::vector < std::vector < real > > qrPts;
		getQuasiRandomPoints( box, qrPts, intN);
		
		
		/* get points from true density */
		std::vector < std::vector < real > > intPts;
		mixMVNptr->prn(intPts, intN);
		
		
		/*get MCMC histogram densities at the remaining integration points and log results*/
		
		std::vector < real > estDensitiesMCMC_IS;
		std::vector < real > estDensitiesMCMC_QR;
		
		{
			PiecewiseConstantFunction pcfSmeared 
			= mcmcPCF.makeSmearZeroValues(1/(1000000.0));
			
			getPCFDensitiesCensor(pcfSmeared, intPts, // intPts is censored in this process
					estDensitiesMCMC_IS, timingIntDensities[rep]);
			intNcensored = intPts.size();
			getPCFDensities(pcfSmeared, qrPts, estDensitiesMCMC_QR);
		}
		
		/*get true densities at the remaining integration points points */
		std::vector < real > trueIntPtDensities_IS;
		getTrueDensities(*mixMVNptr, intPts, trueIntPtDensities_IS, timingIntDensities[rep]);
		/*get true densities at the qr points */
		std::vector < real > trueIntPtDensities_QR;
		getTrueDensities(*mixMVNptr, qrPts, trueIntPtDensities_QR);
		/*get average log true den and pcf den and ratios*/
		real avLogTrueDen = avLogDen(trueIntPtDensities_IS);
		avLogDens[rep].push_back(avLogTrueDen);
		avLogDenRatios[rep].push_back(0.0);
		estL1ErrorsQR[rep].push_back(0.0);
		
		{	
			real avLogEstDen = avLogDen(estDensitiesMCMC_IS);
			avLogDens[rep].push_back(avLogEstDen);
			avLogDenRatios[rep].push_back(avLogTrueDen - avLogEstDen);
			/*approx L1 errors*/
			 real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensitiesMCMC_QR);
			estL1ErrorsQR[rep].push_back(estL1_QR);
		}
		
		subpavings::outputFileStart(repLogFilename);
	
		ofstream os(repLogFilename.c_str(), ios::app);         // append
		if (os.is_open()) {
			os << "Timing\t" << timingMake[rep].back() << endl;
			os << "Leaves\t" << leaves[rep].back() << endl;
			os << "KL\t" << avLogDenRatios[rep].back() << endl;
			os << "L1\t" << estL1ErrorsQR[rep].back() << endl;
					
			os.close();
		}
		else {
			std::cerr << "Error: could not open file named "
				<< repLogFilename << std::endl << std::endl;
		}
		
	}//end reps
	
		
	std::string logFilename;
	{
		ostringstream oss;
		oss << logFilenameBase << "_n" << n 
				<< "_using" << intNcensored << "_" << intN << ".txt";
		logFilename = oss.str();
	}
	
	outputResults(logFilename,
			timingMake,
			timingIntDensities,
			avLogDens,
			avLogDenRatios,
			estL1ErrorsQR,
			leaves);
	
}

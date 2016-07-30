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
using namespace subpavings::kde;
using namespace std;


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

/* overloaded version:
 * - makes its own root box 
 * - use given thinout and samples*/
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

/* overloaded version:
 * - uses given root box
 * - use default thinout and samples */
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
						
/* the actual code that does the MCMC histogram */
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


/*code to do prior selection*/
void doPriorSelect(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t n,
		size_t intN,
		int K, 
		std::vector<double> Temperatures,
		int MaxTempIterations,
		const string& tempFilenameBase,	
		const string& MCMCFilenameBase,	
		const string& logFilenameBase,	
		const string& logMCMCFilenameBase,	
		const string& logCatFilenameBase,	
		MixtureMVN* mixMVNptr,
		bool ComputeCVError,
		unsigned int thinout,
		unsigned int samplesNeeded)
{
	size_t intNcensored = intN;
	size_t intNcensoredMCMC = intN;
	size_t minPoints = 1;
	
	double t_opt = 0.3;//initiate
	const double lowest_double = -std::numeric_limits<double>::max();
	double AvgHeldOutLkls_opt = lowest_double;
	//-numeric_limits<double>::max();//-100000000000000000;
	//cout << "lowest double = " << lowest_double << endl;
				
	/* containers for stuff we will be storing */
	std::vector < std::vector < double> > topt(reps);
	// will be output in logMCMCFilename

	// for temp iterate and RPQ
	std::vector < std::vector < double > > timings(reps);
	std::vector < std::vector < double > > timingIntDensities(reps);
	std::vector < std::vector < real > > avLogDens(reps);
	std::vector < std::vector < real > > avLogDenRatios(reps);
	std::vector < std::vector < real > > estL1ErrorsQR(reps);
	std::vector < std::vector < size_t > > leaves(reps);
	
	// for MCMC
	std::vector < std::vector < double > > timingMakeMCMC(reps);
	std::vector < std::vector < double > > timingDensitiesMCMC(reps);
	std::vector < std::vector < double > > timingIntDensitiesMCMC(reps);
	std::vector < std::vector < real > > avLogDensMCMC(reps);
	std::vector < std::vector < real > > avLogDenRatiosMCMC(reps);
	std::vector < std::vector < real > > estL1ErrorsQRMCMC(reps);
	std::vector < std::vector < size_t > > leavesMCMC(reps);

	/* Start the repetitions */
	for (int rep = 0; rep < reps; ++rep) 
	{
		cout << "\n =========== Rep number " 
			 << (rep+1) << "============" << endl;

		/* a file for each repetition */
		std::string endFilename;
		{
			ostringstream oss;
			oss << "_n" << n << "_r" << (rep+1) << ".txt";
			endFilename = oss.str();
		}
		string tempLogFilename = tempFilenameBase + "RepLog" + endFilename;
	
		long unsigned int dataseed = seed+rep;
		mixMVNptr->resetPRNG(dataseed);
	
		std::vector < std::vector < double > > simdata;
		
		cout << "\nGenerate " << n << " random values:" << endl;
		mixMVNptr->prn(simdata, n);
		//cout << simdata[0] << endl;
		
		// for computations of distances	
		// get points from true density 
		std::vector < std::vector < real > > intPts;
		mixMVNptr->prn(intPts, intN);
		
		/* set up switches */
		bool successfulInsertion = false;
		bool successfulPQSplit = false;
		// set to do Carver PQ + SEBPB for posterior maximization
		bool CarvingMaxPosterior=true;// false means using SEB PQ with minChildPoints
		bool TempIterate=true;//true;
		
		/* set up prior distribution object */
		//LogTemperaturePrior logPrior(10.0);
		LogCatalanTempPrior logPrior(10.0); 
		// t \in (-infinity, infinity)\-2
		//play around with some range
		// n=3000,d=2: 1/100 = undersmoother, 10 is oversmoothed/raggedy, (1/10,1.0) is ok

		/* set up containers for histograms and PCFs of histograms */
		//std::vector< subpavings::AdaptiveHistogram* > hists;
		//std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
		
		{
		/* insert data into an histogram object */
		cout << "\nInsert data and make root box:" << endl; 
		AdaptiveHistogram adhA0; // let it make its own root box
		successfulInsertion = false;
		successfulInsertion = adhA0.insertRvectorsFromVectorOfVecDbls(simdata);
		adhA0.setHoldAllStats(true);
		if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
		size_t n = adhA0.getRootCounter();
		size_t d = adhA0.getDimensions ();
		//cout << n << "\t" << d << endl;
		
		/* begin getting best histogram for first big burst */
		//std::string pqFileNameBase = "pq";
		std::string postFileName = "";
		/*{
			ostringstream oss;
			oss << "Posteriors" << pqFileNameBase << "_d" << d << "_n" << n 
			<< "_r" << (rep+1) << ".txt";
			postFileName = oss.str();
		}*/
   
		std::string checkPostFileNameBase = "";
		/*{
			ostringstream oss;
			oss << "CheckPosterior" << pqFileNameBase << "_d" << d  << "_n" 
			<< n  << "_r" << (rep+1) << ".txt";
			checkPostFileNameBase = oss.str();
		}*/
		
		/* parameters for function findStartingPointsBest */
		int precPQ = 5;
		bool stopOnMaxPosterior = true;//bool for carver PQ
		int keep = 1;
		int chooseStarts = 10;
		clock_t starttime = clock();
		size_t critSEB;
		size_t maxLeavesSEB;
		size_t maxLeavesCarving;
		
		/* containers to hold required results */      
		RealVec AvgHeldOutLkls;
		RealVec AvgHeldOutEmpiricalDeviations;
		RealVec AvgEstL1CV;
		RealVec AvgEstKLCV;
		cxsc::real CVgain=0.0;
  		
		/* optimizing over smoothing parameter Temp */  
		int TempIterations=0;
		unsigned long int seedStarts = seed+TempIterations;  

		if(TempIterate)
		{
				
			/* set up CV parameters */
			const size_t N = adhA0.getRootCounter();
			//size of successfully inserted transformed data
			const size_t KofN = N/K; //cout << KofN << endl; getchar();
			size_t nTrain;
			adhA0.clearAllHistData();//clear the data to make space during CV
			//containers to keep the Training and validation data
            std::vector < std::vector < double > > DataT;
            std::vector < std::vector < double > > DataV;
            
            /* set up for permutations for random shuffling for CV */
			const gsl_rng_type * T;
			gsl_rng * r;
			gsl_permutation * p = gsl_permutation_alloc (N);
			gsl_permutation * q = gsl_permutation_alloc (N);
			gsl_rng_env_setup();
			T = gsl_rng_default;
			r = gsl_rng_alloc (T);
			//printf ("initial permutation:");  
			gsl_permutation_init (p);
			//gsl_permutation_fprintf (stdout, p, " %u"); printf ("\n"); getchar();   
		
			do
			{
				LogCatalanTempPrior logPrior(Temperatures[TempIterations]);
				seedStarts += TempIterations;
				cout << "\n=============================" << endl;
				cout << "seedStarts: " << seedStarts 
					<< "\t Temperature: " 
					<< Temperatures[TempIterations] << endl; 
				//getchar();
          
				real HeldOutLkl = 0.0;
				real HeldOutEmpiricalDeviation = 0.0;
				real EstL1CV = 0.0;
				real EstKLCV = 0.0;

				// a container for our histograms at various temperatures
				AdaptiveHistogram adhA0cv(adhA0.getRootBox()); 
				// make adh for CV with root box from adhA0
				for (int cvI=1; cvI<=K; cvI++)//K-fold CV loop
				{ 
					cout << "----------- CV: " << cvI  
						<< "\t Temperature: " << Temperatures[TempIterations] 
						<< "\t Rep: " << rep + 1 << "----------" << endl;
						
					//getchar(); 
					//randomly shuffles the order of N objects, 
					//each of size size_t, stored in the array base[0..N-1]. 
					gsl_ran_shuffle (r, p->data, N, sizeof(size_t)); 
					//gsl_permutation_fprintf (stdout, p, " %u"); printf ("\n"); getchar();   
					
					successfulInsertion = false;
					
					std::vector< subpavings::AdaptiveHistogram* > histsT;
				
					//insert the training and validation data into Cv containers
					for(size_t i=0; i<KofN; i++) 
					DataV.push_back(simdata[gsl_permutation_get(p,i)]);
					
					for(size_t i=KofN; i<N; i++) 
					DataT.push_back(simdata[gsl_permutation_get(p,i)]);
					//cout << DataT[0] << endl;
			  
					successfulInsertion = adhA0cv.insertRvectorsFromVectorOfVecDbls(DataT);
					//insert training data
					if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
				
					// some guesses for max points in a node to stop posterior queue 
					nTrain = adhA0cv.getRootCounter(); 
					//cout << "nTrain = " << nTrain << endl; getchar();
					critSEB = static_cast<size_t>(std::log(static_cast<double>(nTrain)));
					//can be as low as 1
	    
					// some guesses for maximum leaves we'll let SEB queue go to 
					maxLeavesSEB = nTrain;// / critSEB; // integer division
					maxLeavesCarving = maxLeavesSEB / 2; // integer division
					SPSNodeMeasureVolMassMinus compCarving(nTrain);
					AdaptiveHistogram::PrioritySplitQueueEvaluator 
					evaluatorCarving(compCarving, maxLeavesCarving);
					
					SPSNodeMeasureCount compSEB;
					AdaptiveHistogram::PrioritySplitQueueEvaluator 
					evaluatorSEB(compSEB, critSEB, maxLeavesSEB);
				
					cout << "\n Start findStartingPointsBest..."  << endl;
					CarverSEB::findStartingPointsBest(adhA0cv, histsT, 
						evaluatorCarving, evaluatorSEB, 
						logPrior, minPoints, chooseStarts, 
						keep, stopOnMaxPosterior, 
						postFileName, checkPostFileNameBase, 
						precPQ, seedStarts);
					
					/* make histsT[0] into a PCF object */
					//PiecewiseConstantFunction& pcfT(*histsT[0]);
					PiecewiseConstantFunction* pcfT = new PiecewiseConstantFunction(*histsT[0]);
					pcfT->smearZeroValues(0.0000001);
					//assert (pcfT.getTotalIntegral() == cxsc::real(1.0));
					cout << "Number of leaves (T): " << pcfT->getRootLeaves() << endl;
					
					histsT[0]->clearAllHistData();//clear the data from training big burst
					adhA0cv.clearAllHistData(); //also clear the data from this histogram
					adhA0cv.mergeUp();//Merge the possibly multileaf cv histogram up to just root. 
					cout << "Cleared data from histsT[0], adhA0cv with t data" << endl;
					//getchar(); 
	
					histsT[0]->insertRvectorsFromVectorOfVecDbls(DataV);//insert validation data
					HeldOutLkl += pcfT->getLogLikelihood(*histsT[0]);
					
					PiecewiseConstantFunction* pcfV = new PiecewiseConstantFunction(*histsT[0]);
					HeldOutEmpiricalDeviation += pcfT->getL1Distance(*pcfV);
					//cout << HeldOutLkl << '\t' << HeldOutEmpiricalDeviation << endl; getchar();
				
					histsT[0]->clearAllHistData();//clear the data from validation big burst
					cout << "Cleared data from histsT[0] with v data" << endl;
					//histsT[0]->mergeUp();//merge up to root
					
					cout << "Number of leaves (V): " << pcfV->getRootLeaves() << endl;
					
					if(ComputeCVError) {
					
					cout << "Getting the L1-distance and KL-distance... " << endl;
					/* Get the number of leaves, L1-distance and KL-distance */					
					// get quasi random points in the box 
					ivector box = pcfT->getRootBox();
					std::vector < std::vector < real > > qrPts;
					getQuasiRandomPoints(box, qrPts, intN);
				
					// get points from true density 
					//std::vector < std::vector < real > > intPts;
					//mixMVNptr->prn(intPts, intN);

					// get histogram densities at the remaining integration points and log results
					std::vector < real > estDensities_IS;
					std::vector < real > estDensities_QR;	
					{
						PiecewiseConstantFunction* pcfSmeared = new PiecewiseConstantFunction();
						*pcfSmeared = pcfT->makeSmearZeroValues(1/(1000000.0));
						//PiecewiseConstantFunction pcfSmeared 
						//= pcfT->makeSmearZeroValues(1/(1000000.0));
						
						getPCFDensitiesCensor(*pcfSmeared, intPts, // intPts is censored in this process
							estDensities_IS, timingIntDensities[rep]);
						intNcensored = intPts.size();
						getPCFDensities(*pcfSmeared, qrPts, estDensities_QR);
							
						if (NULL != pcfSmeared) delete pcfSmeared; 
						pcfSmeared = NULL;
					}
		
					/*get true densities at the remaining integration points points */
					std::vector < real > trueIntPtDensities_IS;
					getTrueDensities(*mixMVNptr, intPts, trueIntPtDensities_IS, timingIntDensities[rep]);
					/*get true densities at the qr points */
					std::vector < real > trueIntPtDensities_QR;
					getTrueDensities(*mixMVNptr, qrPts, trueIntPtDensities_QR);
					
					/*get average log true den and pcf den and ratios*/
					real avLogTrueDen = avLogDen(trueIntPtDensities_IS);		
					{	
						real avLogEstDen = avLogDen(estDensities_IS);
						EstKLCV += (avLogTrueDen - avLogEstDen);
						/*approx L1 errors*/
						real boxVol = realVolume(box);
						real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensities_QR);
						EstL1CV += estL1_QR;
			
						cout << "KL: \t L1 : \n" << endl;
						cout << avLogTrueDen - avLogEstDen << "\t" << estL1_QR << endl;
						
						// some cleaning up
						trueIntPtDensities_IS.clear();
						trueIntPtDensities_QR.clear();
						//intPts.clear();
						qrPts.clear();
						estDensities_IS.clear();
						estDensities_QR.clear();
						
					}					
					
					} // end of ComputeCVError

					
					cout << "Cleaning up..." << endl;
					//getchar();
					
					//to free all the contents of histsT, DataV, DataT, pcfV, pcfT
					if (NULL != pcfT) delete pcfT; 
					pcfT = NULL;
					
					if (NULL != pcfV) delete pcfV; 
					pcfV = NULL;
					
					DataV.clear(); DataT.clear();
					
					for (size_t i = 0; i < histsT.size(); ++i)
					{
						if (NULL != histsT[i]) delete histsT[i];
						histsT[i] = NULL;
					} 

					
				} //end of CV for loop

				AvgHeldOutLkls.push_back(HeldOutLkl/double(K));
				AvgHeldOutEmpiricalDeviations.push_back(HeldOutEmpiricalDeviation/double(K));
				
				if(ComputeCVError) 
				{
					AvgEstL1CV.push_back(EstL1CV/double(K));
					AvgEstKLCV.push_back(EstKLCV/double(K));
				}
				else
				{
					AvgEstL1CV.push_back(0);
					AvgEstKLCV.push_back(0);
				}
				
				TempIterations++;
				
			} // end of do
        
			while (TempIterations<MaxTempIterations);// && CVgain>0.1);
				//cout << "Temp: " << Temperatures << endl;
				//cout << "AvgHeldOutLkls:"  << AvgHeldOutLkls << endl 
					//<< "AvgHeldOutEmpiricalDeviations: " << AvgHeldOutEmpiricalDeviations 
					//<< endl;
				//cout << "Temp Iteration Number " << TempIterations << endl; getchar();
				//output the L1-distance and KL distance: plot: get a concave down 
				
			subpavings::outputFileStart(tempLogFilename);
				ofstream os(tempLogFilename.c_str(), ios::app); // append
				if (os.is_open()) {
					os << "Temperature \t AvgHeldOutLkls \t AvgKL \t AvgL1 " << endl;
					for (int i = 0; i <MaxTempIterations; i++)
					{
						os << Temperatures[i] << "\t" 
						   << AvgHeldOutLkls[i] << "\t" 
						   << AvgEstKLCV[i] << "\t" 
						   << AvgEstL1CV[i] << endl;
					}
					//os << "Timing\t" << timingMake[rep].back() << endl;
					//os << "Leaves\t" << leaves[rep].back() << endl;
					//os << "KL\t" << avLogDenRatios[rep].back() << endl;
					//os << "L1\t" << estL1ErrorsQR[rep].back() << endl;
					os.close();
				}
				else {
				std::cerr << "Error: could not open file named "
					<< tempLogFilename << std::endl << std::endl;
				}	

				gsl_permutation_free (p);
				gsl_rng_free (r);		
		} // end of if(TempIterate)
		
		/*=============================================================
		Get the optimal temperature based on the largest AvgHeldOutLkls 
		=============================================================*/
		
		if(TempIterate)
		{
		
			for(int i=0; i < Temperatures.size(); i++)
			{
				if(AvgHeldOutLkls[i]>AvgHeldOutLkls_opt)
				{
					AvgHeldOutLkls_opt = _double(AvgHeldOutLkls[i]);
					t_opt=_double(Temperatures[i]);
				}
			}
			
			topt[rep].push_back(t_opt);
			
			cout << "******************************************** \n"   
				<< "Optimal temperature and AvgHeldOutLkls are : " 
				<< t_opt << '\t' << AvgHeldOutLkls_opt << endl; 
			cout << "Build histogram based on optimal temperature... \n" 
				<< endl;
		}

		/*=============================================================
		 * build a histogram based on the optimal temperature 
		=============================================================*/
		{
		// some guesses for max points in a node to stop posterior queue 
		if(adhA0.getRootCounter()==0)
		{  
			successfulInsertion = false;
			successfulInsertion = adhA0.insertRvectorsFromVectorOfVecDbls(simdata);//insert all data
			if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
		}
	
		critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));//can be as low as 1
		//cout << critSEB << endl;
	
		/* set up containers for histograms and PCFs of histograms */
		std::vector< subpavings::AdaptiveHistogram* > hists;
		std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
		
		// some guesses for maximum leaves we'll let SEB queue go to 
		maxLeavesSEB = n;//*adhA0.getDimensions();// / critSEB; // integer division
		maxLeavesCarving = maxLeavesSEB;//*adhA0.getDimensions();///2; // integer division
		SPSNodeMeasureVolMassMinus compCarving(n);
		AdaptiveHistogram::PrioritySplitQueueEvaluator 
		evaluatorCarving( compCarving, maxLeavesCarving);
		
		SPSNodeMeasureCount compSEB;
		AdaptiveHistogram::PrioritySplitQueueEvaluator 
		evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
		
		LogCatalanTempPrior logPrior(t_opt);
		if(CarvingMaxPosterior)
		{ 
			CarverSEB::findStartingPointsBest(adhA0, hists, 
						evaluatorCarving, evaluatorSEB, logPrior, 
						minPoints, chooseStarts, keep, stopOnMaxPosterior, 
						postFileName, checkPostFileNameBase, precPQ, 
						dataseed);
		}
		else
		{
              // function object to compare nodes on count
              // ie split node with largest count first
              CompCount compCount;
              CompVol compVol;
              //  minimum points to use when splitting.
              // A node will not be splittable if either child would then have
              // < minPoints of data associated with it. 
              size_t minChildPoints = 80;//static_cast<size_t>(std::log(static_cast<double>(n)))
              double minVolume = 0.01;
              adhA0.prioritySplit(compCount, maxLeavesSEB, NOLOG, minChildPoints, minVolume);
              //adhA0.prioritySplit(compVol, maxLeavesSEB, NOLOG, minVolume);
              hists.push_back(& adhA0);
              //AdaptiveHistogram* adhPtr = & adhA0;
              //hists.push_back(adhPtr);
		}
		clock_t endtime = clock();	
		double timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
		timings[rep].push_back(timingStarts);
		cout << "time to get starts = " << timingStarts << endl;

		
		assert(hists.size()==1);
		//hists[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt", 6,true);
		// create a name for the file to output
		// To realize a file output
		//make a piecewise constant function of the best adaptive histogram density estimate
		PiecewiseConstantFunction* pcfPtr = new PiecewiseConstantFunction(*hists[0]);
		adhA0.clearAllHistData(); 
		adhA0.mergeUp();
		cout << "Cleared data in adhA0." << endl; //getchar(); 
		pcfs.push_back(pcfPtr);
		//pcfs[0]->smearZeroValues(0.0000001);
		cxsc::real integral0 = pcfs[0]->getTotalIntegral();
		leaves[rep].push_back(pcfs[0]->getRootLeaves());
		cout << "Number of leaves: " << pcfs[0]->getRootLeaves() << endl;
		//pcfs[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt"+"pcf", 6,true);
		cout << "pcfs[0]->getTotalIntegral() = " << integral0 << endl;
		assert (integral0 == cxsc::real(1.0));
    
		/* Get L1-distance and KL-distance */
		cout << "Getting the L1-distance and KL-distance: " << endl;
		
		// get quasi random points in the box 
		ivector box = pcfs[0]->getRootBox();
		std::vector < std::vector < real > > qrPts;
		intNcensored = intN;
		getQuasiRandomPoints( box, qrPts, intN);
				
		// get points from true density 
		//std::vector < std::vector < real > > intPts;
		//mixMVNptr->prn(intPts, intN);
		//for (int i = 0; i < intPts.size(); i++) {cout << intPts[i] << endl;}
		
		/*get histogram densities at the remaining integration points and log results*/
		std::vector < real > estDensities_IS;
		std::vector < real > estDensities_QR;
		
		{
			PiecewiseConstantFunction pcfSmeared 
			= pcfs[0]->makeSmearZeroValues(1/(1000000.0));			
			getPCFDensitiesCensor(pcfSmeared, intPts, // intPts is censored in this process
					estDensities_IS, timingIntDensities[rep]);
			intNcensored = intPts.size();
			getPCFDensities(pcfSmeared, qrPts, estDensities_QR);
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
			real avLogEstDen = avLogDen(estDensities_IS);
			avLogDens[rep].push_back(avLogEstDen);
			avLogDenRatios[rep].push_back(avLogTrueDen - avLogEstDen);
			/*approx L1 errors*/
		    real boxVol = realVolume(box);
		    real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensities_QR);
			estL1ErrorsQR[rep].push_back(estL1_QR);
			
			cout << "KL: \t L1 : \n" << endl;
			cout << avLogTrueDen - avLogEstDen << "\t" << estL1_QR << endl;
			
			// some cleaning up
			//cout << "before cleaning: " << endl; getchar();
			trueIntPtDensities_IS.clear();
			trueIntPtDensities_QR.clear();
			//intPts.clear();
			qrPts.clear();
			estDensities_IS.clear();
			estDensities_QR.clear();

		} 
		
		//to free all the contents of pcfs at the end
		for (size_t i = 0; i < pcfs.size(); ++i) 
		{
			if (NULL != pcfs[i]) delete pcfs[i];
			pcfs[i] = NULL;
		}
    
		if(CarvingMaxPosterior)
		{ 
			//to free all the contents of hists at the end
			for (size_t i = 0; i < hists.size(); ++i) 
			{
				if (NULL != hists[i]) delete hists[i];
				hists[i] = NULL;
			}
		} 
		
		} // end of getting histogram using optimal temperature
	
		{
		/*=========================================================== 
		  Perform MCMC using optimal temperature 
		=============================================================*/
		//hardcode
		//t_opt = 0.0;
		cout << "==========Begin MCMC for " << t_opt << "=========" << endl;
		
		std::string endFilename;
		{
			ostringstream oss;
			oss << "_n" << n << "_r" << (rep+1) << ".txt";
			endFilename = oss.str();
		}
	
        /* a file for each repetition */
		string repMCMCLogFilename = MCMCFilenameBase + "RepLog" + endFilename;
	
		/* Perform the MCMC */
		cout << "\nGet pcf from mcmc" << endl;
		LogCatalanTempPrior logPrior(t_opt);
		PiecewiseConstantFunction mcmcPCF = doHistMCMC(rep, n, minPoints,	
					MCMCFilenameBase,
					simdata,
					timingMakeMCMC[rep],
					seed,
					thinout,
					samplesNeeded,
					logPrior);
		
		/* Perform distance calculations */
		ivector box = mcmcPCF.getRootBox();
		real boxVol = realVolume(box);
		cout << "\npcf box volume is " << boxVol << endl;
		
		{
			size_t l = mcmcPCF.getRootLeaves();
			cout << "\nMCMC average with " << l << " leaves"<< endl;
			
			leavesMCMC[rep].push_back(l);
					
			/*{
				ostringstream oss;
				oss << MCMCFilenameBase << "_n" << n << "_l" <<l 
				* 	<< "_r" << (rep+1) << ".txt";
				string filename = oss.str();
				mcmcPCF.outputToTxtTabs(filename);
				cout << "Average output to " << filename << endl;
			}*/
		}
		
		/*get quasi random points in the box */
		std::vector < std::vector < real > > qrPts;
		getQuasiRandomPoints( box, qrPts, intN);
		
		/* get points from true density */
		//std::vector < std::vector < real > > intPts;
		//mixMVNptr->prn(intPts, intN);
		//for (int i = 0; i < intPts.size(); i++) {cout << intPts[i] << endl;}
				
		/*get MCMC histogram densities at the remaining integration points and log results*/		
		std::vector < real > estDensitiesMCMC_IS;
		std::vector < real > estDensitiesMCMC_QR;
		
		{
			PiecewiseConstantFunction pcfSmeared 
			= mcmcPCF.makeSmearZeroValues(1/(1000000.0));			
			getPCFDensitiesCensor(pcfSmeared, intPts, // intPts is censored in this process
					estDensitiesMCMC_IS, timingIntDensitiesMCMC[rep]);
			//cout << intPts.size() << endl; getchar();
			intNcensoredMCMC = intPts.size();
			getPCFDensities(pcfSmeared, qrPts, estDensitiesMCMC_QR);
		}
		
		/*get true densities at the remaining integration points points */
		std::vector < real > trueIntPtDensities_IS;
		getTrueDensities(*mixMVNptr, intPts, trueIntPtDensities_IS, 
							timingIntDensitiesMCMC[rep]);
		/*get true densities at the qr points */
		std::vector < real > trueIntPtDensities_QR;
		getTrueDensities(*mixMVNptr, qrPts, trueIntPtDensities_QR);
		/*get average log true den and pcf den and ratios*/
		real avLogTrueDen = avLogDen(trueIntPtDensities_IS);
		avLogDensMCMC[rep].push_back(avLogTrueDen);
		avLogDenRatiosMCMC[rep].push_back(0.0);
		estL1ErrorsQRMCMC[rep].push_back(0.0);
		
		{	
			real avLogEstDen = avLogDen(estDensitiesMCMC_IS);
			avLogDensMCMC[rep].push_back(avLogEstDen);
			avLogDenRatiosMCMC[rep].push_back(avLogTrueDen - avLogEstDen);
			/*approx L1 errors*/
			real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensitiesMCMC_QR);
			estL1ErrorsQRMCMC[rep].push_back(estL1_QR);
		}
		
		subpavings::outputFileStart(repMCMCLogFilename);
	
		ofstream os(repMCMCLogFilename.c_str(), ios::app);         // append
		if (os.is_open()) {
			os << "Timing\t" << timingMakeMCMC[rep].back() << endl;
			os << "Leaves\t" << leavesMCMC[rep].back() << endl;
			os << "KL\t" << avLogDenRatiosMCMC[rep].back() << endl;
			os << "L1\t" << estL1ErrorsQRMCMC[rep].back() << endl;
					
			os.close();
		}
		else {
			std::cerr << "Error: could not open file named "
				<< repMCMCLogFilename << std::endl << std::endl;
		}
		
			// some cleaning up
			//cout << "before cleaning: " << endl; getchar();
			trueIntPtDensities_IS.clear();
			trueIntPtDensities_QR.clear();
			//intPts.clear();
			qrPts.clear();
			estDensitiesMCMC_IS.clear();
			estDensitiesMCMC_QR.clear();
			
		}//end of MCMC
		
		} // end of prior selection, temp iterate		

		simdata.clear(); // clear the simulated data for this rep
		
		cout << "============ end of rep " 
			<< rep+1 << "==============\n" << endl;
			
	} // end of for loop for reps
	 
	 
	// output results
	std::string logFilename;
	{
		ostringstream oss;
		oss << logFilenameBase << "_n" << n 
				<< "_using" << intNcensored << "_" << intN << ".txt";
		logFilename = oss.str();
	}

	cout << "Output results to: " << logFilename << endl;
	outputResults(logFilename,
			timings,
			avLogDens,
			avLogDenRatios,
			estL1ErrorsQR,
			leaves);	


	//output results for MCMC
	std::string logMCMCFilename;
	{
		ostringstream oss;
		oss << logMCMCFilenameBase << "_n" << n 
				<< "_using" << intNcensoredMCMC << "_" << intN << ".txt";
		logMCMCFilename = oss.str();
	}
	
	cout << "Output results to: " << logMCMCFilename << endl;
	outputResults(logMCMCFilename,
			timingMakeMCMC,
			timingIntDensitiesMCMC,
			avLogDensMCMC,
			avLogDenRatiosMCMC,
			estL1ErrorsQRMCMC,
			leavesMCMC,
			topt);
}

/* code to do RPQ histogram using LogCatalanPrior */
void doRPQCat(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t n,
		size_t intN,
		const string& logCatFilenameBase,	
		MixtureMVN* mixMVNptr)
{
	size_t intNcensored = intN;
	size_t minPoints = 1;

				
	/* containers for stuff we will be storing */
	std::vector < std::vector < double > > timings(reps);
	std::vector < std::vector < double > > timingIntDensities(reps);
	std::vector < std::vector < real > > avLogDens(reps);
	std::vector < std::vector < real > > avLogDenRatios(reps);
	std::vector < std::vector < real > > estL1ErrorsQR(reps);
	std::vector < std::vector < size_t > > leaves(reps);

	/* Start the repetitions */
	for (int rep = 0; rep < reps; ++rep) 
	{
		cout << "\n =========== Rep number " 
			 << (rep+1) << "============" << endl;
	
		long unsigned int dataseed = seed+rep;
		mixMVNptr->resetPRNG(dataseed);
	
		std::vector < std::vector < double > > simdata;
		
		cout << "\nGenerate " << n << " random values:" << endl;
		mixMVNptr->prn(simdata, n);
		//cout << simdata[0] << endl;
		
		// for computations of distances	
		// get points from true density 
		std::vector < std::vector < real > > intPts;
		mixMVNptr->prn(intPts, intN);
		
		/* set up switches */
		bool successfulInsertion = false;
		bool successfulPQSplit = false;
		// set to do Carver PQ + SEBPB for posterior maximization
		bool CarvingMaxPosterior=true;// false means using SEB PQ with minChildPoints
		
		/* insert data into an histogram object */
		cout << "\nInsert data and make root box:" << endl; 
		AdaptiveHistogram adhA0; // let it make its own root box
		successfulInsertion = false;
		successfulInsertion = adhA0.insertRvectorsFromVectorOfVecDbls(simdata);
		adhA0.setHoldAllStats(true);
		if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
		size_t n = adhA0.getRootCounter();
		size_t d = adhA0.getDimensions ();
		//cout << n << "\t" << d << endl;
		
		/* begin getting best histogram for first big burst */
		//std::string pqFileNameBase = "pq";
		std::string postFileName = "";
		/*{
			ostringstream oss;
			oss << "Posteriors" << pqFileNameBase << "_d" << d << "_n" << n 
			<< "_r" << (rep+1) << ".txt";
			postFileName = oss.str();
		}*/
   
		std::string checkPostFileNameBase = "";
		/*{
			ostringstream oss;
			oss << "CheckPosterior" << pqFileNameBase << "_d" << d  << "_n" 
			<< n  << "_r" << (rep+1) << ".txt";
			checkPostFileNameBase = oss.str();
		}*/
		
		/* parameters for function findStartingPointsBest */
		int precPQ = 5;
		bool stopOnMaxPosterior = true;//bool for carver PQ
		int keep = 1;
		int chooseStarts = 10;
		clock_t starttime = clock();
		size_t critSEB;
		size_t maxLeavesSEB;
		size_t maxLeavesCarving;
		
		/* set up containers for histograms and PCFs of histograms */
		std::vector< subpavings::AdaptiveHistogram* > hists;
		std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
		
		
		/*=============================================================
		 * build an RPQ histogram based on the LogCatalanPrior 
		=============================================================*/
		{
		// some guesses for max points in a node to stop posterior queue 
		if(adhA0.getRootCounter()==0)
		{  
			successfulInsertion = false;
			successfulInsertion = adhA0.insertRvectorsFromVectorOfVecDbls(simdata);//insert all data
			if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
		}
	
		critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));//can be as low as 1
		//cout << critSEB << endl;
	
		/* set up containers for histograms and PCFs of histograms */
		std::vector< subpavings::AdaptiveHistogram* > hists;
		std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
		
		// some guesses for maximum leaves we'll let SEB queue go to 
		maxLeavesSEB = n;//*adhA0.getDimensions();// / critSEB; // integer division
		maxLeavesCarving = maxLeavesSEB;//*adhA0.getDimensions();///2; // integer division
		SPSNodeMeasureVolMassMinus compCarving(n);
		AdaptiveHistogram::PrioritySplitQueueEvaluator 
		evaluatorCarving( compCarving, maxLeavesCarving);
		
		SPSNodeMeasureCount compSEB;
		AdaptiveHistogram::PrioritySplitQueueEvaluator 
		evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
		
		LogCatalanPrior logPrior;
		if(CarvingMaxPosterior)
		{ 
			CarverSEB::findStartingPointsBest(adhA0, hists, 
						evaluatorCarving, evaluatorSEB, logPrior, 
						minPoints, chooseStarts, keep, stopOnMaxPosterior, 
						postFileName, checkPostFileNameBase, precPQ, 
						dataseed);
		}
		else
		{
              // function object to compare nodes on count
              // ie split node with largest count first
              CompCount compCount;
              CompVol compVol;
              //  minimum points to use when splitting.
              // A node will not be splittable if either child would then have
              // < minPoints of data associated with it. 
              size_t minChildPoints = 80;//static_cast<size_t>(std::log(static_cast<double>(n)))
              double minVolume = 0.01;
              adhA0.prioritySplit(compCount, maxLeavesSEB, NOLOG, minChildPoints, minVolume);
              //adhA0.prioritySplit(compVol, maxLeavesSEB, NOLOG, minVolume);
              hists.push_back(& adhA0);
              //AdaptiveHistogram* adhPtr = & adhA0;
              //hists.push_back(adhPtr);
		}
		clock_t endtime = clock();	
		double timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
		timings[rep].push_back(timingStarts);
		cout << "time to get starts = " << timingStarts << endl;

		
		assert(hists.size()==1);
		//hists[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt", 6,true);
		// create a name for the file to output
		// To realize a file output
		//make a piecewise constant function of the best adaptive histogram density estimate
		PiecewiseConstantFunction* pcfPtr = new PiecewiseConstantFunction(*hists[0]);
		adhA0.clearAllHistData(); 
		adhA0.mergeUp();
		cout << "Cleared data in adhA0." << endl; //getchar(); 
		pcfs.push_back(pcfPtr);
		//pcfs[0]->smearZeroValues(0.0000001);
		cxsc::real integral0 = pcfs[0]->getTotalIntegral();
		leaves[rep].push_back(pcfs[0]->getRootLeaves());
		cout << "Number of leaves: " << pcfs[0]->getRootLeaves() << endl;
		//pcfs[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt"+"pcf", 6,true);
		cout << "pcfs[0]->getTotalIntegral() = " << integral0 << endl;
		assert (integral0 == cxsc::real(1.0));
    
		/* Get L1-distance and KL-distance */
		cout << "Getting the L1-distance and KL-distance: " << endl;
		
		// get quasi random points in the box 
		ivector box = pcfs[0]->getRootBox();
		std::vector < std::vector < real > > qrPts;
		intNcensored = intN;
		getQuasiRandomPoints( box, qrPts, intN);
				
		// get points from true density 
		std::vector < std::vector < real > > intPts;
		mixMVNptr->prn(intPts, intN);
		//for (int i = 0; i < intPts.size(); i++) {cout << intPts[i] << endl;}
		
		/*get histogram densities at the remaining integration points and log results*/
		std::vector < real > estDensities_IS;
		std::vector < real > estDensities_QR;
		
		{
			PiecewiseConstantFunction pcfSmeared 
			= pcfs[0]->makeSmearZeroValues(1/(1000000.0));			
			getPCFDensitiesCensor(pcfSmeared, intPts, // intPts is censored in this process
					estDensities_IS, timingIntDensities[rep]);
			intNcensored = intPts.size();
			getPCFDensities(pcfSmeared, qrPts, estDensities_QR);
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
			real avLogEstDen = avLogDen(estDensities_IS);
			avLogDens[rep].push_back(avLogEstDen);
			avLogDenRatios[rep].push_back(avLogTrueDen - avLogEstDen);
			/*approx L1 errors*/
		    real boxVol = realVolume(box);
		    real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensities_QR);
			estL1ErrorsQR[rep].push_back(estL1_QR);
			
			cout << "KL: \t L1 : \n" << endl;
			cout << avLogTrueDen - avLogEstDen << "\t" << estL1_QR << endl;
			
			// some cleaning up
			//cout << "before cleaning: " << endl; getchar();
			trueIntPtDensities_IS.clear();
			trueIntPtDensities_QR.clear();
			//intPts.clear();
			qrPts.clear();
			estDensities_IS.clear();
			estDensities_QR.clear();

		} 
		
		} // end of getting histogram using LogCatalanPrior
	
		//to free all the contents of pcfs at the end
		for (size_t i = 0; i < pcfs.size(); ++i) 
		{
			if (NULL != pcfs[i]) delete pcfs[i];
			pcfs[i] = NULL;
		}
    
		if(CarvingMaxPosterior)
		{ 
			//to free all the contents of hists at the end
			for (size_t i = 0; i < hists.size(); ++i) 
			{
				if (NULL != hists[i]) delete hists[i];
				hists[i] = NULL;
			}
		} 

		simdata.clear(); // clear the simulated data for this rep
		
		cout << "============ end of rep " 
			<< rep+1 << "==============\n" << endl;
			
	} // end of for loop for reps
	 
	 
	// output results
	std::string logFilename;
	{
		ostringstream oss;
		oss << logCatFilenameBase << "_n" << n 
				<< "_using" << intNcensored << "_" << intN << ".txt";
		logFilename = oss.str();
	}

	cout << "Output results to: " << logFilename << endl;
	outputResults(logFilename,
			timings,
			avLogDens,
			avLogDenRatios,
			estL1ErrorsQR,
			leaves);	

}

/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief Definitions for new version of a class to do MCMC with Gelman-Rubin 
* heuristic automatic sampling rule.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCGRAutoNew.hpp"

#include "realmappedspnode.hpp"
#include "sptools.hpp"
#include "sptypes.hpp"

#include <ctime>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions

#include <gsl/gsl_randist.h>


#define MYDEBUG // extra console output for what is happening in process
//#define MYDEBUG_LOGS // log files for workings collations, averages and diffs to av as chains develop
//#define MYDEBUG_LOG_SEQUENCESTATES // to force logging of sequence states in mcmc
//#define MYDEBUG_XTRA // a few lines of extra console output
//#define MYDEBUG_OUTPUT // extra console output etc for debugging - only use for small examples!
//#define MYDEBUG_CALCS // extra console output for calculations
#define TIMEINFO // for info on loop timing

//#define FORCEFAILMCMCLOOP // debugging flag to force a failure during an MCMC loop


//#define NDEBUG // uncomment this to turn off assertion checking and all for this module only

#ifdef NDEBUG // ie only allow defines for the others if we have not defined NDEBUG for no debugging
	#undef MYDEBUG_LOGS
	#undef MYDEBUG_XTRA
	#undef MYDEBUG_OUTPUT
	#undef MYDEBUG_CALCS
	#undef MYDEBUG
	#undef TIMEINFO

	#undef FORCEFAILINSERTION 

	#undef FORCEFAILMCMCLOOP
#endif

using namespace cxsc;
using namespace subpavings;
using namespace std;

// static data members
const int MCMCGRAutoNew::defaultSamplesMonitored = 1000;
const std::string MCMCGRAutoNew::baseLogpostColName = "logpost_";
const std::string MCMCGRAutoNew::baseSequenceStateFilename = "SequenceStates";
// files for outputing samples
const std::string MCMCGRAutoNew::samplesLogFilename = "LogSamplesFromMCMCGRAuto.txt";

//interval between logging			
const size_t MCMCGRAutoNew::logInterval = 1000;
	
	

MCMCGRAutoNew::MCMCGRAutoNew(const std::vector < MCMCGRAuto::Diagnostic * >& diagPtrs,
							unsigned long int seed)
	: diagObjPtrs(diagPtrs), rgsl(NULL), time(0.0), burntinReachedState(0)
	
{
	
	try {
		
		rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);
	}
	catch (...) {
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (...) {} // catch and swallow
		throw;
	}
}


MCMCGRAutoNew::~MCMCGRAutoNew()
{
	try {
		if (NULL != rgsl) {
			gsl_rng_free(rgsl);
			rgsl = NULL;
		} 
	}
	catch (...) {} // catch and swallow
	
	
}

/* default rhatFlagCounterThreshold to number of diagnostic criteria, ie
all have to be 'true' for convergence to be deemed to have taken place */ 
PiecewiseConstantFunction MCMCGRAutoNew::_doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,  
							size_t minPoints,
							real minVol,
							int runID,
							const std::string& scalarsFileName,   
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	int rhatFlagCounterThreshold = diagObjPtrs.size();
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 
							runID,
							scalarsFileName,   
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41
}
	
	
/* rhatFlagCounterThreshold is how many of the scalar values must have
* diagnostic within limits for sampling to start?
* usually this would probably be the number
* of scalar values being used? */ 
PiecewiseConstantFunction MCMCGRAutoNew::_doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,  
							size_t minPoints,
							real minVol,
							int rhatFlagCounterThreshold, 	
							int runID,
							const std::string& scalarsFileName,   
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{	
	#ifdef MYDEBUG_LOG_SEQUENCESTATES
		loggingInChangeStates = TXT;
	#endif
	
	#ifdef TIMEINFO
		struct tm * timeinfo;
	#endif
	
	// make sure that everything is cleaned up
	clean();
	
	time = 0.0; // reset last runtime
	
	bool volChecking = (minVol > 0.0);
	
	// how many chains are to be run = number starting histograms
	size_t chains = histPtrsVec.size(); 
	
	if (chains < 2) {
		throw std::invalid_argument("Chains < 2");
	}
	
	size_t vecReserve = static_cast<size_t> (3 * samplesNeeded * thinout / chains);
	if (maxLoops/10 > static_cast<int>(vecReserve)) vecReserve = maxLoops/10; 
	if (maxLoops < static_cast<int>(vecReserve)) vecReserve = maxLoops;
		
	ostringstream stmH;
	stmH << runID; 
	std::string runIDstr = stmH.str(); // use for appending to output
	
	cout << "\nInitialising things using the current histogram state: " << endl;
	cout << "Number of diagnostics is : " << diagObjPtrs.size() << endl;
	cout << "rhatFlagCounterThreshold is : " << rhatFlagCounterThreshold << endl;
	cout << "samplesMonitored is : " << samplesMonitored << endl;
	
	// set up proposal distribution object - the stay split merge chain
	double probStay = 1.0E-6;
	UniformSSMProposal proposal(probStay);
	//gat41
	// set up prior distribution object
	//LogCatalanPrior logPrior;
	//LogTemperaturePrior logPrior(0.316474);

	
	// set up containers for the stuff we need pass to the MCMC engine
	vector < ChangeOfStateInformationAutoMCMC > infoObjs;
	
	vector<SPSnodeList> nodeLists(chains);
	Size_tVec numLeavesVec(chains);
	Size_tVec numCherriesVec(chains);
	vector <RealMappedSPnode > shadows(chains);
	vector < RealMappedSPnode::ListPtrs > 
										shadowLeavesVec (chains);
	vector < RealMappedSPnode::ListPtrs > 
										shadowCherriesVec (chains);

	vector<string> sequenceStateFilenames(chains, string(""));
	
	if ( logging == LOGANDGRAPHSAMPLES ) {
		cerr << "Sorry, LOGANDGRAPHSAMPLES not supported: changing to LOGSAMPLES" << endl;
		logging = LOGSAMPLES;
	}
	if ( logging == TXTANDGRAPH ) {
		cerr << "Sorry, TXTANDGRAPH not supported: changing to TXT" << endl;
		logging = TXT;
	}
	if ( (logging == LOGSAMPLES) ) {
		outputFileStart(samplesLogFilename);
	}
	
	std::string filenameGRs;
	{
		ostringstream oss;
		
		// bolt together a list of the diagnostics we are using
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
			oss << "_" << diagObjPtrs[d]->getScalarsName();
		}
		filenameGRs = oss.str();
	}
	
	// a name for a file of the leaves and logposteriors
	std::string filenameScalarsAndLogpost("");
	if ( !scalarsFileName.empty() ) {
		  filenameScalarsAndLogpost = getScalarFilename(
								"ScalarsLogPostTrace"+filenameGRs,
								runIDstr,
								scalarsFileName);
	}
	
	// filenames for diagnostics from each GR object
	std::vector < std::string > filenamesGRResults(diagObjPtrs.size(), "");
	if ( !scalarsFileName.empty() ) {
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
			filenamesGRResults[d] = 
				getScalarFilename(diagObjPtrs[d]->getGRDiagnosticsFilename(),
									runIDstr,
									scalarsFileName);
		}
	}
	bool keepLogs = false;
	#ifdef MYDEBUG_LOGS
		keepLogs = true;
		// filenames for working calculations from each GR object
		std::vector < std::string > filenamesGRWorkingCalcs;
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
			filenamesGRWorkingCalcs.push_back(
				getScalarFilename(diagObjPtrs[d]->getGRWorkingCalcsFilename(),
									runIDstr,
									scalarsFileName) );
		}
	#endif
				
	// one vector of leaves for each chain
	//leavesPtr = new std::vector < Size_tVec >(chains); 
	
	/* one vector of log posteriors for each chain
	 for which we need to keep loglikelihoods and log priors as well*/
	logpost = std::vector < RealVec >(chains);  
	for (size_t ci = 0; ci < chains; ++ci) logpost[ci].reserve(vecReserve);

	
	//get the diagnostic objects to initialise here
	for (size_t d = 0; d < diagObjPtrs.size(); ++d) 
				diagObjPtrs[d]->initialise(chains, vecReserve, keepLogs);

	#ifdef MYDEBUG_LOGS
		/* keep a vector of indicators for whether a state was sampled
		* (not a real, but easier to output it if we treat it like one) */
		sampledInd = RealVec(1, cxsc::real(0.0));
		sampledInd.reserve(vecReserve);
	#endif

	
	bool cancontinue = true;
	

	// this loop is just setting up containers of file names
	// and getting info from the starting histograms that is
	// needed to start the chains
	for (size_t ci = 0; ci < chains; ci++) {
		
		// do not comment these out
		if (loggingInChangeStates != NOLOG) {
			std::ostringstream stm1;
			stm1 << baseSequenceStateFilename << ci << ".txt";
			sequenceStateFilenames[ci] = stm1.str();
		
			outputFileStart(sequenceStateFilenames[ci]);
		}
		
		/* we only need to do this because we are doing a step-by-step change of the
		* histogram states 'from the outside', ie through this example:  we need to
		* collect the stuff the histogram's changeMCMCstate method needs to make one 
		* change.  */
	  
		// set up a container for the leaf children
		SPSnodePtrs leafVec;
		// set up a container for the subleaf children
		SPSnodePtrs cherryVec;

		size_t numLeaves = 0;
		size_t numCherries = 0;

		// fill the container with the leaf children
		histPtrsVec[ci]->getSubPaving()->getLeaves(leafVec);
		// fill the container with the subleaf children
		histPtrsVec[ci]->getSubPaving()->getSubLeaves(cherryVec);

		//check the cherries are all "legal", ie pass isSplittableNode and minVol
		if (!cherryVec.empty()) {
			SPSnodePtrsItr cit;
			for (cit = cherryVec.begin(); cit < cherryVec.end(); ++cit) {
				if  ((volChecking && 
				
					((((*cit)->getParent() == NULL) && ((*cit)->nodeRealVolume() < 2*minVol))
						|| 
					(*cit)->getParent()->nodeRealVolume() < 4*minVol))
					||
					!((*cit)->isSplittableNode(minPoints)))
				{
					cerr << "cherry " << (*cit)->getNodeName() << endl;
					cerr << "cherry count" << (*cit)->getCounter() << endl;
					cerr << "Minvol is " << minVol << endl;
					cerr << "2*Minvol is " << (2*minVol) << endl;
					cerr << "cherry vol is " << (*cit)->nodeRealVolume()  << endl;
					 
					throw std::logic_error(
					"\nIllegal state - cherries do not satisfy minPoints and  minVol for split");
				}
			}
		}
		
		numCherries = cherryVec.size();

		// check if node is still splittable
		if (!leafVec.empty()) {
		/* but only put into the container the leaves which, if split,
			 * would have at least minPoints data points associated with them
			 * or could split with one child getting all the points
			 * and which satisfy minVol for splittability */
			SPSnodePtrsItr lit;
			for (lit = leafVec.begin(); lit < leafVec.end(); ++lit) {
				if  ((!volChecking || 
				
					!((((*lit)->getParent() == NULL) && ((*lit)->nodeRealVolume() < 2*minVol))
						|| 
					(*lit)->getParent()->nodeRealVolume() < 4*minVol))
					&&
					((*lit)->isSplittableNode(minPoints)) )
				{
					// leaf can go into container
					nodeLists[ci].push_back(*lit);
					++numLeaves;
				}
			}
			
		}
		// no need to check on cherries - they can all go in
		if (numCherries > 0)
			 nodeLists[ci].insert(nodeLists[ci].end(), cherryVec.begin(),
											 cherryVec.end());
		if (nodeLists[ci].empty()) {
			 cancontinue = false;
			 break; // break out of the for loop
			 std::cout << "No changeable nodes given minPoints = "
							 << minPoints << " in histogram " << ci
							 << ". Sorry, aborting MCMC." << std::endl;
		}

		numLeavesVec[ci] = numLeaves;
		numCherriesVec[ci] = numCherries;

		/* Make a RMP out of the current histogram and run that
		 * alongside changes to the histogram; when sampling this is 
		 * copied into a PCF ...*/
		shadows[ci] = RealMappedSPnode((*histPtrsVec[ci]->getSubPaving()));
		RealMappedSPnode::Ptrs slVec;
		RealMappedSPnode::Ptrs scVec;
		shadows[ci].getLeaves(slVec);
		shadows[ci].getSubLeaves(scVec);
		
		shadowLeavesVec[ci].assign(slVec.begin(),slVec.end());
		shadowCherriesVec[ci].assign(scVec.begin(),scVec.end());
		
		// initialise things for the collection of data on leaves
		
		// one vector of leaves for each chain
		// record leaves for this first state
		size_t lastStateLeaves = shadowLeavesVec[ci].size();
		//leavesPtr->at(ci).push_back( lastStateLeaves );  
		
		
		//logposterior
		real loglik = histPtrsVec[ci]->getLogLikelihood();
		real logpi = logPrior(lastStateLeaves-1);
		real logposterior = loglik + logpi;
		logpost.at(ci).push_back(logposterior );

		/* use shadows information, and the current hist, to set up
		 * the info object for each chain */
		infoObjs.push_back(
			ChangeOfStateInformationAutoMCMC(
				logposterior,
				lastStateLeaves, shadowCherriesVec[ci].size(),
				*histPtrsVec[ci]) );
		
		
	} // end loop through chains setting up things to be able to start

	// get the diagnostic objects to initialise values for each chain
	for (size_t d = 0; d < diagObjPtrs.size(); ++d) 
					diagObjPtrs[d]->initialiseChainValues(infoObjs);

	bool goodLoop = cancontinue;

	if (cancontinue) cout << "\nAbout to do MCMC" << endl;

	/* set up some variables we need to keep track of states and sampling */
	int samplesSoFar = 0;
	size_t states = 1;  /* keep track of states in the chain = 1 so far,
					since state 1 is the initial shadow state */
	size_t checkedStates = 1;  // states for which convergence is checked 
	// variables for monitoring convergence

	/* a vector of GR flags, one for each GR object
	 * each one is indicator for whether we are burntin on diagnostics  scalar value */
	std::vector<int> rhatFlags(diagObjPtrs.size(), 0);
	//-----------------------

	int burntin = 0; // indicator for whether we consider ourselves burntin yet
	int rhatokay = 0; // indicator for whether all the required rhat criteria are currently met
	burntinReachedState = 0; // keep track of when we (first) reached burnin
	int rhatFlagCounter = 0;
		
										
	// counter to keep track of loops
	int loopCounter = 0;

	/* We also need a rmp to accumulate the samples
	 * Use the same paving box as the hist - use the box from the first one*/
	RealMappedSPnode* samplesRMP 
				= new RealMappedSPnode(histPtrsVec.at(0)->getRootBox());
	PiecewiseConstantFunction samplesAvg;
	
	gsl_rng * r = NULL;
	
	try {
		
		long unsigned int seed = 1 + gsl_rng_uniform_int(rgsl, gsl_rng_max(rgsl) );	
				
		r = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(r, seed);
		
		clock_t startMCMC = clock();	// timing
		
		while (goodLoop && (loopCounter < maxLoops) && (samplesSoFar < samplesNeeded)) 
		{
			#ifdef MYDEBUG_CALCS
				cout << "****** Change from state number " << states << " ******" << endl;
			#endif

			loopCounter++;
			
			// for each histogram and its shadow in turn, change the state
			/* 
			 * this is all a fudge - changeMCMCstate should just be a private
			 * method of the histograms but I think I made it public so that
			 * I could use it here in the example as a first step to being
			 * able to make all of this chain convergence stuff back into
			 * a method of the histograms themselves
			*/

			for (size_t ci = 0; ci < chains; ci++) {
			
				#ifdef MYDEBUG_CALCS
					cout << "--- chain index " << ci << " ---" << endl;
				#endif

			
				/* I refer to the current chain, indexed by ci, as 'the chain
				* in the comments inside this loop */
				
				// infoObj is updated for change in log posterior 
				histPtrsVec[ci]->changeMCMCState(nodeLists[ci],
					 numLeavesVec[ci], numCherriesVec[ci],
					 proposal, logPrior, minPoints,
					 minVol,
					 shadowLeavesVec[ci], shadowCherriesVec[ci],
					 infoObjs[ci], // updated
					 r, loggingInChangeStates,
					 sequenceStateFilenames[ci], states);

				#ifdef MYDEBUG_XTRA
					cout << "completed change of state" << endl;
				#endif
				 
				#ifdef FORCEFAILMCMCLOOP
					// for debugging - force a loop failure and see what happens to program
					if (runID == 2 && states == 5) goodLoop = false;
				#endif 

				if ((numLeavesVec[ci] == 0 && numCherriesVec[ci] == 0)) {
					throw std::runtime_error("No more leaves or cherries in MCMC");
				}
			
				// so assume all is okay if we have not just thrown an exception
				
				// store new logposterior
				/* Only if !burntin or if we are not normal sampling or if samplesSoFar > samplesMonitored*/
				if ( !burntin || (samplingType != NORMAL) || (samplesSoFar < samplesMonitored) ) {
			
					logpost.at(ci).push_back(
									infoObjs[ci].getCurrentLogPosterior() );
				}
				
				
				// update leaves for last histogram's state in the chain
				//size_t lastStateLeaves = infoObjs[ci].getCurrentLeaves();
				//assert (lastStateLeaves == shadowLeavesVec[ci].size());
				//assert (infoObjs[ci].getCurrentCherries() == shadowCherriesVec[ci].size());
				
				//leavesPtr->at(ci).push_back( lastStateLeaves );  
								
			} // end change state for each histogram in turn
			
			// increment number of states histograms have been through
			++states;
			//cout << "------------state " << states << endl;

			
			/* Only do any updates if burntin if we are not normal sampling or if samplesSoFar > samplesMonitored*/
			if ( !burntin || (samplingType != NORMAL) || (samplesSoFar < samplesMonitored) ) {
				
				// increment number of convergence check states histograms have been through
				++checkedStates;
				//cout << "------------state " << states << endl;

			
				// for each diagnostic, do the updates
				for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
						
					diagObjPtrs[d]->updateChainValuesInLoop(infoObjs);

				}	

				/* each chain now has a new state
				 * and info for leaves and logposteriors has been collected
				 * and the sample variance of the diagnostics for each chain 
				 * have been updated,
				 * so we can now work out the convergence diagnostics */


				/* for each diagnostic, get the rhatFlag value and change
				 * the flag counter as necessary 
				 * 
				 * and also check if any required diagnostic is not flagged*/
				bool requiredFlagged = true;
				for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
							
					int newRhatFlag = diagObjPtrs[d]->calcDiagnosticsForLoop();
				
					if (newRhatFlag && !rhatFlags[d]) {
						rhatFlags[d] = 1;
						++rhatFlagCounter;
						
					}
					else if (!newRhatFlag) {
						
						if ( rhatFlags[d]) {
							
							rhatFlags[d] = 0; // update the flag
							--rhatFlagCounter;
							// if this one required, we don't have all the required ones
						}
						if (diagObjPtrs[d]->isRequired()) requiredFlagged = false;
					}
				} // end loop through diagnostics

				if ( !burntin && (rhatFlagCounter >= rhatFlagCounterThreshold) 
							&& requiredFlagged ) {
					rhatokay = 1;
					burntin = 1; 
					// only set burntinReachstate if this is the first time
					if (!burntinReachedState) burntinReachedState = states;
					// except if we are doing cautious sampling
					if (samplingType == CAUTIOUS) burntinReachedState = states;
					
					
					#ifdef MYDEBUG
						// if we have not been burntin, give a message
						 cout << "Burnin convergence test satisfied at state " 
							  << states << endl;
					#endif
				}
				
				// but it may be that we were burntin and no longer are
				else if ( burntin && ((rhatFlagCounter < rhatFlagCounterThreshold) 
									|| !requiredFlagged )) {
					
					// for all sampling types other than NORMAL, say we are not burntin
					if (samplingType != NORMAL) burntin = 0; 
					
					// and if we are doing CAUTIOUS sampling, throw away the existing samples
					if (samplingType == CAUTIOUS) {
						burntinReachedState = 0;
						
						ivector box = samplesRMP->getBox();
						delete samplesRMP; // get rid of the old samples collator
						samplesRMP = NULL;
						samplesRMP = new RealMappedSPnode(box); // and take a new one
						
						samplesSoFar = 0;
						
						// want to change all the 1's in sampledInd so far to 0s
						cxsc::real newVal(0.0);
						
						#ifdef MYDEBUG_LOGS
							std::replace_if (sampledInd.begin(), sampledInd.end(), 
								std::bind2nd(greater< cxsc::real >(),newVal), newVal);
						#endif
					
						// restart the log file if we are logging
						if ( (logging == LOGSAMPLES) ) {
							outputFileStart(samplesLogFilename);
						}
					}	
					
					if (rhatokay) {
						rhatokay = 0;
					
						#ifdef MYDEBUG
						if (!burntin)
							cout << "Burnin convergence test now NOT satisfied at state " 
								  << states << endl;
							
						#endif
					}
				}
				/*output the leaf and logpost traces as we go */
				if (!filenameScalarsAndLogpost.empty() && checkedStates%logInterval == 0) {
					assert (checkedStates >= logInterval);
					int precData = 10;
					outputLeafLogPost(filenameScalarsAndLogpost,
							checkedStates - logInterval, 
							precData);
					
				}
			} // end section that conditionally checks convergence diagnostics	
			
			/* take samples if we are burntin and this is a sampling point according to 
			 * the thinout specified 
			 * note - we will only be in the loop at all if we still need more samples*/
			if (burntin && (( states - burntinReachedState )%thinout == 0)) {
				
				#ifdef MYDEBUG_LOGS
					sampledInd.push_back (cxsc::real(1.0)); 

				#endif

				// take one sample from each chain until we have enough samples
				// and increment samplesSoFar for each one taken
				for (vector < RealMappedSPnode >::iterator rit = shadows.begin(); 
						(rit < shadows.end() && samplesSoFar < samplesNeeded);
						++rit) {
					
					// add the state (it is not normalised at this stage - do all that later
					samplesRMP->operator+=(*rit);

					samplesSoFar++;
					#ifdef MYDEBUG_XTRA
						cout << "samples so far: " << samplesSoFar << endl;
					#endif
					
					if ( (logging == LOGSAMPLES) ) {
						PiecewiseConstantFunction p(*rit);
						p.outputLog(samplesLogFilename, samplesSoFar);
					}
				}
			} // finished taking samples for this loop
			else {
				#ifdef MYDEBUG_LOGS
					sampledInd.push_back (cxsc::real(0.0)); 
				#endif
			}
			
			



			// back into loop
			#if !defined(MYDEBUG_CALCS)
				#ifdef MYDEBUG
					// output a line every now and again so that we know it's still alive
					if (loopCounter%10000 == 0) {
						cout << "\n...Still going: completed change in state number " << states;
						
						if ( !burntin || (samplingType != NORMAL) || (samplesSoFar < samplesMonitored) ) {
			
							for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
								cout << "\n\trhat  for " << diagObjPtrs[d]->getScalarsName()
									<< " = " << diagObjPtrs[d]->getCurrentRhatValue();
							}
							
						}
						cout << "  ...\n" << endl;
						if ( burntin && (samplingType == NORMAL) ) 	cout << "\tburnt-in";
						#ifdef TIMEINFO
							time_t rawtime;
			  				std::time ( &rawtime );
							timeinfo = localtime ( &rawtime );
							cout << "\ttime: " << asctime (timeinfo) << endl;
						#endif
						
					}
				#endif
			#endif
		}    // finished while loop - either loop failed or reached maxLoops or have all our samples

		cancontinue = goodLoop;

		// stop recording time here
		clock_t endMCMC = clock();	
		time = ((static_cast<double>(endMCMC - startMCMC)) / CLOCKS_PER_SEC);
			
		#ifdef MYDEBUG
			cout << "****** finished all loops, states counter is = " << states << " ******" << endl;
			cout << "Computing time for core MCMC process: " << time << " s."<< endl;
		
		#endif

		cout << "\nnumber of samples collected is = " << samplesSoFar << endl;

		// free the random number generator
		try {
			if (r != NULL) gsl_rng_free (r);
			r = NULL;
		}
		catch (...) {} // catch and swallow
	}
	catch (...) {
		
		/*output anything more we have for leaves and log posteriors
		* last scheduled output will have been at last integer multiple 
		* of logInterval loops*/
		if (!filenameScalarsAndLogpost.empty() && (checkedStates % logInterval) ) {
			cout << "\nfinal output of " 
				<< (checkedStates - (checkedStates/logInterval) * logInterval) 
				<< " states" << endl;
			int precData = 5;
			outputLeafLogPost(filenameScalarsAndLogpost,
						(checkedStates/logInterval) * logInterval,
						precData);
		}
		
		// free the random number generator
	
		try {
			if (r != NULL) gsl_rng_free (r);
			r = NULL;
		}
		catch (...) {} // catch and swallow
		try {
			if (samplesRMP != NULL) delete samplesRMP;
			samplesRMP = NULL;
		}
		catch (...) {} // catch and swallow
		
		// clean all the diagnostics
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
			diagObjPtrs[d]->clean();
		}
		
		//clean me
		clean();
		
		throw;
	}
		

	cout << cxsc::RestoreOpt; // undo changes we made to cout printing for cxsc values

	// output stuff whether we succeeded in getting all our samples or not
	{
	
	// output the convergence diagnostics for each diagnostic
		if (!scalarsFileName.empty() ) {
			int precData = 5;
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
				
				#ifdef MYDEBUG_LOGS
					diagObjPtrs[d]->outputResults(filenamesGRResults[d],
											sampledInd, precData);
				#else
					diagObjPtrs[d]->outputResults(filenamesGRResults[d],
												precData);
				#endif
			}
			
		} 
		
		//output leaves and log posteriors together
		/*last scheduled output will have been at last integer multiple of logInterval loops*/
		if (!filenameScalarsAndLogpost.empty() && (checkedStates % logInterval) ) {
			int precData = 5;
			outputLeafLogPost(filenameScalarsAndLogpost,
						(checkedStates/logInterval) * logInterval,
						precData);
						
		}
		
		#ifdef MYDEBUG_LOGS
		{
			// output once for every diagnostic
			int precData = 5;
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
				
					diagObjPtrs[d]->outputCalculations(
												filenamesGRWorkingCalcs[d],
												precData);
			}
		}
		#endif

	} // end of output of everything we want to output whether we have samples or not
		

	/* is all okay with the loop
	 * and we have all our samples */
	if (cancontinue && (samplesSoFar >= samplesNeeded) ) {
		
		// average the rmp over the samples and then normalise
		//samplesRMP->operator/=(cxsc::real(1.0*samplesSoFar));
		samplesAvg = PiecewiseConstantFunction(*samplesRMP);
		samplesAvg.normalise();
		cout << "\nafter normalising, samplesAvg.getTotalIntegral = " << samplesAvg.getTotalIntegral() << endl;
		
		cout << "\n\nFinished MCMC successfully" << endl;
		cout << "Burnt in at " << burntinReachedState << endl;

	
		if (!scalarsFileName.empty() ) {
			cout << "\nFor diagnostics, check output files:" << endl;
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
											
				cout << "\t" << filenamesGRResults[d] << endl;
			}
		
			cout << "and for log posteriors\t" << filenameScalarsAndLogpost << endl;
		}
		if ( (logging == LOGSAMPLES) ) {
			cout << "and for log of samples\t" << samplesLogFilename << endl;
		}

		#ifdef MYDEBUG_LOGS
			cout << "and for working calculations for diagnostics:" << endl;
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {

			cout << "\t" << filenamesGRWorkingCalcs[d] << endl;
		}
		#endif
	
		cout << endl;
	}
	
	
	/* since I throw an exception in the while loop if it is not a good loop,
	 *  really the only reason for failing here is that we did not get the right 
	 * number of samples, but might as well leave it like this - belt & braces*/ 
	if (!cancontinue || (samplesSoFar < samplesNeeded) ) {
		cout << "\nMCMC not successful" << endl;
		if (!scalarsFileName.empty() ) {
			cout << "Output files will not be complete - delete or ignore:" << endl;
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {		
				cout << "\t" << filenamesGRResults[d] << endl;
			}
			cout << "\t" << filenameScalarsAndLogpost << endl;
		}
		#ifdef MYDEBUG_LOGS
			for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
				cout << "\t" << filenamesGRWorkingCalcs[d] << endl;
			}
		#endif
		if ( (logging == LOGSAMPLES) ) {
			cout << "\t" << samplesLogFilename << endl;
		}
		cout << endl;
	}
		
	clean();
	try {
		if (samplesRMP != NULL) delete samplesRMP;
		samplesRMP = NULL;
		
		// clean all the diagnostics
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
			diagObjPtrs[d]->clean();
		}

	}
	catch (...) {} // catch and swallow any exception on cleanup

	if (!cancontinue) {
		throw std::runtime_error("MCMC failed");
	}
	if (samplesSoFar < samplesNeeded) {
		// we have not been able to get the required samples - need to give up
		throw std::runtime_error("Did not get required number of samples");
	}

	return samplesAvg;
	
} 


int MCMCGRAutoNew::getDefaultSamplesMonitored() const
{
	return defaultSamplesMonitored;
}

void MCMCGRAutoNew::outputLeafLogPost(
						const std::string& filenameScalarsAndLogpost,
						int precData) const 
{
	size_t startPos = 0;
	
	outputLeafLogPost(
						filenameScalarsAndLogpost,
						startPos, 
						precData);
	
						
}

void MCMCGRAutoNew::outputLeafLogPost(
						const std::string& filenameScalarsAndLogpost,
						size_t startPos, 
						int precData) const
{	
	//get the names for each diag's scalar
	std::vector < std::string> scalarNames;
	for (size_t d = 0; d < diagObjPtrs.size(); ++d) 
			scalarNames.push_back(diagObjPtrs[d]->getScalarsName());
	
	size_t chains = logpost.size();
	
	std::vector < std::string> colNames;
			
	std::vector < const RealVec* > dataPtrs;
	for (size_t ci = 0; ci < chains; ++ci) {
		
		for (size_t d = 0; d < diagObjPtrs.size(); ++d) {
		
			{
				ostringstream oss;
				oss << scalarNames[d] << "_" << ci;
				colNames.push_back(oss.str());
			}
		
			dataPtrs.push_back(&(diagObjPtrs[d]->getScalarsRef().at(ci)));
			
		}
		
		// then the log post
		{
			ostringstream oss;
			oss << baseLogpostColName << "_" << ci;
			colNames.push_back(oss.str());
		}
		dataPtrs.push_back(&(logpost.at(ci)));
	}
	
	
	if (startPos) { // starting from pos > 0
		
		outputToFileVertical(dataPtrs, 
							filenameScalarsAndLogpost,
							startPos,
							precData);
		
	}
	else { // starting from pos 0
		//output leaves and log posteriors together
		outputToFileVertical(dataPtrs, 
							colNames,
							filenameScalarsAndLogpost, precData);
	}
}

size_t MCMCGRAutoNew::getMaxFlagThreshold() const
{
	return diagObjPtrs.size();
}


double MCMCGRAutoNew::getLastRunTime() const
{
	return time;
}

size_t MCMCGRAutoNew::getFirstBurninTime() const
{
	return burntinReachedState;
}


std::string MCMCGRAutoNew::getScalarFilename(const std::string& scalarType,
								const std::string& runIDstr,
								const std::string& scalarsFileName)
{
	/* Search for the last '/' in the scalarsFilenName, break
	 * it there, join again with scalarType and runIDstr in the middle*/

	std::string retValue;
	std::string rID = "_r" + runIDstr;

	std::string file = scalarsFileName; 
	std::string path("");
	size_t found = file.find_last_of("/");
	if (found!=string::npos) {
		path = file.substr(0,found+1);
		file = file.substr(found+1);
	}
	// file stays as full scalarsFileName
		
	found = file.find_last_of(".");
	if (found!=string::npos) {
		std::string prefix = file.substr(0,found);
		std::string suffix = file.substr(found);
		file = prefix + rID + suffix;
	}
	else file = file + rID;
		
	retValue  = path + scalarType + file;
		
	return retValue;
}

std::string MCMCGRAutoNew::getPath(const std::string& scalarsFileName)
{
	std::string file = scalarsFileName; 
	std::string path("");
	size_t found = file.find_last_of("/");
	if (found!=string::npos) {
		path = file.substr(0,found+1);
	}
		
	return path;
}


void MCMCGRAutoNew::clean() const
{
	
	std::vector < RealVec >().swap( logpost );
	RealVec().swap( sampledInd );
	
}

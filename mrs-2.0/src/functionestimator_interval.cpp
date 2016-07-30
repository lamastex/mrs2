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
\brief FunctionEstimatorInterval definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "functionestimator_interval.hpp"
#include "intervalexpander_estimate.hpp"
#include "realmappedspnode.hpp"
#include "real_estimate.hpp"
#include "interval_estimate.hpp"
#include "intervalmappedspnode_measurers.hpp"
#include "toolz.hpp"
#include "sptools.hpp"



#include "subpaving_exception.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <utility>     // to use stl::pair
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <cassert>

#include <gsl/gsl_randist.h>

//#define MYDEBUGBF
//#define MYDEBUGBFMIN
//#define MYDEBUGPQ
//#define MYDEBUGPQMIN
//#define MYDEBUGPM
//#define MYDEBUGPQGAINCALCS

#ifdef NDEBUG
	#undef MYDEBUGBF
	#undef MYDEBUGBFMIN
	#undef MYDEBUGPQ
	#undef MYDEBUGPQMIN
	#undef MYDEBUGPM
	#undef MYDEBUGPQGAINCALCS


#endif

using namespace subpavings;
using namespace std;

#if defined (MYDEBUGPQMIN) || defined (MYDEBUGBFMIN)
	#include <ctime>
#else
	#define SHOW_TIMEDOTS
	#define TIMEDOTS_NUM 10000 // number of loops for timedots
#endif


// -------------------implementation of FunctionEstimatorInterval class --------------


// ----------- public methods


// initialised constructor 
FunctionEstimatorInterval::FunctionEstimatorInterval(
										const ivector& v,
										const MappedFobj& f,
										int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"FunctionEstimatorInterval::FunctionEstimatorInterval(const ivector&, const MappedFobj&, int lab)");
		}
        rootPaving = new IntervalMappedSPnode(v);
		
		IntervalEstimator estimator(fobj);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

FunctionEstimatorInterval::FunctionEstimatorInterval(const SPnode& spn, 
											const MappedFobj& f, 
											int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check spn has box
		if (spn.isEmpty()) {
			throw subpavings::NoBox_Error(
			"FunctionEstimatorInterval::FunctionEstimatorInterval(const SPnode&, MappedFobj&, int lab");
		}
        rootPaving = new IntervalMappedSPnode(spn);
		
		IntervalEstimator estimator(fobj);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// copy constructor*/
FunctionEstimatorInterval::FunctionEstimatorInterval(
								const FunctionEstimatorInterval& other)
        : rootPaving(NULL), fobj(other.fobj), label(other.label)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new IntervalMappedSPnode(*(other.getSubPaving()));
			
		} // else subpaving is NULL which should be impossible
		else {
			throw NullSubpavingPointer_Error(
				"FunctionEstimatorInterval::FunctionEstimatorInterval(const FunctionEstimatorInterval&)");
		}
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}



//Destructor
FunctionEstimatorInterval::~FunctionEstimatorInterval()
{
	try {
		delete rootPaving;
		rootPaving = NULL;

	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in FunctionEstimatorInterval destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

/*copy assignment operator private and not implemented:
Cannot do assignment because fobj is a reference.*/


const MappedFobj& FunctionEstimatorInterval::getFobjReference() const
{return fobj;}

int FunctionEstimatorInterval::getLabel() const
{return label;}

void FunctionEstimatorInterval::setLabel(int lab)
{label = lab;}

// get whether this has a subpaving.
bool FunctionEstimatorInterval::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

cxsc::ivector FunctionEstimatorInterval::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionEstimatorInterval::getRootBox()");
	}
	return getSubPaving()->getBox();
}

cxsc::real FunctionEstimatorInterval::getRootRangeDiameter() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionEstimatorInterval::getRootRangeDiameter()");
	}
	return getSubPaving()->getRangeDiameter();
}

int FunctionEstimatorInterval::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

cxsc::real FunctionEstimatorInterval::getDomainVolume() const
{
	real retValue(0.0);
	if (hasSubPaving()) {
		retValue = getSubPaving()->nodeRealVolume();
	}
	return retValue;
}
	

// Gets number of leaf nodes in the root paving.
size_t FunctionEstimatorInterval::getRootLeaves() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("FunctionEstimatorInterval::getRootLeaves()");
		
	}
	return getSubPaving()->getNumberLeaves();
}


// Get a string of the leaf node levels.
std::string FunctionEstimatorInterval::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}


subpavings::PiecewiseConstantFunction 
			FunctionEstimatorInterval::makePiecewiseConstantFunction() const
{
	/* make a rmspnode as a copy of this's subpaving and then put the 
	ranges on as the mid-image of the boxes */
	RealMappedSPnode tmp(*getSubPaving());
	
	RealEstimator estimator(fobj);
	tmp.acceptSPValueVisitor(estimator);
	return PiecewiseConstantFunction(tmp, getLabel());
}

void FunctionEstimatorInterval::bruteForceEstimate(cxsc::real tolerance)
{
	string errorMsg(
		"FunctionEstimatorInterval::bruteForceEstimate(cxsc::real)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	IntervalExpanderEstimator estimator(fobj, tolerance);

	getSubPaving()->acceptSPExpandVisitor(estimator);

}

void FunctionEstimatorInterval::breadthFirstBruteForceEstimate(
					size_t maxLeaves,
					long unsigned int seed)
{
	SplittableCheck nodeChecker;
	breadthFirstBruteForceEstimate(nodeChecker, maxLeaves, seed);
	
}
	
void FunctionEstimatorInterval::breadthFirstBruteForceEstimate(
				const FEIEvalObj& fe,
				long unsigned int seed)
{
	SplittableCheck nodeChecker;
	breadthFirstBruteForceEstimate(nodeChecker, fe, seed);
	
}

void FunctionEstimatorInterval::breadthFirstBruteForceEstimate(
					const SPCheckVisitor& nodeChecker,
					size_t maxLeaves,
					long unsigned int seed)
{
	string errorMsg(
		"FunctionEstimatorInterval::breadthFirstBruteForceEstimate(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	gsl_rng * r = NULL;
	
	try {
		
		NodeQueueT nq;
		
		_setupBreadthFirstQueue(nq, nodeChecker);
					
		r = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(r, seed);
		
		// actual number of leaves in subpaving
		size_t numLeaves = getRootLeaves();
			
		#ifdef MYDEBUGBFMIN
			clock_t start = clock();
		#endif
	
		IntervalEstimator estimator(fobj);
		
		/* split until there are no more leaves needing splitting 
		or there are enough leaves */
		while (!nq.empty() && numLeaves < maxLeaves) {
			
			_breadthFirstQueueLoop(
					nodeChecker, nq, estimator, r);
			
			/* actual number of leaves has increased by one 
			irrespective of whether any have actually gone into the queue*/
			numLeaves++;
			
			#ifdef MYDEBUGBFMIN
				if (numLeaves%TIMEDOTS_NUM == 0) {
					double time = static_cast<double>(clock()-start)/CLOCKS_PER_SEC;
					cout << numLeaves << "\t" << time << "\t" << getTotalAreaOfIntervalBand() << endl;
					start = clock();
				}
			#endif
			#ifdef MYDEBUGBF
					cout << "current estimator is " << endl;
					outputToStreamTabs(cout, 10);
					cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << "\n\n" << endl;
			#endif
		}// loop
		
		try {
			if (NULL != r) {
				gsl_rng_free(r);
				r = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
	
		#ifdef MYDEBUGBF
			cout << "\nEnd of queue " << endl;
			if (nq.empty()) cout << "ran out of nodes to expand" << endl;
			else cout << "number of leaves is " << getRootLeaves() << endl;

		#endif
	}
	
    catch (exception const& e) {
        try {
			if (NULL != r) {
				gsl_rng_free(r);
				r = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }

}


void FunctionEstimatorInterval::breadthFirstBruteForceEstimate(
					const SPCheckVisitor& nodeChecker,
					const FEIEvalObj& fe,
					long unsigned int seed)
{
	string errorMsg(
		"FunctionEstimatorInterval::depthFirstBruteForceEstimate(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	gsl_rng * r = NULL;
	
	try {
		
		NodeQueueT nq;
		
		_setupBreadthFirstQueue(nq, nodeChecker);
					
		r = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(r, seed);
		
		// actual number of leaves in subpaving
		size_t numLeaves = getRootLeaves();
			
		#ifdef MYDEBUGBFMIN
			clock_t start = clock();
		#endif
	
		IntervalEstimator estimator(fobj);
	
		/* split until fe satisfied */
		while (!nq.empty() && !(fe(*this)) ) {
			
			_breadthFirstQueueLoop(
					nodeChecker, nq, estimator, r);
			
			/* actual number of leaves has increased by one 
			irrespective of whether any have actually gone into the queue*/
			numLeaves++;
			
			#ifdef MYDEBUGBFMIN
				if (numLeaves%TIMEDOTS_NUM == 0) {
					double time = static_cast<double>(clock()-start)/CLOCKS_PER_SEC;
					cout << numLeaves << "\t" << time << "\t" << getTotalAreaOfIntervalBand() << endl;
					start = clock();
				}
			#endif
			#ifdef MYDEBUGBF
					cout << "current estimator is " << endl;
					outputToStreamTabs(cout, 10);
					cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << "\n\n" << endl;
			#endif
		}// loop
		
		try {
			if (NULL != r) {
				gsl_rng_free(r);
				r = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
	
		#ifdef MYDEBUGBF
			cout << "\nEnd of queue " << endl;
			if (nq.empty()) cout << "ran out of nodes to expand" << endl;
			else cout << "number of leaves is " << getRootLeaves() << endl;

		#endif
	}
	
    catch (exception const& e) {
        try {
			if (NULL != r) {
				gsl_rng_free(r);
				r = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }

}

void FunctionEstimatorInterval::hullPropagation()
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorInterval::hullPropagation(...)");
	}
	getSubPaving()->hullPropagation();
}

bool FunctionEstimatorInterval::prioritySplitOnGain(
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	IntervalBandAreaGainOnSplitMeasurer measure(fobj);
	return prioritySplit(measure, fe, logging, seed);
}

bool FunctionEstimatorInterval::prioritySplitOnGain(
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	IntervalBandAreaGainOnSplitMeasurer measure(fobj);
	return prioritySplit(measure, maxLeaves, logging, seed);
}

bool FunctionEstimatorInterval::prioritySplit(
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	ReimannDiffMeasurer measure;
	return prioritySplit(measure, fe, logging, seed);
}


bool FunctionEstimatorInterval::prioritySplit(
						const IntervalMappedSPnode::Measurer& measure,
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
        string errorMsg("FunctionEstimatorInterval::prioritySplit(...)");
		
		if (!hasSubPaving()) {
			throw NullSubpavingPointer_Error(errorMsg);
		}

		// set up a random number generator for uniform rvs
        rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);

		int i = 0;
		std::string baseFileName("pqOutput");
		std::string s = "";
		if (logging == TXT) {
			// pass to log output to keep track of splits
			s = getUniqueFilename(baseFileName);
			// Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the estimate
			outputLog(s, i);
			
		}
		
		
		std::string baseFileNameStates("pqState");
		if (logging == LOGSAMPLES) {
			ostringstream oss;
			oss << baseFileNameStates << i << ".txt";
			outputToTxtTabs(oss.str());
		}
		
		i++;
				
		PriorityQueueT pq;
		
		pq = _setupPrioritySplitQueue(pq, measure);
		
		bool canContinue = !pq.empty(); // used after loop to stop final logging etc
		
		if(!canContinue) {
			std::cerr << "No splittable leaves to split - aborting" << std::endl;
		}
		
		IntervalEstimator estimator(fobj);
		
		size_t numLeaves = getRootLeaves();
		
		#ifdef MYDEBUGPQMIN
			clock_t start = clock();
		#endif


		// split until the FEIEvalObj fe () operator returns true
		/* check there are nodes in queue at the top of loop rather
		 * than the bottom to ensure that we will stop if fe is 
		 * satisfied before we check queue has nodes in,
		 * ie only return false if fe is not satisfied at end of operation*/ 
		while (canContinue && !fe(*this)) {
			
			canContinue = _prioritySplitQueueLoop(
					measure, pq, estimator, rgsl);
			
			if (logging == TXT) { // logs even if something has gone wrong
				// To add current state of estimator to log file
				outputLog(s, i);
			}
			if (logging == LOGSAMPLES) {
				ostringstream oss;
				oss << baseFileNameStates << i << ".txt";
				outputToTxtTabs(oss.str());
			}
			
			i++;
			
			if (canContinue) numLeaves++;
			
			#ifdef SHOW_TIMEDOTS
			if (numLeaves%TIMEDOTS_NUM == 0) {
					cout << "." << flush;
			}
			#endif
			
			#ifdef MYDEBUGPQ
				cout << "\nat end of loop, number of leaf nodes in subpaving is " << (getRootLeaves()) << endl;
				cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << endl;
					
			#endif
			#ifdef MYDEBUGPQMIN
				
				if (numLeaves%TIMEDOTS_NUM == 0) {
					double time = static_cast<double>(clock()-start)/CLOCKS_PER_SEC;
					cout << numLeaves << "\t" << time << "\t" << getTotalAreaOfIntervalBand() << endl;
					start = clock();
				}
			#endif
			#ifdef MYDEBUGPQ
				cout << "current estimator is " << endl;
				outputToStreamTabs(cout, 10);
				cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << "\n\n" << endl;
				
			#endif
		}// loop
		
		#ifdef SHOW_TIMEDOTS
			if (numLeaves > TIMEDOTS_NUM-1) cout << endl;
		
		#endif
	

		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		
		return canContinue;
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}


bool FunctionEstimatorInterval::prioritySplit(
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	ReimannDiffMeasurer measure;
	return prioritySplit(measure, maxLeaves, logging, seed);
}

 //  using maxLeaves
bool FunctionEstimatorInterval::prioritySplit(
						const IntervalMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	
	gsl_rng * rgsl = NULL;

    try {
        string errorMsg("FunctionEstimatorInterval::prioritySplit(...)");
		
		if (!hasSubPaving()) {
			throw NullSubpavingPointer_Error(errorMsg);
		}

		// set up a random number generator for uniform rvs
        rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);

		
		int i = 0;
		
		std::string baseFileName("pqOutput");
		std::string s = "";
		if (logging == TXT) {
			// pass to log output to keep track of splits
			s = getUniqueFilename(baseFileName);
			// Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the estimate
			outputLog(s, i);
			
		}
		
		std::string baseFileNameStates("pqState");
		if (logging == LOGSAMPLES) {
			ostringstream oss;
			oss << baseFileNameStates << i << ".txt";
			outputToTxtTabs(oss.str());
		}
		
		i++;
	
		size_t numLeaves = getRootLeaves();
		
		PriorityQueueT pq;
		
		pq = _setupPrioritySplitQueue(pq, measure);
		
		bool canContinue = !pq.empty(); // used after loop to stop final logging etc
		
		if(!canContinue) {
			std::cerr << "No splittable leaves to split - aborting" << std::endl;
		}
		
		IntervalEstimator estimator(fobj);
		
		#ifdef MYDEBUGPQMIN
			clock_t start = clock();
		#endif


		// split until we have enough leaves
		/* check there are nodes in queue at the top of loop rather
		 * than the bottom to ensure that we will stop if enough leaves 
		 * before we check queue has nodes in,
		 * ie only return false if max leaves not achieved end of operation*/ 
		while (canContinue && numLeaves < maxLeaves) {
			
			canContinue = _prioritySplitQueueLoop(
					measure, pq, estimator, rgsl);
			
			if (logging == TXT) { // logs even if something has gone wrong
				// To add current state of estimator to log file
				outputLog(s, i);
			}
			if (logging == LOGSAMPLES) {
				ostringstream oss;
				oss << baseFileNameStates << i << ".txt";
				outputToTxtTabs(oss.str());
			}
			i++;
						
			// number of leaves will have increased by 1 for split
			if (canContinue) numLeaves++;
			
			
			#ifdef SHOW_TIMEDOTS
			if (numLeaves%TIMEDOTS_NUM == 0) {
					cout << "." << flush;
			}
			#endif
			#ifdef MYDEBUGPQ
				cout << "\nat end of loop, number of leaf nodes in subpaving is " << (getRootLeaves()) << endl;
				cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << endl;
					
			#endif
			#ifdef MYDEBUGPQMIN
				if (numLeaves%TIMEDOTS_NUM == 0) {
					double time = static_cast<double>(clock()-start)/CLOCKS_PER_SEC;
					cout << numLeaves << "\t" << time << "\t" << getTotalAreaOfIntervalBand() << endl;
					start = clock();
				}
			#endif
			#ifdef MYDEBUGPQ
				cout << "current estimator is " << endl;
				outputToStreamTabs(cout, 10);
				cout << "getTotalAreaOfIntervalBand() = " << (getTotalAreaOfIntervalBand()) << "\n\n" << endl;
				
			#endif
		}// loop
		#ifdef SHOW_TIMEDOTS
			cout << endl;
		#endif

		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		
		return canContinue;
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}

bool FunctionEstimatorInterval::priorityMergeOnLoss(
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
    IntervalBandAreaLossOnMergeMeasurer measure;
	return priorityMerge(measure, fe, logging, seed);
}

bool FunctionEstimatorInterval::priorityMergeOnLoss(
						size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
    IntervalBandAreaLossOnMergeMeasurer measure;
	return priorityMerge(measure, minLeaves, logging, seed);
}

bool FunctionEstimatorInterval::priorityMerge(
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
    ReimannDiffMeasurer measure;
	return priorityMerge(measure, fe, logging, seed);
}



bool FunctionEstimatorInterval::priorityMerge(
						const IntervalMappedSPnode::Measurer& measure,
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
        string errorMsg("FunctionEstimatorInterval::priorityMerge(...)");
		if (!hasSubPaving()) {
			throw NullSubpavingPointer_Error(errorMsg);
		}
		
		if (getSubPaving()->isLeaf()) {
			throw UnfulfillableRequest_Error(errorMsg);
		}
		
		// set up a random number generator for uniform rvs
        rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);

		int i = 0;
		
		std::string baseFileName("pqMergeOutput");
		std::string s = "";
		if (logging == TXT) {
			// pass to log output to keep track of splits
			s = getUniqueFilename(baseFileName);
			// Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the estimate
			outputLog(s, i);
			
		}
		
		std::string baseFileNameStates("pqMergeState");
		if (logging == LOGSAMPLES) {
			ostringstream oss;
			oss << baseFileNameStates << i << ".txt";
			outputToTxtTabs(oss.str());
		}
		
		
		PriorityQueueT pq;
		
		pq = _setupPriorityMergeQueue(pq, measure);
		
		size_t numLeaves = getRootLeaves(); // true number of leaves
		
		bool canContinue = !pq.empty();
		
		if (canContinue) {
			
			bool bigEnough = true;
			
			// merge until the FEIEvalObj fe () operator returns true
			while (bigEnough && !fe(*this)) {
				
				bigEnough = _priorityMergeQueueLoop(measure, pq, rgsl);
							
				if (logging == TXT) { // logs even if something has gone wrong
					// To add current state of estimator to log file
					outputLog(s, i);
				}
				if (logging == LOGSAMPLES) {
					ostringstream oss;
					oss << baseFileNameStates << i << ".txt";
					outputToTxtTabs(oss.str());
				}
				i++;
		
				// number of actual leaves will have reduced by 1 
				if (bigEnough) numLeaves--;
				
			} // loop
		}
		else {
			std::cerr << "No mergable cherries - aborting" << std::endl;
		}
		
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		
		return canContinue;  // true unless we could not even start merge
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}


bool FunctionEstimatorInterval::priorityMerge(
						size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
    ReimannDiffMeasurer measure;
	return priorityMerge(measure, minLeaves, logging, seed);
}



bool FunctionEstimatorInterval::priorityMerge(
						const IntervalMappedSPnode::Measurer& measure,
						size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
        string errorMsg("FunctionEstimatorInterval::priorityMerge(...)");
		if (!hasSubPaving()) {
			throw NullSubpavingPointer_Error(errorMsg);
		}
		
		if (getSubPaving()->isLeaf()) {
			throw UnfulfillableRequest_Error(errorMsg);
		}
		
		// set up a random number generator for uniform rvs
        rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);

		int i = 0;
		
		std::string baseFileName("pqMergeOutput");
		std::string s = "";
		if (logging == TXT) {
			// pass to log output to keep track of splits
			s = getUniqueFilename(baseFileName);
			// Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the estimate
			outputLog(s, i);
			
		}
		
		std::string baseFileNameStates("pqMergeState");
		if (logging == LOGSAMPLES) {
			ostringstream oss;
			oss << baseFileNameStates << i << ".txt";
			outputToTxtTabs(oss.str());
		}
		
		i++;
		
		PriorityQueueT pq;
		
		pq = _setupPriorityMergeQueue(pq, measure);
		
		size_t numLeaves = getRootLeaves(); // true number of leaves
		
		bool canContinue = !pq.empty();
		
		if (canContinue) {
			
			bool bigEnough = true;
			
			// merge until the we have at most minLeaves
			while (bigEnough && numLeaves > minLeaves) {
				
				bigEnough = _priorityMergeQueueLoop(measure, pq, rgsl);
							
				if (logging == TXT) {
					// To add current state of histogram to log file
					outputLog(s, i);
					
				}
				if (logging == LOGSAMPLES) {
					ostringstream oss;
					oss << baseFileNameStates << i << ".txt";
					outputToTxtTabs(oss.str());
				}
				i++;
				
				// number of actual leaves will have reduced by 1 
				if (bigEnough) numLeaves--;
				
			} // loop
		}
		else {
			std::cerr << "No mergable cherries - aborting" << std::endl;
		}
		
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		
		return canContinue;  // true unless we could not even start merge
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}




//splits histogram according to string instruction
//returns true if some splitting was achieved
bool FunctionEstimatorInterval::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"FunctionEstimatorInterval::splitToShape()");
	}
	bool success = false;
	IntervalMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"FunctionEstimatorInterval::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"FunctionEstimatorInterval::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);
		
		
		/* ALSO NEED to set the interval ranges using fobj */
		IntervalEstimator estimator(fobj);
		rootPaving->acceptSPValueVisitor(estimator);

		if (!success) {
			
			handleSplitToShapeError(temp);
	   }
	   
	}
	catch (std::invalid_argument const& ia) {
		cerr << ia.what() << endl;
		handleSplitToShapeError(temp);
		success = false;
	}
	catch (std::logic_error const& le) {
		cerr << le.what() << endl;
		handleSplitToShapeError(temp);
		success = false;
	}
	return success;
	// any other exceptions are unhandled
}




// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec FunctionEstimatorInterval::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}

cxsc::real FunctionEstimatorInterval::getTotalAreaOfIntervalBand() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorInterval::getTotalAreaOfIntervalBand()");
	}
	return getSubPaving()->getTotalLeafIntervalAreaOnBox();
}

cxsc::real FunctionEstimatorInterval::getMaximumAreaOfIntervalBand() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorInterval::getMaximumAreaOfIntervalBand()");
	}
	
	IntervalMappedSPnode::ConstPtrs leaves;
	getSubPaving()->getConstLeaves(leaves);
	
	cxsc::real maxArea(0.0);
	
	for (IntervalMappedSPnode::ConstPtrsItr it = leaves.begin();
			it < leaves.end();
			++it) {
				
		cxsc::real thisArea = (*it)->getIntervalAreaOnBox();
		if (thisArea > maxArea) maxArea = thisArea;
				
	}
	return maxArea;
	
}

cxsc::real FunctionEstimatorInterval::getMaximumIntervalDiameter() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorInterval::getMaximumIntervalDiameter()");
	}
	
	IntervalMappedSPnode::ConstPtrs leaves;
	getSubPaving()->getConstLeaves(leaves);
	
	cxsc::real maxDiam(0.0);
	
	for (IntervalMappedSPnode::ConstPtrsItr it = leaves.begin();
			it < leaves.end();
			++it) {
				
		interval thisRange = (*it)->getRange();
		
		cxsc::real thisDiam = cxsc::diam(thisRange);
		if (thisDiam > maxDiam) maxDiam = thisDiam;
				
	}
	return maxDiam;
	
}


cxsc::real FunctionEstimatorInterval::getMaximumIntervalTolerance() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorInterval::getMaximumIntervalTolerance()");
	}
	
	IntervalMappedSPnode::ConstPtrs leaves;
	getSubPaving()->getConstLeaves(leaves);
	
	cxsc::real maxTol(0.0);
	
	for (IntervalMappedSPnode::ConstPtrsItr it = leaves.begin();
			it < leaves.end();
			++it) {
				
		interval thisRange = (*it)->getRange();
		real thisMidImage = fobj.imageMid((*it)->getBox());
		
		cxsc::real thisTol 
				= cxsc::max(cxsc::abs(Sup(thisRange) - thisMidImage),
					cxsc::abs(thisMidImage - Inf(thisRange)));
		if (thisTol > maxTol) maxTol = thisTol;
				
	}
	return maxTol;
	
}

// Method to output the subpaving, leaves only
std::ostream & FunctionEstimatorInterval::outputToStreamTabs(std::ostream & os,
								int prec) const
{
	if (hasSubPaving()) {
		
		// have to use cxsc io manipulators
		os << cxsc::SaveOpt;
		os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

		getSubPaving()->leavesOutputTabs(os); // the output
		
		os << cxsc::RestoreOpt;
	}
	
    return os;
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void FunctionEstimatorInterval::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void FunctionEstimatorInterval::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the FunctionEstimatorInterval object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabs(os, prec); // the output
			
		}
		if (confirm)
			std::cout << "The output of the FunctionEstimatorInterval "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



void FunctionEstimatorInterval::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void FunctionEstimatorInterval::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the FunctionEstimatorInterval
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
			
		}
		if (confirm)
			std::cout << "Details of the root paving of the FunctionEstimatorInterval "
				<< "has been written to " << s << std::endl << std::endl;
				
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



// Method to output the subpaving
// all nodes, with indentation of children
std::ostream & FunctionEstimatorInterval::outputRootToStreamTabs(
													std::ostream & os,
													int prec) const
{
	if (hasSubPaving()) {
		
		getSubPaving()->nodesAllOutput(os, 1, prec); // the output
		
	}
	
    return os;
}

// Method to add current state of the histogram during splitting to a log file
// Output goes to file named according to argument s
// Output is textToTabs
void FunctionEstimatorInterval::outputLog(const std::string& s, 
											int i, int prec) const
{
    // To add output of the FunctionEstimatorInterval object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
		
		os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        getSubPaving()->leavesOutputTabs(os, prec); // the output
		
		os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}



std::string FunctionEstimatorInterval::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "Reference to function object is  = " << (&fobj) << endl;
	
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}





// --------------------------- private ---------------------------------------

IntervalMappedSPnode* FunctionEstimatorInterval::getSubPaving() const
{return rootPaving;}


// Method to put opening line into a log file
void FunctionEstimatorInterval::outputLogStart(const std::string& s) const
{
    // Make a string with filename and timestamp to start log file
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    ofstream os(s.c_str());         // replace data
    if (os.is_open()) {
        os << "File " << s << " created " <<  asctime (timeinfo) << std::endl;
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}


FunctionEstimatorInterval::PriorityQueueT& 
	FunctionEstimatorInterval::_setupPrioritySplitQueue(
						FunctionEstimatorInterval::PriorityQueueT& pq,
						const IntervalMappedSPnode::Measurer& measure)
{
	// check current state is legal
	if (!(getSubPaving()->checkTreeStateLegal()) ) {
		throw std::logic_error(
			"_setupPrioritySplitQueue(...) : Illegal state");
		
	}

	#ifdef MYDEBUGPQ
		if (getSubPaving()->isLeaf()) {
		
			cout << "root is a leaf" << endl;
		}
		else
			cout << "root has " << getRootLeaves() << " leaves" << endl;
			
	#endif
	
	IntervalMappedSPnode::Ptrs leaves;
	getSubPaving()->getLeaves(leaves);
	
	#ifdef MYDEBUGPQ
		cout << leaves.size() << " leaves" << endl;
			
	#endif
	
	IntervalMappedSPnode::PtrsItr it;
	for (it = leaves.begin(); it < leaves.end(); ++it) {
		// check to insert each of the leaves into the set
		if ( (*it)->isSplittableNode() ) {
			
			real m = measure( (*it) );
			pq.insert( NodePtrMeasurePair((*it), m) );
		}
	}
	#ifdef MYDEBUGPQ
		cout << "pq now has size " << pq.size() << endl;
		
	#endif

	return pq;
}



#if(0) // new equal largest method
bool FunctionEstimatorInterval::_prioritySplitQueueLoop(
					const IntervalMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					const SPValueVisitor<cxsc::interval>& estimator,
					gsl_rng * rgsl)
{
	bool canContinue = true;
	
	if (!pq.empty()) {
		
		NodePtrMeasurePair largest = *(pq.rbegin ()); // the last largest in the set
			
		#ifdef MYDEBUGPQ
			cout << "queue is:" << endl;
			for (PriorityQueueItrT pit = pq.begin(); pit != pq.end(); ++pit) {
				cout << "\t" << pit->toString() << endl;
			}
			cout << "first largest is " << largest.toString() << endl;
		#endif
		
		// find if there are any more equal to largest around
		PriorityQueueItrT equalLargestFirst = pq.end();
		--equalLargestFirst;
		assert((*equalLargestFirst).measure == largest.measure);
		
		size_t numberLargest = 1; // number of =largest
		while (equalLargestFirst != pq.begin()) {
			--equalLargestFirst;
			if ((*equalLargestFirst).measure < largest.measure) {
				++equalLargestFirst;
				break;
			}
			else numberLargest++;
		}
		
		#ifdef MYDEBUGPQ
			cout << "number of equal largest is " << numberLargest << endl;
		#endif

		PriorityQueueItrT mit = equalLargestFirst;
			
		if (numberLargest > 1) {

			// draw a random number in [0,1)
			double rand = gsl_rng_uniform(rgsl);
			
			#ifdef MYDEBUGPQ
				cout << "rand is " << rand << endl;
			#endif
		
			real sum = 0.0;
			bool found = false;
			// random selection of the =largest node to chose
			for (mit=equalLargestFirst; mit!=pq.end(); ++mit) {
				
				sum += 1.0/(1.0*numberLargest);
				#ifdef MYDEBUGPQ
					cout << "considering " << mit->toString();
					cout << ":\t sum = " << sum << endl;
				#endif
				if (rand < sum) {
					#ifdef MYDEBUGPQ
						cout << "choosing this one: " << mit->toString() << endl;
					#endif
					found = true;
					break;
				}
			}
			assert(found);
			
			largest = *(mit); // the chosen largest in the set
			#ifdef MYDEBUGPQ
				cout << "chosen one is " << largest.toString() << endl;;
			#endif
		}
		pq.erase(mit);// take the iterator to chosen largest out of the set

		// split the biggest one and divvie up its data

		try {
			largest.nodePtr->nodeExpand();
			#ifdef MYDEBUGPQ
				cout << "should have now split it" << endl;
			#endif
		}
		catch (UnfulfillableRequest_Error& ure) {
			std::cerr << "Something has gone wrong in the splitting.  Error reported is:"
				<< std::endl;
			std::cerr << ure.what()	<< std::endl;
			std::cerr << "Process aborted with " << getRootLeaves() << " leaves" << std::endl;
			canContinue = false;
		}
		
		#ifdef MYDEBUGPQ
			if (largest.nodePtr->isLeaf()) {
				cout << "*** note that largest is still a leaf" << endl;
			}
		#endif
		assert(!largest.nodePtr->isLeaf());
		
		if (canContinue && !largest.nodePtr->isLeaf()) {

			{
				IntervalMappedSPnode* child 
							= largest.nodePtr->getLeftChild();
				child->acceptSPValueVisitor(estimator); // put range value on child
				
				if (child->isSplittableNode()) {
					
					real m = measure( child );
					#ifdef MYDEBUGPQ
						cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
					#endif
					pq.insert(NodePtrMeasurePair(child, m));
				}
			}
			{
				IntervalMappedSPnode* child 
							= largest.nodePtr->getRightChild();
				child->acceptSPValueVisitor(estimator); // put range value on child
				
				if (child->isSplittableNode()) {
					
					real m = measure( child );
					#ifdef MYDEBUGPQ
						cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
					#endif
					pq.insert(NodePtrMeasurePair(child, m));
				}
			}
		} // end canContinue && largest not a leaf
		
	}
	else { // pq empty
		std::cerr << "Terminated splitting: no splittable nodes left"
			<< std::endl;
		
		canContinue = false;
	}
	return canContinue;
}// loop

#endif


bool FunctionEstimatorInterval::_prioritySplitQueueLoop(
					const IntervalMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					const SPValueVisitor<cxsc::interval>& estimator,
					gsl_rng * rgsl)
{
	bool canContinue = true;
	
	if (!pq.empty()) {
		

		NodePtrMeasurePair largest = *(pq.rbegin ()); // the last largest in the set
		
		#ifdef MYDEBUGPQ
			cout << "queue is:" << endl;
			for (PriorityQueueItrT pit = pq.begin(); pit != pq.end(); ++pit) {
				cout << "\t" << pit->toString() << endl;
			}
			cout << "first largest is " << largest.toString() << endl;
		#endif
		
		// find if there are any more equal to largest around
		pair< PriorityQueueItrT, PriorityQueueItrT > equalLargest;

		equalLargest = pq.equal_range(largest); // everything that = largest
		size_t numberLargest = pq.count(largest); // number of =largest
		assert(numberLargest > 0);

		#ifdef MYDEBUGPQ
			cout << "number of equal largest is " << numberLargest << endl;
			
		#endif

		PriorityQueueItrT mit = equalLargest.first;
			
		if (numberLargest > 1) {

			// draw a random number in [0,1)
			double rand = gsl_rng_uniform(rgsl);
			
			#ifdef MYDEBUGPQ
				cout << "rand is " << rand << endl;
			#endif
		
			real sum = 0.0;
			bool found = false;
			// random selection of the =largest node to chose
			for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
				
				sum += 1.0/(1.0*numberLargest);
				#ifdef MYDEBUGPQ
					cout << "considering " << mit->toString();
					cout << ":\t sum = " << sum << endl;
				#endif
				if (rand < sum) {
					#ifdef MYDEBUGPQ
						cout << "choosing this one: " << mit->toString() << endl;
					#endif
					found = true;
					break;
				}
			}
			assert(found);
			
			largest = *(mit); // the chosen largest in the set
			#ifdef MYDEBUGPQ
				cout << "chosen one is " << largest.toString() << endl;;
			#endif
		}
		pq.erase(mit);// take the iterator to chosen largest out of the set

		// split the biggest one and divvie up its data

		try {
			largest.nodePtr->nodeExpand();
			#ifdef MYDEBUGPQ
				cout << "should have now split it" << endl;
			#endif
		}
		catch (UnfulfillableRequest_Error& ure) {
			std::cerr << "Something has gone wrong in the splitting.  Error reported is:"
				<< std::endl;
			std::cerr << ure.what()	<< std::endl;
			std::cerr << "Process aborted with " << getRootLeaves() << " leaves" << std::endl;
			canContinue = false;
		}
		
		#ifdef MYDEBUGPQ
			if (largest.nodePtr->isLeaf()) {
				cout << "*** note that largest is still a leaf" << endl;
			}
		#endif
		assert(!largest.nodePtr->isLeaf());
		
		if (canContinue && !largest.nodePtr->isLeaf()) {

			{
				IntervalMappedSPnode* child 
							= largest.nodePtr->getLeftChild();
				child->acceptSPValueVisitor(estimator); // put range value on child
				
				if (child->isSplittableNode()) {
					
					real m = measure( child );
					#ifdef MYDEBUGPQ
						cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
					#endif
					pq.insert(NodePtrMeasurePair(child, m));
				}
			}
			{
				IntervalMappedSPnode* child 
							= largest.nodePtr->getRightChild();
				child->acceptSPValueVisitor(estimator); // put range value on child
				
				if (child->isSplittableNode()) {
					
					real m = measure( child );
					#ifdef MYDEBUGPQ
						cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
					#endif
					pq.insert(NodePtrMeasurePair(child, m));
				}
			}
		} // end can continue and largest not a leaf
		
	}
	else { // pq empty
		std::cerr << "Terminated splitting: no splittable nodes left"
			<< std::endl;
		
		canContinue = false;
	}
	return canContinue;
}// loop

/* I experimented with altering the priority merge so that it erases nodes
 * from the back not the start but that was slower.  pq is a set not a 
 * vector of course so presumably the set is not inefficient for erasing
 * like this and something else is making it in fact slower. Makes 
 * me wonder if I should reverse the split pq!) JAH June 2012*/


FunctionEstimatorInterval::PriorityQueueT& 
				FunctionEstimatorInterval::_setupPriorityMergeQueue(
				FunctionEstimatorInterval::PriorityQueueT& pq,
				const IntervalMappedSPnode::Measurer& measure)
{
	IntervalMappedSPnode::Ptrs cherries;
	getSubPaving()->getSubLeaves(cherries);
	#ifdef MYDEBUGPM
		cout << cherries.size() << " cherries" << endl;
			
	#endif
	IntervalMappedSPnode::PtrsItr it;
	for (it = cherries.begin(); it < cherries.end(); it++) {
		
		if ( (*it)->isSplittableNode() ) {
			
			real m = measure( (*it) );
			pq.insert( NodePtrMeasurePair((*it), m) );
		}
		else {
			throw std::logic_error(
				"_setupPriorityMergeQueue(...) : illegal state");
		}
	}

	#ifdef MYDEBUGPM
		cout << "pq now has size " << pq.size() << endl;
		
	#endif
	
	return pq;
}




bool FunctionEstimatorInterval::_priorityMergeQueueLoop(
					const IntervalMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					gsl_rng * rgsl)
{
	bool canContinue = true;
	
	if (!pq.empty()) {
	
		NodePtrMeasurePair smallest = *(pq.begin ()); // the first smallest in the set
		
		#ifdef MYDEBUGPM
			cout << "queue is:" << endl;
			for (PriorityQueueItrT pit = pq.begin(); pit != pq.end(); ++pit) {
				cout << "\t" << pit->toString() << endl;
			}
			cout << "first smallest is " << smallest.toString() << endl;
		#endif
		
		// find if there are any more equal to smallest around
		pair< PriorityQueueItrT, PriorityQueueItrT > equalSmallest;

		equalSmallest = pq.equal_range(smallest); // everything that = largest
		size_t numberSmallest = pq.count(smallest); // number of =largest
		assert(numberSmallest > 0);
		
		#ifdef MYDEBUGPM
			cout << "number of equal smallest is " << numberSmallest << endl;
			
		#endif

		PriorityQueueItrT mit = equalSmallest.first;
			
		if (numberSmallest > 1) {

			// draw a random number in [0,1)
			double rand = gsl_rng_uniform(rgsl);
			
			#ifdef MYDEBUGPM
				cout << "rand is " << rand << endl;
			#endif
		
			real sum = 0.0;
			bool found = false;
			// random selection of the =largest node to chose
			for (mit=equalSmallest.first; mit!=equalSmallest.second; ++mit) {
				
				sum += 1.0/(1.0*numberSmallest);
				#ifdef MYDEBUGPM
					cout << "considering " << mit->toString();
					cout << ":\t sum = " << sum << endl;
				#endif
				if (rand < sum) {
					#ifdef MYDEBUGPM
						cout << "choosing this one: " << mit->toString() << endl;
					#endif
					found = true;
					break;
				}
			}
			assert(found);
			
			smallest = *(mit); // the chosen smallest in the set
			#ifdef MYDEBUGPM
				cout << "chosen one is " << smallest.toString() << endl;;
			#endif
		}
		pq.erase(mit);// take the iterator to chosen smallest out of the set

		// merge 
		smallest.nodePtr->nodeReabsorbChildren();
		
		#ifdef MYDEBUGPM
			if (!smallest.nodePtr->isLeaf()) {
				cout << "*** note that smallest is not a leaf" << endl;
			}
		#endif
		assert(smallest.nodePtr->isLeaf());
		
		if (smallest.nodePtr->isLeaf()) {
			// if smallest has a leaf sibling, smallest's parent is now a cherry
			// and should be inserted into the multiset
			if (smallest.nodePtr->hasLeafSibling()) {

				IntervalMappedSPnode* parent 
							= smallest.nodePtr->getParent();
				real m = measure( parent );
				#ifdef MYDEBUGPM
					cout << "inserting " << parent->getNodeName() << " with measure " << m << " into the set" << endl;
				#endif
				pq.insert(NodePtrMeasurePair(parent, m));
		   }
		}
		
	} 
	else { // pq empty
			std::cerr << "No subleaves left to merge" << std::endl;
			
			canContinue = false;
	}
	
	return canContinue;  
}



FunctionEstimatorInterval::NodeQueueT& FunctionEstimatorInterval::_setupBreadthFirstQueue(
					NodeQueueT& nq,
					const SPCheckVisitor& nodeChecker)
{
	// check current state is legal
	if (!(getSubPaving()->checkTreeStateLegal()) ) {
		throw std::logic_error(
			"_setupBreadthFirstQueue(...) : Illegal state");
	}
	
	#ifdef MYDEBUGBF
		if (getSubPaving()->isLeaf()) {
		
			cout << "root is a leaf" << endl;
		}
		else
			cout << "root has " << getRootLeaves() << " leaves" << endl;
			
	#endif
	
	IntervalMappedSPnode::Ptrs leaves;
	getSubPaving()->getLeaves(leaves);
	
	#ifdef MYDEBUGBF
		cout << leaves.size() << " leaves" << endl;
			
	#endif

	IntervalMappedSPnode::PtrsItr it;
	for (it = leaves.begin(); it < leaves.end(); ++it) {
		// check to insert each of the leaves into the set
		(*it)->acceptSPCheckVisitor(nodeChecker);
		if ( !nodeChecker.getResult() ) { // does not yet meet requirement
			nq.push((*it));
		}
	}
	#ifdef MYDEBUGBF
		cout << "nq now has size " << nq.size() << endl;
		
	#endif
	
	return nq;
}

void FunctionEstimatorInterval::_breadthFirstQueueLoop(
					const SPCheckVisitor& nodeChecker,
					NodeQueueT& nq,
					const SPValueVisitor<cxsc::interval>& estimator,
					gsl_rng * r)
{
	#ifdef MYDEBUGBF
		cout << "queue size is " << nq.size() << endl;
		cout << "front of queue is " << nq.front()->getNodeName() << endl;
		cout << "back of queue is " << nq.back()->getNodeName() << endl;
		
	#endif
		
			
	// expand the one at the front
	nq.front()->nodeExpand();
	
	double rand = gsl_rng_uniform(r);
	//shove its children in at the back, random order
	IntervalMappedSPnode* child1 = NULL;
	IntervalMappedSPnode* child2 = NULL;
	if( rand < 0.5) {
		child1 = nq.front()->getLeftChild();
		child2 = nq.front()->getRightChild();
	}
	else {
		child1 = nq.front()->getRightChild();
		child2 = nq.front()->getLeftChild();
	}
	child1->acceptSPValueVisitor(estimator); // put range value on child1
	child2->acceptSPValueVisitor(estimator); // put range value on child2

	//check to see if we want to put children into the queue
	nodeChecker.visit(child1);
	if (!(nodeChecker.getResult())) {
		nq.push(child1);
		#ifdef MYDEBUGBF
			cout << "added " << child1->getNodeName() << " to the back" << endl;
				
		#endif
	}
	nodeChecker.visit(child2);
	if (!(nodeChecker.getResult())) {
		nq.push(child2);
		#ifdef MYDEBUGBF
			cout << "added " << child2->getNodeName() << " to the back" << endl;
				
		#endif
	}

	//pop the front of the queue
	nq.pop();
	
}
		

//check that the box is okay for the histogram
bool FunctionEstimatorInterval::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void FunctionEstimatorInterval::handleSplitToShapeError(
											IntervalMappedSPnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
			std::cerr << "Your instruction does not describe a proper tree.";
			std::cerr << "  Please check your instruction and try again."
			<< std::endl;
}

// ensure rootPaving is deleted if constructed in failed constructor
void FunctionEstimatorInterval::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}

// implementations for private inner class for nodes and measurements pairs

FunctionEstimatorInterval::NodePtrMeasurePair::NodePtrMeasurePair(
	IntervalMappedSPnode *p, real m) : nodePtr(p), measure(m) {}

bool FunctionEstimatorInterval::NodePtrMeasurePair::operator<(
		const NodePtrMeasurePair& rhs) const
{
	return measure < rhs.measure;
}


std::string 
	FunctionEstimatorInterval::NodePtrMeasurePair::toString() const
{
	ostringstream oss;
	oss << (nodePtr->getNodeName()) << ":\t" << measure;
	return oss.str();
}

// ---------------- implmentation for private inner class for measuring nodes
FunctionEstimatorInterval::
	IntervalBandAreaGainOnSplitMeasurer::
	IntervalBandAreaGainOnSplitMeasurer(const MappedFobj& mf) 
		: measurefobj(mf) {}
			
cxsc::real FunctionEstimatorInterval::IntervalBandAreaGainOnSplitMeasurer::
			operator()(const IntervalMappedSPnode * const imspn) const
{
	/* use spn rather than box because we cannot be sure that 
	 * the interval on spn has been set using measurefobj*/
	cxsc::ivector box = imspn->getBox();
	cxsc::ivector lcBox;
	cxsc::ivector rcBox;
	
	int c = MaxDiamComp(box);
	Lower (box, lcBox, c);
	Upper (box, rcBox, c);
	
	#ifdef MYDEBUGPQGAINCALCS
		cout << "box is " << box << endl;
		cout << "lcBox is " << lcBox << endl;
		cout << "rcBox is " << rcBox << endl;
		cout << "imspn->getRange() is " 
			<< imspn->getRange() << " ";
		cout << "imspn->getRangeDiameter() is " 
			<< imspn->getRangeDiameter() << endl;
		cout << "realVolume(box) is " << (realVolume(box)) << endl;
		cout << "0.5*realVolume(box) is " << (0.5*realVolume(box)) << endl;
		cout << "imspn->getIntervalAreaOnBox() is " 
			<< imspn->getIntervalAreaOnBox() << endl;
		cout << "measurefobj(lcBox) is " << (measurefobj(lcBox)) << endl;
		cout << "measurefobj(rcBox) is " << (measurefobj(rcBox)) << endl;
		cout << "cxsc::diam( measurefobj(lcBox)) is " << (cxsc::diam( measurefobj(lcBox))) << "\t";
		cout << "realVolume(box) * 0.5 * ( cxsc::diam( measurefobj(lcBox) ) is " 
			<< (realVolume(box) * 0.5 * ( cxsc::diam( measurefobj(lcBox) ))) << endl;
		cout << "cxsc::diam( measurefobj(rcBox)) is " << (cxsc::diam( measurefobj(rcBox))) << "\t";
		cout << "realVolume(box) * 0.5 * ( cxsc::diam( measurefobj(rcBox) ) is " 
			<< (realVolume(box) * 0.5 * ( cxsc::diam( measurefobj(rcBox) ))) << endl;
	
	#endif
	
	return ( imspn->getIntervalAreaOnBox() 
		- realVolume(box) * 0.5 * ( cxsc::diam( measurefobj(lcBox) ) + 
									cxsc::diam( measurefobj(rcBox) ) ) );

}

/* ---------------- implementation for private inner class for increase in 
		interval band area if a split node were to be merged.*/
FunctionEstimatorInterval::
	IntervalBandAreaLossOnMergeMeasurer::
	IntervalBandAreaLossOnMergeMeasurer() {}
			
cxsc::real FunctionEstimatorInterval::IntervalBandAreaLossOnMergeMeasurer::
			operator()(const IntervalMappedSPnode * const imspn) const
{
	return imspn->getIntervalAreaDiffToChildren();
}

// ----------------------------- non member functions

//Output all boxes in FunctionEstimatorInterval adh
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::FunctionEstimatorInterval& fei)
{
    fei.outputRootToStreamTabs(os);
    return os;
}










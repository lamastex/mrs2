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
\brief FunctionEstimatorKDE definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "functionestimator_kde.hpp"
#include "real_kde_mid_estimate.hpp"
//#include "realexpander_estimate.hpp"
#include "toolz.hpp"
#include "sptools.hpp"


#include "subpaving_exception.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <cassert>


//#define MYDEBUGPQ
//#define MYDEBUGPM
//#define MYDEBUGPQMIN
//#define MYDEBUGPQGAINCALCS
//#define GET_MEASURE


#ifdef NDEBUG
	#undef MYDEBUGPQ
	#undef MYDEBUGPM
	#undef MYDEBUGPQMIN
	#undef MYDEBUGPQGAINCALCS
	#undef GET_MEASURE

#endif

#if defined (MYDEBUGPQMIN)
	#include <ctime>
#else
	#define SHOW_TIMEDOTS
#endif



using namespace subpavings;
using namespace std;



// -------------------implementation of FunctionEstimatorKDE class --------------


// ----------- public methods


// initialised constructor 
FunctionEstimatorKDE::FunctionEstimatorKDE(
										const ivector& v,
										const RealPointwiseFunctionEstimator& f,
										int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"FunctionEstimatorKDE::FunctionEstimatorKDE(const ivector&, const TypeKDE&, int lab)");
		}
        rootPaving = new RealMappedSPnode(v);
		
		rootPaving->acceptSPValueVisitor(fobj);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

FunctionEstimatorKDE::FunctionEstimatorKDE(const SPnode& spn, 
											const RealPointwiseFunctionEstimator& f, 
											int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check spn has box
		if (spn.isEmpty()) {
			throw subpavings::NoBox_Error(
			"FunctionEstimatorKDE::FunctionEstimatorKDE(const SPnode&, TypeKDE&, int lab");
		}
        rootPaving = new RealMappedSPnode(spn);
		
		rootPaving->acceptSPValueVisitor(fobj);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// copy constructor*/
FunctionEstimatorKDE::FunctionEstimatorKDE(
								const FunctionEstimatorKDE& other)
        : rootPaving(NULL), fobj(other.fobj), label(other.label)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new RealMappedSPnode(*(other.getSubPaving()));
			
		} // else subpaving is NULL which should be impossible
		else {
			throw NullSubpavingPointer_Error(
				"FunctionEstimatorKDE::FunctionEstimatorKDE(const FunctionEstimatorKDE&)");
		}
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}



//Destructor
FunctionEstimatorKDE::~FunctionEstimatorKDE()
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
			std::cerr << "Error in FunctionEstimatorKDE destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

/*copy assignment operator private and not implemented:
Cannot do assignment because fobj is a reference.*/


const RealPointwiseFunctionEstimator& FunctionEstimatorKDE::getFobjReference() const
{return fobj;}

int FunctionEstimatorKDE::getLabel() const
{return label;}

void FunctionEstimatorKDE::setLabel(int lab)
{label = lab;}

// get whether this has a subpaving.
bool FunctionEstimatorKDE::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

cxsc::ivector FunctionEstimatorKDE::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionEstimatorKDE::getRootBox()");
	}
	return getSubPaving()->getBox();
}



int FunctionEstimatorKDE::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

cxsc::real FunctionEstimatorKDE::getDomainVolume() const
{
	real retValue(0.0);
	if (hasSubPaving()) {
		retValue = getSubPaving()->nodeRealVolume();
	}
	return retValue;
}
	
	
// Gets number of leaf nodes in the root paving.
size_t FunctionEstimatorKDE::getRootLeaves() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("FunctionEstimatorKDE::getRootLeaves()");
		
	}
	return getSubPaving()->getNumberLeaves();
}

// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec FunctionEstimatorKDE::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// Get a string of the leaf node levels.
std::string FunctionEstimatorKDE::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}

subpavings::PiecewiseConstantFunction 
		FunctionEstimatorKDE::makePiecewiseConstantFunction() const
{
	return PiecewiseConstantFunction(*getSubPaving(), getLabel());
}



bool FunctionEstimatorKDE::prioritySplit(
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	FunctionEstimatorKDE::
	TotalVariationSplitMeasurer measure(fobj);
	real maxMeasure = cxsc::MaxReal;
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, seed);
}

bool FunctionEstimatorKDE::prioritySplit(
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	FunctionEstimatorKDE::
	TotalVariationSplitMeasurer measure(fobj);
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, seed);
}

bool FunctionEstimatorKDE::prioritySplit(
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	FunctionEstimatorKDE::
	TotalVariationSplitMeasurer measure(fobj);
	real maxMeasure = cxsc::MaxReal;
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, rgsl);
}

bool FunctionEstimatorKDE::prioritySplit(
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	FunctionEstimatorKDE::
	TotalVariationSplitMeasurer measure(fobj);
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, rgsl);
}


bool FunctionEstimatorKDE::prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	real maxMeasure = cxsc::MaxReal;
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, seed);
}

 //  using maxLeaves
bool FunctionEstimatorKDE::prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	
	gsl_rng * rgsl = NULL;

    try {
        string errorMsg("FunctionEstimatorKDE::prioritySplit(...)");
		
		if (!hasSubPaving()) {
			throw NullSubpavingPointer_Error(errorMsg);
		}

		// set up a random number generator for uniform rvs
        rgsl = gsl_rng_alloc (gsl_rng_mt19937);
					
		gsl_rng_set(rgsl, seed);

		bool retValue = prioritySplit(measure,
						maxMeasure,
						maxLeaves,
						logging,
						rgsl);
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		
		return retValue;
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

//  using maxLeaves
bool FunctionEstimatorKDE::prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	real maxMeasure = cxsc::MaxReal;
	return prioritySplit(measure, maxMeasure, maxLeaves, logging, rgsl);
}

//  using maxLeaves
bool FunctionEstimatorKDE::prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	string errorMsg("FunctionEstimatorKDE::prioritySplit(...)");
	
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

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
	
	bool noCheckMeasure = !(maxMeasure < cxsc::MaxReal);
	
	#ifdef MYDEBUGPQMIN
		clock_t start = clock();
	#endif


	// split until we have enough leaves
	/* check there are nodes in queue at the top of loop rather
	 * than the bottom to ensure that we will stop if enough leaves 
	 * before we check queue has nodes in,
	 * ie only return false if max leaves not achieved end of operation*/ 
	while (canContinue && numLeaves < maxLeaves &&
		(noCheckMeasure || pq.empty() || (pq.rbegin()->measure > maxMeasure) )) {
			// note we want to enter the loop even if pq is empty, to get back canContinue = false.
		
		canContinue = _prioritySplitQueueLoop(
				measure, pq, rgsl);
		
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
		if (numLeaves%10000 == 0) {
				cout << "." << flush;
		}
		#endif
		#ifdef MYDEBUGPQ
			cout << "\nat end of loop, number of leaf nodes in subpaving is " << (getRootLeaves()) << endl;
				
		#endif
		#ifdef MYDEBUGPQMIN
			if (numLeaves%10000 == 0) {
				double time = static_cast<double>(clock()-start)/CLOCKS_PER_SEC;
				cout << numLeaves << "\t" << time << endl;
				start = clock();
			}
		#endif
		#ifdef MYDEBUGPQ
			cout << "current estimator is " << endl;
			outputToStreamTabs(cout, 5);
			
		#endif
	}// loop
	#ifdef SHOW_TIMEDOTS
		cout << endl;
	#endif
	
	#ifdef GET_MEASURE
		if (numLeaves == maxLeaves) {
			NodePtrMeasurePair largest = *(pq.rbegin ()); // the last largest in the set
			cout << "With maxLeaves, the largest measure is " << largest.measure << endl;
		}
	#endif

	return canContinue;
}



//splits this according to string instruction
//returns true if some splitting was achieved
bool FunctionEstimatorKDE::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"FunctionEstimatorKDE::splitToShape()");
	}
	bool success = false;
	RealMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"FunctionEstimatorKDE::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"FunctionEstimatorKDE::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);
		
		
		rootPaving->acceptSPValueVisitor(fobj);

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



cxsc::real FunctionEstimatorKDE::getTotalIntegralOfRealEstimate() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorKDE::getTotalIntegralOfRealEstimate)");
	}
	return getSubPaving()->getTotalAbsLeafAreaRangeWithBox();
}


// Method to output the subpaving, leaves only
std::ostream & FunctionEstimatorKDE::outputToStreamTabs(std::ostream & os,
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
void FunctionEstimatorKDE::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void FunctionEstimatorKDE::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the FunctionEstimatorKDE object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabs(os, prec); // the output
			
		}
		if (confirm)
			std::cout << "The output of the FunctionEstimatorKDE "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



void FunctionEstimatorKDE::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void FunctionEstimatorKDE::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the FunctionEstimatorKDE
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
			
		}
		if (confirm)
			std::cout << "Details of the root paving of the FunctionEstimatorKDE "
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
std::ostream & FunctionEstimatorKDE::outputRootToStreamTabs(
													std::ostream & os,
													int prec) const
{
	if (hasSubPaving()) {
		
		getSubPaving()->nodesAllOutput(os, 1, prec); // the output
		
	}
	
    return os;
}

// Method to add current state of the estimator during splitting to a log file
// Output goes to file named according to argument s
// Output is textToTabs
void FunctionEstimatorKDE::outputLog(const std::string& s, 
											int i, int prec) const
{
    // To add output of the FunctionEstimatorKDE object to file
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



std::string FunctionEstimatorKDE::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "Reference to kde function object is  = " << (&fobj) << endl;
	
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}





// --------------------------- private ---------------------------------------


RealMappedSPnode* FunctionEstimatorKDE::getSubPaving() const
{return rootPaving;}



FunctionEstimatorKDE::PriorityQueueT& 
	FunctionEstimatorKDE::_setupPrioritySplitQueue(
						FunctionEstimatorKDE::PriorityQueueT& pq,
						const RealMappedSPnode::Measurer& measure)
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
	
	RealMappedSPnode::Ptrs leaves;
	getSubPaving()->getLeaves(leaves);
	
	#ifdef MYDEBUGPQ
		cout << leaves.size() << " leaves" << endl;
			
	#endif
	
	RealMappedSPnode::PtrsItr it;
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


bool FunctionEstimatorKDE::_prioritySplitQueueLoop(
					const RealMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
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
				RealMappedSPnode* child 
							= largest.nodePtr->getLeftChild();
				child->acceptSPValueVisitor(fobj); // put range value on child
				
				if (child->isSplittableNode()) {
					
					real m = measure( child );
					#ifdef MYDEBUGPQ
						cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
					#endif
					pq.insert(NodePtrMeasurePair(child, m));
				}
			}
			{
				RealMappedSPnode* child 
							= largest.nodePtr->getRightChild();
				child->acceptSPValueVisitor(fobj); // put range value on child
				
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



// Method to put opening line into a log file
void FunctionEstimatorKDE::outputLogStart(const std::string& s) const
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


//check that the box is okay for the estimator
bool FunctionEstimatorKDE::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void FunctionEstimatorKDE::handleSplitToShapeError(
											RealMappedSPnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
			std::cerr << "Your instruction does not describe a proper tree.";
			std::cerr << "  Please check your instruction and try again."
			<< std::endl;
}

// ensure rootPaving is deleted if constructed in failed constructor
void FunctionEstimatorKDE::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}





// implementations for private inner class for nodes and measurements pairs

FunctionEstimatorKDE::NodePtrMeasurePair::NodePtrMeasurePair(
	RealMappedSPnode *p, real m) : nodePtr(p), measure(m) {}

bool FunctionEstimatorKDE::NodePtrMeasurePair::operator<(
		const NodePtrMeasurePair& rhs) const
{
	return measure < rhs.measure;
}


std::string 
	FunctionEstimatorKDE::NodePtrMeasurePair::toString() const
{
	ostringstream oss;
	oss << (nodePtr->getNodeName()) << ":\t" << measure;
	return oss.str();
}

// ---------------- implmentation for private inner class for measuring nodes
FunctionEstimatorKDE::
	TotalVariationSplitMeasurer::
	TotalVariationSplitMeasurer(const RealPointwiseFunctionEstimator& mf) 
		: measurefobj(mf) {}
			
cxsc::real FunctionEstimatorKDE::TotalVariationSplitMeasurer::
			operator()(const RealMappedSPnode * const rmspn) const
{
	
	real thisRange = rmspn->getRange();
	
	
	cxsc::ivector box = rmspn->getBox();
	
	real thisVol = realVolume(box);
	
	cxsc::ivector lcBox;
	cxsc::ivector rcBox;
	
	int c = MaxDiamComp(box);
	Lower (box, lcBox, c);
	Upper (box, rcBox, c);
	
	
	
	real lcMeasure = measurefobj(lcBox);
	real rcMeasure = measurefobj(rcBox);
	
	real totalVar = 0.5 * thisVol * (cxsc::abs(thisRange - lcMeasure) 
				+ cxsc::abs(thisRange - rcMeasure) );
	
	
	
	#ifdef MYDEBUGPQGAINCALCS
		cout << "box is " << box << endl;
		cout << "lcBox is " << lcBox << endl;
		cout << "rcBox is " << rcBox << endl;
		cout << "thisRange = rmspn->getRange() is " 
			<< rmspn->getRange() << " ";
		cout << "realVolume(box) is " << thisVol << endl;
		cout << "0.5*realVolume(box) is " << (0.5*thisVol) << endl;
		cout << "measurefobj(lcBox) is " << lcMeasure << endl;
		cout << "measurefobj(rcBox) is " << rcMeasure << endl;
		cout << "realVolume(box) * 0.5 * ( cxsc::abs(thisRange - lcMeasure) "
			<< " + cxsc::abs(thisRange - rcMeasure) ) is " << totalVar << endl;
	
	#endif
	
	return totalVar;

}

/* ---------------- implementation for private inner class for total variation
 *  if a split node were to be merged.*/
FunctionEstimatorKDE::
	TotalVariationMergeMeasurer::
	TotalVariationMergeMeasurer() {}
			
cxsc::real FunctionEstimatorKDE::TotalVariationMergeMeasurer::
			operator()(const RealMappedSPnode * const rmspn) const
{
	
	real thisRange = rmspn->getRange();
	real lcRange = rmspn->getLeftChild()->getRange();
	real rcRange = rmspn->getRightChild()->getRange();
	
	cxsc::ivector box = rmspn->getBox();
	
	real thisVol = realVolume(box);
	
	real totalVar = 0.5 * thisVol * (cxsc::abs(thisRange - lcRange) 
				+ cxsc::abs(thisRange - rcRange) );
	
	
	
	return totalVar;
}






// ----------------------------- non member functions

//Output all boxes in FunctionEstimatorKDE adh
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::FunctionEstimatorKDE& fek)
{
    fek.outputRootToStreamTabs(os);
    return os;
}










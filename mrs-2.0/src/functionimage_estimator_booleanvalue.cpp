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
\brief FunctionImageEstimatorBooleanValue definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "functionimage_estimator_booleanvalue.hpp"
#include "booleanvalue_image_expander_estimate.hpp"
#include "booleanvalue_image_estimate.hpp"
//#include "BooleanValueMappedSPnode_measurers.hpp"
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
#endif


// -------------------implementation of FunctionImageEstimatorBooleanValue class --------------


// ----------- public methods


// initialised constructor 
FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(
										const ivector& v,
										const MappedFobj& f,
										const cxsc::interval& crit, 
										int lab)
         : rootPaving(NULL), fobj(f), criterion(crit), label(lab)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(const ivector&, const MappedFobj&, int lab)");
		}
        rootPaving = new BooleanValueMappedSPnode(v);
		
		BooleanValueImageEstimator estimator(fobj, crit);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(const SPnode& spn, 
											const MappedFobj& f, 
											const cxsc::interval& crit, 
											int lab)
         : rootPaving(NULL), fobj(f), criterion(crit), label(lab)
{
    try {
        // check spn has box
		if (spn.isEmpty()) {
			throw subpavings::NoBox_Error(
			"FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(const SPnode&, MappedFobj&, int lab");
		}
        rootPaving = new BooleanValueMappedSPnode(spn);
		
		BooleanValueImageEstimator estimator(fobj, crit);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// copy constructor*/
FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(
								const FunctionImageEstimatorBooleanValue& other)
        : rootPaving(NULL), fobj(other.fobj), criterion(other.criterion),
			label(other.label)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new BooleanValueMappedSPnode(*(other.getSubPaving()));
			
		} // else subpaving is NULL which should be impossible
		else {
			throw NullSubpavingPointer_Error(
				"FunctionImageEstimatorBooleanValue::FunctionImageEstimatorBooleanValue(const FunctionImageEstimatorBooleanValue&)");
		}
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}



//Destructor
FunctionImageEstimatorBooleanValue::~FunctionImageEstimatorBooleanValue()
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
			std::cerr << "Error in FunctionImageEstimatorBooleanValue destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

/*copy assignment operator private and not implemented:
Cannot do assignment because fobj is a reference.*/


const MappedFobj& FunctionImageEstimatorBooleanValue::getFobjReference() const
{return fobj;}

cxsc::interval FunctionImageEstimatorBooleanValue::getCriterion() const
{return criterion;}

int FunctionImageEstimatorBooleanValue::getLabel() const
{return label;}

void FunctionImageEstimatorBooleanValue::setLabel(int lab)
{label = lab;}

// get whether this has a subpaving.
bool FunctionImageEstimatorBooleanValue::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

cxsc::ivector FunctionImageEstimatorBooleanValue::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionImageEstimatorBooleanValue::getRootBox()");
	}
	return getSubPaving()->getBox();
}

bool FunctionImageEstimatorBooleanValue::getRootRangeValue() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionImageEstimatorBooleanValue::getRootRangeDiameter()");
	}
	return getSubPaving()->getRange();
}

int FunctionImageEstimatorBooleanValue::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

cxsc::real FunctionImageEstimatorBooleanValue::getDomainVolume() const
{
	real retValue(0.0);
	if (hasSubPaving()) {
		retValue = getSubPaving()->nodeRealVolume();
	}
	return retValue;
}
	

// Gets number of leaf nodes in the root paving.
size_t FunctionImageEstimatorBooleanValue::getRootLeaves() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("FunctionImageEstimatorBooleanValue::getRootLeaves()");
		
	}
	return getSubPaving()->getNumberLeaves();
}


// Get a string of the leaf node levels.
std::string FunctionImageEstimatorBooleanValue::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}


subpavings::SpatialObjectRepresentationBV 
			FunctionImageEstimatorBooleanValue::makeSpatialObjectRepresentationBV() const
{
	return SpatialObjectRepresentationBV(*getSubPaving(), getLabel());
}

void FunctionImageEstimatorBooleanValue::bruteForceEstimate(cxsc::real tolerance)
{
	string errorMsg(
		"FunctionImageEstimatorBooleanValue::bruteForceEstimate(cxsc::real)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	BooleanValueImageExpanderEstimator estimator(fobj, criterion,
									tolerance);

	getSubPaving()->acceptSPExpandVisitor(estimator);

}





//splits function image according to string instruction
//returns true if some splitting was achieved
bool FunctionImageEstimatorBooleanValue::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"FunctionImageEstimatorBooleanValue::splitToShape()");
	}
	bool success = false;
	BooleanValueMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"FunctionImageEstimatorBooleanValue::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"FunctionImageEstimatorBooleanValue::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);
		
		
		/* ALSO NEED to set the interval ranges using fobj */
		BooleanValueImageEstimator estimator(fobj, criterion);
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
IntVec FunctionImageEstimatorBooleanValue::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}

cxsc::real FunctionImageEstimatorBooleanValue::getTotalTrueArea() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionImageEstimatorBooleanValue::getTotalAreaOfIntervalBand()");
	}
	return getSubPaving()->getTotalLeafTrueAreaOnBox();
}


// Method to output the subpaving, leaves only
std::ostream & FunctionImageEstimatorBooleanValue::outputToStreamTabs(std::ostream & os,
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
void FunctionImageEstimatorBooleanValue::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void FunctionImageEstimatorBooleanValue::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the FunctionImageEstimatorBooleanValue object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabs(os, prec); // the output
			
		}
		if (confirm)
			std::cout << "The output of the FunctionImageEstimatorBooleanValue "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



void FunctionImageEstimatorBooleanValue::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void FunctionImageEstimatorBooleanValue::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the FunctionImageEstimatorBooleanValue
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
			
		}
		if (confirm)
			std::cout << "Details of the root paving of the FunctionImageEstimatorBooleanValue "
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
std::ostream & FunctionImageEstimatorBooleanValue::outputRootToStreamTabs(
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
void FunctionImageEstimatorBooleanValue::outputLog(const std::string& s, 
											int i, int prec) const
{
    // To add output of the FunctionImageEstimatorBooleanValue object to file
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



std::string FunctionImageEstimatorBooleanValue::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "Reference to function object is  = " << (&fobj) << endl;
	
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}





// --------------------------- private ---------------------------------------

BooleanValueMappedSPnode* FunctionImageEstimatorBooleanValue::getSubPaving() const
{return rootPaving;}


// Method to put opening line into a log file
void FunctionImageEstimatorBooleanValue::outputLogStart(const std::string& s) const
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

#if(0)
FunctionImageEstimatorBooleanValue::PriorityQueueT& 
	FunctionImageEstimatorBooleanValue::_setupPrioritySplitQueue(
						FunctionImageEstimatorBooleanValue::PriorityQueueT& pq,
						const BooleanValueMappedSPnode::Measurer& measure)
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
	
	BooleanValueMappedSPnode::Ptrs leaves;
	getSubPaving()->getLeaves(leaves);
	
	#ifdef MYDEBUGPQ
		cout << leaves.size() << " leaves" << endl;
			
	#endif
	
	BooleanValueMappedSPnode::PtrsItr it;
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



bool FunctionImageEstimatorBooleanValue::_prioritySplitQueueLoop(
					const BooleanValueMappedSPnode::Measurer& measure,
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
				BooleanValueMappedSPnode* child 
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
				BooleanValueMappedSPnode* child 
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

#endif
		

//check that the box is okay for the histogram
bool FunctionImageEstimatorBooleanValue::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void FunctionImageEstimatorBooleanValue::handleSplitToShapeError(
											BooleanValueMappedSPnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
			std::cerr << "Your instruction does not describe a proper tree.";
			std::cerr << "  Please check your instruction and try again."
			<< std::endl;
}

// ensure rootPaving is deleted if constructed in failed constructor
void FunctionImageEstimatorBooleanValue::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}

#if(0)

// implementations for private inner class for nodes and measurements pairs

FunctionImageEstimatorBooleanValue::NodePtrMeasurePair::NodePtrMeasurePair(
	BooleanValueMappedSPnode *p, real m) : nodePtr(p), measure(m) {}

bool FunctionImageEstimatorBooleanValue::NodePtrMeasurePair::operator<(
		const NodePtrMeasurePair& rhs) const
{
	return measure < rhs.measure;
}


std::string 
	FunctionImageEstimatorBooleanValue::NodePtrMeasurePair::toString() const
{
	ostringstream oss;
	oss << (nodePtr->getNodeName()) << ":\t" << measure;
	return oss.str();
}
#endif

// ----------------------------- non member functions

//Output all boxes in FunctionImageEstimatorBooleanValue adh
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::FunctionImageEstimatorBooleanValue& fei)
{
    fei.outputRootToStreamTabs(os);
    return os;
}










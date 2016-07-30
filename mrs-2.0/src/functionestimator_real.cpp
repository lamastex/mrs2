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
\brief FunctionEstimatorReal definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "functionestimator_real.hpp"
#include "real_estimate.hpp"
#include "realexpander_estimate.hpp"
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

#ifdef NDEBUG
	#undef MYDEBUGPQ
	#undef MYDEBUGPM

#endif

using namespace subpavings;
using namespace std;



// -------------------implementation of FunctionEstimatorReal class --------------


// ----------- public methods


// initialised constructor 
FunctionEstimatorReal::FunctionEstimatorReal(
										const ivector& v,
										const MappedFobj& f,
										int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"FunctionEstimatorReal::FunctionEstimatorReal(const ivector&, const MappedFobj&, int lab)");
		}
        rootPaving = new RealMappedSPnode(v);
		
		RealEstimator estimator(fobj);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

FunctionEstimatorReal::FunctionEstimatorReal(const SPnode& spn, 
											const MappedFobj& f, 
											int lab)
         : rootPaving(NULL), fobj(f), label(lab)
{
    try {
        // check spn has box
		if (spn.isEmpty()) {
			throw subpavings::NoBox_Error(
			"FunctionEstimatorReal::FunctionEstimatorReal(const SPnode&, MappedFobj&, int lab");
		}
        rootPaving = new RealMappedSPnode(spn);
		
		RealEstimator estimator(fobj);
		rootPaving->acceptSPValueVisitor(estimator);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// copy constructor*/
FunctionEstimatorReal::FunctionEstimatorReal(
								const FunctionEstimatorReal& other)
        : rootPaving(NULL), fobj(other.fobj), label(other.label)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new RealMappedSPnode(*(other.getSubPaving()));
			
		} // else subpaving is NULL which should be impossible
		else {
			throw NullSubpavingPointer_Error(
				"FunctionEstimatorReal::FunctionEstimatorReal(const FunctionEstimatorReal&)");
		}
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}



//Destructor
FunctionEstimatorReal::~FunctionEstimatorReal()
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
			std::cerr << "Error in FunctionEstimatorReal destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

/*copy assignment operator private and not implemented:
Cannot do assignment because fobj is a reference.*/


const MappedFobj& FunctionEstimatorReal::getFobjReference() const
{return fobj;}

int FunctionEstimatorReal::getLabel() const
{return label;}

void FunctionEstimatorReal::setLabel(int lab)
{label = lab;}

// get whether this has a subpaving.
bool FunctionEstimatorReal::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

bool FunctionEstimatorReal::hasNegativeFunctionEstimates() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionEstimatorReal::hasNegativeFunctionEstimates()");
	}
	return getSubPaving()->hasNegativeRangeInTree();
}

cxsc::ivector FunctionEstimatorReal::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"FunctionEstimatorReal::getRootBox()");
	}
	return getSubPaving()->getBox();
}



int FunctionEstimatorReal::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

cxsc::real FunctionEstimatorReal::getDomainVolume() const
{
	real retValue(0.0);
	if (hasSubPaving()) {
		retValue = getSubPaving()->nodeRealVolume();
	}
	return retValue;
}
	
	
// Gets number of leaf nodes in the root paving.
size_t FunctionEstimatorReal::getRootLeaves() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("FunctionEstimatorReal::getRootLeaves()");
		
	}
	return getSubPaving()->getNumberLeaves();
}

// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec FunctionEstimatorReal::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// Get a string of the leaf node levels.
std::string FunctionEstimatorReal::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}

subpavings::PiecewiseConstantFunction 
		FunctionEstimatorReal::makePiecewiseConstantFunction() const
{
	return PiecewiseConstantFunction(*getSubPaving(), getLabel());
}

void FunctionEstimatorReal::bruteForceEstimate(cxsc::real tolerance)
{
	string errorMsg(
		"FunctionEstimatorReal::bruteForceEstimate(cxsc::real)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	RealExpanderEstimator estimator(fobj, tolerance);

	getSubPaving()->acceptSPExpandVisitor(estimator);

}



//splits this according to string instruction
//returns true if some splitting was achieved
bool FunctionEstimatorReal::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"FunctionEstimatorReal::splitToShape()");
	}
	bool success = false;
	RealMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"FunctionEstimatorReal::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"FunctionEstimatorReal::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);
		
		
		/* ALSO NEED to set the interval ranges using fobj */
		RealEstimator estimator(fobj);
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



cxsc::real FunctionEstimatorReal::getTotalIntegralOfRealEstimate() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"FunctionEstimatorReal::getTotalIntegralOfRealEstimate)");
	}
	return getSubPaving()->getTotalAbsLeafAreaRangeWithBox();
}


// Method to output the subpaving, leaves only
std::ostream & FunctionEstimatorReal::outputToStreamTabs(std::ostream & os,
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
void FunctionEstimatorReal::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void FunctionEstimatorReal::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the FunctionEstimatorReal object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabs(os, prec); // the output
			
		}
		if (confirm)
			std::cout << "The output of the FunctionEstimatorReal "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



void FunctionEstimatorReal::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void FunctionEstimatorReal::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the FunctionEstimatorReal
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
			
		}
		if (confirm)
			std::cout << "Details of the root paving of the FunctionEstimatorReal "
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
std::ostream & FunctionEstimatorReal::outputRootToStreamTabs(
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
void FunctionEstimatorReal::outputLog(const std::string& s, 
											int i, int prec) const
{
    // To add output of the FunctionEstimatorReal object to file
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



std::string FunctionEstimatorReal::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "Reference to function object is  = " << (&fobj) << endl;
	
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}





// --------------------------- private ---------------------------------------


RealMappedSPnode* FunctionEstimatorReal::getSubPaving() const
{return rootPaving;}


// Method to put opening line into a log file
void FunctionEstimatorReal::outputLogStart(const std::string& s) const
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
bool FunctionEstimatorReal::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void FunctionEstimatorReal::handleSplitToShapeError(
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
void FunctionEstimatorReal::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}


// ----------------------------- non member functions

//Output all boxes in FunctionEstimatorReal adh
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::FunctionEstimatorReal& fei)
{
    fei.outputRootToStreamTabs(os);
    return os;
}










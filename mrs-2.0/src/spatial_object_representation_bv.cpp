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
\brief SpatialObjectRepresentationBV definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "spatial_object_representation_bv.hpp"

#include "subpaving_exception.hpp"

#include <numeric>
#include <functional>
#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <cassert>

#ifdef NDEBUG
	

#endif

using namespace subpavings;
using namespace std;



// -------------------implementation of SpatialObjectRepresentationBV class --------------


// ----------- public methods

// no-argument constructor 
SpatialObjectRepresentationBV::SpatialObjectRepresentationBV()
         : rootPaving(NULL), label(0)
{}

// initialised constructor 
SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(
										const ivector& v,
										int lab)
         : rootPaving(NULL), label(lab)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(const ivector&, int lab)");
		}
        rootPaving = new BooleanValueMappedSPnode(v, false); 
		
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(
								const BooleanValueMappedSPnode& bmspn, 
								int lab)
         : rootPaving(NULL), label(lab)
{
    try {
        // check spn has box
		if (bmspn.isEmpty()) {
			throw subpavings::NoBox_Error(
			"SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(const BooleanValueMappedSPnode&, int lab)");
		}
        
		rootPaving = new BooleanValueMappedSPnode(bmspn);
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}



// copy constructor*/
SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(
								const SpatialObjectRepresentationBV& other)
        : rootPaving(NULL), label(other.label)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new BooleanValueMappedSPnode(*(other.getSubPaving()));
			
		} // else subpaving is NULL
		else {
			throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::SpatialObjectRepresentationBV(const SpatialObjectRepresentationBV&)");
		}
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}


//Destructor
SpatialObjectRepresentationBV::~SpatialObjectRepresentationBV()
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
			std::cerr << "Error in SpatialObjectRepresentationBV destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

// assignment operator
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator=(
    SpatialObjectRepresentationBV rhs)
{
	rhs.swap(*this);
	return *this;

}


int SpatialObjectRepresentationBV::getLabel() const
{return label;}

void SpatialObjectRepresentationBV::setLabel(int lab)
{label = lab;}

const BooleanValueMappedSPnode SpatialObjectRepresentationBV::getCopySubPaving() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"SpatialObjectRepresentationBV::getCopySubpaving()");
	}
	return BooleanValueMappedSPnode(*getSubPaving());
}

// get whether this has a subpaving.
bool SpatialObjectRepresentationBV::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}


cxsc::ivector SpatialObjectRepresentationBV::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"SpatialObjectRepresentationBV::getRootBox()");
	}
	return getSubPaving()->getBox();
}



int SpatialObjectRepresentationBV::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

cxsc::real SpatialObjectRepresentationBV::getDomainVolume() const
{
	real retValue(0.0);
	if (hasSubPaving()) {
		retValue = getSubPaving()->nodeRealVolume();
	}
	return retValue;
}
	
	
// Gets number of leaf nodes in the root paving.
size_t SpatialObjectRepresentationBV::getRootLeaves() const
{
	size_t result = 0;
	if (hasSubPaving()) {
		result = getSubPaving()->getNumberLeaves();
	}
	return result;
}

// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec SpatialObjectRepresentationBV::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// Get a string of the leaf node levels.
std::string SpatialObjectRepresentationBV::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}



//splits this according to string instruction
//returns true if some splitting was achieved
bool SpatialObjectRepresentationBV::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::splitToShape()");
	}
	bool success = false;
	BooleanValueMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"SpatialObjectRepresentationBV::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"SpatialObjectRepresentationBV::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);
		
		if (!success) {
			
			handleSPError(temp);
	   }
	   
	}
	catch (std::invalid_argument const& ia) {
		cerr << ia.what() << endl;
		handleSPError(temp);
		success = false;
	}
	catch (std::logic_error const& le) {
		cerr << le.what() << endl;
		handleSPError(temp);
		success = false;
	}
	return success;
	// any other exceptions are unhandled
}

void SpatialObjectRepresentationBV::allocateRanges(
				const std::vector< bool >& rangesToAllocate)
{
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::allocateRanges(...)");
	}
	vector< BooleanMappedValue > tmp;
	for (size_t i = 0; i < rangesToAllocate.size(); ++i)
		tmp.push_back(BooleanMappedValue(rangesToAllocate[i]));
		
	BooleanValueMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		getSubPaving()->allocateRanges(tmp);
	   
	}
	catch (std::invalid_argument const& ia) {
		cerr << ia.what() << endl;
		handleSPError(temp);
	}
	
	// any other exceptions are unhandled
}

void SpatialObjectRepresentationBV::allocateRanges(
				const std::vector< BooleanMappedValue >& rangesToAllocate)
{
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::allocateRanges(...)");
	}
	BooleanValueMappedSPnode temp(*getSubPaving()); // copy to temp
	try {
		getSubPaving()->allocateRanges(rangesToAllocate);
	   
	}
	catch (std::invalid_argument const& ia) {
		cerr << ia.what() << endl;
		handleSPError(temp);
	}
	
	// any other exceptions are unhandled
}

void SpatialObjectRepresentationBV::reshapeToUnion(
					const SpatialObjectRepresentationBV& other)
{
	if ( !hasSubPaving() || !other.hasSubPaving() ) {
		throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::reshapeToUnion(const SpatialObjectRepresentationBV&)");
	}
	
	getSubPaving()->reshapeToUnion(*other.getSubPaving());
}

SpatialObjectRepresentationBV SpatialObjectRepresentationBV::makeShapeToUnion(
					const SpatialObjectRepresentationBV& other) const
{
	if ( !hasSubPaving() || !other.hasSubPaving() ) {
		throw NullSubpavingPointer_Error(
				"SpatialObjectRepresentationBV::reshapeToUnion(const SpatialObjectRepresentationBV&)");
	}
	
	SpatialObjectRepresentationBV tmp(*this);
	tmp.reshapeToUnion(other);
	return tmp;
}


cxsc::real SpatialObjectRepresentationBV::getTotalVolume() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::getTotalVolume()");
	}
	return getSubPaving()->getTotalLeafTrueAreaOnBox();
}


/* Union to self operator.*/     
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator+= (
					const SpatialObjectRepresentationBV& rhs)
{
	if( !hasSubPaving() || !(rhs.hasSubPaving()) ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator+= (const SpatialObjectRepresentationBV&)");
	}
	
	getSubPaving()->operator+=( *(rhs.getSubPaving()) );
	return *this;
}

/* \brief Union operator.*/
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator+ (
					const SpatialObjectRepresentationBV& rhs) const
{
	SpatialObjectRepresentationBV result =(*this);

	result+= rhs;
	
	if ( label != rhs.getLabel() ) result.setLabel(0);
	
	return result;
}

/* Self-scalar OR operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator+= (
					bool val)
{
	if( !hasSubPaving() ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator+= (bool)");
	}
	
	getSubPaving()->operator+=( val );
	return *this;
}

/* Scalar OR operator*/
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator+ (
					bool val) const
{
	SpatialObjectRepresentationBV result =(*this);

	result+= val;
	
	return result;
}

/*XOR (symmetric set difference) against self operator.*/  
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator-= (
					const SpatialObjectRepresentationBV& rhs)
{
	
	if( !hasSubPaving() || !(rhs.hasSubPaving()) ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator-= (const SpatialObjectRepresentationBV&)");
	}
	
	getSubPaving()->operator-=( *(rhs.getSubPaving()) );
	return *this;
}

/* XOR (symmetric set difference).*/   
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator- (
					const SpatialObjectRepresentationBV& rhs) const
{
	SpatialObjectRepresentationBV result =(*this);

	result-= rhs;

	if ( label != rhs.getLabel() ) result.setLabel(0);
	
	return result;
}

/* Self-scalar XOR (symmetric set difference) operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator-= (
					bool val)
{
	if( !hasSubPaving() ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator-= (bool)");
	}
	
	getSubPaving()->operator-=( val );
	return *this;
}

/* Scalar XOR (symmetric set difference) operator*/
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator- (
					bool val) const
{
	SpatialObjectRepresentationBV result =(*this);

	result-= val;
	
	return result;
}

/* Intersection with self operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator*= (
					const SpatialObjectRepresentationBV& rhs)
{
	
	if( !hasSubPaving() || !(rhs.hasSubPaving()) ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator*= (const SpatialObjectRepresentationBV&)");
	}
	
	getSubPaving()->operator*=( *(rhs.getSubPaving()) );
	return *this;
}

/* Intersection operator.*/   
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator* (
					const SpatialObjectRepresentationBV& rhs) const
{
	SpatialObjectRepresentationBV result =(*this);

	result*= rhs;
	
	if ( label != rhs.getLabel() ) result.setLabel(0);
	
	return result;
}

/* Self-scalar AND operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator*= (
					bool val)
{
	if( !hasSubPaving() ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator*= (bool)");
	}
	
	getSubPaving()->operator*=( val );
	return *this;
}

/* Scalar AND operator*/
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator* (
					bool val) const
{
	SpatialObjectRepresentationBV result =(*this);

	result*= val;
	
	return result;
}


/* Set difference to self operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator/= (
					const SpatialObjectRepresentationBV& rhs)
{
	if( !hasSubPaving() || !(rhs.hasSubPaving()) ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator/= (const SpatialObjectRepresentationBV&)");
	}
	
	getSubPaving()->operator/=( *(rhs.getSubPaving()) );
	return *this;
}

/* Set difference operator.*/   
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator/ (
						const SpatialObjectRepresentationBV& rhs) const
{
	SpatialObjectRepresentationBV result =(*this);

	result/= rhs;

	if ( label != rhs.getLabel() ) result.setLabel(0);
	
	return result;
}

/*Self-scalar set difference operator.*/
SpatialObjectRepresentationBV& SpatialObjectRepresentationBV::operator/= (
						bool val)
{
	if( !hasSubPaving() ) {
		throw subpavings::NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::operator/= (bool val)");
	}
	
	
	getSubPaving()->operator/=( val );
	return *this;
}

/*Scalar set difference operator.*/
const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::operator/ (
						bool val) const
{
	SpatialObjectRepresentationBV result =(*this);

	result/= val;

	return result;
}

const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::makeSlice(
					const std::vector < int >& sliceDims,
					const std::vector < cxsc::real >& slicePts,
					const std::string sliceFilename) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::makeSlice(const std::vector < int >&, const std::vector < cxsc::real >&)");
	}
	BooleanValueMappedSPnode copyRoot(*getSubPaving());
	copyRoot.slice(sliceDims, slicePts, sliceFilename);
	return SpatialObjectRepresentationBV
			(copyRoot, getLabel());
}

const SpatialObjectRepresentationBV SpatialObjectRepresentationBV::makeSlice(
					const std::vector < int >& sliceDims,
					const std::vector < double >& slicePts,
					const std::string sliceFilename) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"SpatialObjectRepresentationBV::makeSlice(const std::vector < int >&, const std::vector < cxsc::real >&)");
	}
	std::vector < cxsc::real > slicePtsReal(slicePts.begin(), slicePts.end());
	
	BooleanValueMappedSPnode copyRoot(*getSubPaving());
	copyRoot.slice(sliceDims, slicePtsReal, sliceFilename);
	return SpatialObjectRepresentationBV
			(copyRoot, getLabel());
}


// get value mapped to leaf node of subpaving containing \a pt
bool SpatialObjectRepresentationBV::pointwiseExtension(
				const rvector& pt) const
{
	
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
					"SpatialObjectRepresentationBV::pointwiseExtension(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
			"SpatialObjectRepresentationBV::pointwiseExtension(const rvector&)");
	}
	bool result = false;
	
	const BooleanValueMappedSPnode * container = 
				getSubPaving()->findContainingNode(pt);
				
	if (container != NULL) {
		result = container->getRange();
		
	}
		
	return result;
}


// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void SpatialObjectRepresentationBV::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void SpatialObjectRepresentationBV::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the SpatialObjectRepresentationBV object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabsTrue(os, prec); // the output
			
		}
		if (confirm)
			std::cout << "The output of the SpatialObjectRepresentationBV "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}



void SpatialObjectRepresentationBV::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void SpatialObjectRepresentationBV::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the SpatialObjectRepresentationBV
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
			
		}
		if (confirm)
			std::cout << "Details of the root paving of the SpatialObjectRepresentationBV "
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
std::ostream & SpatialObjectRepresentationBV::outputRootToStreamTabs(
													std::ostream & os,
													int prec) const
{
	if (hasSubPaving()) {
		
		getSubPaving()->nodesAllOutput(os, 1, prec); // the output
		
	}
	
    return os;
}


// Method to add current state of this during splitting to a log file
// Output goes to file named according to argument s
// Output is textToTabs
void SpatialObjectRepresentationBV::outputLog(const std::string& s, 
											int i, int prec) const
{
    // To add output of the SpatialObjectRepresentationBV object to file
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



std::string SpatialObjectRepresentationBV::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "This label = " << label << endl;
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}


void SpatialObjectRepresentationBV::swap(SpatialObjectRepresentationBV& pcf) // throw()
{
	std::swap(rootPaving, pcf.rootPaving); // just swap the paving pointers
	std::swap(label, pcf.label); // and labels
}


// --------------------------- private ---------------------------------------


BooleanValueMappedSPnode* SpatialObjectRepresentationBV::getSubPaving() const
{return rootPaving;}


// Method to put opening line into a log file
void SpatialObjectRepresentationBV::outputLogStart(const std::string& s) const
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


//check that the box is okay
bool SpatialObjectRepresentationBV::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void SpatialObjectRepresentationBV::handleSPError(
											BooleanValueMappedSPnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
			std::cerr << "Operation aborted:original subpaving restored."
			<< std::endl;
}

// ensure rootPaving is deleted if constructed in failed constructor
void SpatialObjectRepresentationBV::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}


// ----------------------------- non member functions

//Output all boxes in SpatialObjectRepresentationBV pcf
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::SpatialObjectRepresentationBV& pcf)
{
    pcf.outputRootToStreamTabs(os);
    return os;
}


// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap (subpavings::SpatialObjectRepresentationBV & s1, 
		subpavings::SpatialObjectRepresentationBV & s2) // throw ()
{
	s1.swap(s2);
}







/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
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
\brief AdaptiveHistogramCollator definitions
*/

#include "adaptivehistogramcollator.hpp"
#include "subpaving_exception.hpp"

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"

//to use subpavings
#include "sptools.hpp"

// to use histogram penalty function objects
#include "histpenalty.hpp"

#include <string>   // to use the C++ string class
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <functional> // use functionals

#include <math.h> // math library


using namespace subpavings;
using namespace std;


// ---------- implementation of AdaptiveHistogramCollator class -------------

// ---------------- public methods

// default constructor
AdaptiveHistogramCollator::AdaptiveHistogramCollator()
				: rootCollator(NULL)
{
    try {
        rootCollator = new CollatorSPnode();
    }
	catch (exception const& e) {
		constructor_error_handler();
    }

}


// initialised constructor, initialised with an AdaptiveHistogram object
AdaptiveHistogramCollator::AdaptiveHistogramCollator(
	const AdaptiveHistogram& adh)
	: rootCollator(NULL)
{
    try {
		if (!adh.hasSubPaving()) {
			throw NullSubpavingPointer_Error(
			"AdaptiveHistogramCollator::AdaptiveHistogramCollator(const AdaptiveHistogram&)");
		}
		
		rootCollator = new CollatorSPnode(*adh.getSubPaving());
    }
	catch (exception const& e) {
		constructor_error_handler();
    }

}

// copy constructor
AdaptiveHistogramCollator::AdaptiveHistogramCollator(
                const AdaptiveHistogramCollator& other)
				: rootCollator(NULL)
{
    try {
		if (other.getSubPaving() == NULL) {
			// this should not be possible with current constructors
			throw NullSubpavingPointer_Error(
			"AdaptiveHistogramCollator::AdaptiveHistogramCollator(const AdaptiveHistogramCollator&)");
		}
		
        rootCollator = new CollatorSPnode(*(other.getSubPaving()));
    }
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// assignment operator
AdaptiveHistogramCollator& AdaptiveHistogramCollator::operator=(
    AdaptiveHistogramCollator rhs)
{
	rhs.swap(*this);
	return *this;

}


// Destructor.
AdaptiveHistogramCollator::~AdaptiveHistogramCollator()
{ 
	try {
		delete rootCollator;
		rootCollator = NULL;

	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in NewAdaptiveHistogram destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}


AdaptiveHistogramCollator AdaptiveHistogramCollator::importCollator(
					const std::string& s)
{
	//CollatorSPnode* newRoot = NULL;
	try {
		AdaptiveHistogramCollator newCollator;
		
		ifstream ifs(s.c_str());
		
		if ( ifs.is_open() ) {
			
			string nodeLevelsString;
			std::vector < string > strIntervals;
			size_t linesRead = 0;
		
			while(ifs.good() && nodeLevelsString.empty()) {
			
				// read first line - the node levels
				getline (ifs, nodeLevelsString);
				linesRead++;
			}
		
			
			string line;
					
			bool readingBox = true;
			
			string whiteSpace(" \t\r\n");
			string boxStart("[");
			string boxEnd("[");
			while (readingBox && ifs.good()) {

				getline(ifs, line);
				linesRead++;
				if (line.empty()) continue;
				
				size_t start1 = line.find_first_not_of(whiteSpace);
				readingBox = ( (start1 < string::npos)  
							&& (line.find_first_of(boxStart, start1) == 0));
				if (readingBox) {
					readingBox = (line.find_first_of(boxEnd, start1) < string::npos);
				}
				if (readingBox) strIntervals.push_back(line);
			}
			
			if (ifs.good()) linesRead--; // last line read was not a box line
			bool hasRanges = !ifs.eof();
					
			ifs.close();
			// end reading file
			
			if(!nodeLevelsString.empty()) {
			
				if (strIntervals.empty() ) {
					throw std::runtime_error("No box");
				}
			
				int dims = strIntervals.size();
			
				cxsc::ivector box(dims);
				
				for (int i = 1; i <= dims; ++i) {

					strIntervals[i-1] >> box[i];
					
				}
				
				CollatorSPnode newRoot(box);
				
				bool success = newRoot.splitRootToShape(nodeLevelsString);
				
				if (!success) {
					throw std::runtime_error("Could not split root to required shape");
				}
				
				if (hasRanges) {
					vector < RealVec > ranges;
					bool allocate = readVecRealVecFromTxt(ranges, s, linesRead);
					
					if (allocate) newRoot.allocateRanges(ranges);
				}
				
				std::swap(*(newCollator.rootCollator), newRoot);
				
				//newCollator.rootCollator = newRoot;
			}
			// no nodeLevelsString means collator exported was empty
		}
		else {
			
			throw subpavings::IO_Error("Could not open file named " + s);
		}
		return newCollator;

	}
	catch (std::exception const& e) {
		//delete newRoot;
		
		throw subpavings::UnfulfillableRequest_Error(
		"AdaptiveHistogramCollator::importCollator(const std::string&):\n"
						+ string(e.what()) );
	}
}




// Return a pointer to the CollatorPSnode this manages
CollatorSPnode* AdaptiveHistogramCollator::getSubPaving() const
{
	return rootCollator;
} // boost::shared_ptr might be better


bool AdaptiveHistogramCollator::isEmptyCollation() const
{
	return getSubPaving()->isEmptyRangeCollection();
}

cxsc::ivector AdaptiveHistogramCollator::getRootBox() const
{
	if (getSubPaving() == NULL) {
		throw NullSubpavingPointer_Error(
						"AdaptiveHistogramCollator::::getRootBox()");
	}
	if (getSubPaving()->isEmpty()) {
		throw NoBox_Error(
						"AdaptiveHistogramCollator::::getRootBox()");
	}
	return getSubPaving()->getBox();
}

int AdaptiveHistogramCollator::getDimensions() const
{
	int retValue = 0;
	if ((getSubPaving() != NULL) && (!getSubPaving()->isEmpty())) {
		
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

// increment addition operator
AdaptiveHistogramCollator& AdaptiveHistogramCollator::operator+=(
		const AdaptiveHistogramCollator& rhs)
{
	//nothing to add
	if ( rhs.isEmptyCollation() ) {
		
		return *this;
	}
			
	// if this has no subpaving or an empty one, it should just copy the other one
	if ( isEmptyCollation() ) {
		
		*this = rhs;
		return *this;
	}

	// get here only if both have something in their collations
	
		getSubPaving()->addPaving( rhs.getSubPaving() );
		return *this;
}

// addition operator
const AdaptiveHistogramCollator AdaptiveHistogramCollator::operator+(const
    AdaptiveHistogramCollator& rhs) const
{
    AdaptiveHistogramCollator temp(*this);
			
	temp += rhs;
	
	return temp;

}

// increment addition operator
AdaptiveHistogramCollator& AdaptiveHistogramCollator::operator+=(
		const AdaptiveHistogram& rhs)
{

	//nothing to add
	if ( !rhs.hasSubPaving() || rhs.getSubPaving()->isEmpty() ) {
		
		return *this;
	}
			
	// if this has nothing collated
	if ( isEmptyCollation() ) {
		
		*this = AdaptiveHistogramCollator(rhs);
		return *this;
	}

	// get here only if both have subpavings
	
	getSubPaving()->addPaving( rhs.getSubPaving() );
	return *this;
}

// addition operator
const AdaptiveHistogramCollator AdaptiveHistogramCollator::operator+(
	const AdaptiveHistogram& rhs) const
{
	AdaptiveHistogramCollator temp(*this);
			
	temp += rhs;
	
	return temp;
}



// averaging method
const AdaptiveHistogramCollator AdaptiveHistogramCollator::makeAverage() const
{
	
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error("AdaptiveHistogramCollator::makeAverage()");
	}

	AdaptiveHistogramCollator temp(*this);

	temp._average();

	return temp;
}


// normalising method
const AdaptiveHistogramCollator AdaptiveHistogramCollator::makeNormalised() const
{
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error(
			"AdaptiveHistogramCollator::makeNormalised()");
	}
	
	AdaptiveHistogramCollator temp(*this);
	
	temp._normalise();
	return temp;
}

const AdaptiveHistogramCollator AdaptiveHistogramCollator::makeMarginal(
								const std::vector<int>& reqDims) const
{
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error(
			"AdaptiveHistogramCollator::makeMarginal(const std::vector<int>&)");
	}
	
	AdaptiveHistogramCollator temp(*this);
	
	temp._marginalise(reqDims);
	return temp;
}

// do checks and use private coverage method to find coverage
double AdaptiveHistogramCollator::findCoverage(const rvector& pt) const
{
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error("AdaptiveHistogramCollator::findCoverage(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error("AdaptiveHistogramCollator::findCoverage(const rvector&)");
	}
	
	return _coverage(pt);
}


//new AHABC
//find nodes constituting density region with cov coverage
std::vector< const CollatorSPnode* > & AdaptiveHistogramCollator::findDensityRegion(
								std::vector< const subpavings::CollatorSPnode* > & covNodes,
								double cov) const
{
	vector < const CollatorSPnode* > tmp;
	
	// put the leaves into a vector and sort it, smallest to largest
	std::vector< const CollatorSPnode * > leaves;
	getSubPaving()->getConstLeaves(leaves);
	
	sort(leaves.begin(), leaves.end(), nodeCompTotalRangeCollection);
			
	//start iterating from the largest
	vector< const CollatorSPnode* >::const_reverse_iterator rit = leaves.rbegin();
	bool found = false; //found the boxes that gives cov density region
	
	real totalArea = 
				getSubPaving()->getTotalAbsValueTimesVol();
	// the total area of this might not be 1
	real thisCov = totalArea * cov;
	
	cxsc::dotprecision totalCov(0.0);
	
	while (!found && rit < leaves.rend()) {
		// height is box counts/box vol
		// ie count is unnormalised vol of an individual element of histogram
		// box vol * height == box vol * (count / box vol) == count
	
		//accumulate the summary * box vol
		accumulate(totalCov, (*rit)->getTotalAbsValueTimesVol(), 1); 

		//push back the node that fulfill the condition totalCov <= cov 
		//into the container covNodes
		if (!(thisCov < totalCov)) { 
			tmp.push_back((*rit)); 
		} 
		
		// found is true if totalCov >= thisCov
		found = !(totalCov < thisCov) ;
		
		++rit;				
	} // end while
	
	tmp.swap(covNodes);
	return covNodes; 
} 

// do checks and use private density method to find density
double AdaptiveHistogramCollator::findEmpiricalDensity(const rvector& pt) const
{
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error(
			"AdaptiveHistogramCollator::findEmpiricalDensity(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
			"AdaptiveHistogramCollator::findEmpiricalDensity(const rvector&)");
	}
	
	return _empiricalDensity(pt);
}

//new AHABC
bool AdaptiveHistogramCollator::histCollatorBoxContains(const rvector& pt) const
{
	if (isEmptyCollation()) {
		throw UnfulfillableRequest_Error(
				"AdaptiveHistogram::histCollatorBoxContains(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
			"AdaptiveHistogramCollator::histCollatorBoxContains(const rvector&)");
	}
	
	return ( getSubPaving()->findContainingNode(pt) != NULL ); 
}

RealVec AdaptiveHistogramCollator::getL1DistancesToAverage() const
{
	RealVec rvec;
	rvec = getL1DistancesToAverage(rvec);
	return rvec;
}

// take a container and return the same container, which has been
// cleared (if necessary) and re-filled with 
// L1-distances-to-average, one for each histogram in collation
RealVec& AdaptiveHistogramCollator::getL1DistancesToAverage(RealVec& container) const
{
	if (isEmptyCollation()) {
			throw UnfulfillableRequest_Error(
			"AdaptiveHistogramCollator::getL1DistancesToAverage(RealVec&)");
	}

	container = getSubPaving()->getL1DistancesToAverage(container);
   
	return container;

}

RealVec AdaptiveHistogramCollator::getL1DistancesToAverage(
					const AdaptiveHistogramCollator& other) const
{
	RealVec rvec;
	rvec = getL1DistancesToAverage(rvec, other);
	return rvec;
}

// take a container and return the same container, which has been
// cleared (if necessary) and re-filled with 
// L1-distances to average-of-other, one for each histogram in collation
RealVec& AdaptiveHistogramCollator::getL1DistancesToAverage(
						RealVec& container,
						const AdaptiveHistogramCollator& other) const
{
	if (NULL == getSubPaving()) {
		// should never happen
		throw NullSubpavingPointer_Error(
		"AdaptiveHistogramCollator::getL1DistancesToAverage(RealVec&, const AdaptiveHistogramCollator&)");
	}

	if (other.isEmptyCollation()) {
		throw UnfulfillableRequest_Error(
		"AdaptiveHistogramCollator::getL1DistancesToAverage(RealVec&, const AdaptiveHistogramCollator&)");
	}
	
	container = getSubPaving()->getL1DistancesToAverage(container, 
												other.getSubPaving());
   
	return container;
}

RealVec AdaptiveHistogramCollator::getL1Distances(
					const AdaptiveHistogram& adh) const
{
	RealVec rvec;
	rvec = getL1Distances(rvec, adh);
	return rvec;
}

// take a container and return the same container, which has been
// cleared (if necessary) and re-filled with 
// L1-distances to the adaptive histogram, one for each histogram in collation
RealVec& AdaptiveHistogramCollator::getL1Distances(
						RealVec& container,
						const AdaptiveHistogram& adh) const
{
	if (NULL == getSubPaving()) {
		// should never happen
		throw NullSubpavingPointer_Error(
		"AdaptiveHistogramCollator::getL1Distances(RealVec&, const AdaptiveHistogram&)");
	}

	container = getSubPaving()->getL1Distances(container, 
												adh.getSubPaving());
   
	return container;

}

//splits collator histogram according to string instruction
//returns true if some splitting was achieved
//used for testing
bool AdaptiveHistogramCollator::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (NULL == getSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogramCollator::splitToShape()");
	}
	bool success = false;
	CollatorSPnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"AdaptiveHistogramCollator::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"AdaptiveHistogramCollator::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);

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


// Add an AdaptiveHistogram into this collation
void AdaptiveHistogramCollator::addToCollation(const AdaptiveHistogram& adh)
{
    
	(*this) += adh;
	
}

void AdaptiveHistogramCollator::addToCollation(
					const std::vector < AdaptiveHistogram >& samples) 
{
	for (std::vector < AdaptiveHistogram >::const_iterator it = samples.begin();
				it < samples.end();
				++it) {
		(*this) += (*it);
	}
}

// make a .dot file for the histogram
bool AdaptiveHistogramCollator::outputGraphDot() const
{
	bool success = false;
	
    if (isEmptyCollation()) {

        std::cerr << "Sorry, you can't make a graph with nothing collated"
                << std::endl;
    }

    else success = getSubPaving()->outputGraphDot();

    return success;
}

// Get the number of Adaptive Histogram objects collated.
size_t AdaptiveHistogramCollator::getNumberCollated() const
{
	return getSubPaving()->getSizeRangeCollection();
}

// Output the collated normalised histogram heights and bins data to a txt file
void AdaptiveHistogramCollator::outputToTxtTabs(const std::string& s,
                                        int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Output the collated normalised histogram heights and bins data to a txt file
void AdaptiveHistogramCollator::outputToTxtTabs(const std::string& s,
                                        int prec, bool confirm) const
{
    // To generate a file output of the AdaptiveHistogramCollator object
    ofstream os(s.c_str());         // Filename, c-string version

    getSubPaving()->leavesOutputTabs(os, prec);
    if (confirm)
        std::cout << "The output of the AdaptiveHistogramCollator has been "
            << "written to " << s << std::endl << std::endl;

}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
std::ostream & AdaptiveHistogramCollator::outputToStreamTabs(
						std::ostream & os, int prec) const
{
	if (NULL != getSubPaving()) {

		 return getSubPaving()->leavesOutputTabs(os, prec); // the output
	}
	
    else return os;
}

void AdaptiveHistogramCollator::outputAverageToTxtTabs(const
    std::string& s, int prec) const
{
	outputAverageToTxtTabs(s, prec, false);
}


// Output the average data over the collation to a txt file
// this outputs the normalised average histogram heights and bins
void AdaptiveHistogramCollator::outputAverageToTxtTabs(const
    std::string& s, int prec, bool confirm) const
{
    if (NULL != getSubPaving()) {

		// To generate a file output of the AdaptiveHistogramCollator object
		ofstream os(s.c_str());         // Filename, c-string version

		if (NULL != getSubPaving() && !isEmptyCollation()) {
			getSubPaving()->leavesAverageOutputTabs(os, prec);
			if (confirm)
				std::cout << "The output of the average AdaptiveHistogram has been "
					<< "written to " << s << std::endl << std::endl;
		}
		else {
			std::cerr << "Sorry, nothing is in collation to average"
				<< std::endl;
		}
	}
}


//new
// Method to output a description of this collator that can be read
// in again to remake it
void AdaptiveHistogramCollator::exportCollator(
					const std::string& s, 
					int prec) const
{
	ofstream os(s.c_str());
	if (os.is_open()) {
	
		if (!isEmptyCollation() ) {
			vector <RealVec> ranges;
			
			ranges = getSubPaving()->getAllRangeCollections(ranges);
			
			string leafNodeLevels = getSubPaving()->getLeafNodeLevelsString();
			
			ivector box = getRootBox();
		
			os << leafNodeLevels << endl;
			
			os << cxsc::SetPrecision(prec+2, prec);
			
			os << box << endl;
			
			for (vector <RealVec>::const_iterator it = ranges.begin();
					it < ranges.end();
					++it) {
				os << (*it) << endl;
			}
		}
		os.close();
		
	}
	else {
			std::cerr << "Error: could not open file named "
				<< s << std::endl << std::endl;
	}
}

// Method to add current state of the histogram collator to a log file
// Output goes to file named according to argument s
void AdaptiveHistogramCollator::outputLog(const std::string& s, 
					const int i, int prec) const
{
    // To add output of the AdaptiveHistogramCollator object to file
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



void AdaptiveHistogramCollator::swap(AdaptiveHistogramCollator& adh) // throw()
{
	std::swap(rootCollator, adh.rootCollator); // just swap the paving pointers
}

// private methods



// average this
// assumes a non-null rootCollator
void AdaptiveHistogramCollator::_average()
{
	getSubPaving()->average();
}

// normalise this
// assumes a non-null rootCollator
void AdaptiveHistogramCollator::_normalise()
{
	getSubPaving()->normalise();
}

//marginalise this
// assumes a non-null rootCollator
void AdaptiveHistogramCollator::_marginalise(
						const std::vector<int>& reqDims)
{
	getSubPaving()->marginalise(reqDims);
}

// internal method to find coverage
double AdaptiveHistogramCollator::_coverage(const rvector& pt) const
{
	double cov = 0.0;
	
	real totalArea = 
				getSubPaving()->getTotalAbsValueTimesVol();
	
	//if total value times vol is 0, coverage will always be 0
	if (totalArea > cxsc::real(0.0) ) {
	
		const CollatorSPnode * container = 
						getSubPaving()->findContainingNode(pt);

		if (container != NULL) {
			
			real culmArea = totalArea;
		
			// put the leaves into a vector and sort it, smallest to largest
			std::vector< const CollatorSPnode * > leaves;
			getSubPaving()->getConstLeaves(leaves);
			
			sort(leaves.begin(), leaves.end(), nodeCompTotalRangeCollection);
			
			bool found = FALSE;
			cxsc::real containerHeightNonNorm = container -> getTotalRangeCollection();
			
			std::vector< const CollatorSPnode * >::const_reverse_iterator
										rit = leaves.rbegin();
			
			while (!found && rit < leaves.rend()) {
				
				// height is box counts/box vol
				// ie count is unnormalised vol of an individual element of histogram
				// box vol * height == box vol * (count / box vol) == count
			
				// check the boxes
				// stop at the first box not taller than this one
				cxsc::real thisHeightNonNorm = (*rit)->getTotalRangeCollection();
				
				if ( thisHeightNonNorm > containerHeightNonNorm ) {
					// decrement cumulative counts (unnormalised histogram volume)
					culmArea -= (*rit)->getTotalAbsValueTimesVol();
				}
				
				/*
				if ( fabs( thisHeightNonNorm 
						- containerHeightNonNorm ) > 
						DBL_EPSILON*cxsc::max(thisHeightNonNorm, containerHeightNonNorm) ) {
					// decrement cumulative counts (unnormalised histgram volume)
					culmArea -= (*rit)->getTotalAbsValueTimesVol();
					
				}
				*/
				
				else { 
					found = TRUE;	// break out of loop
				}
				++rit;
				
			} // end while
				
			// if we have not found we have a problem,
			// since findContainingNode said it was here somewhere
		
			if (!found) {
				throw std::logic_error(
				"AdaptiveHistogramCollator::_coverage(const rvector&) : lost container");
			}	
			
			cov = cxsc::_double(culmArea/totalArea);
		}
	}
	return cov;
}
 
// internal method to find empiricalDensity
double AdaptiveHistogramCollator::_empiricalDensity(const rvector& pt) const
{
	double den = 0.0;
	
	real totalArea = getSubPaving()->getTotalAbsValueTimesVol();
	// this would be 1 if the thing was normalised, but it probably isn't
	
	// if the total area is 0, density will always be 0
	if (totalArea > cxsc::real(0.0) ) {			
		
		const CollatorSPnode * container = 
				getSubPaving()->findContainingNode(pt);
				
		if (container != NULL) {
				
			cxsc::real total = 	container->getTotalRangeCollection();
		
			if (total < 0.0) {
				throw ( std::logic_error(
				"AdaptiveHistogramCollator::_empiricalDensity(const rvector&) : Negative density") );
			}
			
			den = cxsc::_double(total/totalArea);
		}
	}
	return den;
}

void AdaptiveHistogramCollator::handleSplitToShapeError(CollatorSPnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
		std::cerr << "Your instruction does not describe a proper tree.";
		std::cerr << "  Please check your instruction and try again."
		<< std::endl;
}

// ensure rootCollator is deleted if constructed in failed constructor
void AdaptiveHistogramCollator::constructor_error_handler() 
{
	try {
		
			delete rootCollator;
			rootCollator = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}


// --------------- end private methods


// ---------- end implementation of AdaptiveHistogramCollators -----------

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap (subpavings::AdaptiveHistogramCollator & a1, 
		subpavings::AdaptiveHistogramCollator & a2) // throw ()
{
	a1.swap(a2);
}


//Output all boxes in AdaptiveHistogramCollator adhc
std::ostream & subpavings::operator<<(std::ostream &os,
                    const AdaptiveHistogramCollator& adhc)
{
    if (NULL != adhc.getSubPaving()) {
        adhc.getSubPaving()->nodesAllOutput(os, 1);
		os << std::endl;
    }

    return os;
}





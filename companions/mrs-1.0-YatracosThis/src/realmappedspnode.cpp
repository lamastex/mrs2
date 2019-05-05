/*
* Copyright (C) 2010, 2011, 2012 Jennifer Harlow
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

/*!/ \file
\brief RealMappedSPnode definitions.
*/

#include "realmappedspnode.hpp"

#include <gsl/gsl_math.h> //gsl_isnan(), isinf

//#define MARG_OUTPUT // for console output from marginalisation
//#define SLICE_OUTPUT // for console output from slice
//#define DEBUG_L1

using namespace subpavings;
using namespace std;


// ------------------------ public member functions ----------------

// Destructor.
RealMappedSPnode::~RealMappedSPnode()  {}


// no-argument constructor,
RealMappedSPnode::RealMappedSPnode()
	: MappedSPnode<cxsc::real>() {}   // uses the base SPnode class default constructor


// initialised constructor
// initialised with a box
RealMappedSPnode::RealMappedSPnode(
						const ivector& v)
	: MappedSPnode<cxsc::real>(v, real(0.0)) {}


 // initialised constructor
// initialised with a labeled box
RealMappedSPnode::RealMappedSPnode(
						const LabBox& lb)
	: MappedSPnode<cxsc::real>(lb, real(0.0)) {}



// initialised constructor
// initialised with a box, a value for the range,
RealMappedSPnode::RealMappedSPnode(
						const ivector& v, const cxsc::real& range)
	: MappedSPnode<cxsc::real>(v, range) {}



// initialised constructor
// initialised with a labeled box and a value for the range
RealMappedSPnode::RealMappedSPnode(
						const LabBox& lb, const cxsc::real& range)
	: MappedSPnode<cxsc::real>(lb, range) {}


// Copy constructor
// copies from given node downwards
RealMappedSPnode::RealMappedSPnode(
						const SPnode& spn)
{
	if (!spn.isEmpty()) { 
		theBox = new ivector( spn.getBox() );
	}
	nodeName = spn.getNodeName();
	range = 0.0;
	
	//recursion on the children
	if (spn.hasLCwithBox()) {
		nodeAddLeft(new RealMappedSPnode(
			*(spn.getLeftChild())));
	}
	else leftChild=NULL;

	if (spn.hasRCwithBox()) {
		nodeAddRight(new RealMappedSPnode(
			*(spn.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
// adjusts for volume but does not "normalise" 
//ie does not divide by total counter (which it does not know)
RealMappedSPnode::RealMappedSPnode(
						const SPSnode& spn)
{
	if (!spn.isEmpty()) { 
		theBox = new ivector( spn.getBox() );
	}
	nodeName = spn.getNodeName();
	range = (1.0*spn.getCounter()/nodeRealVolume());
	
	//recursion on the children
	if (spn.hasLCwithBox()) {
		nodeAddLeft(new RealMappedSPnode(
			*(spn.getLeftChild())));
	}
	else leftChild=NULL;

	if (spn.hasRCwithBox()) {
		nodeAddRight(new RealMappedSPnode(
			*(spn.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
RealMappedSPnode::RealMappedSPnode(
						const RealMappedSPnode& other)
{
	
	if (other.theBox != NULL) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new RealMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new RealMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
RealMappedSPnode::RealMappedSPnode(
						const MappedSPnode<cxsc::real>& other)
{
	if (!other.isEmpty()) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new RealMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new RealMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;
}


// Copy assignment operator
// copies from given node downwards
RealMappedSPnode& RealMappedSPnode::operator=(
						RealMappedSPnode rhs)
{
	rhs.swapRMSPSR(*this); // make sure we use our version of swap
	return(*this);
}

// Copy assignment operator
// copies from given node downwards
RealMappedSPnode& RealMappedSPnode::operator=(
						MappedSPnode<cxsc::real> rhs)
{
	rhs.swapMSPSR(*this); // make sure we use our version of swap
	return(*this);
}



// parent and child accessors have to hide the base class implementation
// this is not good but otherwise we get the base class return type
// I've asked around and I can't find a way around it ...

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
RealMappedSPnode* RealMappedSPnode::getParent() const
{ return (RealMappedSPnode*) parent; }


// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node, cast to this type
RealMappedSPnode* RealMappedSPnode::getLeftChild() const
{ return (RealMappedSPnode*) leftChild; }


// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node, cast this type
RealMappedSPnode* RealMappedSPnode::getRightChild() const
{ return (RealMappedSPnode*) rightChild; }

// return true if there is a negative range in the tree rooted at this.
bool RealMappedSPnode::hasNegativeRangeInTree() const
{
	bool retValue = false;
	if (range < real(0.0)) {
		retValue = true;
	}
	else if (!isLeaf()) {
		retValue = getLeftChild()->hasNegativeRangeInTree();
		if (!retValue) retValue = getRightChild()->hasNegativeRangeInTree();
	}
	return retValue;
	
}

// return true if there is an infinite range in the tree rooted at this.
bool RealMappedSPnode::hasInfiniteRangeInTree() const
{
	bool retValue = false;
	if (range == cxsc::Infinity) {
		retValue = true;
	}
	else if (!isLeaf()) {
		retValue = getLeftChild()->hasInfiniteRangeInTree();
		if (!retValue) retValue = getRightChild()->hasInfiniteRangeInTree();
	}
	return retValue;
	
}

bool RealMappedSPnode::operator<(
		const RealMappedSPnode& rhs) const
{
	return (getRange() < rhs.getRange());
}

// Find where data should be 
// childInd is an indicator for which child is being checked
// throws exception if node has no box
const RealMappedSPnode* RealMappedSPnode::findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd) const
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"RealMappedSPnode::findContainingNode(const cxsc::rvector&, OPERATIONS_ON)");
	}
	// start at the top
	const RealMappedSPnode* retObj = NULL;
	
	if(nodeContains(pt, childInd)) {
		
		if(isLeaf()) {

			// give this node as return value
			retObj = this;

		} // end of isLeaf

		// if not a leaf and contains data
		// recurse on the children if any
		else {

			if(hasRCwithBox()){
				
				retObj =
				(getRightChild())->findContainingNode(
					pt, ON_RIGHT);
				
			}
			// only try left if we did not find on the right
			if(retObj == NULL && hasLCwithBox()) {
				retObj =
				(getLeftChild())->findContainingNode(
					pt, ON_LEFT);
				
			}
		}

	} // end if node contains

	// will return null if does not contain the data
	
	return retObj;
}


// add two sibling nodes to this provided this is a leaf
// comp argument is passed to Upper() and Lower()
// split the box in half normal to dimension set by comp
void RealMappedSPnode::nodeExpand(int comp)
{
	// can only expand if there is a box
	if (NULL == theBox) {
		throw NoBox_Error("RealMappedSPnode::nodeExpand(int)");
	}

	// only do something if this node is a leaf
	if (isLeaf()) {
		
		RealMappedSPnode* newLC = NULL;
		RealMappedSPnode* newRC = NULL;
		
		try {

			// ivectors to become boxes for new children
			ivector lC, rC;
			// Call Lower() and Upper() to put split boxes
			// into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);

			// make and add the new children
			newLC = new RealMappedSPnode(lC, range);
			newRC = new RealMappedSPnode(rC, range);
			
			nodeAddLeft(newLC);
			nodeAddRight(newRC);
			// both children get the same range as this
			
			//name the new children
			getLeftChild()->setNodeName(nodeName + "L");
			getRightChild()->setNodeName(nodeName + "R");

			// new children have summary from this
		}
		catch(std::exception& e) {
		// overkill, but try to deal with all eventualities...
			try {
				if (getLeftChild() != NULL) {
					delete (getLeftChild());
					leftChild = NULL;
				}
				
				if (getRightChild() != NULL) {
					delete (getRightChild());
					rightChild = NULL;
				}
				if (newLC != NULL) {
					delete newLC;
					newLC = NULL;
				}
				if (newRC != NULL) {
					delete newRC;
					newRC = NULL;
				}
			}
			catch(std::exception& e) {} //catch and swallow
			
			throw; // rethrow original exception
		}
	}
}



// add two sibling nodes to this provided this is a leaf
// finds its own comp argument then calls nodeExpand(int comp)
void RealMappedSPnode::nodeExpand()
{
	int maxdiamcomp; // variable to hold first longest dimension
	double temp = ::MaxDiam(getBox(), maxdiamcomp);
	nodeExpand(maxdiamcomp); // complete nodeExpand
}

// maximum range in the tree (range of this if this is a leaf)
cxsc::real RealMappedSPnode::getMaxRangeForLeavesInTree() const
{
	if (isLeaf()) {
		return getRange();
	}
	else {
		cxsc::real rLC = getLeftChild()->getMaxRangeForLeavesInTree();
		cxsc::real rRC = getRightChild()->getMaxRangeForLeavesInTree();
		return (rLC > rRC? rLC : rRC);
	}
}


// Marginalise
// marginalise from given node downwards
void RealMappedSPnode::marginalise(
		const std::vector<int>& reqDims)
{
	if ( getParent() != NULL ) {
		throw NonRootNode_Error(
			"RealMappedSPnode::marginalise(const std::vector<int>&)");
	}
	
	_start_marginalise(reqDims);
}

const RealMappedSPnode RealMappedSPnode::makeMarginalised(
							const std::vector<int>& reqDims) const
{
	RealMappedSPnode result = *this;

	result._start_marginalise(reqDims);
	
	return result;
}

void RealMappedSPnode::normalise()
{
	if ( getParent() != NULL ) {
		throw NonRootNode_Error(
			"RealMappedSPnode::normalise()");
	}
	_normalise();
}
		
const RealMappedSPnode RealMappedSPnode::makeNormalised() const
{
	RealMappedSPnode result = *this;

	result._normalise();
	
	return result;
}

//20160904
std::pair<size_t, cxsc::real> RealMappedSPnode::getNonZeroBoxSummary() const
{
	if(isEmpty()) 
		throw NoBox_Error("RealMappedSPnode::getNonZeroBoxSummary");
		
	size_t nNonZeroBoxes = 0;
	real vNonZeroBoxVolumes(0.0);
	
	accumulateNonZeroBoxSummary(nNonZeroBoxes, vNonZeroBoxVolumes);
	
	return std::pair <size_t, real >(nNonZeroBoxes, 
								vNonZeroBoxVolumes/nodeRealVolume());
	
}

void RealMappedSPnode::smearRanges(cxsc::real zeroRange,
							cxsc::real ratioRange)
{
	if(parent != NULL) 
		throw NonRootNode_Error("RealMappedSPnode::smearRanges(cxsc::real, cxsc::real)");
	
	if(isEmpty()) 
		throw NoBox_Error("RealMappedSPnode::smearRanges(cxsc::real, cxsc::real)");
		
	if(!(zeroRange > 0.0))  
		throw std::invalid_argument("RealMappedSPnode::smearRanges(cxsc::real, cxsc::real)");
	
	if(!((ratioRange > 0.0) && (ratioRange < 1.0)) ) 
		throw std::invalid_argument("RealMappedSPnode::smearRanges(cxsc::real, cxsc::real)");
	
	_smearRanges(zeroRange, ratioRange);

}


//accumluate number of non-empty boxes and total non-empty box vol
void RealMappedSPnode::accumulateNonZeroBoxSummary(size_t& nNonZeroBoxes, 
								real& vNonZeroBoxVolumes) const
{
	if (isLeaf()) {
		if ((range > 0.0) || (range < 0.0)) {
			++nNonZeroBoxes;
			vNonZeroBoxVolumes+=nodeRealVolume();
			
		}
	}
	else { // recurse
		getLeftChild()->accumulateNonZeroBoxSummary(nNonZeroBoxes, 
												vNonZeroBoxVolumes);
		getRightChild()->accumulateNonZeroBoxSummary(nNonZeroBoxes, 
												vNonZeroBoxVolumes);
	}
	// we want to count all empty leaves even if parent is empty too
}

void RealMappedSPnode::_smearRanges(cxsc::real zeroRange,
									cxsc::real ratioRange)
{
	if (isLeaf()) {
		if ((range > 0.0) || (range < 0.0)) {
			range = range*(1.0-ratioRange);
		}
		else range = zeroRange;
	}
	else {
		
		//recurse first
		real childSumBefore = getLeftChild()->getRange() + 
								getRightChild()->getRange();
		getLeftChild()->_smearRanges(zeroRange, ratioRange);
		getRightChild()->_smearRanges(zeroRange, ratioRange);
		range *= (getLeftChild()->getRange() + 
						getRightChild()->getRange())/childSumBefore;
	}
}



// overide base class
void RealMappedSPnode::slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts)
{
	
	#ifdef SLICE_OUTPUT
		std::cout << "In RealMappedSPnode::slice, I am " << getNodeName() << std::endl;
	#endif

	std::vector<cxsc::real> fullSlicePts = sliceCheck(sliceDims, slicePts);
	
	RealMappedSPnode temp(*this);
	temp._slice(sliceDims, fullSlicePts);
	swapRMSPSR(temp);

}


const RealMappedSPnode RealMappedSPnode::makeSlice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts) const
{
	
	std::vector<cxsc::real> fullSlicePts = sliceCheck(sliceDims, slicePts);
	
	RealMappedSPnode result(*this);
	result._slice(sliceDims, fullSlicePts);
	return result;

}

// get distance between two pavings
// checks boxes match
// throws exception if this is empty, or if rhs empty
// if both have no data, distance is 0
// if one has data and the other no data, distance is integral over the 
// one with data
cxsc::real RealMappedSPnode::getL1Distance(
					const RealMappedSPnode& other) const
{
	if ( isEmpty() || other.isEmpty() ) {
		throw NoBox_Error(
			"RealMappedSPnode::getL1Distance(const RealMappedSPnode&)");
	}
	
	//std::cout << getBox() << "\t" << other.getBox() << std::endl;
	
	if ( getBox() != other.getBox() ) {
		throw IncompatibleDimensions_Error(
			"RealMappedSPnode::getL1Distance(const RealMappedSPnode&)");
	}
	
	
	
	cxsc::dotprecision retDP(0.0);
	
	if (this != & other) retDP = _getL1distance(retDP, &other);
	
	return cxsc::rnd(retDP);
		
}

// NEW JUNE 2012 for log posteriors
cxsc::real RealMappedSPnode::getLogLikelihood(const SPSnode& spn) const
{
	if ( isEmpty() || spn.isEmpty() ) {
		throw NoBox_Error(
			"RealMappedSPnode::getLogLikelihood(const SPSnode&)");
	}
	
	if ( getBox() != spn.getBox() ) {
		throw IncompatibleDimensions_Error(
			"RealMappedSPnode::getLogLikelihood(const SPSnode&)");
	}
	cxsc::real result(0.0);
	
	if (spn.getCounter() > 0) {
		
		cxsc::dotprecision loglik(0.0);
		
		int isnan = 0;
		int isposinf =0;
		int isneginf =0; // all these 3 changed by ref
		
		_getLogLikelihood(loglik, isnan, isposinf, isneginf, &spn);
		
		/*nan if nan flag is set*/
		if (isnan )	result = cxsc::SignalingNaN;
		else result = cxsc::rnd(loglik);
		/* result could be -inf if there are 0 ranges in tree where spn 
		 * has points or inf if there are inf ranges in tree where spn has points*/ 
	}

	return result;
}
	
//-------------------	
/* Get the the "area" of the range and the box of an element
of a Scheffe set based on a given box*/ 
cxsc::real RealMappedSPnode::getIntegralForScheffeElement
																	(ivector& box, cxsc::real vol, bool split) const
{
	cxsc::real result = 0.0;
	
	//cout << box << "vs" << getBox() << endl;

	if ( isLeaf() ) {
		if (box <= getBox() && split == false ) {
			result += getRange() * vol;
		}
		else if ( isLeaf() && split == true ) {
			result += getRange() * nodeRealVolume();
		}
		return result;
	}  	
	else if ( box == getBox() ) {
		split = true;
		return (
		getLeftChild()->getIntegralForScheffeElement(box, vol, split)
		+ getRightChild()->getIntegralForScheffeElement(box, vol, split));			
	}	
	else {
		return (
		getLeftChild()->getIntegralForScheffeElement(box, vol, split)
		+ getRightChild()->getIntegralForScheffeElement(box, vol, split));
	}
}	
//---------------	


/*Volume of box represented by this multiplied by
real range of this.*/
cxsc::real RealMappedSPnode::getRealAreaRangeWithBox() const
{
	return getRange() * nodeRealVolume();
}

/*Volume of box represented by this multiplied by
real range of this.*/
cxsc::dotprecision 
	RealMappedSPnode::getDotPrecisionAreaRangeWithBox() const
{
	cxsc::dotprecision result(0.0);
	cxsc::accumulate(result, getRange(), nodeRealVolume());
	return result;
}

cxsc::real 
	RealMappedSPnode::getTotalLeafAreaRangeWithBox() const
{
	cxsc::real total = cxsc::Infinity;
	
	if (!hasInfiniteRangeInTree()) {
		total = rnd(getTotalDotPrecisionLeafAreaRangeWithBox());
	}
	
	return total;

}

cxsc::dotprecision 
	RealMappedSPnode::getTotalDotPrecisionLeafAreaRangeWithBox()
																		const
{
	if (isLeaf()) return getDotPrecisionAreaRangeWithBox();
	else {
		return ( getLeftChild()->getTotalDotPrecisionLeafAreaRangeWithBox()
			+ getRightChild()->getTotalDotPrecisionLeafAreaRangeWithBox() );
	}
}

cxsc::real 
	RealMappedSPnode::getTotalAbsLeafAreaRangeWithBox() const
{
	cxsc::real total = cxsc::Infinity;
	
	if (!hasInfiniteRangeInTree()) {
		total = rnd(getTotalDotPrecisionAbsLeafAreaRangeWithBox());
	}
	
	return total;
}

cxsc::dotprecision 
	RealMappedSPnode::getTotalDotPrecisionAbsLeafAreaRangeWithBox()
																		const
{
	if (isLeaf()) return cxsc::abs(getDotPrecisionAreaRangeWithBox());
	else {
		return ( getLeftChild()->getTotalDotPrecisionAbsLeafAreaRangeWithBox()
			+ getRightChild()->getTotalDotPrecisionAbsLeafAreaRangeWithBox() );
	}
}

cxsc::real 
	RealMappedSPnode::getTotalAbsDiffLeafAreaRangeWithBox(
						const RealMappedSPnode& rmsp) const
{
	cxsc::real total = cxsc::Infinity;
	
	if (!hasInfiniteRangeInTree() && !rmsp.hasInfiniteRangeInTree()) {
		RealMappedSPnode diff = *this - rmsp;
		total = diff.getTotalAbsLeafAreaRangeWithBox();
		
	}
	return total;
}

cxsc::dotprecision 
	RealMappedSPnode::
		getTotalDotPrecisionAbsDiffLeafAreaRangeWithBox(
						const RealMappedSPnode& rmsp) const
{
	RealMappedSPnode diff = *this - rmsp;
	return diff.getTotalDotPrecisionAbsLeafAreaRangeWithBox();
}

/* Return a reference to a container of nodes.
 
Contents of container are the leaves descended from this, 
or this if this is a leaf, left to right order. */
RealMappedSPnode::Ptrs& RealMappedSPnode::getLeaves(
			RealMappedSPnode::Ptrs& leaves)
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getLeaves(leaves);
	}

	if (hasRCwithBox()) {
		getRightChild()->getLeaves(leaves);
	}

	if ( isLeaf() ) { // this is a leaf
		leaves.push_back(this);
	}
	return leaves;
}


	
/* Return a reference to a container of const nodes.
 
Contents of container are the leaves descended from this, 
or this if this is a leaf, left to right order. */
RealMappedSPnode::ConstPtrs& RealMappedSPnode::getConstLeaves(
			RealMappedSPnode::ConstPtrs& leaves) const
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getConstLeaves(leaves);
	}

	if (hasRCwithBox()) {
		getRightChild()->getConstLeaves(leaves);
	}

	if ( isLeaf() ) { // this is a leaf
		leaves.push_back(this);
	}
	return leaves;
}

/* Return a reference to a container of nodes.
 
Contents of container are the sub-leaves descended from this, 
or this if this is a sub-leaf, left to right order. 

A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
RealMappedSPnode::Ptrs& RealMappedSPnode::getSubLeaves(
			RealMappedSPnode::Ptrs& subleaves)
{
	if (isSubLeaf()) { // this is a subleaf
		subleaves.push_back(this);
	}
	
	//if children, recurse on the children
	else if (!isLeaf()) {
		getLeftChild()->getSubLeaves(subleaves);
		getRightChild()->getSubLeaves(subleaves);
	}

	return subleaves;
	
}

/* Return a reference to a container of const nodes.
 
Contents of container are the sub-leaves descended from this, 
or this if this is a sub-leaf, left to right order. 

A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
RealMappedSPnode::ConstPtrs& RealMappedSPnode::getConstSubLeaves(
			RealMappedSPnode::ConstPtrs& subleaves) const
{
	if (isSubLeaf()) { // this is a subleaf
		subleaves.push_back(this);
	}
	//if children, recurse on the children
	else if (!isLeaf()) {
		getLeftChild()->getConstSubLeaves(subleaves);
		getRightChild()->getConstSubLeaves(subleaves);
	}

	
	return subleaves;
	
}


void RealMappedSPnode::swapRMSPSR(RealMappedSPnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
	rightChild and nodeName are inherited from base class */
	MappedSPnode < cxsc::real >::swapMSPSR(spn); // use the MSP version
}


// Marginalise
// sort out outDims from reqDims
void RealMappedSPnode::_start_marginalise(
		const std::vector<int>& reqDims)
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"RealMappedSPnode::_start_marginalise(const std::vector<int>&)");
	}
	
	
	if (reqDims.empty()) {
		throw std::invalid_argument(
			"RealMappedSPnode::_start_marginalise(const std::vector<int>&) : reqDims.empty()");
	}
	
	// want to find what dimensions to take out given required dimensions
	ivector box = getBox();
	int dim = VecLen(box);
	int boxLB = Lb(box);
	int boxUB = Ub(box);
	//each of the required dims must be there
	std::vector<int> sorted = reqDims;
	sort(sorted.begin(), sorted.end());
	
	// remove any duplicates
	vector<int>::iterator it = unique (sorted.begin(), sorted.end());
	sorted.resize( it - sorted.begin() );

	if ( (*(sorted.begin()) < 1)) {
		throw std::invalid_argument(
		"RealMappedSPnode::_start_marginalise(const std::vector<int>&) : Dimensions < 1");
	}
	
	if (*(sorted.rbegin()) > boxUB - boxLB + 1)  {
		throw std::invalid_argument(
		"RealMappedSPnode::_start_marginalise(const std::vector<int>&) : Dimension too large for box");
	}
	// could use min and max, but we want the not-req dims anyway
	std::vector<int> outDims;
	
	for (int i = 1; i <= dim; i++) {
		if (!(find(reqDims.begin(), reqDims.end(), i) < reqDims.end())) {
			// dim of box was not in reqDims 
			outDims.push_back(i);
		}
	}
	
	// now we have the dimensions to take out
	_marginalise(outDims); // use internal version
}

// Marginalise, internal version
// marginalise this from the given node down
void RealMappedSPnode::_marginalise(
		const std::vector<int>& outDims)
{
	#ifdef MARG_OUTPUT
		std::cout << "I am " << getNodeName() << std::endl;
	#endif
	
	ivector box = getBox();
	
	// find if this node split on one of the given dimensions
	// returns -1 if no split
	// make sure we do this before fiddling with the children!
	int splitDim = getSplitDim();
	
	// deal with children first
	if (!isLeaf()) {
		
		#ifdef MARG_OUTPUT
			std::cout << "not a leaf - marginalising children\n" << std::endl;
		#endif
	
		getLeftChild()->_marginalise(outDims);		
		getRightChild()->_marginalise(outDims);	
	}								
					
	// now deal with this node itself
	#ifdef MARG_OUTPUT
		if (!isLeaf()) std::cout << "\nback in " << getNodeName() << std::endl;
		std::cout << "\tmy split dimension is " << splitDim << std::endl;
	#endif
	
	std::vector<int>::const_iterator found 
			= find (outDims.begin(), outDims.end(), splitDim);
	
	/* if we split on any of the given dimensions
	 then we drop this from the tree entirely and
	 replace it with the result
	 of adding the two new children together */
	if (found < outDims.end()) { // split on one of the outDims
		// so this will become the result of adding together
		// the two marginlised children
		// note - can never be in here if this is a leaf
		
		#ifdef MARG_OUTPUT
			std::cout << "\tmy split dim is one of the dimensions to be removed" << std::endl;
		#endif
		
		// save who our parent is
		RealMappedSPnode* savedParent = getParent();
		
		// and detach ourselves from parent
		parent = NULL;
		
		// make a copy of the addition left and
		// right children - use the copy constructor to
		// make a temporary, do addition and then just swap 
		
		RealMappedSPnode temp(*getLeftChild());
		
		// keep the same number of elements in the range collection
		temp += (*getRightChild());
		
		swapRMSPSR(temp); //swap me and temp
		
		#ifdef MARG_OUTPUT
			std::cout << "\tre-made me out of parallel addition collation of my children:" 
									<< std::endl;
			std::cout << "\tmy range collection is " << range << std::endl;
			std::cout << "\tand my box is " << getBox() << std::endl;
			if (!isLeaf() ) {
				std::cout << "\tand my children (before renaming) are:" << std::endl;
				getLeftChild()->oneLineOutput(std::cout, 2);
				getRightChild()->oneLineOutput(std::cout, 2);
			}
		#endif
		
		// restore relationship to parent
		parent = savedParent;
		
	}
	else { // did not split on an outdim or is a leaf
		// have to contract this
		// marginalised children (if any) will still be attached
		
		#ifdef MARG_OUTPUT
			std::cout << "\tI did not split on a dimension to take out, or I am a leaf, so need to contract box" << std::endl;
			std::cout << "\tmy box is " << box << std::endl;
		#endif
		
		int dim = VecLen(box);
		int boxLB = Lb(box);
	
		int newDims = dim - outDims.size();
		ivector newBox = ivector(newDims); 
		int index = Lb(newBox);
		int oldindex = boxLB;
	
		// put in the upper and lower bounds for the new box
		// for each dimension that stays	
		for (; oldindex <= Ub(box); oldindex++) {
			std::vector<int>::const_iterator fit 
			= find (outDims.begin(), outDims.end(), (oldindex - boxLB + 1));
			if (!(fit < outDims.end())) { // keep this one
				newBox[index] = box[oldindex];
				index++;
			}
		}
		
		// find the volume we missed
		cxsc::real missingVol = realVolume(box)/realVolume(newBox);
		
		#ifdef MARG_OUTPUT
			std::cout << "\tnew box is " << newBox << std::endl;
			std::cout << "\tmissing volume is " << missingVol << std::endl;
		#endif
		
		//store the child node locations and then temporarily detach them
		RealMappedSPnode* savedLC = getLeftChild();
		RealMappedSPnode* savedRC = getRightChild();
		leftChild = NULL;
		rightChild = NULL;
		
		// also need to store parent
		RealMappedSPnode* savedParent = getParent();
		parent = NULL;
		
		#ifdef MARG_OUTPUT
			std::cout << "\tcurrent range is " << range << std::endl;
		#endif
	
		cxsc::real temp = getRange();
		
		temp *= missingVol;
		
		#ifdef MARG_OUTPUT
			std::cout << "scaled up range is " << temp << std::endl;
		#endif
		
		// replace contents of this with contents of a newly made node	
		RealMappedSPnode 
				tempNode( RealMappedSPnode(newBox, temp) );
		this->swapRMSPSR(tempNode);
		
		#ifdef MARG_OUTPUT
			std::cout << "\tafter contracting, my box is " << getBox() << std::endl;
			std::cout << "\tand after scaling up my range ,";
			std::cout << "\n\tmy range is " << range << std::endl;
		#endif
		
		// put the child pointers back, and reattach to parent
		leftChild = savedLC;
		rightChild = savedRC;
		parent = savedParent;
		
	} // finished else
	
	//if we are the root, recursively rename everything
	if (getParent() == NULL) {
		recursiveRename();
		#ifdef MARG_OUTPUT
			std::cout << "\nNow recursively rename everything from me down\n\n" << std::endl;
		#endif
	}
}


void RealMappedSPnode::_normalise()
{
	if ( isEmpty() ) {
		throw NoBox_Error(
			"RealMappedSPnode::_normalise()");
	}
	
	cxsc::real normaliser = getTotalAbsLeafAreaRangeWithBox();
	if (normaliser == cxsc::Infinity) {
		throw std::runtime_error(
			"RealMappedSPnode::_normalise() : Normalising constant is Infinity");
	}
	if (normaliser <= 0.0) {
		throw std::runtime_error(
			"RealMappedSPnode::_normalise() : Normalising constant is <= 0.0");
	}
	(*this) /= normaliser;
}

// L1 distance between this and another node, with the root counters for normalisers
// no checks on boxes since this should done by another function that calls this one
// no checks null, empty extra - will just explode if there is a problem
cxsc::dotprecision& RealMappedSPnode::_getL1distance(
				cxsc::dotprecision& disL1,
				const RealMappedSPnode * const other) const
{
	#ifdef DEBUG_L1
		std::cout << "\nIn _getL1distance, I am " << getNodeName() << std::endl;
	#endif
	
	// other is not a leaf
	if (!other->isLeaf()) {
		
		#ifdef DEBUG_L1
			std::cout << "Other is not a leaf, other name is " << other->getNodeName() << std::endl;
		#endif
		
		if ( isLeaf() ) {
			
			#ifdef DEBUG_L1
				std::cout << "I am a leaf" << std::endl;
			#endif
		
			// turn it around and use nodeL1Distance with other
			other->nodeL1Distance(disL1, getRange());
		}
		
		else { // I am not a leaf, so recurse
		
			#ifdef DEBUG_L1
				std::cout << "I am not a leaf: recursing" << std::endl;
			#endif
			
			getLeftChild()->_getL1distance(disL1, other->getLeftChild());
			getRightChild()->_getL1distance(disL1, other->getRightChild());
		}
	}
	
	else { // other is a leaf
	
		#ifdef DEBUG_L1
			std::cout << "Other is a leaf, other name is " << other->getNodeName() << std::endl;
		#endif
		
		nodeL1Distance(disL1, other->getRange());
	
	}
	
	return disL1;
}

cxsc::dotprecision& RealMappedSPnode::nodeL1Distance(
								cxsc::dotprecision& disL1,
								cxsc::real other_v) const
{
	#ifdef DEBUG_L1
		std::cout << "\nIn nodeL1distance, I am " << getNodeName() 
		<< ", and disL1 is " << (rnd(disL1)) << endl;
		std::cout << "\nother_v is " << other_v << std::endl;
	#endif
		
	// this is not a leaf
	if (!isLeaf() ) {
		
		#ifdef DEBUG_L1
			std::cout << "I am not a leaf" << std::endl;
		#endif
		
		disL1 = getLeftChild()->nodeL1Distance(disL1,  other_v );
		disL1 = getRightChild()->nodeL1Distance(disL1, other_v );
	}
	else { // this is a leaf
	
		// for me, calculate difference to other_v
		accumulate(disL1, 
				cxsc::abs(cxsc::abs(range) - cxsc::abs(other_v)),
				nodeRealVolume() );
			
		#ifdef DEBUG_L1
			std::cout << "I am a leaf, my range value is " << getRange() << endl;
			std::cout << ", and other value is " << other_v << std::endl;
			std::cout << "my volume is " << nodeRealVolume() << endl;
			std::cout << "adding the following to disL1: ";
			std::cout << (cxsc::abs(cxsc::abs(range) - cxsc::abs(other_v))*nodeRealVolume()) << std::endl;
		#endif
		
	}
	
	#ifdef DEBUG_L1
		std::cout << "disL1 is now: " << (rnd(disL1)) << std::endl;
	#endif
	
	return disL1;
}

// NEW JUNE 2012 for log posteriors
cxsc::dotprecision& RealMappedSPnode::_getLogLikelihood(
								cxsc::dotprecision& loglik,
								int& isnan, int& isposinf, int& isneginf,
								const SPSnode * const spn) const
{
	#ifdef DEBUG_LL
		std::cout << "\nIn _getLogLikelihood, I am " << getNodeName() 
		<< ", and loglik is " << (rnd(loglik)) << endl;
		<< ", and isnan is " << isnan << ", and isinf is " << isinf <<  << endl;
		std::cout << "\nspn is " << spn->getNodeName() << std::endl;
	#endif
	
	// not already nan , and neither this nor spn is a leaf, so we recurse
	if (!isnan && !isLeaf() && !(spn->isLeaf())) {
		
		#ifdef DEBUG_LL
			std::cout << "All okay and neither I nor spn is a leaf" << std::endl;
		#endif
		
		loglik = getLeftChild()->_getLogLikelihood(loglik, 
								isnan, isposinf, isneginf,
								spn->getLeftChild() );
		if (!isnan) { // check again on isnan
			loglik = getRightChild()->_getLogLikelihood(loglik, 
								isnan, isposinf, isneginf,
								spn->getRightChild() );
		}
	}
	else if (!isnan) { // this is a leaf OR spn is a leaf - can already be inf
	
		// add in nj*log(hj) for me
		size_t n = spn->getCounter();
		
		#ifdef DEBUG_LL
			std::cout << "I am a leaf, my range value is " << getRange() << endl;
			std::cout << ", and spn->getCounter is " << n << std::endl;
		#endif
			
		if (n > 0) {
			
			/* if already +infinity we can reset this to nan (- range) or -infinity (0 range)
			* if already -infinity we can reset this to nan (-range) */
		
			real rng = getRange();
			
			if (rng < 0.0) { // nan whenever we find a negative range even if n ==0;
				isnan = 1; // change by ref
				#ifdef DEBUG_LL
					std::cout << "rng < 0.0: adding nothing to log lik, ";
					std::cout << "and isnan = " << isnan << std::endl;
					std::cout << "loglik is now: " << (rnd(loglik)) << std::endl;
				#endif
			}
			
			// -inf will be < 0 and also IsInfinity so the else if gives +infs
			else if (!isneginf) {
				if ( gsl_isinf(_double(rng)) ) {
					if (!isposinf ) {
						// only reset to +ve infinity if not already infinity (+ve or n-ve)
						isposinf = 1; // change by ref
						loglik = cxsc::dotprecision(cxsc::Infinity);
						#ifdef DEBUG_LL
							std::cout << "rng is Infinity: log lik = inf, ";
							std::cout << "and isposinf = " << isposinf << std::endl;
							std::cout << "loglik is now: " << (rnd(loglik)) << std::endl;
						#endif
					}
				}
				else if (rng > 0.0) { // not +ve infinity
					if (!isposinf) {
						accumulate(loglik, (1.0*n), cxsc::ln(getRange()));
						#ifdef DEBUG_LL
							std::cout << "adding the following to log lik: ";
							std::cout << (n * cxsc::ln(rng)) << std::endl;
						#endif
					}
				}
				else { // range must be 0.0 
					isneginf = 1; // change by ref
					loglik = cxsc::dotprecision(-cxsc::Infinity);
					#ifdef DEBUG_LL
						std::cout << "rng is 0.0: log lik = -inf, ";
						std::cout << "and isneginf = " << isneginf << std::endl;
						std::cout << "loglik is now: " << (rnd(loglik)) << std::endl;
					#endif
				}
			}
		}
			
	}
	
	return loglik;
}


std::ostream& RealMappedSPnode::oneLineOutput(std::ostream& os, int level) const
{
	// do me
	for (int i = 0; i < level; ++i) { os << "\t"; }
	os << getNodeName() << "\tRange: " << range;
	os << "\tbox: " << getBox() << std::endl;
	
	// do children
	if (hasLCwithBox()) getLeftChild()->oneLineOutput(os, level+1);
	if (hasRCwithBox()) getRightChild()->oneLineOutput(os, level+1);
	
	return os;
}

// uses comparison operator for nodes
bool subpavings::nodePtrCompare(const subpavings::RealMappedSPnode* lhs,
			const subpavings::RealMappedSPnode* rhs)
{
	return ((*lhs) < (*rhs));
} 

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::RealMappedSPnode & s1, 
			subpavings::RealMappedSPnode & s2) // throw ()
	{
		s1.swapRMSPSR(s2);
	}





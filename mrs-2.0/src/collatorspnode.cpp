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

/*!/ \file         
\brief CollatorSPnode definitions.
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "collatorspnode.hpp"

#include "toolz.hpp"

#include "sptools.hpp"

#include "SmallClasses.hpp"

#include "spnodevisitor.hpp"

#include "subpaving_exception.hpp"

// to use std input/output
#include <iostream>
// format manipulation on streams
#include <iomanip>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <functional>

//#define MARG_OUTPUT // for console output from marginalisation
//#define DEBUG_L1 // comment out to turn off debugging output for L1 distances

#ifdef NDEBUG
	#undef MARG_OUTPUT
	#undef DEBUG_L1
#endif

using namespace subpavings;
using namespace std;




// ------------------------ public member functions ----------------

// Destructor.
CollatorSPnode::~CollatorSPnode()  {}



// no-argument constructor,
CollatorSPnode::CollatorSPnode() {}   // uses the base SPnode class default constructor


// initialised constructor
// initialised with a box
CollatorSPnode::CollatorSPnode(ivector& v)
	: SPnode(v) {}


// Copy constructor
// copies from given node downwards
CollatorSPnode::CollatorSPnode(const CollatorSPnode& other)
	: SPnode(), rangeCollection(other.rangeCollection)
	 
{
	if (other.theBox != NULL) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.nodeName;
	
	
	//recursion on the children
	if (other.leftChild) {
		nodeAddLeft(new CollatorSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.rightChild) {
		nodeAddRight(new CollatorSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;
	

}

// constructor initialised with an SPSnode
// and a normalising constant, eg sum of counts in each node for a histogram
// the summary becomes count /(normalising constant * vol) for the SPSnode
CollatorSPnode::CollatorSPnode(const SPSnode& spn, size_t bigN)
{
	if ( spn.isEmpty() ) {
		throw NoBox_Error(
			"CollatorSPnode::CollatorSPnode(const SPSnode&, size_t)");
	}
	
	/* recursion means that this is repeated if bigN == 0, but 
	 * it would be very unusual for a root with bigN == 0 to have 
	 * children*/
	if ( bigN == 0) {
		bigN = spn.getRootCounter();
	}
	
	theBox = NULL;	
	nodeName = spn.getNodeName();
	
	cxsc::real nodeValue(0.0);
	
	if (bigN == 0 && spn.getCounter() != 0) { // should never happen
		
		throw std::logic_error(
		"CollatorSPnode::CollatorSPnode(const SPSnode&, size_t) : error with counter");
	}
	// allow 0 value for bigN if the spn counter is also 0 -> node value is 0
	if (bigN != 0 && spn.getCounter() != 0) {
		
		nodeValue = (1.0 * spn.getCounter())/
						(1.0 * bigN * spn.nodeRealVolume());
		
	}
	
	rangeCollection = RangeCollectionHist( nodeValue );
	theBox = new ivector(spn.getBox());
	
	//recursion on the children
	if (spn.hasLCwithBox()) {
		nodeAddLeft(new CollatorSPnode(*spn.getLeftChild(), bigN));
	}
	else leftChild=NULL;

	if (spn.hasRCwithBox()) {
		nodeAddRight(new CollatorSPnode(*spn.getRightChild(), bigN));
	}
	else rightChild=NULL;

}

// Copy assignment operator
// copies from given node downwards
CollatorSPnode& CollatorSPnode::operator=(CollatorSPnode rhs)
{
	rhs.swapCollator(*this);
	return *this;
}



// parent and child accessors have to hide the base class implementation
// this is not good but otherwise we get the base class return type
// I've asked around and I can't find a way around it ...

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
CollatorSPnode* CollatorSPnode::getParent() const
{ return (CollatorSPnode*) parent; }


// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node, cast to this type
CollatorSPnode* CollatorSPnode::getLeftChild() const
{ return (CollatorSPnode*) leftChild; }


// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node, cast this type
CollatorSPnode* CollatorSPnode::getRightChild() const
{ return (CollatorSPnode*) rightChild; }


// Get the total of the range collection
cxsc::real CollatorSPnode::getTotalRangeCollection() const
{
	return rangeCollection.getTotal();
}

// Get number of items in the rangeCollection.
size_t CollatorSPnode::getSizeRangeCollection() const
{ return rangeCollection.getSize(); }


bool CollatorSPnode::isEmptyRangeCollection() const
{ return rangeCollection.isEmpty(); }

// Get the total 'area' of the histogram represented by this
// ie sum of abs values x vols for leaf nodes for which this is ancestor
cxsc::real CollatorSPnode::getTotalAbsValueTimesVol() const
{
	// if everything is positive, value x vol for this is equal
	// to sum of leaf node values x vols
	cxsc::real result = getNodeAbsValueTimesVol();
	
	// but if not, and not a leaf, we have to accumulate the absolute values
	if (!isLeaf() && !checkNoNegativeTotalValues() ) {
		// start with an empty dot precision
		// as many copies of empty reals as there are histograms in collation
		cxsc::dotprecision accs(0.0);
		
		// use to accumulate leaf node absolute areas
		accs = getLeafNodeAbsAreaAccumulations(accs);
		
		result = cxsc::rnd( accs );
	}
	return result;
}

//new
std::vector < RealVec >& CollatorSPnode::getAllRangeCollections(
					std::vector < RealVec >& container) const
{
	// fill in me and then recurse, left then right
	container.push_back( rangeCollection.getRangeCollectionValues() );
	
	if (hasLCwithBox()) getLeftChild()->getAllRangeCollections(container);
	if (hasRCwithBox()) getRightChild()->getAllRangeCollections(container);
	
	return container;
}


//new AHABC
// return a reference to a container of reals
// contents being the total normalised values of the leaves descended from this,
// or of this if this is a leaf
// left to right order
RealVec& CollatorSPnode::getLeafNodeFHatValues(RealVec& vals) const
{
	
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getLeafNodeFHatValues(vals);
	}

	if (hasRCwithBox()) {
		getRightChild()->getLeafNodeFHatValues(vals);
	}

	if (isLeaf()) { // this is a leaf
		vals.push_back(getTotalRangeCollection());
	}
	
	if (getParent() == NULL) {
		// this is the parent node
		// normalise the values 
		cxsc::real normalisingConstant = getTotalAbsValueTimesVol();
		if (normalisingConstant <= 0.0) {
			throw std::logic_error(
				"CollatorSPnode::getLeafNodeFHatValues(RealVec&) : Normalising constant is <= 0.0");
		}
		if (normalisingConstant != cxsc::real(1.0)) {
			transform (vals.begin(), vals.end(), vals.begin(),
					bind2nd(divides< cxsc::real>(), normalisingConstant) );
		}
	}
	return vals;
}

RealVec CollatorSPnode::getL1DistancesToAverage() const
{
	RealVec retValue;
	retValue = getL1DistancesToAverage(retValue);
	return retValue;
}

// take a container and return the same container, which has been
// cleared (if necessary) and re-filled with 
// L1-distances-to-average, one for each histogram in collation
RealVec& CollatorSPnode::getL1DistancesToAverage(RealVec& container) const
{
		if ( isEmptyRangeCollection() ) {
			throw UnfulfillableRequest_Error(
				"CollatorSPnode::getL1DistancesToAverage(RealVec&)");
		}
		
		std::size_t n = getSizeRangeCollection();
		
		// start with a container of empty dot precisions
		// as many copies of empty dot prec as there are histograms in collation
		VecDotPrec accs( n, cxsc::dotprecision(0.0) );
		
		// if there is just one collated, distance will be 0
		if (n > 1) {
		
			// fill this with dot precisions
			// each one accumulating (leaf node diff to average) x vol for a particular histogram  
			accs = getLeafNodeAbsDiffToAvAreaAccumulations(accs);
		}
		
		RealVec temp( n );
		
		for (size_t i = 0; i < n; ++i) {
			// put rounded L1 diff into the container
			temp.at(i) = rnd( accs.at(i) ); // round to nearest
		}
		temp.swap(container); 
		return container;
}

RealVec CollatorSPnode::getL1DistancesToAverage(
							const CollatorSPnode * const other) const
{
	RealVec retValue;
	retValue = getL1DistancesToAverage(retValue, other);
	return retValue;
}

// get distances between two pavings
// checks boxes match
// throws exception if this is empty, or if rhs is null or empty
// throws exception if other has nothing collated
// if this has nothing collated, returns an empty RealVec
RealVec& CollatorSPnode::getL1DistancesToAverage(RealVec& container, 
							const CollatorSPnode * const other) const
{
	if ( other == NULL ) {
		throw NullSubpavingPointer_Error(
		"CollatorSPnode::getL1DistancesToAverage(RealVec&, const CollatorSPnode * const)");
	}
	
	if (other->isEmptyRangeCollection()) 
	{
		throw UnfulfillableRequest_Error(
		"CollatorSPnode::getL1DistancesToAverage(RealVec&, const CollatorSPnode * const)");
	}
	
	if ( getBox() != other->getBox() ) {
		throw IncompatibleDimensions_Error(
		"CollatorSPnode::getL1DistancesToAverage(RealVec&, const CollatorSPnode * const)");
	}
	
	

	RealVec temp;
	
	if ( !isEmptyRangeCollection() ) {
		
		std::size_t n = getSizeRangeCollection();
	
		VecDotPrec vecDP( n, cxsc::dotprecision(0.0) );
	
		// take the average of other
		CollatorSPnode av = other->makeAverage();
		
		vecDP = _getL1distances(vecDP, av);
		
		temp.reserve(n);
		
		for (VecDotPrecIt it = vecDP.begin(); it < vecDP.end(); ++it) {
			temp.push_back ( cxsc::rnd(*it) );
		}
	}
	
	temp.swap(container); // clear container and minimize memory
	
	return container;	
}

RealVec CollatorSPnode::getL1Distances(
							const SPSnode * const spn) const
{
	RealVec retValue;
	retValue = getL1Distances(retValue, spn);
	return retValue;
}

RealVec& CollatorSPnode::getL1Distances(RealVec& container, 
							const SPSnode * const spn) const
{
	if ( spn == NULL ) {
		throw NullSubpavingPointer_Error(
		"CollatorSPnode::getL1Distances(RealVec&, const SPSnode * const)");
	}
	
	// turn spn into a collator
	CollatorSPnode temp(*spn);
	
	container = getL1DistancesToAverage(container, &temp);
	
	return container;
}

// return a reference to a container of SPSnodes
// contents being the leaves descended from this, or this if this is a leaf
// left to right order

std::vector< const CollatorSPnode * >& 
		CollatorSPnode::getConstLeaves(
			std::vector< const CollatorSPnode * >& leaves) const
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getConstLeaves(leaves);
	}

	if (hasRCwithBox()) {
		getRightChild()->getConstLeaves(leaves);
	}

	// this is a leaf
	if (!hasLCwithBox() && !hasRCwithBox()) { 
		leaves.push_back(this);
	}
	return leaves;
}

//Output for averages for all the  leaf boxes in this, using tab delimiters
std::ostream& CollatorSPnode::leavesAverageOutputTabs(
						std::ostream &os, int prec) const
{
	// uses  member function leafAverageOutputTabs for node output
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		leafAverageOutputTabs(os, prec);
		os << "\n";
	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesAverageOutputTabs(os, prec);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesAverageOutputTabs(os, prec);
	}
    return os;

}

// Print the details of a specific node in a subpaving
std::ostream& CollatorSPnode::nodePrint(
								std::ostream &os) const
{
	// output for box in form:
	// box, volume, summary data

	if( !isEmpty() ) { // do nothing if there is no box
		
		ivector thisBox = *theBox; // copy theBox

		os << nodeName << "\tBox is:\t";

		for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
			// c-xsc default output for intervals
			os << "\t" << thisBox[i];   }

		os << "\tBox volume is:\t" << nodeRealVolume();
		os << "\trange data:\t" ;

		rangeCollection.outputTabs(os);

	}
	return os;
}

//Output for all the boxes in this
std::ostream& CollatorSPnode::nodesAllPrint(std::ostream &os,
									int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
	nodesAllPrint(os);
	
	os << cxsc::RestoreOpt;

	return os;
}

//Output for all the boxes in this
std::ostream& CollatorSPnode::nodesAllPrint(std::ostream &os) const
{
	if (!isEmpty()) { 
		nodePrint(os);

		//recurse on the children
		if (hasLCwithBox()) {
			getLeftChild()->nodesAllPrint(os);
		}

		if (hasLCwithBox()) {
			getRightChild()->nodesAllPrint(os);
		}
	}

	return os;
}

//Output the range collection
std::ostream& CollatorSPnode::outputNodeRangeCollection(std::ostream &os,
									int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
	rangeCollection.outputTabs(os);

	os << cxsc::RestoreOpt;
	
	return os;
}

//Output the range collection
std::ostream& CollatorSPnode::outputNodeRangeCollection(std::ostream &os) const
{
	rangeCollection.outputTabs(os);
	return os;
}


// recursively allocate a rangeCollection to this and children
// allocation order is this, left child with remainder of allocation, right child with remainder
std::size_t CollatorSPnode::allocateRanges(
			const vector< RealVec >& rangesToAllocate, std::size_t index)
{
	if (index >= rangesToAllocate.size()) {
		throw std::invalid_argument(
		"CollatorSPnode::allocateRanges(	const vector< RealVec >&, size_t) : range allocations too short");
	}

	rangeCollection = RangeCollectionHist(rangesToAllocate[index]);

	size_t newIndex = index+1;
	if (hasLCwithBox()) {
		newIndex = getLeftChild()->allocateRanges(rangesToAllocate, newIndex);
	}
	if (hasRCwithBox()) {
		newIndex = getRightChild()->allocateRanges(rangesToAllocate, newIndex);
	}

	// have done all the children
	// if this is the root, want to check we have used all the allocations
	if (NULL == getParent()) {
		if (newIndex < rangesToAllocate.size()) {
			std::cerr << "Warning - there are more ranges in the ranges to ";
			std::cerr << "allocate than there are nodes in the tree" << std::endl;
		}
	}

	return newIndex;
}

// Find where data should be 
// childInd is an indicator for which child is being checked
// throws exception if node has no box
const CollatorSPnode* CollatorSPnode::findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd) const
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"CollatorSPnode::findContainingNode(const cxsc::rvector&, OPERATIONS_ON)");
	}
	// start at the top
	const CollatorSPnode* retObj = NULL;
	
	if(nodeContains(pt, childInd)) {
		
		if(isLeaf()) {

			// give this node as return value
			retObj = this;

		} // end of isLeaf

		// if not a leaf and contains data
		// recurse on the children if any
		else {

			if(rightChild!=NULL && !rightChild->isEmpty()){
				
				retObj =
				(getRightChild())->findContainingNode(
					pt, ON_RIGHT);
				
			}
			// only try left if we did not find on the right
			if(retObj == NULL && leftChild!=NULL &&
								!leftChild->isEmpty()) {
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
void CollatorSPnode::nodeExpand(int comp)
{
	
	// can only expand if there is a box
	if ( isEmpty() ) {
		throw NoBox_Error("CollatorSPnode::nodeExpand(int)");
	}

	// only do something if this node is a leaf
	if (isLeaf()) {
		
		CollatorSPnode* newLC = NULL;
		CollatorSPnode* newRC = NULL;
		
		try {
			// ivectors to become boxes for new children
			ivector lC, rC;
			// Call Lower() and Upper() to put split boxes
			// into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);
			
			// make and add the new children
			newLC = new CollatorSPnode(lC, rangeCollection);
			newRC = new CollatorSPnode(rC, rangeCollection);
			
			nodeAddLeft(newLC);
			nodeAddRight(newRC);
			// both children get the same rangeCollection as this

			//name the new children
			getLeftChild()->setNodeName(nodeName + "L");
			getRightChild()->setNodeName(nodeName + "R");

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



// as base class
void CollatorSPnode::nodeExpand()
{
	SPnode::nodeExpand();
}


// lazy collation:  just string together the range collections:  Rc = <Ra, Rb>
void CollatorSPnode::addPaving(
				const CollatorSPnode * const rhs)
{
	// if rhs is  NULL or rhs is empty, do nothing
	if (rhs != NULL && !rhs->isEmpty()) {
		
		// this is not empty, boxes don't match
		if (!isEmpty() && getBox() != rhs->getBox() ) {
			throw IncompatibleDimensions_Error(
			"CollatorSPnode::addPaving(const CollatorSPnode * const)");
		}
		
		// gtee that rhs is non null and not empty, but this could be empty
		_lazyCollationNonMinimalUnion(rhs);
	}
}

// lazy collation:  just string together the range collections:  Rc = <Ra, Rb>
void CollatorSPnode::addPaving(
				const SPSnode * const rhs)
{
	// if rhs is  NULL or rhs is empty, do nothing
	if (rhs != NULL && !rhs->isEmpty()) {
	
		// make the sps node into a collator
		CollatorSPnode temp(*rhs);
		
		// add temp to this
		addPaving(&temp);
	}
}


// completes addition/union operation on this tree
// takes a tree with lazy collation and converts to one with
// rangeCollections with one value equal to the sum of the
// present range collection elements
void CollatorSPnode::reduceCollation()
{
	// reduce the range collection on this node
	nodeRangeCollectionReduce();
	// recurse on children
	if (hasLCwithBox())
		getLeftChild()->reduceCollation();
	if (hasRCwithBox())
		getRightChild()->reduceCollation();
}


// addition with lazy addition of collation
CollatorSPnode& CollatorSPnode::operator+= (
									const CollatorSPnode& add)
{
	this->addPaving(&add);
	return *this;
	
}

const CollatorSPnode CollatorSPnode::operator+ (
								const CollatorSPnode& add) const
{
	CollatorSPnode result =(*this);
	
	result+= add;
	
	return result;
}

CollatorSPnode& CollatorSPnode::operator+= (const SPSnode& add)
{
	this->addPaving(&add);
	return *this;
}
			
const CollatorSPnode CollatorSPnode::operator+ (
											const SPSnode& add) const
{
	CollatorSPnode result =(*this);
	
	result+= add;
	
	return result;
}


CollatorSPnode& CollatorSPnode::operator*= (cxsc::real& mult)
{
	selfScalarMult(mult);
	return *this;
}

const CollatorSPnode CollatorSPnode::operator* (cxsc::real& mult) const
{
	CollatorSPnode result = *this;     
	result *= mult;            
	return result;
}


CollatorSPnode& CollatorSPnode::operator/= (cxsc::real& div)
{
		selfScalarDiv(div);
		return *this;
}

const CollatorSPnode CollatorSPnode::operator/ (cxsc::real& div) const
{
	CollatorSPnode result = *this;     
	result /= div;            
	return result;
}


// Marginalise
// marginalise from given node downwards
void CollatorSPnode::marginalise(
		const std::vector<int>& reqDims)
{
	if ( getParent() != NULL ) {
		throw NonRootNode_Error(
			"CollatorSPnode::marginalise(const std::vector<int>&)");
	}
	
	_start_marginalise(reqDims);
}

const CollatorSPnode CollatorSPnode::makeMarginalised(
							const std::vector<int>& reqDims) const
{
	CollatorSPnode result = *this;

	result._start_marginalise(reqDims);
	
	return result;
}

void CollatorSPnode::normalise()
{
	if ( getParent() != NULL ) {
		throw NonRootNode_Error(
			"CollatorSPnode::normalise()");
	}
	
	_normalise();
	
}

const CollatorSPnode CollatorSPnode::makeNormalised() const
{
	
	CollatorSPnode result = *this;

	result._normalise();
	
	return result;
}


// average this
void CollatorSPnode::average()
{
		if ( getParent() != NULL ) {
			throw NonRootNode_Error(
			"CollatorSPnode::average()");
		}
		
		_average();
}

// make a CollatorSPnode which represents an average of the summary
const CollatorSPnode CollatorSPnode::makeAverage() const
{
	CollatorSPnode result = (*this);
	result._average();
	return result;

}


void CollatorSPnode::swapCollator(CollatorSPnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
        rightChild and nodeName are inherited from base class */
		SPnode::swap(spn); // use the base version
		
		std::swap(rangeCollection, spn.rangeCollection);
        
}

std::string CollatorSPnode::nodeStringSummary() const
{
	std::ostringstream oss;
	
	oss << "I am " << getNodeName() << "(address " << this << "),\n";
	oss << "Dimension is " << getDimension() << ", address of box is " << theBox << "\n"; 
	oss << "address of rangeCollection is " << &rangeCollection << ", rangeCollection is ";
	outputNodeRangeCollection(oss) << "\n";
	oss << "my parent is ";
	if (getParent() != NULL) oss << getParent()->getNodeName();
	else oss << "NULL"; 
	oss << ", my left child is ";
	if (getLeftChild() != NULL) oss << getLeftChild()->getNodeName();
	else oss << "NULL"; 
	oss << ", my right child is ";
	if (getRightChild() != NULL) oss << getRightChild()->getNodeName();
	else oss << "NULL";
	
	return oss.str();
	
}


// -------------------------- protected member functions -------------

// initialised constructor
// initialised with a box and a collection for the rangeCollection
CollatorSPnode::CollatorSPnode(ivector& v, 
						const RangeCollectionHist& rangeColl)
	: SPnode(v), rangeCollection(rangeColl)
{}


// Print the details of a single leaf node, using tab delimiters
std::ostream& CollatorSPnode::leafOutputTabs(std::ostream &os) const
{
	if( !isEmpty() ) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox

		// output the nodeName, nodeVolume
		os << nodeName;
		
		os << "\t" << nodeRealVolume() << "\t";

		// followed by the rangeCollection
		rangeCollection.outputTabs(os);

		// followed by intervals making up box using Inf & Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
	}
    return os;
}

// Print the average of the summary of a single leaf node, using tab delimiters
std::ostream& CollatorSPnode::leafAverageOutputTabs(
									std::ostream &os, int prec) const
{
	if( !isEmpty() ) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox

		// output the nodeName, nodeVolume
		os << nodeName;

		os << cxsc::SaveOpt;
		os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
		
		os << "\t" << nodeRealVolume();
		// followed by the average
		os << "\t";
		rangeCollection.outputAverage(os, prec);

		// followed by intervals making up box using Inf & Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
		os << cxsc::RestoreOpt;
	}
    return os;
}

// Accessor for the rangeCollection.
// return a copy of the rangeCollection
RangeCollectionHist  CollatorSPnode::getRangeCollection() const
{   
	return rangeCollection;
}

// Get the average value of the range collection
cxsc::real CollatorSPnode::getAverageRangeCollectionValue() const
{
	return rangeCollection.getAverageValue();
}

// Get the 'area' of the histogram element represented by this
// ie value x vol for this node
cxsc::real CollatorSPnode::getNodeAbsValueTimesVol() const
{
	cxsc::real vol = nodeRealVolume();
	cxsc::real tot = getTotalRangeCollection();
	return cxsc::abs(vol * tot);
}


bool CollatorSPnode::checkNoNegativeTotalValues(bool checkSoFar) const
{
	checkSoFar = checkSoFar && rangeCollection.checkRangeCollectionAllPositive();
	
	if (checkSoFar && hasLCwithBox()) {
		checkSoFar = getLeftChild()->checkNoNegativeTotalValues(checkSoFar);
	}
	if (checkSoFar&& hasRCwithBox()) {
		checkSoFar = getRightChild()->checkNoNegativeTotalValues(checkSoFar);
	}
	return checkSoFar;
}



// Get a collection of containing average of the range collection
RangeCollectionHist CollatorSPnode::getAverageRangeCollection() const
{
	return rangeCollection.makeAverageRangeCollection();
}

// reduce the rangeCollection one value, the accumulation of the present values
void CollatorSPnode::nodeRangeCollectionReduce()
{
		rangeCollection.reduce();
}

// element by element multiplication of this's range collection elements by mult
void CollatorSPnode::nodeScalarMult(cxsc::real& mult)
{
	if (!isEmpty()) {
		rangeCollection*=(mult);
		
	}
}


// scalar multiplier using scalar of type double
// recursive application of nodeScalarMult to multiply whole tree by mult
void CollatorSPnode::selfScalarMult(cxsc::real& mult)
{
	nodeScalarMult(mult);
	// recurse on children
	if (hasLCwithBox())
		getLeftChild()->selfScalarMult(mult);
	if (hasRCwithBox())
		getRightChild()->selfScalarMult(mult);
}


// element by element divison of this's range collection elements by div
void CollatorSPnode::nodeScalarDiv(cxsc::real& div)
{
	if (!isEmpty()) {
		rangeCollection/=(div);
	}
}


// scalar divison 
// recursive application of nodeScalarDiv to divide whole tree by div
void CollatorSPnode::selfScalarDiv(cxsc::real& div)
{
		nodeScalarDiv(div);
		// recurse on children
		if (hasLCwithBox())
			getLeftChild()->selfScalarDiv(div);
		if (hasRCwithBox())
			getRightChild()->selfScalarDiv(div);
}

void CollatorSPnode::_normalise()
{
	if ( isEmpty() ) {
		throw NoBox_Error(
			"CollatorSPnode::_normalise()");
	}
	
	if ( isEmptyRangeCollection() ) {
		throw UnfulfillableRequest_Error(
			"CollatorSPnode::_normalise()");
	}
	
	cxsc::real normalisingConstant = getTotalAbsValueTimesVol();
	if (normalisingConstant <= 0.0) {
		throw std::logic_error("CollatorSPnode::_normalise() : Normalising constant is <= 0.0");
	}
	_normalise(normalisingConstant);
}

void CollatorSPnode::_normalise(const cxsc::real& normalisingConstant)
{
	// if this is a root we use its value x vol as normaliser
	
	// normalise ourselves
	// just need to divide range collection by normalising constant
	rangeCollection.reduce();  // reduce to one value
	rangeCollection /= normalisingConstant;
	
	// recurse on left child
	if ( hasLCwithBox() ) {
		getLeftChild()->_normalise(normalisingConstant);
	}
	
	// recurse on right child
	if ( hasRCwithBox() ) {
		getRightChild()->_normalise(normalisingConstant);
	}
	
}

// Marginalise
// sort out outDims from reqDims
void CollatorSPnode::_start_marginalise(
		const std::vector<int>& reqDims)
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"CollatorSPnode::_start_marginalise(const std::vector<int>&)");
	}
	
	
	if (reqDims.empty()) {
		throw std::invalid_argument(
			"CollatorSPnode::_start_marginalise(const std::vector<int>&) : reqDims.empty()");
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
		"CollatorSPnode::_start_marginalise(const std::vector<int>&) : Dimensions < 1");
	}
	
	if (*(sorted.rbegin()) > boxUB - boxLB + 1)  {
		throw std::invalid_argument(
		"CollatorSPnode::_start_marginalise(const std::vector<int>&) : Dimension too large for box");
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
void CollatorSPnode::_marginalise(
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
		CollatorSPnode* savedParent = getParent();
		parent = NULL;
		
		// make a copy of the lazy collation of left and
		// right children - use the copy constructor to
		// make a temporary, do union and then just swap 
		
		CollatorSPnode temp(*getLeftChild());
		//temp._lazyCollationNonMinimalUnion(getRightChild());
		// change so that we keep the collation all separate
		temp._parallelAddNonMinimalUnion(getRightChild());
		
		swapCollator(temp); //swap me and temp
		
		#ifdef MARG_OUTPUT
			//std::cout << "\tre-made me out of lazy collation of my children:" << std::endl;
			std::cout << "\tre-made me out of parallel addition collation of my children:" 
									<< std::endl;
			std::cout << "\tmy range collection is ";
			outputNodeRangeCollection(std::cout);
			std::cout << std::endl;
			std::cout << "\tand my box is " << getBox() << std::endl;
			if (!isLeaf() ) {
				std::cout << "\tand my children (before renaming) are:" << std::endl;
				getLeftChild()->oneLineOutput(std::cout, 2);
				getRightChild()->oneLineOutput(std::cout, 2);
			}
		#endif
		
		// restore relationship to parent
		parent = savedParent;
		
		// reduce marginal summaries to just one value
		// no - keep all the marginal summaries
		//reduceCollation();
		/*
		#ifdef MARG_OUTPUT
			std::cout << "\tAfter reducing range collection to one value, range collection is ";
			outputNodeRangeCollection(std::cout);
			std::cout << std::endl;
			if (!isLeaf() ) {
				std::cout << "\tand my children (before renaming) are:" << std::endl;
				getLeftChild()->oneLineOutput(std::cout, 2);
				getRightChild()->oneLineOutput(std::cout, 2);
			}
			else { std::cout << std::endl; }
		#endif
		*/
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
		CollatorSPnode* savedLC = getLeftChild();
		CollatorSPnode* savedRC = getRightChild();
		leftChild = NULL;
		rightChild = NULL;
		
		// also need to store parent
		CollatorSPnode* savedParent = getParent();
		parent = NULL;
		
		#ifdef MARG_OUTPUT
			std::cout << "\tcurrent range collection is ";
			outputNodeRangeCollection(std::cout);
			std::cout << std::endl;
		#endif
	
		RangeCollectionHist temp = getRangeCollection();
		//temp.reduce();
		temp*=missingVol;
		
		// replace contents of this with contents of a newly made node	
		CollatorSPnode tempNode( CollatorSPnode(newBox, temp) );
		this->swapCollator(tempNode);
		
		#ifdef MARG_OUTPUT
			std::cout << "\tafter contracting, my box is " << getBox() << std::endl;
			std::cout << "\tand after scaling up my range collection,";
			//std::cout << " and reducing to it one figure,\n\tmy range collection is ";
			std::cout << "\n\tmy range collection is ";
			outputNodeRangeCollection(std::cout);
			std::cout << std::endl;
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


// average this
void CollatorSPnode::_average()
{
	// replace my range collection with average
	RangeCollectionHist tmp = getAverageRangeCollection();
	std::swap(tmp, rangeCollection);
	
	// recurse on children
	if (getLeftChild() != NULL) {
		getLeftChild()->_average();
	}
	if (getRightChild() != NULL) {
		getRightChild()->_average();
	}
}

std::ostream& CollatorSPnode::oneLineOutput(std::ostream& os, int level) const
{
	// do me
	for (int i = 0; i < level; ++i) { os << "\t"; }
	os << getNodeName() << "\tRangeCollection: ";
	outputNodeRangeCollection(os);
	os << "\tbox: " << getBox() << std::endl;
	
	// do children
	if (hasLCwithBox()) getLeftChild()->oneLineOutput(os, level+1);
	if (hasRCwithBox()) getRightChild()->oneLineOutput(os, level+1);
	
	return os;
}


// accumulates absolute value of difference of leaf node areas to average,
// one dot precision in vector for each histogram in collation
// this method accumulates results for over all the leaves recursively
VecDotPrec& CollatorSPnode::getLeafNodeAbsDiffToAvAreaAccumulations(
									VecDotPrec& areaAcc) const
{
	if (getLeftChild() != NULL) {
		areaAcc = getLeftChild()->getLeafNodeAbsDiffToAvAreaAccumulations(areaAcc);
	}
	if (getRightChild() != NULL) {
		areaAcc = getRightChild()->getLeafNodeAbsDiffToAvAreaAccumulations(areaAcc);
	}
	// add on our absolute differences to av x vol if we are a leaf
	if ( isLeaf() ) {

		cxsc::real vol = nodeRealVolume();
		
		RealVec ranges = rangeCollection.getAbsDiffsToAverageRangeCollection();
		
		RealVecItr it = ranges.begin();
		VecDotPrecIt dpit;

		for (dpit = areaAcc.begin(); dpit < areaAcc.end(); dpit++) {

			cxsc::accumulate((*dpit), (*it), vol);

			it++;
		}
	}

	return areaAcc;
}

// accumulates absolute value of leaf node areas,
// one dot precision in vector for each histogram in collation
// this method accumulates results for over all the leaves recursively
cxsc::dotprecision& CollatorSPnode::getLeafNodeAbsAreaAccumulations(
									cxsc::dotprecision& areaAcc) const
{
	if (getLeftChild() != NULL) {
		areaAcc = getLeftChild()->getLeafNodeAbsAreaAccumulations(areaAcc);
	}
	if (getRightChild() != NULL) {
		areaAcc = getRightChild()->getLeafNodeAbsAreaAccumulations(areaAcc);
	}
	// add on our absolute value x vol if we are a leaf
	if ( isLeaf() ) {

		cxsc::real vol = nodeRealVolume();
		
		//RealVec ranges = rangeCollection.getRangeCollectionValues();
		
		cxsc::accumulate( areaAcc, cxsc::abs( getTotalRangeCollection() ), 
									nodeRealVolume() );
	}
		
	return areaAcc;
}


// ------------------ private methods ---------------------------------



// gtee that rhs is non empty and not empty, but this could be empty
void CollatorSPnode::_lazyCollationNonMinimalUnion(
				const CollatorSPnode * const rhs)
{
	// should not need to check on rhs in this condition
	if (isEmpty() && rhs != NULL && !rhs->isEmpty()) {

		*this = CollatorSPnode(*rhs);
		return;
	}
	// if (rhs == NULL || rhs->isEmpty()) do nothing

	// both not null and not empty
	if (!isEmpty() && rhs != NULL && !rhs->isEmpty()) {

		// make copy of rhs to reshape

		CollatorSPnode rhsTemp(*rhs);

		// reshape this and the temporary copy
		_reshapeTreesToUnion(this, &rhsTemp);

		// do the range collections
		_lazyCollationRangeCollections(&rhsTemp);
		return;

	}
}

// gtee that rhs is non empty and not empty, but this could be empty
// range collections combined by addition of elements pair by pair
// ie ith element of this range collection += ith element of rhs' range collection
void CollatorSPnode::_parallelAddNonMinimalUnion(
				const CollatorSPnode * const rhs)
{
	// should not need to check on rhs in this condition
	if (isEmpty() && rhs != NULL && !rhs->isEmpty()) {

		*this = CollatorSPnode(*rhs);
		return;
	}
	// if (rhs == NULL || rhs->isEmpty()) do nothing

	// both not null and not empty
	if (!isEmpty() && rhs != NULL && !rhs->isEmpty()) {

		// make copy of rhs to reshape

		CollatorSPnode rhsTemp(*rhs);

		// reshape this and the temporary copy
		_reshapeTreesToUnion(this, &rhsTemp);

		// do the range collections
		_parallelAddRangeCollections(&rhsTemp);
		return;

	}
}

// reshape two trees so that both have the shape that is the union of the two
void CollatorSPnode::_reshapeTreesToUnion(CollatorSPnode * const lhs,
							CollatorSPnode * const rhs)
{
	if (lhs == NULL || rhs == NULL) {

		throw NullSubpavingPointer_Error(
		"CollatorSPnode::_reshapeTreesToUnion(CollatorSPnode * const, CollatorSPnode * const)");
	}

	// lhs is not a leaf, rhs is a leaf
	if (!lhs->isLeaf() && rhs->isLeaf()) {

		//we need to expand rhs
		rhs->nodeExpand();
	}

	// lhs is a leaf, rhs is not a leaf
	if (lhs->isLeaf() && !rhs->isLeaf()) {

		//we need to expand lhs
		lhs->nodeExpand();
	}

	// we have made sure that if at least one of them was not a leaf,
	// then neither are now leaves
	// now recurse on the children
	if (!lhs->isLeaf() && !rhs->isLeaf()) {
		_reshapeTreesToUnion(lhs->getLeftChild(), rhs->getLeftChild());
		_reshapeTreesToUnion(lhs->getRightChild(), rhs->getRightChild());

	}
}



// assumes that lhs, rhs are the same as this up to the depth of this
// ie are roots of trees of the same shape from these nodes downwards
void CollatorSPnode::_lazyCollationRangeCollections(
				const CollatorSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"CollatorSPnode::_lazyCollationRangeCollections(const CollatorSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_lazyCollationRangeCollections(
												other->getLeftChild());
		getRightChild()->_lazyCollationRangeCollections(
												other->getRightChild());

	}

	// make this range collection have the contents of this and other combined
	rangeCollection += other->rangeCollection;
}

// parallel add range collections from these two nodes down
// assumes that this and other have the same shape from here down
void CollatorSPnode::_parallelAddRangeCollections(
				const CollatorSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"CollatorSPnode::_parallelAddRangeCollections(const CollatorSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_parallelAddRangeCollections(
												other->getLeftChild());
		getRightChild()->_parallelAddRangeCollections(
												other->getRightChild());
	}

	// make this range collection into parallel addition of this and other combined
	rangeCollection.parallelAdd(other->rangeCollection);
}

// L1 distance between this and another node
// no checks on boxes since this should be done by another function that calls this one...
// no checks empty extra - will just explode if there is a problem
// disL1 will have as many dot precisions in it as there are elements in the collator
VecDotPrec& CollatorSPnode::_getL1distances(VecDotPrec& disL1,
										const CollatorSPnode& other) const 
{
	#ifdef DEBUG_L1
		std::cout << "\nIn _getL1distances, I am " << getNodeName() << std::endl;
	#endif
	
	// other is not a leaf
	if (!other.isLeaf()) {
		
		#ifdef DEBUG_L1
			std::cout << "Other is not a leaf, other name is " << other.getNodeName() << std::endl;
		#endif
		
		if ( isLeaf() ) {
			
			#ifdef DEBUG_L1
				std::cout << "I am a leaf" << std::endl;
			#endif
		
			// turn it around and use nodeL1Distance with other
			// have to do this for every element in this, getting back a dot precision each time
			VecDotPrecIt dpit = disL1.begin();
			for (RangeCollectionHist::RangeCollectionHistItr it = rangeCollection.begin();
				it < rangeCollection.end();
				++it) {
					VecDotPrec tmp(1, *dpit);
					tmp = other.nodeL1Distances(tmp, *it);
					*dpit = tmp.back();
					
				}
		}
		
		else { // I am not a leaf, so recurse
		
			#ifdef DEBUG_L1
				std::cout << "I am not a leaf: recursing" << std::endl;
			#endif
			
			getLeftChild()->_getL1distances(disL1, *other.getLeftChild());
			getRightChild()->_getL1distances(disL1, *other.getRightChild());
		}
		
		return disL1;
	}
	
	else { // other is a leaf
	
		#ifdef DEBUG_L1
			std::cout << "Other is a leaf, other name is " << other.getNodeName() << std::endl;
		#endif
		
		cxsc::real other_h = other.getTotalRangeCollection();
		nodeL1Distances(disL1, other_h);
	
	}
	
	return disL1;
}

VecDotPrec& CollatorSPnode::nodeL1Distances(VecDotPrec& disL1,
								cxsc::real other_h) const
{
	#ifdef DEBUG_L1
		std::cout << "\nIn nodeL1distance, I am " << getNodeName() << ", and disL1 is ";
		for (VecDotPrecIt it = disL1.begin(); it < disL1.end(); ++it) {
			cout << cxsc::rnd(*it) << "\t";
		}
		std::cout << "\nother_h is " << other_h << std::endl;
	#endif
		
	// this is not a leaf
	if (!isLeaf() ) {
		
		#ifdef DEBUG_L1
			std::cout << "I am not a leaf" << std::endl;
		#endif
		
		disL1 = getLeftChild()->nodeL1Distances(disL1,  other_h );
		disL1 = getRightChild()->nodeL1Distances(disL1, other_h );
	}
	else { // this is a leaf
	
		// for each hist in me, calculate difference to other_h
		RangeCollectionHist::RangeCollectionHistItr hit = rangeCollection.begin();
		for (VecDotPrecIt dpit= disL1.begin(); 
			dpit < disL1.end(); 
			++ dpit) {
				
			accumulate(*dpit, cxsc::abs(*hit - other_h), nodeRealVolume() );
			++hit;
				
		}
			
		#ifdef DEBUG_L1
			std::cout << "I am a leaf, my values are ";
			for (RangeCollectionHist::RangeCollectionHistItr rit = rangeCollection.begin(); 
					rit < rangeCollection.end();
					++rit) {
					cout << (*rit) << "\t";
			}
			std::cout << ", and other value is " << other_h << std::endl;
			std::cout << "my volume is " << nodeRealVolume() << endl;
			std::cout << "adding the following to disL1: ";
			for (RangeCollectionHist::RangeCollectionHistItr rit = rangeCollection.begin(); 
					rit < rangeCollection.end();
					++rit) {
					cout << ( cxsc::abs(*rit - other_h) * nodeRealVolume() ) << "\t";
			}
			std::cout << std::endl;
		#endif
		
	}
	
	#ifdef DEBUG_L1
		std::cout << "disL1 is now: ";
		for (VecDotPrecIt it = disL1.begin(); it < disL1.end(); ++it) {
			cout << cxsc::rnd(*it) << "\t";
		}
		std::cout << std::endl;
	#endif
	
	return disL1;
}



// ----------------- non member tools functions ----------------------

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::CollatorSPnode & s1, 
			subpavings::CollatorSPnode & s2) // throw ()
	{
		s1.swapCollator(s2);
	}


//Output for all boxes in collator
std::ostream & subpavings::operator<<(std::ostream &os,
					const CollatorSPnode& spn)
{
	spn.nodesAllOutput(os, 1);
	os << std::endl;
	return os;
}

bool subpavings::nodeCompTotalRangeCollection(
						const CollatorSPnode * const lhs,
						const CollatorSPnode * const rhs)
{
	return (lhs->getTotalRangeCollection() 
							< rhs->getTotalRangeCollection());
}

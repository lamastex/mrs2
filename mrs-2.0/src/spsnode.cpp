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
\brief SPSnode (StatsSubPaving) definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "spsnode.hpp"
#include "realmappedspnode.hpp"

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
#include "sptools.hpp"
// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

// to use spsnode splitting classes
#include "splitdecisionobj.hpp"

#include "subpaving_exception.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <stdexcept>

// include fstream so as to be able to output a file
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>

// format manipulation on streams
#include <iomanip>

#include <cassert>

//#define DEBUG_L1 // comment out to turn off debugging output for L1 distances
//#define EMPSDEBUG // debug output for EMPS
//#define DEBUG_NEWMCMC


//#define DEBUG_CHECK_NODE_COUNT
//#define DEBUG_MCMC_SPLIT
//#define DEBUG_MCMC_SPLIT_FAIL

#ifdef NDEBUG
	#undef DEBUG_L1
	#undef EMPSDEBUG
	#undef DEBUG_NEWMCMC
	#undef DEBUG_CHECK_NODE_COUNT
	#undef DEBUG_MCMC_SPLIT
	#undef DEBUG_MCMC_SPLIT_FAIL

#endif

using namespace subpavings;
using namespace std;

// ------------------------ public member functions -----------------

// Default constructor
SPSnode::SPSnode() :  counter(0), splitDim(-1), splitValue(0.0),
							countsOnly(false)
{
	//invokes the base class default constructor
	// then does additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	// reserve space
	spaceIndication = static_cast<size_t>(defaultMaxPts);
	// not sure whether to do this or not - leave for the moment
	dataItrs.reserve(spaceIndication);
	
}

// initialised constructor, initialised with one ivector for the box
// countsOnly will default to false
SPSnode::SPSnode(const ivector& v) : SPnode(v),
	counter(0), splitDim(-1), splitValue(0.0), countsOnly(false)
{
	//invokes the base class constructor with ivector
	// and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	spaceIndication = static_cast<size_t>(defaultMaxPts);
	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication);
	
}

// initialised constructor, initialised with one ivector for the box
// and value for countsOnly
SPSnode::SPSnode(const ivector& v, bool cntOnly) : SPnode(v),
	counter(0), splitDim(-1), splitValue(0.0), countsOnly(cntOnly)
{
	//invokes the base class constructor with ivector
	// and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised
	
	spaceIndication = static_cast<size_t>(defaultMaxPts);
	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication);
	
}

// initialised constructor, initialised with one ivector for the box
// and max for spaceIndication
// countsOnly defaults to false
SPSnode::SPSnode(const ivector& v, size_t max) :
	SPnode(v),
	spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
	countsOnly(false)
{
	//invokes the base class constructor with ivector argument
	// and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication+1);
}


// initialised constructor, initialised with one ivector for the box
// and max for spaceIndication and value for countsOnly
SPSnode::SPSnode(const ivector& v, size_t max, bool cntOnly) :
	SPnode(v),
	spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
	countsOnly(cntOnly)
{
	//invokes the base class constructor with ivector argument
	// and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication+1);
}

// initialised constructor, initialised with a LabBox (labeled box)
// and a max for spaceIndication
// and optionally with cntOnly for countsOnly, defaults to false
SPSnode::SPSnode(const LabBox& lb, size_t max, bool cntOnly) : SPnode(lb),
	spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
	countsOnly(cntOnly)
{
	//invokes the base class constructor with LabBox argument
	//and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication+1);
}

// initialised constructor, initialised with a LabBox (labeled box)
// and optionally with cntOnly for countsOnly, defaults to false
SPSnode::SPSnode(const LabBox& lb, bool cntOnly) : SPnode(lb), counter(0),
	splitDim(-1), splitValue(0.0), countsOnly(cntOnly)
{
	//invokes the base class constructor with LabBox argument
	// and then initialises additional data members

	//dpSums, a vector of dotprecision terms, is not initialised
	//dpSumProducts, similarly not initialised

	spaceIndication = static_cast<size_t>(defaultMaxPts);
	//reserve space - not sure if important - leave for moment
	dataItrs.reserve(spaceIndication);
}

//Copy constructor
// copies from given node downwards
SPSnode::SPSnode(const SPSnode& other) : SPnode(),
	spaceIndication(other.spaceIndication),
	counter(other.counter), dpSums(other.dpSums),
		dpSumProducts(other.dpSumProducts), splitDim(other.splitDim),
	splitValue(other.splitValue), countsOnly(other.countsOnly)
{
	if (other.theBox != NULL) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.nodeName;
	
	//reserve space
	dataItrs.reserve((other.dataItrs).size());
	//copy dataItrs from other to this
	dataItrs = other.dataItrs;
	//recursion on the children
	if (other.leftChild) {
		nodeAddLeft(new SPSnode(*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.rightChild) {
		nodeAddRight(new SPSnode(*(other.getRightChild())));
	}
	else rightChild=NULL;

}


//copy assignment operator
//copies from this node downwards
// copy ellision version, rhs passed by value, no check for self assignment
SPSnode& SPSnode::operator=(SPSnode rhs)
{
	rhs.swapSPS(*this); // make sure we use our version of swap
	return(*this);
}

// Accessor for the counter.
size_t SPSnode::getCounter() const
{ return counter; }

// Accessor for the split dimension.
// Overrides base class
int SPSnode::getSplitDim() const
{ return splitDim; }

// Accessor for the split value.
// Overrides base class
real SPSnode::getSplitValue() const
{ return splitValue; }

// Accessor for the countsOnly value.
real SPSnode::getCountsOnly() const
{ return countsOnly; }

// Accessor for the node's data collection.
// Returns a copy of the node's collection of iterators to the big data set.
NodeData SPSnode::getData() const
{ return dataItrs; }

// Accessor for the parent of a node
// Returns a copy of the pointer to parent node.
SPSnode* SPSnode::getParent() const
{ return (SPSnode*) parent; }

// Accessor for the left child.
// Returns a copy of the pointer to leftChild node cast to this node type
SPSnode* SPSnode::getLeftChild() const
{ return (SPSnode*) leftChild; }

// Accessor for the right child
//Returns a copy of the pointer to rightChild node cast to this node type
SPSnode* SPSnode::getRightChild() const
{ return (SPSnode*) rightChild; }

bool SPSnode::checkTreeStateLegal() const 
{
	return SPnode::checkTreeStateLegal();
}

bool SPSnode::checkTreeStateLegal(size_t minChildPoints, 
								double minVol) const
{
	// check current state is legal by looking at everything not a leaf
	bool legal = true;
	if ( !isLeaf() ) {
		// to be splittable, need at least 2*min vol in this node, so each child has min vol
		legal = isSplittableNode(minChildPoints, minVol);
		if (legal && hasLCwithBox() ) {
				legal = 
					getLeftChild()->checkTreeStateLegal(minChildPoints,
														minVol);
		}
		if (legal && hasRCwithBox() ) {
				legal = 
					getRightChild()->checkTreeStateLegal(minChildPoints,
														minVol);
		}
	}
	
	return legal;
}

bool SPSnode::checkTreeStateLegal(size_t minChildPoints)
{
	// check current state is legal by looking at everything not a leaf
	bool legal = true;
	if ( !isLeaf() ) {
		legal = isSplittableNode(minChildPoints);
		if (legal && hasLCwithBox() ) {
				legal = 
					getLeftChild()->checkTreeStateLegal(minChildPoints);
		}
		if (legal && hasRCwithBox() ) {
				legal = 
					getRightChild()->checkTreeStateLegal(minChildPoints);
		}
	}
	return legal;
}


bool SPSnode::isSplittableNode() const
{
	return SPnode::isSplittableNode();
} 


bool SPSnode::isSplittableNode(size_t minChildPoints, 
								double minVol) const 
{
    bool retValue = true;
	
	if (minVol > 0.0) retValue = (!(nodeVolume() < 2*minVol));
	
	if (retValue) {
		retValue = isSplittableNode(minChildPoints);
	}
	else {
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode: node failed vol test" << endl;
		#endif
	}
	return retValue;
}

bool SPSnode::isSplittableNode(size_t minChildPoints) const
{
    bool retValue = isSplittableNode();  //basic check
	if (!retValue) {
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode: node failed basic is splittable test" << endl;
		#endif
		#ifdef DEBUG_MCMC_SPLIT_FAIL
			cout << "Failed isSplittableNode: I am " << nodeName << endl;
			{
				ivector box = getBox();
				interval maxD = box[MaxDiamComp(box)];
	
				cout << cxsc::SaveOpt;
				cout << Scientific << SetPrecision(35,30);
				cout << "interval to be split is " << maxD << endl;
				cout << cxsc::RestoreOpt;
			}
		#endif
	}
	#ifdef DEBUG_CHECK_NODE_COUNT
		cout << "isSplittableNode minChildPoints = " << minChildPoints << endl;
	#endif
	if (retValue && minChildPoints > 0) {
		retValue = false; // need to retest
	
		size_t  minChildCount = getMinChildCountIfSplitNEW();
		
				
		if ( (counter >= minChildPoints) &&
			((minChildCount == 0) || (minChildCount >= minChildPoints)) ) {
				retValue = true;
			}
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode minChildCount = " << minChildCount << endl;
			cout << "(minChildCount >= minChildPoints) = " << (minChildCount >= minChildPoints) << endl;
	
			cout << "isSplittable = " << retValue << endl;
		#endif
	}
    return retValue;
}

// get the number of datapoints currently associated with this which would
// be associated with the new left child if this node were to split
// remember that the left child's box is open at the split
size_t SPSnode::getLeftCountIfSplit() const
{
	size_t leftCount = 0;
	
	if (isLeaf()) {
		// first find what the dimension for the split would be 
		// if the split were made
		// right hand child's box would be if that child
		// were to be created
		cxsc::ivector box = getBox();
		
		int split = MaxDiamComp(box);
		
		cxsc::real midSplit = cxsc::mid(box[split]);

		// left child would have everything up to but not including
		// midSplit, on the split dimension
		NodeDataItr it;

		for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
			// DataItrs is a container of iterators to a BigDataCollection
			// increment rightCount if the point is in rC
			if(  (**it)[split] < midSplit ) leftCount++;
		}
	}
	else { // already split
		leftCount = getLeftChild()->getCounter();
	
	}

	return leftCount;
}

// The count the right child would have if this node was split.
// Does not split the nodes, just calculates how many of the data points
// currently associated with this node would go to the right child
// if the node were to be split.
size_t SPSnode::getRightCountIfSplit() const
{
	size_t rightCount = 0;
	
	if (isLeaf()) {
		rightCount = counter - getLeftCountIfSplit();
	}
	else {
		rightCount = getRightChild()->getCounter();
	}
    return rightCount;
	
}

// Smallest number of points in either child if this was split.
size_t SPSnode::getMinChildCountIfSplit() const
{
	size_t min = getLeftCountIfSplit();
	if ((counter - min) < min) min = counter - min;
	return min;
	
}

// Smallest number of points in either child if this was split.
size_t SPSnode::getMinChildCountIfSplitNEW() const
{
	size_t min = 0;
	
	if (isLeaf()) {
		min = getLeftCountIfSplit();
	}
	else {
		min = getLeftChild()->getCounter();
	}
	
	if ((counter - min) < min) min = counter - min;
	
	return min;
	
}


// get the number of datapoints currently associated with this which would
// be associated with the left and right children of new children
// if this node were split
// will return a container of the number of points the children
// of each child of target might have, in order
// [0] = left child's left child count, [1] = left child's rght child count,
// [2] = rght child's left child count, [3] = rght child's rght child count,
Size_tVec& SPSnode::getChildrensLeftAndRightCountsIfSplit
					(Size_tVec& grandchildCounts) const
{
	// first find what the children's boxes would be would be
	int splitMe; // variable to hold first longest dimension
	ivector box = getBox();
	MaxDiam(box, splitMe);

	// ivectors to be new boxes for new children
	ivector rCBox;
	ivector lCBox;
	
	// Call Upper() to get what would be the right hand child box
	Upper(box, rCBox, splitMe);
	// Call Lower() to get what would be the left hand child box
	Lower(box, lCBox, splitMe);
	
	// mid point of my box on first longest dimension
	cxsc::real midSplit = cxsc::mid(box[splitMe]);

	// and if those children were split
	// left Child 
	int splitChildren = MaxDiamComp(lCBox);
	
	cxsc::real midSplitLC = cxsc::mid(lCBox[splitChildren]);
	
	// right child 
	// will split on the same dimension as LC
	
	cxsc::real midSplitRC = cxsc::mid(rCBox[splitChildren]);
	
	// now find how many of this node's data points would go right
	// and left children of left and right children
	size_t rightRightCount = 0;
	size_t rightLeftCount = 0;
	size_t leftRightCount = 0;
	size_t leftLeftCount = 0;
	NodeDataItr it;

	for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
		// DataItrs is a container of iterators to a BigDataCollection
		rvector p = **it;
		// increment left child?
		if ( p[splitMe] < midSplit ) {
			if ( p[splitChildren] < midSplitLC ) leftLeftCount++;
			else rightLeftCount++;
		}
		else { // on right of me
			if ( p[splitChildren] < midSplitRC) leftRightCount++;
			else rightRightCount++;
		}
	}

	grandchildCounts.push_back(leftLeftCount);
	grandchildCounts.push_back(rightLeftCount);
	grandchildCounts.push_back(leftRightCount);
	grandchildCounts.push_back(rightRightCount);

	return grandchildCounts;
	
}

// get the number of datapoints currently associated with this which would
// be associated with the left and right children of new children
// if this node were split
// will return a container of the number of points the children
// of each child of target might have, in order
// [0] = left child's left child count, [1] = left child's rght child count,
// [2] = rght child's left child count, [3] = rght child's rght child count,
Size_tVec& SPSnode::getChildrensLeftAndRightCountsIfSplitNEW
					(Size_tVec& grandchildCounts) const
{
	size_t rightRightCount = 0;
	size_t rightLeftCount = 0;
	size_t leftRightCount = 0;
	size_t leftLeftCount = 0;
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn getChildrensLeftAndRightCountsIfSplitNEW, isLeaf = " << isLeaf() << endl;
	
	#endif
	
	if (isLeaf()) {
		// first find what the children's boxes would be would be
		int splitMe; // variable to hold first longest dimension
		ivector box = getBox();
		MaxDiam(box, splitMe);

		// ivectors to be new boxes for new children
		ivector rCBox;
		ivector lCBox;
		
		// Call Upper() to get what would be the right hand child box
		Upper(box, rCBox, splitMe);
		// Call Lower() to get what would be the left hand child box
		Lower(box, lCBox, splitMe);
		
		// mid point of my box on first longest dimension
		cxsc::real midSplit = cxsc::mid(box[splitMe]);

		// and if those children were split
		// left Child 
		int splitChildren = MaxDiamComp(lCBox);
		
		cxsc::real midSplitLC = cxsc::mid(lCBox[splitChildren]);
		
		// right child 
		// will split on the same dimension as LC
		
		cxsc::real midSplitRC = cxsc::mid(rCBox[splitChildren]);
		
		// now find how many of this node's data points would go right
		// and left children of left and right children
		NodeDataItr it;

		for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
			// DataItrs is a container of iterators to a BigDataCollection
			rvector p = **it;
			// increment left child?
			if ( p[splitMe] < midSplit ) {
				if ( p[splitChildren] < midSplitLC ) leftLeftCount++;
				else rightLeftCount++;
			}
			else { // on right of me
				if ( p[splitChildren] < midSplitRC) leftRightCount++;
				else rightRightCount++;
			}
		}
	}
	else { // is not a leaf
	
		// left child
		{
			#ifdef DEBUG_NEWMCMC
				cout << "getLeftChild()->isLeaf = " << getLeftChild()->isLeaf() << endl;
			#endif
			size_t leftCount = getLeftChild()->getCounter();
			if (getLeftChild()->isLeaf()) {
				leftLeftCount = getLeftChild()->getLeftCountIfSplit();
			}
			else { // left child is not a leaf
				leftLeftCount = getLeftChild()->getLeftChild()->getCounter();
			}
			rightLeftCount = leftCount - leftLeftCount;
		}
		// right child
		{
			#ifdef DEBUG_NEWMCMC
				cout << "getRightChild()->isLeaf = " << getRightChild()->isLeaf() << endl;
			#endif
			size_t rightCount = getRightChild()->getCounter();
			if (getRightChild()->isLeaf()) {
				leftRightCount = getRightChild()->getLeftCountIfSplit();
			}
			else { // right child is not a leaf
				leftRightCount = getRightChild()->getLeftChild()->getCounter();
			}
			rightRightCount = rightCount - leftRightCount;
		}
		
	}

	grandchildCounts.push_back(leftLeftCount);
	grandchildCounts.push_back(rightLeftCount);
	grandchildCounts.push_back(leftRightCount);
	grandchildCounts.push_back(rightRightCount);

	return grandchildCounts;
	
}


// fills in container of leaf counts, left to right
Size_tVec& SPSnode::getLeafNodeCounts(Size_tVec& counts) const
{
	if (getLeftChild() != NULL) {
		getLeftChild()->getLeafNodeCounts(counts);
	}
	if (getRightChild() != NULL) {
		getRightChild()->getLeafNodeCounts(counts);
	}
	if (isLeaf()) {

		counts.push_back(counter);
	}
	return counts;
}


// return a reference to a container of SPSnodes
// contents being the leaves descended from this, or this if this is a leaf
// left to right order
SPSnodePtrs& SPSnode::getLeaves(SPSnodePtrs& leaves)
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

// return a reference to a container of const SPnodes
// contents being the leaves descended from this, or this if this is a leaf
// left to right order
SPSnodeConstPtrs& SPSnode::getConstLeaves(SPSnodeConstPtrs& leaves) const
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

// return a reference to a container of SPSnodes
// contents being the sub-leaf children of the given node
// sub-leaf nodes are the parents of leaf nodes and only have leaf nodes
// as children
// left to right order
SPSnodePtrs& SPSnode::getSubLeaves(SPSnodePtrs& subleaves)
{
	if (isSubLeaf()) { // this is a subleaf
		subleaves.push_back(this);
	}
	
	//else if children, recurse on the children
	else if (!isLeaf()) {
		getLeftChild()->getSubLeaves(subleaves);

		getRightChild()->getSubLeaves(subleaves);
	}

	return subleaves;
	
}

// return a reference to a container of SPSnodes
// contents being the sub-leaf children of the given node
// sub-leaf nodes are the parents of leaf nodes and only have leaf nodes
// as children
// left to right order
SPSnodeConstPtrs& SPSnode::getConstSubLeaves(SPSnodeConstPtrs& subleaves) const
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

/* return a reference to a container of SPSnodes
 * contents being the leaves in the intersection of this and spn
 * left to right order */
SPSnodePtrs& SPSnode::getLeavesInIntersection(
		const SPnode * const spn,
		SPSnodePtrs& leaves)
{
	//if both have children, recurse on the children
	if (hasLCwithBox() && spn->hasLCwithBox()) {
		getLeftChild()->getLeavesInIntersection(
						spn->getLeftChild(),
						leaves);
	}

	if (hasRCwithBox()  && spn->hasRCwithBox()) {
		getRightChild()->getLeavesInIntersection(
						spn->getRightChild(),
						leaves);
	}

	if ( isLeaf() || spn->isLeaf() ) { // either is a leaf
		leaves.push_back(this);
	}
	return leaves;
}

/* return a reference to a container of const SPSnodes
 * contents being the leaves in the intersection of this and \a spn,
 * left to right order */
SPSnodeConstPtrs& SPSnode::getConstLeavesInIntersection(
		const SPnode * const spn,
		SPSnodeConstPtrs& leaves) const
{
	//if both have children, recurse on the children
	if (hasLCwithBox() && spn->hasLCwithBox()) {
		getLeftChild()->getConstLeavesInIntersection(
							spn->getLeftChild(),
							leaves);
	}

	if (hasRCwithBox()  && spn->hasRCwithBox()) {
		getRightChild()->getConstLeavesInIntersection(
							spn->getRightChild(),
							leaves);
	}

	if ( isLeaf() || spn->isLeaf() ) { // either is a leaf
		leaves.push_back(this);
	}
	
	// if only one is a leaf
	return leaves;
}

/* return a reference to a container of SPSnodes
 * contents being the sub-leaf children of the intersection 
 * of this and \a spn,
 * left to right order */
SPSnodePtrs& SPSnode::getSubLeavesInIntersection(
									const SPnode * const spn,
									SPSnodePtrs& subleaves)
{
	if ( !isLeaf() && !(spn->isLeaf()) ) {
	
		if ( isSubLeaf() || spn->isSubLeaf() ) { // either is a subleaf
			subleaves.push_back(this);
		}
		
		//else (must be children) recurse on the children
		else {
			getLeftChild()->getSubLeavesInIntersection(
								spn->getLeftChild(),
								subleaves);

			getRightChild()->getSubLeavesInIntersection(
								spn->getRightChild(),
								subleaves);
		}
	}

	return subleaves;
	
}

/* return a reference to a container of SPSnodes
 * contents being the sub-leaf children of the intersection 
 * of this and \a spn,
 * left to right order */
SPSnodeConstPtrs& SPSnode::getConstSubLeavesInIntersection(
									const SPnode * const spn,
									SPSnodeConstPtrs& subleaves) const
{
	if ( !isLeaf() && !(spn->isLeaf()) ) {
	
		if ( isSubLeaf() || spn->isSubLeaf() ) { // either is a subleaf
			subleaves.push_back(this);
		}
		
		//else (must be children) recurse on the children
		else {
			getLeftChild()->getConstSubLeavesInIntersection(
								spn->getLeftChild(),
								subleaves);

			getRightChild()->getConstSubLeavesInIntersection(
								spn->getRightChild(),
								subleaves);
		}
	}

	return subleaves;
	
}


std::pair<size_t, cxsc::real> SPSnode::getNonEmptyBoxSummary() const
{
	if(isEmpty()) 
		throw NoBox_Error("SPSnode::getNonEmptyBoxSummary");
		
	size_t nNonEmptyBoxes = 0;
	real vNonEmptyBoxVolumes(0.0);
	
	accumulateNonEmptyBoxSummary(nNonEmptyBoxes, vNonEmptyBoxVolumes);
	
	return std::pair <size_t, real >(nNonEmptyBoxes, 
								vNonEmptyBoxVolumes/nodeRealVolume());
	
}

//Returns the sum of the count over volume in the leaf nodes
real SPSnode::getSumLeafCountOverVol() const
{
	dotprecision sum(0.0);
	sum = accumulateLeafCountOverVol(sum);
	
	return rnd(sum);
}

//Returns the count in the smallest (by count) leaf node.
size_t SPSnode::getSmallestLeafCount() const
{
	if (isLeaf() ) {
		return getCounter();
	}
	else {
		size_t myL = getLeftChild()->getSmallestLeafCount();
		size_t myR = getRightChild()->getSmallestLeafCount();
		return (myL < myR ? myL : myR);
	}
}

// Returns the count in the largest (by count) leaf node.
size_t SPSnode::getLargestLeafCount() const
{
	if (isLeaf() ) {
		return getCounter();
	}
	else {
		size_t myL = getLeftChild()->getLargestLeafCount();
		size_t myR = getRightChild()->getLargestLeafCount();
		return (myL > myR ? myL : myR);
	}
}


// the 'unnormalised' height of the histogram element represented by this
double SPSnode::getCountOverVolume() const
{
	return (getCounter()/nodeVolume());
	
}


// get the count in the ultimate root node ancestor of this node
size_t SPSnode::getRootCounter() const
{
	size_t retValue = 0;
	if (getParent() == NULL) { // this is root
		retValue = counter;
	}
	else {
		// recurse upwards
		retValue = getParent()->getRootCounter();
	}
	return retValue;
}


// Get the mean of the data covered by the box of a node
// rvector of signaling nans if countsOnly or counter == 0
rvector SPSnode::getMean() const
{
	int dimension = getDimension();
	// set up an rvector retMean of the correct dimensions
	rvector retMean(dimension);
	// loop through the elements in the dpSums vector
	for (int i = 0; i< dimension; i++) {

		// if no data elements each element or if only counts are held,
		// that element of the mean is cxsc::SignalingNaN
		if (countsOnly || (counter == 0)) {
			// cxsc::rvector is indexed 1 to n
			retMean[i+1] = cxsc::SignalingNaN;
		}
		// if data elements, find the element-by-element mean
		else {
			// default cxsc rounding dotprecision rnd_next
			retMean[i+1] = rnd(dpSums[i])/(1.0*counter);
		}
	}// end loop through the elements in dpSums

	return retMean;
}

// Get the variance-covariance vector of the data covered
// by the box of a node
RealVec& SPSnode::getVarCovar(RealVec& varCovar) const
{
	RealVec temp;
	int dimension = getDimension();
	temp.reserve(dimension*dimension);

	// loop through the elements in the dpSumProducts vector
	for (int k = 0; k < dimension*dimension; k++) {

		// counts only held or if 0 or 1 data points
		// each element of the var-covar is a cxsc::SignalingNaN
		if (countsOnly || (counter <= 1)) {
			temp.push_back(cxsc::SignalingNaN);
		}
		// if >1 data points find element-by-element var-covar

		/*the var-covar is the sample var-covar
		which is
		[sumproduct(i,j)-sum(i)sum(j)/counter]/(counter-1)

		element k in the vector of dotprecison sumproducts
		corresponds to row k/n, (row 0 to n-1)
		and column k-row*n (col 0 to n-1)
		in a matrix view of the sumproducts */

		else {
			int i = k/dimension; // row  (int/int = int)
			int j = k - i*dimension; // column

			// make another dotprecision variable
			dotprecision temp1 = dpSumProducts[k];

			dotprecision temp2(0.0);
			// sum(i) x sum(j)
			// default cxsc rounding dotprecision rnd_next
			accumulate(temp2,  rnd(dpSums[i]),
					rnd(dpSums[j]));

			real div = -1.0/counter;

			// sumproduct(i,j) - sum(i)(sum(j)/counter
			// default cxsc rounding
			accumulate(temp1, rnd(temp2), div);
			// calculate the variance covariance element
			temp.push_back(rnd(temp1)/(1.0*(counter-1)));

		}
	}// end loop through the elements in dpSumProducts

	temp.swap(varCovar); // std swap
	return varCovar;
	
}


RealVec SPSnode::getVarCovar() const
{
	RealVec retVarCovar;
	retVarCovar = getVarCovar(retVarCovar);
	return retVarCovar;
}



// Print the details of a of a specific node in a subpaving
std::ostream& SPSnode::nodePrint(std::ostream &os) const
{
	// output for box in form:
	// box, volume, counter, mean, variance covariance, and data
	
	if(theBox != NULL) { // do nothing if there is no box
	
		ivector thisBox = *theBox; // copy theBox

		os << "Box is :";

		for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
			// c-xsc default output for intervals
			os << "  " << thisBox[i];
		}

		os << std::endl;
		os << "Box volume is " << nodeRealVolume() << std::endl;
		os << "Counter = " << counter << std::endl;

		os << "Mean is ";
		nodeMeanPrint(os);
		os << "\nVariance-covariance matrix is" << std::endl;
		nodeVarCovarPrint(os);
		os << "Data is" << std::endl;
		nodeDataPrint(os);
		os << std::endl;

		os << std::endl;
	}

	return os;

}


// Get this node's contribution to loglikelihood
real SPSnode::getLogLik(const size_t n) const
{
	// contribution to loglikelihood is counter*log(counter/(n * vol))

	dotprecision contribution(0.0);
	if ((n > 0) && (counter > 0)) {

		accumulate(contribution, 1.0*counter, log(1.0*counter));
		accumulate(contribution, -1.0*counter, log(1.0*n));
		real nv = nodeRealVolume();
		if (nv < cxsc::MinReal) {
			accumulate(contribution, -1.0*counter, log(nodeVolume()));
		}
		else {
			
			accumulate(contribution, -1.0*counter, ln(nv));
		}
	}
	// default cxsc rounding to nearest
	return rnd(contribution);
	
}

// used for MCMC log posteriors(JUNE 2012 changes)
// Get unscaled log likelihood over tree rooted at this 
cxsc::real SPSnode::getUnscaledTreeLogLik() const
{
	if(isEmpty()) 
		throw NoBox_Error("SPSnode::getUnscaledTreeLogLik()");
		
	dotprecision nlogn(0.0);
	int ndepth = 0.0;
	int depth = 0;
	
	_getUnscaledTreeLogLik(nlogn, ndepth, depth);
	// this accumulates values in nlogn, ndepth
	
	dotprecision ndlog2(0.0);
	accumulate (ndlog2, cxsc::Ln2_real, 1.0*ndepth);
	
	return (cxsc::rnd(nlogn + ndlog2));
	

}


// Get change in scaled contribution to log likelihood split
dotprecision SPSnode::getSplitChangeLogLik() const
{
	
	dotprecision change(0.0);
	
	// if counter is 0 there can be no change on splitting
	if (counter > 0) {

		// first find what the left hand child's counter would be if
		// that child were to be created
		size_t leftCount = getLeftCountIfSplit();

		// current number of data points associated to node is counter
		size_t rightCount = counter-leftCount;

		// current node volume from nodeVolume; each child will have half

		// change is
		//      (lc_count*ln(lc_count) + rc_count*ln(rc_count) +
		//          count*ln(2)) - (count*ln(count)
		// if we split and lc_count, rc_count were the new counts in
		// left and right children respectively
		// note that the terms involving the total count in the histogram
		// and the volume of this node cancel so this change
		// is effectively scaled and does not need to use n

		dotprecision currentEMP(0.0);
		dotprecision childEMP(0.0);

		if (leftCount > 0) accumulate(childEMP, 1.0*leftCount,
						 log(1.0*leftCount));

		if (rightCount > 0) accumulate(childEMP, 1.0*rightCount,
						 log(1.0*rightCount));

		accumulate(childEMP, 1.0*counter, log(2.0));

		accumulate(currentEMP, 1.0*counter, log(1.0*counter));

		change = childEMP - currentEMP;
	}
	return change;
}

// Get change in contribution to log likelihood on merge of leaf children
dotprecision SPSnode::getMergeChangeLogLik() const
{
	if (isLeaf() ) {
		throw subpavings::UnfulfillableRequest_Error(
					"SPSnode::getMergeChangeLogLik()");
	}
	dotprecision change(0.0);
	
	// if counter is 0 there can be no change on merging
	if (counter > 0) {

		// first find what the left hand child's counter is
		size_t leftCount = getLeftChild()->getCounter();

		// and right child
		size_t rightCount = getRightChild()->getCounter();

		// change is (count*ln(count)
		//      - (lc_count*ln(lc_count) + rc_count*ln(rc_count) +
		//          count*ln(2))
		// note that the terms involving the total count in the histogram
		// and the volume of this node cancel so this change
		// is effectively scaled and does not need to use n

		dotprecision currentEMP(0.0);
		dotprecision childEMP(0.0);
		if (leftCount > 0) accumulate(childEMP, 1.0*leftCount,
						 log(1.0*leftCount));
	
		if (rightCount > 0) accumulate(childEMP, 1.0*rightCount,
						 log(1.0*rightCount));
	
		accumulate(childEMP, 1.0*counter, log(2.0));

		accumulate(currentEMP, 1.0*counter, log(1.0*counter));

		change = currentEMP - childEMP;
	}
	return change;
}


// Returns the best (smallest positive or most negative) change in EMP
// from splitting any leaf node under COPERR
// n is total points, used for scaling
dotprecision SPSnode::getBestSplitChangeEMPCOPERR(const size_t n) const
{
	if (isLeaf() ) {
		return getSplitChangeEMPCOPERR(n);
	}
	else {
		dotprecision myChangeL = getLeftChild()->getBestSplitChangeEMPCOPERR(n);
		dotprecision myChangeR = getRightChild()->getBestSplitChangeEMPCOPERR(n);
		return (myChangeL < myChangeR ? myChangeL : myChangeR);
	}
}

// Returns the best (smallest positive or most negative) change in EMP
// from splitting any leaf node under AIC
// term involving n cancels out of change
dotprecision SPSnode::getBestSplitChangeEMPAIC() const
{
	if (isLeaf() ) {
		return getSplitChangeEMPAIC();
	}
	else {
		dotprecision myChangeL = getLeftChild()->getBestSplitChangeEMPAIC();
		dotprecision myChangeR = getRightChild()->getBestSplitChangeEMPAIC();
		return (myChangeL < myChangeR ? myChangeL : myChangeR);
	}
}

// Returns the best (smallest positive or most negative) change in EMP
// from merging any subleaf node under COPERR
// n is total points, used for scaling
dotprecision SPSnode::getBestMergeChangeEMPCOPERR(const size_t n) const
{
	
	if (isLeaf() ) {
		throw subpavings::UnfulfillableRequest_Error(
				"SPSnode::getBestMergeChangeEMPCOPERR(const size_t)");
	}
	if (isSubLeaf() ) {
		return getMergeChangeEMPCOPERR(n);
	}
	else {
		if (getLeftChild()->isLeaf()) {
			return getRightChild()->getBestMergeChangeEMPCOPERR(n);
		}
		else if (getRightChild()->isLeaf()) {
			return getLeftChild()->getBestMergeChangeEMPCOPERR(n);
		}
		else {
			dotprecision myChangeL = getLeftChild()->getBestMergeChangeEMPCOPERR(n);
			dotprecision myChangeR = getRightChild()->getBestMergeChangeEMPCOPERR(n);
			return (myChangeL < myChangeR ? myChangeL : myChangeR);
		}
	}
	
}

// Returns the best (smallest positive or most negative) change in EMP
// from merging any subleaf node under AIC
// term involving n cancels out of change
dotprecision SPSnode::getBestMergeChangeEMPAIC() const
{
	if (isLeaf() ) {
		throw subpavings::UnfulfillableRequest_Error(
					"SPSnode::getBestMergeChangeEMPAIC()");
	}
	
	if (isSubLeaf() ) {
		return getMergeChangeEMPAIC();
	}
	else {
		
		if (getLeftChild()->isLeaf()) {
			return getRightChild()->getBestMergeChangeEMPAIC();
		}
		else if (getRightChild()->isLeaf()) {
			return getLeftChild()->getBestMergeChangeEMPAIC();
		}
		else {
			dotprecision myChangeL = getLeftChild()->getBestMergeChangeEMPAIC();
			dotprecision myChangeR = getRightChild()->getBestMergeChangeEMPAIC();
			return (myChangeL < myChangeR ? myChangeL : myChangeR);
		}
	}

}

// Get this node's scaled contribution to sum term
real SPSnode::getEMPContributionCOPERR(const size_t n) const
{
	// current number of data points associated to node is counter
	// current node volume from nodeVolume, and each child will have half

	#ifdef EMPSDEBUG
		cout << "in getEMPContributionCOPERR, n = " << n << endl;
	#endif

	dotprecision contribution(0.0);
	if ((n > 0) && (counter > 0)) {
		accumulate(contribution, -(1.0*counter)/(1.0*n),
				(1.0*counter)/((n*1.0)*nodeRealVolume()));
	}

	// contribution is -counter^2/(n^2 * vol)
	// default cxsc rounding to nearest

	#ifdef EMPSDEBUG
		cout << "EMPContributionCOPERR = " << rnd(contribution) << endl;
	#endif
	
	//return contribution;
	return rnd(contribution);
}

// Get this node's scaled contribution to EMP under AIC.
real SPSnode::getEMPContributionAIC(const size_t n) const
{
	#ifdef EMPSDEBUG
		cout << "in getEMPContributionAIC, n = " << n << endl;
		cout << "EMPContributionAIC = " << -getLogLik(n) << endl;
	#endif
	return -getLogLik(n);
}

// Get change in scaled contribution to EMP under COPERR on split
dotprecision SPSnode::getSplitChangeEMPCOPERR(const size_t n) const
{
	#ifdef EMPSDEBUG
		cout << "in getSplitChangeEMPCOPERR, n = " << n << endl;
	#endif
	
	// change is 1/(n^2 * vol) x (counter^2 - 2(lc_count^2 + rc_count^2))
	// if we split and lc_count, rc_count were the new counts in
	// left and right children respectively
	// Change is scaled by n, total points in histogram
	dotprecision change;
	change = 0.0;
	
	#ifdef EMPSDEBUG
		cout << "counter = = " << counter << endl;
	#endif
	
	if ((n > 0) && (counter > 0)) {
		// first find what the left hand child's counter would be if that child
		// were to be created
		size_t leftCount = getLeftCountIfSplit();

		#ifdef EMPSDEBUG
			cout << "leftCount = " << leftCount << endl;
		#endif
		// current number of data points associated to node is counter
		// current node volume from nodeVolume, and each child will have half

		real nv = nodeRealVolume();
		#ifdef EMPSDEBUG
			cout << "node volume = " << nv << endl;
			cout << "1.0/((n*1.0)*nv) = " << (10.0/((n*1.0)*nv)) << endl;
			cout << "1.0/((n*1.0)*MinReal) = " << (10.0/((n*1.0)*MinReal)) << endl;
		#endif
			
		if (nv < cxsc::MinReal) {
			change = -cxsc::Infinity;
		}
		else {
			accumulate(change, (1.0*counter)/((n*1.0)*nv),
								(1.0*counter)/(1.0*n));
			#ifdef EMPSDEBUG
				cout << "change stage 1 rnd(change) = " << rnd(change) << endl;
			#endif
			accumulate(change, (1.0*leftCount)/(1.0*n),
								-(2.0*leftCount)/((n*1.0)*nv));
			#ifdef EMPSDEBUG
				cout << "change stage 2 = " << rnd(change) << endl;
			#endif
			accumulate(change, (1.0*(counter - leftCount))/(1.0*n),
								-(2.0*(counter - leftCount))/((1.0*n)*nv));
			#ifdef EMPSDEBUG
				cout << "change stage 3 = " << rnd(change) << endl;
			#endif
		}
	}
	#ifdef EMPSDEBUG
		cout << "returning rnd(change) = " << rnd(change) << endl;
	#endif
	return change;
}

// Get change in sum term in EMP under AIC on split.
dotprecision SPSnode::getSplitChangeEMPAIC() const
{
	return -getSplitChangeLogLik();

}

// Get change in scaled contribution to EMP under COPERR on merge
dotprecision SPSnode::getMergeChangeEMPCOPERR(const size_t n) const
{
	if (isLeaf() ) {
		throw subpavings::UnfulfillableRequest_Error(
					"SPSnode::getMergeChangeEMPCOPERR(const size_t)");
	}
	
	// change is 1/(n^2 * vol) x (2(lc_count^2 + rc_count^2) - counter^2)
	// Change is scaled by n, total points in histogram
	dotprecision change;
	change = 0.0;
	
	if ((n > 0) && (counter > 0)) {
		
		real nv = nodeRealVolume();
					
		if (nv < cxsc::MinReal) {
			change = cxsc::Infinity;
		}
		else {
			
			// first find what the left hand child's counter is
			size_t leftCount = getLeftChild()->getCounter();

			// and right child
			size_t rightCount = counter - leftCount;

			// current number of data points associated to node is counter
			// current node volume from nodeVolume, and each child will have half

			accumulate(change, (1.0*leftCount)/(1.0*n),
								(2.0*leftCount)/(n*nodeVolume()));
			accumulate(change, (1.0*(rightCount))/(1.0*n),
								(2.0*(rightCount))/(n*nodeVolume()));
			accumulate(change, -(1.0*counter)/(n*nodeVolume()),
									(1.0*counter)/(1.0*n));
		}
	}
	return change;
}



// Get change in sum term in EMP under AIC on merge.
dotprecision SPSnode::getMergeChangeEMPAIC() const
{
	if (isLeaf() ) {
		throw subpavings::UnfulfillableRequest_Error(
					"SPSnode::getMergeChangeEMPAIC()");
	}
	return -getMergeChangeLogLik();
}

void SPSnode::reshapeToUnion(const SPnode& other)
{
	SPnode::reshapeToUnion(other);
	
}

void SPSnode::reshapeToUnion(const SPnode& other,
						size_t minChildPoints)
{
	if (isEmpty() || other.isEmpty()) {
		throw NoBox_Error(
		"SPnode::reshapeToUnion(const SPnode&, size_t)");
	}
	if ( getBox() != other.getBox() )  {
		throw IncompatibleDimensions_Error(
		"SPnode::reshapeToUnion(const SPnode&, size_t)");
	}
	if ( !other.checkTreeStateLegal() )
		throw runtime_error(
		"SPnode::reshapeToUnion(const SPnode&, size_t) : other has illegal tree state");
	
	std::string baseErrorFilename("ReshapeErrors");
	std::string errorFilename = getUniqueFilename(baseErrorFilename);
	
	bool success = this->_reshapeToUnion(
					&other, minChildPoints, errorFilename);
	
	// if we returned success there should be no file with that name
	if(!success) {
		std::cerr << "\nCould not exactly reshape this to the union:"
			<< " check " << errorFilename << " for errors\n" << endl;
	}
	
}

//Output for all the  leaf boxes in this, using tab delimiters
// including unscaled EMP contributions and changes if split
std::ostream& SPSnode::leavesOutputTabsWithEMPs(const size_t bigN,
						std::ostream &os, int prec) const
{
	if (parent == NULL) { // root
		std::string headers = "node \t vol \t count \t EMP COPERR ";
		headers += "\t &change \t EMP AIC \t &change \t dimensions \n";
		os << headers;
	}

	// uses  member function leafOutputTabsWithEMPs to generate node output
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		leafOutputTabsWithEMPs(bigN, os, prec);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabsWithEMPs(bigN, os, prec);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabsWithEMPs(bigN, os, prec);
	}
    return os;

}


//Output for all the  leaf boxes in this, using tab delimiters
//including output for the height of histogram bins for a
// normalised histogram based on this tree with this as root
std::ostream& SPSnode::leavesOutputTabsWithHistHeight(
					std::ostream &os, int prec) const
{
	leavesOutputTabsWithHistHeight(counter, os, prec);
	return (os);
}

//Output for all the  leaf boxes in this, using tab delimiters
//including output for the height of histogram bins for a
// normalised histogram based on tree with total number of data points bigN.
std::ostream& SPSnode::leavesOutputTabsWithHistHeight(const size_t bigN,
					std::ostream &os, int prec) const
{
	// uses  member function leafOutputTabs to generate node output
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		leafOutputTabsWithHistHeight(bigN, os, prec);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabsWithHistHeight(bigN, os, prec);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabsWithHistHeight(bigN, os, prec);
	}
    return os;

}

//Output for all the  leaf boxes in this, using tab delimiters
// including unscaled EMP contributions and changes if split
std::ostream& SPSnode::leavesOutputTabsWithHistHeightAndEMPs(
				const size_t bigN, std::ostream &os, int prec) const
{
	if (parent == NULL) { // root
		std::string headers = "node \t vol \t count \t height ";
		headers += "\t EMP COPERR \t &change \t EMP AIC \t &change ";
		headers += "\t dimensions \n";
		os << headers;
	}

	// uses  member function leafOutputTabsWithEMPs to generate node output
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		leafOutputTabsWithHistHeightAndEMPs(bigN, os, prec);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabsWithHistHeightAndEMPs(bigN,
																os, prec);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabsWithHistHeightAndEMPs(bigN,
																os, prec);
	}
    return os;

}

// Get the scaled EMP sum term under COPERR for the tree rooted at this
dotprecision SPSnode::getEMPSumCOPERR(const size_t n) const
{
	dotprecision retValue;
	retValue = 0.0;

	// uses  member function getEMPContributionCOPERR for leaf value
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		retValue = getEMPContributionCOPERR(n);
	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		retValue = retValue + getLeftChild()->getEMPSumCOPERR(n);
	}

	if (getRightChild()!=NULL) {
		retValue = retValue + getRightChild()->getEMPSumCOPERR(n);
	}
	return retValue;

}

// Get the unscaled EMP sum term under AIC for the tree rooted at this
dotprecision SPSnode::getEMPSumAIC(const size_t n) const
{
	dotprecision retValue;
	retValue = 0.0;

	// uses  member function getEMPContributionAIC for leaf result
	if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
		retValue = getEMPContributionAIC(n);
	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		retValue+=getLeftChild()->getEMPSumAIC(n);
	}

	if (getRightChild()!=NULL) {
		retValue+=getRightChild()->getEMPSumAIC(n);
	}
	return retValue;
}


// Expand and split data
void SPSnode::nodeExpand(int comp)
{
	nodeExpansionOnly(comp);    // expand the node
	SplitNever sn;              // dummy split decision maker
	splitData(sn);            // split the data with no further splitting


}

// Expand and split the data with further splitting
void SPSnode::nodeExpand(const SplitDecisionObj& boolTest, int comp)
{
	nodeExpansionOnly(comp);    // expand the node
	// split the data, allowing for further splitting

	splitData(boolTest);
}


// as base class
void SPSnode::nodeExpand()
{
	SPnode::nodeExpand();
}


// Expand and split with further splitting
// finds its own comp argument
void SPSnode::nodeExpand(const SplitDecisionObj& boolTest)
{
	int maxdiamcomp; // variable to hold first longest dimension
	MaxDiam(getBox(), maxdiamcomp);
	nodeExpand(boolTest, maxdiamcomp); // complete nodeExpand
}

// reabsorb both the children of this node
// the stats in this node will be correct so all we have to do is to
// associate the child data with this and delete the children
// This would work even if children are not leaves since gatherData()
// gets data from descendents of the node it operates on
void SPSnode::nodeReabsorbChildren()
{
	// first recursively deal with the children of the children
	if (hasLCwithBox())
		getLeftChild()->nodeReabsorbChildren();
	if (hasRCwithBox())
		getRightChild()->nodeReabsorbChildren();

	if (hasLCwithBox()) {
		getLeftChild()->gatherData(dataItrs);
		delete getLeftChild();
		
	}

	if (hasRCwithBox()) {
		getRightChild()->gatherData(dataItrs);
		delete getRightChild();
		
	}

	// reset splitDim and splitValue to their defaults
	splitDim = -1;
	splitValue = 0.0;

	leftChild = NULL;
	rightChild = NULL;
}

bool SPSnode::splitRootAtLeastToShapeSPS(std::vector < size_t > reqDepths,
								size_t minPoints)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPSnode::splitRootAtLeastToShape(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::splitRootAtLeastToShape(...)");
	}
	
	size_t myDepth = 0;
	
	bool success = _splitAtLeastToShapeSPS(reqDepths,
								myDepth,
								minPoints);
	
	if (success && !reqDepths.empty()) 
		throw std::invalid_argument("SPSnode::splitRootAtLeastToShape(...)");
	
	return success;
	
}

// recursively split children according to instruction
/* reqDepths is our instruction, index is where we should look
 * ensures this has at least shape of instruction, but 
 * this can also be more split at some or all points.
 * Splits from the right*/
bool SPSnode::randomMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						bool saveInstructions)
{
	size_t minPoints = 0;
	
	return randomMCMCSplitRootAtLeastSPS(numLeaves, partitioner,
					rmsp, nlogn, ndepth, minPoints, 
					saveInstructions);
						
}

bool SPSnode::randomMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						bool saveInstructions)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPSnode::randomMCMCSplitRootAtLeastSPS(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::randomMCMCSplitRootAtLeastSPS(...)");
	}
	if (rmsp == NULL) {
		throw NullSubpavingPointer_Error(
		"SPSnode::randomMCMCSplitRootAtLeastSPS(...)");
	}
	
	if(saveInstructions) partitioner.initialiseInstructions(numLeaves);
	
	size_t myDepth = 0;
	
	bool success = _randomMCMCSplitAtLeastSPS(numLeaves,
						partitioner,
						rmsp,
						myDepth,
						nlogn,
						ndepth,
						minPoints);
	
	// clear the instructions if not successful?
	
	#ifdef DEBUG_MCMC_SPLIT
			cout << "randomMCMCSplitRootAtLeastSPS" << endl;
			if (success) cout << "returning sucessful"  << endl;
			else cout << "returning unsucessfully"  << endl;
	#endif
	return success;
	
}

bool SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						bool saveInstructions,
						const std::string& failureLogFilename)
{

	size_t minPoints = 0;
						
	return randomKnuthMCMCSplitRootAtLeastSPS(
						numLeaves,
						partitioner,
						rmsp,
						nlogn, 
						ndepth,
						minPoints,
						saveInstructions,
						failureLogFilename);
}


bool SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						bool saveInstructions,
						const std::string& failureLogFilename)
{
	if ( isEmpty() ) {
		throw NoBox_Error("randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	if (rmsp == NULL) {
		throw NullSubpavingPointer_Error(
		"SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	
	if(saveInstructions) partitioner.initialiseInstructions(numLeaves);
	
	size_t myDepth = 0;
	int p = numLeaves-1;
	int q = numLeaves-1;
	std:: stack < size_t > lastPair;
	
	bool success = _randomKnuthMCMCSplitAtLeastSPS(p, q,
						lastPair,
						partitioner,
						rmsp,
						myDepth,
						nlogn,
						ndepth,
						minPoints,
						failureLogFilename);
	
	// clear the instructions if not successful?
	
	#ifdef DEBUG_MCMC_SPLIT
			cout << "randomKnuthMCMCSplitRootAtLeastSPS" << endl;
			if (success) cout << "returning sucessful"  << endl;
			else cout << "returning unsucessfully"  << endl;
	#endif
	return success;
	
}

#if(0)
bool SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						bool saveInstructions)
{
	if ( isEmpty() ) {
		throw NoBox_Error("randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	if (rmsp == NULL) {
		throw NullSubpavingPointer_Error(
		"SPSnode::randomKnuthMCMCSplitRootAtLeastSPS(...)");
	}
	
	if(saveInstructions) partitioner.initialiseInstructions(numLeaves);
	
	size_t myDepth = 0;
	int p = numLeaves-1;
	int q = numLeaves-1;
		
	bool success = _randomKnuthMCMCSplitAtLeastSPS(p, q,
						partitioner,
						rmsp,
						myDepth,
						nlogn,
						ndepth,
						minPoints);
	
	// clear the instructions if not successful?
	
	#ifdef DEBUG_MCMC_SPLIT
			cout << "randomKnuthMCMCSplitRootAtLeastSPS" << endl;
			if (success) cout << "returning sucessful"  << endl;
			else cout << "returning unsucessfully"  << endl;
	#endif
	return success;
	
}
#endif

// Inserts data into this node
// we are actually inserting an iterator to the data
// childInd is an indicator for which child is being checked
SPSnode* SPSnode::insertOneFind(BigDataItr newItr,
								OPERATIONS_ON childInd,
								const SplitDecisionObj& boolTest)
{
	rvector newData = *newItr;

	// start at the top

	SPSnode* retObj = NULL;

	if(nodeContains(newData, childInd)) {
		
		recalculateStats(newData);

		bool wasLeaf = (isLeaf());
		
		// if it is a leaf, add the data and return this object
		if(wasLeaf) {

			dataItrs.push_back(newItr);

			// give this node as return value
			retObj = this;

			// split if we need to
			if (boolTest(this)) {
				// expand and split data to children

				nodeExpand(boolTest);

			} // end if we need to split
			
		} // end of wasLeaf

		// if not a leaf before we had split, and contains data
		// recurse on the children if any
		if (!wasLeaf) {
			
			if( rightChild != NULL && !rightChild->isEmpty() ){
				
				retObj =
				(getRightChild())->insertOneFind(
					newItr, ON_RIGHT, boolTest);
			}
			// only try left if we did not find on the right
			if(retObj == NULL && leftChild!=NULL &&
								!leftChild->isEmpty()) {
									
				retObj =
				(getLeftChild())->insertOneFind(newItr,
				ON_LEFT, boolTest);
			}
		}

	} // end if node contains

	// will return null if does not contain the data

	return retObj;
}


// Find where data should be 
// childInd is an indicator for which child is being checked
const SPSnode* SPSnode::findContainingNode(const rvector& pt,
								OPERATIONS_ON childInd) const
{
	// start at the top
	if (isEmpty()) {
		throw NoBox_Error(
			"SPSnode::findContainingNode(const rvector&, OPERATIONS_ON)");
	}

	const SPSnode* retObj = NULL;

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

// add two non-minimal pavings in a union operation,
// but with no data attached to it - up to the manager to add data
void SPSnode::unionTreeStructure(const SPSnode * const rhs)
{
	if ( isEmpty() ) {
		throw NoBox_Error(
			"SPSnode::unionTreeStructure(const SPSnode * const)");
	}
	if ( rhs != NULL && !rhs->isEmpty() && (getBox() != rhs->getBox()) )
		{
			throw IncompatibleDimensions_Error(
				"SPSnode::unionTreeStructure(const SPSnode * const)");
		}
	// if rhs is null or empty we carry on
	
	unionNoData(rhs);
	
}

// get distance between two pavings
// checks boxes match
// throws exception if this is empty, or if rhs is null or empty
// if both have no data, distance is 0
// if one has data and the other no data, distance is 1
cxsc::real SPSnode::getL1Distance(const SPSnode * const other) const
{
	cxsc::real retValue(0.0);
	
	if ( other == NULL )
	{
		throw NullSubpavingPointer_Error(
				"SPSnode::getL1Distance(const SPSnode * const)");
	}
	
	if ( isEmpty() || other->isEmpty() ) {
		throw NoBox_Error(
				"SPSnode::getL1Distance(const SPSnode * const)");
	}
	
	if ( getBox() != other->getBox() )
		{
		throw IncompatibleDimensions_Error(
				"SPSnode::getL1Distance(const SPSnode * const)");
		}
	
	std::size_t thisBigN = getCounter();
	std::size_t otherBigN = other->getCounter();
	
	// if both BigN's are zero, L1 distance is 0
	// else if thisBigN is 0 or other bigN is 0 , L1 distance will be 1
	
	if( (thisBigN + otherBigN)) { // at least one of them has some data
	
		if (thisBigN * otherBigN) { // both have some data
		
			cxsc::dotprecision retDP(retValue);
			retDP = _getL1distance(retDP, other, thisBigN, otherBigN);
			retValue = cxsc::rnd(retDP);
		
		}
		else { // only one has some data
		
			retValue += 1.0;
		
		}
	}
	
	return retValue;	
		
}

// clear all the data held in the tree rooted at this, setting all counters and
// sums to zero and emptying data containers
// no change to countsOnly
// can only be called on a root node
void SPSnode::clearAllDataHeld() const
{
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::clearAllDataHeld()");
	}
	stripData(); // clears recursively on this and children
}

// if setTo is true, clear all the optional stats held 
// and set countsOnly to true
// if setTo is false, recalcualte the optional stats held and 
// set countsOnly to false
// can only be called on a root node
void SPSnode::setCountsOnly(bool setTo)
{
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPSnode::setCountsOnly(bool)");
	}
	
	if (setTo) { // change to keeping counts only
		stripOptionalStatsOnly(); // clears recursively on this and children
	}
	else { // change to not keeping counts only, ie to keeping all stats
	
		addOptionalStatsOnly(); // adds recursively on this and children
	}
}


void SPSnode::swapSPS(SPSnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
	rightChild and nodeName are inherited from base class */
	SPnode::swap(spn); // use the base version
	
	std::swap(spaceIndication, spn.spaceIndication);
	std::swap(counter, spn.counter);
	std::swap(dpSums, spn.dpSums);
	std::swap(dpSumProducts, spn.dpSumProducts);
	std::swap(dataItrs, spn.dataItrs);
	std::swap(splitDim, spn.splitDim);
	std::swap(splitValue, spn.splitValue);
	std::swap(countsOnly, spn.countsOnly);
	
}

std::string SPSnode::nodeStringSummary() const
{
	std::ostringstream oss;
	
	oss << "I am " << getNodeName() << "(address " << this << "),\n";
	oss << "Dimension is " << getDimension(); 
	if (isEmpty() ) {
		oss << ", the box is NULL\n"; 
	}
	else {
		oss << ", address of box is " << theBox << "\n"; 
	}
	oss << "spaceIndication is " << spaceIndication << ", splitDim is " << splitDim << ", splitValue is " << splitValue << "\n"; 
	oss << "address of dataItrs is " << &dataItrs << ", countsOnly is " << countsOnly << ", counter is " << counter << "\n";
	oss << "getMean is ";
	nodeMeanPrint(oss);
	oss << "\ngetVarCovar is\n";
	nodeVarCovarPrint(oss);
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

// for debugging
std::string SPSnode::doubleCheckStringSummary() const
{
	std::ostringstream oss;
	
	oss << "I am " << getNodeName() << "(address " << this << "),\n";
	oss << "Dimension is " << getDimension() << ", address of box is " << theBox << "\n"; 
	oss << "spaceIndication is " << spaceIndication << ", splitDim is " << splitDim << ", splitValue is " << splitValue << "\n"; 
	oss << "countsOnly is " << countsOnly << ", counter is " << counter << "\n";
	oss << "address of dataItrs is " << &dataItrs;
	if (isLeaf() ) {
		
		oss << " and dataItrs (giving addresses pointed to by iterators) is " << endl;
	
		for (NodeDataConstItr cit = dataItrs.begin(); cit < dataItrs.end(); ++cit) {
			
			oss << &(**cit) << "\t";
		}
	}
	oss << "\n" << endl;
	
	if ( getLeftChild() ) oss << getLeftChild()->doubleCheckStringSummary() << endl;
	if ( getRightChild() ) oss << getRightChild()->doubleCheckStringSummary() << endl;
	
	return oss.str();
	
}
// ---------------------- protected member functions -------------------

// strip the data from this node and, recursively, its children
// note makes no change to counts only
void SPSnode::stripData() const
{
	//strip me
	counter = 0;
	VecDotPrec().swap (dpSums);
	VecDotPrec().swap (dpSumProducts);
	NodeData().swap (dataItrs);

	// and children
	if (getLeftChild() != NULL) {
		getLeftChild()->stripData();
	}
	if (getRightChild() != NULL) {
		getRightChild()->stripData();
	}
}

// strip the totals held in the sums and sum products from this node 
// and set counts only to true
// and do this recursively on the children
void SPSnode::stripOptionalStatsOnly()
{
	if (!countsOnly) {
		//strip me of the sum collectors only
		VecDotPrec().swap (dpSums);
		VecDotPrec().swap (dpSumProducts);
		countsOnly = true;
	}
	
	// and children
	if (getLeftChild() != NULL) {
		getLeftChild()->stripOptionalStatsOnly();
	}
	if (getRightChild() != NULL) {
		getRightChild()->stripOptionalStatsOnly();
	}
}

// add the totals held in the sums and sum products from this node 
// and set counts only to false
// and do this recursively on the children
void SPSnode::addOptionalStatsOnly()
{
	if (countsOnly) {
		NodeData container;
		recalculateOptionalStatsAndGatherData(container);
	}
}


// split data between two new children
// using an SplitDecisionObj to see if the children should be further split
void SPSnode::splitData(const SplitDecisionObj& boolTest)
{

	// check that both children exist
	if (!hasLCwithBox() || !hasRCwithBox()) {
		throw UnfulfillableRequest_Error(
						"SPSnode::splitData(const SplitDecisionObj&)");
	}

	NodeDataItr dataItr; // iterator

	//divvie the data up amongst the children
	for (dataItr = dataItrs.begin();
		dataItr!= dataItrs.end(); dataItr++) {
		BigDataItr newItr = *dataItr;

		//calls insertOneFind on the children of this node
		// so stats are not recalculated for this node itself
		SPSnode* reinsertedInto = NULL;

		if(rightChild!=NULL && !rightChild->isEmpty()) {

			reinsertedInto =
				(getRightChild())->insertOneFind(
				newItr, ON_RIGHT, boolTest);
		}

		// only try the left if it's not on the right
		if(reinsertedInto==NULL && leftChild!=NULL
		&& !leftChild->isEmpty()) {

			reinsertedInto =
				(getLeftChild())->insertOneFind(
				newItr, ON_LEFT, boolTest);
		}

	}

	clearNodeData();         //clear the data in this node
}



// Only expand the node - no reallocation of data
// add two sibling nodes to this provided that this is a leaf
// comp argument is passed to Upper() and Lower()
// these functions split box in half normal to dimension set by comp
void SPSnode::nodeExpansionOnly(int comp)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPSnode::nodeExpansionOnly(int )");
	}
	// only do something if this SPSnode is a leaf
	if (isLeaf()) {
		
		SPSnode* newLC = NULL;
		SPSnode* newRC = NULL;
		
		try {
			// ivectors to become boxes for new children
			ivector lC, rC;

			// Call Lower() and Upper() to put the split
			// boxes into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);
			
			// when making new children, use constructor
			// that will give space indication (for data)
			// of the size of this node's dataItrs
			size_t space = dataItrs.size();

			newLC = new SPSnode(lC, space, countsOnly);
			newRC = new SPSnode(rC, space, countsOnly);
			
			nodeAddLeft(newLC);
			nodeAddRight(newRC);

			//name the new children
			getLeftChild()->setNodeName(nodeName + "L");
			getRightChild()->setNodeName(nodeName + "R");

			// store the split dimension in this
			splitDim = comp;

			// store the split value in this
			// the split value is the infinum of interval
			// of right child box for dimension split on
			splitValue = _double(Inf(
				((getRightChild())->getBox())[comp]));
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


// recursively split children according to instruction
/* reqDepths is our instruction, index is where we should look
 * ensures this has at least shape of instruction, but 
 * this can also be more split at some or all points.
 * Splits from the right*/
bool SPSnode::_splitAtLeastToShapeSPS(
						std::vector < size_t >& reqDepths,
						size_t myDepth,
						size_t minPoints)
{
	if(reqDepths.empty()) 
		throw std::invalid_argument(
			"SPSnode::_splitAtLeastToShapeSPS(...): : reqDepths.empty()");
		
	bool success = true;
	
	size_t depth = reqDepths.back();
	
	if (myDepth < depth) { // need to try to go down more

		if ( isLeaf() ) {
			
			if (isSplittableNode(minPoints)) {
				// split, using the nodeExpand for this subtype if not base
				nodeExpand();
			}
			else success = false;
		}

		// if we are okay, send instruction down RIGHT SIDE FIRST
		if (success) success = getRightChild()->_splitAtLeastToShapeSPS(reqDepths,
						myDepth+1,
						minPoints);
		if (success) success = getLeftChild()->_splitAtLeastToShapeSPS(reqDepths,
						myDepth+1,
						minPoints);
		
	}
	else if (myDepth == depth) {     // split enough

		/* I am not necessarily a leaf though - allow for mcmc situation
		 * where I want to be split to at least this much */
	
		// knock the last element out
		reqDepths.erase(reqDepths.begin()+reqDepths.size()-1);
	}
	else { // myDepth is > depth
		throw std::invalid_argument(
			"SPSnode::_splitAtLeastToShapeSPS(...): node depth > instruction depth");
	}

	return success;
}


/* random partitioning of this and children until tree rooted at 
this has numleaves or a node cannot be split. rmsp mirrors this,
* and components for calculating likelihood are accumulated in nlogn
* and n depth
* 
* This is my old method for uniform partitioning, before the Knuth
* based version below.*/
bool SPSnode::_randomMCMCSplitAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						size_t myDepth,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints)
{
	#ifdef DEBUG_MCMC_SPLIT
		cout << "_randomMCMCSplitAtLeastSPS, I am " << nodeName << endl;
		cout << "numLeaves = " << numLeaves << " and myDepth = " << myDepth << endl;
	
	#endif
	
	bool success = true;
	if (numLeaves == 1) {
		/* and collect the likelihood from me  */
		// accumulate counter * log counter to nlogn for each child
		if (counter > 0) accumulate(nlogn, 1.0*counter, log(1.0*counter));
		/* accumulate counter*myDepth to ndepth*/
		ndepth += (counter * myDepth);
	}
	else if (numLeaves == 2) {
		if (isSplittableNode(minPoints)) {
			
			/* split, using the nodeExpand for this subtype if not base
			 * we don't need to check it's a leaf for nodeExpand, but 
			 * we only need to replace rmsp if we have split here*/
			if (isLeaf()) {
				nodeExpand();
				// and do the rmsp if it also needs expanding
				rmsp->replaceMe(RealMappedSPnode(*this));
				#ifdef DEBUG_MCMC_SPLIT
					cout << "expanded " << nodeName << endl;
				#endif
			}
			else {
				#ifdef DEBUG_MCMC_SPLIT
					cout << nodeName << " was already expanded" << endl;
				#endif
			}
			
			/* and collect the likelihoods from the immediate children if it's not a leaf*/
			size_t counterLC = getLeftChild()->getCounter();
			size_t counterRC = counter - counterLC;
			
			// accumulate counter * log counter to nlogn for each child
			if (counterLC > 0) accumulate(nlogn, 1.0*counterLC, log(1.0*counterLC));
			if (counterRC > 0) accumulate(nlogn, 1.0*counterRC, log(1.0*counterRC));
			/* accumulate counterLC*(myDepth+1) +  counterRC*(myDepth+1) to ndepth
			 * which is the same as accumulating counter * (myDepth+1)*/
			ndepth += (counter * (myDepth+1));
				
		}
		else {
			success = false;
			
			#if defined (DEBUG_MCMC_SPLIT) || defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "could not expand " << nodeName << endl;
			#endif
			#if defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "numLeaves = " << numLeaves << " and myDepth = " << myDepth << endl;
			#endif
		}
	}
	else {
		
		unsigned long int left = 
					partitioner.generateStatePartition(numLeaves);
		
		#ifdef DEBUG_MCMC_SPLIT
			if (!isLeaf()) cout << "left = " << left << endl;
		#endif
		
		if (isSplittableNode(minPoints)) {
			/* split, using the nodeExpand for this subtype if not base
			 * we don't need to check it's a leaf for nodeExpand, but 
			 * we only need to replace rmsp if we have split here*/
			if (isLeaf()) {
				nodeExpand();
				// and do the rmsp if it also needs expanding
				rmsp->replaceMe(RealMappedSPnode(*this));
				#ifdef DEBUG_MCMC_SPLIT
					cout << "expanded " << nodeName << endl;
				#endif
			}
			else {
				#ifdef DEBUG_MCMC_SPLIT
					cout << nodeName << " was already expanded" << endl;
				#endif
			}
			
			if (!isLeaf()) {
				
				#ifdef DEBUG_MCMC_SPLIT
					cout << "recursing\n" << endl;
				#endif
				
				// left side
				success = getLeftChild()->_randomMCMCSplitAtLeastSPS(
						left,
						partitioner,
						rmsp->getLeftChild(),
						myDepth + 1,
						nlogn,
						ndepth,
						minPoints);
				// right side
				if (success) success = getRightChild()->_randomMCMCSplitAtLeastSPS(
						numLeaves - left,
						partitioner,
						rmsp->getRightChild(),
						myDepth + 1,
						nlogn,
						ndepth,
						minPoints);
			}
			else throw std::logic_error("SPSnode::_randomMCMCSplitAtLeastSPS(...)");
		}
		else {
			success = false;
			#if defined (DEBUG_MCMC_SPLIT) || defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "could not expand " << nodeName << endl;
			#endif
			#if defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "numLeaves = " << numLeaves << " and myDepth = " << myDepth << endl;
			#endif
		}
	}
	
	
	return success;
}


int SPSnode::levelsUpToNextRight() const
{
	#ifdef DEBUG_MCMC_SPLIT
		cout << "\nlevelsUpToNextRight, I am " << nodeName << endl;
	#endif
	assert (getParent() != NULL);
	
	int returnValue = 0;
	
	if (this == getParent()->getRightChild()) {
		
		returnValue = 1 + getParent()->levelsUpToNextRight();
	}
	#ifdef DEBUG_MCMC_SPLIT
		cout << "\nlevelsUpToNextRight, I am " << nodeName;
		cout << " returning levelsUp = " << returnValue << endl;
	#endif
	return returnValue;
}


/* random partitioning of this and children until tree rooted at 
this has numleaves or a node cannot be split. rmsp mirrors this,
* and components for calculating likelihood are accumulated in nlogn
* and n depth*/
bool SPSnode::_randomKnuthMCMCSplitAtLeastSPS(
						unsigned long int p,
						unsigned long int q,
						std::stack< size_t>& lastPair,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						size_t myDepth,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						const std::string& failureLogFilename,
						bool across)
{
	#ifdef DEBUG_MCMC_SPLIT
		cout << "\n_randomKnuthMCMCSplitAtLeastSPS, I am " << nodeName << endl;
		cout << "p = " << p << " q = " << q << " and myDepth = " << myDepth << endl;
	
	#endif
		
	bool success = true;
	
	int dec = 0;
	
	//generate Knuth decision
	dec = partitioner.generateKnuthDecision(p, q, across);
	#ifdef DEBUG_MCMC_SPLIT
		cout << "dec = " << dec << endl;
	#endif
	/* dec = 1 is equivalent to '(', 0 to ')' */
	
	if (dec) { //'('
		assert( p > 0);
		--p;
		/* Expand, then call algorithm on left child */
		
		if (isSplittableNode(minPoints)) {
		
			/* split, using the nodeExpand for this subtype if not base
			 * we don't need to check it's a leaf for nodeExpand, but 
			 * we only need to replace rmsp if we have split here*/
			if (isLeaf()) {
				nodeExpand();
				// and do the rmsp if it also needs expanding
				rmsp->replaceMe(RealMappedSPnode(*this));
				#ifdef DEBUG_MCMC_SPLIT
					cout << "expanded " << nodeName << endl;
				#endif
			}
			else {
				#ifdef DEBUG_MCMC_SPLIT
					cout << nodeName << " was already expanded" << endl;
				#endif
			}
			
			lastPair.push(myDepth+1);
			
			success = getLeftChild()->_randomKnuthMCMCSplitAtLeastSPS(
						p, q,
						lastPair,
						partitioner,
						rmsp->getLeftChild(),
						myDepth+1,
						nlogn, 
						ndepth,
						minPoints,
						failureLogFilename);
		}
		else {
			success = false;
			
			if (!failureLogFilename.empty()) { // log
				ostringstream oss;
				oss << nodeName << "\t" << counter << "\t" << nodeRealVolume() << "\t";
				prettyPrint(oss, getBox());
				outputFile(failureLogFilename, oss.str());
			}
			
			#if defined (DEBUG_MCMC_SPLIT) || defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "could not expand " << nodeName << endl;
			#endif
			#if defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "p = " << p << " q = " << q << " and myDepth = " << myDepth << endl;
			#endif
		}
		
	}
	else { //')'
	
		/* collect the likelihood from me  */
		// accumulate counter * log counter to nlogn for each child
		if (counter > 0) accumulate(nlogn, 1.0*counter, log(1.0*counter));
		/* accumulate counter*myDepth to ndepth*/
		ndepth += (counter * myDepth);
		
		/* only keep choosing if q > 0 */ 
		if ( q > 0) {
		
			--q;
			
			/* move to the last pair  */
			size_t moveToLevel = lastPair.top();
			assert (moveToLevel <= myDepth); 
			lastPair.pop();
			
			int levelsUp = myDepth - moveToLevel;
			
			#ifdef DEBUG_MCMC_SPLIT
				cout << "levelsUp = " << levelsUp << endl;
			#endif
			
			SPSnode* nextRight = getParent();
			#if defined (DEBUG_MCMC_SPLIT)
				cout << "nextRight = " << nextRight->getNodeName() << endl;
			#endif
			RealMappedSPnode* rmspNextRight = rmsp->getParent();
			#if defined (DEBUG_MCMC_SPLIT)
				cout << "rmspNextRight ="  << rmspNextRight->getNodeName() << endl;
			#endif
			
			for (int i = levelsUp; i > 0; --i) {
				nextRight = nextRight->getParent();
				#if defined (DEBUG_MCMC_SPLIT)
					cout << "nextRight = " << nextRight->getNodeName() << endl;
				#endif
				rmspNextRight = rmspNextRight->getParent();
				#if defined (DEBUG_MCMC_SPLIT)
					cout << "rmspNextRight ="  << rmspNextRight->getNodeName() << endl;
				#endif
			}
			nextRight = nextRight->getRightChild();
			#if defined (DEBUG_MCMC_SPLIT)
				cout << "nextRight = " << nextRight->getNodeName() << endl;
			#endif
			rmspNextRight = rmspNextRight->getRightChild();
			#if defined (DEBUG_MCMC_SPLIT)
				cout << "rmspNextRight ="  << rmspNextRight->getNodeName() << endl;
			#endif
			
			#if defined (DEBUG_MCMC_SPLIT)
				cout << "heading right to " << nextRight->getNodeName();
				cout << " (" << rmspNextRight->getNodeName() << ")" << endl;
			#endif
			
			/* flag if we are going across from left sibling to right */	
			across = (levelsUp == 0);
			
			success = nextRight->_randomKnuthMCMCSplitAtLeastSPS(
							p, q,
							lastPair,
							partitioner,
							rmspNextRight,
							myDepth-levelsUp,
							nlogn, 
							ndepth,
							minPoints,
							failureLogFilename,
							across);
		}
	}	
	
	return success;
}


dotprecision& SPSnode::accumulateLeafCountOverVol(dotprecision& sum) const
{
	if (isLeaf()) {  // this is a leaf
		accumulate(sum, 1.0*counter, (1.0/nodeRealVolume()));
	}

	else { // this is not a leaf

		sum = getLeftChild()->accumulateLeafCountOverVol(sum);
		sum = getRightChild()->accumulateLeafCountOverVol(sum); 
	}
	return sum;
	
}

//accumluate number of non-empty boxes and total non-empty box vol
void SPSnode::accumulateNonEmptyBoxSummary(size_t& nNonEmptyBoxes, 
								real& vNonEmptyBoxVolumes) const
{
	if (isLeaf()) {
		if (counter > 0) {
			++nNonEmptyBoxes;
			vNonEmptyBoxVolumes+=nodeRealVolume();
			
		}
	}
	else { // recurse
		getLeftChild()->accumulateNonEmptyBoxSummary(nNonEmptyBoxes, 
												vNonEmptyBoxVolumes);
		getRightChild()->accumulateNonEmptyBoxSummary(nNonEmptyBoxes, 
												vNonEmptyBoxVolumes);
	}
	// we want to count all empty leaves even if parent is empty too
}


// Print the data in a node if any
std::ostream& SPSnode::nodeDataPrint(std::ostream &os) const
{
	if (!dataItrs.empty()) {

		NodeDataItr dataItr;
		
		int dimension = getDimension();
		
		for (dataItr = dataItrs.begin();
			dataItr!= dataItrs.end(); dataItr++) {

			BigDataItr bigIt = *dataItr;
			rvector theData = *bigIt;

			for (int i = 1; i < dimension + 1; i++) {
				os << theData[i] << "  "; // print data
			}   // end loop through data elements

		} // end loop through data container
	} // end if counter > 0
	// if no data, ie counter = 0, then just return os unaltered

	return os;
}

// Print the mean of the data in a node
std::ostream& SPSnode::nodeMeanPrint(std::ostream &os) const
{

	if ((counter > 0) && !countsOnly) {
		int dimension = getDimension();
		// loop through the elements in the dpSums vector
		for (int i = 0; i< dimension; i++) {
			// default cxsc rounding of dotprecision
			// to rnd_next
			os << "  " << (rnd(dpSums[i])/(1.0*counter));

		}// end loop through the elements in dpSums
		
	} // end if
	// if no data, ie counter = 0, or if we are only keeping counts
	// then just return os unaltered

	return os;

}

// Print the variance covariance matrix of the data in a node
std::ostream& SPSnode::nodeVarCovarPrint(std::ostream &os) const
{
	int dimension = getDimension();
		
	if ((counter > 0) && !countsOnly) {

		RealVec varCovar;
		varCovar = getVarCovar(varCovar);

		/* element k in the vector representing the
		variance-covariance matrix corresponds to
		row k/n, (row 0 to n-1) and column k-row*n (col 0 to n-1)
		in a matrix view variance-covariance */
		
		// loop through the elements and print as matrix
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				os << "  " << varCovar[(i*dimension)+j];
			}
			os << std::endl;
		}
	}
	else { // just put in new lines
		for (int i = 0; i < dimension; i++) {
			os << std::endl;
		}
	}
	return os;

}

// Print the details of a single leaf node, using tab delimiters
std::ostream& SPSnode::leafOutputTabs(std::ostream &os) const
{
	if(theBox != NULL) { // do nothing if there is no box
	
		ivector thisBox = *theBox; // copy of theBox
	
		// output the node name, nodeVolume, counter
		os << nodeName;
		
		os << "\t" << nodeRealVolume();
		os << "\t" << counter;
		// followed by the intervals of box using Inf and Sup
		// ie unlike cxsc output, there is no [  ] around them

		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}

	}
    return os;
}



// Print the details of a single leaf node, using tab delimiters
// including EMP contributions and changes if split
std::ostream& SPSnode::leafOutputTabsWithEMPs(const size_t bigN,
						std::ostream &os, int prec) const
{
	if(theBox != NULL) { // do nothing if there is no box

		
		ivector thisBox = *theBox; // copy of theBox

		// output the name, nodeVolume, counter
		os << nodeName;
		
		os << cxsc::SaveOpt;
		os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
		
		os << "\t" << nodeRealVolume();
		os << "\t" << counter;
		// EMP contributions and changes if split
		os << "\t" << getEMPContributionCOPERR(bigN);
		os << "\t" << rnd(getSplitChangeEMPCOPERR(bigN));
		os << "\t" << getEMPContributionAIC(bigN);
		os << "\t" << rnd(getSplitChangeEMPAIC());
		
		// followed by the intervals of box using Inf and Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
		os << cxsc::RestoreOpt;

	}
    return os;
}

// Print the details of a single leaf node, using tab delimiters
// includes the height = n/(N*vol) where n is count in this leaf node,
// N is count over whole histogram, vol is volume of this leaf node
std::ostream& SPSnode::leafOutputTabsWithHistHeight(const size_t bigN,
						std::ostream &os, int prec) const
{
	if(theBox != NULL) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy of theBox
		
		// output the node name, nodeVolume, counter, counter/(bigN * vol)
		os << cxsc::SaveOpt;
		os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
		
		cxsc::real vol = nodeRealVolume();
		
		os << setprecision(prec);
		os << "\t" << vol;
		os << "\t" << counter;
		os << "\t" << (1.0 * counter)/(vol * (1.0* bigN));
		// followed by the intervals of box using Inf and Sup
		// ie unlike cxsc output, there is no [  ] around them
		
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; ++i) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
		os << cxsc::RestoreOpt;

	}
    return os;
}

// Print the details of a single leaf node, using tab delimiters
// including EMP contributions and changes if split
std::ostream& SPSnode::leafOutputTabsWithHistHeightAndEMPs(
								const size_t bigN, std::ostream &os,
								int prec) const
{
	if(theBox != NULL) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy of theBox
		
		// output the node name, nodeVolume, counter, counter/(n * vol)
		os << nodeName;
		
		os << cxsc::SaveOpt;
		os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
		
		cxsc::real vol = nodeRealVolume();
		
		os << "\t" << vol;
		os << "\t" << counter;
		os << "\t" << (1.0 * counter)/(vol * (1.0*bigN));
		// EMP contributions and changes if split
		os << "\t" << getEMPContributionCOPERR(bigN);
		os << "\t" << rnd(getSplitChangeEMPCOPERR(bigN));
		os << "\t" << getEMPContributionAIC(bigN);
		os << "\t" << rnd(getSplitChangeEMPAIC());

		// followed by the intervals of box using Inf and Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
		os << cxsc::RestoreOpt;

	}
    return os;
}

// Clears the node's data collection.
void SPSnode::clearNodeData() const
{ NodeData().swap (dataItrs); }

// recalculate the counter and accumulated sum
// and accumulated sumproducts
void SPSnode::recalculateStats(const rvector& newdata) const
{
	counter++;  // update the counter

	if (!countsOnly) {

		recalculateSums(newdata); // update the sums

		recalculateSumProducts(newdata); // update the sumproducts
	}
}

// recalculate the accumulated sum
void SPSnode::recalculateSums(const rvector& newdata) const
{
	int dimension = getDimension();
		
	if (dpSums.empty()) {   //nothing in the sums yet
		// reserve space in dpSums for all elements of the mean
		
		dpSums.reserve(dimension);

		// for each dimnsn of data, initialise element
		for (int i = 0; i< dimension; i++) {
			dotprecision dp;
			dp = 0.0;
			dpSums.push_back(dp);
		}
		
	}

	// make a dot precision variable out of the ith element
	// of the rvector of new data and store in dpSums
	for (int i = 1; i< dimension + 1; i++) {
		// rvectors indexed 1 to n, vectors indexed 0 to n-1
		
		accumulate(dpSums[i-1], newdata[i], 1.0);
	}

}


// recalculate the accumulated sumproducts
void SPSnode::recalculateSumProducts(const rvector& newdata) const
{
	/* the sumproducts can be thought of as an nxn matrix,
	which is implemented here as a nxn element vector of
	dotprecision variables, using row-major order.
	Ie the m-th element (m = 0, . . . nxn-1) is in row floor(m/n)
	and column m-rowxn in the matrix configuration.
	Or, the sumproduct of elements i and j in an rvector,
	i,j = 0,...,n-1, is element m=(ixn+j) of the sumproducts
	vector. */
	int dimension = getDimension();
		
	if (dpSumProducts.empty()) {    //nothing there yet
		// reserve space for all elements
		dpSumProducts.reserve(dimension*dimension);

		// for each dimnsn^2 of data, initialise element
		for (int i = 0; i< (dimension*dimension); i++) {
			dotprecision dp;
			dp = 0.0;
			dpSumProducts.push_back(dp);
		}
	}

	// make a dot precision variable out of the ith element
	// and jth element of the of the rvector of new data and
	// store in dpSumProducts.
	for (int i = 1; i < dimension + 1; i++) {
		// only need to do columns 1 to i because of symmetry
		for (int j = 1; j< i + 1; j++) {

			int index = (i-1)*dimension + (j-1);
			// rvectors indexed 1 to n
			accumulate(dpSumProducts[index],
					newdata[i], newdata[j]);

			//if not on the diagonal of the matrix,
			// we can also fill in the symmetric element
			if (i!=j) {
				int sym_index = (j-1)*dimension
					+ (i-1);
				dpSumProducts[sym_index] =
					dpSumProducts[index];
			} // end if
		}// end j-loop
	}// end i-loop

	// sumproducts has been updated for new datapoint
}


// gather up all the data in children of a node
void SPSnode::recalculateOptionalStatsAndGatherData(NodeData& container)
{
	countsOnly = false;
	
	// make sure the current containers are empty
	VecDotPrec().swap (dpSums);
	VecDotPrec().swap (dpSumProducts);
	
	NodeData myContainer;

	if (!isLeaf()) {
		
		if (hasLCwithBox()) {
			getLeftChild()->recalculateOptionalStatsAndGatherData(myContainer);
		}
		if (hasRCwithBox()) {
			getRightChild()->recalculateOptionalStatsAndGatherData(myContainer);
		
		}
	}
	else { // is a leaf
		// just my data
		
		myContainer = dataItrs;
		
	}
	
	// now do my stats
	assert(myContainer.size() == getCounter() );
	recalculateOptionalStatsForData(myContainer);
	
	// add my data to the container
	container.insert(container.end(),
					myContainer.begin(),
					myContainer.end());
	
}


void SPSnode::recalculateOptionalStatsForData(const NodeData& nodedata) const
{
	if (!countsOnly) {
		for (NodeDataConstItr it = nodedata.begin(); 
										it < nodedata.end(); it++) {
			// nodedata is a container of iterators to a BigDataCollection
			rvector rv = *(*it);
			
			recalculateSums(rv);
			recalculateSumProducts(rv);
		}
	}
}


// ---------------------- private member functions -------------------


// gather up all the data in children of a node
NodeData& SPSnode::gatherData(NodeData& container) const
{
	if (!isLeaf()) {
		if (hasLCwithBox()) {
			container =
				getLeftChild()->gatherData(container);
		}
		if (hasRCwithBox()) {
			container =
				getRightChild()->gatherData(container);
		}
	}
	else { // is a leaf
		// copy data from dataItrs into temp container
		container.insert(container.end(),
						dataItrs.begin(),
						dataItrs.end());
	}

	return container;
}


// this is guaranteed to be non-empty, rhs can be empty or null
// if rhs is not null or empty, it must have the same box as this
 void SPSnode:: unionNoData(const SPSnode * const rhs)
 {
	// strip the data off me
	stripData();
	
	// both not null and not empty
	if (rhs != NULL && !rhs->isEmpty()) {

		// reshape using rhs as the other node
		_reshapeToUnion(rhs);
		
	}
}

void SPSnode::_reshapeToUnion(const SPnode * const other)
{
	SPnode::_reshapeToUnion(other);
}

/* reshape this tree to have union of this shape and shape of other
 no checks on boxes since this should be redundant if used by unionNoData...
*/
bool SPSnode::_reshapeToUnion(const SPnode * const other,
							size_t minChildPoints,
							const std::string& errorFilename)
{
	// indictator for being able to do union exactly
	bool success = true;
	
	if ( other != NULL && !(other->isEmpty()) ) {

		// this is not a leaf, other is a leaf
		if (!isLeaf() && other->isLeaf()) {

			// no need to do anything
		}

		// this is a leaf, other is not a leaf
		if (isLeaf() && !other->isLeaf()) {

			//we need to expand this
			if (isSplittableNode(minChildPoints)) nodeExpand();
			else {
				success = false;
				// log file
				std::string line = "Could not split " + getNodeName() 
				 + " because of minChildPoints";
				
				outputFile(errorFilename, line); 
				
			}
			
		}

		// now recurse on the children if both have children
		// note - it won't go here is !success because still isLeaf()
		if (!isLeaf() && !other->isLeaf()) {
			success = getLeftChild()->_reshapeToUnion(
					other->getLeftChild(), minChildPoints, errorFilename);
			success = getRightChild()->_reshapeToUnion(
					other->getRightChild(), minChildPoints, errorFilename)
					&& success;
		}
		
	}
	return success;
}

// used for MCMC log posteriors(JUNE 2012 changes)
// Get unscaled log likelihood over tree rooted at this 
void SPSnode::_getUnscaledTreeLogLik(dotprecision& nlogn, 
										int& ndepth, int depth) const
{
	
	// if counter is 0 there is no log likelihood
	if (counter > 0) {

		if (isLeaf()) {
			
			// accumulate counter * log counter to nlogn
			accumulate(nlogn, 1.0*counter, log(1.0*counter));
			// accumulate counter * depth to ndepth
			ndepth += counter * depth;
			
		}
		
		else {
			
			getLeftChild()->_getUnscaledTreeLogLik(nlogn, ndepth, depth+1);
			getRightChild()->_getUnscaledTreeLogLik(nlogn, ndepth, depth+1);
			
		}
	}

}


// L1 distance between this and another node, with the root counters for normalisers
// no checks on boxes since this should done by another function that calls this one
// no checks null, empty extra - will just explode if there is a problem
cxsc::dotprecision& SPSnode::_getL1distance(cxsc::dotprecision& disL1,
										const SPSnode * const other,
										std::size_t thisBigN,
										std::size_t otherBigN) const
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
			std::size_t my_n = getCounter();
			other->nodeL1Distance(disL1, otherBigN, my_n, thisBigN);
		}
		
		else { // I am not a leaf, so recurse
		
			#ifdef DEBUG_L1
				std::cout << "I am not a leaf: recursing" << std::endl;
			#endif
			
			getLeftChild()->_getL1distance(disL1, other->getLeftChild(),
												thisBigN, otherBigN);
			getRightChild()->_getL1distance(disL1, other->getRightChild(),
											thisBigN, otherBigN);
		}
	}
	
	else { // other is a leaf
	
		#ifdef DEBUG_L1
			std::cout << "Other is a leaf, other name is " << other->getNodeName() << std::endl;
		#endif
		
		std::size_t other_n = other->getCounter();
		nodeL1Distance(disL1, thisBigN, other_n, otherBigN);
	
	}
	
	return disL1;
}


cxsc::dotprecision& SPSnode::nodeL1Distance(cxsc::dotprecision& disL1,
								std::size_t thisBigN, 
								std::size_t other_n, std::size_t otherBigN) const
{
	#ifdef DEBUG_L1
		std::cout << "\nIn nodeL1distance, I am " << getNodeName() << ", and disL1 is " << cxsc::rnd(disL1) << std::endl;
		std::cout << "thisBigN is " << thisBigN << ", other_n is " << other_n << ", otherBigN " << otherBigN << std::endl;
	#endif
		
	// this is not a leaf
	if (!isLeaf() ) {
		
		#ifdef DEBUG_L1
			std::cout << "I am not a leaf" << std::endl;
		#endif
		
		disL1 = getLeftChild()->nodeL1Distance(disL1, thisBigN, 
											other_n, (otherBigN * 2) );
		disL1 = getRightChild()->nodeL1Distance(disL1, thisBigN, 
											other_n, (otherBigN * 2) );
	}
	else { // this is a leaf
	
		cxsc::dotprecision other(0.0);
		cxsc::dotprecision me(0.0);
		accumulate(other, (1.0 * other_n), (1.0/otherBigN) );
		accumulate(me, (1.0 * getCounter()), (1.0/thisBigN) );
		
		#ifdef DEBUG_L1
			std::cout << "I am a leaf, my normalised count is " << cxsc::rnd(me) << ", and other normalised count is " << cxsc::rnd(other) << std::endl;
			std::cout << "adding " << cxsc::rnd( cxsc::abs(me - other) ) << " to disL1" << std::endl;
		#endif
		
		disL1 += cxsc::abs(me - other);
		
		#ifdef DEBUG_L1
			std::cout << "disL1 is now " << cxsc::rnd( disL1 ) << std::endl;
		#endif
		
	}
	
	#ifdef DEBUG_L1
			std::cout << "returning from nodeL1distance, disL1 = " << cxsc::rnd( disL1 ) << std::endl;
	#endif
	
	return disL1;
}
	
// end of definitions for SPSnode class

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::SPSnode & s1, 
			subpavings::SPSnode & s2) // throw ()
{
	s1.swapSPS(s2);
}


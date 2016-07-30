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

/*!/ \file
\brief BooleanValueMappedSPnode definitions.
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "booleanvaluemappedspnode.hpp"

using namespace subpavings;
using namespace std;

				
// ------------------------ public member functions ----------------

// Destructor.
BooleanValueMappedSPnode::~BooleanValueMappedSPnode()  {}


// no-argument constructor,
BooleanValueMappedSPnode::BooleanValueMappedSPnode() 
	: MappedSPnode<BooleanMappedValue>() {}   
	// uses the base SPnode class default constructor


// initialised constructor
// initialised with a box
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const ivector& v)
	: 	MappedSPnode<BooleanMappedValue>(v, BooleanMappedValue(false))
	 {}


 // initialised constructor
// initialised with a labeled box
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const LabBox& lb)
	: 	MappedSPnode<BooleanMappedValue>(lb, BooleanMappedValue(false))
	 {}



// initialised constructor
// initialised with a box, a value for the range,
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const ivector& v, 
											bool range)
	: MappedSPnode<BooleanMappedValue>(v, BooleanMappedValue(range)) {}



// initialised constructor
// initialised with a labeled box and a value for range
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const LabBox& lb, 
											bool range)
	: MappedSPnode<BooleanMappedValue>(lb, BooleanMappedValue(range)) {}

// initialised constructor
// initialised with a box, a bv value for the range,
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const ivector& v, 
						const BooleanMappedValue& range)
	: MappedSPnode<BooleanMappedValue>(v, range) {}



// initialised constructor
// initialised with a labeled box and a bv value for range
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const LabBox& lb, 
						const BooleanMappedValue& range)
	: MappedSPnode<BooleanMappedValue>(lb, range) {}



// Copy constructor
// copies from given node downwards
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const SPnode& other)
	: MappedSPnode<BooleanMappedValue>()
{
	
	if (!other.isEmpty()) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = BooleanMappedValue(false);
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new BooleanValueMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new BooleanValueMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
BooleanValueMappedSPnode::BooleanValueMappedSPnode(const BooleanValueMappedSPnode& other)
	: MappedSPnode<BooleanMappedValue>()
{
	
	if (other.theBox != NULL) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new BooleanValueMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new BooleanValueMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
BooleanValueMappedSPnode::BooleanValueMappedSPnode(
					const MappedSPnode<BooleanMappedValue>& other)
	: MappedSPnode<BooleanMappedValue>()
{
	if (!other.isEmpty()) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new BooleanValueMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new BooleanValueMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;
}


// Copy assignment operator
// copies from given node downwards
BooleanValueMappedSPnode& BooleanValueMappedSPnode::operator=(
							BooleanValueMappedSPnode rhs)
{
	rhs.swapBMSP(*this); // make sure we use our version of swap
	return(*this);
}

// Copy assignment operator
// copies from given node downwards
BooleanValueMappedSPnode& BooleanValueMappedSPnode::operator=(
							MappedSPnode<BooleanMappedValue> rhs)
{
	rhs.swapMSPSR(*this); // make sure we use our version of swap
	return(*this);
}

void BooleanValueMappedSPnode::replaceMe(MappedSPnode<BooleanMappedValue> newNode)
{
	
	BooleanValueMappedSPnode* pNode = getParent();
	newNode.swapMSPSR(*this);
	
	if (NULL != pNode) {
		
		this->parent = pNode;
	}
}

void BooleanValueMappedSPnode::replaceMe(BooleanValueMappedSPnode newNode)
{
	
	BooleanValueMappedSPnode* pNode = getParent();
	newNode.swapBMSP(*this);
	
	if (NULL != pNode) {
		
		this->parent = pNode;
	}
}


// parent and child accessors have to hide the base class implementation
// this is not good but otherwise we get the base class return type
// I've asked around and I can't find a way around it ...

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
BooleanValueMappedSPnode* BooleanValueMappedSPnode::getParent() const
{ return (BooleanValueMappedSPnode*) parent; }


// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node, cast to this type
BooleanValueMappedSPnode* BooleanValueMappedSPnode::getLeftChild() const
{ return (BooleanValueMappedSPnode*) leftChild; }


// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node, cast this type
BooleanValueMappedSPnode* BooleanValueMappedSPnode::getRightChild() const
{ return (BooleanValueMappedSPnode*) rightChild; }

bool BooleanValueMappedSPnode::getRange() const
{
	return range.getValue();
}

// Find where data should be 
// childInd is an indicator for which child is being checked
// throws exception if node has no box
const BooleanValueMappedSPnode* BooleanValueMappedSPnode::findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd) const
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"IntervalMappedSPnode::findContainingNode(const cxsc::rvector&, OPERATIONS_ON)");
	}
	// start at the top
	const BooleanValueMappedSPnode* retObj = NULL;
	
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
void BooleanValueMappedSPnode::nodeExpand(int comp)
{
	// can only expand if there is a box
	if (NULL == theBox) {
		throw NoBox_Error("BooleanValueMappedSPnode::nodeExpand(int)");
	}

	// only do something if this node is a leaf
	if (isLeaf()) {
		
		BooleanValueMappedSPnode* newLC = NULL;
		BooleanValueMappedSPnode* newRC = NULL;
		
		try {

			// ivectors to become boxes for new children
			ivector lC, rC;
			// Call Lower() and Upper() to put split boxes
			// into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);

			// make and add the new children
			newLC = new BooleanValueMappedSPnode(lC, range);
			newRC = new BooleanValueMappedSPnode(rC, range);
			
			nodeAddLeft(newLC);
			nodeAddRight(newRC);
			// both children get the same range as this
			
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
void BooleanValueMappedSPnode::nodeExpand()
{
	SPnode::nodeExpand();
}

// nodeReabsorbChildren() can use the base class implementation
// (the range for this should be correct so just delete the children)


cxsc::real BooleanValueMappedSPnode::getTotalLeafTrueAreaOnBox() const
{
	if (isLeaf()) return getTrueAreaOnBox();
	else {
		return ( getLeftChild()->getTotalLeafTrueAreaOnBox()
			+ getRightChild()->getTotalLeafTrueAreaOnBox() );
	}
}

/*Volume of box represented by this multiplied by
diameter of interval range of this.*/
cxsc::real BooleanValueMappedSPnode::getTrueAreaOnBox() const
{
	cxsc::real retr(0.0);
	
	if (range.getValue()) {
		retr = nodeRealVolume();
	}
	return retr;
}


/* Return a reference to a container of nodes.
 
Contents of container are the leaves descended from this, 
or this if this is a leaf, left to right order. */
BooleanValueMappedSPnode::Ptrs& BooleanValueMappedSPnode::getLeaves(
			BooleanValueMappedSPnode::Ptrs& leaves)
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
BooleanValueMappedSPnode::ConstPtrs& BooleanValueMappedSPnode::getConstLeaves(
			BooleanValueMappedSPnode::ConstPtrs& leaves) const
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
BooleanValueMappedSPnode::Ptrs& BooleanValueMappedSPnode::getSubLeaves(
			BooleanValueMappedSPnode::Ptrs& subleaves)
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
BooleanValueMappedSPnode::ConstPtrs& BooleanValueMappedSPnode::getConstSubLeaves(
			BooleanValueMappedSPnode::ConstPtrs& subleaves) const
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


void BooleanValueMappedSPnode::swapBMSP(BooleanValueMappedSPnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
	rightChild and nodeName are inherited from base class */
	SPnode::swap(spn); // use the base version
	std::swap(range, spn.range); 
}



//Output for all the  leaf boxes in this with true values, using tab delimiters
std::ostream& BooleanValueMappedSPnode::leavesOutputTabsTrue(std::ostream &os) const
{
	// uses  member function leafOutputTabs to generate node output
	if (range.getValue() && !(isEmpty()) && isLeaf() ) { // this is a non-empty leaf
		leafOutputTabs(os);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabsTrue(os);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabsTrue(os);
	}
    return os;

}

//Output for all the  leaf boxes with true values in this, using tab delimiters
std::ostream& BooleanValueMappedSPnode::leavesOutputTabsTrue(
									std::ostream &os, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	leavesOutputTabsTrue(os);
	os << cxsc::RestoreOpt;
	
	return os;

}

#if(0)
/* Give this a range equal to range of this AND range of other.*/
void BooleanValueMappedSPnode::_addRanges(
						const BooleanValueMappedSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"BooleanValueMappedSPnode::_addRanges(const BooleanValueMappedSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_addRanges(
									other->getLeftChild());
		getRightChild()->_addRanges(
									other->getRightChild());

	}

	// deal with this range collection
	range = (range || other->range);
}

/* Give this a range equal to range of this XOR (symmetric set differnce)
 * range of other,
 * ie set difference. */
void BooleanValueMappedSPnode::_subtractRanges(
						const BooleanValueMappedSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"BooleanValueMappedSPnode::_subtractRanges(const BooleanValueMappedSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_subtractRanges(other->getLeftChild());
		getRightChild()->_subtractRanges(other->getRightChild());

	}

	// deal with this range collection
	if (range && other->range) range = false;
	else range = (range || other->range);
}


/* Give this a range equal to range of this AND range of other.*/
void BooleanValueMappedSPnode::_multRanges(
							const BooleanValueMappedSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"BooleanValueMappedSPnode::_multRanges(const BooleanValueMappedSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_multRanges(other->getLeftChild());
		getRightChild()->_multRanges(other->getRightChild());

	}

	// deal with this range collection
	range = (range && other->range);
}

/* Give this a range equal to range of this set difference range of other.*/
void BooleanValueMappedSPnode::_divRanges(
						const BooleanValueMappedSPnode * const other)
{
	if (other == NULL) {
		throw NullSubpavingPointer_Error(
			"BooleanValueMappedSPnode::_divRanges(const BooleanValueMappedSPnode * const)");
	}

	// recurse on the children if any first
	if (!isLeaf()) {
		getLeftChild()->_divRanges(	other->getLeftChild() );
		getRightChild()->_divRanges( other->getRightChild() );

	}

	// deal with this range collection
	if (!range || (other->range)) range = false;
}




/* Make a non-minimal union of subpavings using union (OR) of ranges.
*/
void BooleanValueMappedSPnode::_addNonMinimalUnion(
			   const BooleanValueMappedSPnode& rhs)
{
	
	//if (rhs == NULL || rhs->isEmpty()) do nothing
	
	// should not need to check on rhs here 
	if ((isEmpty()) && !rhs.isEmpty()) {

		*this = BooleanValueMappedSPnode(rhs);
	}
			
	if (!isEmpty() && !rhs.isEmpty()) {

		// make copies of rhs to reshape

		BooleanValueMappedSPnode rhsTemp(rhs);

		// reshape
		_reshapeTreesToUnion(&rhsTemp);

		_addRanges(&rhsTemp);

	}
}

/* Make a non-minimal union of subpavings using
  XOR = symmetric difference of ranges.*/
void BooleanValueMappedSPnode::_subtractNonMinimalUnion(
			   const BooleanValueMappedSPnode& rhs)
{
	// should not need to check on rhs not null here 
	if ((isEmpty()) || rhs.isEmpty()) {

		throw IncompatibleDimensions_Error(
			"BooleanValueMappedSPnode::_subtractNonMinimalUnion(const BooleanValueMappedSPnode * const)");
	}
	
	// make copies of rhs to reshape

	BooleanValueMappedSPnode rhsTemp(rhs);

	// reshape
	_reshapeTreesToUnion(&rhsTemp);

	_subtractRanges(&rhsTemp);

}

/* Make a non-minimal union of subpavings using
intersection of ranges.*/
void BooleanValueMappedSPnode::_multiplyNonMinimalUnion(
			   const BooleanValueMappedSPnode& rhs)
{
	// should not need to check on rhs not null here 
	if ((isEmpty()) || rhs.isEmpty()) {

		throw IncompatibleDimensions_Error(
			"BooleanValueMappedSPnode::_multiplyNonMinimalUnion(const BooleanValueMappedSPnode * const)");
	}
	
	// make copies of rhs to reshape

	BooleanValueMappedSPnode rhsTemp(rhs);

	// reshape
	_reshapeTreesToUnion(&rhsTemp);

	_multRanges(&rhsTemp);

}

/* Make a non-minimal union of subpavings using
set difference of ranges.*/
void BooleanValueMappedSPnode::_divideNonMinimalUnion(
			   const BooleanValueMappedSPnode& rhs)
{
	if ((isEmpty()) || rhs.isEmpty()) {

		throw IncompatibleDimensions_Error(
			"BooleanValueMappedSPnode::_divideNonMinimalUnion(const BooleanValueMappedSPnode * const)");
	}
	
	// make copies of rhs to reshape

	BooleanValueMappedSPnode rhsTemp(rhs);

	// reshape
	_reshapeTreesToUnion(&rhsTemp);

	_divRanges(&rhsTemp);

}



/*Union range collection of this.
*/
void BooleanValueMappedSPnode::_scalarAdd(bool add)
{
	
	if (!isEmpty()) {
		range = (range || add);
		// recurse on children
		if (hasLCwithBox())
			getLeftChild()->_scalarAdd(add);
		if (hasRCwithBox())
			getRightChild()->_scalarAdd(add);
	}
}

/*XOR = symmetric set difference range collection of this.*/
void BooleanValueMappedSPnode::_scalarSubtract(bool sub)
{
	
	if (!isEmpty()) {
		if ((range && sub)) range = false;
		else range = (range || sub);
		// recurse on children
		if (hasLCwithBox())
			getLeftChild()->_scalarSubtract(sub);
		if (hasRCwithBox())
			getRightChild()->_scalarSubtract(sub);
	}
}


/*Intersect range of this.*/
void BooleanValueMappedSPnode::_scalarMult(bool mult)
{
	
	if (!isEmpty()) {
		range = range && mult;
		// recurse on children
		if (hasLCwithBox())
			getLeftChild()->_scalarMult(mult);
		if (hasRCwithBox())
			getRightChild()->_scalarMult(mult);
	}
}

/*set difference */
void BooleanValueMappedSPnode::_scalarDiv(bool div)

{
	
	if (!isEmpty()) {
		if (!range || div) range = false;
		// recurse on children
		if (hasLCwithBox())
			getLeftChild()->_scalarDiv(div);
		if (hasRCwithBox())
			getRightChild()->_scalarDiv(div);
	}
}

#endif

// Print the details of a single leaf node, using tab delimiters.*/
std::ostream& BooleanValueMappedSPnode::leafOutputTabs(std::ostream &os) const
{
	if(theBox != NULL) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox

		// output the nodeName, nodeVolume
		os << nodeName;
		real vol = nodeRealVolume();
		os << "\t" << vol;

		// followed by the range
		os << "\t" << range;

		// followed by intervals making up box using Inf & Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i])
				<< "\t" << Sup(thisBox[i]);
		}
		os << std::flush;
		
	}
	return os;
}

// Print the details of a specific node in a subpaving.*/
std::ostream& BooleanValueMappedSPnode::nodePrint(std::ostream &os) const
{
	// output for box in form:
	// box, volume, summary data

	if(theBox != NULL) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox

		os << nodeName << "\tBox is:\t";

		for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
			// c-xsc default output for intervals
			os << "\t" << thisBox[i];   }

		os << "\tBox volume is:\t" << nodeVolume();
		os << "\trange data:\t" << range << std::endl;
		
	}
	return os;
}


std::ostream& BooleanValueMappedSPnode::oneLineOutput(std::ostream& os, int level) const
{
	// do me
	for (int i = 0; i < level; ++i) { os << "\t"; }
	os << getNodeName() << "\tRange: " << range 
								<< "\tbox: " << getBox() << std::endl;
	
	// do children
	if (hasLCwithBox()) getLeftChild()->oneLineOutput(os, level+1);
	if (hasRCwithBox()) getRightChild()->oneLineOutput(os, level+1);
	
	return os;
}


// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::BooleanValueMappedSPnode & s1, 
			subpavings::BooleanValueMappedSPnode & s2) // throw ()
	{
		s1.swapBMSP(s2);
	}


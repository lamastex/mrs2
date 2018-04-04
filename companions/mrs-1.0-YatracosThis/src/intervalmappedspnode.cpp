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
\brief IntervalMappedSPnode definitions.
*/


#include "intervalmappedspnode.hpp"

using namespace subpavings;
using namespace std;

//#define MYDEBUGHULL 


// ------------------------ public member functions ----------------

// Destructor.
IntervalMappedSPnode::~IntervalMappedSPnode()  {}


// no-argument constructor,
IntervalMappedSPnode::IntervalMappedSPnode() 
	: MappedSPnode<cxsc::interval>() {}   
	// uses the base SPnode class default constructor


// initialised constructor
// initialised with a box
IntervalMappedSPnode::IntervalMappedSPnode(const ivector& v)
	: 	MappedSPnode<cxsc::interval>(v, cxsc::interval(0.0,0.0))
	 {}


 // initialised constructor
// initialised with a labeled box
IntervalMappedSPnode::IntervalMappedSPnode(const LabBox& lb)
	: 	MappedSPnode<cxsc::interval>(lb, cxsc::interval(0.0,0.0))
	 {}



// initialised constructor
// initialised with a box, a value for the range,
IntervalMappedSPnode::IntervalMappedSPnode(const ivector& v, 
											const cxsc::interval& range)
	: MappedSPnode<cxsc::interval>(v, range) {}



// initialised constructor
// initialised with a labeled box and a value for range
IntervalMappedSPnode::IntervalMappedSPnode(const LabBox& lb, 
											const cxsc::interval& range)
	: MappedSPnode<cxsc::interval>(lb, range) {}


// Copy constructor
// copies from given node downwards
IntervalMappedSPnode::IntervalMappedSPnode(const SPnode& other)
{
	
	if (!other.isEmpty()) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = cxsc::interval(0.0,0.0);
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new IntervalMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new IntervalMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
IntervalMappedSPnode::IntervalMappedSPnode(const IntervalMappedSPnode& other)
{
	
	if (other.theBox != NULL) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new IntervalMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new IntervalMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;

}

// Copy constructor
// copies from given node downwards
IntervalMappedSPnode::IntervalMappedSPnode(
					const MappedSPnode<cxsc::interval>& other)

{
	if (!other.isEmpty()) {
		theBox = new ivector( other.getBox() );
	}
	nodeName = other.getNodeName();
	range = other.getRange();
	
	//recursion on the children
	if (other.hasLCwithBox()) {
		nodeAddLeft(new IntervalMappedSPnode(
			*(other.getLeftChild())));
	}
	else leftChild=NULL;

	if (other.hasRCwithBox()) {
		nodeAddRight(new IntervalMappedSPnode(
			*(other.getRightChild())));
	}
	else rightChild=NULL;
}


// Copy assignment operator
// copies from given node downwards
IntervalMappedSPnode& IntervalMappedSPnode::operator=(
							IntervalMappedSPnode rhs)
{
	rhs.swapIMSP(*this); // make sure we use our version of swap
	return(*this);
}

// Copy assignment operator
// copies from given node downwards
IntervalMappedSPnode& IntervalMappedSPnode::operator=(
							MappedSPnode<cxsc::interval> rhs)
{
	rhs.swapMSPSR(*this); // make sure we use our version of swap
	return(*this);
}



// parent and child accessors have to hide the base class implementation
// this is not good but otherwise we get the base class return type
// I've asked around and I can't find a way around it ...

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
IntervalMappedSPnode* IntervalMappedSPnode::getParent() const
{ return (IntervalMappedSPnode*) parent; }


// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node, cast to this type
IntervalMappedSPnode* IntervalMappedSPnode::getLeftChild() const
{ return (IntervalMappedSPnode*) leftChild; }


// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node, cast this type
IntervalMappedSPnode* IntervalMappedSPnode::getRightChild() const
{ return (IntervalMappedSPnode*) rightChild; }

// comparison operator, returns true iff range this < rhs.range
bool IntervalMappedSPnode::operator<(
		const IntervalMappedSPnode& rhs) const
{
	return (getRange() < rhs.getRange());
}

/*Self-scalar addition operator*/
IntervalMappedSPnode& IntervalMappedSPnode::operator+= (const real& val)
{
	MappedSPnode<cxsc::interval>::operator+=(cxsc::interval(val));
	return *this;
}

/*Scalar addition operator.*/
const IntervalMappedSPnode IntervalMappedSPnode::operator+ (const real& val) const
{
	return this->MappedSPnode<cxsc::interval>::operator+(cxsc::interval(val));;
}

/*Self-scalar subtraction operator*/
IntervalMappedSPnode& IntervalMappedSPnode::operator-= (const real& val)
{
	MappedSPnode<cxsc::interval>::operator-=(cxsc::interval(val));
	return *this;
}

/*Scalar subtraction operator.*/
const IntervalMappedSPnode IntervalMappedSPnode::operator- (const real& val) const
{
	return this->MappedSPnode<cxsc::interval>::operator-(cxsc::interval(val));;
}

/*Self-scalar multiplication operator*/
IntervalMappedSPnode& IntervalMappedSPnode::operator*= (const real& val)
{
	MappedSPnode<cxsc::interval>::operator*=(cxsc::interval(val));
	return *this;
}

/*Scalar multiplication operator.*/
const IntervalMappedSPnode IntervalMappedSPnode::operator* (const real& val) const
{
	return this->MappedSPnode<cxsc::interval>::operator*(cxsc::interval(val));;
}

/*Self-scalar division operator*/
IntervalMappedSPnode& IntervalMappedSPnode::operator/= (const real& val)
{
	MappedSPnode<cxsc::interval>::operator/=(cxsc::interval(val));
	return *this;
}

/*Scalar division operator.*/
const IntervalMappedSPnode IntervalMappedSPnode::operator/ (const real& val) const
{
	return this->MappedSPnode<cxsc::interval>::operator/(cxsc::interval(val));;
}

// Find where data should be 
// childInd is an indicator for which child is being checked
// throws exception if node has no box
const IntervalMappedSPnode* IntervalMappedSPnode::findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd) const
{
	if ( isEmpty() ) {
		throw NoBox_Error(
		"RealMappedSPnode::findContainingNode(const cxsc::rvector&, OPERATIONS_ON)");
	}
	// start at the top
	const IntervalMappedSPnode* retObj = NULL;
	
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
void IntervalMappedSPnode::nodeExpand(int comp)
{
	// can only expand if there is a box
	if (NULL == theBox) {
		throw NoBox_Error("IntervalMappedSPnode::nodeExpand(int)");
	}

	// only do something if this node is a leaf
	if (isLeaf()) {
		
		IntervalMappedSPnode* newLC = NULL;
		IntervalMappedSPnode* newRC = NULL;
		
		try {
             
			// ivectors to become boxes for new children
			ivector lC, rC;
			// Call Lower() and Upper() to put split boxes
			// into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);
	
			// make and add the new children
			newLC = new IntervalMappedSPnode(lC, range);
			newRC = new IntervalMappedSPnode(rC, range);

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
void IntervalMappedSPnode::nodeExpand()
{
	int maxdiamcomp; // variable to hold first longest dimension
	double temp = ::MaxDiam(getBox(), maxdiamcomp);
	nodeExpand(maxdiamcomp); // complete nodeExpand
}

// nodeReabsorbChildren() can use the base class implementation
// (the range for this will be correct so just delete the children)



/* works bottom up to calculate the 
  range for this node as the interval hull of the ranges of the children.*/
void IntervalMappedSPnode::hullPropagation()
{
	#ifdef MYDEBUGHULL
		std::cout << "In hullPropagation, I am " << nodeName << std::endl;
	#endif 
	// first recursively deal with the children of the children
	if (hasLCwithBox())
		getLeftChild()->hullPropagation();
	if (hasRCwithBox())
		getRightChild()->hullPropagation();

	// now deal with this
	if (hasLCwithBox() && hasRCwithBox()) {
		
		#ifdef MYDEBUGHULL
			std::cout << "Back in hullPropagation for " << nodeName << std::endl;
			std::cout << "my range is " << range << std::endl;
			std::cout << "getLeftChild()->range is " << (getLeftChild()->range) << std::endl;
			std::cout << "getRightChild()->range is " << (getRightChild()->range) << std::endl;
			std::cout << "(getLeftChild()->getRange() | getRightChild()->getRange()) is " << ((getLeftChild()->getRange() | getRightChild()->getRange())) << std::endl;
		
		#endif 
			
		// hull
		range = (getLeftChild()->getRange() | getRightChild()->getRange());

	}
}



// overide base class
void IntervalMappedSPnode::slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts)
{
	
	#ifdef SLICE_OUTPUT
		std::cout << "In IntervalMappedSPnode::slice, I am " << getNodeName() << std::endl;
	#endif

	std::vector<cxsc::real> fullSlicePts = sliceCheck(sliceDims, slicePts);
	
	IntervalMappedSPnode temp(*this);
	temp._slice(sliceDims, fullSlicePts);
	swapIMSP(temp);

}


/*Diameter of interval range of this.*/
cxsc::real IntervalMappedSPnode::getRangeDiameter() const
{
	cxsc::real retr(0.0);
	
	if (!cxsc::IsEmpty(range)) {
		retr = cxsc::diam(range);
	}
	return retr;
}

cxsc::real IntervalMappedSPnode::getTotalLeafIntervalAreaOnBox() const
{
	if (isLeaf()) return getIntervalAreaOnBox();
	else {
		return ( getLeftChild()->getTotalLeafIntervalAreaOnBox()
			+ getRightChild()->getTotalLeafIntervalAreaOnBox() );
	}
}

/*Volume of box represented by this multiplied by
diameter of interval range of this.*/
cxsc::real IntervalMappedSPnode::getIntervalAreaOnBox() const
{
	cxsc::real retr(0.0);
	
	if (!cxsc::IsEmpty(range)) {
		retr = getRangeDiameter() * nodeRealVolume();
	}
	return retr;
}

/*Volume of box represented by this multiplied by
diameter of interval range of this, less sum of the same for the children.
Returns interval area of this if no children.*/
cxsc::real IntervalMappedSPnode::getIntervalAreaDiffToChildren() const
{
	cxsc::real retr(0.0);
	
	if (isLeaf() ) {
		retr = getIntervalAreaOnBox();
	}
	else {
		retr = nodeRealVolume() * (getRangeDiameter() 
				- 0.5*(getLeftChild()->getRangeDiameter() 
						+ getRightChild()->getRangeDiameter() ) );
	}
	
	return retr;
}

/* Return a reference to a container of nodes.
 
Contents of container are the leaves descended from this, 
or this if this is a leaf, left to right order. */
IntervalMappedSPnode::Ptrs& IntervalMappedSPnode::getLeaves(
			IntervalMappedSPnode::Ptrs& leaves)
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
IntervalMappedSPnode::ConstPtrs& IntervalMappedSPnode::getConstLeaves(
			IntervalMappedSPnode::ConstPtrs& leaves) const
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
IntervalMappedSPnode::Ptrs& IntervalMappedSPnode::getSubLeaves(
			IntervalMappedSPnode::Ptrs& subleaves)
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
IntervalMappedSPnode::ConstPtrs& IntervalMappedSPnode::getConstSubLeaves(
			IntervalMappedSPnode::ConstPtrs& subleaves) const
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


void IntervalMappedSPnode::swapIMSP(IntervalMappedSPnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
	rightChild and nodeName are inherited from base class */
	SPnode::swap(spn); // use the base version
	std::swap(range, spn.range); // use our specialisation of swap
}



// Print the details of a single leaf node, using tab delimiters.*/
std::ostream& IntervalMappedSPnode::leafOutputTabs(std::ostream &os) const
{
	if(theBox != NULL) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox

		// output the nodeName, nodeVolume
		os << nodeName;
		real vol = nodeRealVolume();
		os << "\t" << vol;

		// followed by the range
		os << "\t" << cxsc::Inf(range) << "\t" << cxsc::Sup(range);

		// followed by intervals making up box using Inf & Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i]) << "\t" << Sup(thisBox[i]);
		}
		os << std::flush;
		
	}
	return os;
}





// Print the details of a specific node in a subpaving.*/
std::ostream& IntervalMappedSPnode::nodePrint(std::ostream &os) const
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
		os << "\trange data:\t" << cxsc::Inf(range) << "\t" 
								<< cxsc::Sup(range) 
								<< std::endl;
		
	}
	return os;
}


std::ostream& IntervalMappedSPnode::oneLineOutput(std::ostream& os, int level) const
{
	// do me
	for (int i = 0; i < level; ++i) { os << "\t"; }
	os << getNodeName() << "\tRange: " << cxsc::Inf(range) << "\t" 
								<< cxsc::Sup(range) 
								<< "\tbox: " << getBox() << std::endl;
	
	// do children
	if (hasLCwithBox()) getLeftChild()->oneLineOutput(os, level+1);
	if (hasRCwithBox()) getRightChild()->oneLineOutput(os, level+1);
	
	return os;
}


// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::IntervalMappedSPnode & s1, 
			subpavings::IntervalMappedSPnode & s2) // throw ()
	{
		s1.swapIMSP(s2);
	}

template <>
void std::swap(cxsc::interval & i1, 
			cxsc::interval & i2) // throw ()
{
	cxsc::real saveInf1 = cxsc::Inf(i1);
	cxsc::real saveSup1 = cxsc::Sup(i1);
	SetInf(i1, cxsc::Inf(i2));
	SetSup(i1, cxsc::Sup(i2));
	SetInf(i2, saveInf1);
	SetSup(i2, saveSup1);
}


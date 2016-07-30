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
\brief SPnode (SubPaving) and associated non-member function definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "spnode.hpp"

#include "subpaving_exception.hpp"
#include "toolz.hpp"
#include "sptools.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <stdexcept>

// include fstream so as to be able to output a file from spImage
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>

#include <algorithm> // for algorithms wth stl containers

#include <iterator> // ostream iterator 


#include <cassert> // assertions

//#define DEBUG_RESHAPE
//#define DEBUG_PARTITIONS

#ifdef NDEBUG
	#undef DEBUG_RESHAPE
	#undef DEBUG_PARTITIONS
#endif

using namespace std;
using namespace subpavings;


// ----------------------- SPnode class definitions ------------------

// ---------------- public methods

// default constructor
SPnode::SPnode() :  theBox(NULL),
		nodeName("X"), parent(NULL), leftChild(NULL), rightChild(NULL)
	{}


// initialised constructor, initialised with one ivector for the box
SPnode::SPnode(const ivector& v) : theBox(NULL),
		nodeName("X"), parent(NULL), leftChild(NULL), rightChild(NULL)
{
	try {
		int d = Ub(v) - Lb(v) + 1; 
		if (d < 1) {
			throw MalconstructedBox_Error(
								"SPnode::SPnode(const ivector&, int)");
		}
		theBox = new ivector(v);
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}

// initialised constructor, initialised with a LabBox (labeled box)
SPnode::SPnode(const LabBox& lb) : theBox(NULL), 
	nodeName("X"), parent(NULL), leftChild(NULL), rightChild(NULL) 
{
	try {
		ivector v = lb.Box;
		int d = Ub(v) - Lb(v) + 1; 
		if (d < 1) {
			throw MalconstructedBox_Error(
						"SPnode::SPnode(const LabBox&)");
		}
		theBox = new ivector(v);
	}
	catch (exception const& e) {
		constructor_error_handler();
	}

}


//Copy constructor
//copies from the given node downwards
SPnode::SPnode(const SPnode& other) : 
		theBox(NULL), nodeName(other.nodeName),
		parent(NULL), 
		leftChild(NULL), rightChild(NULL)
		
{
	
	try {
		
		if (other.theBox != NULL) {
			
			theBox=new ivector(*other.theBox);
		}

		//recursion on the children
		if (other.leftChild) {
			leftChild=new SPnode(*other.leftChild);
			leftChild->parent = this;
		}
		else leftChild=NULL;

		if (other.rightChild) {
			rightChild=new SPnode(*other.rightChild);
			rightChild->parent = this;
		}
		else rightChild=NULL;
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}


// Destructor.
SPnode::~SPnode()
{
	try {
		
		delete theBox;
		theBox = NULL;
		delete leftChild;
		leftChild = NULL;
		delete rightChild;
		rightChild = NULL;
	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in SPnode destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}


//copy assignment operator
//copies from this node downwards
// no check for self-assignment - copy elision and copy-and-swap idiom
SPnode& SPnode::operator=(SPnode rhs)
{
	rhs.swap(*this);
	return *this;
}

// recursively rename from the top down
void SPnode::recursiveRename()
{
	if(getLeftChild() != NULL) {
		getLeftChild()->setNodeName(nodeName + "L");
		getLeftChild()->recursiveRename();
	}
	if(getRightChild() != NULL) {
		getRightChild()->setNodeName(nodeName + "R");
		getRightChild()->recursiveRename();
	}
}

/* Return boolean to indicate if node is splittable.
 * 
A node is splittable if 
* there is at least one value in the basic number screen 
* between the inf and sup of the interval on the coordinate
* with maximum diameter (ie the box can be split on that interval)
* and the node volume is >= 2 * cxsc::MinReal (the
smallest representable real number).*/
bool SPnode::isSplittableNode() const
{
	ivector box = getBox();
	interval maxD = box[MaxDiamComp(box)];
	return ( (succ(succ(Inf(maxD))) <= Sup(maxD))
		&& (realVolume(box) >= 2*cxsc::MinReal) );
}

// Accessor for the dimension of theBox of a node.
size_t SPnode::getDimension() const
{
	int d = 0;
	if (!isEmpty()) {
		ivector box = getBox(); 
		d = Ub(*theBox) - Lb(*theBox) + 1; 
	}
	return d;
}

// Accessor for theBox of a node.
// Returns a copy of the object pointed to by theBox of a node.
ivector SPnode::getBox() const
{
	
	// check there is a box to get
	if (NULL == theBox) {
		throw NoBox_Error("SPnode::getBox()");
	}
	return *theBox;
	
}

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
SPnode* SPnode::getParent() const
{ return parent; }

// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node.
SPnode* SPnode::getLeftChild() const
{ return leftChild; }

// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node.
SPnode* SPnode::getRightChild() const
{ return rightChild; }



// Get the node name.
std::string SPnode::getNodeName() const
{
	return nodeName;
}

// Set the node name.
void SPnode::setNodeName(std::string newname)
{
	nodeName = newname;
}


// Return the volume of the box as a double.
double SPnode::nodeVolume() const
{ 
	if (isEmpty() ) {
		throw NoBox_Error("SPnode::nodeVolume()");
	}
	
	return Volume(getBox());
}
	
// Return the volume of the box as a real.
real SPnode::nodeRealVolume() const
{ 
	if (isEmpty() ) {
		throw NoBox_Error("SPnode::nodeRealVolume()");
	}
	return realVolume(getBox()); }

// find if this node is a subleaf node
// a sub-leaf node is the parent of leaf nodes and only have leaf nodes
// as children
bool SPnode::isSubLeaf() const
{

	bool retValue = false;
	if (hasLCwithBox() && hasRCwithBox()) { // has two children
		if (getLeftChild()->isLeaf() && getRightChild()->isLeaf()) {
			retValue = true;
		}
	}
	else if (hasLCwithBox()) { // only has left child
		if (getLeftChild()->isLeaf()) retValue = true;
	}
	else if (hasRCwithBox()) { // only has right child
			if (getRightChild()->isLeaf()) retValue = true;
	}
	// default return value false

	return retValue;
}

// Check if this SPnode is empty in the sense of having no box
bool SPnode::isEmpty() const
{return ( NULL == theBox) ; }

// Check if this SPnode is a leaf.
bool SPnode::isLeaf() const
{return ( (NULL == leftChild) && (NULL == rightChild)); }

// Check if this has a non-empty left child.
bool SPnode::hasLCwithBox() const
{return ( (leftChild != NULL) && ((leftChild->theBox) != NULL)); }

// Check if this has a non-empty right child.
bool SPnode::hasRCwithBox() const
{return ( (rightChild != NULL) &&
	((rightChild->theBox) != NULL)); }


// find node's split dimension
// note - if I took this out and had spnodes record split dimension I'd have
// to make sure split dimension was adjusted if we marginalised or sliced
int SPnode::getSplitDim() const
{
	int splitDim = -1;	
	
	if ( !isLeaf() ) {
	
		ivector box = getBox();
		int dim = VecLen(box);
		int boxLB = Lb(box);
		
		ivector boxChild;
		
		if (hasRCwithBox()) boxChild = getRightChild()->getBox();
		else if (hasLCwithBox()) boxChild = getLeftChild()->getBox();
		else throw std::logic_error(
			"SPnode::getSplitDim() :Cannot find split dimension");						
		
		int index = 0;
		int boxChildLB = Lb(boxChild);
		while (splitDim < 1 && index < dim) {
			if ((Inf(box[boxLB + index]) 
					!= Inf(boxChild[boxChildLB + index]))
				||
				(Sup(box[boxLB + index]) != 
					Sup(boxChild[boxChildLB + index]))) {
						// found splitDim
						splitDim = boxLB + index;
			}
			index ++;
		} // end while
		
		if (splitDim < 0) {
			throw std::logic_error(
			"SPnode::getSplitDim() :Cannot find split dimension");
		}
	} // end !isLeaf
						
	return splitDim;
}

real SPnode::getSplitValue() const
{
	real splitValue(0.0);
	
	if (!isLeaf()) {
	
		int splitDim = getSplitDim();
		
		if ( hasLCwithBox() ) {
	
			splitValue = Sup(getLeftChild()->getBox()[splitDim]);
		}
		else if ( hasRCwithBox() ) {
	
			splitValue = Inf(getRightChild()->getBox()[splitDim]);
		}
		else throw std::logic_error("SPnode::getSplitValue()");	
	}
	return splitValue;	
}

// get the total node depth in the tree rooted at this
unsigned long int SPnode::getTotalLeafDepth() const
{
	unsigned long int thisDepth = getNodeDepth();
	
	return _getTotalLeafDepth(thisDepth);
}




// get number of leaves
size_t SPnode::getNumberLeaves() const
{
	size_t retVal=0;
	
	if (isLeaf()) retVal = 1; // leaf

	else {

		if (hasLCwithBox()) 
			retVal += getLeftChild()->getNumberLeaves();
		if (hasRCwithBox()) 
			retVal += getRightChild()->getNumberLeaves();

	}

	return retVal;
}

//get number of cherries
size_t SPnode::getNumberCherries() const
{
	size_t retVal = 0;
	
	if (isSubLeaf()) retVal = 1;

	//if children, recurse on the children
	else if (!isLeaf()) {
		if (hasLCwithBox()) 
			retVal += getLeftChild()->getNumberCherries();
		if (hasRCwithBox()) 
			retVal += getRightChild()->getNumberCherries();

	}
	
	return retVal;
}

// return a reference to a container of SPnodes
// contents being the leaves descended from this, or this if this is a leaf
// left to right order
SPnodePtrs& SPnode::getSPnodeLeaves(SPnodePtrs& leaves)
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getSPnodeLeaves(leaves);
	}

	if (hasRCwithBox()) {
		getRightChild()->getSPnodeLeaves(leaves);
	}

	if ( isLeaf() ) { // this is a leaf
		leaves.push_back(this);
	}
	return leaves;
}

// return a reference to a container of const SPnodes
// contents being the leaves descended from this, or this if this is a leaf
// left to right order
SPnodeConstPtrs& SPnode::getConstSPnodeLeaves(SPnodeConstPtrs& leaves) const
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getConstSPnodeLeaves(leaves);
	}

	if (hasRCwithBox()) {
		getRightChild()->getConstSPnodeLeaves(leaves);
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
SPnodePtrs& SPnode::getSPnodeSubLeaves(SPnodePtrs& subleaves)
{
	if (isSubLeaf()) { // this is a subleaf
		subleaves.push_back(this);
	}
	
	//else if children, recurse on the children
	else if (!isLeaf()) {
		getLeftChild()->getSPnodeSubLeaves(subleaves);

		getRightChild()->getSPnodeSubLeaves(subleaves);
	}

	return subleaves;
	
}

// return a reference to a container of SPSnodes
// contents being the sub-leaf children of the given node
// sub-leaf nodes are the parents of leaf nodes and only have leaf nodes
// as children
// left to right order
SPnodeConstPtrs& SPnode::getConstSPnodeSubLeaves(
						SPnodeConstPtrs& subleaves) const
{
	if (isSubLeaf()) { // this is a subleaf
		subleaves.push_back(this);
	}
	//if children, recurse on the children
	else if (!isLeaf()) {
		getLeftChild()->getConstSPnodeSubLeaves(subleaves);

		getRightChild()->getConstSPnodeSubLeaves(subleaves);

	}
	
	return subleaves;
	
}

// return a reference to a container of reals
// contents being the volumes of the leaves descended from this,
// or the volume of this if this is a leaf
// left to right order
RealVec& SPnode::getLeafNodeVolumes(RealVec& vols) const
{
	//if children, recurse on the children
	if (hasLCwithBox()) {
		getLeftChild()->getLeafNodeVolumes(vols);
	}

	if (hasRCwithBox()) {
		getRightChild()->getLeafNodeVolumes(vols);
	}

	if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
		vols.push_back(nodeVolume());
	}
	return vols;
}

// fills in the leaf node levels, left to right
IntVec& SPnode::getLeafNodeLevels(IntVec& levels, int level) const
{
	
	if (getLeftChild()!=NULL) {
		getLeftChild()->getLeafNodeLevels(levels, level+1);
	}
	if (getRightChild()!=NULL) {
		getRightChild()->getLeafNodeLevels(levels, level+1);
	}
	if (getLeftChild()==NULL && getRightChild()==NULL) {

		levels.push_back(level);
	}
	return levels;
}


// returns a string of the leaf levels
// left to right, 0 is root
std::string SPnode::getLeafNodeLevelsString() const
{
	IntVec levels;
	getLeafNodeLevels(levels); // fill container of leaf levels
	
	bool compact = true;
	return toString(levels, compact);
}


// Get the node depth from the root
int SPnode::getNodeDepth() const
{
	int count = 0;
	if (getParent() != NULL) {
		count ++;
		count += getParent()->getNodeDepth();
	}
	
	return count;
}

//Returns the height of the tree rooted at this node
int SPnode::getTreeHeight() const
{
	int myDep = getNodeDepth();
	int maxChildDepth = myDep;
	
	if ( !isLeaf() ) {
		IntVec levels;
		levels = getLeafNodeLevels(levels, myDep);
		
		maxChildDepth = *( std::max_element (levels.begin(), levels.end()) );
	}
	return maxChildDepth - myDep;
}


//Returns the volume of the smallest (by vol) leaf node.
double SPnode::getSmallestLeafVol() const
{
	if (isLeaf() ) {
		return nodeVolume();
	}
	else {
		
		if ( hasLCwithBox() && hasRCwithBox() ) {
			double myL = getLeftChild()->getSmallestLeafVol();
			double myR = getRightChild()->getSmallestLeafVol();
		
			return (myL < myR ? myL : myR);
		}
		else if ( hasLCwithBox() ) {
			return getLeftChild()->getSmallestLeafVol();
		}
		else if ( hasRCwithBox() ) {
			return getRightChild()->getSmallestLeafVol();
		}
		else throw std::runtime_error("SPnode::getSmallestLeafVol()");
	}
	
}

// Returns the volume of the largest (by vol) leaf node.
double SPnode::getLargestLeafVol() const
{
	if (isLeaf() ) {
		return nodeVolume();
	}
	else {
		double myL = 0.0;
		double myR = 0.0;
		
		if ( hasLCwithBox() ) myL = getLeftChild()->getLargestLeafVol();
		if ( hasRCwithBox() ) myR = getRightChild()->getLargestLeafVol();
		
		return (myL > myR ? myL : myR);
	}
}

bool SPnode::checkTreeStateLegal() const
{
	// check current state is legal by looking at everything not a leaf
	bool legal = true;
	if ( !isLeaf() ) {
		legal = isSplittableNode();
		if (legal && hasLCwithBox() ) {
				legal = getLeftChild()->checkTreeStateLegal();
		}
		if (legal && hasRCwithBox() ) {
				legal = getRightChild()->checkTreeStateLegal();
		}
	}
	
	return legal;
}

 // used by nodes to find if this has a leaf sibling
bool SPnode::hasLeafSibling() const
{
	bool retValue = false;
	if (parent != NULL && getParent()->hasRCwithBox() && getParent()->hasLCwithBox()) {
		if (getParent()->getRightChild() == this) {
			retValue = getParent()->getLeftChild()->isLeaf();
		}
		else {
			retValue = getParent()->getRightChild()->isLeaf();
		}
	}
	return retValue;
}

cxsc::real SPnode::getRootVolume() const

{
	if (parent == NULL) return nodeRealVolume();
	else return getParent()->getRootVolume();
}


void SPnode::acceptSPnodeVisitor(const SPnodeVisitor& visitor)
{
	visitor.visit(this);
	
	if (hasLCwithBox()) (getLeftChild())->acceptSPnodeVisitor(visitor);
	if (hasRCwithBox()) (getRightChild())->acceptSPnodeVisitor(visitor);

}

void SPnode::acceptSPCheckVisitor(const SPCheckVisitor& visitor) const
{
	visitor.visit(this);
	
}

// Print the box of a specific node in a subpaving
std::ostream& SPnode::nodePrint(std::ostream &os) const
{
// output in form: box

	if( !isEmpty() ) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox
		
		for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
			os << "[ " << Inf(thisBox[i]) << " , "
				<< Sup(thisBox[i]) << " ]";
			if (i < Ub(thisBox) ) {
				os << " , ";
			}
		}
	}
	return os;

}

// Print the box of a specific node in a subpaving
std::ostream& SPnode::nodePrint(std::ostream &os, int prec) const
{

	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	nodePrint(os);
	os << cxsc::RestoreOpt;
	
	return os;
}


//Output for all the  leaf boxes in this, using tab delimiters
std::ostream& SPnode::leavesOutputTabs(std::ostream &os) const
{
	// uses  member function leafOutputTabs to generate node output
	if ( !(isEmpty()) && isLeaf() ) { // this is a non-empty leaf
		leafOutputTabs(os);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabs(os);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabs(os);
	}
    return os;

}

//Output for all the  leaf boxes in this, using tab delimiters
std::ostream& SPnode::leavesOutputTabs(std::ostream &os, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	leavesOutputTabs(os);
	os << cxsc::RestoreOpt;
	
	return os;

}

//Output for all the nodes in this, using tab delimiters
// uses LeafOutputTabs, which is a virtual function
// and can be redefined in derived classes
std::ostream& SPnode::nodesAllOutput(std::ostream &os,
									int level) const
{

	if ( !isEmpty() ) {
		leafOutputTabs(os);
		os << "\n";
	}
	if(getLeftChild()!=NULL) {
		for (int i = 0; i < level; i++) {
			os << "\t";
		}
		getLeftChild()->nodesAllOutput(os, (level+1));
	}
	
	if(getRightChild()!=NULL) {
		for (int i = 0; i < level; i++) {
			os << "\t";
		}
		getRightChild()->nodesAllOutput(os, (level+1));
	}
	
	return os;
}

//Output for all the nodes in this, using tab delimiters
// uses LeafOutputTabs, which is a virtual function
// and can be redefined in derived classes
std::ostream& SPnode::nodesAllOutput(std::ostream &os,
									int level, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	nodesAllOutput(os, level);
	os << cxsc::RestoreOpt;
	
	return os;

}

// Method to add current state to a log file
// Output goes to file named according to argument s
void SPnode::outputLog(const std::string& s, int i, int prec) const
{
    // To add output to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        leavesOutputTabs(os, prec); // the output
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

// Method to add current state to a log file
// Output goes to file named according to argument s
void SPnode::outputLog(const std::string& s, int i) const
{
    // To add output to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        leavesOutputTabs(os); // the output
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

// make a .dot file for the subpaving
bool SPnode::outputGraphDot() const
{
	bool success = false;

	std::string toParse = getChildNodeNames();

	std::string baseFileName = "graph";
	std::string suffix = ".dot";
	std::string s = getUniqueFilename(baseFileName, suffix);
	outputFile(s, "digraph G {"); // opening line

	if (toParse.length() > 0) {

		success = parseForGraphDot(s, toParse);

	}
	else if (nodeName.length() > 0) {
		// only the root node
		std::string line = "\t " + nodeName + ";";
		outputFile(s, line);
		success = true;

   }

	outputFile(s, "}"); // closing line

	// make the image from the dot file
	makeDotImage(s);

	return success;
}

// Check if a node contains a datapoint
// it is assumed that the node will have a box
// childInd is an indicator for which child is being checked
bool SPnode::nodeContains(const rvector& p,
						OPERATIONS_ON childInd) const
{
	bool retValue = false;
	
	// only check for total containment if this is a parent node
	// or to be treated as such
	if ((parent == NULL) || (childInd == ON_PARENT)) {
	
		// cast p to an ivector
		ivector pvec = _ivector(p);
	
		//find if p is in the box
		if (pvec <= getBox()) {
			retValue = true;
		}
	}
	// if not to be treated as a whole box, we assume it was in the parent
	// and only check ourselves with respect to the split dimension
	// and if this is a right child it can be anywhere 
	// but  this is a left child,
	// we need to check the split
	// find what dimension the parent was split on and what
	// the split value was
	// pvector must be strictly less than parentSplitValue
	// on the split dimension
	else if (parent != NULL) { // truly not a parent node
	
		int parentSplitDim = (getParent())->getSplitDim();
		
		if (childInd == ON_LEFT) { 
			retValue = (p[parentSplitDim] < Sup(getBox()[parentSplitDim]));
		}
		if (childInd == ON_RIGHT) { 
			retValue = !(p[parentSplitDim] < Inf(getBox()[parentSplitDim]));
		}
	}
	
	return retValue;
}

// add two sibling nodes to this provided that this node is a leaf
// comp argument is passed to Upper() and Lower()
// split the box in half normal to dimension set by comp
void SPnode::nodeExpand(int comp)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::nodeExpand(int )");
	}
	
	// only do something if this SPnode is a leaf
	if (isLeaf()) {
		
		SPnode* newLC = NULL;
		SPnode* newRC = NULL;
		
		try {

			// ivectors to be new boxes for new children
			ivector lC, rC;

			// Call Lower() and Upper() to put split boxes
			// into lC and rC respectively
			Lower(getBox(), lC, comp);
			Upper(getBox(), rC, comp);
			
			newLC = new SPnode(lC);
			newRC = new SPnode(rC);
			
			nodeAddLeft(newLC);
			nodeAddRight(newRC);
			
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

// finds the dimension to split on and then calls nodeExpand(int comp)
void SPnode::nodeExpand()
{
	int maxdiamcomp; // variable to hold first longest dimension
	MaxDiam(getBox(), maxdiamcomp);
	nodeExpand(maxdiamcomp); // complete nodeExpand

}

// reabsorb children into this
// needs to recursively delete the children
void SPnode::nodeReabsorbChildren()
{
	// first recursively deal with the children of the children
	if (hasLCwithBox())
		getLeftChild()->nodeReabsorbChildren();
	if (hasRCwithBox())
		getRightChild()->nodeReabsorbChildren();

	// now deal with itself
	if (leftChild != NULL) delete getLeftChild();
	leftChild = NULL;
	if (rightChild != NULL) delete getRightChild();
	rightChild = NULL;

}


// add lChild onto this node
void SPnode::nodeAddLeft(SPnode *lChild)
{
	if (lChild != NULL) {
		if ( isEmpty() ) {
			throw NoBox_Error("SPnode::nodeAddLeft(SPnode *)");
		}
		if ( lChild->isEmpty() ) {
			throw NoBox_Error("SPnode::nodeAddLeft(SPnode *)");
		}
		if ( lChild->getDimension() != getDimension() ) {
			throw IncompatibleDimensions_Error("SPnode::nodeAddLeft(SPnode *)");
		}

		// check no existing left child
		if (hasLCwithBox()) {
			throw std::logic_error(
				"SPnode::nodeAddLeft(SPnode *)");
		}
		// check box to be added
		if (*(lChild->theBox) == *theBox) {

			throw UnfulfillableRequest_Error(
				"SPnode::nodeAddLeft(SPnode *) : Child box same as this's box.");
		}
		if (hasRCwithBox()) {
			// check that this new box fits with the current box
			if (*theBox !=
				(*(lChild->theBox) | *(getRightChild()->theBox))) {
					throw IncompatibleDimensions_Error(
										"SPnode::nodeAddLeft(SPnode *)");
			}

		}

		this->leftChild = lChild;
		leftChild->parent = this;
	}
}

// add Child onto this node
void SPnode::nodeAddRight(SPnode *rChild)
{
	if (rChild != NULL) {
		if ( isEmpty() ) {
			throw NoBox_Error("SPnode::nodeAddRight(SPnode *)");
		}
		if ( rChild->isEmpty() ) {
			throw NoBox_Error("SPnode::nodeAddRight(SPnode *)");
		}
		if ( rChild->getDimension() != getDimension() ) {
			throw IncompatibleDimensions_Error("SPnode::nodeAddRight(SPnode *)");
		}

		// check no existing right child
		if (hasRCwithBox()) {
			throw std::logic_error(
				"SPnode::nodeAddRight(SPnode *)");
		}
		
		// check box to be added
		if (*(rChild->theBox) == *theBox) {

			throw UnfulfillableRequest_Error(
				"SPnode::nodeAddRight(SPnode *) : Child box same as this's box");
		}

		if (hasLCwithBox()) {
			// check that this new box fits with the current box
			if (*theBox !=
				(*(rChild->theBox) | *(getLeftChild()->theBox))) {
				throw IncompatibleDimensions_Error(
										"SPnode::nodeAddRight(SPnode *)");
			}
		}

		this->rightChild = rChild;
		rightChild->parent = this;
	}
}


// Makes a string of child names left to right order
std::string SPnode::getChildNodeNames() const
{
	std::string retStr = "";
	SPnode * LChild =  getLeftChild();
	SPnode * RChild =  getRightChild();

	if (LChild != NULL) {
		retStr += " ";
		retStr += LChild->getNodeName();
	}
	if (RChild != NULL) {
		retStr += " ";
		retStr += RChild->getNodeName();
	}
	if (LChild != NULL) {
		retStr += LChild->getChildNodeNames();
	}
	if (RChild != NULL) {
		retStr += RChild->getChildNodeNames();
	}
	return retStr;
}

/* \brief Reshape so that the tree rooted at this has shape that
is the union of this shape and the shape of another tree.

Throws a NoBox_Error if this has no box or if \a other has no box. 
Throws an IncompatibleDimensions_Error if boxes of this and \a other
are not the same.

\param other is the tree to make the union against.
\pre This has a box and that box is identical to the box of \a other. 
\post the tree rooted at this has shape that is the
union of the shape of this before the operation and the shape of 
\a other.  \a other is unchanged.      */
void SPnode::reshapeToUnion(const SPnode& other)
{
	if (isEmpty() || other.isEmpty()) {
		throw NoBox_Error(
		"SPnode::reshapeToUnion(const SPnode&)");
	}
	if ( getBox() != other.getBox() )  {
		throw IncompatibleDimensions_Error(
		"SPnode::reshapeToUnion(const SPnode&)");
	}
	if ( !other.checkTreeStateLegal() )
		throw runtime_error(
		"SPnode::reshapeToUnion(const SPnode&) : other has illegal tree state");
	
	this->_reshapeToUnion(&other);
	
}

bool SPnode::randomSplitRootAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						bool saveInstructions)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::randomSplitRootAtLeast(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::randomSplitRootAtLeast(...)");
	}
		
	if(saveInstructions) partitioner.initialiseInstructions(numLeaves);
	
	
	bool success = _randomSplitAtLeast(numLeaves,
						partitioner);
	
	// clear the instructions if not successful?
	
	#ifdef DEBUG_MCMC_SPLIT
			cout << "randomSplitRootAtLeast" << endl;
			if (success) cout << "returning sucessful"  << endl;
			else cout << "returning unsucessfully"  << endl;
	#endif
	return success;
	
}

bool SPnode::randomNaturalSplitRootAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						bool saveInstructions)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::randomNaturalSplitRootAtLeast(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::randomNaturalSplitRootAtLeast(...)");
	}
		
	if(saveInstructions) partitioner.initialiseInstructions(numLeaves);
	
	
	bool success = _randomNaturalSplitAtLeast(numLeaves,
						partitioner);
	
	// clear the instructions if not successful?
	
	#ifdef DEBUG_MCMC_SPLIT
			cout << "randomNaturalSplitRootAtLeast" << endl;
			if (success) cout << "returning sucessful"  << endl;
			else cout << "returning unsucessfully"  << endl;
	#endif
	return success;
	
}


// split a root box to a shape specified by the instruction string
bool SPnode::splitRootToShape(std::string instruction)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::splitRootToShape(std::string)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::splitRootToShape(std::string)");
	}
	
	bool success = false;
	std::string res = splitLeft(instruction);
	if (res == "") success = true;
	
	return success;
	
}

bool SPnode::splitRootAtLeastToShape(std::vector < size_t > reqDepths)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::splitRootAtLeastToShape(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::splitRootAtLeastToShape(...)");
	}
	
	size_t myDepth = 0;
	
	bool success = _splitAtLeastToShape(reqDepths,
								myDepth);
	
	if (success && !reqDepths.empty()) 
		throw std::invalid_argument("SPnode::splitRootAtLeastToShape(...)");
	
	return success;
	
}

/* assumes a binary tree type instruction, ie any non-root node has a sibling
 * will expand or prune as required but if it expands, will use nodeExpand for this 
 * type of node.
 * Make bool return type to allow subtypes to override if necessary */
bool SPnode::reshapeToMirrorValidSplitPartitions(
						const std::vector < unsigned long int >& partitions,
						unsigned long int numLeaves)
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::reshapeToMirrorValidSplitPartitions(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::reshapeToMirrorValidSplitPartitions(...)");
	}
	if ( (numLeaves > 2) && partitions.empty() ) {
		throw std::invalid_argument(
			"SPnode::reshapeToMirrorValidSplitPartitions(...) : instructions empty()");
	}
	if (!numLeaves) {
		throw std::invalid_argument(
			"SPnode::reshapeToMirrorValidSplitPartitions(...) : numLeaves = 0");
	}
	
	size_t currentIndex = 0;
	
	currentIndex = _reshapeToMirrorValidSplitPartitions(
						partitions,
						numLeaves,
						currentIndex);
	
	assert(currentIndex >= partitions.size());
    
    return (currentIndex >= partitions.size());
}

// return total number of leaves
unsigned long int SPnode::getPartitions(
				std::vector < unsigned long int >& partitions) const
{
	if ( isEmpty() ) {
		throw NoBox_Error("SPnode::getMyPartitions(...)");
	}
	if (getParent() != NULL) {
		throw NonRootNode_Error("SPnode::getMyPartitions(...)");
	}
	
	std::vector < unsigned long int >tmp;
	partitions = tmp; // clears partitions just in case
	
	tmp.reserve(1000000); // arbitrary large amount
		
	unsigned long int numLeaves = _getMyPartitions(tmp);
	
	#ifdef DEBUG_PARTITIONS
		cout << "\ntmp size on return is = " << tmp.size() << endl; 
	#endif
	partitions.resize(tmp.size());
	
	reverse_copy (tmp.begin(), tmp.end(), partitions.begin());

	#ifdef DEBUG_PARTITIONS
		cout << "\nafter reverse copy, partitions size is = " << partitions.size() << endl; 
	#endif

	return numLeaves;
}

void SPnode::swap (SPnode& spn) //throw()
{
	std::swap(theBox, spn.theBox); // theBox is a pointer
	std::swap(nodeName, spn.nodeName);
	std::swap(parent, spn.parent);
	// can just swap child pointers
	std::swap(leftChild, spn.leftChild);
	std::swap(rightChild, spn.rightChild);
	// children have to be repointed to the new swapped parents
	if (leftChild != NULL) leftChild->parent = this;
	if (rightChild != NULL) rightChild->parent = this;
	if (spn.leftChild != NULL) spn.leftChild->parent = &spn;
	if (spn.rightChild != NULL) spn.rightChild->parent = &spn;
	
}

std::string SPnode::nodeStringSummary() const
{
	std::ostringstream oss;
	
	oss << "I am " << getNodeName() << "(address " << this << "),\n";
	oss << "Dimension is " << getDimension() << ", address of box is " << theBox << "\n"; 
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


// -----------------protected methods



// recursively split children according to instruction string
std::string SPnode::splitLeft(std::string instruction)
{
	// parse string to get first level out
	std::string comma = ", ";
	size_t startpos = instruction.find_first_not_of(comma);
	size_t endpos = std::string::npos;
	size_t newstartpos = std::string::npos;
	int depth = -1;
	std::string str = "";
	if (startpos != std::string::npos) {
		endpos = instruction.find_first_of(comma, startpos);
		//not the last digit
		if (endpos != std::string::npos) {
			str = instruction.substr(startpos, endpos-startpos);
			newstartpos = instruction.find_first_not_of(comma,
													endpos + 1);
		}
		//last digit
		else {
			str = instruction.substr(startpos);
		}
	}
	if (str.empty()) {
		throw std::invalid_argument(
				"SPnode::splitLeft(std::string): missing instruction depth");
	}
	depth = atoi(str.c_str()); // 0 if not valid integer
	
	int myDepth = getNodeDepth();
	if (myDepth < depth) { // need to go down more

		if (isLeaf()) {
			// split, using the nodeExpand for this subtype if not base
			nodeExpand();
		}

		// send splitLeft instruction down to the left
		instruction = getLeftChild()->splitLeft(instruction);
		instruction = getRightChild()->splitLeft(instruction);
		
	}
	else if (myDepth == depth) {     // split enough
		if (!isLeaf()) {
			throw std::logic_error(
					"SPnode::splitLeft(std::string): not a leaf here");
		}
		else if (newstartpos != std::string::npos) {
			instruction = instruction.substr(newstartpos);
		}
		else instruction = "";
	}
	else { // myDepth is > depth
		throw std::logic_error(
			"SPnode::splitLeft(std::string): node depth > instruction depth");
	}

	return instruction;
}

// recursively split children according to instruction
/* reqDepths is our instruction, index is where we should look
 * ensures this has at least shape of instruction, but 
 * this can also be more split at some or all points.
 * Splits from the right*/
bool SPnode::_splitAtLeastToShape(
						std::vector < size_t >& reqDepths,
						size_t myDepth)
{
	if(reqDepths.empty()) 
		throw std::invalid_argument(
			"SPnode::_splitRight(...): : reqDepths.empty()");
		
	bool success = true;
	
	size_t depth = reqDepths.back();
	
	if (myDepth < depth) { // need to try to go down more

		if ( isLeaf() ) {
			
			if (isSplittableNode()) {
				// split, using the nodeExpand for this subtype if not base
				nodeExpand();
			}
			else success = false;
		}

		// if we are okay, send instruction down RIGHT SIDE FIRST
		if (success) success = getRightChild()->_splitAtLeastToShape(reqDepths,
						myDepth+1);
		if (success) success = getLeftChild()->_splitAtLeastToShape(reqDepths,
						myDepth+1);
		
	}
	else if (myDepth == depth) {     // split enough

		/* I am not necessarily a leaf though - allow for situation
		 * where I want to be split to at least this much */
	
		// knock the last element out
		reqDepths.erase(reqDepths.begin()+reqDepths.size()-1);
	}
	else { // myDepth is > depth
		throw std::invalid_argument(
			"SPnode::_splitRight(...): node depth > instruction depth");
	}

	return success;
}

/* assumes a binary tree type instruction, ie any non-root node has a sibling
 * will expand or prune as required but if it expands, will use nodeExpand for this 
 * type of node */
size_t SPnode::_reshapeToMirrorValidSplitPartitions(
						const std::vector < unsigned long int >& partitions,
						unsigned long int numLeaves,
						size_t currentIndex)
{
	#ifdef DEBUG_RESHAPE
        size_t n = partitions.size();
		cout << "\nIn _reshapeToMirrorValidSplitPartitions. I am " << nodeName << endl;
		cout << "partitions.size() = " << n << " and numLeaves = " << numLeaves << " and currentIndex = " << currentIndex << endl;
	#endif

	if(numLeaves == 1) {
		// chop off children
		nodeReabsorbChildren();
		#ifdef DEBUG_RESHAPE
			cout << "chopped off children: isLeaf() = " << (isLeaf()) << endl;
		#endif
		// increment current index
		//currentIndex++;
		#ifdef DEBUG_RESHAPE
			cout << "chopped off children: isLeaf() = " << (isLeaf()) << endl;
			cout << "currentIndex = " << currentIndex << endl;
		#endif
	}
	else if(numLeaves == 2) {
		if(isLeaf()) nodeExpand();
		else {
			// chop off children's children
			getLeftChild()->nodeReabsorbChildren();
			getRightChild()->nodeReabsorbChildren();
			#ifdef DEBUG_RESHAPE
				cout << "chopped off children's children" << endl;
				cout << "getLeftChild()->isLeaf() = " << (getLeftChild()->isLeaf()) << endl;
				cout << "getRightChild()->isLeaf() = " << (getRightChild()->isLeaf()) << endl;
			#endif
		}
		//++currentIndex;
		#ifdef DEBUG_RESHAPE
			cout << "currentIndex = " << currentIndex << endl;
		#endif
	}
		
	else {
		assert(currentIndex < partitions.size());
		unsigned long int thisInstruction = partitions[currentIndex];
		#ifdef DEBUG_RESHAPE
			cout << "partitions[currentIndex] = " << thisInstruction << endl;
		#endif
		
        if(isLeaf()) nodeExpand();
        
        currentIndex = getLeftChild()->_reshapeToMirrorValidSplitPartitions(
                    partitions, thisInstruction, currentIndex+1);
        #ifdef DEBUG_RESHAPE
            cout << "\nafter going Left, back in _reshapeToMirrorValidSplitPartitions for " << nodeName << endl;
            cout << "currentIndex = " << currentIndex << endl;
        #endif
        currentIndex = getRightChild()->_reshapeToMirrorValidSplitPartitions(
                    partitions, numLeaves - thisInstruction, currentIndex);
        #ifdef DEBUG_RESHAPE
            cout << "\nafter going right, back in _reshapeToMirrorValidSplitPartitions for " << nodeName << endl;
            cout << "currentIndex = " << currentIndex << endl;
        #endif

	}
	return currentIndex;
	
}

/* not sure how much sense this will make if not a binary tree */
unsigned long int SPnode::_getMyPartitions(
			std::vector < unsigned long int >& partitions) const
{
	#ifdef DEBUG_PARTITIONS
		cout << "\nI am " << nodeName << endl;
		cout << "Partitions.size = " << partitions.size() << endl;
	#endif
	unsigned long int leaves = 0;
	
	//if I am a leaf do nothing
	if (isLeaf()) {
		leaves = 1;

	}
	else if (isSubLeaf()) { // I have two leaf children
		
		leaves = 2;
		// no need to record anything
	
	}
	else {
		
		#ifdef DEBUG_PARTITIONS
			cout << "recursing" << endl;
		#endif 
		
		if (NULL != rightChild) {
			unsigned long int rightLeaves 
					= getRightChild()->_getMyPartitions(partitions);
			
			leaves += rightLeaves;
			#ifdef DEBUG_PARTITIONS
				cout << "\nI am " << nodeName << ", after going right, added " << rightLeaves << " leaves now is " << leaves << endl;
			#endif 
		}
		
		if (NULL != leftChild) {
			unsigned long int leftLeaves 
					= getLeftChild()->_getMyPartitions(partitions);
			
			leaves += leftLeaves;
			partitions.push_back(leftLeaves);
			#ifdef DEBUG_PARTITIONS
				cout << "\nI am " << nodeName << ", after going left, added " << leftLeaves << " leaves now is " << leaves << endl;
			#endif 
		}
		assert (leaves > 2);
		
			
	}
	#ifdef DEBUG_PARTITIONS
		cout << "\n I am " << nodeName << " returning leaves = " << leaves << endl;
		cout << "Partitions.size = " << partitions.size() << endl;
	#endif 
	return leaves;
	
}



// reshape this tree to have union of this shape and shape of other
// no checks on boxes since this should be redundant if used by unionNoData...
//
void SPnode::_reshapeToUnion(const SPnode * const other)
{
		if ( other != NULL && !(other->isEmpty()) ) {

		// this is not a leaf, other is a leaf
		if (!isLeaf() && other->isLeaf()) {

			// no need to do anything
		}

		// this is a leaf, other is not a leaf
		if (isLeaf() && !other->isLeaf()) {

			/* we need to expand this, */
			nodeExpand();
			
		}

		// now recurse on the children if both have children
		if (!isLeaf() && !other->isLeaf()) {
			getLeftChild()->_reshapeToUnion(other->getLeftChild());
			getRightChild()->_reshapeToUnion(other->getRightChild());
		}
	}
}

/* random partitioning of this and children until tree rooted at 
this has at least numleaves or a node cannot be split.  This can
* already be more split than the random partitioner requires in which
* case this will be unchanged*/
bool SPnode::_randomSplitAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner)
{
	#ifdef DEBUG_MCMC_SPLIT
		cout << "_randomSplitAtLeast, I am " << nodeName << endl;
		cout << "numLeaves = " << numLeaves << endl;
	
	#endif
	
	bool success = true;
	/* if (numLeaves == 1)  nothing to be done  */
	
	if (numLeaves == 2) {
		if (isSplittableNode()) {
			
			nodeExpand();
				
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
	else if (numLeaves > 2) {
		
		unsigned long int left = 
					partitioner.generateStatePartition(numLeaves);
		
		#ifdef DEBUG_MCMC_SPLIT
			if (!isLeaf()) cout << "left = " << left << endl;
		#endif
		
		if (isSplittableNode()) {
			
			nodeExpand();
				
			if (!isLeaf()) {
				
				#ifdef DEBUG_MCMC_SPLIT
					cout << "recursing\n" << endl;
				#endif
				
				// left side
				success = getLeftChild()->_randomSplitAtLeast(
						left,
						partitioner);
				// right side
				if (success) success = getRightChild()->_randomSplitAtLeast(
						numLeaves - left,
						partitioner);
			}
			else throw std::logic_error("SPnode::_randomSplitAtLeast(...)");
		}
		else {
			success = false;
			#if defined (DEBUG_MCMC_SPLIT) || defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "could not expand " << nodeName << endl;
			#endif
			#if defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "numLeaves = " << numLeaves << endl;
			#endif
		}
	}
	
	
	return success;
}

/* random partitioning of this and children until tree rooted at 
this has at least numleaves or a node cannot be split.  This can
* already be more split than the random partitioner requires in which
* case this will be unchanged*/
bool SPnode::_randomNaturalSplitAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner)
{
	#ifdef DEBUG_MCMC_SPLIT
		cout << "_randomNaturalSplitAtLeast, I am " << nodeName << endl;
		cout << "numLeaves = " << numLeaves << endl;
	
	#endif
	
	bool success = true;
	/* if (numLeaves == 1)  nothing to be done  */
	
	if (numLeaves == 2) {
		if (isSplittableNode()) {
			
			nodeExpand();
				
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
	else if (numLeaves > 2) {
		
		unsigned long int left = 
					partitioner.generateNaturalStatePartition(numLeaves);
		
		#ifdef DEBUG_MCMC_SPLIT
			if (!isLeaf()) cout << "left = " << left << endl;
		#endif
		
		if (isSplittableNode()) {
			
			nodeExpand();
				
			if (!isLeaf()) {
				
				#ifdef DEBUG_MCMC_SPLIT
					cout << "recursing\n" << endl;
				#endif
				
				// left side
				success = getLeftChild()->_randomNaturalSplitAtLeast(
						left,
						partitioner);
				// right side
				if (success) success = getRightChild()->_randomNaturalSplitAtLeast(
						numLeaves - left,
						partitioner);
			}
			else throw std::logic_error("SPnode::_randomNaturalSplitAtLeast(...)");
		}
		else {
			success = false;
			#if defined (DEBUG_MCMC_SPLIT) || defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "could not expand " << nodeName << endl;
			#endif
			#if defined (DEBUG_MCMC_SPLIT_FAIL)
				cout << "numLeaves = " << numLeaves << endl;
			#endif
		}
	}
	
	
	return success;
}

// get the total node depth in the tree rooted at this
unsigned long int SPnode::_getTotalLeafDepth(
		unsigned long int thisDepth) const
{
	unsigned long int retVal = 0;
	
	if (isLeaf()) retVal = thisDepth;

	else {

		if (hasLCwithBox()) 
			retVal += getLeftChild()->_getTotalLeafDepth(thisDepth+1);
		if (hasRCwithBox()) 
			retVal += getRightChild()->_getTotalLeafDepth(thisDepth+1);

	}

	return retVal;
}

// Print the details of a single leaf node, using tab delimiters
std::ostream& SPnode::leafOutputTabs(std::ostream &os) const
{

	if( !isEmpty() ) { // do nothing if there is no box

		ivector thisBox = *theBox; // copy theBox
		
		// output the name, nodeVolume
		os << nodeName;
		os << "\t" << nodeRealVolume();

		// followed by intervals making up box using Inf & Sup
		// ie unlike cxsc output, there is no [  ] around them
		for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

			os << "\t" << Inf(thisBox[i]) << "\t"
				<< Sup(thisBox[i]);
		}
		
	}
    return os;

}

// Print the details of a single leaf node, using tab delimiters
std::ostream& SPnode::leafOutputTabs(std::ostream &os, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	leafOutputTabs(os);
		
	os << cxsc::RestoreOpt;
    return os;

}

// ------------------- private methods ------------------------------

// ensure theBox and children (if any) are deleted if constructed in failed constructor
void SPnode::constructor_error_handler() 
{
	try {
			delete theBox;
			theBox = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	try {
			delete leftChild;
			leftChild = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	try {
			delete rightChild;
			rightChild = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}




// ------------------- end of SPnode class definitions ---------------




// -------------------- start of SPnode non-member functions



//Output the leaf boxes in SubPaving spn
std::ostream & subpavings::operator<<(std::ostream &os, const SPnode& spn)
{
	// uses nodePrint to generate node output

	if (!(spn.isEmpty()) && spn.isLeaf()) { // spn is non-empty leaf
		spn.nodePrint(os);
	}

	//recurse on the children
	if (spn.hasLCwithBox()) {
		os << (*(spn.getLeftChild()));
	}
	if (spn.hasRCwithBox()) {
		os << (*(spn.getRightChild()));
	}

	//in the case where spn is empty we just return os

	return os;
}

// check if a SubPaving is a leaf
bool subpavings::isLeaf(const SPnode * const spn)
{
	// FALSE if spn is a null pointer, true if spn is not NULL
	bool retVal = (spn!=NULL);


	if (retVal) { // if spn points to a non-empty node
		retVal = spn->isLeaf();
	}

	return retVal;
}

// check if a SubPaving is empty
bool subpavings::isEmpty(const SPnode *const spn)
{
	// return true if spn is a null pointer or
	// node spn points to is empty
	return ((spn==NULL) || (spn->isEmpty()));
}

// get volume of a SubPaving
double subpavings::spVolume(const SPnode * const spn)
{
	double retVal =0.0;

	if (!(isEmpty(spn)) && isLeaf(spn)) {

		retVal = spn->nodeVolume();
	}

	// recurse on children
	if (!(isEmpty(spn)) && !(isLeaf(spn))) {

		retVal += spVolume(spn->getLeftChild());
		retVal += spVolume(spn->getRightChild());
	}

	// case of isEmpty(spn) retValue = 0 by default

	return retVal;
}

// get number of leaves of a SubPaving
size_t subpavings::spLeaves(const SPnode * const spn)
{
	size_t retVal=0;

	if (!(isEmpty(spn))) retVal = spn->getNumberLeaves();
	
	return retVal;
}



// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::SPnode & s1, 
			subpavings::SPnode & s2) // throw ()
{
	s1.swap(s2);
}

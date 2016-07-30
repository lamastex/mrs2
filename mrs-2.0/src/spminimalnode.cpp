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
\brief SPMinimalnode (minimal SubPaving) 
and associated non-member functions definitions.
*/

#include "spminimalnode.hpp"

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
#include "sptools.hpp"
#include "voxel_tools.hpp"
#include "sptemplates.hpp"
#include "subpaving_exception.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <stdexcept>



using namespace std;
using namespace subpavings;

// ----------------------- SPMinimalnode class definitions ------------------


// ---------------- public methods

// default constructor
SPMinimalnode::SPMinimalnode() 
	//invokes the base class default constructor
	{}


// initialised constructor, initialised with one ivector for the box
SPMinimalnode::SPMinimalnode(ivector& v) : SPnode(v)
{ }

// initialised constructor, initialised with a LabBox (labeled box)
SPMinimalnode::SPMinimalnode(LabBox& lb) : SPnode(lb)
{ }


//Copy constructor
//copies from the given node downwards
SPMinimalnode::SPMinimalnode(const SPMinimalnode& other) :
	SPnode(*(other.theBox))
{
		nodeName = other.nodeName;
		
		//recursion on the children
		if (other.leftChild) {
			nodeAddLeft(new SPMinimalnode(*(other.getLeftChild())));
		}
		else leftChild=NULL;

		if (other.rightChild) {
			nodeAddRight(new SPMinimalnode(*(other.getRightChild())));
		}
		else rightChild=NULL;

}


// Destructor.
SPMinimalnode::~SPMinimalnode()
{//will call base class destructor
}


//copy assignment operator
//copies from this node downwards
SPMinimalnode& SPMinimalnode::operator=(SPMinimalnode rhs)
{
	rhs.swapMin(*this); // make sure we use our version of swap
	return(*this);
	
}

// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
SPMinimalnode* SPMinimalnode::getParent() const
{ return (SPMinimalnode*) parent; }

// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node.
SPMinimalnode* SPMinimalnode::getLeftChild() const
{ return (SPMinimalnode*) leftChild; }

// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node.
SPMinimalnode* SPMinimalnode::getRightChild() const
{ return (SPMinimalnode*) rightChild; }

// computes a minimal subpaving from two sibling subpavings
// a subpaving is minimal if it has no sibling leaves
// a minimal subpaving is created by discarding sibling leaves
// lChild and rChild are the two subpavings to be reunited
void SPMinimalnode::nodeReunite(SPMinimalnode *lChild, 
								SPMinimalnode *rChild)
{
	// *this is the node which will become the parent

	// if both subpavings are leaves and hull of boxes is x,
	// discard them: *this is a leaf
	if (lChild->isLeaf() && rChild->isLeaf()) {
		if (*theBox !=
			(*(lChild->theBox) | *(rChild->theBox))) {
			throw subpavings::IncompatibleDimensions_Error(
						"nodeReunite(SPMinimalnode*, SPMinimalnode*)");
		}

		//discard the two subpavings given
		delete lChild;
		lChild = NULL;
		delete rChild;
		rChild = NULL;
	}

	else {  // at least one of the children is not a leaf
		// this has to adopt them rather than reuniting them
		nodeAdoptLeft(lChild);
		nodeAdoptRight(rChild);
		recursiveRename(); // recursively rename child branches
	}
}

// graft lChild onto this node
// lChild could be a leaf or a non-leaf
// takes care of the data associated with lChild/its descendents
// used when we are building a subpaving from the leaf nodes upwards
void SPMinimalnode::nodeAdoptLeft(SPMinimalnode *lChild)
{
	// *this is the node which will become the parent

	// point parent and child pointers in the right directions
	// nodeAddLeft() checks labels, hull size, present children
	nodeAddLeft(lChild);
}

// graft lChild onto this node
// lChild could be a leaf or a non-leaf
// takes care of the data associated with lChild/its descendents
// used when we are building a subpaving from the leaf nodes upwards
void SPMinimalnode::nodeAdoptRight(SPMinimalnode *rChild)
{
	// *this is the node which will become the parent

	// point parent and child pointers in the right directions
	// nodeAddRight() checks labels, hull size, present children
	nodeAddRight(rChild);
}

/* this is totally overcomplicated because it is designed for minimal
 * subpavings, ie we can't just check that something is in this box!*/

// test for ivector z contained in the subpaving represented by this
BOOL_INTERVAL SPMinimalnode::spContains(const ivector& z) const
{
	// z is assumed not to be empty
	// nb Intersection() gives error if unequal index sets passed

	BOOL_INTERVAL retValue = BI_FALSE; // for the return value

	// case of a non-empty leaf
	if (!isEmpty() && isLeaf()) {

		ivector r; // temporary,to be passed to Intersection

		if (z<=getBox()) {
			retValue = BI_TRUE;
		}

		// result is indeterminate if there is an
		// intersection but z is not wholly in theBox
		else if (Intersection(r, z, getBox())) {
			retValue = BI_INDET;
		}

		// Case that there is no intersection
		else retValue = BI_FALSE;

	} // end (!isEmpty() && isLeaf())

	//case of an non-empty non-leaf
	if (!isEmpty() && !isLeaf()) {
	//

		ivector Lz, Rz; // ivectors passed to Intersection()
		// will contain intersection after Intersection() call

		// to hold results of tests on left and right children
		BOOL_INTERVAL Ltest = BI_FALSE;
		BOOL_INTERVAL Rtest = BI_FALSE;

		// indicators for tested left and right sides
		bool LtestSuccess = false;
		bool RtestSuccess = false;

		// Find if there is a leftChild with a box
		if (hasLCwithBox() &&
			Intersection(Lz, z, leftChild->getBox())) {
			// Lz contains intersctn of z & leftChild box

			// test Lz and left child node
			Ltest = (getLeftChild()->spContains(Lz));
			LtestSuccess = true;
		}


		// Find if there is a rightChild with a box
		if (hasRCwithBox() &&
			Intersection(Rz, z, rightChild->getBox())) {
			// Rz contains intersctn of z & rightChild box

			// test Rz and right child node
			Rtest = (getRightChild()->spContains(Rz));
			RtestSuccess = true;
		}

		// if both children tested
		if (LtestSuccess && RtestSuccess) {
			//return value is the result of both tests
				// if the same or BI_INDET if diff
			Ltest==Rtest ?
				retValue = Ltest : retValue=BI_INDET;
		}

		// if has two children but neither was tested
		// ie neither Intersection() returned true
		if (hasRCwithBox() && hasLCwithBox()
			&& !LtestSuccess && !RtestSuccess) {
			retValue = BI_FALSE;
			// note that the AIA book has BI_TRUE here
			// but this can't be correct
		}

		// if has two children but only right was tested
		// ie left Intersection() returned false
		if (hasRCwithBox() && hasLCwithBox()
			&& !LtestSuccess && RtestSuccess) {
			retValue = Rtest;
			// return value result of test of right side
		}

		// if has two children but only left was tested
		// ie right Intersection() returned false
		if (hasRCwithBox() && hasLCwithBox()
			&& LtestSuccess && !RtestSuccess) {
			retValue = Ltest;
			// return value result of test of left side
		}

		// if has right child only and that child was tested
		// ie Intersection() returned true
		if (hasRCwithBox()
			&& !hasLCwithBox() && RtestSuccess) {
			// if all of z contained in right child's box
			if (Rz==z) {
				retValue = Rtest;
			}
			// return false if Rtest false, else INDET
			else {
				Rtest==BI_FALSE ? retValue = BI_FALSE
					: retValue = BI_INDET;
			}
		}

		// if has right child only and that child not tested
		// ie Intersection() returned false
		if (hasRCwithBox()
			&& !hasLCwithBox() && !RtestSuccess) {
			retValue = BI_FALSE;
		}

		// if has left child only and that child was tested
		// ie Intersection() returned true
		if (!hasRCwithBox() && hasLCwithBox()
			&& LtestSuccess) {
			// if whole of z contained in left child's box
			if (Lz==z) {
				retValue = Ltest;
			}
			// return false if Ltest false, otherwise INDET
			else {
				Ltest==BI_FALSE ?
					retValue = BI_FALSE :
					retValue = BI_INDET;
			}
		}

		// if has left child only & that child was not tested
		// ie Intersection() returned false
		if (!hasRCwithBox() && hasLCwithBox()
			&& !LtestSuccess) {
			retValue = BI_FALSE;
		}

		// case no children covered by isLeaf() block above

	} // end of (!isEmpty() && !isLeaf())

	// case of isEmpty() being true is take care of
	// by default return value being BI_FALSE

	return retValue;

} // end of spContains for ivector

/* this is totally overcomplicated because it is designed for minimal
 * subpavings, ie we can't just check that something is in this box!*/

// test for rvector p contained in subpaving represented by this node
BOOL_INTERVAL SPMinimalnode::spContains(const rvector& p) const
{
	// p is assumed not to be empty
	// nb Intersection() gives error if unequal index sets passed

	BOOL_INTERVAL retValue = BI_FALSE; // for the return value

	//cast p to an ivector
	ivector pvector = _ivector(p);

	// case of a non-empty leaf
	if (!isEmpty() && isLeaf()) {


		//find if p is in the box
		if (pvector <= getBox()) {
			retValue = BI_TRUE;
		}
		//else retValue keeps default value of BI_FALSE


	} // end (!isEmpty() && isLeaf())

	//case of an non-empty non-leaf
	if (!isEmpty() && !isLeaf()) {

		// to hold results of tests on left and right children
		BOOL_INTERVAL Ltest = BI_FALSE;
		BOOL_INTERVAL Rtest = BI_FALSE;

		// Find if there is a leftChild with a box
		if (hasLCwithBox()) {
			// test left child node
			Ltest = (getLeftChild()->spContains(p));
		}


		// Find if there is a rightChild with a box
		if (hasRCwithBox()) {
			// test Rz and right child node
			Rtest = (getRightChild()->spContains(p));
		}

		if ((Ltest==BI_TRUE) || (Rtest==BI_TRUE)) {
			retValue = BI_TRUE;
		}
		//else retValue keeps default value of BI_FALSE

		// case no children taken care of by isLeaf() above

	} // end of (!isEmpty() && !isLeaf())

	// case isEmpty() true covered by default retValue = BI_FALSE

	return retValue;

} // end of spContains for rvector




// Return a container of boxes represented by the finest common
// level of nodes between two subpavings.
// ie the 'outer jacket' that is the collection of smallest boxes that fits
// both of the inner subpavings
BoxVec& SPMinimalnode::vecLeafBoxOuterJacket(BoxVec& boxes,
				const SPMinimalnode * const spn1, 
				const SPMinimalnode * const spn2)
{
	// only do something if both nodes are non-null and boxes match
	if (spn1 != NULL && spn2 != NULL &&
					(spn1->getBox() == spn2-> getBox())) {

		// if both are leaves or
		// one of them does not have a child that the other has
		// then push back the box of the node we are in now
		if (spn1->isLeaf() && spn2->isLeaf()) {
			boxes.push_back(spn1->getBox());
		}
		else if ((spn1->getLeftChild() != NULL
					&& spn2->getLeftChild() == NULL)
					|| (spn1->getLeftChild() == NULL
					&& spn2->getLeftChild() != NULL)
					|| (spn1->getRightChild() != NULL
					&& spn2->getRightChild() == NULL)
					|| (spn1->getRightChild() == NULL
					&& spn2->getRightChild() != NULL)) {
			boxes.push_back(spn1->getBox());
		}

		else {
			// we recurse on the children if both have both children
			if (spn1->getLeftChild() != NULL
						&& spn2->getLeftChild() != NULL) {
				boxes = vecLeafBoxOuterJacket(boxes,
							spn1->getLeftChild(), spn2->getLeftChild());
			}
			if (spn1->getRightChild() != NULL
						&& spn2->getRightChild() != NULL) {

				boxes = vecLeafBoxOuterJacket(boxes,
							spn1->getRightChild(), spn2->getRightChild());
			}
		}
	}
	return boxes;
}

// Return a minimal subpaving representing the finest common
// level of nodes between two subpavings.
// ie the 'outer jacket' that is the finest subpaving that fits
// both of the inner subpavings
SPMinimalnode* SPMinimalnode::spLeafBoxOuterJacket(const SPMinimalnode * const spn1,
								const SPMinimalnode * const spn2)
{
	SPMinimalnode* jacketSP = NULL;
	try {
		
		BoxVec jacketBoxes;
		jacketBoxes = vecLeafBoxOuterJacket(jacketBoxes,
						spn1, spn2);
		if (!jacketBoxes.empty()) {

			ImageList listBoxes;
			listBoxes.insert(listBoxes.end(), jacketBoxes.begin(),
												jacketBoxes.end());
			ivector root = spn1->getBox();

			jacketSP = makeTreeFromLeaves(root, listBoxes);

			if (jacketSP != NULL) {
				jacketSP->setNodeName("X");
				jacketSP->recursiveRename();
			}

		}

		return jacketSP;
	}
	catch(std::exception const& e) {
		delete jacketSP;
		jacketSP = NULL;
		throw;
	}
}

// Return the sum of the volume of the outer jacket around two subpavings
double SPMinimalnode::volOuterJacket(const SPMinimalnode * const spn1,
							const SPMinimalnode * const spn2)
{
	BoxVec jacketBoxes;
	jacketBoxes = vecLeafBoxOuterJacket(jacketBoxes,
					spn1, spn2);
	BoxVecItr bit;
	double jacketVol = 0;
	for (bit = jacketBoxes.begin(); bit < jacketBoxes.end(); bit++) {
		jacketVol += Volume (*bit);
	}
	return jacketVol;
}

//Return a container of boxes represented by intersection of subpavings.
BoxVec& SPMinimalnode::vecLeafBoxIntersection(BoxVec& boxes,
				const SPMinimalnode * const spn1, 
				const SPMinimalnode * const spn2)
{
	// only do something if both nodes are non-null and boxes match
	if (spn1 != NULL && spn2 != NULL &&
					(spn1->getBox() == spn2-> getBox())) {

		// if both are leaves (and boxes match), push back the matching box
		if (spn1->isLeaf() && spn2->isLeaf()) boxes.push_back(spn1->getBox());

		// if one is a leaf and one is not, the children of the non-leaf
		// are all in the intersection
		if (!(spn1->isLeaf()) && spn2->isLeaf()) {
			SPnodeConstPtrs leaves1;
			spn1->getConstSPnodeLeaves(leaves1);
			SPnodeConstPtrsItr it;
			for (it = leaves1.begin(); it < leaves1.end(); it++) {
				boxes.push_back((*it)->getBox());
			}
		}

		if (spn1->isLeaf() && !(spn2->isLeaf())) {
			SPnodeConstPtrs leaves2;
			spn2->getConstSPnodeLeaves(leaves2);
			SPnodeConstPtrsItr it;
			for (it = leaves2.begin(); it < leaves2.end(); it++) {
				boxes.push_back((*it)->getBox());
			}
		}

		// if both have children, recurse on the children
		if (spn1->hasRCwithBox() && spn2->hasRCwithBox())
			boxes = vecLeafBoxIntersection(boxes, spn1->getRightChild(),
													spn2->getRightChild());
		if (spn1->hasLCwithBox() && spn2->hasLCwithBox())
			boxes = vecLeafBoxIntersection(boxes, spn1->getLeftChild(),
													spn2->getLeftChild());
	}

	return boxes;
}

//Return a subpaving representing intersection of subpavings.
SPMinimalnode* SPMinimalnode::spLeafBoxIntersection(const SPMinimalnode * const spn1,
								const SPMinimalnode * const spn2)
{
	SPMinimalnode* interSP = NULL;
	
	try {
		BoxVec interBoxes;
		interBoxes = vecLeafBoxIntersection(interBoxes,
						spn1, spn2);
		if (!interBoxes.empty()) {

			ImageList listBoxes;
			listBoxes.insert(listBoxes.end(), interBoxes.begin(),
												interBoxes.end());
			ivector root = spn1->getBox();

			interSP = makeTreeFromLeaves(root, listBoxes);

			if (interSP != NULL) {
				interSP->setNodeName("X");
				interSP->recursiveRename();
			}
		}

		return interSP;
	}
	catch(std::exception const& e) {
		delete interSP;
		interSP = NULL;
		throw;
	}
	
}

// Return the sum of the volume of intersection between two subpavings
double SPMinimalnode::volIntersection(const SPMinimalnode * const spn1,
							const SPMinimalnode * const spn2)
{
	BoxVec interBoxes;
	interBoxes = vecLeafBoxIntersection(interBoxes,
					spn1, spn2);
	BoxVecItr bit;
	double interVol = 0;
	for (bit = interBoxes.begin(); bit < interBoxes.end(); bit++) {
		interVol += Volume (*bit);
	}
	return interVol;
}

// Return a container of boxes represented by difference between subpavings,
// all the space in spn1 that is not in spn2
BoxVec& SPMinimalnode::vecLeafBoxDifference(BoxVec& boxes,
				const SPMinimalnode * const spn1, const SPMinimalnode * const spn2)
{
	if (spn1 != NULL && spn2 == NULL) {

		SPnodeConstPtrs leaves1;
		spn1->getConstSPnodeLeaves(leaves1);
		SPnodeConstPtrsItr it;
		for (it = leaves1.begin(); it < leaves1.end(); it++) {
			boxes.push_back((*it)->getBox());
		}
	}

	if (spn1 != NULL && spn2 != NULL
			&& (spn1->getBox() == spn2->getBox())) {

		if (spn1->isLeaf() && !(spn2->isLeaf())) {
			// now we want to get any part of the box of spn1 that is
			// not represented by whatever boxes hang off spn2
			boxes = vecBoxNodeDifference(boxes, spn1->getBox(), spn2);
		}

		// if spn1 is not a leaf but spn2 is, then there is nothing
		// represented in spn1 that is not already represented in spn2

		if (!(spn1->isLeaf()) && !(spn2->isLeaf())) {
			// recurse on the children
			boxes = vecLeafBoxDifference(boxes, spn1->getLeftChild(),
						spn2->getLeftChild());
			boxes = vecLeafBoxDifference(boxes, spn1->getRightChild(),
						spn2->getRightChild());
		}
	}

	// if spn1 is null and spn2 is not null we don't do anything
	// if both not null but the boxes don't match we don't do anything
	return boxes;
}

// Return a container of boxes represented by difference between box1 and
// the boxes represented by the leaves of spn2
BoxVec& SPMinimalnode::vecBoxNodeDifference(BoxVec& boxes,
	ivector box1, const SPMinimalnode * const spn2)
{
	if (spn2 == NULL) boxes.push_back(box1);

	else if (!(spn2->isLeaf())) {
		// bisect box1
		int maxdiamcomp = MaxDiamComp(box1);
		// new ivectors from splitting root along its biggest dimension
		ivector leftbox = Lower(box1, maxdiamcomp);
		ivector rightbox = Upper(box1, maxdiamcomp);
		boxes = vecBoxNodeDifference(boxes, leftbox, spn2->getLeftChild());
		boxes = vecBoxNodeDifference(boxes, rightbox, spn2->getRightChild());
	}
	// else spn2 must be a leaf so there is no difference

	return boxes;
}

//Return a subpaving representing the difference between subpavings.
SPMinimalnode* SPMinimalnode::spLeafBoxDifference(const SPMinimalnode * const spn1,
								const SPMinimalnode * const spn2)
{
	SPMinimalnode* diffSP = NULL;
	
	try {
		BoxVec diffBoxes;
		diffBoxes = vecLeafBoxDifference(diffBoxes, spn1, spn2);
		if (!diffBoxes.empty()) {

			ImageList listBoxes;
			listBoxes.insert(listBoxes.end(), diffBoxes.begin(),
												diffBoxes.end());
			ivector root = spn1->getBox();

			diffSP = makeTreeFromLeaves(root, listBoxes);

			if (diffSP != NULL) {
				diffSP->setNodeName("X");
				diffSP->recursiveRename();
			}
		}

		return diffSP;
	}
	catch(std::exception const& e) {
		delete diffSP;
		diffSP = NULL;
		throw;
	}
}

// Return the sum of the volume of difference between two subpavings
double SPMinimalnode::volDifference(const SPMinimalnode * const spn1,
							const SPMinimalnode * const spn2)
{
	BoxVec diffBoxes;
	diffBoxes = vecLeafBoxDifference(diffBoxes,
					spn1, spn2);
	BoxVecItr bit;
	double diffVol = 0;
	for (bit = diffBoxes.begin(); bit < diffBoxes.end(); bit++) {
		diffVol += Volume (*bit);
	}
	return diffVol;
}


// Forms a minimal image subpaving from boxes based on voxels
/*
Make a minimal subpaving tree from a list of interval vectors which approx
to the leaves of the tree.  The root of the subpaving tree will have
Box = root, and the boxes in the list will have the same width in each
dimension and be formable by a series of bisections of the given root.
i.e. each leaf will be the same size as each other leaf and will be a
'square' hypercube.
*/
SPMinimalnode* SPMinimalnode::makeTreeFromVoxels(ivector& root, ImageList& leafList,
				double spacing, size_t dim)
{
	SPMinimalnode* newNode = NULL;  // for return value
	
	try {

		if (!leafList.empty()) {
			//sort using the volCompare function
			leafList.sort(volCompare);   // sorts smallest to largest

		   // test if root is equal to the largest image element, ie the last
			bool isRootEqual = (root == *leafList.rbegin());

			// with given spacing, max width in any should be 1/spacing
			// so as soon as max root dimension is below 2/spacing, it should
			// be in the subpaving if it has any image data in it
			// so be conservative and take 1.5/spacing (between 1/ and 2/)
			double eps = 1.5/spacing;
			int maxdiamcomp = 0;  // to take value calculated from MaxDiam
			// find the maximum diameter, put the max dimension into maxdiamcomp
			double maxDiamRoot = MaxDiam(root, maxdiamcomp);

			// test if maximum root width is smaller than eps
			bool isRootSmall = (maxDiamRoot < eps);

			// if the list has some images in it
			// and either if the root is equal to the largest box in the list
			// or if the root max diameter is < eps
			// return a new node based on root
			if (isRootEqual || isRootSmall) {

			   newNode = new SPMinimalnode(root);
			   
			}
			// if the list has some images in it
			// and the root is not equal to the largest box in the list
			// and the root is not small
			// bisect the root, divide up the list, and recurse
			else {

				// new ivectors from splitting root along its biggest dimension
				ivector leftbox = Lower(root, maxdiamcomp);
				ivector rightbox = Upper(root, maxdiamcomp);

				// create two empty lists for the left and right side
				ImageList leftlist, rightlist;

				ImageListItr it; // iterator to for the list

				// find a conservative minimum expected volume for a leaf as
				// half of the product of 1/spacing over all dimensions
				// we expect leaves to be the result of successive bisections so
				// if the volume of the intersection is not what we expect
				// then the intersection will be just the boundary
				real minVol = 0.5;
				for (size_t i = 0; i < dim; i++) minVol /= spacing;

				// iterate through the current list, put the intersection of any
				// element with the leftbox into new left list, & intersection
				// of any element with the new right box into new right list
				for (it=leafList.begin(); it!=leafList.end(); it++) {
					ivector interLeft;  // intersection with left hull
					ivector interRight;  // intersection with right hull

					if (Intersection(interLeft, *it, leftbox)) {
						if (Volume(interLeft) > minVol)
							leftlist.push_back(interLeft);
					}

					if (Intersection(interRight, *it, rightbox)) {
						if (Volume(interRight) > minVol)
							rightlist.push_back(interRight);
					}

				} // end of iteration through list elements

				// recursively call makeTreeFromVoxels with leftbox, leftlist
				// and rightbox, rightlist
				// reunite the results using root as the box for parent node
				// makeTreeFromVoxels creates a minimal subpaving
				// (no sibling child nodes) on the root

				newNode = Reunite<SPMinimalnode>(makeTreeFromVoxels(leftbox,
												leftlist, spacing, dim),
									makeTreeFromVoxels(rightbox,
												rightlist, spacing, dim), root);

			} // end of is list has elements and first box does not contain root
		}

		// if there is nothing in the list we return the default
			// initialisation value of NULL

		return newNode;
	}
	catch(std::exception const& e) {
		delete newNode;
		newNode = NULL;
		throw;
	}
	

}

// makes a new paving from a vtk file
// expects 3d, structured point data in file
SPMinimalnode* SPMinimalnode::vtkPaving(const std::string filename)
{
	SPMinimalnode* newTree = NULL;
	
	try {

		size_t expectedDims = 3;

		IntVec Xs;
		IntVec Ys;
		IntVec Zs;

		// get the spacing and fill in the coordinates vectors using
		// getCoordinatesFromVtk.
		// getCoordinatesFromVtk expects 10 header lines with spacings on 4th
		IntVec spacing = getCoordinatesFromVtk(Xs, Ys, Zs, filename);
		bool success = ((spacing.size() == expectedDims) && (spacing[0] > 0)
								&& (spacing[0] == spacing[1])
								&& (spacing[0] == spacing[2]));

		if (success) {

			ivector rootbox(expectedDims);

			double maxXYZ = spacing[0] * 1.0;

			ImageList boxes;

			double totalListVol = 0;

			IntVecItr xIt = Xs.begin();
			IntVecItr yIt = Ys.begin();
			IntVecItr zIt = Zs.begin();

			for (xIt = Xs.begin(); xIt < Xs.end(); xIt++, yIt++, zIt++) {
				ivector box(expectedDims);
				// make the box so that its boundaries are on the inner edges
				// of the voxel
				interval xdim(Sup(_interval(*xIt/maxXYZ)),
							Inf(_interval(*xIt + 1.0)/maxXYZ));
				interval ydim(Sup(_interval(*yIt/maxXYZ)),
							Inf(_interval(*yIt + 1.0)/maxXYZ));
				interval zdim(Sup(_interval(*zIt/maxXYZ)),
							Inf(_interval(*zIt + 1.0)/maxXYZ));
				box[1] = xdim;
				box[2] = ydim;
				box[3] = zdim;
				boxes.push_back(box);
				totalListVol += Volume(box);

			}

			rootbox[1] = interval(0.0, 1.0);
			rootbox[2] = interval(0.0, 1.0);
			rootbox[3] = interval(0.0, 1.0);


			newTree = makeTreeFromVoxels(rootbox, boxes,
										maxXYZ, expectedDims);
		}

		if (newTree != NULL) {
			newTree->setNodeName("X");
			newTree->recursiveRename();
		}
	
		return newTree;
	}
	catch(std::exception const& e) {
		delete newTree;
		newTree = NULL;
		throw;
	}
	
}

// Forms a minimal image subpaving from leaf boxes
/*
Make a minimal subpaving tree from a list of interval vectors which are
the leaves of the tree.  The root of the subpaving tree will have
Box = root, and the boxes in the list will be formable by a series of
bisections of the given root.
*/
SPMinimalnode* SPMinimalnode::makeTreeFromLeaves(ivector& root, 
												ImageList& leafList)
{
	SPMinimalnode* newNode = NULL;  // for return value
	
	try {

		if (!leafList.empty())
		{
			//sort using the volCompare function
			leafList.sort(volCompare);   // sorts smallest to largest

			// test if root is equal to the largest image element, ie the last
			bool isRootEqual = (root == *leafList.rbegin());

			double smallestVol = Volume(*leafList.begin());

			//is root smaller than twice the size of the smallest box in it
			bool isRootSmall = (Volume(root) < 2*smallestVol);

			// if the list has some images in it
			// and if the root is equal to the largest box in the list
			// or the root is close enough
			// return a new node based on root
			if (isRootEqual || isRootSmall) {

				newNode = new SPMinimalnode(root);
				
			}
			// if the list has some images in it
			// and the root is not equal to the largest box in the list
			// and the root is not small
			// bisect the root, divide up the list, and recurse

			if (!isRootEqual && !isRootSmall) {

				int maxdiamcomp = 0;
				MaxDiam(root, maxdiamcomp);
				// new ivectors from splitting root along its biggest dimension
				ivector leftbox = Lower(root, maxdiamcomp);
				ivector rightbox = Upper(root, maxdiamcomp);

				// create two empty lists for the left and right side
				ImageList leftlist, rightlist;

				ImageListItr it; // iterator to for the list

				// iterate through the current list, put the intersection of any
				// element with the leftbox into new left list, & intersection
				// of any element with the new right box into new right list
				// but only put in whole leaves, ie with vol >= smallest vol
				for (it=leafList.begin(); it!=leafList.end(); it++) {
					ivector interLeft;  // intersection with left hull
					ivector interRight;  // intersection with right hull

					if (Intersection(interLeft, *it, leftbox)) {
						if (Volume(interLeft) > 0.5*smallestVol)
							leftlist.push_back(interLeft);
					}

					if (Intersection(interRight, *it, rightbox)) {
						if (Volume(interRight) > 0.5*smallestVol)
							rightlist.push_back(interRight);
					}

				} // end of iteration through list elements

				// recursively call makeTreeFromLeaves with leftbox, leftlist
				// and rightbox, rightlist
				// reunite the results using root as the box for parent node
				// makeTreeFromLeavess creates a minimal subpaving
				// (no sibling child nodes) on the root

				newNode = Reunite<SPMinimalnode>(makeTreeFromLeaves(leftbox,leftlist),
									makeTreeFromLeaves(rightbox, rightlist),
																root);

			} // end of is list has elements and first box not root
		}

		// if there is nothing in the list we return the default
			// initialisation value of NULL

		return newNode;
	}
	catch(std::exception const& e) {
		delete newNode;
		newNode = NULL;
		throw;
	}

}

void SPMinimalnode::swapMin(SPMinimalnode& spn) //throw() // don't hide base class version
{
	/* theBox, parent, leftChild,
	rightChild and nodeName are inherited from base class */
	SPnode::swap(spn); // use the base version
	
}

// protected methods


// ------------------- end of SPMinimalnode class definitions ---------------

// check for containment of ivector or box in the SubPaving
BOOL_INTERVAL subpavings::operator<=(const ivector& z, const SPMinimalnode * const spmn)
{
	BOOL_INTERVAL retValue = BI_FALSE;

	if (spmn!=NULL)
	{
		if (VecLen(z) != static_cast<int>(spmn->getDimension())) {
			throw IncompatibleDimensions_Error(
				"subpavings::operator<=(const ivector& , const SPMinimalnode * const )");
		}

		retValue = (*spmn).spContains(z);
	}
	return retValue;
	
}


// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::SPMinimalnode & s1, 
			subpavings::SPMinimalnode & s2) // throw ()
	{
		s1.swapMin(s2);
	}


/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
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

/*!/ \file:     spnode.cpp
\brief SPnode (SubPaving) and associated non-member function definitions
*/

#include "spnode.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <exception>

// include fstream so as to be able to output a file from spImage
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
#include "sptools.hpp"
#include "sptemplates.hpp"

// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

#include "spnodevisitor.hpp"

//src_trunk_0701
#include "subpaving_exception.hpp"

using namespace std;

namespace subpavings {

    // ----------------------- SPnode class definitions ------------------

//src_trunk_0701
// ------------------- private methods ------------------------------

//src_trunk_0701
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

	// -----------------protected methods
	
	// recursively split children according to instruction string
    std::string SPnode::splitLeft(std::string instruction)
    {
        std::string errorCode = "XX";
        // parse string to get first level out
        std::string comma = ", ";
        size_t startpos = instruction.find_first_not_of(comma);
        size_t endpos = std::string::npos;
        size_t newstartpos = std::string::npos;
        int depth = 0;
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
        depth = atoi(str.c_str()); // 0 if not valid integer
        if (depth == 0) {
            std::cerr << "Problem in parsing string" << endl;
            instruction = errorCode;
        }
        else { // depth should be > 0
            int myDepth = nodeName.length() - 1;
            if (myDepth < depth) { // need to go down more

                if (isLeaf()) {
                    // split, using the nodeExpand for this subtype if not base
                    nodeExpand();
                }

                // send splitLeft instruction down to the left
                instruction = getLeftChild()->splitLeft(instruction);
                if (instruction != errorCode) {
                    // send splitLeft instruction to right child
                    instruction = getRightChild()->splitLeft(instruction);
                }
            }
            else if (myDepth == depth) {     // split enough
                if (!isLeaf()) {
                    std::cerr << "Problem with instruction string:" << std::endl;
                    std::cerr << "cannot use with current histogram" << std::endl;
                    instruction = errorCode;
                }
                else if (newstartpos != std::string::npos) {
                    instruction = instruction.substr(newstartpos);
                }
                else instruction = "";
            }
            else { // myDepth is > depth
                std::cerr << "Problem with instruction string:" << std::endl;
                std::cerr << "cannot use with current histogram" << std::endl;
                instruction = errorCode;
            }
        }

        return instruction;
    }

	//src_trunk_0701
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



    // ---------------- public methods

    // default constructor
    SPnode::SPnode() :  theBox(NULL), dimension(0), label(0),
            parent(NULL), leftChild(NULL), rightChild(NULL), nodeName("X")
        {}


	//src_trunk_0701
	// initialised constructor, initialised with one ivector for the box
SPnode::SPnode(const ivector& v) : theBox(NULL),
	parent(NULL), leftChild(NULL), rightChild(NULL), nodeName("X")
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


//src_trunk_0701
// initialised constructor, initialised with a LabBox (labeled box)
SPnode::SPnode(const LabBox& lb) : theBox(NULL), parent(NULL),
	leftChild(NULL), rightChild(NULL), nodeName("X")
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



    // initialised constructor, initialised with one ivector for the box
    // and optionally initialised with lab for label
    // label defaults to 0 (see declaration)
    SPnode::SPnode(ivector& v, int lab) : label(lab),
        parent(NULL), leftChild(NULL), rightChild(NULL), nodeName("X")
    {
        try {
            theBox = new ivector(v);
            dimension = Ub(v) - Lb(v) + 1;
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with a LabBox (labeled box)
    SPnode::SPnode(LabBox& lb) : label(lb.L), parent(NULL),
        leftChild(NULL), rightChild(NULL), nodeName("X")
    {
        try {
            theBox = new ivector(lb.Box);
            dimension = Ub(lb.Box) - Lb(lb.Box) + 1;
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << std::endl;
            throw;
        }

    }


    //Copy constructor
    //copies from the given node downwards
    SPnode::SPnode(const SPnode& other) : dimension(other.dimension),
            label(other.label), parent(NULL), nodeName(other.nodeName)
    {
        try {
            theBox=new ivector(*other.theBox);

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
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << std::endl;
            throw;
        }
}


    // Destructor.
    SPnode::~SPnode()
    {delete theBox; delete leftChild; delete rightChild;}


    //copy assignment operator
    //copies from this node downwards
    SPnode& SPnode::operator=(const SPnode& rhs)
    {

        try {

            // delete the current children (deletes their children as well)
            if (leftChild != NULL) {
                delete leftChild;
                leftChild = NULL;
            }

            if (rightChild != NULL) {
                delete rightChild;
                rightChild = NULL;
            }
            // and delete the current box
            if (theBox != NULL) {
                delete theBox;
                theBox = NULL;
            }
            parent=NULL;

            theBox=new ivector(*rhs.theBox);
            dimension = rhs.dimension;
            label = rhs.label;
            nodeName = rhs.nodeName;

            //recursion on the children
            if (rhs.leftChild) {
                leftChild=new SPnode(*rhs.leftChild);
                leftChild->parent = this;
            }
            else leftChild=NULL;

            if (rhs.rightChild) {
                rightChild=new SPnode(*rhs.rightChild);
                rightChild->parent = this;
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << std::endl;
            throw;
        }

        return *this;

    }

        // recursively rename from the top down
    void SPnode::recursiveRename()
    {
        if(getLeftChild()!=NULL) {
            getLeftChild()->setNodeName(nodeName + "L");
            getLeftChild()->recursiveRename();
        }
        if(getRightChild()!=NULL) {
            getRightChild()->setNodeName(nodeName + "R");
            getRightChild()->recursiveRename();
        }
    }

	//from src_trunk_0701
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
    //{ return dimension; }

    // Accessor for theBox of a node.
    // Returns a copy of the object pointed to by theBox of a node.
    ivector SPnode::getBox() const
    { return *theBox; }

    // Accessor for label of a node.
    int SPnode::getLabel() const
    { return label; }

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

    // set the label
    //  will only set label for child-less root nodes
    void SPnode::setLabel(int lab)
    {
        if (isLeaf() && (getParent()==NULL)) {
            label = lab;
        }
        else {
            std::cout << "Cannot set label for node which "
                << "is not a childless root node"
                << std::endl;
        }
    }

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
        { return Volume(getBox()); }

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

	//src_trunk_0701
	// Return the volume of the box as a real.
	real SPnode::nodeRealVolume() const
	{ 
		if (isEmpty() ) {
			throw NoBox_Error("SPnode::nodeRealVolume()");
		}
		return realVolume(getBox()); 
	}

    // Check if this SPnode is empty.
    // Can only check if an actual subpaving object is empty, not if it is null.
    bool SPnode::isEmpty() const
    {return ( theBox==NULL) ; }

    // Check if this SPnode is a leaf.
    bool SPnode::isLeaf() const
    {return ( (leftChild==NULL) && (rightChild==NULL)); }

    // Check if this has a non-empty left child.
    bool SPnode::hasLCwithBox() const
    {return ( (leftChild!=NULL) && ((leftChild->theBox)!=NULL)); }

    // Check if this has a non-empty right child.
    bool SPnode::hasRCwithBox() const
    {return ( (rightChild!=NULL) &&
        ((rightChild->theBox)!=NULL)); }

//src_trunk_0701
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

	//src_trunk_0701
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



    // return a reference to a container of SPnodes
    // contents being the leaves descended from this, or this if this is a leaf
    // left to right order
    SPnodePtrs& SPnode::getSPnodeLeaves(SPnodePtrs& leaves) const
    {
        //if children, recurse on the children
        if (hasLCwithBox()) {
            getLeftChild()->getSPnodeLeaves(leaves);
        }

        if (hasRCwithBox()) {
            getRightChild()->getSPnodeLeaves(leaves);
        }

        if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
            leaves.push_back(const_cast<SPnode*>(this));
        }
        return leaves;
    }
    
    //src_trunk_0701
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
    
    // return a reference to a container of doubles
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

    //gat41 commented this out
    // fills in the leaf node levels, left to right
/*    IntVec& SPnode::getLeafNodeLevels(IntVec& levels) const
    {
        if (nodeName.length() > 0) {

            if (getLeftChild()!=NULL) {
                getLeftChild()->getLeafNodeLevels(levels);
            }
            if (getRightChild()!=NULL) {
                getRightChild()->getLeafNodeLevels(levels);
            }
            if (getLeftChild()==NULL && getRightChild()==NULL) {

                levels.push_back(nodeName.length()-1);
            }
        }
        else { // nodeName.length() == 0
            int level = 0;

            if (getLeftChild()!=NULL) {
                getLeftChild()->getLeafNodeLevels(level+1, levels);
            }
            if (getRightChild()!=NULL) {
                getRightChild()->getLeafNodeLevels(level+1, levels);
            }
            if (getLeftChild()==NULL && getRightChild()==NULL) {

                levels.push_back(level);
            }
        }
        return levels;
    }


    // fills in the leaf node levels, left to right, based on root node level
    // this version used when nodeNames not set, ie nodeName.length() == 0
    IntVec& SPnode::getLeafNodeLevels(const int level, IntVec& levels) const
    {

        if (getLeftChild()!=NULL) {
            getLeftChild()->getLeafNodeLevels(level+1, levels);
        }
        if (getRightChild()!=NULL) {
            getRightChild()->getLeafNodeLevels(level+1, levels);
        }
        if (getLeftChild()==NULL && getRightChild()==NULL) {

            levels.push_back(level);
        }
        return levels;
    }
*/

//src_trunk_0701
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
        IntVecItr it;

        std::string str = "";

        for (it = levels.begin(); it < levels.end(); it++) {
            stringstream out;
            out << (*it);
            str = str + out.str() + ",";
        };
        // trim off any lead or trailing ", "
        size_t startpos = str.find_first_not_of(" ,");
        size_t endpos = str.find_last_not_of(" ,");
        // if all commas or empty return an empty string
        if(( std::string::npos == startpos ) || ( std::string::npos == endpos)) {
            str = "";
        }
        else {
            str = str.substr( startpos, endpos-startpos+1 );
        }
        return str;
    }


    // Get the node depth.
    int SPnode::getNodeDepth() const
    {
        int depth = 0;
        if (nodeName.length() > 1)
            depth = nodeName.length() - 1;
        return depth;
    }


    //Returns the depth of the tree descending from this node
    int SPnode::getDepth() const
    {
        int depth = 0;

        // if this is a leaf the depth is 0
        if (!isLeaf()) {  // this is not a leaf
            // set up a container for the leaf children
            SPnodePtrs leaves;
            // fill the container with the leaf children
            getSPnodeLeaves(leaves);

            // find the deepest leaf child
            SPnodePtrsItr it;

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->getNodeDepth() > depth) {

                    depth = (*it)->getNodeDepth();
                }
            }
        } // end if not a leaf

        return depth;
    }


    //Returns the volume of the smallest (by vol) leaf node.
    double SPnode::getSmallestLeafVol() const
    {
        double smallestVol = 0.0;

        if (isLeaf()) {  // this is a leaf
            smallestVol = nodeVolume();
        }

        else { // this is not a leaf
            // set up a container for the leaf children
            SPnodePtrs leaves;
            // fill the container with the leaf children
            getSPnodeLeaves(leaves);

            // find the smallest child by volume
            SPnodePtrsItr it;
            SPnode* smallest = *(leaves.begin());

            smallestVol = smallest->nodeVolume();

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->nodeVolume() < smallestVol) {

                    smallestVol = (*it)->nodeVolume();
                }
            }
        } // end else not a leaf

        return smallestVol;
    }

    // Returns the volume of the largest (by vol) leaf node.
    double SPnode::getLargestLeafVol() const
    {
        double largestVol = 0.0;

        if (isLeaf()) {  // this is a leaf
            largestVol = nodeVolume();
        }

        else { // this is not a leaf

            // set up a container for the leaf children
            SPnodePtrs leaves;
            // fill the container with the leaf children
            // could be just this if no children
            getSPnodeLeaves(leaves);

            // find the largest child by volume
            SPnodePtrsItr it;
            largestVol = (*(leaves.begin()))->nodeVolume();

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->nodeVolume() > largestVol) {
                    largestVol = (*it)->nodeVolume();
                }
            }
        }

        return largestVol;
    }

//src_trunk_0701
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
        if (parent != NULL) {
            // check the final letter on my name
            std::string mySide = nodeName.substr(nodeName.length()-1);

            // siblings cannot be assumed to exist
            if((mySide == "R" && getParent()->hasLCwithBox()
                    && getParent()->getLeftChild()->isLeaf())
                ||
                (mySide == "L" && getParent()->hasRCwithBox()
                    && getParent()->getRightChild()->isLeaf())) {

                retValue = true;
            }
        }
        return retValue;
    }

	void SPnode::accept(SPnodeVisitor& visitor)
    {
        visitor.visit(this);
    }

	//gat41
	void SPnode::collectRange(SPnodeVisitor& visitor)
    {
        visitor.tellMe(this);
    }

	//src_trunk_0701
	void SPnode::acceptSPCheckVisitor(const SPCheckVisitor& visitor) const
{
	visitor.visit(this);
	
}

    // Print the box of a specific node in a subpaving
    std::ostream& SPnode::nodePrint(std::ostream &os) const
    {
    // output in form: box

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
                os << "[ " << Inf(thisBox[i]) << " , "
                    << Sup(thisBox[i]) << " ]";
                if (i < Ub(thisBox) ) {
                    os << " , ";
                }
            }
            os << std::endl;
        }
        return os;

    }

    // Print the details of a single leaf node, using tab delimiters
    std::ostream& SPnode::leafOutputTabs(std::ostream &os) const
    {

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            // output the label, nodeVolume
            os << label;
            os << "\t" << nodeVolume();

            // followed by intervals making up box using Inf & Sup
            // ie unlike cxsc output, there is no [  ] around them
            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i]) << "\t"
                    << Sup(thisBox[i]);
            }

        }

    }

    //Output for all the  leaf boxes in this, using tab delimiters
    std::ostream& SPnode::leavesOutputTabs(std::ostream &os) const
    {
        // uses  member function leafOutputTabs to generate node output
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafOutputTabs(os);
            return (os << "\n");

        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesOutputTabs(os);
        }

        if (getRightChild()!=NULL) {
            getRightChild()->leavesOutputTabs(os);
        }

    }
    
    //src_trunk_0701
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

        if (!isEmpty()) {
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

//src_trunk_0701
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

//src_trunk_0701


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




    // test for ivector z contained in the subpaving represented by this
    BOOL_INTERVAL SPnode::spContains(const ivector& z) const
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
                Ltest = ((*leftChild).spContains(Lz));
                LtestSuccess = true;
            }


            // Find if there is a rightChild with a box
            if (hasRCwithBox() &&
                Intersection(Rz, z, rightChild->getBox())) {
                // Rz contains intersctn of z & rightChild box

                // test Rz and right child node
                Rtest = ((*rightChild).spContains(Rz));
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


    // test for rvector p contained in subpaving represented by this node
    BOOL_INTERVAL SPnode::spContains(const rvector& p) const
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
                Ltest = ((*leftChild).spContains(p));
            }


            // Find if there is a rightChild with a box
            if (hasRCwithBox()) {
                // test Rz and right child node
                Rtest = ((*rightChild).spContains(p));
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



    // add two sibling nodes to this provided that this node is a leaf
    // comp argument is passed to Upper() and Lower()
    // split the box in half normal to dimension set by comp
    void SPnode::nodeExpand(int comp)
    {
        try
        {
            // only do something if this SPnode is a leaf
            if (isLeaf()) {

                // ivectors to be new boxes for new children
                ivector lC, rC;

                // Call Lower() and Upper() to put split boxes
                // into lC and rC respectively
                Lower(getBox(), lC, comp);
                Upper(getBox(), rC, comp);

                // make children and make this their parent
                leftChild = new SPnode(lC, label);
                leftChild->parent = this;
                rightChild = new SPnode(rC, label);
                rightChild->parent = this;

                //name the new children
                getLeftChild()->setNodeName(nodeName + "L");
                getRightChild()->setNodeName(nodeName + "R");

            }
        }

        catch (bad_alloc&)
        {

            std::cout << "Error allocating memory "
                << "in SPnode::nodeExpand()" << std::endl;
            throw;
        }

    }

    // finds the dimension to split on and then calls nodeExpand(int comp)
    void SPnode::nodeExpand()
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(maxdiamcomp); // complete nodeExpand

    }
    
     // finds the dimension to split on and then calls nodeExpand(int comp)
    // will also bring validation data down
    void SPnode::nodeExpand(bool boolVal)
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(maxdiamcomp);
    }

    // reabsorb children into this
    // just needs to delete the children
    void SPnode::nodeReabsorbChildren()
    {
        delete leftChild;
        delete rightChild;
        leftChild = NULL;
        rightChild = NULL;

    }


    // computes a minimal subpaving from two sibling subpavings
    // a subpaving is minimal if it has no sibling leaves
    // a minimal subpaving is created by discarding sibling leaves
    // lChild and rChild are the two subpavings to be reunited
    void SPnode::nodeReunite(SPnode *lChild, SPnode *rChild)

    {
        // *this is the node which will become the parent

        // check that the labels match and exit if not
        if ((lChild->label != label ) || (rChild->label != label)) {
             throw SPnodeException("Box labels do not match");
        }

        // if both subpavings are leaves and hull of boxes is x,
        // discard them: *this is a leaf
        if (lChild->isLeaf() && rChild->isLeaf()) {
            if (*theBox !=
                (*(lChild->theBox) | *(rChild->theBox))) {
                throw SPnodeException("Boxes to reunite not = this box");
            }

            //discard the two subpavings given
            delete lChild;
            delete rChild;

        }

        else {  // at least one of the children is not a leaf
            // this has to adopt them rather than reuniting them
            nodeAdoptLeft(lChild);
            nodeAdoptRight(rChild);
            recursiveRename(); // recursively rename child branches
        }

    }

 /*
    // add lChild onto this node
    void SPnode::nodeAddLeft(SPnode *lChild)
    {
        // *this is the node which will become the parent

        // check that the labels match and exit if not
        if (lChild->label != label) {
            throw SPnodeException("Labels on boxes do not match");

        }

        // check  box for this is union of boxes of proposed children
        if (hasLCwithBox()) {
            throw SPnodeException("Cannot adopt left with existing left child");

        }
        if (hasRCwithBox()) {
            // check that this new box fits with the current box
            if (*theBox !=
                (*(lChild->theBox) | *(getRightChild()->theBox))) {
                    throw SPnodeException("Child boxes do not fit together");
            }

        }

        this->leftChild = lChild;
        leftChild->parent = this;
    }

    // add Child onto this node
    void SPnode::nodeAddRight(SPnode *rChild)
    {
        // *this is the node which will become the parent

        // check that the labels match and exit if not
        if (rChild->label != label) {
            throw SPnodeException("Box labels do not match");

        }

        // check  box for this is union of boxes of proposed children
        if (hasRCwithBox()) {
            throw SPnodeException("Cannot adopt right with existing right child");

        }
        if (hasLCwithBox()) {
            // check that this new box fits with the current box
            if (*theBox !=
                (*(rChild->theBox) | *(getLeftChild()->theBox))) {
                throw SPnodeException("Child boxes do not fit together");
            }

        }

        this->rightChild = rChild;
        rightChild->parent = this;

    }
*/

//src_trunk_0701
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

//src_trunk_0701
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


    // graft lChild onto this node
    // lChild could be a leaf or a non-leaf
    // takes care of the data associated with lChild/its descendents
    // used when we are building a subpaving from the leaf nodes upwards
    void SPnode::nodeAdoptLeft(SPnode *lChild)
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
    void SPnode::nodeAdoptRight(SPnode *rChild)
    {
        // *this is the node which will become the parent

        // point parent and child pointers in the right directions
        // nodeAddRight() checks labels, hull size, present children
        nodeAddRight(rChild);
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

	//src_trunk_0701
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

	// split a root box to a shape specified by the instruction string
    bool SPnode::splitRootToShape(std::string instruction)
    {
        bool success = false;
        if (parent == NULL) {
            std::string res = splitLeft(instruction);
            if (res == "") success = true;
        }
        else {
            std::cerr << "Sorry, not a root paving" << std::endl;
        }
        return success;
    }

    // Return a container of boxes represented by the finest common
    // level of nodes between two subpavings.
    // ie the 'outer jacket' that is the collection of smallest boxes that fits
    // both of the inner subpavings
    BoxVec& SPnode::vecLeafBoxOuterJacket(BoxVec& boxes,
                    const SPnode * const spn1, const SPnode * const spn2)
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

    // Return a container of boxes represented by the finest common
    // level of nodes between two subpavings.
    // ie the 'outer jacket' that is the finest subpaving that fits
    // both of the inner subpavings
    SPnode* SPnode::spLeafBoxOuterJacket(const SPnode * const spn1,
                                    const SPnode * const spn2)
    {
        SPnode* jacketSP = NULL;
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

    // Return the sum of the volume of the outer jacket around two subpavings
    double SPnode::volOuterJacket(const SPnode * const spn1,
                                const SPnode * const spn2)
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
    BoxVec& SPnode::vecLeafBoxIntersection(BoxVec& boxes,
                    const SPnode * const spn1, const SPnode * const spn2)
    {
        // only do something if both nodes are non-null and boxes match
        if (spn1 != NULL && spn2 != NULL &&
                        (spn1->getBox() == spn2-> getBox())) {

            // if both are leaves (and boxes match), push back the matching box
            if (spn1->isLeaf() && spn2->isLeaf()) boxes.push_back(spn1->getBox());

            // if one is a leaf and one is not, the children of the non-leaf
            // are all in the intersection
            if (!(spn1->isLeaf()) && spn2->isLeaf()) {
                SPnodePtrs leaves1;
                spn1->getSPnodeLeaves(leaves1);
                SPnodePtrsItr it;
                for (it = leaves1.begin(); it < leaves1.end(); it++) {
                    boxes.push_back((*it)->getBox());
                }
            }

            if (spn1->isLeaf() && !(spn2->isLeaf())) {
                SPnodePtrs leaves2;
                spn2->getSPnodeLeaves(leaves2);
                SPnodePtrsItr it;
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
    SPnode* SPnode::spLeafBoxIntersection(const SPnode * const spn1,
                                    const SPnode * const spn2)
    {
        SPnode* interSP = NULL;
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

    // Return the sum of the volume of intersection between two subpavings
    double SPnode::volIntersection(const SPnode * const spn1,
                                const SPnode * const spn2)
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
    BoxVec& SPnode::vecLeafBoxDifference(BoxVec& boxes,
                    const SPnode * const spn1, const SPnode * const spn2)
    {
        if (spn1 != NULL && spn2 == NULL) {

            SPnodePtrs leaves1;
            spn1->getSPnodeLeaves(leaves1);
            SPnodePtrsItr it;
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
    BoxVec& SPnode::vecBoxNodeDifference(BoxVec& boxes,
        ivector box1, const SPnode * const spn2)
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
    SPnode* SPnode::spLeafBoxDifference(const SPnode * const spn1,
                                    const SPnode * const spn2)
    {
        SPnode* diffSP = NULL;
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

    // Return the sum of the volume of difference between two subpavings
    double SPnode::volDifference(const SPnode * const spn1,
                                const SPnode * const spn2)
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
    SPnode* SPnode::makeTreeFromVoxels(ivector& root, ImageList& leafList,
                    double spacing, size_t dim)
    {
        SPnode* newNode = NULL;  // for return value

        if (!leafList.empty())
        {
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

               try {
                    newNode = new SPnode(root);
               }
                catch (bad_alloc&)
                {
                    std::cout << "Error allocating memory "
                        << "in makeTreeFromVoxels(...)" << std::endl;
                    throw;
                }
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

                newNode = Reunite<SPnode>(makeTreeFromVoxels(leftbox,
                                                leftlist, spacing, dim),
                                    makeTreeFromVoxels(rightbox,
                                                rightlist, spacing, dim), root);

            } // end of is list has elements and first box does not contain root
        }

        // if there is nothing in the list we return the default
            // initialisation value of NULL

        return newNode;

    }

    // makes a new paving from a vtk file
    // expects 3d, structured point data in file
    SPnode* SPnode::vtkPaving(const std::string filename)
    {
        SPnode* newTree = NULL;

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

    // Forms a minimal image subpaving from leaf boxes
    /*
    Make a minimal subpaving tree from a list of interval vectors which are
    the leaves of the tree.  The root of the subpaving tree will have
    Box = root, and the boxes in the list will be formable by a series of
    bisections of the given root.
    */
    SPnode* SPnode::makeTreeFromLeaves(ivector& root, ImageList& leafList)
    {
        SPnode* newNode = NULL;  // for return value

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

                try {
                    newNode = new SPnode(root);
                }
                catch (bad_alloc&)
                {
                    std::cout << "Error allocating memory "
                        << "in makeTreeFromVoxels(...)" << std::endl;
                    throw;
                }
            }
            // if the list has some images in it
            // and the root is not equal to the largest box in the list
            // and the root is not small
            // bisect the root, divide up the list, and recurse

            if (!isRootEqual && !isRootSmall) {

                int maxdiamcomp = 0;
                double rootDiam = MaxDiam(root, maxdiamcomp);
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

                newNode = Reunite<SPnode>(makeTreeFromLeaves(leftbox,leftlist),
                                    makeTreeFromLeaves(rightbox, rightlist),
                                                                root);

            } // end of is list has elements and first box not root
        }

        // if there is nothing in the list we return the default
            // initialisation value of NULL

        return newNode;

    }


	//src_trunk_0701
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


    // ------------------- end of SPnode class definitions ---------------

    // ----------------------------- SPnode exceptions definitions

    SPnodeException::SPnodeException(std::string ss) : s(ss) {}
    SPnodeException::~SPnodeException () throw () {}
    const char* SPnodeException::what() const throw() { return s.c_str(); }



    // -------------------- start of SPnode non-member functions

    //Output the leaf boxes in SubPaving spn
    std::ostream & operator<<(std::ostream &os, const SPnode * const spn)
    {
        // uses nodePrint to generate node output

        if (!(isEmpty(spn)) && isLeaf(spn)) { // spn is non-empty leaf
            spn->nodePrint(os);
        }

        //recurse on the children
        if (!(isEmpty(spn)) && !(isLeaf(spn))) {
            os << (spn->getLeftChild());
            os << (spn->getRightChild());
        }

        //in the case where spn is empty we just return os

        return os;
    }

    // check for containment of ivector or box in the SubPaving
    BOOL_INTERVAL operator<=(const ivector& z, const SPnode * const spn)
    {
        BOOL_INTERVAL retValue = BI_FALSE;

        if (spn!=NULL)
        {
            if (VecLen(z) != spn->getDimension()) {
                throw SPnodeException("Dimension do not match");
            }

            retValue = (*spn).spContains(z);
        }
        return retValue;
    }


    // check if a SubPaving is a leaf
    bool isLeaf(const SPnode * const spn)
    {
        // FALSE if spn is a null pointer, true if spn is not NULL
        bool retVal = (spn!=NULL);


        if (retVal) { // if spn points to a non-empty node
            retVal = spn->isLeaf();
        }

        return retVal;
    }

    // check if a SubPaving is empty
    bool isEmpty(const SPnode *const spn)
    {
        // return true if spn is a null pointer or
        // node spn points to is empty
        return ((spn==NULL) || (spn->isEmpty()));
    }

    // get volume of a SubPaving
    double spVolume(const SPnode * const spn)
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
    size_t spLeaves(const SPnode * const spn)
    {
        size_t retVal=0;

        if (!(isEmpty(spn)) && isLeaf(spn)) retVal = 1; // leaf

        // recurse on children
        if (!(isEmpty(spn)) && !(isLeaf(spn))) {

            retVal += spLeaves(spn->getLeftChild());
            retVal += spLeaves(spn->getRightChild());

        }

        // if spn is empty we will return retVal = 0 by default

        return retVal;
    }

    //gloria's additions=======================================
	     // get number of nodes of a SubPaving
    size_t spTotalNodes(const SPnode * const spn)
    {
        size_t retVal=0;
        if (!(isEmpty(spn))) retVal = 1; // leaf
        // recurse on children
        if (!(isEmpty(spn))) {
            retVal += spTotalNodes(spn->getLeftChild());
            retVal += spTotalNodes(spn->getRightChild());
        }
        // if spn is empty we will return retVal = 0 by default
        return retVal;
    }
	 //======end of gloria's additions==============================

} // end namespace subpavings

	//src_trunk_0701
	// Full specializations of the templates in std namespace can be added in std namespace.
	template <>
	void std::swap(subpavings::SPnode & s1, 
				subpavings::SPnode & s2) // throw ()
	{
		s1.swap(s2);
	}

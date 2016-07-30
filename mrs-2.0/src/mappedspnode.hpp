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
\brief MappedSPnode definitions.
*/

#ifndef __MAPPEDSP_HPP__
#define __MAPPEDSP_HPP__

// put it all in the header for the moment and sort out the template issues later

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "toolz.hpp"

#include "sptools.hpp"

#include "spnode.hpp"

#include "SmallClasses.hpp"

#include "sp_expand_visitor.hpp"

#include "sp_value_visitor.hpp"

#include "subpaving_exception.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <set>

#include <ctime>
#include <gsl/gsl_rng.h>

//#define MYDEBUGVISITOR
//#define MYDEBUGVISITORE
//#define SLICE_OUTPUT

// capturing a leaf box that a slice goes through
//#define CAPTURE_SLICE_BOX

//#define LOGARITHMETIC // log state from addition operation

#ifdef NDEBUG
	#undef MYDEBUGVISITOR
	#undef MYDEBUGVISITORE
	#undef SLICE_OUTPUT

#endif


#ifdef LOGARITHMETIC
	#include <fstream>
	#include <sstream>
	
#endif


namespace subpavings {
	
	
	/*! \brief A templated derived class based on SPnode.

    The base class SPnode is a node in the representation of a regular
    subpaving as a binary tree.  A node represents a box (interval vector).
    SPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.
     
    The MappedSPnode template is parameterised on the type T
	of the range data mapped to (associated with) each node 
    (i.e. with each element of the
    subpaving).  Examples of range data types include reals, or 
    intervals, or colours ...
    
	A MappedSPnode node has a single value of type T mapped
	onto it.  This value is referred to as the "range value" or "range".
	
    This templatised class provides ways to perform some basic
    arithmetic-like operations on the subpaving and the data mapped
    to it.
	
	It is probably best to regard this class as a class for non-minimal
	node-trees, ie a node either has two children or no children.  The
	expand method will automatically make two children and there is 
	no way of then creating a minimal tree by cutting off an 'unwanted'
	child.  This means that the whole root paving is regarded as 
	the 'support' of the mapping. 

    */
	
	
    template <typename T>
    class MappedSPnode : public SPnode {

	

	public:


    // ------------------------ public member functions ----------------

    /*! \brief Destructor.*/
    virtual ~MappedSPnode()  {}



    /*! \brief No-argument constructor. */
    MappedSPnode() {}   // uses the base SPnode class default constructor


    /*! \brief Initialised constructor.
    
    Initialised with a box.*/
    explicit MappedSPnode(const ivector& v)
        : SPnode(v) {}


    /*! \brief Initialised constructor
    
    Initialised with a labeled box.*/
    explicit MappedSPnode(const LabBox& lb)
        : SPnode(lb) {}



    /*! \brief Initialised constructor
    
    Initialised with a box and a value of T for the range.*/
    MappedSPnode(const ivector& v, const T& r)
        : SPnode(v), range(r)
    {}


    /*! \brief Initialised constructor.
    
    Initialised with a labeled box and a value of T for the range.*/
    MappedSPnode(const LabBox& lb, const T& r)
        : SPnode(lb), range(r)
    {}


    /*! \brief Copy constructor.
	
	Copies from a plain SPnode.
    
    Copies from given node downwards.*/
    explicit MappedSPnode(const SPnode& other)
        : SPnode()
    {
		if (other.theBox != NULL) {
			theBox = new ivector( other.getBox() );
		}
	
		nodeName = other.nodeName;

		//recursion on the children
		if (other.leftChild) {
			nodeAddLeft(new MappedSPnode(
				*(other.getLeftChild())));
		}
		else leftChild=NULL;

		if (other.rightChild) {
			nodeAddRight(new MappedSPnode(
				*(other.getRightChild())));
		}
		else rightChild=NULL;
        
    }
	
	/*! \brief Copy constructor.
	
	Copies from a %MappedSPnode. 
    
    Copies from given node downwards.*/
    MappedSPnode(const MappedSPnode& other)
        : SPnode()
    {
		if (other.theBox != NULL) {
			theBox = new ivector( other.getBox() );
		}
	
		range = other.range;
		nodeName = other.nodeName;

		//recursion on the children
		if (other.leftChild) {
			nodeAddLeft(new MappedSPnode(
				*(other.getLeftChild())));
		}
		else leftChild=NULL;

		if (other.rightChild) {
			nodeAddRight(new MappedSPnode(
				*(other.getRightChild())));
		}
		else rightChild=NULL;
        
    }
	
	/*! \brief Copy assignment operator.
    
    Copies from given node downwards.*/
    MappedSPnode<T>& operator=(MappedSPnode<T> rhs)
    {
        rhs.swapMSPSR(*this); // make sure we use our version of swap
		return(*this);
    }
	
	/*! \brief Replace the properties of this node and its descendents
	with the properties of another node and its descendents.
	
	Copies \a newNode node downwards into
	this, but keeps the relationship of this with 
	its parent, ie if this is a node in a tree, the part of the tree
	rooted at this is made identical to \a newNode but the relationship
	to the rest of the tree,through the link between this and its
	parent, is retained. */
	void replaceMe(MappedSPnode<T> newNode)
	{
		
		MappedSPnode<T>* pNode = getParent();
		newNode.swapMSPSR(*this);
		
		if (NULL != pNode) {
			
			this->parent = pNode;
		}
	}



    // parent and child accessors have to hide the base class implementation
    // this is not good but otherwise we get the base class return type
    // I've asked around and I can't find a way around it ...

    /*! \brief Accessor for the parent of a node.
     
    \return a copy of the pointer to parent node, cast to this type.*/
    MappedSPnode<T>* getParent() const
    { return (MappedSPnode<T>*) parent; }


    /*! \brief Accessor for the left child of a node.
    
    \return a copy of the pointer to leftChild node, cast to this type.*/
    MappedSPnode<T>* getLeftChild() const
    { return (MappedSPnode<T>*) leftChild; }


    /*! \brief Accessor for the right child of a node.
    
    \return a copy of the pointer to rightChild node, cast this type.*/
    MappedSPnode<T>* getRightChild() const
    { return (MappedSPnode<T>*) rightChild; }


    /*! \brief Accessor for the range.
    
    \return a copy of the range.*/
    T  getRange() const
    {   
        return T(range);
    }

	/*! \brief Set the range.
    
    \param r the value for the range.*/
    void  setRange(T r)
    {   
        range = r;
    }



   
	/*! \brief Get a string summary of this.
	
	Used for debugging.*/
	std::string nodeStringSummary() const
	{
		std::ostringstream oss;
		
		oss << "I am " << getNodeName() << "(address " << this << "),\n";
		oss << "Dimension is " << getDimension() << ", address of box is " << theBox << "\n"; 
		oss << "range is " << range << "\n";
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
	
	
    /*! \brief Recursively allocate a collection of ranges to this and children.
	
	Throws a std::invalid_argument if there are fewer values in 
	\a rangesToAllocate than there are nodes in this.
    
    Allocation order is this, left child with remainder of allocation, 
    right child with remainder.
	
	\param rangesToAllocate a collection of range values to apply to the 
	nodes of the tree rooted at this.
	\param index is the index of the value in \a rangesToAllocate to
	become the value of the range for this.
	\pre rangesToAllocate.size() - index >= number of nodes in tree
	rooted at this.
	\post The nodes in the tree rooted at this have values taken from 
	the values in \a rangesToAllocate, starting at the value indexed
	\a index which is applied to this, then allocating the following
	values to the left child of this and its descendants 
	and then the right child of this and its descendants.*/
    size_t allocateRanges(const std::vector< T >& rangesToAllocate, size_t index = 0)
    {
		if (index >= rangesToAllocate.size()) {
			throw std::invalid_argument("Range allocations too short");

		}

		range = rangesToAllocate[index];

		std::size_t newIndex = index+1;
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
				throw std::invalid_argument("More ranges than nodes in tree");
			}

		}

		return newIndex;
	
    }


    /*! \brief Splits paving according to string instruction.
	
	The instruction specifies the required shape in terms of the depth of
	the leaf nodes, in left to right order.  The depth of a leaf node is
	equivalent to the number of bisections of the root box required to make
	the box represented by that leaf.  i.e., the root has depth 0 and if
	that were bisected, the root node would have two child nodes each of
	level 1.  Bisection of any one of the boxes represented by a
	child would give two more children, each of level 2 (2 bisections), etc.

	Leaf levels in \a instruction can be separated by commas, spaces, or 
	both.
	
	For example, an instruction string  "3, 3, 2, 1" would give an SPSnode
	tree with 4 leaves, 2 of which would be level 3 (i.e. representing
	boxes resulting from 3 bisections of the root box, each of which would
	have volume 1/8 the size of the root box).  Another leaf would represent
	a box resulting from 2 bisections of the root box (volume 1/4 that of
	the root box) and the 'right-most' leaf (in a drawing of the tree) would
	be the result of a single bisection of the root box and would have half
	the volume of the root box.  This is a valid instruction string because
	it is possible to get leaves of those levels by a series of successive
	bisections of the root box and the volume of the leaves adds up to the
	volume of the root box.
	
	Throws the following exceptions:
	<ul>
	<li>Throws a NoBox_Error if the box is NULL.</li>
	<li>Throws a NonRootNode_Error if this is not a root node
	(ie if this has a parent node).</li>
	<li>Throws an std::invalid_argument exception if the \a
	instruction constain invalid characters (anything other
	than digits, whitespace or commas).</li> 
	<li>Throws an std::logic_error if the instruction is valid but does
	not describe an achievable tree shape.</li>
	</ul>

	\param instruction specifies the required shape, eg "3, 3, 2, 1".
	\return true if the instruction could be successfully carried out,
	false if the instruction could not be carried out successfully.
	\pre The instruction must be valid and describe an achievable
	tree shape.  This must be a non-empty root node.  
	\post getLeafNodeLevelsString() == instruction */
	bool splitToShape(std::string instruction)
    {
		bool success = false;

		// checks:  is this the root?
		if (NULL != parent) {
			throw NonRootNode_Error("MappedSPNodeSingleRange<T>::splitToShape(std::string)");
		}

		// checks:  is there a root box
		if (NULL == theBox) {
			throw NoBox_Error("MappedSPNodeSingleRange<T>::splitToShape(std::string)");
		}


		// checks: is the string properly formed?
		if (instruction.length() == 0) {
			throw std::invalid_argument(
			"MappedSPNodeSingleRange<T>::splitToShape(std::string) : No instruction");
		}
		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
			"MappedSPNodeSingleRange<T>::splitToShape(std::string) : Illegal character");
		}

		// all seems to be okay, we can start spliting the root
		// specify what to look for as numbers or decimal point or + or -

		success = splitRootToShape(instruction);

		if (!success) {
			throw std::logic_error(
			"MappedSPNodeSingleRange<T>::splitToShape(std::string) : instruction not a proper tree");
			
	   }
	   
		return success;
   }


    /*! \brief Add two sibling nodes to this provided this is a leaf.
    
    Creates two children each with half of the box of this.    
    Split the box in half normal to dimension set by comp.
    comp argument is passed to Upper() and Lower()
    
    \param comp the dimension on which to bisect this to create children.*/
    virtual void nodeExpand(int comp)
    {
        // can only expand if there is a box
		if (isEmpty()) {
			throw NoBox_Error("MappedSPnode<T>::nodeExpand(int)");
		}
		
		// only do something if this node is a leaf
		if (isLeaf()) {
			
			MappedSPnode<T>* newLC = NULL;
			MappedSPnode<T>* newRC = NULL;
			
			try {
				// ivectors to become boxes for new children
				ivector lC, rC;
				// Call Lower() and Upper() to put split boxes
				// into lC and rC respectively
				Lower(getBox(), lC, comp);
				Upper(getBox(), rC, comp);

				// make and add the new children
				newLC = new MappedSPnode<T>(lC, range);
				newRC = new MappedSPnode<T>(rC, range);
				
				nodeAddLeft(newLC);
				nodeAddRight(newRC);
				// both children get the same range as this
				
				//name the new children
				getLeftChild()->setNodeName(nodeName + "L");
				getRightChild()->setNodeName(nodeName + "R");

				// new children have range collection from this
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


    /*! \brief Add two sibling nodes to this provided this is a leaf.
     
    Creates two children each with half of the box of this.*/
    // as base class
    void nodeExpand()
    {
        SPnode::nodeExpand();
    }

	/*! \brief Slice this.
	
	Slice this on at the point jointly specified by \a sliceDims
	and \a slicePts.  For example, if this has a 3-dimensional 
	root box [-1,1]x[-1,1]x[-1,1], ie dimensions {1,2,3}
	and \a sliceDims = {2,3} and 
	\a slicePts is (0.0,0.5) then we are slicing at point 0.0 on 
	dimension 2, point 0.5 on dimension 3.  After the opertation this
	would then have only one-dimensional boxes on dimensions 
	{1,2,3}\{2,3} = {1},
	each box being one that did contain point 0.0 on 
	dimension 2, point 0.5 on dimension 3.  The range will be 
	unchanged).  Any boxes 
	associated with the original this that did not contain 
	the point 0.0 on dimension 2 and 0.5 on dimension 3 will not be
	represented in its new form.  
	 
	\param sliceDims is a vector of dimensions to slice on, indexed 
	from 1 onwards.
	\param slicePts is a vector of points to slice on, assumed to 
	correspond to the dimensions in \a sliceDims, ie the ith value
	in \a sliceDims gives the dimension for the ith point in
	\a slicePts.
	\param sliceFilename is the name of file to use to capture
	the boxes of this that are used in the slice by outputting 
	them to the file named \a sliceFilename.  Defaults to 
	the empty string "".  If \a sliceFilename is the empty string
	"" no boxes will be captured to file.*/
	virtual void slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts,
			const std::string& sliceFilename = "")
	{
		
		#ifdef SLICE_OUTPUT
			std::cout << "In MappedSPnode<T>::slice, I am " << getNodeName() << std::endl;
		#endif
		
		std::vector<cxsc::real> fullSlicePts = sliceCheck(sliceDims, slicePts);
		
		//start the file
		if (!sliceFilename.empty()) {
			std::ofstream os(sliceFilename.c_str());
			if (os.is_open()) os.close();
		}
		
		_slice(sliceDims, fullSlicePts, sliceFilename);
	}
	

    // nodeReabsorbChildren() can use the base class implementation
    // (the range for this will be correct so just delete the children)

	

	/*! \brief  Accept a SPExpandVisitor.
		
	If this is a leaf, and if the node is willing to split,
	accepts the visitor and sets the range equal to
	the value returned by the visitor's visit operation, then 
	recursively asks children to accept visitor too.
	
	The node's willingness to split is based on the results of the
	isSplittableNode() method.
	
	\param visitor a visitor capable of
	expanding this according to its own criteria
	and also of returning a 
	value to used as the range of this. */
	void acceptSPExpandVisitor(const SPExpandVisitor<T>& visitor)
	{
		#ifdef MYDEBUGVISITOR
			std::cout << "Using mappedspnode_sr expander accept" << std::endl;
		#endif
		
		// only accept the visit if this is a leaf and the node is splittable
		if (isLeaf() && isSplittableNode()) {
			range = visitor.visit(this);
			#ifdef MYDEBUGVISITOR
				std::cout << "after visit, range value = " << range << std::endl;
			#endif
		}
		
		
		#ifdef MYDEBUGVISITOR
			if (isLeaf() && !isSplittableNode()) {
				std::cout << "box is too small to split: volume is " << nodeRealVolume() << std::endl;
			}
		
			if(!isLeaf()) std::cout << "now visit children" << std::endl;
		#endif
		
		if (hasLCwithBox()) (getLeftChild())->acceptSPExpandVisitor(visitor);
		if (hasRCwithBox()) (getRightChild())->acceptSPExpandVisitor(visitor);
	}
	
	/*! \brief  Accept a SPValueVisitor visitor.
		
	Accepts the visitor and sets the range equal to
	the value returned by the vistor's visit operation.
	
	Recursively asks children to accept visitor too.
	
	\param visitor a visitor capable of returning a 
	value to be used as the range of this. */
	void acceptSPValueVisitor(const SPValueVisitor<T>& visitor)
	{
		#ifdef MYDEBUGVISITORE
			std::cout << "Using mappedspnode_sr valuer accept " << getNodeName() << std::endl;
		#endif
		
		range = visitor.visit(this);
		#ifdef MYDEBUGVISITORE
			std::cout << "after visit, range value = " << range << std::endl;
		#endif
		
		#ifdef MYDEBUGVISITORE
			if(!isLeaf()) std::cout << "now visit children" << std::endl;
		#endif
		
		if (hasLCwithBox()) (getLeftChild())->acceptSPValueVisitor(visitor);
		if (hasRCwithBox()) (getRightChild())->acceptSPValueVisitor(visitor);
	}
	
	/*! \brief Addition to self operator.
	\param add the object to add to this.
	\pre either both this and \a add are empty or both are
	non-empty and the boxes match.
	\post this is a non-minimal union of this and \a add
	with a range equal to the 
	sum of the values in ranges of this before the operation 
	and of \a add. If both this and \a add were empty before
	the operation, this will be empty post the operation.*/     
    MappedSPnode<T>& operator+= (const MappedSPnode<T>& add)
    {
		// if both empty, do nothing
		if (!isEmpty() || !add.isEmpty()) {
			
			// just one empty or boxes don't match
			if ( isEmpty() || add.isEmpty() ||
				(getBox() != add.getBox() ) ) {
				throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::operator+=(const MappedSPnode<T>& const)");
			}
			
			#ifdef LOGARITHMETIC
				int recordNumber = 0;
				std::string filename("MappedAddition");
				
				this->_addNonMinimalUnionAdditionSpecial(add,
										filename, recordNumber);
			#else
				this->_addNonMinimalUnion(add);
			#endif
		
		}
		return *this;
    }
    

	/*! \brief Addition operator.
    \param add the object to add to this to give the object returned.
	\return a non-minimal union of this and \a add.
	with a range equal to the 
	sum of the values in ranges of this before the operation 
	and of \a add.  If both this and \a add were empty before
	the operation, the object returned will be empty.
    \pre either both this and \a add are empty or both are
	non-empty and the boxes match.*/
	const MappedSPnode<T> operator+ (const MappedSPnode<T>& add) const
    {
		MappedSPnode<T> result =(*this);
	
		result+= add;
		
		return result;
    }
	
	/*! \brief Self-scalar addition operator.
    \param add the value to add to this.
	\post this has the same tree structure as before the operation 
	with a range equal to the range of this before 
	the operation plus \a add.  */
    MappedSPnode<T>& operator+= (const T& add)
    {
		_scalarAdd(add);
		return *this;
	}
	
	/*! \brief Scalar addition operator.
    
	\param add the value to add to this to give the object returned.
	\return a mapped paving with the same tree structure as this 
	before the operation with range equal to the range of this before 
	the operation plus \a add.  */
    const MappedSPnode<T> operator+ (const T& add) const
    {
        MappedSPnode<T> result =(*this);
	
		result+= add;
	
		return result;
    }
	

	/*! \brief Subtraction from self operator.
	\param sub the object to subtract from this.
	\pre either both this and \a sub are empty or both are
	non-empty and the boxes match.
	\post this is a non-minimal union of this and \a sub
	with a range equal to
	the range of this before
	the operation less the range of \a sub.
	If both this and \a sub were empty before
	the operation, this will be empty post the operation.*/  
    MappedSPnode<T>& operator-= (const MappedSPnode<T>& sub)
    {
		
		// if both empty, do nothing
		if (!isEmpty() || !sub.isEmpty()) {
			
			// just one empty or boxes don't match
			if ( isEmpty() || sub.isEmpty() ||
				(getBox() != sub.getBox() ) ) {
				throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::operator-=(const MappedSPnode<T>& const)");
			}
			
			this->_subtractNonMinimalUnion(sub);
		
		}
		return *this;
	}
    
	/*! \brief Subtraction operator.
	\param sub the object to subtract from this to give the object returned.
	\return a non-minimal union of this and \a sub
	with a range equal to the 
	value of the range of this before
	the operation less the range of of \a sub.  
	If both this and \a sub were empty before
	the operation, the returned object will be empty post the operation.
	\pre either both this and \a sub are empty or both are
	non-empty and the boxes match.*/   
    const MappedSPnode<T> operator- (const MappedSPnode<T>& sub) const
    {
        MappedSPnode<T> result =(*this);
	
		result-= sub;
	
		return result;
    }

	/*! \brief Self-scalar subtraction operator.
    \param sub the value to subtract from this.
	\post this has the same tree structure as before the operation 
	with a range equal to the range of this before 
	the operation minus \a sub.  */
    MappedSPnode<T>& operator-= (const T& sub)
    {
		_scalarSubtract(sub);
		return *this;
	}
	
	/*! \brief Scalar subtraction operator.
    
	\param sub the value to subtract from this to give the object returned.
	\return a mapped paving with the same tree structure as this 
	before the operation with range equal to the range of this before 
	the operation minus \a sub.  */
    const MappedSPnode<T> operator- (const T& sub) const
    {
        MappedSPnode<T> result =(*this);
	
		result-= sub;
	
		return result;
    }
	
	/*! \brief Multiplication of self operator.
	\param mult the object to multiply this by.
	\pre either both this and \a mult are empty or both are
	non-empty and the boxes match.
	\post this is a non-minimal union of this and \a mult
	with a range equal to
	the range of this before
	the operation multiplied by the range
	of \a mult.  If both this and \a mult were empty before
	the operation, this will be empty post the operation.*/
    MappedSPnode<T>& operator*= (const MappedSPnode<T>& mult)
    {
		
		// if both empty, do nothing
		if (!isEmpty() || !mult.isEmpty()) {
			
			// just one empty or boxes don't match
			if ( isEmpty() || mult.isEmpty() ||
				(getBox() != mult.getBox() ) ) {
				throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::operator*=(const MappedSPnode<T>& const)");
			}
			// both must be non-empty
			
			this->_multiplyNonMinimalUnion(mult);

			
		}
		
		return *this;
	}
	
    /*! \brief Multiplication operator.
    \param mult the object to multiply this by to give the object returned.
	\return a non-minimal union of this and \a mult
	with a range equal to
	the range of this before
	the operation multiplied by the range of \a mult.
	If both this and \a mult were empty before
	the operation, the returned object will be empty post the operation.
	\pre either both this and \a mult are empty or both are
	non-empty and the boxes match.*/   
    const MappedSPnode<T> operator* (const MappedSPnode<T>& mult) const
    {
        MappedSPnode<T> result =(*this);
	
		result*= mult;
		
		return result;
    }
	
	/*! \brief Self-scalar multiplication operator.
    \param mult the value to multiply this by.
	\post this has the same tree structure as before the operation 
	with a range equal to the range of this before 
	the operation multiplied by \a mult.  */
    MappedSPnode<T>& operator*= (const T& mult)
    {
		_scalarMult(mult);
		return *this;
	}
	
	/*! \brief Scalar multiplication operator.
    
	\param mult the value to multiply this by to give the object returned.
	\return a mapped paving with the same tree structure as this 
	before the operation with range equal to the range of this before 
	the operation multiplied by \a mult.  */
    const MappedSPnode<T> operator* (const T& mult) const
    {
        MappedSPnode<T> result =(*this);
	
		result*= mult;
	
		return result;
    }
	
	
	/*! \brief Division of self operator.
	\param div the object to divide this by.
	\pre either both this and \a div are empty or both are
	non-empty and the boxes match.
	\post this is a non-minimal union of this and \a div
	with a range equal to
	the range of this before
	the operation divided by the range of
	of \a div.  If both this and \a div were empty before
	the operation, this will be empty post the operation.*/
    MappedSPnode<T>& operator/= (const MappedSPnode<T>& div)
    {
		// if both empty, do nothing
		if (!isEmpty() || !div.isEmpty()) {
			
			// just one empty or boxes don't match
			if ( isEmpty() || div.isEmpty()  ||
				(getBox() != div.getBox() ) ) {
				throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::operator/=(const MappedSPnode<T>& const)");
			}
			// both must be non-empty
			
			this->_divideNonMinimalUnion(div);

			
		}
		return *this;
	}
	
	/*! \brief Division operator.
    \param div the object to divide this by to give the object returned.
	\return a non-minimal union of this and \a div
	with a range equal to
	the range of this before
	the operation divided by the range
	of \a div.  If both this and \a div were empty before
	the operation, the returned object will be empty post the operation.
	\pre either both this and \a div are empty or both are
	non-empty and the boxes match.*/   
    const MappedSPnode<T> operator/ (const MappedSPnode<T>& div) const
    {
        MappedSPnode<T> result =(*this);
	
		result/= div;
	
		return result;
    }
	
	/*! \brief Self-scalar division operator.
    \param div the value to divide this by.
	\post this has the same tree structure as before the operation.
	with a range equal to the range of this before 
	the operation divided by \a div.   */
    MappedSPnode<T>& operator/= (const T& div)
    {
		_scalarDiv(div);
		return *this;
	}
	
	/*! \brief Scalar division operator.
    \param div the value to divide this by to give the object returned.
	\return an object with the same tree structure as this 
	before the operation with a range 
	equal to the range of this before 
	the operation divided by \a div.  */
    const MappedSPnode<T> operator/ (const T& div) const
    {
        MappedSPnode<T> result =(*this);
	
		result/= div;
	
		return result;
    }
	
	/*! \brief Change this so that it
	has the minimum number of leaves necessary to represent
	the same overall 'shape' (combination of subpaving volumes
	and values) as before, by recombining any pair of sibling
	leaves with identical ranges so that their 
	former parent becomes a leaf with the same range
	as formerly belonged to each child. 
	
	Eg a node with two leafs where both children have the same value
	\f$ r \f$ mapped to them can be represented as a single leaf 
	with the value \f$ r \f$ mapped to it.  minimiseLeaves can 
	be used after arithmetic operations on nodes
	to deal with cases where sibling leaf nodes in the result
	are equivalent in terms of their range values.
	
	The minimisation
	process works from the leaves upwards, so that a tree could
	potentially be reduced to a single node if every pair 
	of sibling leaves encountered in the process can be 
	recombined. */
	void minimiseLeaves()
	{
		if (!isLeaf()) { // can't do anything if this is a leaf
			
			// recurse first
			if ( hasLCwithBox() ) getLeftChild()->minimiseLeaves();
			if ( hasRCwithBox() ) getRightChild()->minimiseLeaves();
			
			// now do me
			if ( isSubLeaf() && 
				(getLeftChild()->getRange() == getRightChild()->getRange()) ) {
				
					// make this range collection one of the childrens
					range = getLeftChild()->getRange();
					delete leftChild;
					leftChild = NULL;
					delete rightChild;
					rightChild = NULL;
			}
		
		}
	}
	
	/*! \brief Swap the properties of this and another.
	 * 
	The properties swapped include all the node properties like name,
	box, etc, and children, and the range.
	
    \param spn the node to swap with.
	\post this has the properties that \a spn had before the operation,
	and spn has the properties this had before the operation.  */
    void swapMSPSR(MappedSPnode& spn) //throw() // don't hide base class version
	{
		/* theBox, parent, leftChild,
		rightChild and nodeName are inherited from base class */
		SPnode::swap(spn); // use the base version
		
		std::swap(range, spn.range);
	}
	
	/*! \brief Print the details of a specific node in a subpaving.
	
	Prints:
	nodeName [tab] "Box is" [tab]
	... each interval in the box in cxsc output format ...
	[tab] "Box volume is" [tab] volume
	[tab] "range data is" [tab] range
	[newline]
	* 
	\param os the stream to print to.*/
    virtual std::ostream& nodePrint(std::ostream &os) const
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
	
	


    protected:
    // -------------------------- protected member functions -------------
    
    /*! \brief Print the details of a single leaf node, using tab delimiters.*/
    virtual std::ostream& leafOutputTabs(std::ostream &os) const
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

		
	/*! \brief Reshape this and \rhs so that 
	both this and rhs have identical shapes after the operation 
	and that shape is the non-minimal union of shape the of this
	shape and the shape of rhs before the operation.
    */
 	void _reshapeTreesToUnion(MappedSPnode<T> * const rhs)
	{
		if (rhs == NULL) {

			throw NullSubpavingPointer_Error(
			"MappedSPnode<T>::_reshapeTreesToUnion(MappedSPnode<T> * const)");
		}

		// this is not a leaf, rhs is a leaf
		if (!isLeaf() && rhs->isLeaf()) {

			//we need to expand rhs
			rhs->nodeExpand();
		}

		// this is a leaf, rhs is not a leaf
		if (isLeaf() && !rhs->isLeaf()) {

			//we need to expand this
			nodeExpand();
		}

		// we have made sure that if at least one of them was not a leaf,
		// then neither are now leaves
		// now recurse on the children
		if (!isLeaf() && !rhs->isLeaf()) {
			getLeftChild()->_reshapeTreesToUnion(rhs->getLeftChild());
			getRightChild()->_reshapeTreesToUnion(rhs->getRightChild());

		}
	}
	
	#ifdef LOGARITHMETIC
		// only used for generating log file output from addition 
		int _reshapeTreesToUnionAdditionSpecial(
											MappedSPnode<T> * const rhs,
											const std::string& filename,
											int recordNumber)
		{
			if (rhs == NULL) {

				throw NullSubpavingPointer_Error(
				"MappedSPnode<T>::_reshapeTreesToUnion(MappedSPnode<T> * const)");
			}
			
			// this is not a leaf, rhs is a leaf
			if (!isLeaf() && rhs->isLeaf()) {

				//we need to expand rhs
				rhs->nodeExpand();
			}

			// this is a leaf, rhs is not a leaf
			if (isLeaf() && !rhs->isLeaf()) {
				//we need to expand this
				nodeExpand();
			}
			// add the ranges only after expanding
			range = range + rhs->range;
			
			//log it
			if (isLeaf()) {
				logFromRoot(filename, recordNumber);
				++recordNumber;
			}
			
			// we have made sure that if at least one of them was not a leaf,
			// then neither are now leaves
			// now recurse on the children
			if (!isLeaf() && !rhs->isLeaf()) {
				recordNumber = 
					getLeftChild()->_reshapeTreesToUnionAdditionSpecial(
						rhs->getLeftChild(),
						filename,
						recordNumber);
				recordNumber = 
					getRightChild()->_reshapeTreesToUnionAdditionSpecial(
						rhs->getRightChild(),
						filename,
						recordNumber);

			}
			return recordNumber;
		}
		
		//log it
		void logFromRoot(const std::string& filename,
						int recordNumber, int prec = 5)
		{
			if (NULL == parent) {
				std::ostringstream oss;
				oss << filename << recordNumber << ".txt";
				std::ofstream os;
				os.open((oss.str()).c_str());
				if (os.is_open()) {
					leavesOutputTabs(os, prec);
					os.close();
				}
				else {
					std::cerr << "Error: could not open log file"
							<< std::endl << std::endl;
				}
			}
			else getParent()->logFromRoot(filename, recordNumber, prec);
						
		}
			
	#endif
	
	/*! \brief Give this a range equal to range of this and range of other.
	
    Assumes that this and other are the same from leaves up to this and other,
    i.e. are roots of trees of the same shape up to the depth of this
    Note that other could be a taller trees than this: range collections
    are only combined up to the level of this. */
    virtual void _addRanges(	const MappedSPnode<T> * const other)
	{
		if (other == NULL) {
			throw NullSubpavingPointer_Error(
				"MappedSPnode<T>::_addRanges(const MappedSPnode<T> * const)");
		}

		// recurse on the children if any first
		if (!isLeaf()) {
			getLeftChild()->_addRanges(
													other->getLeftChild());
			getRightChild()->_addRanges(
													other->getRightChild());

		}

		// deal with this range collection
		range = range + other->range;
	}
	
	/*! \brief Give this a range equal to range of this - range of other.
	
    Assumes that this and other are the same from leaves up to this and other,
    i.e. are roots of trees of the same shape up to the depth of this
    Note that other could be a taller trees than this: range collections
    are only combined up to the level of this. */
    virtual void _subtractRanges(	const MappedSPnode<T> * const other)
	{
		if (other == NULL) {
			throw NullSubpavingPointer_Error(
				"MappedSPnode<T>::_subtractRanges(const MappedSPnode<T> * const)");
		}

		// recurse on the children if any first
		if (!isLeaf()) {
			getLeftChild()->_subtractRanges(
													other->getLeftChild());
			getRightChild()->_subtractRanges(
													other->getRightChild());

		}

		// deal with this range collection
		range = range - other->range;
	}
	
	/*! \brief Give this a range equal to range of this * range of other.
	
    Assumes that this and other are the same from leaves up to this and other,
    i.e. are roots of trees of the same shape up to the depth of this
    Note that other could be a taller trees than this: range collections
    are only combined up to the level of this. */
    virtual void _multRanges(	const MappedSPnode<T> * const other)
	{
		if (other == NULL) {
			throw NullSubpavingPointer_Error(
				"MappedSPnode<T>::_multRanges(const MappedSPnode<T> * const)");
		}

		// recurse on the children if any first
		if (!isLeaf()) {
			getLeftChild()->_multRanges(
													other->getLeftChild());
			getRightChild()->_multRanges(
													other->getRightChild());

		}

		// deal with this range collection
		range = range * other->range;
	}
	
	/*! \brief Give this a range equal to range of this divided by range of other.
	
    Assumes that this and other are the same from leaves up to this and other,
    i.e. are roots of trees of the same shape up to the depth of this
    Note that other could be a taller trees than this: range collections
    are only combined up to the level of this. */
    virtual void _divRanges(	const MappedSPnode<T> * const other)
	{
		if (other == NULL) {
			throw NullSubpavingPointer_Error(
				"MappedSPnode<T>::_divRanges(const MappedSPnode<T> * const)");
		}

		// recurse on the children if any first
		if (!isLeaf()) {
			getLeftChild()->_divRanges(
													other->getLeftChild());
			getRightChild()->_divRanges(
													other->getRightChild());

		}

		// deal with this range collection
		range = range / other->range;
	}
	
	


	/*! \brief Make a non-minimal union of subpavings using
	  addition of ranges.
   */
   virtual void _addNonMinimalUnion(
                   const MappedSPnode<T>& rhs)
    {
        
		//if (rhs == NULL || rhs->isEmpty()) do nothing
		
		// should not need to check on rhs here 
		if ((isEmpty()) && !rhs.isEmpty()) {

			*this = MappedSPnode<T>(rhs);
		}
				
		if (!isEmpty() && !rhs.isEmpty()) {

			// make copies of rhs to reshape

			MappedSPnode<T> rhsTemp(rhs);

			// reshape
			_reshapeTreesToUnion(&rhsTemp);

			_addRanges(&rhsTemp);

		}
	}
	
	// only used for generating log file output from addition 
	#ifdef LOGARITHMETIC
		virtual void _addNonMinimalUnionAdditionSpecial(
					   const MappedSPnode<T>& rhs,
					   const std::string& filename,
					   int recordNumber)
		{
			
			//if (rhs == NULL || rhs->isEmpty()) do nothing
			
			// should not need to check on rhs here 
			if ((isEmpty()) && !rhs.isEmpty()) {

				*this = MappedSPnode<T>(rhs);
			}
					
			if (!isEmpty() && !rhs.isEmpty()) {

				// make copies of rhs to reshape

				MappedSPnode<T> rhsTemp(rhs);

				// do the whole lot
				_reshapeTreesToUnionAdditionSpecial(&rhsTemp,
									filename, recordNumber);

			}
		}
	#endif
	
	/*! \brief Make a non-minimal union of subpavings using
	  subtraction of ranges.
    
    \pre this and *rhs are non-empty.*/
    virtual void _subtractNonMinimalUnion(
                   const MappedSPnode<T>& rhs)
    {
        // should not need to check on rhs not null here 
		if ((isEmpty()) || rhs.isEmpty()) {

			throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::_subtractNonMinimalUnion(const MappedSPnode<T> * const)");
		}
		
		// make copies of rhs to reshape

		MappedSPnode<T> rhsTemp(rhs);

		// reshape
		_reshapeTreesToUnion(&rhsTemp);

		_subtractRanges(&rhsTemp);

	}

	/*! \brief Make a non-minimal union of subpavings using
	multiplication of ranges.
    
    \pre this and *rhs are non-empty.*/
    virtual void _multiplyNonMinimalUnion(
                   const MappedSPnode<T>& rhs)
    {
		// should not need to check on rhs not null here 
		if ((isEmpty()) || rhs.isEmpty()) {

			throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::_multiplyNonMinimalUnion(const MappedSPnode<T> * const)");
		}
		
        // make copies of rhs to reshape

		MappedSPnode<T> rhsTemp(rhs);

		// reshape
		_reshapeTreesToUnion(&rhsTemp);

		_multRanges(&rhsTemp);

	}

	/*! \brief Make a non-minimal union of subpavings using
	division of ranges.
    
    \pre this and *rhs are non-empty.*/
	virtual void _divideNonMinimalUnion(
                   const MappedSPnode<T>& rhs)
    {
		if ((isEmpty()) || rhs.isEmpty()) {

			throw IncompatibleDimensions_Error(
				"MappedSPnode<T>::_divideNonMinimalUnion(const MappedSPnode<T> * const)");
		}
		
        // make copies of rhs to reshape

		MappedSPnode<T> rhsTemp(rhs);

		// reshape
		_reshapeTreesToUnion(&rhsTemp);

		_divRanges(&rhsTemp);

	}
	
	/*! \brief Increment range collection of this.
    
    Increment range by type T scalar. 
    \param add the scalar increment, of same type as range, to use.
    */
    virtual void _scalarAdd(const T& add)
    {
		
		if (!isEmpty()) {
			range = range + add;
			// recurse on children
			if (hasLCwithBox())
				getLeftChild()->_scalarAdd(add);
			if (hasRCwithBox())
				getRightChild()->_scalarAdd(add);
		}
    }
	
	/*! \brief Increment range collection of this.
    
    Increment range by type T scalar. 
    \param sub the scalar decrement, of same type as range, to use.
    */
    virtual void _scalarSubtract(const T& sub)
    {
		
		if (!isEmpty()) {
			range = range - sub;
			// recurse on children
			if (hasLCwithBox())
				getLeftChild()->_scalarSubtract(sub);
			if (hasRCwithBox())
				getRightChild()->_scalarSubtract(sub);
		}
    }
	
	
	/*! \brief Scale up range of this.
    
    Multiplication this range by type T scalar. 
    \param mult the scalar multiplier, of same type as range, to use.
    */
    virtual void _scalarMult(const T& mult)
    {
		
		if (!isEmpty()) {
			range = range * mult;
			// recurse on children
			if (hasLCwithBox())
				getLeftChild()->_scalarMult(mult);
			if (hasRCwithBox())
				getRightChild()->_scalarMult(mult);
		}
    }

	/*! \brief Scale down range of this.
    
    Division this range by type T scalar.
    \param div the scalar divisor, of same type as range, to use.
    */
    virtual void _scalarDiv(const T& div)
    {
		
		if (!isEmpty()) {
			range = range / div;
			// recurse on children
			if (hasLCwithBox())
				getLeftChild()->_scalarDiv(div);
			if (hasRCwithBox())
				getRightChild()->_scalarDiv(div);
		}
    }

	// assumes the slice pts are inside the box of this node somewhere
	virtual void _slice(
				const std::vector < int >& sliceDims,
				const std::vector < cxsc::real >& slicePts,
				unsigned randID = 0)
	{
		#ifdef SLICE_OUTPUT
			std::cout << "In MappedSPnodeSR<T>::slice, I am " << getNodeName() << std::endl;
		#endif
		
		#ifdef CAPTURE_SLICE_BOX
			if (randID == 0) { 
				for (time_t t = time(NULL) + 1; time(NULL) < t; ) {}
				gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937);
				long seed = time (NULL) * getpid();    
				gsl_rng_set (rgsl, seed);
				randID = gsl_rng_get (rgsl);
				gsl_rng_free(rgsl);
				std::cout << "Slicing:  sliceBoxesFilename = sliceBoxesFilename" 
					<< randID << "...\n" << std::endl;
			}
			std::ostringstream oss;
			
			oss << "sliceBoxes" << randID;
			for (size_t i = 0; i < sliceDims.size(); ++i) {
					oss << "_" << _double(slicePts[sliceDims[i]-1]);
			}
			oss << ".txt";
			std::string sliceBoxesFilename = oss.str();
				
		
			if (NULL == parent) { // start outputfile
				std::ofstream os(sliceBoxesFilename.c_str());
				if (os.is_open()) os.close();
			}
		#endif
		
		bool splitOnSliceDim = false;
		
		int splitDim = -1;
		
		if (!isLeaf()) {
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tI am not a leaf: get split dim" << std::endl;
			#endif
			
			/* Check whether we split on a split dim
			  before we do anything to the children*/
			splitDim = getSplitDim();
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tsplitDim = " << splitDim << std::endl;
			#endif
			
			std::vector< int >::const_iterator found 
				= find (sliceDims.begin(), sliceDims.end(), splitDim);
		
			// can we find splitDim in the sliceDims?
			splitOnSliceDim = (found < sliceDims.end()); 
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tsplitOnSliceDim = " << splitOnSliceDim << std::endl;
			#endif 
		}
		
		/* if I am a leaf, or did not split on a slice dimension,
		 * adjust the box */ 
		if (!splitOnSliceDim) {
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tI did not split on a dimension to take out, or I am a leaf: first deal with children ..." << std::endl;
			#endif
			
			if (!isLeaf()) {
				// Slice on the children first
				getLeftChild()->_slice(sliceDims, slicePts, randID);
				getRightChild()->_slice(sliceDims, slicePts, randID);
			}
			
			#ifdef CAPTURE_SLICE_BOX
				if (isLeaf()) {
					std::ofstream os(sliceBoxesFilename.c_str(), std::ios::app);
					if (os.is_open()) {
						leafOutputTabs(os);
						os << "\n";
						os.close();
					}
				}
			#endif
			
			// now deal with this node itself
			#ifdef SLICE_OUTPUT
				if (!isLeaf()) {
					std::cout << "\nback in " << getNodeName() << std::endl;
					std::cout << "\tI did not split on a dimension to take out, or I am a leaf, so need to contract box" << std::endl;
				}
				std::cout << "\tmy split dimension is " << getSplitDim() << std::endl;
				std::cout << "\tmy box is " << getBox() << std::endl;
				std::cout << "\tmy range is " << range << std::endl;
			#endif
			
			ivector box = getBox();
			int dim = VecLen(box);
			int boxLB = Lb(box);
		
			int newDims = dim - sliceDims.size();
			ivector newBox = ivector(newDims); 
			int index = Lb(newBox);
			int oldindex = boxLB;
		
			// put in the upper and lower bounds for the new box
			// for each dimension that stays	
			for (; oldindex <= Ub(box); oldindex++) {
				std::vector<int>::const_iterator fit 
				= find (sliceDims.begin(), sliceDims.end(), (oldindex - boxLB + 1));
				if (!(fit < sliceDims.end())) { // keep this one
					newBox[index] = box[oldindex];
					index++;
				}
			}
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tnew box is " << newBox << std::endl;
			#endif
			
			//store the child node locations and then temporarily detach them
			MappedSPnode <T>* savedLC = getLeftChild();
			MappedSPnode <T>* savedRC = getRightChild();
			leftChild = NULL;
			rightChild = NULL;
			
			// also need to store parent
			MappedSPnode <T>* savedParent = getParent();
			parent = NULL;
			
			T temp = getRange();
			
			// replace contents of this with contents of a newly made node	
			MappedSPnode <T> tempNode(newBox, temp);
			this->swapMSPSR(tempNode);
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tafter contracting, my box is " << getBox() << std::endl;
				std::cout << "\n\tmy range is " << range << std::endl;
			#endif
			
			// put the child pointers back, and reattach to parent
			leftChild = savedLC;
			rightChild = savedRC;
			parent = savedParent;
			
		}
		else {
		
			/* if we split on any of the slice dimensions
				 then we drop this from the tree entirely and
				 replace it with the child who actually contains the slice */
				
			/*to do this, we just have to know where we split and look at the slice
			 * points on the same dimension and see if they are above or below our
			 * split. We know we are not a leaf if we are in here at all.  */
			#ifdef SLICE_OUTPUT
				std::cout << "I DID split on a dimension to take out (so I can't be a leaf)" << std::endl;
			#endif 
			
			real splitValue = Sup(getLeftChild()->getBox()[splitDim]);	
			
			real slPt = slicePts.at(splitDim-1); // vectors indexed 0 - n-1
			
			#ifdef SLICE_OUTPUT
				std::cout << "split dimension is " << splitDim << " and split value is " << splitValue << std::endl;
				std::cout << "slPt = slicePts.at(splitDim-1) is " << slPt << std::endl;
			#endif
			
			MappedSPnode <T> temp;
			
			// save parent
			MappedSPnode <T>* savedParent = getParent();
			// detach from parent temporarily
			parent = NULL;
				
			if (slPt < splitValue) { // crucial slice point is in left child
			
				#ifdef SLICE_OUTPUT
					std::cout << "slPt is in left child" << getBox() << std::endl;
					std::cout << "need to do getLeftChild()->_slice(sliceDims, slicePts) and then copy it " << std::endl;
				#endif
				
				// Slice on the left child
				getLeftChild()->_slice(sliceDims, slicePts,randID);
				
				// copy left child
				temp = *getLeftChild();
			}
			else { // crucial slice point is in right child
			
				#ifdef SLICE_OUTPUT
					std::cout << "slPt is in right child" << getBox() << std::endl;
					std::cout << "need to do getRightChild()->_slice(sliceDims, slicePts) and then copy it " << std::endl;
				#endif
				
				// Slice on the right child
				getRightChild()->_slice(sliceDims, slicePts, randID);
				
				// copy left child
				temp = *getRightChild();
			}
			
			
			
			#ifdef SLICE_OUTPUT
				std::cout << "\nback in " << getNodeName() << std::endl;
				std::cout << "\tNow copy the temp into me" << std::endl;
				std::cout << "\tmy split dimension is " << getSplitDim() << std::endl;
				std::cout << "\tmy box is " << getBox() << std::endl;
				std::cout << "\tmy range is " << range << std::endl;
			#endif
			
			swapMSPSR(temp); //swap me and temp
			
				
			#ifdef SLICE_OUTPUT
				std::cout << "\tre-made me out of copy of child that is sliced:" << std::endl;
				std::cout << "\tmy name (before renaming) was " << temp.getNodeName() << std::endl;
				std::cout << "\tmy range is now " << range << std::endl;
				std::cout << "\tand my box is " << getBox() << std::endl;
				if (!isLeaf() ) {
					std::cout << "\tand my children (before renaming) are:" << std::endl;
					getLeftChild()->oneLineOutput(std::cout, 2);
					getRightChild()->oneLineOutput(std::cout, 2);
				}
			#endif
			
			// restore relationship to parent
			parent = savedParent;
			
				
		}// end else
		
		//if we are the root, recursively rename everything
		if (getParent() == NULL) {
			setNodeName("X"); // my name might be wrong
			recursiveRename();
			#ifdef MARG_OUTPUT
				std::cout << "\nNow recursively rename everything from me down\n\n" << std::endl;
			#endif
		}
	}
	
		// assumes the slice pts are inside the box of this node somewhere
	virtual void _slice(
				const std::vector < int >& sliceDims,
				const std::vector < cxsc::real >& slicePts,
				const string& sliceFilename)
	{
		#ifdef SLICE_OUTPUT
			std::cout << "In MappedSPnodeSR<T>::slice, I am " << getNodeName() << std::endl;
		#endif
		
		bool splitOnSliceDim = false;
		
		int splitDim = -1;
		
		if (!isLeaf()) {
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tI am not a leaf: get split dim" << std::endl;
			#endif
			
			/* Check whether we split on a split dim
			  before we do anything to the children*/
			splitDim = getSplitDim();
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tsplitDim = " << splitDim << std::endl;
			#endif
			
			std::vector< int >::const_iterator found 
				= find (sliceDims.begin(), sliceDims.end(), splitDim);
		
			// can we find splitDim in the sliceDims?
			splitOnSliceDim = (found < sliceDims.end()); 
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tsplitOnSliceDim = " << splitOnSliceDim << std::endl;
			#endif 
		}
		
		/* if I am a leaf, or did not split on a slice dimension,
		 * adjust the box */ 
		if (!splitOnSliceDim) {
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tI did not split on a dimension to take out, or I am a leaf: first deal with children ..." << std::endl;
			#endif
			
			if (!isLeaf()) {
				// Slice on the children first
				getLeftChild()->_slice(sliceDims, slicePts, sliceFilename);
				getRightChild()->_slice(sliceDims, slicePts, sliceFilename);
			}
			else if ( !sliceFilename.empty()) {
				
				std::ofstream os(sliceFilename.c_str(), std::ios::app);
				if (os.is_open()) {
					leafOutputTabs(os);
					os << "\n";
					os.close();
				}
			}
						
			// now deal with this node itself
			#ifdef SLICE_OUTPUT
				if (!isLeaf()) {
					std::cout << "\nback in " << getNodeName() << std::endl;
					std::cout << "\tI did not split on a dimension to take out, or I am a leaf, so need to contract box" << std::endl;
				}
				std::cout << "\tmy split dimension is " << getSplitDim() << std::endl;
				std::cout << "\tmy box is " << getBox() << std::endl;
				std::cout << "\tmy range is " << range << std::endl;
			#endif
			
			ivector box = getBox();
			int dim = VecLen(box);
			int boxLB = Lb(box);
		
			int newDims = dim - sliceDims.size();
			ivector newBox = ivector(newDims); 
			int index = Lb(newBox);
			int oldindex = boxLB;
		
			// put in the upper and lower bounds for the new box
			// for each dimension that stays	
			for (; oldindex <= Ub(box); oldindex++) {
				std::vector<int>::const_iterator fit 
				= find (sliceDims.begin(), sliceDims.end(), (oldindex - boxLB + 1));
				if (!(fit < sliceDims.end())) { // keep this one
					newBox[index] = box[oldindex];
					index++;
				}
			}
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tnew box is " << newBox << std::endl;
			#endif
			
			//store the child node locations and then temporarily detach them
			MappedSPnode <T>* savedLC = getLeftChild();
			MappedSPnode <T>* savedRC = getRightChild();
			leftChild = NULL;
			rightChild = NULL;
			
			// also need to store parent
			MappedSPnode <T>* savedParent = getParent();
			parent = NULL;
			
			T temp = getRange();
			
			// replace contents of this with contents of a newly made node	
			MappedSPnode <T> tempNode(newBox, temp);
			this->swapMSPSR(tempNode);
			
			#ifdef SLICE_OUTPUT
				std::cout << "\tafter contracting, my box is " << getBox() << std::endl;
				std::cout << "\n\tmy range is " << range << std::endl;
			#endif
			
			// put the child pointers back, and reattach to parent
			leftChild = savedLC;
			rightChild = savedRC;
			parent = savedParent;
			
		}
		else {
		
			/* if we split on any of the slice dimensions
				 then we drop this from the tree entirely and
				 replace it with the child who actually contains the slice */
				
			/*to do this, we just have to know where we split and look at the slice
			 * points on the same dimension and see if they are above or below our
			 * split. We know we are not a leaf if we are in here at all.  */
			#ifdef SLICE_OUTPUT
				std::cout << "I DID split on a dimension to take out (so I can't be a leaf)" << std::endl;
			#endif 
			
			real splitValue = Sup(getLeftChild()->getBox()[splitDim]);	
			
			real slPt = slicePts.at(splitDim-1); // vectors indexed 0 - n-1
			
			#ifdef SLICE_OUTPUT
				std::cout << "split dimension is " << splitDim << " and split value is " << splitValue << std::endl;
				std::cout << "slPt = slicePts.at(splitDim-1) is " << slPt << std::endl;
			#endif
			
			MappedSPnode <T> temp;
			
			// save parent
			MappedSPnode <T>* savedParent = getParent();
			// detach from parent temporarily
			parent = NULL;
				
			if (slPt < splitValue) { // crucial slice point is in left child
			
				#ifdef SLICE_OUTPUT
					std::cout << "slPt is in left child" << getBox() << std::endl;
					std::cout << "need to do getLeftChild()->_slice(sliceDims, slicePts) and then copy it " << std::endl;
				#endif
				
				// Slice on the left child
				getLeftChild()->_slice(sliceDims, slicePts, sliceFilename);
				
				// copy left child
				temp = *getLeftChild();
			}
			else { // crucial slice point is in right child
			
				#ifdef SLICE_OUTPUT
					std::cout << "slPt is in right child" << getBox() << std::endl;
					std::cout << "need to do getRightChild()->_slice(sliceDims, slicePts) and then copy it " << std::endl;
				#endif
				
				// Slice on the right child
				getRightChild()->_slice(sliceDims, slicePts, sliceFilename);
				
				// copy left child
				temp = *getRightChild();
			}
			
			
			
			#ifdef SLICE_OUTPUT
				std::cout << "\nback in " << getNodeName() << std::endl;
				std::cout << "\tNow copy the temp into me" << std::endl;
				std::cout << "\tmy split dimension is " << getSplitDim() << std::endl;
				std::cout << "\tmy box is " << getBox() << std::endl;
				std::cout << "\tmy range is " << range << std::endl;
			#endif
			
			swapMSPSR(temp); //swap me and temp
			
				
			#ifdef SLICE_OUTPUT
				std::cout << "\tre-made me out of copy of child that is sliced:" << std::endl;
				std::cout << "\tmy name (before renaming) was " << temp.getNodeName() << std::endl;
				std::cout << "\tmy range is now " << range << std::endl;
				std::cout << "\tand my box is " << getBox() << std::endl;
				if (!isLeaf() ) {
					std::cout << "\tand my children (before renaming) are:" << std::endl;
					getLeftChild()->oneLineOutput(std::cout, 2);
					getRightChild()->oneLineOutput(std::cout, 2);
				}
			#endif
			
			// restore relationship to parent
			parent = savedParent;
			
				
		}// end else
		
		//if we are the root, recursively rename everything
		if (getParent() == NULL) {
			setNodeName("X"); // my name might be wrong
			recursiveRename();
			#ifdef MARG_OUTPUT
				std::cout << "\nNow recursively rename everything from me down\n\n" << std::endl;
			#endif
		}
	}
	
	/*! \brief Check slice parameters and return 
	a full vector of slice points.
	* 
	* Checks that:
	* <ul>
	* <li>there is a box,</li>
	* <li>that dimensions of box are compatible
	* with sliceDims given,</li>
	* <li>that slicePts is of same size as sliceDims,</li>
	* <li>that both sliceDims and slicePts are non-empty,</li>
	* <li>that sliceDims does not contain every dimension in this,</li>
	* <li>that for the ith dimension in sliceDims, the i-th value in slicePts
	* is inside the interval of the box on the ith dimension in sliceDims.</li>
	* </ul>
	* 
	* Throws exceptions if any problems are encountered.
	* 
	* \param sliceDims is a vector of dimensions to slice on.
	* \param slicePts is a vector of points to slice on, assumed to 
	* correspond to the dimensions in \a sliceDims, ie the ith value
	* in \a sliceDims gives the dimension for the ith point in
	* \a slicePts.
	* \return a vector of slice points with as many elements in as 
	* the dimensions of this, with the values from \a slicePts 
	* as the positions given by \a sliceDims and 0's (dummy values)
	* otherwise.	*/
	virtual std::vector < cxsc::real > sliceCheck(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts) const
	{
		std::vector < cxsc::real > fullSlicePts;
		
		if ( getParent() != NULL ) {
			throw NonRootNode_Error(
				"sliceCheck error found");
		}
		if ( !sliceDims.empty() || !slicePts.empty()) {
			
		
			if ( sliceDims.size() != slicePts.size()) {
				throw std::invalid_argument(
					"sliceCheck error found: sliceDims.size() != slicePts.size()");
			}
			
			/* get ordered unique values in sliceDims*/
			std::set<int> dimsSet(sliceDims.begin(), sliceDims.end());
			
			if ( dimsSet.size() < sliceDims.size()) {
				throw std::invalid_argument(
					"sliceCheck error found: duplicate dimensions in sliceDims");
			}
			
			if ( dimsSet.size() >= getDimension()) {
				throw std::invalid_argument(
					"sliceCheck error found: number of dimensions to slice on >= dimensions of this");
			}
			
			
			/* check that each of the values in sliceDims is a dimension of this*/
			ivector box = getBox();
			int lb = Lb(box);
			int ub = Ub(box);  // ub-lb+1 = dim of box
			
			// make a full vector of dim reals, 0.0's in all positions
			fullSlicePts =  std::vector< cxsc::real >(ub-lb+1, real(0.0));
			
			#ifdef SLICE_OUTPUT
				std::cout << "\nin slice" << getNodeName() << std::endl;
				std::cout << "lb = " << lb << " and ub = " << ub << std::endl;
			#endif
			
			std::vector < real >::const_iterator rit = slicePts.begin();
			for ( std::vector < int >::const_iterator cit = sliceDims.begin();
					cit < sliceDims.end();
					++cit, ++rit) {
				int d = *cit;
				#ifdef SLICE_OUTPUT
					std::cout << "checking d =" << d << std::endl;
					std::cout << "*rit =" << (*rit) << std::endl;
				#endif
					
				if (d < lb || d > ub) {
					throw std::invalid_argument(
					"sliceCheck error found: illegal dimension in sliceDims");
				}
				#ifdef SLICE_OUTPUT
					std::cout << "box[d] =" << box[d] << std::endl;
					std::cout << "(*rit) <= box[d] =" << ((*rit) <= box[d]) << std::endl;
				#endif
				real r = *rit;
				if (!(r <= box[d])) {
					throw std::invalid_argument(
					"sliceCheck error found: illegal pt in slicePts");
				}
				#ifdef SLICE_OUTPUT
					std::cout << "fullSlicePts[" << (d-lb) << "] = " << r << std::endl;
				#endif
				
				fullSlicePts[d-lb] = r;
			}
			
		} // end check at least one is non-empty
		else  { // both empty
			throw std::invalid_argument(
				"sliceCheck error found: sliceDims and slicePts both empty");
		}
		return fullSlicePts;
	}
	

	virtual std::ostream& oneLineOutput(std::ostream& os, int level) const
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

	// data members
	
	/*! \brief A range of type T.
    */
    T range;



	
	private:

    /* theBox, parent, leftChild,
    rightChild and nodeName are inherited from base class.
    */

    


 }; // end MappedSPnode<T> class



    // ----------------- non member tools functions ----------------------


} // end namespace subpavings

/*! Specialisation of std::swap for MappedSPnode type.*/
namespace std
{
	template <typename T>
	void swap(subpavings::MappedSPnode<T> & s1, 
			subpavings::MappedSPnode<T> & s2) // throw ()
	{
		s1.swapMSP(s2);
	}
}


#endif

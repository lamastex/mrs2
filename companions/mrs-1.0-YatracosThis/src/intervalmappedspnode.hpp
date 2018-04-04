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
\brief IntervalMappedSPnode declarations.
*/

#ifndef __INTERVALMAPPEDSP_HPP__
#define __INTERVALMAPPEDSP_HPP__


#include "mappedspnode.hpp"


#include <vector>
#include <list>

namespace subpavings {

    /*! \brief A derived class based on 
	 MappedSPnode < cxsc::interval >.

    The base class MappedSPnode is a node in the representation of a 
	mapped regular subpaving as a binary tree, where the type mapped 
	to the nodes is a real.  A node represents a box (interval vector).
    MappedSPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.  Each node will have
	a single interval value mapped to it.
	* 
	A IntervalMappedSPnode provides the functionality of the 
	MappedSPnode <interval> 
	and extended functionality appropriate for
	interval mapped nodes.
	
	@todo Could add the the arithmetic operations provided for intervals
	- see cxsc imath files.
	*/
	class IntervalMappedSPnode : 
						public MappedSPnode<cxsc::interval> {



    public:

		//typedefs
		typedef std::vector < IntervalMappedSPnode* >
			Ptrs;
		typedef Ptrs::iterator
			PtrsItr;
		typedef std::vector < const IntervalMappedSPnode* >
			ConstPtrs;
		typedef ConstPtrs::const_iterator
			ConstPtrsItr;
		typedef std::list <IntervalMappedSPnode* >
			ListPtrs;
		typedef ListPtrs::iterator
			ListPtrsItr;
		typedef std::list < const IntervalMappedSPnode* >
			ListConstPtrs;
		typedef ListConstPtrs::const_iterator
			ListConstPtrsItr;

		/*! \brief Inner interface for types measuring nodes. */
		class Measurer {
			public:
			virtual cxsc::real operator()(
				const IntervalMappedSPnode * const imspn) const = 0;
			virtual ~Measurer(){};
		};

		// ------------------------ public member functions ----------------

		virtual ~IntervalMappedSPnode();

		/*! \brief No-argument constructor. */
		IntervalMappedSPnode();


		/*! \brief initialised constructor.
		 
		Initialised with a box.
		Range is set to cxsc::interval(0.0,0.0). */
		explicit IntervalMappedSPnode(const ivector& v);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box. 
		Range is set to cxsc::interval(0.0,0.0).*/
		explicit IntervalMappedSPnode(const LabBox& lb);


		/*! \brief Initialised constructor.
		
		Initialised with a box and an interval value for the rangeCollection. */
		explicit IntervalMappedSPnode(const ivector& v, 
										const cxsc::interval& range);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box 
		and an interval value for the rangeCollection. */
		explicit IntervalMappedSPnode(const LabBox& lb,
										const cxsc::interval& range);

		/*! \brief Copy constructor.
		 * 
		Range is set to csxc::interval(0.0, 0.0).
		
		Copies from given SPnode downwards. */
		explicit IntervalMappedSPnode(const SPnode& other);

		/*! \brief Copy constructor.
		
		Copies from given %IntervalMappedSPnode downwards. */
		IntervalMappedSPnode(const IntervalMappedSPnode& other);

		/*! \brief Copy constructor.
		
		Copies from a given %MappedSPnode<cxsc::interval> 
		node downwards. */
		IntervalMappedSPnode(
					const MappedSPnode<cxsc::interval>& other);

		/*! \brief Copy assignment operator.
		
		Copies from given %IntervalMappedSPnode downwards. */
		IntervalMappedSPnode& operator=(IntervalMappedSPnode rhs);

		/*! \brief Copy assignment operator.
		
		Copies from a given %MappedSPnode<cxsc::interval> node downwards. */
		IntervalMappedSPnode& operator=(
							MappedSPnode<cxsc::interval> rhs);

		// parent and child accessors have to hide the base class implementation
		// this is not good but otherwise we get the base class return type
		// I've asked around and I can't find a way around it ...

		/*! \brief Accessor for the parent of a node.
		
		Returns a copy of the pointer to parent node.*/
		IntervalMappedSPnode* getParent() const;


		/*! \brief Accessor for the left child of a node.
		
		Returns a copy of the pointer to leftChild node. */
		IntervalMappedSPnode* getLeftChild() const;


		/*! \brief Accessor for the right child of a node.
		
		Returns a copy of the pointer to rightChild node. */
		IntervalMappedSPnode* getRightChild() const;

		/*! \brief Less-than operator.
		
		\return true iff the range of this is < range of \a rhs. */
		bool operator<(const IntervalMappedSPnode& rhs) const;

		/*! \brief Self-scalar addition operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to add to this.
		\post this has the same tree structure as before the operation 
		with a range equal to the range of this before 
		the operation plus cxsc::interval(\a val).  */
		IntervalMappedSPnode& operator+= (const real& val);
		
		/*! \brief Scalar addition operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to add to this to give the object returned.
		\return a mapped paving with the same tree structure as this 
		before the operation with range equal to the range of this before 
		the operation plus cxsc::interval(\a val).  */
		const IntervalMappedSPnode operator+ (const real& val) const;
		
		/*! \brief Self-scalar subraction operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to subtract from this.
		\post this has the same tree structure as before the operation 
		with a range equal to the range of this before 
		the operation minus cxsc::interval(\a val).  */
		IntervalMappedSPnode& operator-= (const real& val);
		
		/*! \brief Scalar subtraction operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to subtract from this to give the object returned.
		\return a mapped paving with the same tree structure as this 
		before the operation with range equal to the range of this before 
		the operation minus cxsc::interval(\a val).  */
		const IntervalMappedSPnode operator- (const real& val) const;
		
		/*! \brief Self-scalar multiplication operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to multiply this by.
		\post this has the same tree structure as before the operation 
		with a range equal to the range of this before 
		the operation multiplied by cxsc::interval(\a val).  */
		IntervalMappedSPnode& operator*= (const real& val);
		
		/*! \brief Scalar multiplication operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to multiply this by to give the object returned.
		\return a mapped paving with the same tree structure as this 
		before the operation with range equal to the range of this before 
		the operation multiplied by cxsc::interval(\a val).  */
		const IntervalMappedSPnode operator* (const real& val) const;

		/*! \brief Self-scalar division operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to divide this by.
		\post this has the same tree structure as before the operation 
		with a range equal to the range of this before 
		the operation divided by cxsc::interval(\a val).  
		\pre \a val != 0.0.*/
		IntervalMappedSPnode& operator/= (const real& val);
		
		/*! \brief Scalar division operator with real scalar.
		
		\note that because there are implicit conversions from 
		types integer and double to cxsc::real in the cxsc::library, 
		this operation also provides equivalent functionality for
		operand \a val of type integer and double. 
		
		\param val the value to divide this by to give the object returned.
		\return a mapped paving with the same tree structure as this 
		before the operation with range equal to the range of this before 
		the operation divided by cxsc::interval(\a val).  
		\pre \a val != 0.0.*/
		const IntervalMappedSPnode operator/ (const real& val) const;

		/*! Get a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.
		 * 
		\return NULL if no leaf node descendent of this contains
		\a pt, otherwise a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.*/
		const IntervalMappedSPnode* findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd  = ON_PARENT) const;

		/*! \brief  Add two sibling child nodes to this provided this is a leaf.
		
		Each new child gets half of the box associated with this, 
		splitting the box in half normal to dimension set by comp. */
		virtual void nodeExpand(int comp);

		/*! \brief  Add two sibling child nodes to this provided this is a leaf.
		
		Each new child gets half of the box associated with this, 
		splitting the box in half normal to the first longest dimension
		of the box associated with this. */
		virtual void nodeExpand();

		/*! \brief Propagate the interval hull of the children upwards.

		hullPropagates children then recalculates the range for this
		node as the interval hull of the ranges of the children.*/
		virtual void hullPropagation();
		
		/*! \brief Slice this.
		 * 
		* \param sliceDims is a vector of dimensions to slice on.
		* \param slicePts is a vector of points to slice on, assumed to 
		* correspond to the dimensions in \a sliceDims, ie the ith value
		* in \a sliceDims gives the dimension for the ith point in
		* \a slicePts.	*/
		virtual void slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts);


		/*! \brief  Get the diameter of interval range of this.*/
		cxsc::real getRangeDiameter() const;


		/*! \brief  Get the "area" of the interval range over the box
		of this.
		
		The "area" returned is calculated as the volume of the box
		represented by this multiplied by the diameter of the 
		interval range of this (which can also be seen as the volume
		of the box represented by this with the interval
		range of this added as an extra dimension).
		
		\return volume of box represented by this multiplied by
		diameter of interval range of this.*/
		virtual cxsc::real getIntervalAreaOnBox() const;
		
		/*! \brief  Get the total "area" of the interval ranges 
		over the boxes for the leaves of this.
		
		For each leaf descendent of this 
		the "area" returned is calculated as the volume of the box
		represented  multiplied by the diameter of the 
		interval range of the node (which can also be seen as the volume
		of the box represented with the interval
		range of the node added as an extra dimension).
		
		\return total over all leaf descendents of this of (volume
		of box represented multiplied by diameter of interval range).*/
		virtual cxsc::real getTotalLeafIntervalAreaOnBox() const;

		/*! \brief  Get the difference between the "area" of the 
		interval range over the box of this and the sum of the 
		"areas" of the interval ranges of the children over their
		boxes.
		 
		\return The volume of box represented by this multiplied by
		diameter of interval range of this, less sum of the same 
		for the children (returns interval area of this if no children.)*/
		cxsc::real getIntervalAreaDiffToChildren() const;

		/*! \brief Return a reference to a container of nodes.
 
		Contents of container are the leaves descended from this, 
		or this if this is a leaf, left to right order. */
		Ptrs& 
			getLeaves(Ptrs& leaves);

		/*! \brief Return a reference to a container of const nodes.
 
		Contents of container are the leaves descended from this, 
		or this if this is a leaf, left to right order. */
		ConstPtrs& 
			getConstLeaves(ConstPtrs& leaves) const;

		/*! \brief Return a reference to a container of nodes.
 
		Contents of container are the sub-leaves descended from this, 
		or this if this is a sub-leaf, left to right order. 

		A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
		Ptrs& 
			getSubLeaves(Ptrs& subleaves);

		/*! \brief Return a reference to a container of const nodes.
 
		Contents of container are the sub-leaves descended from this, 
		or this if this is a sub-leaf, left to right order. 

		A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
		ConstPtrs& 
			getConstSubLeaves(ConstPtrs& subleaves) const;


		void swapIMSP(IntervalMappedSPnode& spn); //throw()

		/*! \brief Print the details of a specific node in a subpaving.*/
		virtual std::ostream& nodePrint(std::ostream &os) const;

	private:
		
		/*! \brief Print the details of a single leaf node, using tab delimiters.*/
		virtual std::ostream& leafOutputTabs(std::ostream &os) const;

		
		/*! \brief A quick one-line summary of a node. 
		
		\param level Controls number of tabs indenting output.*/
		std::ostream& oneLineOutput(std::ostream& os,
											int level = 0) const;
		
			
	
}; // end IntervalMappedSPnode class



    // ----------------- non member tools functions ----------------------


} // end namespace subpavings

// Full specializations of the templates in std namespace can be added in std namespace.
namespace std
{
	template <>
	void swap(subpavings::IntervalMappedSPnode & s1, 
			subpavings::IntervalMappedSPnode & s2); // throw ()
	
	template <>
	void swap(cxsc::interval & i1, 
			cxsc::interval & i2); // throw ()

	
}

#endif

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
\brief BooleanValueMappedSPnode declarations.
*/

#ifndef __BOOLEANVALUEMAPPEDSP_HPP__
#define __BOOLEANVALUEMAPPEDSP_HPP__


#include "mappedspnode.hpp"
#include "booleanvalue.hpp"

#include <vector>
#include <list>

namespace subpavings {
	
	
    /*! \brief A derived class based on 
	 MappedSPnode < bool >.

    The base class MappedSPnode is a node in the representation of a 
	mapped regular subpaving as a binary tree, where the type mapped 
	to the nodes is a real.  A node represents a box (interval vector).
    MappedSPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.  Each node will have
	a single boolean value mapped to it.
	* 
	A BooleanValueMappedSPnode provides the functionality of the 
	MappedSPnode <bool> 
	and extended functionality appropriate for
	boolean mapped nodes.
	
	*/
	class BooleanValueMappedSPnode : 
						public MappedSPnode<BooleanMappedValue> {



    public:

		/*! @ name Definitions for collection types for
		 *  BooleanValueMappedSPnodes */
		typedef std::vector < BooleanValueMappedSPnode* >
			Ptrs;
		typedef Ptrs::iterator
			PtrsItr;
		typedef std::vector < const BooleanValueMappedSPnode* >
			ConstPtrs;
		typedef ConstPtrs::const_iterator
			ConstPtrsItr;
		typedef std::list <BooleanValueMappedSPnode* >
			ListPtrs;
		typedef ListPtrs::iterator
			ListPtrsItr;
		typedef std::list < const BooleanValueMappedSPnode* >
			ListConstPtrs;
		typedef ListConstPtrs::const_iterator
			ListConstPtrsItr;
			
		

		#if(0)
		/*! \brief Inner interface for types measuring nodes. */
		class Measurer {
			public:
			virtual cxsc::real operator()(
				const BooleanValueMappedSPnode * const imspn) const = 0;
			virtual ~Measurer(){};
		};
		#endif

		// ------------------------ public member functions ----------------

		/*! \brief Destructor. */
		virtual ~BooleanValueMappedSPnode();

		/*! \brief No-argument constructor. */
		BooleanValueMappedSPnode();


		/*! \brief initialised constructor.
		 
		Initialised with a box.
		Range is set to false. */
		explicit BooleanValueMappedSPnode(const ivector& v);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box. 
		Range is set to false.*/
		explicit BooleanValueMappedSPnode(const LabBox& lb);


		/*! \brief Initialised constructor.
		
		Initialised with a box and a boolean value for the range. */
		BooleanValueMappedSPnode(const ivector& v, 
										bool range);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box 
		and a boolean value for the range. */
		BooleanValueMappedSPnode(const LabBox& lb,
										bool range);
										
		/*! \brief Initialised constructor.
		
		Initialised with a box and a BooleanMappedValue for the range. */
		BooleanValueMappedSPnode(const ivector& v, 
							const BooleanMappedValue& range);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box 
		and a BooleanMappedValue for the range. */
		 BooleanValueMappedSPnode(const LabBox& lb,
							const BooleanMappedValue& range);

		/*! \brief Copy constructor.
		 * 
		Range is set to false.
		
		Copies from given SPnode downwards. */
		explicit BooleanValueMappedSPnode(const SPnode& other);

		/*! \brief Copy constructor.
		
		Copies from given %BooleanValueMappedSPnode downwards. */
		BooleanValueMappedSPnode(const BooleanValueMappedSPnode& other);

		/*! \brief Copy constructor.
		
		Copies from a given %MappedSPnode<BooleanMappedValue> 
		node downwards. */
		BooleanValueMappedSPnode(
					const MappedSPnode<BooleanMappedValue>& other);

		/*! \brief Copy assignment operator.
		
		Copies from given %BooleanValueMappedSPnode downwards. */
		BooleanValueMappedSPnode& operator=(BooleanValueMappedSPnode rhs);

		/*! \brief Copy assignment operator.
		
		Copies from a given %MappedSPnode<BooleanMappedValue> node downwards. */
		BooleanValueMappedSPnode& operator=(
							MappedSPnode<BooleanMappedValue> rhs);

		/*! \brief Replace the properties of this node and its descendents
		with the properties of another node and its descendents.
		
		Copies \a newNode node downwards into
		this, but keeps the relationship of this with 
		its parent, ie if this is a node in a tree, the part of the tree
		rooted at this is made identical to \a newNode but the relationship
		to the rest of the tree,through the link between this and its
		parent, is retained. */
		void replaceMe(BooleanValueMappedSPnode newNode);
		
		/*! \brief Replace the properties of this node and its descendents
		with the properties of another node and its descendents.
		
		Copies \a newNode node downwards into
		this, but keeps the relationship of this with 
		its parent, ie if this is a node in a tree, the part of the tree
		rooted at this is made identical to \a newNode but the relationship
		to the rest of the tree,through the link between this and its
		parent, is retained. */
		void replaceMe(MappedSPnode<BooleanMappedValue> newNode);

		// parent and child accessors have to hide the base class implementation
		// this is not good but otherwise we get the base class return type
		// I've asked around and I can't find a way around it ...

		/*! \brief Accessor for the parent of a node.
		
		Returns a copy of the pointer to parent node.*/
		BooleanValueMappedSPnode* getParent() const;


		/*! \brief Accessor for the left child of a node.
		
		Returns a copy of the pointer to leftChild node. */
		BooleanValueMappedSPnode* getLeftChild() const;


		/*! \brief Accessor for the right child of a node.
		
		Returns a copy of the pointer to rightChild node. */
		BooleanValueMappedSPnode* getRightChild() const;

		/*! \brief Get the value represented by this, as 
		 * a boolean.*/
		bool getRange() const;

		/*! Get a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.
		 * 
		\return NULL if no leaf node descendent of this contains
		\a pt, otherwise a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.*/
		const BooleanValueMappedSPnode* findContainingNode(
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

		#if(0)
		/*! \brief Slice this.
		 * 
		* \param sliceDims is a vector of dimensions to slice on.
		* \param slicePts is a vector of points to slice on, assumed to 
		* correspond to the dimensions in \a sliceDims, ie the ith value
		* in \a sliceDims gives the dimension for the ith point in
		* \a slicePts.	\param sliceFilename is the name of file to use to capture
		* the boxes of this that are used in the slice by outputting 
		* them to the file named \a sliceFilename.  Defaults to 
		* the empty string "".  If \a sliceFilename is the empty string
		* "" no boxes will be captured to file.*/
		virtual void slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts,
			const std::string& sliceFilename = "");

		#endif
		
		/*! \brief  Get the "area" of the true range over the box
		of this.
		
		The "area" returned is calculated as the volume of the box
		represented by this if the range is true, 0.0 otherwise.
		
		\return volume of box represented by this if range is true
		0.0 otherwise.*/
		virtual cxsc::real getTrueAreaOnBox() const;

		/*! \brief  Get the total "area" of the true ranges 
		over the boxes for the leaves of this.
		
		For each leaf descendent of this 
		the "area" returned is calculated as the volume of the box
		if the range it true, else 0.0.
		
		\return total over all leaf descendents of this of (volume
		of box represented for all true leaves).*/
		virtual cxsc::real getTotalLeafTrueAreaOnBox() const;

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


		void swapBMSP(BooleanValueMappedSPnode& spn); //throw()

		/*! \brief Print the details of a specific node in a subpaving.*/
		virtual std::ostream& nodePrint(std::ostream &os) const;
		
		/*! \brief Output all leaf nodes with true values mapped to them.*/
		virtual std::ostream& leavesOutputTabsTrue(std::ostream &os) const;
		
		/*! \brief Output all leaf nodes with true values mapped to them.
		 * 
		 * \param prec the precision for the output*/
		virtual std::ostream& leavesOutputTabsTrue(std::ostream &os, 
					int prec) const;

	private:
		#if(0)
			/*! \brief Give this a range equal to range of this AND range of other.

			Assumes that this and other are the same from leaves up to this and other,
			i.e. are roots of trees of the same shape up to the depth of this
			Note that other could be a taller trees than this: range collections
			are only combined up to the level of this. */
			void _addRanges(	const BooleanValueMappedSPnode * const other);

			/*! \brief Give this a range equal to range of this 
			 * XOR = symmetric set difference range of other.

			Assumes that this and other are the same from leaves up to this and other,
			i.e. are roots of trees of the same shape up to the depth of this
			Note that other could be a taller trees than this: range collections
			are only combined up to the level of this. */
			void _subtractRanges(	const BooleanValueMappedSPnode * const other);

			/*! \brief Give this a range equal to range of this 
			AND range of other.

			Assumes that this and other are the same from leaves up to this and other,
			i.e. are roots of trees of the same shape up to the depth of this
			Note that other could be a taller trees than this: range collections
			are only combined up to the level of this. */
			void _multRanges(	const BooleanValueMappedSPnode * const other);

			/*! \brief Give this a range equal to range of this 
			set difference range of other.

			Assumes that this and other are the same from leaves up to this and other,
			i.e. are roots of trees of the same shape up to the depth of this
			Note that other could be a taller trees than this: range collections
			are only combined up to the level of this. */
			void _divRanges(	const BooleanValueMappedSPnode * const other);



			/*! \brief Make a non-minimal union of subpavings using
			  union (OR) of ranges.
			*/
			void _addNonMinimalUnion(
						   const BooleanValueMappedSPnode& rhs);
			
			/*! \brief Make a non-minimal union of subpavings using
			  XOR = symmetric set difference of ranges.

			\pre this and *rhs are non-empty.*/
			void _subtractNonMinimalUnion(
						   const BooleanValueMappedSPnode& rhs);
			
			/*! \brief Make a non-minimal union of subpavings using
			intersection of ranges.

			\pre this and *rhs are non-empty.*/
			void _multiplyNonMinimalUnion(
						   const BooleanValueMappedSPnode& rhs);

			/*! \brief Make a non-minimal union of subpavings using
			set difference of ranges.

			\pre this and *rhs are non-empty.*/
			void _divideNonMinimalUnion(
						   const BooleanValueMappedSPnode& rhs);


		
			/*! \brief Union range collection of this.
    
			\param add the value to OR with the range of this.
			*/
			void _scalarAdd(bool add);
			
			/*! \brief XOR = symmetric set difference range collection of this.
			
			\param sub the value to XOR with range of this.
			*/
			void _scalarSubtract(bool sub);
			
			/*! \brief Intersect range of this.
			
			\param mult the value to AND with range of this.
			*/
			void _scalarMult(bool mult);
			

			/*! \brief Set difference range of this.
			
			\param div the value to set difference with the range of this.
			*/
			void _scalarDiv(bool div);
			
		#endif
		
		/*! \brief Print the details of a single leaf node, using tab delimiters.*/
		virtual std::ostream& leafOutputTabs(std::ostream &os) const;

		
		/*! \brief A quick one-line summary of a node. 
		
		\param level Controls number of tabs indenting output.*/
		std::ostream& oneLineOutput(std::ostream& os,
											int level = 0) const;
		
			
	
}; // end BooleanValueMappedSPnode class



    // ----------------- non member tools functions ----------------------


} // end namespace subpavings



// Full specializations of the templates in std namespace can be added in std namespace.
namespace std
{
	template <>
	void swap(subpavings::BooleanValueMappedSPnode & s1, 
			subpavings::BooleanValueMappedSPnode & s2); // throw ()
	
}
#endif

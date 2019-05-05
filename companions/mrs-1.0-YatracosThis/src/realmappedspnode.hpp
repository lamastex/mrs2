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
\brief RealMappedSPnode declarations.
*/

#ifndef __REALMAPPEDSP_HPP__
#define __REALMAPPEDSP_HPP__


#include "spsnode.hpp"
#include "mappedspnode.hpp"

#include <vector>
#include <list>

namespace subpavings {

    /*! \brief A derived class based on MappedSPnode < cxsc::real >.

    The base class MappedSPnode is a node in the representation of a 
	mapped regular subpaving as a binary tree, where the type mapped 
	to the nodes is a real.  A node represents a box (interval vector).
    MappedSPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.  Each node will have
	a single real value mapped to it.
	* 
	A RealMappedSPnode provides the functionality of the 
	MappedSPNode <real> and extended functionality appropriate for
	real mapped nodes, such as marginalisation.
	* 
	* Note that arithmetical operations with non-real operands can
	* take place provided that there is an implicit conversion from
	* the operand type to a cxsc::real, eg the operations this += 2 
	* or this + 2 will give the same results as this += cxsc::real(2.0)
	* or this + cxsc::real(2.0) respectively, because there is a 
	* implicit conversion from an integer to a cxsc::real in the cxsc 
	* library.  Similarly, operations with operands of type double 
	* will also succeed.
	* 
	@todo  Really it seems like we should not have both the RealMappedSPnode
	and the RealMappedSPnode classes.  Should look at the collator stuff
	and see if that can be the RealMappedSPnode and just replace all this 
	code by that and rename it.  That would mean that we lose the
	connection with the templatised MappedSPnode class which would be a 
	good thing really - the RealMappedSPnode is adding functionality
	to the extent that it is not appropriate to see it as a derived class
	from the mapped nodes any more. Combining collators and real mapped
	spnodes would deal with all the duplication of code that there is
	at present. */
	class RealMappedSPnode : public MappedSPnode<cxsc::real> {



    public:


		//typedefs
		typedef std::vector <RealMappedSPnode* >
			Ptrs;
		typedef Ptrs::iterator
			PtrsItr;
		typedef std::vector < const RealMappedSPnode* >
			ConstPtrs;
		typedef ConstPtrs::const_iterator
			ConstPtrsItr;
		typedef std::list <RealMappedSPnode* >
			ListPtrs;
		typedef ListPtrs::iterator
			ListPtrsItr;
		typedef std::list < const RealMappedSPnode* >
			ListConstPtrs;
		typedef ListConstPtrs::const_iterator
			ListConstPtrsItr;

		/*! \brief Inner interface for types measuring nodes. */
		class Measurer {
			public:
			virtual cxsc::real operator()(
				const RealMappedSPnode * const rmspn) const = 0;
			virtual ~Measurer(){};
		};
		
		
		// ------------------------ public member functions ----------------

		virtual ~RealMappedSPnode();

		/*! \brief No-argument constructor.*/
		RealMappedSPnode();


		/*! \brief initialised constructor.
		 
		Initialised with a box.
		Range is set to 0.0. */
		explicit RealMappedSPnode(const ivector& v);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box. 
		Range is set to 0.0.*/
		explicit RealMappedSPnode(const LabBox& lb);


		/*! \brief Initialised constructor.
		
		Initialised with a box and and a real value for the range.*/
		explicit RealMappedSPnode(
						const ivector& v, const cxsc::real& range);

		/*! \brief Initialised constructor.
		
		Initialised with a labeled box 
		and a real value for the rangeCollection. */
		explicit RealMappedSPnode(const LabBox& lb, 
											const cxsc::real& range);

		/*! \brief Constructor initialised with an SPnode.
		
		Range is set to csxc::real(0.0);
		
		\param spn An SPnode to copy.
		\pre \a spn must have a box.
		\post This has a subpaving identical to \a spn and the value
		mapped onto each node of this is 0.0.*/
		explicit RealMappedSPnode(const SPnode& spn);

		/*! Constructor initialised with an SPSnode.
		
		A NoBox_Error is thrown if \a other has no box. 
		
		\param spn An SPSnode to copy.
		\pre \a spn must have a box.
		\post This has a subpaving identical to \a spn and the value
		mapped onto each node of this is the value of count/nodeVolume
		from the equivalent node in \a spn.*/
		explicit RealMappedSPnode(const SPSnode& spn);


		/*! \brief Copy constructor.
		
		Copies from given %RealMappedSPnode downwards. */
		RealMappedSPnode(const RealMappedSPnode& other);

		/*! \brief Copy constructor.
		
		Copies from a given %MappedSPnode<cxsc::real> node downwards. */
		RealMappedSPnode(const MappedSPnode<cxsc::real>& other);

		/*! \brief Copy assignment operator.
		
		Copies from given %RealMappedSPnode downwards. */
		RealMappedSPnode& operator=(RealMappedSPnode rhs);

		/*! \brief Copy assignment operator.
		
		Copies from a given %MappedSPnode<cxsc::real> node downwards. */
		RealMappedSPnode& operator=(MappedSPnode<cxsc::real> rhs);

		// parent and child accessors have to hide the base class implementation
		// this is not good but otherwise we get the base class return type
		// I've asked around and I can't find a way around it ...

		/*! \brief Accessor for the parent of a node.
		
		Returns a copy of the pointer to parent node.*/
		RealMappedSPnode* getParent() const;


		/*! \brief Accessor for the left child of a node.
		
		Returns a copy of the pointer to leftChild node. */
		RealMappedSPnode* getLeftChild() const;


		/*! \brief Accessor for the right child of a node.
		
		Returns a copy of the pointer to rightChild node. */
		RealMappedSPnode* getRightChild() const;

		/*! \brief Less-than operator.
		
		\return true iff the range of this is < range of \a rhs. */
		bool operator<(const RealMappedSPnode& rhs) const;


		/*! Get a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.
		 * 
		\return NULL if no leaf node descendent of this contains
		\a pt, otherwise a pointer to the leaf node descendent of this 
		 whose box contains the point \a pt.*/
		const RealMappedSPnode* findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd  = ON_PARENT) const;

		/*! Get whether there is a negative range in the tree 
		rooted at this. The tree rooted at this includes this
		 node, subterminal nodes and leaf nodes).
		
		\return true if the tree rooted at this has a node with
		a negative range.*/
		bool hasNegativeRangeInTree() const;

		/*! Get whether there is an infinite range in the tree 
		rooted at this. The tree rooted at this includes this
		 node, subterminal nodes and leaf nodes).
		
		\return true if the tree rooted at this has a node with
		an infinite range (ie range cxsc::Infinity) .*/
		bool hasInfiniteRangeInTree() const;


		/*! \brief  Add two sibling child nodes to this provided this is a leaf.
		
		Each new child gets half of the box associated with this, 
		splitting the box in half normal to dimension set by comp. */
		virtual void nodeExpand(int comp);

		/*! \brief  Add two sibling child nodes to this provided this is a leaf.
		
		Each new child gets half of the box associated with this, 
		splitting the box in half normal to the first longest dimension
		of the box associated with this. */
		virtual void nodeExpand();

		/*! \brief  Get maximum value of the range for any of the leaves
		of the tree rooted at this.
		* 
		Returns the range of this if this is a leaf.
	
		\return maximum value of the range for any of the leaves of
		the tree rooted at this.	*/
		virtual cxsc::real getMaxRangeForLeavesInTree() const;
		
				
		/*! \brief Marginalise this.

		Marginalises to take out the unwanted dimensions and 
		adjust the range so that node vol x range 
		is the same as before marginalisation, and does
		this recursively for all children of this so that 
		after the marginalisation the sum over all leaves of (node vol x 
		range)
		is the same as before marginalisation.
		
		\note allowed dimensions start at 1, ie dimensions to
		marginalise on can include 1, 2, ... d 
		where d = dimensions of this.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NonRootNode_Error if this has
		a parent (ie. is not a root node).</li>
		<li>Throws a NoBox_Error if this is empty (has no box).</li>
		<li>Throws an std::invalid_argument if the required dimensions
		\a reqDim is empty or contains dimensions outside the 
		range of the dimensions of this.</li>
		</ul>
		
		\param reqDims is a vector of the dimensions to include in 
		marginal.
		\pre \a reqDims must be compatible with current dimensions
		and this must be a root node (have no parent) with a box.
		\post This will have a tree structure that is the 
		result of the marginalisation of the original subpaving
		to include only the required dimensions and will have
		a range adjusted to preserve the total range x volume
		given the change in
		the volume of the subpaving box that the node represents.*/
		virtual void marginalise(const std::vector<int>& reqDims);
		
		/*! \brief Make a marginalised version of this.

		Makes a new %RealMappedSPnode tree that s the marginal of this, 
		marginalising to take out the unwanted dimensions and 
		adjust the range so that node vol x range 
		is the same as before, and does
		this recursively for all children so that 
		after the marginalisation the sum over all leaves of (node vol x 
		range) of the new tree 
		is the same as for this.
		
		\note allowed dimensions start at 1, ie dimensions to
		marginalise on can include 1, 2, ... d 
		where d = dimensions of this.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if this is empty (has no box).</li>
		<li>Throws an std::invalid_argument if the required dimensions
		\a reqDim is empty or contains dimensions outside the 
		range of the dimensions of this.</li>
		</ul>
		
		\param reqDims is a vector of the dimensions to include in 
		marginal.
		\return  A RealMappedSPnode with a subpaving that is the 
		result of the marginalisation of the subpaving of this
		to include only the required dimensions and 
		a range adjusted to preserve the range x node volume
		given the change in
		the volume of the subpaving box that the node represents
		compared to that in this.
		\pre \a reqDims must be compatible with current dimensions
		and this must have a box.*/
		const RealMappedSPnode makeMarginalised(
							const std::vector<int>& reqDims) const;
		
		/*! \brief Normalise this.
		
		Normalising means changing the range of each node
		so that the total over all leaf nodes of the sum of
		(box volumes xabs(range)) is 1.0.   
	 
		Throws the following exceptions:
		<ul>
		<li>Throws a NonRootNode_Error if this has
		a parent (ie. is not a root node).</li>
		<li>Throws a NoBox_Error if this is empty (has no box).</li>
		<li>Throws an std::runtime_error if 
		this has no 'area', ie getTotalAbsLeafAreaRangeWithBox() <= 0.</li>
		<li>Throws an std::runtime_error if 
		this has no 'area', ie getTotalAbsLeafAreaRangeWithBox() 
		== Infinity.</li>
		</ul>
		
		\pre this must be a root node (have no parent) with a box and
		getTotalAbsLeafAreaRangeWithBox() != 0 and
		getTotalAbsLeafAreaRangeWithBox() != Infinity.
		\post this will have the same tree structure as before
		and the ranges on the nodes will be adjusted so that 
		getTotalAbsLeafAreaRangeWithBox() = 1.0.  */
		virtual void normalise();
		
		/*! \brief Make a normalised version of this.

		Makes a new %RealMappedSPnode tree that is the
		normalised version of this, 
		ie has the same subpaving but the total over all leaf nodes 
		of the sum of (box volumes x range) is 1.0.   
	
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if this is empty (has no box).</li>
		<li>Throws an std::runtime_error if 
		this has no 'area', ie getTotalAbsLeafAreaRangeWithBox() <= 0.</li>
		<li>Throws an std::runtime_error if 
		this has no 'area', ie getTotalAbsLeafAreaRangeWithBox() 
		== Infinity.</li>
		</ul>
		
		\return  A RealMappedSPnode with a subpaving that is the 
		the same as the subpaving of this this with ranges such that
		the total over all leaf nodes 
		of the sum of (box volumes x range) is 1.0. 
		\pre this must have a box and
		getTotalAbsLeafAreaRangeWithBox() != 0 and
		getTotalAbsLeafAreaRangeWithBox() != Infinity.*/
		const RealMappedSPnode makeNormalised() const;
				
		
		//20160904 - for KL computations
		std::pair<size_t, cxsc::real> getNonZeroBoxSummary() const;
		void smearRanges(cxsc::real zeroRange, cxsc::real ratioRange);

		
		/*! \brief Slice this.
		 * 
		* \param sliceDims is a vector of dimensions to slice on.
		* \param slicePts is a vector of points to slice on, assumed to 
		* correspond to the dimensions in \a sliceDims, ie the ith value
		* in \a sliceDims gives the dimension for the ith point in
		* \a slicePts.	*/
		void slice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts);
			
		/*! \brief Make a %RealMappedSPnode that is a slice of this.
		 * 
		* \param sliceDims is a vector of dimensions to slice on.
		* \param slicePts is a vector of points to slice on, assumed to 
		* correspond to the dimensions in \a sliceDims, ie the ith value
		* in \a sliceDims gives the dimension for the ith point in
		* \a slicePts.	
		\return A %RealMappedSPnode that is a slice of this
		on dimensions \a sliceDims at points \a slicePts.*/
		const RealMappedSPnode makeSlice(
			const std::vector < int >& sliceDims,
			const std::vector < cxsc::real >& slicePts) const;
		
		/*! Gets the L1 distance between this and another 
		%RealMappedSPnode.

        The L1 distance is defined as the sum of the absolute values
		of the differences in 'area' represented by the leaf nodes
		of this and the other paving.  The 'area' represented
		by a leaf node is the absolute value mapped to that leaf node 
		multiplied by
		the volume of the box associated with the node.  
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if either this or 
		by \a other have no box.</li>
		<li>Throws a IncompatibleDimensions_Error if the dimensions
		and sizes of the boxes of this and \a other are not the same. 
		</ul>
		 
		\note This will not attempt to adjust for any difference in
		total integral between this and \a other: the L1 distance 
		is simply taken as the difference between the 'areas' of
		the leaf boxes.  
		
		\note If this or \a other has leaf nodes with infinite ranges
		then the L1 distance between them will be infinite (cxsc::Infinity).  

        \param other the root of a 
		%RealMappedSPnode tree to calculate
		the L1 distance against.
        \pre Both this and \a other must have boxes and those
		boxes must be the same.
		\post this will be unchanged.       */
		cxsc::real getL1Distance(
					const RealMappedSPnode& other) const;

		/*! \brief Get a 'log likelihood' using 
		positive values from this and counts from \a spn.
		  
		Treats the rages mapped onto the nodes of this like 
		function values and uses the counts from \a spn together
		with these values to calculate a 'log-likelihood'.  
		* 
		If range < 0.0 (including -infinity) is encountered for any
		leaf where \a spn has points, the method will return 
		cxsc::SignalingNaN.
		Otherwise if range == 0.0 is encountered for any
		leaf where \a spn has points, the method will return
		-cxsc::Infinity.
		Otherwise if range == infinity is encountered for any
		leaf where \a spn has points, the method will return
		-cxsc::Infinity.
		 
		\param spn is the root of an SPSnode tree which 
		provides information about counts.
		\return the sum over the leaf nodes with positive values
		in the intersection of 
		the tree rooted at this and the tree rooted at \a spn of
		the product of the 
		log of the value on the leaf descendent of this and the
		count on the corresponding leaf descendent of \a spn. Returns
		cxsc::SignalingNaN if negative ranges are encountered, else
		returns -cxsc::Infinity if zero ranges are encountered, else
		returns cxsc::Infinity if any infinite ranges are encountered. 
		\pre this and \a spn are both non-empty and 
		both have the same box.*/
		cxsc::real getLogLikelihood(const SPSnode& spn) const;
				
		/*! \brief Get the the "area" of the range and the box of an element
		 *  of a Scheffe set based on a given box
		*/ 
		virtual cxsc::real
			getIntegralForScheffeElement(ivector& box, cxsc::real vol, bool split) const;

		
		/*! \brief  Get the "area" of the range and the box
		of this.
		
		The "area" returned is calculated as the volume of the box
		represented by this multiplied by the real range of this 
		(which can also be seen as the volume
		of the box represented by this with another dimension
		equal to an interval that has Inf 0 and Sup equal to the
		range of this).
		
		\return volume of box represented by this multiplied by
		real range of this.		
		\pre this must have a boxe.*/
		virtual cxsc::real getRealAreaRangeWithBox() const;
		
		/*! \brief  Get the "area" of the range and the box
		of this as a dotprecision type.
		
		The "area" returned should be the dotprecision product
		of the volume of the box represented by this and 
		the real range of this.  \b But note that if the range is
		cxsc::Infinity, the dotprecision result is the same
		as dotprecision(0.0).  This is due to the working of the cxsc
		library.  The results can therefore be highly misleading
		if the range is cxsc::Infinity 
		
		\return dotprecision product of volume of box represented by this
		and real range of this.
		\pre Range of this is  not cxsc::Infinity.
		\pre this must have a boxe.*/
		virtual cxsc::dotprecision getDotPrecisionAreaRangeWithBox() const;
		
		
		/*! \brief  Get the total "area" for the leaves of this
		of the real ranges and boxes.
		
		For each leaf descendent of this 
		the "area" returned is calculated as the volume of the box
		represented multiplied by the real
		range of the node (which can also be seen as the volume
		of the box represented with another dimension equal to an 
		interval that has Inf 0 and Sup equal to the range of the node).
		
		If there are an infinite ranges in the tree rooted at this, 
		the method will return cxsc::Infinity.
		
		\return total over all leaf descendents of this of (volume
		of box represented multiplied by real range) (returns 
		cxsc::Infinity if hasInfiniteRangeInTree()).
		\pre all leaf descendents of this must have boxes.*/
		virtual cxsc::real getTotalLeafAreaRangeWithBox() const;
		
		/*! \brief  Get the total "area" of the range and the box
		of this as a dotprecision type.
		
		The "area" returned should be the sum over all the 
		leaves of this of the dotprecision product
		of the volume of the box represented by the leaf node and 
		the real range of the leaf note.  \b But note that if the range
		of a leaf node is cxsc::Infinity, the dotprecision result for
		that leaf is the same
		as dotprecision(0.0).  This is due to the working of the cxsc
		library.  The results can therefore be highly misleading
		if any leaf node has range cxsc::Infinity. 
		
		\return dotprecision sum over leaf nodes of dotprecision 
		product of volume of box represented by the leaf and the 
		range of the leaf.
		\pre !hasInfiniteRangeInTree().		
		\pre all leaf descendents of this must have boxes.*/
		virtual cxsc::dotprecision 
				getTotalDotPrecisionLeafAreaRangeWithBox() const;
		
		/*! \brief  Get the total absolute "area" for the leaves of this
		of the real ranges and boxes.
		
		For each leaf descendent of this 
		the "area" returned is calculated as the volume of the box
		represented multiplied by the absolute value of the real
		range of the node (which can also be seen as the volume
		of the box represented with another dimension equal to an 
		interval that has Inf 0 and Sup equal to the range of the node).
		
		If there are an infinite ranges in the tree rooted at this, 
		the method will return cxsc::Infinity.
		
		\return total over all leaf descendents of this of abs(volume
		of box represented multiplied by real range) (returns 
		cxsc::Infinity if hasInfiniteRangeInTree()).
		\pre all leaf descendents of this must have boxes.*/
		virtual cxsc::real getTotalAbsLeafAreaRangeWithBox() const;
		
		/*! \brief  Get the total absolute "area" of the range and the box
		of this as a dotprecision type.
		
		The "area" returned should be the sum over all the 
		leaves of this of the absolute value of the dotprecision product
		of the volume of the box represented by the leaf node and 
		the real range of the leaf note.  \b But note that if the range
		of a leaf node is cxsc::Infinity, the dotprecision result for
		that leaf is the same
		as dotprecision(0.0).  This is due to the working of the cxsc
		library.  The results can therefore be highly misleading
		if any leaf node has range cxsc::Infinity. 
		
		\return dotprecision sum over leaf nodes of absolute value of 
		the dotprecision 
		product of volume of box represented by the leaf and the 
		range of the leaf.
		\pre !hasInfiniteRangeInTree().
		\pre all leaf descendents of this must have boxes.*/
		virtual cxsc::dotprecision 
				getTotalDotPrecisionAbsLeafAreaRangeWithBox() const;

		/*! \brief  Get the total over the leaves of the absolute
		value of the difference between the leaf "area" for this 
		and for \a rmsp.
		
		For each leaf descendent of the subpaving union of this and 
		\a rmsp 
		the "area" returned is calculated as the volume of the box
		represented multiplied by the absolute value of the difference
		between the range on that leaf for this and for \a rmsp real.
		
		This method should return the same as if the method 
		getTotalAbsLeafAreaRangeWithBox() was called on a node 
		which is the difference between this and \a rmsp 
		(ie *this - rmsp).
		
		If there are an infinite ranges in the tree rooted at this, 
		or \a rmps, the method will return cxsc::Infinity.
		
		\param rmsp the %RealMappedSPnode to calculate 
		the difference against.
		\return total over the leaves of the union of this and \a rmsp
		of the absolute
		value of the difference between the leaf "area" for this 
		and for \a rmsp (returns 
		cxsc::Infinity if hasInfiniteRangeInTree() or if 
		rmsp.hasInfiniteRangeInTree()).
		\pre this and \a rmsp have the same box.
		\pre all leaf descendents of this and \a rmsp must have boxes.*/
		virtual cxsc::real getTotalAbsDiffLeafAreaRangeWithBox(
						const RealMappedSPnode& rmsp) const;
		
		/*! \brief  Get the total over the leaves of the absolute
		value of the difference between the leaf "area" for this 
		and for \a rmsp as a dotprecision type.
		
		For each leaf descendent of the subpaving union of this and 
		\a rmsp 
		the "area" returned is calculated as the volume of the box
		represented multiplied by the absolute value of the difference
		between the range on that leaf for this and for \a rmsp real.
		\b But note that if the range
		of a leaf node of this or \a rmsp
		is cxsc::Infinity, the dotprecision result for
		that leaf is the same
		as dotprecision(0.0).  This is due to the working of the cxsc
		library.  The results can therefore be highly misleading
		if any leaf node in this or \a rmsp has range cxsc::Infinity. 
		
		This method should return the same as if the method 
		getTotalDotPrecisionAbsLeafAreaRangeWithBox() was called on 
		a node which is the difference between this and \a rmsp 
		(ie *this - rmsp).
		
		\param rmsp the %RealMappedSPnode to calculate 
		the difference against.
		\return dotprecision sum over the leaves of the union of this and \a rmsp
		of the absolute
		value of the difference between the leaf "area" for this 
		and for \a rmsp.
		\pre this and \a rmsp have the same box.
		\pre !hasInfiniteRangeInTree() and !rmsp.hasInfiniteRangeInTree().
		\pre all leaf descendents of this and \a rmsp must have boxes.*/
		virtual cxsc::dotprecision 
				getTotalDotPrecisionAbsDiffLeafAreaRangeWithBox(
						const RealMappedSPnode& rmsp) const;

		/*! \brief Return a reference to a container of nodes.
 
		Contents of container are the leaves descended from this, 
		or this if this is a leaf, left to right order. */
		Ptrs& getLeaves(Ptrs& leaves);

		/*! \brief Return a reference to a container of const nodes.
 
		Contents of container are the leaves descended from this, 
		or this if this is a leaf, left to right order. */
		ConstPtrs& getConstLeaves(ConstPtrs& leaves) const;

		/*! \brief Return a reference to a container of nodes.
 
		Contents of container are the sub-leaves descended from this, 
		or this if this is a sub-leaf, left to right order. 

		A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
		Ptrs& getSubLeaves(Ptrs& subleaves);

		/*! \brief Return a reference to a container of const nodes.
 
		Contents of container are the sub-leaves descended from this, 
		or this if this is a sub-leaf, left to right order. 

		A sub-leaf (aka "cherry") is a node with two leaf child nodes.*/
		ConstPtrs& getConstSubLeaves(ConstPtrs& subleaves) const;

		void swapRMSPSR(RealMappedSPnode& spn); //throw()

	private:
		
		//20160904
		void accumulateNonZeroBoxSummary(size_t& nNonZeroBoxes, 
								real& vNonZeroBoxVolumes) const;
		void _smearRanges(cxsc::real zeroRange, cxsc::real ratioRange);

		
		/*! \brief Non-public version of marginalisation to
		have only the required dimensions.*/
		void _start_marginalise(
				const std::vector<int>& reqDims);
		
		/*! \brief Non-public version of marginalisation to
		take out the unwanted dimensions.*/
		void _marginalise(const std::vector<int>& outDims);
		
		/*! \brief Non-public version of normalisation.*/
		void _normalise();
		
		/*! \brief Accumulate L1 distance between leaf
		 * descendents of this and leaf descendents of 
		 * node pointed to by \a other into \a disL1.*/
		cxsc::dotprecision& _getL1distance(
				cxsc::dotprecision& disL1,
				const RealMappedSPnode * const other) const;
		
		/*! \brief Accumulate L1 distance between this 
		 * and another identical node with value \a other_v
		 * mapped onto it, ie add nodevolume*(range - other_v)
		 * onto \a disL1.*/
		cxsc::dotprecision& nodeL1Distance(
								cxsc::dotprecision& disL1,
								cxsc::real other_v) const;
		
		/*! \brief Accumulate log likelihood using 
		 * positive values from this and counts from \a spn.
		 * 
		 * For all leaf nodes with positive values (ie > 0.0)
		 * in the intersection of the tree rooted
		 * at this and the tree rooted at \a spn, including this
		 * if this is a leaf, accumulates n_j * ln(v_j) into \a loglik
		 * where v_j is the value on the leaf node of this and 
		 * n_j is the count from the corresponding leaf node of
		 * \a spn. 
		 * 
		 * If v_j < 0.0 (including -infinity) is encountered for any
		 * leaf where \a spn has points, the isnan flag is set to 1.
		 * Otherwise if v_j == 0.0 is encountered for any
		 * leaf where \a spn has points, the isneginf flag is set to 1 
		 * and the log likelihood becomes -infinity.
		 * Otherwise if v_j == infinity is encountered for any
		 * leaf where \a spn has points, the isposinf flag is set to 1 
		 * and log likelihood becomes +infinity.*/
		cxsc::dotprecision& _getLogLikelihood(
								cxsc::dotprecision& loglik,
								int& isnan, int& isposinf, int& isneginf,
								const SPSnode * const spn) const;
														
		/*! \brief A quick one-line summary of a node. 
		
		\param level Controls number of tabs indenting output.*/
		std::ostream& oneLineOutput(std::ostream& os,
											int level = 0) const;
		
			
	
}; // end RealMappedSPnode class




    // ----------------- non member tools functions ----------------------

	/*! \brief Less-than operator using pointers.
		
	\return true iff *lhs < *rhs. */
	bool nodePtrCompare(const subpavings::RealMappedSPnode* lhs,
			const subpavings::RealMappedSPnode* rhs);


} // end namespace subpavings

// Full specializations of the templates in std namespace can be added in std namespace.
namespace std
{
	template <>
	void swap(subpavings::RealMappedSPnode & s1, 
			subpavings::RealMappedSPnode & s2); // throw ()
	
}

#endif

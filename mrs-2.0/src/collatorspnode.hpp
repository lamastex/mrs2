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
\brief CollatorSPnode declarations.
*/

#ifndef __NEW_COLLATORSP_HPP__
#define __NEW_COLLATORSP_HPP__

#include "rangecollection_hist.hpp"

#include "spsnode.hpp"

#include "spnode.hpp"

#include "SmallClasses.hpp"

#include <iostream>

namespace subpavings {
	
    /*! \brief A derived class based on SPnode.

    The base class SPnode is a node in the representation of a regular
    subpaving as a binary tree.  A node represents a box (interval vector).
    SPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.

	The CollatorSPnode class collates data from a collection of SPSnodes
    or objects from classes derived from SPSnode.

    An entire tree structure represents a non-minimal union over 
    the trees rooted at each node collated, ie
    the boxes of the leaves of the collated tree make a subpaving that is the
    union of the subpavings represented by each tree (tree rooted 
    at each node) collated.

    Each node has a container (the rangeCollection) that 
    represents the ordered collection of normalised 'heights' for
    the SPSnodes collated, where height is defined as
    (1/total of counters in entire tree) * node counter/node volume.
    The contents of this container are represented in the order in
    which they are added to the collation, the value for the last one added
    being at the end of the collection.
	
	In this documentation the SPSnodes collated by a CollatorSPnode 
	* are referred to as the 'elements' of the collator. The information
	* held in each node for each element gives the
	* 'height' of the equivalent node of that element.  The 'area' 
	* in a node represented for a single element in the collation is the 
	* element 'height' multipled by the volume of the box that the
	* node represents.  The area represented for a single element in 
	* the collation over the whole tree is the sum over all
	* the leaf nodes of the leaf node areas for that element.  If 
	* no modifying operations have been applied to the 
	* rangeCollection after the collation of an SPnode then the area 
	* of the element representing that SPnode in the collection will be
	* 1.0.  Operations which modify the rangeCollection so that the 
	* total area of an element is no longer 1.0 include 
	* reduceRangeCollection(), *=, and /=.  The total area represented
	* by a collation is the sum over all its elements of the areas 
	* represented by by each element.  
    */
    class CollatorSPnode : public SPnode {


		public:

			/*! \brief No-argument constructor. */
			CollatorSPnode();


			/*! \brief Initialised constructor with box.
		 
			Throws a MalconstructedBox_Error if \a v has no practical
			(0 or negative) dimensions.

			\param v interval vector for the box this represents.
			\pre v should be a proper ivector, with length >= 1.
			*/
			explicit CollatorSPnode(ivector& v);


			 /*! \brief Copy constructor.       */
			CollatorSPnode(const CollatorSPnode& other);

			/*! Constructor initialised with an SPSnode and 
			a normalising constant.
			
			The normalising constant is used to construct a 
			NewCollatorNode from an SPSnode.  The SPSnode is 
			represented in the collation by its normalised
			'height' 
			SPS node count /(normalising constant * SPS node vol).
			
			A NoBox_Error is thrown if \a spn has no box. 
			
			The given normaliser \a bigN defaults to 0.  If a value
			of \a bigN = 0 is passed to the method, it is replaced with
			the counter for the ultimate ancestor of this.
			
			If both the counter for the ultimate ancestor of this and
			\a bigN are 0, a std::logic_error exception is thrown.

			\param spn An SPSnode to collate into this.
			\param bigN the normalising constant, defaulting to 0.
			\pre \a spn must have a box.*/
			explicit CollatorSPnode(const SPSnode& spn,
									size_t bigN = 0);

			/*! Destructor.*/
			virtual ~CollatorSPnode();

			/*! Copy assignment operator. */
			CollatorSPnode& operator=(CollatorSPnode rhs);


			// parent and child accessors have to hide the base class implementation
			// this is not good but otherwise we get the base class return type
			// I've asked around and I can't find a way around it ...
			
			/*
			 * These accessor methods shadow equivalent methods in the base
			class.  Thus the method used is determined  at compile time,
			not run time as would be the case if virtual methods were used.
			Because the pointers to parents and children are part of the
			base class definition, the methods have to cast the base class
			form to the derived class form in order for the pointer
			returned to be able to be used with derived class members.*/

			// Accessor for the parent of a node.
			//Returns a copy of the pointer to parent node.
			/*! @name Accessors for links between the nodes.
			Note that pointers for parent, leftChild, and rightChild are
			not reference counted so there could potentially be problems
			with the use of returned pointers (for instance, being used to
			delete nodes). These pointers might be better implemented with
			boost::shared_ptr .
			*/

			// Hide the base class versions of these methods.
			//@{
			/*! \brief Accessor for the parent of a node.

			Returns a pointer to parent node.    */
			CollatorSPnode* getParent() const;


			/*! \brief Accessor for the left child of a node.

			Returns a pointer to leftChild node.        */
			CollatorSPnode* getLeftChild() const;


			/*! \brief Accessor for the right child of a node.

			Returns a pointer to rightChild node.        */
			CollatorSPnode* getRightChild() const;

			/*! \brief Get the total of the values held for all 
			the nodes collated by this.*/
			virtual cxsc::real getTotalRangeCollection() const;
			
			/*! \brief Get whether this node has nothing collated.
			
			\return true if this has nothing collated, false otherwise.*/
			virtual bool isEmptyRangeCollection() const;

			/*! \brief Get the number of subpavings collated into this.*/
			virtual std::size_t getSizeRangeCollection() const;
			
			/*! \brief Get the total 'area' represented by this collation.
			
			The total area represented by the collation is defined as the
			the sum of (sum of [absolute] height values in collation) 
			x vols for leaf nodes for which this is ancestor.
			
			\return The total 'area' represented by this.*/
			virtual cxsc::real getTotalAbsValueTimesVol() const;
			
			//new
			/*! \brief Fill in a container of height values for this
			and its descendants.
			* 
			Fills in order this, left child and descendents, right
			child and descendants.
			
			\return The container, filled in for this and its
			descendants.*/
			std::vector < RealVec >& getAllRangeCollections(
					std::vector < RealVec >& container) const;
					
			//new AHABC
			// get a container of Fhats (normalised heights, left to right)
			virtual RealVec& getLeafNodeFHatValues(RealVec& vals) const;

			/*! @name Get a collection of the L1 distance values
			for each element of the collation against the average
			element in the collation.
			 
			The L1 distance for a element against the average is 
			defined over all the leaves 
			of the entire collation.  The L1 distance 
			for an element in this is the sum,
			over these leaves, of the absolute differences 
			between the 'height' value for the element for 
			the leaf and the 'height' value for the average for the
			leaf, multiplied by the volume
			of the leaf.
			
			Throws a std::runtime_error if this has nothing collated, ie
			if isEmptyRangeCollection() == true
			
			\param container a reference to a container to use to 
			store the L1 distance values.  Any contents of the 
			given container will be discarded before new values
			are added.  
			\return An ordered collection of the L1 distance values
			for each element of the collation against the average
			element in the collation, in the same order as the 
			elements are added to this collation .
			\pre this should have at least one element.*/
			//@{
			virtual RealVec getL1DistancesToAverage() const;

			virtual RealVec& getL1DistancesToAverage(RealVec& container) const;
			//@}
			
			// this is not in general symmetric:  it is symmetric in the
			// special case that this and other both have just a single value
			// in their range collections, otherwise other is averaged
			// before we find the distances
			/*! @name Get a collection of the L1 distance values
			for each element of the collation against the average
			element over another collation.
			 
			The L1 distance for a element of this against the average 
			for another collation is 
			defined over all the leaves of a non-minimal union
			of the entire collations of this and the other.
			The L1 distance for an element in this is the sum,
			over these leaves, of the absolute differences 
			between the 'height' value for the element for
			the leaf and the 'height' value for the average of 
			the other for the leaf, multiplied by the volume
			of the leaf.
			* 
			Throws the following exceptions:
			<ul>
			<li>Throws a NullSubpavingPointer_Error if \a other is NULL.</li>
			<li>Throws a NoBox_Error if either this or the subpaving
			pointed to by \a other are empty.</li>
			<li>Throws an IncompatibleDimensions_Error if the boxes of 
			this and the subpaving pointed to by \a other do not have
			the same dimensions and size.</li>\
			<li>Throws a std::runtime_error if the collation of
			 the subpaving pointed to by \a other  has nothing
			 collated.</li>
			</ul>
			
			\param other a pointer to another collator node:  L1 distances
			are calculated for all the elements of this
			against the average of the colallation pointed to by \a other. 
			\param container a reference to container to use to 
			store the L1 distance values.  Any contents of the 
			given container will be discarded before new values
			are added.  
			\return An ordered collection of the L1 distance values
			for each element of the collation against the average
			element in \a other, in the same order as the 
			elements are added to this collation .
			\pre \a other should not be a NULL pointer, both this and 
			the collated pointed to by \a other should have boxes 
			and those boxes should be of equal sizes and dimensions, and
			the subpaving pointed to by \a other should have at least one
			element.*/
			//@{
			virtual RealVec getL1DistancesToAverage(
							const CollatorSPnode * const other) const;
							
			virtual RealVec& getL1DistancesToAverage(RealVec& container, 
							const CollatorSPnode * const other) const;
			//@}
			
			/*! @name Get a collection of the L1 distance values
			for each element of the collation against a statistical
			subpaving.
			 
			The L1 distance for a element of this against a
			statistical subpaving is 
			defined over all the leaves of a non-minimal union
			of this collation and the given statistical
			subpaving.  The L1 distance for an element in this is the sum,
			over these leaves, of the absolute differences 
			between the 'height' value for that element for
			the leaf and the 'height' value of the statistical subpaving
			for the leaf (ie counter/volume normalised by total count 
			in the whole paving), multiplied by the volume
			of the leaf.
			
			Throws a NullSubpavingPointer_Error if \a adh is NULL.
			
			\param adh A pointer to a statistical subpaving to 
			calculate distances against.
			\param container a reference to container to use to 
			store the L1 distance values.  Any contents of the 
			given container will be discarded before new values
			are added.  
			\return An ordered collection of the L1 distance values
			for each element of the collation against \a adh,
			in the same order as the 
			elements are added to this collation.
			\pre \a spn should be a non-NULL pointer.*/
			//@{
			virtual RealVec getL1Distances(
							const SPSnode * const spn) const;
							
			virtual RealVec& getL1Distances(RealVec& container, 
							const SPSnode * const spn) const;
			//@}

			/*! \brief Get a container of all descendent leaf nodes.

			Will contain just this if this is a leaf.

			\param leaves a reference to a container of node pointers to fill in.
			\return a reference to the container \a leaves 
			filled with pointers to leaf nodes.   */
			std::vector< const CollatorSPnode * >& getConstLeaves(
				std::vector< const CollatorSPnode * >& leaves) const;
			
			/*! \brief Output for for <b>all leaves</b> of a binary tree

			Output intended for a txt file, in numeric form only.

			Outputs the average over the summary for each leaf.
			\param os is the stream to send to.
			\param prec is the precision used for printing, defaulting to 5.
			*/
			virtual std::ostream& leavesAverageOutputTabs(std::ostream &os,
							int prec = 5) const;

			/*! Print the details of a specific node in a subpaving
			\param os is the stream to send to.*/
			virtual std::ostream& nodePrint(std::ostream &os) const;

			/*!@name Output for all the boxes in this
			
			But, unlike nodesAllOutput, does not do tab indentation
			
			\param os is the stream to send to.
			\param prec is the precision used for printing.	*/
			//@{
			virtual std::ostream& nodesAllPrint(std::ostream &os) const;
			virtual std::ostream& nodesAllPrint(std::ostream &os,
							int prec) const;
			//@}
			
			/*! @name Output the range collection for this
			\param os is the stream to send to.
			\param prec is the precision used for printing.*/
			//@{
			virtual std::ostream& outputNodeRangeCollection(
												std::ostream &os) const;
			virtual std::ostream& outputNodeRangeCollection(
								std::ostream &os, int prec) const;
			//@}
			
			/*! \brief Allocate collections of values to become 
			rangeCollections to this and children.
			
			Allocation order is this , left child and its decsendents, 
			then right child and its descendents.
			
			Used for testing.
			
			Throws a std::invalid_argument exception if the number of elements
			in \a rangesToAllocate is less than the number of elements collated
			in this. 
			
			\param rangesToAllocate is a collection of ranges to allocate.
			\param index is the index of the range in 
			\a rangesToAllocate to this.
			\return the index of the next range in the collection to
			be allocated.
			\pre \a rangesToAllocate should have the same number of ranges
			in as there are elements in this collation.*/
			std::size_t allocateRanges(
						const std::vector< RealVec >& rangesToAllocate,
						std::size_t index = 0);
			
			/*! \brief Expand a leaf node.

			Expand a leaf node to have two children and pass
			rangeCollection down to children.

			Throws a NoBox_Error if this is empty.

			\param comp is the dimension on which to to bisect theBox.
			\pre this should be non-empty.*/
			virtual void nodeExpand(int comp);


			/*! \brief Expand a leaf node.

			Expand a leaf node to have two children and pass
			rangeCollection down to children.

			Finds the splitting dimension.
			*/
			virtual void nodeExpand();

			/*! \brief Find the leaf node which would contain a data point.

			\note for efficiency, there is no check that the data point 
			given has the same dimension as this;  the results of trying this
			operation with a point with incompatible dimensions are undefined. 

			Throws a NoBox_Error if this is empty.

			\param pt a data point.
			\param childInd an indicator for whether the current node is a
			treated as a left or right child or a root.  
			Defaults to ON_PARENT.
			\return a pointer to the node which would contain the point,
			or NULL if no node contains the point.
			\pre This should have a box to check for containment against.*/
			const CollatorSPnode* findContainingNode(
								const cxsc::rvector& pt,
								OPERATIONS_ON childInd = ON_PARENT) const;
  

			// nodeReabsorbChildren() can use the base class implementation
			// (the rangeCollection for this will be correct so just delete the children)


			/*! \brief Add the contents of another collator to this.

			Tree structures are combined using a non-minimal union
			and the rangeCollection for the subpaving pointed to by \a 
			rhs is appended to the end 
			of this's rangeCollection.  The effect is as though all 
			the statistical subpavings collated by the subpaving pointed
			to by \a rhs had been collated into this. 

			If \a rhs points to NULL, or the collation pointed to by
			\a rhs has no box, this is unchanged by the operation. 
			
			Throws an IncompatibleDimensions_Error if both this and 
			the subpaving pointed to by \a rhs have boxes and the 
			dimensions or sizes of those boxes are not the same.

			\param rhs pointer to the root of the collator to be added
			to this.
			\pre If both this and the subpaving pointed 
			to by \a rhs have boxes, those boxes should have the same 
			dimensions and sizes.  
			\post The tree structure (subpaving structure) for this becomes 
			the result of the non-minimal union of the tree structure
			(subpaving) for this and the tree structure (subpaving) 
			for the subpaving pointed to by \a rhs.  The rangeCollection 
			for the subpaving pointed to by \a rhs is 
			appended to the end of the original rangeCollection
			for this.  */
			void addPaving(const CollatorSPnode * const rhs);

			/*! \brief Collate a statistical subpaving into this.

			Tree structures are combined using a non-minimal union
			and information giving a  'height' value for the equivalent
			node in the subpaving pointed to by \a rhs is appended to 
			the end of this's rangeCollection.

			If \a rhs points to NULL or has no box, 
			this is unchanged by the operation. 
			
			Throws an IncompatibleDimensions_Error if both this and 
			the subpaving pointed to by \a rhs have boxes and the 
			dimensions or sizes of those
			boxes are not the same.

			\param rhs pointer to the root of the statistical
			subpaving to be added to this.
			\pre If both this and the subpaving pointed 
			to by \a rhs have boxes, those boxes should have the same 
			dimensions and sizes.  
			\post The tree structure (subpaving structure) for this becomes 
			the result of the non-minimal union of the tree structure
			(subpaving) for this and the tree structure (subpaving) 
			for the subpaving pointed to by \a rhs.  Information 
			giving 'height'  for the subpaving pointed to by \a rhs is 
			appended to the end of the original rangeCollection
			for this.  */
			void addPaving(const SPSnode * const rhs);

			/*! Reduce the rangeCollection of this to represent the sum of the
			heights of the statistical subpavings collated.
			 
			\post If the before the operation the rangeCollection had
			> 1 values, then afterwars the rangeCollection will contain
			one value representing the sum of the heights of the 
			statistical subpavings collated.  If before the operation 
			the rangeCollection had <= 1 elements, then is unchanged
			by the operation.*/
			virtual void reduceCollation();
			
			/*! \brief Add the contents of another collator to this.

			Tree structures are combined using a non-minimal union
			and the rangeCollection for \a add is appended to the end 
			of this's rangeCollection.  The effect is as though all 
			the statistical subpavings collated by \a add had been 
			collated into this. 

			Throws an IncompatibleDimensions_Error if both this and 
			\a add have boxes and the dimensions or sizes of those
			boxes are not the same.

			\param add the root of the collator to be added
			to this.
			\pre If both this and \a rhs have boxes, 
			those boxes should have the same 
			dimensions and sizes.  
			\post The tree structure (subpaving structure) for this becomes 
			the result of the non-minimal union of the tree structure
			(subpaving) for this and the tree structure (subpaving) 
			for \a add.  The rangeCollection for \a add is 
			appended to the end of the original rangeCollection
			for this.  */
			CollatorSPnode& operator+= (const CollatorSPnode& add);
			
			/*! \brief Add the contents of two collators together.

			The collator returned has a tree structures (subpaving) 
			that is teh non-minimal union of the trees of this and \a add.
			The collator returned has a rangeCollection 
			that is the result of appending the rangeCollection for 
			\a add to the end of this's rangeCollection.

			Throws an IncompatibleDimensions_Error if both this and 
			\a add have boxes and the dimensions or sizes of those
			boxes are not the same.

			\param add the root of the collator to be added
			to this.
			\return a collator with a tree structure (subpaving 
			structure) that is the result of the non-minimal union of
			the tree structure (subpaving) for this and the tree 
			structure (subpaving) for \a add.  The rangeCollection for
			the collator returned is the result of appending the 
			rangeCollection of \a add to the end of the rangeCollection
			for this.  
			\pre If both this and \a rhs have boxes, 
			those boxes should have the same 
			dimensions and sizes.  */
			const CollatorSPnode operator+ (const CollatorSPnode& add) const;

			/*! \brief Add a representation of a statistical
			subpaving into this collation.

			Tree structures are combined using a non-minimal union
			and information about the 'height' of the 
			equivalent node to this in \a add is appended to the end 
			of this's rangeCollection.  

			Throws an IncompatibleDimensions_Error if both this and 
			\a add have boxes and the dimensions or sizes of those
			boxes are not the same.

			\param add the root of the statistical subpaving to be added
			to this.
			
			\pre If both this and \a rhs have boxes, 
			those boxes should have the same 
			dimensions and sizes.  
			\post The tree structure (subpaving structure) for this becomes 
			the result of the non-minimal union of the tree structure
			(subpaving) for this and the tree structure (subpaving) 
			for \a add.  Information about the 'height' of the 
			equivalent node to this in \a add is appended to the end 
			of this's rangeCollection.  */
			CollatorSPnode& operator+= (const SPSnode& add);
			
			/*! \brief Add the contents of a collator and a 
			statistical subpaving together.

			The tree structures of the result is the non-minimal union
			of the tree structures of this and \a add.  The 
			rangeCollection of the result contains information about 
			the 'height' of the 
			equivalent node to this in \a add appended to the end 
			of this's rangeCollection.  

			Throws an IncompatibleDimensions_Error if both this and 
			\a add have boxes and the dimensions or sizes of those
			boxes are not the same.

			\param add the root of the statistical subpaving to be added
			to this.
			\return the result of adding a representation of the
			statistical subpaving \a add to this collation.
			\pre If both this and \a rhs have boxes, 
			those boxes should have the same 
			dimensions and sizes.  */
			const CollatorSPnode operator+ (const SPSnode& add) const;


			// subtraction makes no sense for collators - we cannot do the opposite of concatenate
			
			// element by element multiplication or division by a scalar
			
			/*! \brief Multiply this collation by a scalar.

			Each element in the rangeCollection for this is
			changed to represent a 'height' \a mult x the 'height' it 
			represented before the operation.  

			\param mult the scalar multiplier to apply to this.
			
			\post The tree structure (subpaving structure) for this 
			is unchanged. Each element in the rangeCollection for this is
			changed to represent a 'height' \a mult x the 'height' it 
			represented before the operation.  */
			CollatorSPnode& operator*= (cxsc::real& mult);

			/*! \brief Get the results of multiplying this collation
			by a scalar.

			Each element in the rangeCollection for of the collation
			returned is the result of changing the rangeCollection of
			this so that each element represents a 'height' 
			\a mult x the 'height' it represented before the operation.  

			\param mult the scalar multiplier to apply.
			
			\return A collation which has the same tree structure 
			(subpaving structure) as this. Each element in the 
			rangeCollection of the collation
			returned is the result of changing the rangeCollection of
			this so that each element represents a 'height' 
			\a mult x the 'height' it represented before the operation. */
			const CollatorSPnode operator* (cxsc::real& mult) const;

			/*! \brief Divide this collation by a scalar.

			Each element in the rangeCollection for this is
			changed to represent the 'height' it 
			represented before the operation divided by \a div. 
			* 
			Throws a std::invalid_argument exception if \a div == 0.0.

			\param div the scalar divisor to apply to this.
			\pre \a div != 0.0.
			\post The tree structure (subpaving structure) for this 
			is unchanged. Each element in the rangeCollection for this is
			changed to represent a 'height' = the 'height' it 
			represented before the operation divided by \a div.  */
			CollatorSPnode& operator/= (cxsc::real& div);

			/*! \brief Get the results of dividing this collation
			by a scalar.

			Each element in the rangeCollection for of the collation
			returned is the result of changing the rangeCollection of
			this so that each element represents a 'height' 
			= the 'height' it represented before the operation divided
			by \a div.  

			Throws a std::invalid_argument exception if \a div == 0.0.

			\param div the scalar divisor to apply.
			\return A collation which has the same tree structure 
			(subpaving structure) as this. Each element in the 
			rangeCollection of the collation
			returned is the result of changing the rangeCollection of
			this so that each element represents a 'height' = the 'height'
			it represented before the operation divided by \a div.  
			\pre \a div != 0.0.*/
			const CollatorSPnode operator/ (cxsc::real& div) const;
			
			/*! \brief Marginalised this.

			Marginalises to take out the unwanted dimensions and 
			adjust the rangeCollection so that for each element in the
			collation, sum of (node vol x 'height' collated) 
			is the same as before marginalisation, and hence that the 
			overall sum over all leaves of (node vol x 
			accumulated 'heights' collated)
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
			a rangeCollection with the same number of elements
			as before marginalisation but where each 'height' is
			adjusted to preserve the area given the change in
			the volume of the subpaving box that the node represents.*/
			virtual void marginalise(const std::vector<int>& reqDims);
			
			/*! \brief Make a marginalised version of this.

			Makes a new collator that this the marginal of this, 
			marginalising to take out the unwanted dimensions and 
			adjust the rangeCollection for this so that for each element
			in the collation, sum of (node vol x 'height' collated) 
			is the same as for this, and hence that the 
			overall sum over all leaves of (node vol x 
			accumulated 'heights' collated)
			is the same as this.
			
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
			\return  A collator with a subpaving that is the 
			result of the marginalisation of the subpaving of this
			to include only the required dimensions and 
			a rangeCollection with the same number of elements
			as this but where each 'height' is
			adjusted to preserve the area given the change in
			the volume of the subpaving box that the node represents
			compared to that in this.
			\pre \a reqDims must be compatible with current dimensions
			and this must have a box.*/
			const CollatorSPnode makeMarginalised(
								const std::vector<int>& reqDims) const;
			
			
			/*! \brief Normalise this.
			
			Normalising means changing the rangeCollection to represent
			the (normalised) 'height' of a single statistical subpaving.
			ie. the rangeCollection for a node becomes a single value 
			such that the sum over all leaf nodes of the leaf node areas
			(node volume x 'height') is 1.0.   
		 
			\note normalise() will give the same result as average 
			if applied to a node where each element has total
			area = 1 (i.e the total area represented by all elements in 
			the collation is equal to the number of elements in the
			collation), but not if those elements do not each have area
			 1.0

			Throws the following exceptions:
			<ul>
			<li>Throws a NonRootNode_Error if this has
			a parent (ie. is not a root node).</li>
			<li>Throws a NoBox_Error if this is empty (has no box).</li>
			<li>Throws an UnfulfillableRequest_Error if this has nothing
			collated.</li>
			<li>Throws an std::logic_error if 
			this has no 'area', ie getTotalAbsValueTimesVol() == 0.</li>
			</ul>
			
			\pre this must be a root node (have no parent) with a box and
			at least one element collated and 
			getTotalAbsValueTimesVol() != 0
			\post this will have the same tree structure as before
			and the collation will include one element and the area 
			represented by that element will be 1.0.  */
			virtual void normalise();
			
			/*! \brief Make the normalised version of this.
			
			Normalising means changing the rangeCollection to represent
			the (normalised) 'height' of a single statistical subpaving.
			ie. the rangeCollection for a node becomes a single value 
			such that the sum over all leaf nodes of the leaf node areas
			(node volume x 'height') is 1.0.   
		 
			\note makeNormalised() will give the same result as 
			makeAverage() if applied to a node where each element has total
			area = 1 (i.e the total area represented by all elements in 
			the collation is equal to the number of elements in the
			collation), but not if those elements do not each have area
			 1.0
			
			Throws the following exceptions:
			<ul>
			<li>Throws a NoBox_Error if this is empty (has no box).</li>
			<li>Throws an UnfulfillableRequest_Error if this has nothing
			collated.</li>
			<li>Throws an std::logic_error if 
			this has no 'area', ie getTotalAbsValueTimesVol() == 0.</li>
			</ul>
			
			\return A collator with the same tree structure (subpaving)
			as this and one element in its collation such that the area 
			represented by that element is 1.0.
			\pre this must have a box and
			at least one element collated and 
			getTotalAbsValueTimesVol() != 0.*/
			const CollatorSPnode makeNormalised() const;
			
			/*! \brief Average this.
			
			Averaging means changing the rangeCollection to represent
			the average 'height' over the elements previously in
			the collation.
			ie. the rangeCollection for a node becomes a single value 
			representing a height that is the average height over
			all the elements previously in the collation.
		 
			\note average() will give the same result as normalise() 
			if applied to a node where each element has total
			area = 1 (i.e the total area represented by all elements in 
			the collation is equal to the number of elements in the
			collation), but not if those elements do not each have area
			 1.0.
			 
			Throws the following exceptions:
			<ul>
			<li>Throws a NonRootNode_Error if this has
			a parent (ie. is not a root node).</li>
			<li>Throws an UnfulfillableRequest_Error if this has nothing
			collated.</li>
			</ul>
			
			\pre this is a root node with at least one element collated.
			\post this will have the same tree structure as before
			and the collation will include one element representing
			the average of the heights of the elements previously
			in the collation.  */
			virtual void average();

			/*! \brief Make the average of this.
			
			Averaging means making a rangeCollection to represent
			the average 'height' over the elements averaged.
			ie. the rangeCollection for a node is a single value 
			representing a height that is the average height over
			all the elements in the collation to be averaged.
		 
			Throws an UnfulfillableRequest_Error if this has nothing
			collated.
			
			\return a collator with the same tree structure as this
			and a rangeCollection with one element representing
			the average of the heights of the elements in the 
			equivalent node of this.  
			\pre this has at least one element collated.*/
			const CollatorSPnode makeAverage() const;

			/*! \brief Swap this and another node.
		
			Swaps all the data members of this with the other node. 
					
			\param spn a reference to the node to swap with
			\post this is identical,in terms of its data members, 
			to spn before the swap, and spn is
			identical to this after the swap.*/
			void swapCollator(CollatorSPnode& spn); //throw()
			
			/*! \brief Get a string summary of this node's properties.
		
			Just this node, not this and its descendents.
			
			\return the string summary.		*/
			std::string nodeStringSummary() const;

		protected:
		
			// -------------------------- protected member functions -------------

			/*! \brief Constructor initialised with a box
			and a collection for the rangeCollection. */
			explicit CollatorSPnode(
								ivector& v, 
								const RangeCollectionHist& rangeColl);


			/*! Print the details of a single leaf node, using tab delimiters
			\param os is the stream to send to.		*/
			virtual std::ostream& leafOutputTabs(std::ostream &os) const;

			 /*! \brief Output for a node in a binary tree, tab-delimited

			Output intended for a txt file, in numeric form only.
			Outputs the average over the summary.

			Replaces the format that that the cxsc::<< operator produces
			for interval vectors.   The format used here 
			produces alpha-numeric tab-delimited data.  The format
			for an n-dimensional interval vector is

			nodeName [tab]  volume [tabl] average summary [tab]
			Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
			[tab] Inf(ivector[n]) [tab] Sup(ivector[n]
			
			\param os is the stream to send to
			\param prec is the precision used for printing, defaulting to 5
			*/
			virtual std::ostream& leafAverageOutputTabs(std::ostream &os,
											int prec = 5) const;

			/*! \brief A quick one-line summary of a collator node. */
			virtual std::ostream& oneLineOutput(std::ostream& os,
													int level = 0) const;

			/*! \brief Get a copy of the rangeCollection. */
			RangeCollectionHist getRangeCollection() const;

			/*! \brief Get the total area of the elements collated in this. */
			virtual cxsc::real getNodeAbsValueTimesVol() const;

			/*! \brief Check that there are no negative values
			in the 'heights' represented in this collation.  Return true
			if no negative heights found. */
			virtual bool checkNoNegativeTotalValues(
										bool checkSoFar = true) const;
			
			/*! \brief Get the average of the heights of the elments
			in this. */
			virtual cxsc::real getAverageRangeCollectionValue() const;

			/*! \brief Get a Range Collection representing the average
			of the heights of the elments in this. */
			virtual RangeCollectionHist getAverageRangeCollection() const;

			/*! \brief Reduce the rangeCollection for this node 
			to a single value representing the heights of the 
			elements originally represented in this's rangeCollection. 
			* 
			* Non-recursive: operates on this node only, not descendants.*/
			virtual void nodeRangeCollectionReduce();
			
			
			/*! \brief Multiply the rangeCollection for this node 
			by a scalar value.
			* 
			* Non-recursive: operates on this node only, not descendants.*/
			virtual void nodeScalarMult(cxsc::real& mult);
			
			/*! \brief Multiply the rangeCollections for this node 
			and its descendants by a scalar value.*/
			virtual void selfScalarMult(cxsc::real& mult);
			
			/*! \brief Divide the rangeCollection for this node 
			by a scalar value.
			* 
			* Non-recursive: operates on this node only, not descendants.*/
			virtual void nodeScalarDiv(cxsc::real& div);
			
			/*! \brief Divide the rangeCollections for this node 
			and its descendants by a scalar value.*/
			virtual void selfScalarDiv(cxsc::real& div);
			
			//internal versions
			/*! \brief Non-public version of normalise().*/
			virtual void _normalise();
			
			/*! \brief Normalise this using given normalising constant.
		 
			 \param normalisingConstant to normalise this and its descendants.			 */
			virtual void _normalise(const cxsc::real& normalisingConstant);
			
			/*! \brief Non-public version of marginalisation to
			have only the required dimensions.*/
			virtual void _start_marginalise(
					const std::vector<int>& reqDims);
			
			/*! \brief Non-public version of marginalisation to
			take out the unwanted dimensions.*/
			virtual void _marginalise(const std::vector<int>& outDims);
			
			/*! \brief Non-public version of average.*/
			virtual void _average();

			/*! Accumulates absolute value of differences
			of leaf node areas to the leaf node area of the average
			over the collation.
			
			One dot precision in vector for each histogram in collation.
			Accumulate results for over all the leaf 
			descendents of this.*/
			virtual VecDotPrec& getLeafNodeAbsDiffToAvAreaAccumulations(
									VecDotPrec& areaAcc) const;
									
			/*! Accumulates absolute value of leaf node areas.
			
			One dot precision in vector for each histogram in collation.
			Accumulate results for over all the leaf 
			descendents of this.*/
			virtual cxsc::dotprecision& getLeafNodeAbsAreaAccumulations(
									cxsc::dotprecision& areaAcc) const;
			
			
			
			// data members
			/*! \brief A collection of ranges.
			*/
			RangeCollectionHist rangeCollection;


		private:
		
			/*! \brief Internal version of concatention of collators.
			
			Make this's tree structure into non-minimal union of tree 
			structure of this and \a rhs, and append from rangeCollection
			of \a rhs onto this's rangeCollection.*/
			void _lazyCollationNonMinimalUnion(
				const CollatorSPnode * const rhs);

			/*! \brief Alternative internal version of addition of collators.
			
			Range collections are combined by 'parallel addition' rather
			than appending one to the other.  
			ie ith element of this rangeCollection += 
			ith element of \a rhs's rangeCollection.
			
			Makes this's tree structure into non-minimal union of tree 
			structure of this and \a rhs, rangeCollection into 
			paralled addition of rangeCollections of this and \a rhs. 
			\pre rangeCollections of this and \a rhs should be the
			same size.*/
			void _parallelAddNonMinimalUnion(
							const CollatorSPnode * const rhs);


			/*! \brief Reshape two trees so that both have the shape 
			that is the non-minimal union of the two.
			
			Results undefined if boxes are not the same.
			\pre \a lhs and \a rhs both non-NULL pointers.
			\post both \a lhs and \a rhs have the same tree shape and 
			that shape is the non-minimal union of the two.*/
			static void _reshapeTreesToUnion(CollatorSPnode * const lhs,
										CollatorSPnode * const rhs);

			/*! \brief Append the rangeCollections of other to the rangeCollection
			of this, for this node and its descendants.		 */
			void _lazyCollationRangeCollections(
				const CollatorSPnode * const other);
			
			// parallel add range collections from these two nodes down
			// assumes that this and other have the same shape from here down
			/*! \brief Combine the rangeCollections of this and another
			collator by 'parallel addition' 
			
			ie ith element of this rangeCollection += 
			ith element of \a rhs's rangeCollection. */
			void _parallelAddRangeCollections(
				const CollatorSPnode * const other);

			
			/*! @name Internal methods to get L1 distances.*/
			//@{
			VecDotPrec& _getL1distances(VecDotPrec& disL1,
										const CollatorSPnode& other) const;
											
			VecDotPrec& nodeL1Distances(VecDotPrec& disL1,
								cxsc::real other_h) const;
			//@}
			
			
			
	}; // end CollatorSPnode class



    // ----------------- non member tools functions ----------------------

/*! \brief Output operator for CollatorSPnodes.
*/
std::ostream & operator<<(std::ostream &os,
                        const CollatorSPnode& spn);

/*! \brief Comparison of CollatorSPnodes using total of range collection.
*/
bool nodeCompTotalRangeCollection(const CollatorSPnode * const lhs,
                            const CollatorSPnode * const rhs);


} // end namespace subpavings

/*! A specialisation of std::swap for NewCollatorNode types.*/
namespace std
{
	template <>
	void swap (subpavings::CollatorSPnode & s1, 
			subpavings::CollatorSPnode & s2);
}

#endif

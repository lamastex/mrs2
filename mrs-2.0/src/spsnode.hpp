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
\brief SPSnode (StatsSubPaving) declarations

*/

#ifndef ___SPSNODE_HPP__
#define ___SPSNODE_HPP__

#include "spnode.hpp" 

#include "MCMCPartitionGenerator.hpp"

#include <utility> // for pair
#include <stack> // for pair

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    //! Forward class declarations
    class SplitDecisionObj;
	class RealMappedSPnode;

    /*! \brief A derived class based on SPnode for processing sample data.

    The base class SPnode is a node in the representation of a regular
    subpaving as a binary tree.  A node represents a box (interval
    vector). SPnodes are linked together to form the tree.  The initial
    box of the subpaving is the box represented by the root node of the
    tree.  A box which has been split will be represented as node with 1
    or 2 children.

    A subpaving of [<b>x</b>] (union of non-overlapping subboxes of
    [<b>x</b>]) is represented by the leaves (degenerate/child-less)
    nodes in the tree.

    The SPSnode class has additional data members for statistical data
    analysis.  The SPSnode class is used to form a regular subpaving
    representing containers of sample data where some data-related
    criteria is used to determine when the subpaving should be bisected.
    For example, keeping the number of data points associated with a node
    below a specified maximum -- "a maximally statistically equivalent
    blocks" criterion.

    Data points are of type csxc::rvector.

    Leaves of the SPSnode class have data associated with them in the form
    of pointers to some big collection of sample data.  If an SPSnode is
    bisected the data associated with it descends to its children, so that
    only leaf SPSnodes have data associated with them.  However,
    "recursively-computable statistical summaries", such as, count, sum,
    etc, of the data which would be contained in the box an SPSnode
    represents are kept for all SPSnodes and continue to be updated
    when the SPSnode has children and data reaching the node is passed on
    to be finally associated with a leaf.

    For an algebraic statistical formalisation of Fisher's ideas on
    recursively computable statistics see "Notions of Sufficiency" by
    S.L. Lauritzen, Contributed Paper, 44th Session of the International
    Statistical Institute, Madrid, Spain, September 12th--22nd, 1983.

    By default, all recursively computable statistics provided are maintained
    in each SPSnode.  However, since this uses memory and is not always needed,
    an SPSnode can be constructed to only maintain count statistics.

    When a box is split along some dimension d, only the right (upper)
    subbox is closed on the actual split value (the midpoint of the
    interval which is the dth element in the interval vector/box to be
    bisected) whereas the left (lower) box has an open interval on this
    value.   Thus [1 5] is subdivided to [1 3) and [3 5].  The parent
    SPnode class can ignore this point but it must be addressed if we
    are to be able to decide which subbox (left or  right) a datapoint
    sitting exactly on the split value in the split dimension should
    descend to.  This is an extension of the empirical distribution to
    a partition of the root box with the leaves of the subpaving.

    The class also needs to know which dimension the box represented by
    a node was split on and what the split value was, so that data reaching
    the node can be correctly (ie, in cognisance of the open and closed
    intervals described above) associated with either the right or left
    child.  These values therefore become data members of the class,
    with default values for leaf nodes.
    */
    class SPSnode : public SPnode {

	public:
        /*! \brief Default constructor.
        */
        SPSnode();

        /*! \brief Initialised constructor.

		Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.
        
        \param v interval vector for the box this represents.
		\param cntOnly an indicator for whether all available stats
        are maintained (cntOnly = false) or just counts (cntOnly = true).
		\pre v should be a proper ivector, with length >= 1.
		\post this will have box and countsOnly as specified.
        */
        SPSnode(const ivector& v, bool cntOnly);

        /*! \brief Initialised constructor.
		
		Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

        \param v interval vector for the box this represents.
		\pre v should be a proper ivector, with length >= 1.
		\post this will have box as specified and
		countsOnly = false.        */
        explicit SPSnode(const ivector& v);


        /*! \brief Initialised constructor.

        Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

        \param v interval vector for the box this represents.
		\param max an indication of the maximum number of datapoints
		that may eventually be associated with a leaf.  This is used for 
		efficiency only and does not constitute a maximum count for a
		node.
		\param cntOnly an indicator for whether all available stats
        are maintained (cntOnly = false) or just counts (cntOnly = true).
		\pre v should be a proper ivector, with length >= 1.
		\post this will have box and countsOnly as specified.*/
        SPSnode(const ivector& v, size_t max, bool cntOnly);

        /*! \brief Initialised constructor.

        Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

        \param v interval vector for the box this represents.
		\param max an indication of the maximum number of datapoints
		that may eventually be associated with a leaf.  This is used for 
		efficiency only and does not constitute a maximum count for a
		node.
		\pre v should be a proper ivector, with length >= 1.
		\post this will have box as specified and
		countsOnly = false.        */
		SPSnode(const ivector& v, size_t max);

        /*! \brief Initialised constructor.

		Throws a MalconstructedBox_Error if \a lb 
		has a box with no practical dimensions.

        \param lb a labeled box which provides the box
		to be used for this.
		\param cntOnly an indicator for whether all available stats
        are maintained (cntOnly = false) or just counts (cntOnly = true).
		\pre lb should represent a proper box.
		\post this will have box and countsOnly as specified. */
        explicit SPSnode(const LabBox& lb, bool cntOnly = false);

        /*! \brief Initialised constructor.

        Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

        \param lb a labeled box which provides the box
		to be used for this.
		\param max an indication of the maximum number of datapoints
		that may eventually be associated with a leaf.  This is used for 
		efficiency only and does not constitute a maximum count for a
		node.
		\param cntOnly an indicator for whether all available stats
        are maintained (cntOnly = false) or just counts (cntOnly = true).
		\pre lb should represent a proper box.
		\post this will have box and countsOnly as specified. */
        SPSnode(const LabBox& lb, size_t max, bool cntOnly = false);


        /*! \brief Copy constructor.
         * 
         * Note that if the node to be copied is a leaf, and therefore
         * has a container of iterators to some data collection, then
         * the node created by a copy constructor will simply copy that
         * container of iterators - ie, both the original and copy
         * are associated with the same data collection.  If the node
         * is contained within by some other object that also manages
         * the data collection, it is advisable to ensure that copying
         * that other object creates a different data collection and 
         * that the leaf nodes of the copied subpaving tree reference
         * that new data collection. 
         
		 \param other the node to be copied to construct this. 
        */
        SPSnode(const SPSnode& other);

        // Use base class destructor

        /*! \brief Copy assignment operator.
        */
        SPSnode& operator=(SPSnode rhs);

        /*! \brief Accessor for the counter.
        */
        size_t getCounter() const;

        /*! \brief Get the split dimension.
        */
        virtual int getSplitDim() const;

        /*! \brief Get the split value.
        */
        virtual real getSplitValue() const;

        /*! \brief Accessor for the countsOnly value.
        */
        real getCountsOnly() const;

        /*! \brief Accessor for the node's data collection.

        Returns a copy of the node's collection of iterators to the
        big data set.
        */
        NodeData getData() const;

				/*
		 * These accessor methods shadow equivalent methods in the base
		class.  Thus the method used is determined  at compile time,
		not run time as would be the case if virtual methods were used.
		Because the pointers to parents and children are part of the
		base class definition, the methods have to cast the base class
		form to the derived class form in order for the pointer
		returned to be able to be used with derived class members.*/
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
        SPSnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Returns a pointer to leftChild node.        */
        SPSnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

        Returns a pointer to rightChild node.       */
        SPSnode* getRightChild() const;
        
		//@}
		
		/*! \brief Check tree rooted at this is legal with respect
		to isSplittableNode().
		* 
		'Legal' means that all non-leaf nodes in the tree are 
		splittable, ie return isSplittableNode() = true; 
        
		\return true if all non-leaf nodes in the tree rooted at
		this are splittable, false if any non-leaf node in the tree 
		rooted at this is not splittable.*/
        bool checkTreeStateLegal() const ;
		
		/*! \brief Check tree rooted at this is legal with respect
		to isSplittableNode(size_t minChildPoints, double minVol).
		 
		'Legal' means that all non-leaf nodes in the tree are 
		splittable, ie return 
		isSplittableNode(size_t minChildPoints, double minVol) = true; 
        
		\param minChildPoints is the minimum number of points that 
		there would be
		in the children if the node were to be split.
		\param minVol is the minimum node volume to be tested for.
		\return true if all non-leaf nodes in the tree rooted at
		this are splittable with respect to \a minChildPoints and
		\a minVol, false if any non-leaf node in the tree 
		rooted at this is not splittable according to these criteria.
		eg, a tree with a cherry node with volume < 2*minVol is not legal.*/
        bool checkTreeStateLegal(size_t minChildPoints, 
								double minVol) const;
		
		/*! \brief Check tree rooted at this is legal with respect
		to isSplittableNode(size_t minChildPoints).
		 
		'Legal' means that all non-leaf nodes in the tree are 
		splittable, ie return 
		isSplittableNode(size_T minChildPoints) = true; 
        
		\param minChildPoints is the minimum number of points that
		there would be
		in the children if the node were to be split.
		\return true if all non-leaf nodes in the tree rooted at
		this are splittable with respect to \a minChildPoints,
		false otherwise.*/
        bool checkTreeStateLegal(size_t minChildPoints);

		
		/*! Return boolean to indicate if node is splittable.
		* 
		A node is splittable if the node volume is >= 2 * cxsc::MinReal (the
		smallest representable real number).*/
		virtual bool isSplittableNode() const;
		
		/*! \brief Method to check whether a node is splittable.

		Decides whether a node is splittable based on checking volume 
		and number of points that would result in child nodes on split.

		Node must satisfy the basic isSplittableNode() test \b and
		volume must be >=2*minVol to split \b and if 
		\a minChildPoints > 0, then
		<ul>
		<li>either the node must have at least minChildPoints and all the points go
			to one of the children (the other getting none)</li>
		<li>or the smallest number of points which would go to the either of the
		prospective new children must be >= minChildPoints</li>
		</ul>

		\note if \a minVol is < minimum for passing 
		isSplittableNode() then node will fail test even if node
		volume > minVol, ie isSplittableNode() overrides all other 
		criteria.

		Thus in general the method will only return true if the given node satisfies
		both the minVol test and, if it were to be split, both children would have
		at least minChildPoints data points, \b but if all the data points would go
		to \b one child (none to the other), this will also satisfy the
		minChildPoints test.
		
		If the node has already been split, the test will use the actual numbers
		of points in the children; if the node is a leaf (ie not split) then
		the test will consider the number of points that would go to the 
		each child if it were to be split.

		\param minChildPoints is the minimum number of points that 
		there would be
		in the children if the node were to be split.
		\param minVol is the minimum node volume to be tested for.  This
		is the minimum volume that each \b child is permitted to have, i.e.
		to be splittable the node volume itself must be >= 2*minVol.
		\return true if has been a test conditions satisfied, false otherwise.
		*/
    	bool isSplittableNode(size_t minChildPoints, 
								double minVol) const;
		
		/*! \brief Method to check whether a node is splittable.

		Decides whether a node is splittable based 
		number of points that would result in child nodes on split.

		Node must satisfy the basic isSplittableNode() test 
		\b and if \a minChildPoints > 0, then
		<ul>
		<li>either the node must have at least minChildPoints and all the points go
			to one of the children (the other getting none)</li>
		<li>or the smallest number of points which would go to the either of the
		prospective new children must be >= minChildPoints</li>
		</ul>

		Thus in general the method will only return true if, 
		if it were to be split, both children would have
		at least minChildPoints data points, \b but if all the data points would go
		to \b one child (none to the other), this will also satisfy the
		minChildPoints test.
		
		If the node has already been split, the test will use the actual numbers
		of points in the children; if the node is a leaf (ie not split) then
		the test will consider the number of points that would go to the 
		each child if it were to be split.

		\param minChildPoints is the minimum number of points that there would be
		in the children if the node were to be split.
		\return true if has been a test conditions satisfied, false otherwise.
		*/
		bool isSplittableNode(size_t minChildPoints) const;

        /*! \brief The count the left child would have if this node was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the left child
        if the node were to be split.
		
		Note that the left child's interval on the split dimension would
		be an open interval. 
        */
        size_t getLeftCountIfSplit() const;

		/*! \brief The count the right child would have if this node was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the right child
        if the node were to be split.
		
		Note that the left child's interval on the split dimension would
		be a closed interval. 
        */
        size_t getRightCountIfSplit() const;

        
        /*! \brief Smallest number of points in either child if  this was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the left and right child
        if the node were to be split.
        */
        size_t getMinChildCountIfSplit() const;
		
		size_t getMinChildCountIfSplitNEW() const;

        /*! \brief return a container of counts for prospective grandchildren.

        Should be called only on leaf nodes.

        returns an indexable container of the number of points the prospective
        children of each prospective child (ie all four prospective
        grandchildren) would be associated with, indexed like this
        [0] = left child's left child count, [1] = left child's rght child count,
        [2] = rght child's left child count, [3] = rght child's rght child count,
        \param grandchildCounts a reference to a container to be filled with
        the prospective grandchild counts
        \return grandchildCounts filled with the prospective grandchild counts.
        */
        Size_tVec& getChildrensLeftAndRightCountsIfSplit
                    (Size_tVec& grandchildCounts) const;
		
		Size_tVec& getChildrensLeftAndRightCountsIfSplitNEW
                    (Size_tVec& grandchildCounts) const;

        /*! Fills in container of leaf counts, left to right.

        Traverses the leaves left to right, puts the leaf counts into container.

        \param counts is reference to the container to fill in.
        \return the container filled in with leaf counts.
        */
        Size_tVec& getLeafNodeCounts(Size_tVec& counts) const;


        /*! @name Get a container of all descendent leaf nodes.

        Will contain just this if this is a leaf.

        \param leaves a reference to a container of node pointers to fill in.
        \return a reference to the container \a leaves 
		filled with pointers to leaf nodes.   */
		
		//@{
		
		/*! \brief Returns pointers to non-const nodes.
		 * \todo This is bad bad bad and should not happen.  Change when 
		need for it for testing is over... 	*/ 
        SPSnodePtrs& getLeaves(SPSnodePtrs& leaves);

		/*! \brief Returns pointers to const nodes.*/ 
		SPSnodeConstPtrs& getConstLeaves(SPSnodeConstPtrs& leaves) const;
		
		//@}
		
        /*! @name Get a container of all sub-leaf nodes.

        Sub-leaf nodes (aka 'cherry nodes) have at least one child
		but any child must be a leaf,
        ie sub-leaves are the parents only of leaf nodes.

		Will just contain this if this is a sub-leaf.

        \param subleaves a reference to a container of node pointers to fill in.
        \return a reference to the container \a subleaves 
			filled with pointers to sub-leaf nodes.   		  */
		
		//@{
		
		/*! \brief Returns pointers to non-const nodes.
		\todo This is bad bad bad and should not happen.  Change when 
		need for it for testing is over... */ 
        SPSnodePtrs& getSubLeaves(SPSnodePtrs& subleaves);

		/*! \brief Returns pointers to const nodes.*/ 
		SPSnodeConstPtrs& getConstSubLeaves(SPSnodeConstPtrs& subleaves) const;
		
		//@}
		
		SPSnodePtrs& getLeavesInIntersection(
						const SPnode * const spn,
						SPSnodePtrs& leaves);
		
		SPSnodeConstPtrs& getConstLeavesInIntersection(
						const SPnode * const spn,
						SPSnodeConstPtrs& leaves) const;
		
		SPSnodePtrs& getSubLeavesInIntersection(
						const SPnode * const spn,
						SPSnodePtrs& subleaves);
		
		SPSnodeConstPtrs& getConstSubLeavesInIntersection(
						const SPnode * const spn,
						SPSnodeConstPtrs& subleaves) const;						
				
		/*! \brief The count in the node's ultimate ancestor root.
        */
        size_t getRootCounter() const;

        /*! \brief The count divided by the node volume.
        */
        double getCountOverVolume() const;

        /*! \brief Get the sample mean.

        \return If means are held (see countsOnly) and there 
		is at least one data point (see getCounter()) then return
		the sample mean. Otherwise return an rvector of 
		cxsc::SignalingNaN values.
        */
        rvector getMean() const;

        /** @name Get the sample variance-covariance matrix.

        This calculates the sample variance-covariance matrix for the
		sample and returns it as a d*d-dimensional vector of reals, 
		where d is the dimension of the data.

        cov(i,j) = [sumproduct(i,j)-sum(i)xsum(j)/counter]/(counter-1)
		
		cov(i,j) is at index i*d+j in the returned RealVec (indices 
		from 0 to d*d-1, where d is the dimension of the data).

        \return If variance-covariances are held (see countsOnly) 
		and there 
		are at least two data points (see getCounter()) then return
		the sample variance-covariance matrix as a vector.
		Otherwise return a RealVec of cxsc::SignalingNaN values.
        */
        
		//@{

        RealVec getVarCovar() const;

        RealVec& getVarCovar(RealVec& varCovar) const;

        //@}

		/*! \brief Get summary information on non-empty
		leaf box numbers and volumes.

		\return A pair, where the first value in the pair is the
		total number of non-empty leaf boxes in the tree rooted at this, and 
		the second value is the proportion of the total volume
		of this that is in the non-empty leaf boxes in 
		the tree rooted at this.*/
		std::pair<size_t, cxsc::real> getNonEmptyBoxSummary() const;

		/*! \brief Get the sum of the count over volume in the leaf nodes.
        */
        real getSumLeafCountOverVol() const;

        /*! \brief Get the count of the leaf with the smallest count.

        Returns the count in the smallest (by count) leaf node.
        */
        size_t getSmallestLeafCount() const;

        /*! \brief Get the count in the leaf with the smallest count.

        Returns the count of the largest (by count) leaf node.
        */
        size_t getLargestLeafCount() const;


        /*! \brief Get this node's contribution to loglikelihood.

        A leaf node's contribution to the log likelihood of overall state given
        the data is (count in leaf x ln(count in leaf / (n * vol of leaf)))
        where n is the total number of data points in the histogram.
		
		The loglikelihood of an entire tree of nodes comes only from the
		contributions of the current leaf nodes, but this method
		will not throw an exception if called on a non-leaf node. This
		allows 'what-if' type calculations to be carried out (eg what
		if the node reabsorbed its children)...
		
        \param n the value to use for scaling, the total number of points.
        \return counter * (ln(counter/(n*volume))) for this node.
        */
        virtual real getLogLik(const size_t n) const;

        /*! \brief Get change in log likelihood on split of this node.

        log likelihood is sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

        The split change loglikelihood only makes sense for a leaf
		node, but the calculation can be carried out for a non-leaf
		and so this method will not throw an exception if called on 
		a non-leaf node.
		
        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were split.
        */
        virtual dotprecision getSplitChangeLogLik() const;

		/*! \brief Get the unscaled log likihood for the tree
		rooted at this.
		 
		The unscaled log likelihood is sum over leaves of ( count in leaf
        x ( ln(count in leaf) + depth of leaf * ln2) )
        
		To get the scaled log likelihood we would add
		ln (1/(n x vol)^n ) = -n x ln (n x vol) where n is 
		the counter for this and vol is the volume of the box associated
		with this.  If we are taking ratios of likelihoods 
		for trees rooted at nodes with the same vol and n (eg different
		states of a tree rooted at this) we can 
		often ignore this scaling factor.
		
        \return the unscaled log likihood for the tree
		rooted at this.
        \pre this has a box.*/
        virtual cxsc::real getUnscaledTreeLogLik() const;

		/*! \brief Get change in log likelihood on merge of this' leaf chidren.

        log likelihood is sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

        The merge change only makes sense for a node 
		that is not already a leaf: this method throws an 
		UnfulfillableRequest_Error if called on a leaf node.
		
        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node's leaf chilren
        were merged back into this.
        */
        virtual dotprecision getMergeChangeLogLik() const;
		
		
        /*! \brief Get best change in EMP under COPERR from splitting any leaf.

        \param n the value to use for scaling, the total number of points.
        \return best (lowest or most negative) scaled change in COPERR from
        splitting any of the leaves.
        */

        dotprecision getBestSplitChangeEMPCOPERR(const size_t n) const;


        /*! \brief Get best change in EMP under AIC from splitting any leaf.

        n, the value to use for scaling, cancels out of change.
		
        \return best (lowest or most negative) scaled change in AIC from
        splitting any of the leaves.
        */
        dotprecision getBestSplitChangeEMPAIC() const;

        /*! \brief Get best change in EMP under COPERR from merging any subleaf.
		
		Finding the best merge change only makes sense for a node 
		that is not already a leaf: this method throws an 
		UnfulfillableRequest_Error if called on a leaf node.

        \param n the value to use for scaling, the total number of points.
        \return best (lowest or most negative) scaled change in COPERR from
        merging any of the subleaves.
        */
        dotprecision getBestMergeChangeEMPCOPERR(const size_t n) const;

        /*! \brief Get best change in EMP under AIC from merging any subleaf.

        n, the value to use for scaling, cancels out of change.

        Finding the best merge change only makes sense for a node 
		that is not already a leaf: this method throws an 
		UnfulfillableRequest_Error if called on a leaf node.

		\return best (lowest or most negative) scaled change in AIC from
        merging any of the subleaves.
        */
        dotprecision getBestMergeChangeEMPAIC() const;

        /*! \brief Get this node's scaled contribution to EMP under COPERR.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / volume of leaf)
        where n is the total number of data points in the histogram
        
		The EMP of an entire tree of nodes comes only from the
		contributions of the current leaf nodes, but this method
		will not throw an exception if called on a non-leaf node. This
		allows 'what-if' type calculations to be carried out (eg what
		if the node reabsorbed its children)...
        
		\param n the value to use for scaling, the total number of points.
        \return -counter^2/(n*volume) for this node.
        */
        real getEMPContributionCOPERR(const size_t n) const;

        /*! \brief Get this node's scaled contribution to EMP under AIC.

        Under AIC, EMP is -1 x sum over leaves of (counts in leaf
        x ln(count in leaf / (n * vol of leaf)))
        where n is the total number of data points in the histogram.
        And this is -1 * the node's contribution to the loglikelihood of the
        data given the current state.

		The EMP of an entire tree of nodes comes only from the
		contributions of the current leaf nodes, but this method
		will not throw an exception if called on a non-leaf node. This
		allows 'what-if' type calculations to be carried out (eg what
		if the node reabsorbed its children)...
        
		\param n the value to use for scaling, the total number of points.
        \return counter * (-ln(counter/(n*volume))) for this node.
        */
        real getEMPContributionAIC(const size_t n) const;

        /*! \brief Get scaled change in sum term in EMP under COPERR on split.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / (n * volume of leaf))
        where n is the total number of data points in the histogram.

		The split change only makes sense for a leaf
		node, but the calculation can be carried out for a non-leaf
		and so this method will not throw an exception if called on 
		a non-leaf node.

        \param n the value to use for scaling, the total number of points.
        \return the change in the sum over leaves of
        (counts in leaf squared over (n * volume of leaf)
        which would result if this node were split.
        */
        dotprecision getSplitChangeEMPCOPERR(const size_t n) const;

        /*! \brief Get change in sum term in EMP under AIC on split.

        Under AIC, EMP is -1 x sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

		The split change only makes sense for a leaf
		node, but the calculation can be carried out for a non-leaf
		and so this method will not throw an exception if called on 
		a non-leaf node.

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were split.
        */
        dotprecision getSplitChangeEMPAIC() const;

        /*! \brief Get scaled change in sum term in EMP under COPERR on merge.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / (n * volume of leaf))
        where n is the total number of data points in the histogram.

        The merge change only makes sense for a node 
		that is not already a leaf: this method throws an 
		UnfulfillableRequest_Error if called on a leaf node.

        \param n the value to use for scaling, the total number of points.
        \return the change in the sum over leaves of
        (counts in leaf squared over (n * volume of leaf)
        which would result if this node were merged.
        */
        dotprecision getMergeChangeEMPCOPERR(const size_t n) const;

        /*! \brief Get change in sum term in EMP under AIC on merge.

        Under AIC, EMP is -1 x sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

        The merge change only makes sense for a node 
		that is not already a leaf: this method throws an 
		UnfulfillableRequest_Error if called on a leaf node.

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were merged.
        */
        dotprecision getMergeChangeEMPAIC() const;

		/*! \brief Reshape so that the tree rooted at this has shape that
		is the union of this shape and the shape of another tree.
		
		Throws a NoBox_Error if this has no box or if \a other has no box. 
		
		Throws an IncompatibleDimensions_Error if boxes of this and \a other
		are not the same.
		
		Throws an std::runtime_error if \a other has an illegal state
		(see checkTreeStateLegal()).
		* 
		\param other is the tree to make the union against.
		\pre This has a box and that box is identical to the box of \a other. 
		\post the tree rooted at this has shape that is the
		union of the shape of this before the operation and the shape of 
		\a other.  \a other is unchanged.      */
		void reshapeToUnion(const SPnode& other);
		
		/*! \brief Reshape so that the tree rooted at this has shape that
		is as close as possible to the union of this shape 
		and the shape of another tree.
		
		If \a other is more split than this, this will
		 not exactly follow the shape of \a other if the resulting
		 nodes would not splittable according to 
		 isSplittableNode(size_t minChildPoints).  If any node 
		 cannot be split to follow the shape of \a other due to 
		 \a minChildPoints, a message will be printed to std::cerr.
				
		Throws a NoBox_Error if this has no box or if \a other has no box. 
		
		Throws an IncompatibleDimensions_Error if boxes of this and \a other
		are not the same.
		
		Throws an std::runtime_error if \a other has an illegal state
		(see checkTreeStateLegal()).
		 
		\param other is the tree to make the union against.
		\param minChildPoints is the minumum child points to use
		to check if this can be split in order to follow \a other.
		\pre This has a box and that box is identical to the box of \a other. 
		\post the tree rooted at this has shape that is the
		union of the shape of this before the operation and the shape of 
		\a other.  \a other is unchanged.      */
		void reshapeToUnion(const SPnode& other,
						size_t minChildPoints);

		
        /*! \brief Output details of a specific node.

        This is intended for console output or output
        to a mixed alpha and numeric file.
		\param os is the stream to send to
        */
        virtual std::ostream& nodePrint(std::ostream &os) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Includes scaled COPERR and AIC contributions and changes if split.
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leavesOutputTabsWithEMPs(const size_t bigN,
                        std::ostream &os, int prec = 5) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Includes output for the height of histogram bins for a normalised
        histogram based on tree with total number of data points of this root.
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leavesOutputTabsWithHistHeight(std::ostream &os,
                                    int prec = 5) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

         Includes output for the height of histogram bins for a normalised
        histogram based on tree with total number of data points of this root.
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leavesOutputTabsWithHistHeight(const size_t bigN,
                        std::ostream &os, int prec = 5) const;


        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabsWithHistHeightAndEMPs() to output information
        for each leaf node.
         Includes output for the height of histogram bins for a normalised
        histogram based on tree with total number of data points n.
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leavesOutputTabsWithHistHeightAndEMPs(const size_t bigN,
                            std::ostream &os, int prec = 5) const;


        /*! \brief Get scaled  EMP sum under COPERR for tree rooted at this.

        \param n the total number of points in the histogram, for scaling
        */
        dotprecision getEMPSumCOPERR(const size_t n) const;

        /*! \brief Get the unscaled EMP sum under AIC for tree rooted at this.

        \param n the total number of points in the histogram, for scaling
        */
        dotprecision getEMPSumAIC(const size_t n) const;

        /*! \brief Expand a leaf node.

        Expand a leaf node to have two children and pass
        data  down to the children with no further splitting.

        Equivalent to bisecting a box in a regular subpaving.
        Makes two new sibling child nodes of this.
        \param comp is the dimension on which to to bisect theBox.
        */
        virtual void nodeExpand(int comp);

        /*! \brief Expand a leaf node.

        Expand a leaf node to have two children and pass
        data down to the children, allowing for further splitting.

        \param boolTest is a reference to an object providing a function
        operator determining whether to split a node when a data point arrives.
        \param comp is the dimension on which to to bisect theBox.
        */
        virtual void nodeExpand(const SplitDecisionObj& boolTest, int comp);

        /*! \brief Expand a leaf node.

        Expand a leaf node to have two children and pass
        data  down to the children with no further splitting.

        Finds the splitting dimension.
        */
        virtual void nodeExpand();


        /*! \brief Expand a leaf node

        Expand the leaf node to have two children and pass
        data down to the children, allowing for further splitting.

        Finds the dimension to split on.

        \param boolTest is a reference to an object providing a function
        operator determining whether to split a node when a data point arrives.
        */
        virtual void nodeExpand(const SplitDecisionObj& boolTest);

        /*! \brief Reabsorbs both children of the node.

        Effectively reverses any split of the node.

        Data associated with the children is pushed back up to this
        and the splitDim and splitValue reset to leaf defaults.

        Works even if the children are not leaves.
        */
        virtual void nodeReabsorbChildren();
		
		/*! \brief Try to get tree rooted at this into shape at least as 
		deep as described by instruction \a reqDepths.
		
		Tries to ensure that this has at least shape of \a reqDepths, but 
		allows the tree rooted at this to also be more split at some 
		or all nodes.
		
		This can only follow the instruction if the nodes it is required
		to split are splittable nodes, with respect to both width 
		and \a minPoints (to be split a node 
		must return isSplittableNode(minPoints) = true).
		
		\param reqDepths a collection of leaf levels to try to split
		to.
		\param minPoints the minimum child points used to determine
		whether a node is splittable (defaults to 0). 
		\return true if this was able to follow the instruction, false
		otherwise (ie if instruction could not be followed because
		nodes were not splittable). 
		\pre this has a box.
		\pre this is a root node.
		\pre \a reqDepths describes a valid string (including is not empty).
		\post If returns true, the tree rooted at this has at least the
		shape described by \a reqDepths; if returns false, the tree rooted
		at this will have been split according to \a reqDepths, working
		from the right, until a node that cannot be split needs to 
		be split (all attempts to split will then have ceased).*/
		bool splitRootAtLeastToShapeSPS(std::vector < size_t > reqDepths,
								size_t minPoints = 0);
		
		/* docs */
		bool randomMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						bool saveInstructions = false);
						
		bool randomMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						bool saveInstructions = false);
		
		bool randomKnuthMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						bool saveInstructions,
						const std::string& failureLogFilename = "");
						
		bool randomKnuthMCMCSplitRootAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						bool saveInstructions,
						const std::string& failureLogFilename = "");
						
		
		/*! \brief Tries to insert data into the leaves of the tree
		 rooted at this.

        Passes down the tree, seeking the
        leaf node whose box contains the data point.
		
		If the box that this represents does not contain the data point
		then the return value is NULL.
		  
		If the box that this represents contains the data point, 
		statistics held by this are updated (subject to the value
		of countsOnly for this).  
		
		If this also a leaf node then 
		the data is associated with this leaf and the address
		of this becomes the return value.  If the SplitDecisionObj
		\a boolTest determines that this should expand following the
		insertion of the data, this is expanded.
		
		If the box that this
		represents contains the data point but this is not a leaf node,
		the return value will be the address of the ultimate leaf
		descendent of this that contains the data point. 
		
		\note For efficiency, there is no check that this has
		a box to check for containment against.  The results of trying this
		operation on a node with no box are undefined.  There is also no 
		check that the data point given has the same dimension as this;
		the results of trying this operation with a point with
		incompatible dimensions are undefined. 

        \param newItr an iterator to the data in big data collection.
        \param childInd an indicator for whether the current node is a
        treated as a left or right child or a root.
        \param boolTest is a reference to an object 
		determining whether to split a node when a data point arrives.
        \return a pointer to the node the data was 'inserted' into,
        (before it was split, if insertion of the data was followed
		by expansion of the node), or NULL if no insert.
		*/
        SPSnode* insertOneFind(BigDataItr newItr,
                            OPERATIONS_ON childInd,
                            const SplitDecisionObj& boolTest);


        /*! \brief Find the leaf node which would contain a data point.

        No data is actually inserted.
		
		Throws a NoBox_Error if this has no box.
		
		\note for efficiency, there is no check that the data point 
		given has the same dimension as this;  the results of trying this
		operation with a point with incompatible dimensions are undefined. 

        \param pt a data point.
        \param childInd an indicator for whether the current node is a
        treated as a left or right child or a root.  
		Defaults to ON_PARENT.
        \return a pointer to the node which would contain the point,
        or NULL if no node contains the point.
        \pre This should have a box to check for containment against.*/
	    const SPSnode* findContainingNode(const rvector& pt,
								OPERATIONS_ON childInd = ON_PARENT) const;
        
        /*! Makes the non-minimal union of the tree structures
		of this and another node.

        Discards the actual data.  The resulting structure returned
		has no data associated with it.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if this has no box.</li>
		<li>Throws a IncompatibleDimensions_Error if the node pointed
		to by \a rhs has a box but the dimensions
		and sizes of the boxes of this and the node pointed to by
		\a rhs are not the same (
		however, it is acceptable for \a rhs to be a NULL pointer
		or to point to a node having no box).</li>
		</ul>

        \param rhs a pointer to the root of SPSnode tree to make a union
		against, which can be a NULL pointer.  \a rhs is unchanged by
		the operation.
        \pre This must have a box.  If the node pointed to by \a rhs
		has a box, the size
		and dimensions of that box must be the same as those for 
		this's box.
		\post the tree rooted at this will have the the non-minimal 
		union of its original paving structure and that of the tree tooted
		at the node pointed to by \a rhs but no data:  all 
		counts will be 0, means and variance-covariances will be 
		undefined, and all data containers will be empty.
        */
		void unionTreeStructure(const SPSnode * const rhs);
		
		/*! Gets the L1 distance between this and another subpaving.

        The L1 distance is defined as the sum of the absolute values
		of the differences in 'area' represented by the leaf nodes
		of this and the other paving.  The 'area' represented
		by a leaf node is the proportion of the total count of the
		tree the leaf is part of that is in that leaf node.  
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NullSubpavingPointer_Error if the pointer
		to the other paving is NULL.</li>
		<li>Throws a NoBox_Error if either this or the node pointed to
		by \a other have no box.</li>
		<li>Throws a IncompatibleDimensions_Error if the dimensions
		and sizes of the boxes of this and the box of the node pointed
		to by \a other are not the same. 
		</ul>

        \param other a pointer to the root of SPSnode tree to calculate
		the L1 distance against.
        \pre \a other must be a non-NULL pointer.  Both this and the 
		node pointed to by \a other must have boxes and those
		boxes must be the same.
		\post this will be unchanged.       */
		cxsc::real getL1Distance(const SPSnode * const other) const;

		/*! \brief Clear all the data associated with the tree
		rooted at this.
		
		Calling this method resets counters to 0, discards values held
		to enable the computation of optional statistics, and discards
		data associated with root nodes.  The structure of the tree
		rooted at this is unchanged. 
		
		Throws a NonRootNode_Error if this is not a root node.
		
		\pre this is a root node.	
		\post For all nodes in the tree, getCounter() = 0 and the 
		optional statistcs are not available. Leaf nodes of the tree
		have no data associated with them.  The indicator countsOnly
		is unchanged for all nodes in the tree and the tree structure
		is unchanged (ie each node previously in the tree is still
		there).*/
		virtual void clearAllDataHeld() const;
		
		/*! \brief Set an indicator controlling whether this
		maintains all available statistics or not.
		
		If the node is changed from holding all available statistics to
		not holding all available statistics, the existing values
		held to enable calculation of the unwanted statistics are discarded.  
		
		If the node is changed from not holding all available statistics 
		to holding all available statistics, values enabling the 
		calculation of the optional statistics are recalculated.  
		
		Throws a NonRootNode_Error if this is not a root node.
		
		\param setTo the indicator for whether only counts are
		maintained (\a setTo = true) or whether all available
		statistics are maintained (\a setTo = false).
		\pre this is a root node.	
		\post countsOnly = setTo.  If \a setTo is true, values
		enabling the calculation of optional statistics are not
		held.  If \a setTo is false,  values
		enabling the calculation of optional statistics are
		held.*/
		virtual void setCountsOnly(bool setTo);
		
		/*! \brief Swap this and another node.
		
		Swaps all the data members of this with the other node. 
				
		\param spn a reference to the node to swap with
		\post this is identical,in terms of its data members, 
		to spn before the swap, and spn is
		identical to this after the swap.*/
		void swapSPS(SPSnode& spn); // throw()
		
		/*! \brief Get a string summary of this node's properties.
		
		Just this node, not this and its descendents.
		
		\return the string summary.		*/
		virtual std::string nodeStringSummary() const;
		
		/*! \brief Get a string summary of the properties of nodes
		in the tree rooted at this which can be used for checking
		and debugging.
		
		Includes the descendents of this node.  
		
		Shows the addresses of the values pointed to by the dataItrs
		associated with leaf nodes, rather than values themselves.
		
		\return the string summary.		*/
		virtual std::string doubleCheckStringSummary() const;

    protected:
	
		/*! Strip the data from this node and its children, seting 
		counters to 0 and discarding other statistics held, and 
		discarding data held by leaf nodes.  */
		virtual void stripData() const;
		
		/*! Discard the values held to enable the calculation of 
		optional statistics, and set countsOnly to true, in this
		and the children of this. */
		virtual void stripOptionalStatsOnly();

		/*! Recalculate the values held to enable the calculation of 
		optional statistics, and set countsOnly to false, in this
		and the children of this. */
		virtual void addOptionalStatsOnly();
		
		/*! \brief Send the data associated with this down to children.

        Children may then be resplit using boolTest.
        */
        virtual void splitData(const SplitDecisionObj& boolTest);


		/*! \brief Expand the node with no reallocation of data.

        Bisect box, make two new nodes (one for each half box) and graft
        onto this node provided that this node is a leaf. Equivalent
        to bisecting a box in a regular subpaving.
        */
        virtual void nodeExpansionOnly(int comp);
		
		/*! \brief Try to get tree rooted at this into shape at least as 
		deep as described by instruction \a reqDepths.
		
		Tries to ensures that this has at least shape of \a reqDepths, but 
		allows the tree rooted at this to also be more split at some 
		or all nodes.
		
		This can only follow the instruction if the nodes it is required
		to split are splittable nodes, with respect to both width 
		and \a minPoints (to be split a node 
		must return isSplittableNode(minPoints) = true).
		
		\param reqDepths a collection of leaf levels to try to split
		to.
		\param myDepth the depth this node is in the tree.
		\param minPoints the minimum child points used to determine
		whether this is splittable. 
		\return true if this was able to follow the instruction, false
		otherwise (ie if instruction could not be followed because
		nodes were not splittable). 
		\pre \a reqDepths is not empty*/
		bool _splitAtLeastToShapeSPS(
								std::vector < size_t >& reqDepths,
								size_t myDepth,
								size_t minPoints);

		bool _randomMCMCSplitAtLeastSPS(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						size_t myDepth,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints);
		
		
		
		int levelsUpToNextRight() const;
		
		bool _randomKnuthMCMCSplitAtLeastSPS(
						unsigned long int p,
						unsigned long int q,
						std::stack< size_t>& lastPair,
						const MCMCPartitionGenerator& partitioner,
						RealMappedSPnode* rmsp,
						size_t myDepth,
						dotprecision& nlogn, 
						int& ndepth,
						size_t minPoints,
						const std::string& failureLogFilename,
						bool across = false);

		
		/*! \brief Internal method to reshape this to a union.*/
		void _reshapeToUnion(const SPnode * const other);
		
		/*! \brief Internal method to reshape this to a union
		with a restriction of minChildPoints.
		
		\param minChildPoints is the minimum child points to allow 
		for a split (as usual, if one child gets all the points and the 
		number of points in this is >= minChildPoints, the split is 
		allowed).
		\param errorFileName is a file to which to log messages 
		about nodes which could not be split because of \a minChildPoints. 
		\return true if this was able to reshape to exactly the union
		of this and \a other, false if the extent to which descendents
		of this could split to mimic other was limted by \a minChildPoints.*/
		bool _reshapeToUnion(const SPnode * const other,
							size_t minChildPoints,
							const std::string& errorFilename);

		/*! \brief Internal method to accumulate the components of the 
		unscaled log likelihood for a tree rooted at this.*/
		virtual void _getUnscaledTreeLogLik(dotprecision& nlogn, 
										int& ndepth, int depth) const;

		
		/*! \brief Accumulate a sum of node counts divided by node volumes.

        \param sum a reference to an accumulation of the sum so far.
        \return a reference to the accumulation including the leaf nodes 
        of this.*/
		dotprecision& accumulateLeafCountOverVol(dotprecision& sum) const;

		/*! \brief Accumulate summary information on non-empty
		box numbers and volumes.

        \param nNonEmptyBoxes a reference to an accumulation the
		number of non-empty boxes so far.
		\param vNonEmptyBoxVolumes a reference to an accumulation the
		volume of non-empty boxes so far.*/
		void accumulateNonEmptyBoxSummary(size_t& nNonEmptyBoxes, 
								real& vNonEmptyBoxVolumes) const;

        /*! \brief Print the data in a specified format.

        Replaces the format that the cxsc::<< operator produces for
        vectors.  The format used here
        produces numeric tab-delimited data.  The format for an
        n-dimensional real vector data point is:

        rvector[1] [tab] . . . [tab] rvector[n]
        */
        virtual std::ostream& nodeDataPrint(std::ostream &os) const;

        /*! \brief Print the mean in a specified format.
        */
        virtual std::ostream& nodeMeanPrint(std::ostream &os) const;

        /*! \brief Print the variance-covariance in a specified format.
        */
        virtual std::ostream& nodeVarCovarPrint(std::ostream &os) const;

        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.  The format used here
        produces alpha-numeric tab-delimited data.  The format
        for an n-dimensional interval vector is:

        nodeName [tab] volume [tabl] counter [tab] 
		Inf(ivector[1]) [tab] Sup(ivector[1]) [tab] ...
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n])
		\param os is the stream to send to.        */
        virtual std::ostream& leafOutputTabs(std::ostream &os) const;

		/*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here produces 
		alpha-numeric tab-delimited data.

        The format for a d-dimensional interval vector is

        nodeName [tab] volume [tab] counter [tab] 
        scaled EMP contribution COPERR [tab]
        change in scaled EMP contribution COPERR if split [tab]
        scaled EMP contribution AIC [tab]
        change in scaled EMP contribution AIC if split [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[d]) [tab] Sup(ivector[d]
        \param bigN is total datapoints, used by the emps calculation
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leafOutputTabsWithEMPs(const size_t bigN,
                        std::ostream &os, int prec = 5) const;

            
        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Includes the height of a bin represented by this leaf node for a
        normalised histogram, ie counter/(node volume * total count in tree)

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.  The format used here
		produces alpha-numeric tab-delimited data.  The format
        for a d-dimensional interval vector is:

        nodeName [tab] volume [tab] counter [tab] counter/(volume*total count)
        [tab] Inf(ivector[1]) [tab] Sup(ivector[1]) [tab] ...
        [tab] Inf(ivector[d]) [tab] Sup(ivector[d])
        \param bigN is total data points, used by the height calculation
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leafOutputTabsWithHistHeight(
                                const size_t bigN, std::ostream &os,
                                int prec = 5) const;


        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here
		produces alpha-numeric tab-delimited data.

        Includes the height of a bin represented by this leaf node for a
        normalised histogram, ie counter/(node volume * total count in tree)

        The format for a d-dimensional interval vector is

        nodeName [tab] volume [tabl] counter [tab] counter/(volume*total count)
        scaled EMP contribution COPERR [tab]
        change in scaled EMP contribution COPERR if split [tab]
        scaled EMP contribution AIC [tab]
        change in scaled EMP contribution AIC if split [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[d]) [tab] Sup(ivector[d]
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        std::ostream& leafOutputTabsWithHistHeightAndEMPs(const size_t bigN,
                        std::ostream &os, int prec = 5) const;

		/*! \brief Clears the node's data collection.
		
		Note - only clears the data collection in this node,
		and does nothing to counters or stats - used when expanding
		a node.
        */
        virtual void clearNodeData() const;

		/*! \brief Recalculate summary statistics associated with node.

        Recalculates counter and sums (used for mean) and
        sumproducts (used for variance-covariance).
        */
        virtual void recalculateStats(const rvector& newdata) const;

        /*! \brief Recalculate summary statistics associated with node.

        Recalculates sums (used for mean).
        */
        virtual void recalculateSums(const rvector& newdata) const;

        /*! \brief Recalculate summary statistics associated with node.

        Recalculates sumproducts (used for variance-covariance).
        */
        virtual void recalculateSumProducts(const rvector& newdata) const;
		
		/*! \brief Recalculates optional summary statistics associated
		 with node from scratch, using all the data currently held.

        Recalculates sums and sumproducts.
		
		Uses only the data currently associated with the node in the
		the recalculation.
		* 
		Does not recalculate non-optional summaries, like counter, 
		because these should be right already
        */
        
		virtual void recalculateOptionalStatsAndGatherData(
											NodeData& container);

		
		virtual void recalculateOptionalStatsForData(
										const NodeData& nodedata) const;
		// data members
		/* theBox, parent, leftChild,
        rightChild and nodeName are inherited from base class */

        /*! \brief An indication of the maximum number of data points
        a node needs to carry.

        This is used for efficiency only to reserve vector space and a
        node can have more than this maximum number of data points
        associated with it.  Defaults to defaultMaxPts.
        */
        size_t spaceIndication;

        /*! @name mutable data members.
        These data members are mutable to allow them to be modified by
        const functions as data passed to or through the node.

        Only leaf nodes have data associated with them but the
        recursively computable statistics, such as counter and sum, are
        maintained for all nodes.  Thus when a data point is sent to
        the root node and progresses down the tree to find which
        leaf node it should be associated with, the counter is
        incremented and the data sum increased for each non-leaf node
        it passes through (ie, where it is contained in the box of
        that node but that node is not a leaf node so the box has been
        sub- divided and the datapoint continues on to one of the
        children).
        */
        
		//@{
        
		//! A counter for how many data points are covered by theBox
        mutable size_t counter;

        /*! \brief A container representing the sum of the data points
        covered by theBox.

        The sums are used for calculating the mean and also the sample
        variance-covariance matrix for the data associated with a node.

        cxsc::dotprecision accumulators are used to maintain the sum
        of the data in each dimension of the data because floating
        point arithmetic can result in inaccuracies during summation,
        especially in large boxes.

        We could use Kahan summation instead with a lot more work.
        Kahan summation relies on adding a number of points in a
        sequence and recovering data lost in one summation during the
        next one.  When we simply add two numbers, Kahan summation has
        no chance to recover the lost part.  We would have to implement
        this by keeping the lost part, say having a vector of pairs,
        and trying to re-add the lost part each time.
        See http://en.wikipedia.org/wiki/Kahan_summation_algorithm .
        The same may apply to arguments for using gsl_mean using a
        simpler reccurence relation.  Speed comparisons have not been
        performed on the three alternative possible implementations
        of the recursively computable sample sum or sample mean.
        However we get most relable and accurate sums using
        cxsc::dotprecision accumulators.
        */
        mutable VecDotPrec dpSums;

        /*! \brief A container representing the sumproduct matrix of the
        data points covered by theBox.

        The sumproducts matrix is used to obtain the sample
        variance-covariance matrix.

        The for n-dimensional data the sample variance-covariance matrix is
        an nxn matrix where the element in row i, column j is the sample
        covariance between the ith-dimension and jth-dimension of the data,
        which is [sumproduct(i,j)-sum(i)sum(j)/counter]/(counter-1).

        So by keeping the sum product and sums up to date, we can calculate a
        covariance on demand.

        The sumproducts can be thought of as an nxn matrix where the element
        in row i, column j is the sum over all the datapoints associated with
        that box of the products of the ith element and jth element in the
        datapoints.  ie for each datapoint, we take the product of the ith
        and jth elements and then sum the products over all the datapoints.

        Data points are rvectors so each element is a real, and the the
        accumulation (sum) of products of reals is implemented here with a
        dotprecision accumulator.

        The sumproduct matrix is stored here as a nxn element vector of
        dotprecision variables (where n is the dimensions of the rvectors or
        data points), using row-major order.

        Ie the m-th element (m = 0, . . . nxn-1) is in row floor(m/n)
        and column m-rowxn in the matrix configuration.

        Or, the sumproduct of elements i and j in an rvector,
        i,j = 0,...,n-1, is element m=(ixn+j) of the sumproducts
        vector.
        */
        mutable VecDotPrec dpSumProducts;

        /*! \brief A container for the association of data with a node.

        Data is associated with a node via this container of iterators.
        The iterators can, very loosely, in the sense in which they are
        used here, be thought of as pointers to a big data collection
        of all data points.  Only leaf nodes can have anything in this
        container.  However, not all leaf nodes will necessarily have
        something in this container:  the container will be empty if no
        data points are covered by the box represented by a leaf node.
        */
        mutable NodeData dataItrs;
        
		//@}

        /*! \brief Dimension the node's box has been split along.

        This is needed to accurately divide data between a node's left
        and right children.  Defaults to -1 if the node is a leaf.
        */
        int splitDim;

        /*! \brief The value, on split dimension, where node's box was split.

        This is needed to accurately divide data between a node's left
        and right children.  Defaults to 0.0 if the node is a leaf.
        */
        real splitValue;

        /*! \brief Determines the amount of statistical summary data in node.

        If true, only counts are maintained in the node.  If false, counts and
        other statistics (sums, sumproducts) are also maintained.

        Setting countsOnly to true can help to conserve memory when
        means (from sums) and covariances (from sumproducts) will not be needed.
        */
        bool countsOnly;
		

    private:
        
        /*! \brief To define the default maximum number of datapoints
        the node is expected to have associated with it.

        This is used for efficiency only to reserve vector space and a
        node can have more than this default maximum number of data
        points associated with it.
        */
        enum {defaultMaxPts = 1000};

        /*! \brief Gather data associated with a node and its descendents.

        \return a reference to a container of data associated with leaf
		children of the node.
        */
        NodeData& gatherData(NodeData& container) const;

		/*! Makes the non-minimal union of this and the tree
		rooted at the node pointed to by \a rhs. */					
		void unionNoData(const SPSnode * const rhs);

		
		/*! Calculate the L1 distance between this and another subpaving.*/
		cxsc::dotprecision& _getL1distance(cxsc::dotprecision& disL1,
							const SPSnode * const other,
							std::size_t thisBigN,
							std::size_t otherBigN) const;
		
		/*! Calculate the L1 distance between this and another node.*/								
		cxsc::dotprecision& nodeL1Distance(cxsc::dotprecision& disL1,
							std::size_t thisBigN, 
							std::size_t other_n, std::size_t otherBigN) const;
		

    };    // end of SPSnode class derived from SPnode class


} // end namespace subpavings

/*! A specialisation of std::swap for SPSnode types.*/
namespace std
{
	template <>
	void swap(subpavings::SPSnode & s1, 
			subpavings::SPSnode & s2); // throw ()
	
}

#endif


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

#ifndef ___SPSNODE_HPP__
#define ___SPSNODE_HPP__

// for spnode class
#include "spnode.hpp" // includes sptypes.hpp

/*! \file spsnode.hpp
\brief SPSnode (StatsSubPaving) declarations

*/


/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    //! Forward class declarations
    class SplitDecisionObj;

    /*! \brief StatsSubPaving is an alias for a pointer to an SPSnode.
    */
    typedef SPSnode* StatsSubPaving;


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
    private:
        /* theBox, dimension, label, parent, leftChild,
        rightChild and nodeName are inherited from base class */

        /*! \brief To define the default maximum number of datapoints
        the node is expected to have associated with it.

        This is used for efficiency only to reserve vector space and a
        node can have more than this default maximum number of data
        points associated with it.
        */
        enum {defaultMaxPts = 1000};

        /*! \brief Recalculate summary statistics associated with node.

        Recalculates counter and sums (used for mean) and
        sumproducts (used for variance-covariance).
        */
        void recalculateStats(rvector& newdata) const;

        /*! \brief Recalculate summary statistics associated with node.

        Recalculates sums (used for mean).
        */
        void recalculateSums(rvector& newdata) const;

        /*! \brief Recalculate summary statistics associated with node.

        Recalculates sumproducts (used for variance-covariance).
        */
        void recalculateSumProducts(rvector& newdata) const;

        /*! \brief Expand the node with no reallocation of data.

        Bisect box, make two new nodes (one for each half box) and graft
        onto this node provided that this node is a leaf. Equivalent
        to bisecting a box in a regular subpaving.
        */
        void nodeExpansionOnly(int comp);

        /*! \brief Send the data associated with this down to children.

        Children may then be resplit using boolTest.
        */
        void splitData(const SplitDecisionObj& boolTest);

        /*! \brief Print the data in a specified format.

        Replaces the format that the cxsc::<< operator produces for
        vectors.  The format used here includes the box label and
        produces numeric tab-delimited data.  The format for an
        n-dimensional real vector data point is:

        label [tab] rvector[1] [tab] . . . [tab] rvector[n]
        */
        std::ostream& nodeDataPrint(std::ostream &os) const;

        /*! \brief Print the mean in a specified format.
        */
        std::ostream& nodeMeanPrint(std::ostream &os) const;

        /*! \brief Print the variance-covariance in a specified format.
        */
        std::ostream& nodeVarCovarPrint(std::ostream &os) const;

        /*! \brief Return a reference to the node data.

        \return a reference to a container of data associated with the
        node <b> and its descendents </b>.
        */
        static NodeData& gatherData(NodeData& container, SPSnode * spn);


/*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.

        The format for a d-dimensional interval vector is

        nodeName [tab] counter [tab] volume [tab]
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
        std::ostream& leafOutputTabsWithEMPs(const size_t bigN,
                        std::ostream &os, const int prec = 5) const;


        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Includes the height of a bin represented by this leaf node for a
        normalised histogram, ie counter/(node volume * total count in tree)

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.  The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for a d-dimensional interval vector is:

        nodeName [tab] counter [tab] volume [tabl] counter/(volume*total count)
        [tab] Inf(ivector[1]) [tab] Sup(ivector[1]) [tab] ...
        [tab] Inf(ivector[d]) [tab] Sup(ivector[d])
        \param bigN is total data points, used by the height calculation
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        virtual std::ostream& leafOutputTabsWithHistHeight(
                                const size_t bigN, std::ostream &os,
                                const int prec = 5) const;


        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.

        Includes the height of a bin represented by this leaf node for a
        normalised histogram, ie counter/(node volume * total count in tree)

        The format for a d-dimensional interval vector is

        nodeName [tab] counter [tab] volume [tabl] counter/(volume*total count)
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
                        std::ostream &os, const int prec = 5) const;



        /*! \brief Set the splitDimension and SplitValue when children grafted.

        Sets the splitDimension and SplitValue for this node when children
        are grafted on.
        Called by nodeAdoptLeft or nodeAdoptRight.

        \post this has splitDimension and SplitValue corresponding to children.
        */
        void setSplits();


        /*! Makes the non-minimal union of nodes with no data.

        Calls itself recursively to adds two pavings together as the union of
        the two but discards the actual data.

        Does not rename the nodes from root downwards.

        \param lhs pointer to root of first SPSnode tree operand.
        \param rhs pointer to root of second SPSnode tree operand.
        \return a pointer to root of a new SPSnode tree whose leaves are the
        union of the leaves of lhs, rhs and which has no data.
        */
        static SPSnode* unionNoData(const SPSnode * const lhs,
                            const SPSnode * const rhs);


    protected:

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

        Setting countsOnly to false can help to conserve memory when
        means (from sums) and covariances (from sumproducts) will not be needed.
        */
        bool countsOnly;

		//src_trunk_0701
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
		//src_trunk_0701

    public:
        /*! \brief Default constructor.
        */
        explicit SPSnode();

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents
        and a value for countsOnly which controls whether all available stats
        are maintained (false) or just counts (true), and optionally
        initialised with a label which defaults to 0 if not provided.
        */
        explicit SPSnode(ivector& v, bool cntOnly, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents
        and optionally initialised with a label which defaults to 0
        if not provided.  The value for countsOnly will default to false
        (i.e., all stats maintained by default).
        */
        explicit SPSnode(ivector& v, int lab = 0);


        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents,
        a space indication, and a value for countsOnly which controls wether
        all available stats are maintained (false) or just counts (true), and
        optionally initialised with a label for the model which defaults
        to 0 if not provided.
        */
        explicit SPSnode(ivector& v, size_t max, bool cntOnly, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents,
        a space indication, and optionally initialised with a label for the
        model which defaults to 0 if not provided. The value for countsOnly
        will default to false (i.e., all stats maintained by default).
        */
        explicit SPSnode(ivector& v, size_t max, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with a LabBox for the labeled box it represents.  Also
        optionally initialised with a value for countsOnly, defaults to
        false (i.e., all stats maintained by default).
        */
        explicit SPSnode(LabBox& lb, bool cntOnly = false);

        /*! \brief Initialised constructor.

        Initialised with a LabBox for the labeled box it represents,
        and a space indication.  Also optionally initialised with a value for
        countsOnly, defaults to false (i.e., all stats maintained by default).
        */
        explicit SPSnode(LabBox& lb, size_t max, bool cntOnly = false);


        /*! \brief Copy constructor.
        */
        explicit SPSnode(const SPSnode& other);

        // Use base class destructor

        /*! \brief Copy assignment operator.
        */
        SPSnode& operator=(const SPSnode& rhs);

        static SPSnode* strippedConstructor(const SPSnode * const other);  
        
        /*! \brief Accessor for the counter.
        */
        size_t getCounter() const;
       

        /*! \brief Accessor for the split dimension.
        */
        int getSplitDim() const;

        /*! \brief Accessor for the split value.
        */
        real getSplitValue() const;

        /*! \brief Accessor for the countsOnly value.
        */
        real getCountsOnly() const;

        /*! \brief Accessor for the node's data collection.

        Returns a copy of the node's collection of iterators to the
        big data set.
        */
        NodeData getData() const;

        /*! \brief Clears the node's data collection.
        */
        void clearData() const;
		  
		  /*! \brief Clears the node's data collection and counter.
		  */
		  void makeEmptyNode();

        /*! @name Accessors for links between the nodes.
        These accessor methods shadow equivalent methods in the base
        class.  Thus the method used is determined  at compile time,
        not run time as would be the case if virtual methods were used.
        Because the pointers to parents and children are part of the
        base class definition, the methods have to cast the base class
        form to the derived class form in order for the pointer
        returned to be able to be used with derived class members.

        Note that pointers for parent, leftChild, and rightChild are
        not reference counted so there could potentially be problems
        with the use of returned pointers (for instance, being used to
        delete nodes). These pointers might be better implemented with
        boost::shared_ptr .
        */

        //@{
        /*! \brief Accessor for the parent of a node.

        Hides the base class version of this method.

        Returns a copy of the pointer to parent node.
        */
        SPSnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Hides the base class version of this method.

        Returns a copy of the pointer to leftChild node.
        */
        SPSnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

        Hides the base class version of this method.

        Returns a copy of the pointer to rightChild node.
        */
        SPSnode* getRightChild() const;
        //@}

		//--src_trunk_0701
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
		rooted at this is not splittable according to these criteria.*/
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
        //--src_trunk_0701


		//--src_trunk_0701
				
		/*! Return boolean to indicate if node is splittable.
		* 
		A node is splittable if the node volume is >= 2 * cxsc::MinReal (the
		smallest representable real number).*/
		virtual bool isSplittableNode() const;
		
		/*! \brief Method to check whether a node is splittable.

		Decides whether a node is splittable based on checking volume 
		and number of points that would result in child nodes on split.

		Node must satisfy the basic isSplittableNode() test \b and
		volume must be >=minVol to split \b and if 
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
		\param minVol is the minimum node volume to be tested for.
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
		//--src_trunk_0701

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
        //src_trunk_0701
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
                    
        /*! \brief Smallest volume of either child if this was split.

        Does not split the nodes, just returns the volume of the left and right
		  child if the node were to be split.
        */
        double getMinChildVolIfSplit() const;

        /*! Fills in container of leaf counts, left to right.

        Traverses the leaves left to right, puts the leaf counts into container.

        \param counts is reference to the container to fill in.
        \return the container filled in with leaf counts.
        */
        Size_tVec& getLeafNodeCounts(Size_tVec& counts) const;


        /*! \brief Return a reference to all descendent leaf nodes.

        Will be just this if this is a leaf.

        \return a reference to a container of node pointers.
        */
        SPSnodePtrs& getLeaves(SPSnodePtrs& leaves) const;

        /*! \brief Return a reference to all nodes.

        \return a reference to a container of node pointers.
        */
        SPSnodePtrs& getAllNodes(SPSnodePtrs& allNodes) const;

        /*! \brief Return a reference to all sub-leaf descendent nodes.

        Sub-leaf nodes have at least one child but any child must be a leaf,
        ie sub-leaves are the parents of leaf nodes.

        Will be just this if this is a subleaf.

        \return a reference to a container of node pointers.
        */
        SPSnodePtrs& getSubLeaves(SPSnodePtrs& subleaves) const;

        

        /*! \brief The count in the node's ultimate ancestor root.
        */
        size_t getRootCounter() const;


        /*! \brief Get the sample mean.

        This calculates the sample mean from the accumulators for the
            sums of data point elements.
        */
        rvector getMean() const;
        
        //gat41
        /*! \brief Get the uniform mean vector where each element is the midpoint 
         * of the coordinate.
         */
         rvector getUniformMean() const;
         
         //gat41
         /*! \brief Get the Chebyshev distance for the mean.
         */
         real getChebDistMean() const; 

			//gat41
			/*! \brief get the empirical mass of the node
			 */
			 double getEmpMass() const;
			 
			 //gat41
         /*! \brief Get the Chebyshev distance for the var-covar.
         */
         real getChebDistCovar() const; 
			
			//gat41
			/*! \brief Get the Bhattarchaya coefficient.
			*/
			real getHellingerDist() const;
			real getHellingerDist1D() const;
			 
        /** @name Get the sample variance-covariance matrix.

        This calculates the sample variance-covariance matrix from
        the accumulators for the sumproducts and sums of data point
        elements.

        cov(i,j) = [sumproduct(i,j)-sum(i)xsum(j)/counter]/(counter-1)

        \return a RealVec or reference to a RealVec
        representing the sample variance-covariance matrix in
        row-major order.
        */
        //@{

        RealVec getVarCovar() const;

        RealVec& getVarCovar(RealVec& varCovar) const;
	
			//gat41
			RealVec getUniformVarCovar() const;
			RealVec& getUniformVarCovar(RealVec& varCovar) const;
        //@}

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


        /*! \brief Get this leaf node's contribution to loglikelihood.

        A leaf node's contribution to the log likelihood of overall state given
        the data is (count in leaf x ln(count in leaf / (n * vol of leaf)))
        where n is the total number of data points in the histogram.

        Should only be called on leaf nodes.
        \param n the value to use for scaling, the total number of points.
        \return counter * (ln(counter/(n*volume))) for this node.
        */
        real getLogLik(const size_t n) const;

        /*! \brief Get change in log likelihood on split of this node.

        log likelihood is sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

        Should only be called on leaf nodes.

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were split.
        */
        dotprecision getSplitChangeLogLik() const;

        /*! \brief Get change in log likelihood on merge of this' leaf chidren.

        log likelihood is sum over leaves of (counts in leaf
        x ln(count in leaf / (n * volume of leaf)))
        where n is the total number of data points in the histogram.

        Should only be called on cherry leaf nodes (two leaf children).

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node's leaf chilren
        were merged back into this.
        */
        dotprecision getMergeChangeLogLik() const;


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

        \param n the value to use for scaling, the total number of points.
        \return best (lowest or most negative) scaled change in COPERR from
        merging any of the subleaves.
        */
        dotprecision getBestMergeChangeEMPCOPERR(const size_t n) const;

        /*! \brief Get best change in EMP under AIC from merging any subleaf.

        n, the value to use for scaling, cancels out of change.

        \return best (lowest or most negative) scaled change in AIC from
        merging any of the subleaves.
        */
        dotprecision getBestMergeChangeEMPAIC() const;

        /*! \brief Get this node's scaled contribution to EMP under COPERR.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / volume of leaf)
        where n is the total number of data points in the histogram
        Should only be called on leaf nodes
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

        Should only be called on leaf nodes.
        \param n the value to use for scaling, the total number of points.
        \return counter * (-ln(counter/(n*volume))) for this node.
        */
        real getEMPContributionAIC(const size_t n) const;

        /*! \brief Get scaled change in sum term in EMP under COPERR on split.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / (n * volume of leaf))
        where n is the total number of data points in the histogram.

        Should only be called on leaf nodes.

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

        Should only be called on leaf nodes.

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were split.
        */
        dotprecision getSplitChangeEMPAIC() const;

        /*! \brief Get scaled change in sum term in EMP under COPERR on merge.

        Under COPERR, EMP is -1/n^2 x sum over leaves of (counts in leaf
        squared / (n * volume of leaf))
        where n is the total number of data points in the histogram.

        Should only be called on subleaf nodes (two leaf children).

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

        Should only be called on subleaf nodes (two leaf children)

        \return the change in the sum over leaves of (counts in leaf squared
        over volume of leaf) which would result if this node were merged.
        */
        dotprecision getMergeChangeEMPAIC() const;


		//src_trunk_0701
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
		//--src_trunk_0701

        /*! \brief Output details of a specific node.

        This is intended for console output or output
        to a mixed alpha and numeric file.
        */
        virtual std::ostream& nodePrint(std::ostream &os) const;

        /*! \brief Output for a node in a binary tree, tab-delimited.

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.  The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for an n-dimensional interval vector is:

        label [tab] counter [tab] volume [tabl] Inf(ivector[1]) [tab]
        Sup(ivector[1]) [tab] ...
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n])
        */
        virtual std::ostream& leafOutputTabs(std::ostream &os) const;


        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabs() to output information
        for each leaf node.
        */
        virtual std::ostream& leavesOutputTabs(std::ostream &os) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabs() to output information
        for each leaf node.
        Includes scaled COPERR and AIC contributions and changes if split.
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        std::ostream& leavesOutputTabsWithEMPs(const size_t bigN,
                        std::ostream &os, const int prec = 5) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabsWithHistHeight() to output information
        for each leaf node.
        Includes output for the height of histogram bins for a normalised
        histogram based on tree with total number of data points of this root.
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        std::ostream& leavesOutputTabsWithHistHeight(std::ostream &os,
                                    const int prec) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabsWithHistHeight() to output information
        for each leaf node.
        Includes output for the height of histogram bins for a normalised
        histogram based on tree with total number of data points of this root.
        \param bigN is total data points, used by emps and height calculations
        \param os is the stream to send to
        \param prec is the precision used for printing, defaulting to 5
        */
        std::ostream& leavesOutputTabsWithHistHeight(const size_t bigN,
                        std::ostream &os, const int prec = 5) const;


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
        std::ostream& leavesOutputTabsWithHistHeightAndEMPs(const size_t bigN,
                            std::ostream &os, const int prec = 5) const;


        /*! \brief Get scaled  EMP sum under COPERR for tree rooted at this.

        \param n the total number of points in the histogram, for scaling
        */
        dotprecision getEMPSumCOPERR(const size_t n) const;

        /*! \brief Get the unscaled EMP sum under AIC for tree rooted at this.

        \param n the total number of points in the histogram, for scaling
        */
        dotprecision getEMPSumAIC(const size_t n) const;

        /*! \note virtual BOOL_INTERVAL spContains(const ivector& z) const
        is not redeclared so use base. */

        /*! \note virtual BOOL_INTERVAL spContains(const rvector& p) const
        is not redeclared so use base. */

               /*! \brief Check if the box a node represents contains a datapoint p.

        \param p the value of the data point being tested for containment in
        the box represented by this node.
        \param childInd indicates whether this should be considered
        to be a left child or a right child (ie where we need to take
        splitting dimension and value into account) or as a parent
        node (default).
		
		If a node is being considered as a child node (ie childInd is 
		ON_RIGHT or ON_LEFT) and actaully has a parent node, 
		then is is assumed that the data 
		point would have been in the box associated with the parent
		node, and the question now is just whether it falls into 
		the left or right child. Thus no test for full containment of
		the data point in the box is carried out, and this method
		just checks whether the data point would have fallen 
		from the parent into the left or right child.
		 
		If childInd is ON_PARENT, or if the node has no parent node,
		then the method checks for full containment of the data point
		within the box associated with this node.

        Note that the interval on the parent's split dimension
        of the right child's box is closed at the split value and
        the interval of the left child's box is open.  A datapoint
        whose element in the split dimension is exactly the split
        value should be assessed to be in the right child's box but
        not the left child's box.
        */
        bool nodeContains(const rvector& p,
                        OPERATIONS_ON childInd = ON_PARENT) const;
								
			/*! \brief Get the number of points in any box.

        \param z query box
        \param countBox
		  \param countInBox 

        Later.
        */
        int spsContains(ivector & z, int countBox, int countInBox) const;							
								

        /*! \brief Expand a leaf node.

        Expand a leaf node to have two children and pass
        data  down to the children with no further splitting.

        Uses nodeExpansion() and splitData().

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
        void nodeExpand(const SplitDecisionObj& boolTest, int comp);

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
        void nodeExpand(const SplitDecisionObj& boolTest);

        /*! \brief Reabsorbs both children of the node.

        Effectively reverses any split of the node.

        Data associated with the children is pushed back up to this
        and the splitDim and splitValue reset to leaf defaults.

        Works even if the children are not leaves.
        */
        virtual void nodeReabsorbChildren();


        /*! \brief Try to reunite nodes to form one leaf.

        Note that the nodes provided, lChild and rChild, are not
        the actual children of this, they are potential children
        which we are trying to either totally bring into
        this (if there are two of them) or to graft onto this
        if there is only one of them.  This is typically a new, part-formed
        node whose formation can be completed by reuniting already two
        already-formed nodes into it or by adding on one child if only
        one is available.  nodeReunite is used in building a tree upwards
        (rather than in pruning leaves of formed tree from the bottom up).

        If two potential children are provided and they are both
        leaves, it combines the two leaf siblings into this.  If the
        potential children are not leaves or if only one potential
        child is provided, it grafts the potential child/children
        onto this as its child/children.

        Data associated with the children is
        moved to the new parent and statistics recalculated.

        \warning nodeReunite would not normally be used with SPSnodes
        but is in the base class and is reimplemented to try do it
        appropriately for this derived class should it be needed.
        This function is untested.
        */
        virtual void nodeReunite(SPSnode *lChild, SPSnode* rChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a left child rather than attempting
        to reunite two children into this.
        \warning not thoroughly tested */
        virtual void nodeAdoptLeft(SPSnode *lChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a right child rather than attempting
        to reunite two children into this.
        \warning Not thoroughly tested */
        virtual void nodeAdoptRight(SPSnode *rChild);


        /*! \brief Inserts data into this node.

        Called recursively from the root and through the tree, seeking
        leaf node whose box contains the data point.  If data is
        inserted, this  method also tests whether the node should be
        expanded following the addition of the data.  Following an
        expansion, insertOneFind is used again to pass the the node's
        data down to its new children.

        \param newItr an iterator to the data in big data collection.
        \param childInd an indicator for whether the current node is a
        treated as a left or right child or a root.  This is passed to
        nodeContains() when testing whether the node contains the data.
        \param boolTest is a reference to an object providing a function
        operator determining whether to split a node when a data point arrives.
        This object can a dummy which never splits the node.
        \return a pointer to the node the data is 'inserted' into,
        before it is split, or NULL if no insert.
        */
        SPSnode* insertOneFind(BigDataItr newItr,
                            OPERATIONS_ON childInd,
                            const SplitDecisionObj& boolTest);


        /*! Makes the non-minimal union of two tree as a new tree with no data.

        Adds two pavings together as the union of the two but
        discards the actual data.  The tree manager should provide new data.

        Renames nodes from new root downwards.

        \param lhs pointer to root of first SPSnode tree operand.
        \param rhs pointer to root of second SPSnode tree operand.
        \return a pointer to root of a new SPSnode tree whose leaves are the
        union of the leaves of lhs, rhs and which has no data.
        */
        static SPSnode* unionTreeStructure(const SPSnode * const lhs,
                            const SPSnode * const rhs);

			//src_trunk_0701
			/*! \brief Swap this and another node.
			Swaps all the data members of this with the other node. 

			\param spn a reference to the node to swap with
			\post this is identical,in terms of its data members, 
			to spn before the swap, and spn is
			identical to this after the swap.*/
			void swapSPS(SPSnode& spn); // throw()


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

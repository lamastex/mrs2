/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
* Copyright (C) 2011 Gloria Teng
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

#ifndef ___SPSVnode_HPP__
#define ___SPSVnode_HPP__

// for spnode class
#include "spnode.hpp" // includes sptypes.hpp
#include "spsnode.hpp" 

/*! \file SPSVnode.hpp
\brief SPSVnode (StatsSubPavingVal) declarations

*/

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    //! Forward class declarations
    class SplitDecisionObj;

    /*! \brief StatsSubPavingVal is an alias for a pointer to an SPSVnode.
    */
    typedef SPSVnode* StatsSubPavingVal;


    /*! \brief A derived class based on SPSnode for processing sample data split
               to training and validation sets.

   The SPSVnode class is used to form a regular subpaving representing 
   containers of sample data using data-splitting methods. Here we implement the
   minimum distance estimator based on Devroye and Lugosi (200.)..clearer citation here.
    
    Leaves of the SPSVnode class have data associated with them in the form
    of pointers to some big collection of sample data for both the training and
    validation sets. If an SPSVnode is bisected the data associated with it descends to its children, so that only leaf SPSVnodes have data associated with them.  However, for the training set,  "recursively-computable statistical summaries", such as, count, sum, etc, of the data which would 
    be contained in the box an SPSVnode represents are kept for all SPSVnodes 
    and continue to be updated when the SPSVnode has children and data reaching the node is passed on to be finally associated with a leaf. For the validation set, only the count is maintained.
    
    By default, all recursively computable statistics provided except counts are not maintained in each SPSVnode, since this uses memory and is always 
    needed in the data splitting method implemented here.   
    
    Each node also has a boolean flag to indicate if this node has been split
    or not.
    */
    
    //---------------declaration of the SPSVnode class----------------------//
    class SPSVnode : public SPSnode {
    //-----------private methods
    private:
        /* theBox, dimension, label, parent, leftChild,
        rightChild, nodeName and recursively computable statistcs are inherited from base class */

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
        void recalculateStats(rvector& newdata, bool boolVal) const;

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
        void splitData(const SplitDecisionObj& boolTest, bool boolVal);

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
        static NodeData& gatherData(NodeData& container, SPSVnode * spn);


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

        \param lhs pointer to root of first SPSVnode tree operand.
        \param rhs pointer to root of second SPSVnode tree operand.
        \return a pointer to root of a new SPSVnode tree whose leaves are the
        union of the leaves of lhs, rhs and which has no data.
        */
        static SPSVnode* unionNoData(const SPSVnode * const lhs,
                            const SPSVnode * const rhs);


    protected:

        /*! \brief An indication of the maximum number of data points
        a node needs to carry.

        This is used for efficiency only to reserve vector space and a
        node can have more than this maximum number of data points
        associated with it.  Defaults to defaultMaxPts.
        */
        size_t spaceIndication;

        /*! \brief A boolean flag to know if this node was being split or node.
         
         False (not split) by default.
        */
        bool justSplit;
        
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
        //! A counter for how many data points from the validation set that are covered by theBox
        mutable size_t Vcounter;
         
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
        mutable VecDotPrec dpVSums;

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
        mutable VecDotPrec dpVSumProducts;

        /*! \brief A container for the association of validation data with a 
						 node.

        Data is associated with a node via this container of iterators.
        The iterators can, very loosely, in the sense in which they are
        used here, be thought of as pointers to a big data collection
        of all data points.  Only leaf nodes can have anything in this
        container.  However, not all leaf nodes will necessarily have
        something in this container:  the container will be empty if no
        data points are covered by the box represented by a leaf node.
        */
        mutable NodeData VdataItrs;
        //@}
        
	//-----------end of private methods
   
   //-----------public methods   
      public:
        /*! \brief Default constructor.
        */
        explicit SPSVnode();

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents
        and a value for countsOnly which controls whether all available stats
        are maintained (false) or just counts (true), and optionally
        initialised with a label which defaults to 0 if not provided.
        */
        explicit SPSVnode(ivector& v, bool cntOnly, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents
        and optionally initialised with a label which defaults to 0
        if not provided.  The value for countsOnly will default to false
        (i.e., all stats maintained by default).
        */
        explicit SPSVnode(ivector& v, int lab = 0);


        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents,
        a space indication, and a value for countsOnly which controls wether
        all available stats are maintained (false) or just counts (true), and
        optionally initialised with a label for the model which defaults
        to 0 if not provided.
        */
        explicit SPSVnode(ivector& v, size_t max, bool cntOnly, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with an interval vector for the box it represents,
        a space indication, and optionally initialised with a label for the
        model which defaults to 0 if not provided. The value for countsOnly
        will default to false (i.e., all stats maintained by default).
        */
        explicit SPSVnode(ivector& v, size_t max, int lab = 0);

        /*! \brief Initialised constructor.

        Initialised with a LabBox for the labeled box it represents.  Also
        optionally initialised with a value for countsOnly, defaults to
        false (i.e., all stats maintained by default).
        */
        explicit SPSVnode(LabBox& lb, bool cntOnly = true);

        /*! \brief Initialised constructor.

        Initialised with a LabBox for the labeled box it represents,
        and a space indication.  Also optionally initialised with a value for
        countsOnly, defaults to false (i.e., all stats maintained by default).
        */
        explicit SPSVnode(LabBox& lb, size_t max, bool cntOnly = true);

        /*! \brief Copy constructor.
        */
        explicit SPSVnode(const SPSVnode& other);

        // Use base class destructor

        /*! \brief Copy assignment operator.
        */
        SPSVnode& operator=(const SPSVnode& rhs);

        static SPSVnode* strippedConstructor(const SPSVnode * const other);
                     
        /*! \brief Accessor for the node's validation data counter.

        Returns Vcounter.
        */
        size_t getVcounter() const;

	    /*! \brief Accessor for the node's justSplit boolean flag.

        Returns justSplit.
        */
        bool getJustSplit() const;

       /*! \brief Accessor for the node's validation data collection.

        Returns a copy of the node's collection of iterators to the
        big data set.
        */
        NodeData getVdata() const;
         
         /*! \brief Clears the node's data collection for both training and
                    validation set.
         */          
         void clearData() const;

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
        SPSVnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Hides the base class version of this method.

        Returns a copy of the pointer to leftChild node.
        */
        SPSVnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

        Hides the base class version of this method.

        Returns a copy of the pointer to rightChild node.
        */
        SPSVnode* getRightChild() const;
        //@}

        /*! \brief The count the right child would have if this node was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the right child
        if the node were to be split.
        */
        size_t getRightCountIfSplit() const;

        /*! \brief The count the left child would have if this node was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the left child
        if the node were to be split.
        */
        size_t getLeftCountIfSplit() const;

        /*! \brief Smallest number of points in either child if  this was split.

        Does not split the nodes, just calculates how many of the data points
        currently associated with this node would go to the left and right child
        if the node were to be split.
        */
        size_t getMinChildCountIfSplit() const;

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
        SPSVnodePtrs& getLeaves(SPSVnodePtrs& leaves) const;

        /*! \brief Return a reference to all sub-leaf descendent nodes.

        Sub-leaf nodes have at least one child but any child must be a leaf,
        ie sub-leaves are the parents of leaf nodes.

        Will be just this if this is a subleaf.

        \return a reference to a container of node pointers.
        */
        SPSVnodePtrs& getSubLeaves(SPSVnodePtrs& subleaves) const;

		  /*! \brief Return a reference to all nodes.

        \return a reference to a container of node pointers.
        */
        SPSVnodePtrs& getAllNodes(SPSVnodePtrs& allNodes) const;


        /*! \brief The count in the node's ultimate ancestor root.
        */
        size_t getRootCounter() const;
        
        /*! \brief The count in the node's ultimate ancestor root.
        */
        size_t getRootVcounter() const;

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
			
			//gat41
			RealVec getUniformVarCovar() const;
			RealVec& getUniformVarCovar(RealVec& varCovar) const;

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

        label [tab] counter [tab] Vcounter [tab] volume [tabl] 
        Inf(ivector[1]) [tab] Sup(ivector[1]) [tab] ...
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
        void nodeExpand(const SplitDecisionObj& boolTest, bool boolVal);
        
        void nodeExpand(bool boolVal);

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

        \warning nodeReunite would not normally be used with SPSVnodes
        but is in the base class and is reimplemented to try do it
        appropriately for this derived class should it be needed.
        This function is untested.
        */
        virtual void nodeReunite(SPSVnode *lChild, SPSVnode* rChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a left child rather than attempting
        to reunite two children into this.
        \warning not thoroughly tested */
        virtual void nodeAdoptLeft(SPSVnode *lChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a right child rather than attempting
        to reunite two children into this.
        \warning Not thoroughly tested */
        virtual void nodeAdoptRight(SPSVnode *rChild);


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
        SPSVnode* insertOneFind(BigDataItr newItr,
                            OPERATIONS_ON childInd,
                            const SplitDecisionObj& boolTest, bool boolVal);
                            
                             /*! \brief Check if the box a node represents contains a datapoint p.

        \param p the value of the data point being tested for containment in
        the box represented by this node.
        \param childInd indicates whether this should be considered
        to be a left child or a right child (ie where we need to take
        splitting dimension and value into account) or as a parent
        node.

        childInd, together with the splitValue and splitDimension of
        the parent, is used to to make sure that the containment
        assessment takes notice of the open and closed intervals at
        the split value on the split dimension that result from
        splitting a box.  The interval in the split dimension
        of the right child's box is closed at the split value and
        the interval of the left child's box is open.  A datapoint
        whose element in the split dimension is exactly the split
        value should be assessed to be in the right child's box but
        not the left child's box.
        */
        bool nodeContains(const rvector& p,
                        OPERATIONS_ON childInd) const;


        /*! Makes the non-minimal union of two tree as a new tree with no data.

        Adds two pavings together as the union of the two but
        discards the actual data.  The tree manager should provide new data.

        Renames nodes from new root downwards.

        \param lhs pointer to root of first SPSVnode tree operand.
        \param rhs pointer to root of second SPSVnode tree operand.
        \return a pointer to root of a new SPSVnode tree whose leaves are the
        union of the leaves of lhs, rhs and which has no data.
        */
        static SPSVnode* unionTreeStructure(const SPSVnode * const lhs,
                            const SPSVnode * const rhs);

         //--------end of public methods
    };    // end of SPSVnode class derived from SPnode class


} // end namespace subpavings

#endif

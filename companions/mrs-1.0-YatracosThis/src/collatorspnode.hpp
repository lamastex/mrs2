/*
* Copyright (C) 2007, 2008, 2009, 2010, 2011 Raazesh Sainudiin
* Copyright (C) 2009, 2010 Jennifer Harlow
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

#ifndef ___COLLATORSPNODE_HPP__
#define ___COLLATORSPNODE_HPP__

/*! \file collatorspnode.hpp
\brief CollatorSPnode declarations

*/

#include "rangecollection_hist.hpp"

#include "spnode.hpp"

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    /*! \brief A derived class based on SPnode for creating summaries.

    Creates summaries from other nodes from the SP node family.

    The base class SPnode is a node in the representation of a regular
    subpaving as a binary tree.  A node represents a box (interval vector).
    SPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.

    The CollatorSPnode class collates data from a collection of SPnodes
    or objects from classes derived from SPnode.

    An entire tree structure represents a union over each tree collated, ie
    the boxes of the leaves of the collated tree make a subpaving that is the
    union of the subpavings represented by each tree collated.

    Each node has a container structure (the summary) holding one value for
    each corresponding collated.

    */
    class CollatorSPnode : public SPnode {
    private:
        /* theBox, dimension, label, parent, leftChild,
        rightChild and nodeName are inherited from base class.
        */

        /*! \brief A container of summary values from the collated subpavings.
        */
        mutable VecDbl summary;

        /*! Private initialised constructor.

        Initialised with a pointer to an SPSnode and a normalising constant,
        eg sum of counts in each node for a histogram
        The summary becomes count /(normalising constant * vol) for the SPSnode.

        \param bigN the normalising constant.
        */
        CollatorSPnode(const SPSnode * const spn, size_t bigN);


        /*! \brief Negates the summary for every node in tree rooted at this.
        */
        void nodeNegate(double c);


         /*! \brief Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.
        Outputs the average over the summary.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for an n-dimensional interval vector is

        nodeName [tab]  volume [tabl] average summary [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]
        */
        std::ostream& leafAverageOutputTabs(
                            std::ostream &os) const;


        /*! \brief Add to the accumulation of absolute areas for leaf nodes.

        If this node is a leaf, add its accumulation of summary areas to the
        overall accumulation. Area for a value s in the summary is s * volume
        of node.
        */
        VecDotPrec& getLeafNodeAbsAreaAccumulations(VecDotPrec& areaAcc) const;

			/*! \brief Get the total accumulated value of absolute areas
		  for leaf nodes.

        Area for a value s in a summary is s * volume of node.
        Areas for a leaf node are the areas for each subpaving collated by
        that node.  Accumulated areas over all the leaf nodes are, for each
        subpaving collated, the accumulated areas for that subpaving over all
        the leaf nodes.

        \return total of accumulated areas of leaf nodes as a real.
        */
        real getLeafNodeAbsAreaAccumulationTotal() const;
    

        /*! \brief Add to accumulation of summary values for leaf nodes.

        If this node is a leaf, add its summary values to the
        overall accumulation. If a value in the summary for a particular
        collated subpaving, s, s is added to the accumulation for that subpaving.
        */
        VecDotPrec& getLeafNodeSummaryAccumulations(VecDotPrec& summAcc) const;
        
        
        /*! Replace lhs data with the the dot difference of lhs and rhs data.

        If lhsSummary is originally <h1..hn> and rhsSummary is <H1..Hm>
        then the lhsSummary becomes <h1-H1, .. Hn-H1, .. , hn-H1, .. ,Hn-Hm>
        */
        static VecDbl& dotDifferenceSummary(VecDbl& lhsSummary,
                                VecDbl& rhsSummary);


        /*! \brief make this the root of a tree representing summary differences

        The tree will be the union of this tree and the spn tree
        The summary is the element-by-element difference between this's summary
        and the spn's summary.
        ie if this node has summary <h1..hn> and spn's equivalent node summary
        is <H1..Hm>
        then this's summary becomes <h1-H1, .. Hn-H1, .. , hn-H1, .. ,Hn-Hm>
        */
        bool dotDiffPaving(CollatorSPnode * const spn);
        
        
    public:

        /*! \brief default constructor,
        */
        explicit CollatorSPnode();

        /*! \brief Initialised constructor.

        Initialised with a pointer to an SPSnode.
        The collatorSPSnode summary information is the k/vol of
        the SPSnode, ie effectively the height of a histogram bin
        formed by the box associated with that node.
        */
        explicit CollatorSPnode(const SPSnode * const spn);

        /*! \brief Initialised constructor,

        Initialised with a box, a label, and a vector summary,
        */
        explicit CollatorSPnode(ivector& v, int lab, VecDbl summ);

        /*! \brief Copy constructor.
        */
        explicit CollatorSPnode(const CollatorSPnode& other);

        /*! \brief Copy assignment operator.
        */
        CollatorSPnode& operator=(const CollatorSPnode& rhs);


        /*! make a CollatorSPnode which represents an average.
        */
        CollatorSPnode* makeAverageCollation() const;


        // Use base class destructor.


        /*! @name Accessors for links between the nodes.
        These accessor methods shadow equivalent methods in the base
        class.  Thus the method used is determined  at compile time,
        not run time as would be the case if virtual methods were used.
        Because the pointers to parents and children are part of the
        base class definition, the methods have to cast the base class
        from to the derived class form in order for the pointer
        returned to be able to be used with derived class members.

        Note that pointers for parent, leftChild, and rightChild are
        not reference counted so there could potentially be problems
        with the use of returned pointers (for instance, being used to
        delete nodes). These pointers might be better implemented with
        boost::shared_ptr .
        */

        //@{
        /*! \brief Accessor for the parent of a node.

        Returns a copy of the pointer to parent node.
        */
        CollatorSPnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Returns a copy of the pointer to leftChild node.
        */
        CollatorSPnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

        Returns a copy of the pointer to rightChild node.
        */
        CollatorSPnode* getRightChild() const;
        //@}

        /*! \brief Accessor for the summary.
        */
        VecDbl getSummary() const;

        /*! \brief Get number of subpavings summarised.
        */
        size_t getNumberSummarised() const;


        /*! \brief Get the accumulations of absolute areas for leaf nodes.

        Area for a value s in a summary is s * volume of node.
        Areas for a leaf node are the areas for each subpaving collated by
        that node.  Accumulated areas over all the leaf nodes are, for each
        subpaving collated, the accumulated areas for that subpaving over all
        the leaf nodes.

        \return container of the accumulated leaf areas for all collated
        subpavings, in the same order as the subpavings are in the collator.
        */
        RealVec getLeafNodeAbsAreaAccumulations() const;

		/*! \brief Get the accumulations of summary values for leaf nodes.

        \return container of the accumulated summary values for all collated
        subpavings, in the same order as the subpavings are in the collator.
        */
        RealVec getLeafNodeSummaryAccumulations() const;
        
		
		/*! \brief Output the details of a specific node.

        This is intended for console output or output
        to a mixed alpha and numeric file.
        */
        virtual std::ostream& nodePrint(std::ostream &os) const;


        /*! \brief Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.
        Outputs summary for the node.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for an n-dimensional interval vector with m elements in the summary is

        label [tab] volume [tab] summary[1] [tab] ... summary[m] [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]
        */
        virtual std::ostream& leafOutputTabs(
                            std::ostream &os) const;



        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabs() to output summary information
        for each leaf node.
        */
        virtual std::ostream& leavesOutputTabs(
        std::ostream &os) const;



        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafAverageOutputTabs() to output information
        for each leaf node.  Outputs the average over the summary for each leaf.
        */
        std::ostream& leavesAverageOutputTabs(std::ostream &os) const;


        /*! \brief Expand leaf node to make two more leaves as children
        and copy summary down to the children.

        Equivalent to bisecting a box in a regular subpaving.  Makes
        two new sibling child nodes of this one and grafts them on.
        \param comp is the dimension on which to bisect theBox.
        */
        virtual void nodeExpand(int comp);

        /*! \brief Expand leaf node to make two more leaves as children
        and copy summary down to the children.

        Finds the dimension to split on and passes to
        NodeExpand(int comp) to carry out the bisection.
        */
        virtual void nodeExpand();

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
        leaves, it combines the two leaf siblings into this.
        If the potential children are not leaves or if only one
        potential child is provided, it grafts the potential
        child/children onto this as its child/children.

        Summary data associated with the children is related to
        the new parent.

        \warning nodeReunite would not normally be used with
        CollatorSPnodes but is in the base class and is
        reimplemented to try do it appropriately for this derived
        class should it be needed.  This function is untested.
        */
        virtual void nodeReunite(CollatorSPnode *lChild,
                                CollatorSPnode* rChild);

        /*! \brief Builds a higher level of a tree from existing nodes.
        This adopts a left child rather than attempting
        to reunite two children into this.
        */
        virtual void nodeAdoptLeft(CollatorSPnode *lChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a right child rather than attempting
        to reunite two children into this.
        */
        virtual void nodeAdoptRight(CollatorSPnode *rChild);


        /*! \brief Add two collator subpavings together.

        Collator returned has summary with contents of both operands' summaries.

        \param lhs pointer to the root of one collator to be added.
        \param lhs pointer to the root of the other collator to be added.
        \return a pointer to the root of a new collator subpaving where the
        summmary is the combined summary of the operands (lhs first, then rhs).
        */
        static CollatorSPnode* addPavings(
                    const CollatorSPnode * const lhs,
                    const CollatorSPnode * const rhs);


        /*! \brief Subtract one collator subpavings from another together.

        Collator returned has tree which is the union of lhs and rhs trees and
        summary with lhs summary and negative of rhs summary values.

        \param lhs pointer to the root of one collator.
        \param lhs pointer to the root of the other collator to be subtracted.
        \return a pointer to the root of a new collator subpaving where the
        summmary is the summary of lhs and the negative values from the summary
        of rhs.
        */
        static CollatorSPnode* subtractPavings(
                    const CollatorSPnode * const lhs,
                    const CollatorSPnode * const rhs, double c);


        /*! \brief Make a tree which holds differences of this to avg over this.

        The tree has the same structure as this tree and the summary values are
        the difference between the summary values for this
        and the average summary value over the summary for this
        ie if this has summary <h1..hn> and the mean is h. then
        then the new summary becomes <h1-h., .. Hn-h.>
        */
        CollatorSPnode* makeDifferencesToAveragePaving() const;

		
		/*! \brief Sum of variances of a scalar value.

        The variance for this scalar value is squared sum of 'area' of difference.

        Gives the sum of the variances over the collation where the variance
        of one of the collated subpaving trees is taken as the
        square of the sum of the 'areas' of difference between that collated
        subpaving tree and the average subpaving tree over the collation,
        where 'area' is absolute value of summary value * node volume.

        \return the sum of variances over all the collated subpaving trees.
        */
        real getSumVarsAreaScalar() const;
        
        
        /*! \brief Sum of variances of a scalar value.

        This scalar value is the total summary values.

        Gives the sum of the variances over the collation where the variance
        of one of the collated subpaving trees is taken as the
        square of the difference between the total summary values for that
        tree and the total summary values for the average tree across
        the collation.

        \return the sum of variances over all the collated subpaving trees.
        */
        real getSumVarsTotalSummarisedValueScalar() const;


		/*! \brief Incorporate a Collator subpaving to this summmary.

        Adjusts the tree rooted at this node for the tree being added.
        The tree rooted at this node will expand if necessary to have at least
        all the leaves of the tree of the tree being added.  The summary for
        this node will increase to include the summary of the node being added.

        \param spn pointer to the Collator subpaving to be incorporated into
        the summary.
        \return true if the paving could be added, false otherwise.
        \note that the contents of the Collator subpaving added will be
         altered: it is assumed that the Collator subpaving to be added is a
         temporary object made and passed in from a wrapper class
        */
        bool addPaving(CollatorSPnode * const spn);



        /*! \brief Incorporate negation of a Collator subpaving to this summmary.

        Adjusts the tree rooted at this node for the tree being added.
        The tree rooted at this node will expand if necessary to have at least
        all the leaves of the tree of the tree being added.  The summary for
        this node will increase to include the negation of the summary
        of the node being added.

        \param spn pointer to the Collator subpaving to be incorporated into
        the summary.
        \note that the contents of the Collator subpaving added will be
         altered: it is assumed that the Collator subpaving to be added is a
         temporary object made and passed in with a wrapper class
        */
        void addNegatedPaving(const CollatorSPnode * const spn, double c);
		  
		 //--removed then re-inserted
		 /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafAccumulationOutputTabs() to output information
        for each leaf node.  Outputs a accumulated summary for each leaf
        node, ie the sum of the summary.
        */
        std::ostream& leavesAccumulationOutputTabs(
                                            std::ostream &os) const;
														  
														  
      
        /*! \brief Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.
        Outputs the accumulated sum of the summary for the node,
        ie the summary added together.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for an n-dimensional interval vector with m elements in the summary is

        label [tab] volume [tab] empirical measure [tab]
		  sum(summary[1]... summary[m]) [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]
        */
        std::ostream& leafAccumulationOutputTabs(
                            std::ostream &os) const;
		
		    /*! Find the accumulated summary for a node.

        \return the sum over the summary for this node.
        */
        real nodeAccumulation() const;			
		  
		  /*! Change the current summary in the leaf node by adding mass
		  */
		 void leafMakeNewFhat(double wt, std::vector<double> & fhatNew);
		 
		 /*! Change the current summary in the leaf nodes by adding mass
		 */
		  void leavesMakeNewFhat(double wt, std::vector<double>  &fhatNew);	 
		
		 /*! return a reference to a container of node pointers.
       */ 
		std::vector<CollatorSPnode*>& getLeaves(std::vector<CollatorSPnode*>& leaves) const;
		
		 /*! return a reference to a container of node pointers.
       */ 
		std::vector<CollatorSPnode*>& getAllNodes(std::vector<CollatorSPnode*>& allNodes) const;
		
		
		/*! \brief Total the summaries from this node downwards.

        Works recursively to make each summary into a single element
        */
		void totaliseSummaries();
		
		/*! \brief Make a marginalised version of a given node.

        Private version, works recursively to marginalise from this node down.
        
        Marginalises to take out the given dimensions and adjust summaries
        so that overall sum of (node vol x accumulated summaries) 
        is the same as for \a rhs.
        
        \param rhs is the root of the node to marginalise.
        \outDims is a vector of the dimensions to take out.
        \return root of a tree of marginalised nodes.
        */
		static CollatorSPnode* _marginalise(
			const CollatorSPnode * const rhs,
			const std::vector<int>& outDims);
			
			        /*! \brief Make a marginalised version of subpaving with root
        node \a rhs.

        Marginalises to take out the given dimensions and adjust summaries
        so that overall sum of (node vol x accumulated summaries) 
        is the same as for \a rhs.
        
        \param rhs is the root of the node to marginalise.
        \reqDims is a vector of the dimensions to include in marginal.
        \return root of a tree of marginalised nodes.
        \pre \a rhs must be non-NULL.
        \pre \a reqDims must be compatible with current dimensions.
        \note allowed dimensions start at 1, ie dimensions to
        marginalise on can include 1, 2, ... #dimensions of \rhs
        \post returned tree will have one summary value for each node and 
        have sum of 
        (node vol x accumulated summaries) = that for \a rhs.
        */
		static CollatorSPnode* marginalise(
				const CollatorSPnode * const rhs,
				const std::vector<int>& reqDims);
		
		
		/*! \brief Find the dimension a node split on.

        split dimension is -1 if node has not split.
        
        \return dimension this node split on if it is not a leaf,
        -1 it if is a leaf.
        */
		int getSplitDim() const;
   
		/*! Get the total of the summary in this.
        */
		double getTotalSummary() const;
		
		/*! Get the total of the summary in this.
      */
		double getTotalSummaryAv() const;


		// Jenny addition for Gloria's convergence work
		// take a container and return the same container, which has been
		// cleared (if necessary) and re-filled with 
		// L1-distances-to-average, one for each histogram in collation
		RealVec& getL1DistancesToAverage(RealVec& container) const;
		
		//gat41
		 /*! \brief Find the node that fulfills the Yatracos condition by 
		     comparing the rows of the growing Yatracos matrix against the 
			  columns.
			  \param theta current split number
			  \param k split number in the set 0:(theta-1)
       */
       bool nodeCheckRowSummary(int theta, int k);

       /*! \brief Find the node that fulfills the Yatracos condition by 
		     comparing the columns of the growing Yatracos matrix against the 
			  rows.
		     \param theta current split number
			  \param k split number in the set 0:(theta-1)
       */
       bool nodeCheckColSummary(int theta, int k);

		
		//gat41
		  /*! \brief Get the delta value for a specific theta
      
				Get the delta value for a specific theta.
				\param k split number in the set 0:(theta-1) 
				\param theta the current split number
				\return the delta value
      */
      double getNodeDelta(int thisTheta, size_t sizeColl);
      
      //gat41
      /*! \brief Get the Yatracos set for a particular pair. */
   	void getYatSet(
	       std::set<CollatorSPnode*, std::less<CollatorSPnode*> > & YatSetRow,
          std::set<CollatorSPnode*, std::less<CollatorSPnode*> > & YatSetCol,
          size_t cand1, size_t cand2);
		
		//gat41
		 /*! \brief Find the node that fulfills the Scheffe set, 
		  * 			f_theta1 > f_theta2 for candidates f_theta1 and f_theta2, 
		  * 			where theta1 < theta2 
			  \param theta1 position of candidate f_theta1
			  \param theta2 position of candidate f_theta2
       */
       bool getScheffeNode(int theta1, int theta2);
       
		//gat41
		/*! \brief Get the Scheffe set for a particular pair. */
   	void getScheffeSet(
	       std::set<CollatorSPnode*, std::less<CollatorSPnode*> > & ScheffeSet,
          size_t cand1, size_t cand2);
          
       //src_trunk_0701
       	/*! \brief Swap this and another node.
		
			Swaps all the data members of this with the other node. 
					
			\param spn a reference to the node to swap with
			\post this is identical,in terms of its data members, 
			to spn before the swap, and spn is
			identical to this after the swap.*/
			void swapCollator(CollatorSPnode& spn); //throw()

};
    // end of CollatorSPnode class derived from SPnode class

/*! \brief Output operator for CollatorSPnodes.
*/
std::ostream & operator<<(std::ostream &os,
                        const CollatorSPnode& spn);


/*! \brief Comparison of CollatorSPnodes using total of summaries.
*/
bool nodeCompTotalSummary(const CollatorSPnode * const lhs,
                            const CollatorSPnode * const rhs);
	
/*! \brief Comparison of CollatorSPnodes using the average of the total 
 * 			of summaries.
*/
bool nodeCompTotalSummaryAv(const CollatorSPnode * const lhs,
                            const CollatorSPnode * const rhs);
	


} // end namespace subpavings

/*! A specialisation of std::swap for NewCollatorNode types.*/
/*namespace std
{
	template <>
	void swap (CollatorSPnode & s1, 
			CollatorSPnode & s2);
}*/

#endif


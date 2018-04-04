/*
* Copyright (C) 2007, 2008, 2009, 2010, 2011 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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

#ifndef ___CollatorSPVnode_HPP__
#define ___CollatorSPVnode_HPP__

/*! \file CollatorSPVnode.hpp
\brief CollatorSPVnode declarations

*/

#include "spnode.hpp"

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/

//=============================subpavings namespace=============================
namespace subpavings {

    /*! \brief A derived class based on SPnode for creating summaries.

    Creates summaries and empirical measures from other nodes from the SPnode
	 family.

    The base class SPnode is a node in the representation of a regular
    subpaving as a binary tree.  A node represents a box (interval vector).
    SPnodes are linked together to form the tree.  The initial box of
    the subpaving is the box represented by the root node of the tree.
    A box which has been split will be represented as node with one or
    two children.  A subpaving of [<b>x</b>] (union of non-overlapping sub-
    boxes of [<b>x</b>]) is represented by the leaves
    (degenerate/ child-less) nodes in the tree.

    The CollatorSPVnode class collates data from a collection of SPnodes
    or objects from classes derived from SPnode.

    An entire tree structure represents a union over each tree collated, ie
    the boxes of the leaves of the collated tree make a subpaving that is the
    union of the subpavings represented by each tree collated.

    Each node has a container structure (the summary) holding one value for
    each corresponding collated, and also the empirical measure of the data 
	 contained inside the node. 

    */
    class CollatorSPVnode : public SPnode {
    private:
        /* theBox, dimension, label, parent, leftChild,
        rightChild and nodeName are inherited from base class.
        */

        /*! \brief A container of summary values from the collated subpavings.
        */
        mutable VecDbl summary;
       
        /*! \brief Empirical measure for the validation data from
        the collated subpavings.
        */
        mutable double Vemp;

        /*! Private initialised constructor
            whatSum indicats what type of summary is to be tracked
            1: histogram density estimate
            2: counts only
        */
        CollatorSPVnode(const SPSVnode * const spn, size_t bigN, size_t bigM,
                        int whatSum);
        
        /*! \brief Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.
        Outputs the average over the summary.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here includes the box
        label and produces numeric tab-delimited data.  The format
        for an n-dimensional interval vector is

        nodeName [tab]  volume [tab] empirical measure [tab] 
		  average summary [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]
        */
        std::ostream& leafAverageOutputTabs(
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
        std::ostream& leafDifferenceOutputTabs(
                            std::ostream &os) const;

        /*! Find the accumulated summary for a node.

        \return the sum over the summary for this node.
        */
        double nodeAccumulation() const;
        double nodeDifference() const;
		  
		   /*! return a reference to a container of node pointers.
       */ 
		std::vector<CollatorSPVnode*>& getLeaves(std::vector<CollatorSPVnode*>& leaves) const;
		
		 /*! return a reference to a container of node pointers.
       */ 
		std::vector<CollatorSPVnode*>& getAllNodes(std::vector<CollatorSPVnode*>& allNodes) const;
		  
		  

        /*! Find the accumulation of absolute values in the summary for a node.

        \return the sum over the absolute values in the summary for this node.
        */
        double nodeAbsAccumulation() const;

        /*! Find the accumulated summary multiplied by volume for a node.

        \return the sum over the summary for this node x node volume.
        */
        double nodeAccumulationMultVol() const;

        /*! Find accumulation of absolute values in summmary x volume in node.

        \return the sum over the absolute values in the summary, multiplied by
        the volume, for this node.
        */
        double nodeAbsAccumulationMultVol() const;

    public:

        /*! \brief default constructor.
        */
        explicit CollatorSPVnode();

        /*! \brief Initalised constructor.
   
        Initialised with a pointer to an SPSVnode.
        The CollatorSPSVnode summary information is the k/vol of
        the SPSVnode, ie effectively the height of a histogram bin
        formed by the box associated with that node.
		  The empirical measure is k/m of the SPSVnode, where m is the total 
		  number of validation points. 
        */
        explicit CollatorSPVnode(const SPSVnode * const spn, int whatSum);

        /*! \brief Initialised constructor.
     
        Initialised with a box, a label, a vector summary, and an empirical 
		  measure.
        */
        explicit CollatorSPVnode(ivector& v, int lab, VecDbl summ, 
                                double Vsumm);
        

        /*! \brief Copy constructor.
        */
        explicit CollatorSPVnode(const CollatorSPVnode& other);
		  
        /*! \brief Copy constructor for the subtracted ADHVC.
        */
        explicit CollatorSPVnode(const CollatorSPVnode& other, int toSubtract);
        
        /*! \brief Copy assignment operator.
        */
        CollatorSPVnode& operator=(const CollatorSPVnode& rhs);

        /*! make a CollatorSPVnode which represents an average.
        */
        CollatorSPVnode* makeAverageCollation() const;

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
        CollatorSPVnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Returns a copy of the pointer to leftChild node.
        */
        CollatorSPVnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

        Returns a copy of the pointer to rightChild node.
        */
        CollatorSPVnode* getRightChild() const;
        //@}

        /*! \brief Accessor for the summary.
        */
        VecDbl getSummary() const;

        /*! \brief Accessor for the validation summary.
        */
        double getVemp() const;
        
        /*! \brief Get number of subpavings summarised.
        */
        size_t getNumberSummarised() const;

        /*! \brief sum over leaf nodes of absolute accumulated summary x volume.

        \return the sum over the descendent leaf nodes of this of the absolule
        accumulated summary x volume for each leaf.
        */
        double leavesAbsAccumulationMultVol() const;

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

        label [tab] volume [tab] empirical measure [tab] \
		  summary[1] [tab] ... summary[m] [tab]
        Inf(ivector[1]) [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]
        */
        virtual std::ostream& leafOutputTabs(std::ostream &os) const;
        virtual std::ostream& leafOutputTabs(std::ostream &os, 
																int whichColl) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafOutputTabs() to output summary information
        for each leaf node.
        */
        virtual std::ostream& leavesOutputTabs(
        std::ostream &os) const;
        virtual std::ostream& leavesOutputTabs(
        std::ostream &os, int whichColl) const;	

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafAverageOutputTabs() to output information
        for each leaf node.  Outputs the average over the summary for each leaf.
        */
        std::ostream& leavesAverageOutputTabs(std::ostream &os) const;

        /*! \brief Output for for <b>all leaves</b> of a binary tree

        Output intended for a txt file, in numeric form only.

        Recursively uses leafAccumulationOutputTabs() to output information
        for each leaf node.  Outputs a accumulated summary for each leaf
        node, ie the sum of the summary.
        */
        std::ostream& leavesAccumulationOutputTabs(
                                            std::ostream &os) const;
        std::ostream& leavesDifferenceOutputTabs(
                                            std::ostream &os) const;


        /*! \brief Expand leaf node to make two more leaves as children
        and copy summary and Vemp down to the children.

        Equivalent to bisecting a box in a regular subpaving.  Makes
        two new sibling child nodes of this one and grafts them on.
        \param comp is the dimension on which to bisect theBox.
        */
        virtual void nodeExpand(int comp);

        /*! \brief Expand leaf node to make two more leaves as children
        and copy summary and Vemp down to the children.

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
       // virtual void nodeReunite(CollatorSPVnode *lChild,
       //                         CollatorSPVnode* rChild);
			
		  virtual void nodesReunite();

        /*! \brief Builds a higher level of a tree from existing nodes.
        This adopts a left child rather than attempting
        to reunite two children into this.
        */
        virtual void nodeAdoptLeft(CollatorSPVnode *lChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a right child rather than attempting
        to reunite two children into this.
        */
        virtual void nodeAdoptRight(CollatorSPVnode *rChild);
		  
        /*! \brief Add paving with Vemp to collation
        */
        bool addPavingWithVal(CollatorSPVnode * const spn);
        
         /*! \brief Add two collator subpavings together.

        Collator returned has summary with contents of both operands' summaries.

        \param lhs pointer to the root of one collator to be added.
        \param lhs pointer to the root of the other collator to be added.
        \return a pointer to the root of a new collator subpaving where the
        summmary is the combined summary of the operands (lhs first, then rhs).
        */
        static CollatorSPVnode* addPavings(
                    const CollatorSPVnode * const lhs,
                    const CollatorSPVnode * const rhs);


        /*! \brief Subtract one collator subpavings from another together.

        Collator returned has tree which is the union of lhs and rhs trees and
        summary with lhs summary and negative of rhs summary values.

        \param lhs pointer to the root of one collator.
        \param lhs pointer to the root of the other collator to be subtracted.
        \return a pointer to the root of a new collator subpaving where the
        summmary is the summary of lhs and the negative values from the summary
        of rhs.
        */
        static CollatorSPVnode* subtractPavings(
                    const CollatorSPVnode * const lhs,
                    const CollatorSPVnode * const rhs, double c);

        /*! \brief Negates the summary for every node in tree rooted at this.
        */
        void nodeNegate(double c);
        
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
        void addNegatedPaving(const CollatorSPVnode * const spn, double c);
		  
   
		 /*! \brief  Get a CollatorSPVnode pointer to the corresponding SPSVnode that was split 
	   
			\param chosenLargest the current splitted node
			\return splitCollNode the current splitted node as a CollatorSPVnode	  	
      */
	   bool getSplitNodePtrCSPV(CollatorSPVnode * &splitCollNode, SPSVnode * const spn);
	 
       /*! \brief Find the node that fulfills the Scheffe set, 
		  * 			f_theta1 > f_theta2 for candidates f_theta1 and f_theta2, 
		  * 			where theta1 < theta2 
			  \param theta1 position of candidate f_theta1
			  \param theta2 position of candidate f_theta2
       */
       bool getScheffeNode(int theta1, int theta2);
		 	 
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

       /*! \brief Get the delta value for a specific theta
      
				Get the delta value for a specific theta.
				\param k split number in the set 0:(theta-1) 
				\param theta the current split number
				\return the delta value
      */
      double getNodeDelta(int thisTheta);
		
		/*! \brief Get the Yatracos set for a particular pair. */
   	void getYatSet(
	       std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > & YatSetRow,
          std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > & YatSetCol,
          size_t cand1, size_t cand2);
			 
		/*! \brief Get the Scheffe set for a particular pair. */
   	void getScheffeSet(
	       std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > & ScheffeSet,
          size_t cand1, size_t cand2);
     
}; 
//=========end of CollatorSPVnode class derived from SPnode class================

//================non-member function declarations==============================
/*! \brief Output operator for CollatorSPVnodes.
*/
std::ostream & operator<<(std::ostream &os,
                        const CollatorSPVnode& spn);

/*! \brief Negate a double.
*/
double opNegate(double d);

} // end namespace subpavings

#endif

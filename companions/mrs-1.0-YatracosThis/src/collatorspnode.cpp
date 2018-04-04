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

/*!/ \file:     collatorspnode.cpp
\brief CollatorSPnode definitions
*/

#include "collatorspnode.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <exception>

//to use algorithms
#include <algorithm>
//to use functionals
#include <functional>
//to use numeric accumulate
#include <numeric>

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"

#include "spsnode.hpp"

using namespace std;

namespace subpavings {


    // -------------------------- private member functions -------------

    // private initialised constructor, initialised with a pointer to an SPSnode
    // and a normalising constant, eg sum of counts in each node for a histogram
    // the summary becomes count /(normalising constant * vol) for the SPSnode
    CollatorSPnode::CollatorSPnode(const SPSnode * const spn, size_t bigN)
    {
        try {
            theBox = new ivector(spn->getBox());
            dimension = spn->getDimension();
            label = spn->getLabel();
            nodeName = spn->getNodeName();

            // add the summary to the vector summary
            summary.push_back((spn->getCounter())/
                                (1.0 * bigN * spn->nodeVolume()));

            //recursion on the children
            if (spn->hasLCwithBox()) {
                nodeAddLeft(new CollatorSPnode(spn->getLeftChild(), bigN));
            }
            else leftChild=NULL;

            if (spn->hasRCwithBox()) {
                nodeAddRight(new CollatorSPnode(spn->getRightChild(), bigN));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

    }


    // negates the summary for every node in tree rooted at this
    void CollatorSPnode::nodeNegate(double c) 
	 {
      // transform(summary.begin(), summary.end(), summary.begin(),
      //                                negate<double>());
      for (int i=0; i < summary.size(); i++) {	
			double temp = c * summary[i];
			if ( temp == -0) { temp = 0; }
			summary[i] = temp;		    
		}

        // recurse on children
        if (hasLCwithBox()) getLeftChild()->nodeNegate(c);
        if (hasRCwithBox()) getRightChild()->nodeNegate(c);
    }


    // Print the average of the summary of a single leaf node, using tab delimiters
    std::ostream& CollatorSPnode::leafAverageOutputTabs(
                                        std::ostream &os) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            double summ = 0;
            VecDblIt it;

            summ = accumulate(summary.begin(), summary.end(), summ);

            double vol = nodeVolume();
            double av =  summ/(1.0*summary.size());

            // output the nodeName, nodeVolume
            os << nodeName;
            os << "\t" << vol;
            // followed by the average
            os << "\t" << av;

            // followed by intervals making up box using Inf & Sup
            // ie unlike cxsc output, there is no [  ] around them
            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
        }
    }


    // make the dot difference of two summaries
    VecDbl& CollatorSPnode::dotDifferenceSummary(VecDbl& lhsSummary,
                                VecDbl& rhsSummary)
    {
    // replace lhs data with the the dot difference of lhs and rhs data
        VecDbl originalSummary = lhsSummary;
        lhsSummary.clear();
        VecDblIt it;
        VecDblIt sit;
        for (sit = rhsSummary.begin(); sit < rhsSummary.end(); sit++) {
            for (it = originalSummary.begin(); it < originalSummary.end();
                                                        it++) {
                lhsSummary.push_back(*it - *sit); // push back the diff
            }
        }
        return lhsSummary;
    }


    // accumulates absolute value of leaf node areas,
    // left to right
    VecDotPrec& CollatorSPnode::getLeafNodeAbsAreaAccumulations(
                                        VecDotPrec& areaAcc) const
    {
        if (getLeftChild()!=NULL) {
            areaAcc = getLeftChild()->getLeafNodeAbsAreaAccumulations(areaAcc);
        }
        if (getRightChild()!=NULL) {
            areaAcc = getRightChild()->getLeafNodeAbsAreaAccumulations(areaAcc);
        }
        if (getLeftChild()==NULL && getRightChild()==NULL) {

            double vol = nodeVolume();

            VecDblIt it = summary.begin();
            VecDotPrecIt dpit;

            for (dpit = areaAcc.begin(); dpit < areaAcc.end(); dpit++) {

                cxsc::accumulate((*dpit), abs(*it), vol);

                it++;
            }
        }

        return areaAcc;
    }

	 // double total value of accumulated absolute leaf node areas
    real CollatorSPnode::getLeafNodeAbsAreaAccumulationTotal() const
    {
        size_t n = summary.size();
        dotprecision emptyDP(0.0);
        VecDotPrec areaAcc(n, emptyDP); // n copies of empty dot prec
		dotprecision tot(0.0);
		
        areaAcc = getLeafNodeAbsAreaAccumulations(areaAcc);
        RealVec retvalues;

		VecDotPrecIt it;
        for (it = areaAcc.begin(); it < areaAcc.end(); it++) {
            tot += *it; 
        }

        return rnd(tot); // round to nearest;
    }


    // accumulates the summary values over leaf nodes
    // left to right
    VecDotPrec& CollatorSPnode::getLeafNodeSummaryAccumulations(
                                        VecDotPrec& summAcc) const
    {
        if (getLeftChild()!=NULL) {
            summAcc =
                getLeftChild()->getLeafNodeSummaryAccumulations(summAcc);
        }
        if (getRightChild()!=NULL) {
            summAcc =
                getRightChild()->getLeafNodeSummaryAccumulations(summAcc);
        }
        if (getLeftChild()==NULL && getRightChild()==NULL) {

            VecDblIt it = summary.begin();
            VecDotPrecIt dpit;

            for (dpit = summAcc.begin(); dpit < summAcc.end(); dpit++) {
                cxsc::accumulate((*dpit), (*it), 1.0);
                it++;
            }
        }

        return summAcc;
    }

    // turn this into the root of a CollatorSPnode tree which has tree which is
    // the union of this tree and the spn tree and summary which is the dot
    // difference between this's summary and the spn's summary.
    // ie if this node has summary <h1..hn> and spn's equivalent node summary
    // is <H1..Hm>
    // then this's summary becomes <h1-H1, .. Hn-H1, .. , hn-H1, .. ,Hn-Hm>
    // have not specifed const data for the CollatorSPnode pointer,
    // because if we do that we can't expand it
    // but note that the CollatorSPnode passed in CAN BE ALTERED
    // private because this should never be used in external code
    bool CollatorSPnode::dotDiffPaving(CollatorSPnode * const spn)
    {

        bool retValue = false;
        bool done = false;  // indicator for done adding

        try {

            if (spn == NULL) {
                done = true;

            }

            // if the boxes are not the same we can't do anything
            if (!done && (theBox != NULL) && (*theBox != spn->getBox())) {
                throw SPnodeException("Boxes do not match");
            }

            // if this has no box yet it has not incorporated anything
            // and so we just use spn to construct this
            if (!done && (theBox == NULL)) {

                ivector v = spn->getBox();

                theBox = new ivector(v);
                dimension = Ub(v) - Lb(v) + 1;
                label = spn->getLabel();

                summary = spn->summary;

                //recursion on the children
                if (spn->leftChild) {
                    nodeAddLeft(new CollatorSPnode(
                        *(spn->getLeftChild())));

                }
                else leftChild=NULL;

                if (spn->rightChild) {
                    nodeAddRight(new CollatorSPnode(
                        *(spn->getRightChild())));
                 }
                else rightChild=NULL;

                done = true;
                retValue = true;

            } // end if theBox==NULL

            // do the rest only if done is not true

            // if this is a leaf and the paving to be added is a leaf
            // this summary becomes the dot difference of the summaries
            if (!done && !done && isLeaf() && spn->isLeaf()) {

                // make this summary into the dot difference of this summary
                // and spn summary
                VecDbl spnSummary = spn->getSummary();
                summary = dotDifferenceSummary(summary, spnSummary);
            }

            // else not done and not both leaves,
            // if this is not a leaf or the paving to be added
            // is not a leaf, we may need to split
            // and we will need to recurse further
            else if (!done && (!isLeaf() || !(spn->isLeaf()))) {

                // if this is leaf and spn not we need to split this
                if (isLeaf()) { // so spn can't be a leaf

                    nodeExpand();
                }

                // if spn is leaf and this is not we need to split spn
                // THIS WILL CHANGE the CollatorSPnode pointed to by spn

                if (spn->isLeaf()) { // so this can't be a leaf

                    spn->nodeExpand();
                }

                // make this summary into the dot difference of this summary
                // and spn summary
                VecDbl spnSummary = spn->getSummary();
                summary = dotDifferenceSummary(summary, spnSummary);


                // if they are were neither leaves originally
                // we go straight on to recursing with the children
                // otherwise expansions above are followed by recursion

                // recurse with children
                retValue = getLeftChild()->dotDiffPaving(
                        spn->getLeftChild());
                retValue = getRightChild()->dotDiffPaving(
                        spn->getRightChild());

            } // end of dealing with case where at least one is not a leaf
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory in dotDiffPaving" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return retValue;
    }



    // ------------------------ public member functions ----------------

    //default constructor,
     CollatorSPnode::CollatorSPnode() {}

    // initialised constructor, initialised with a pointer to an SPSnode
    // the summary becomes the k/(N*vol) for the SPSnode
    CollatorSPnode::CollatorSPnode(const SPSnode * const spn)
    {
        try {
            theBox = new ivector(spn->getBox());
            dimension = spn->getDimension();
            label = spn->getLabel();
            nodeName = spn->getNodeName();

            // add the normalised count/volume to the vector summary

            size_t rootCounter = spn->getRootCounter();

            summary.push_back(spn->getCounter()/
                                (rootCounter * spn->nodeVolume()));

            //recursion on the children using constructor which normalises counts
            if (spn->hasLCwithBox()) {
                nodeAddLeft(new CollatorSPnode(spn->getLeftChild(),
                                                    rootCounter));
            }
            else leftChild=NULL;

            if (spn->hasRCwithBox()) {
                nodeAddRight(new CollatorSPnode(spn->getRightChild(),
                                                    rootCounter));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

    }

    // initialised constructor
    // initialised with a box, a label, and a vector summary
    CollatorSPnode::CollatorSPnode(ivector& v, int lab, VecDbl summ)
        : SPnode(v, lab)
    {
        try {
            // copy the vector summary
            summary = summ;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }


    // Copy constructor
    // copies from given node downwards
    CollatorSPnode::CollatorSPnode(const CollatorSPnode& other)
        : SPnode(*(other.theBox), other.label)
    {
        try {
            summary = other.summary;
            nodeName = other.nodeName;

            //recursion on the children
            if (other.leftChild) {
                nodeAddLeft(new CollatorSPnode(
                    *(other.getLeftChild())));
            }
            else leftChild=NULL;

            if (other.rightChild) {
                nodeAddRight(new CollatorSPnode(
                    *(other.getRightChild())));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

    }

    //Copy assignment operator
    // copies from given node downwards
    CollatorSPnode& CollatorSPnode::operator=(const CollatorSPnode& rhs)
    {
        try {

            // delete the current children (deletes their children as well)
            if (leftChild != NULL) {
                delete getLeftChild();
                leftChild = NULL;
            }

            if (rightChild != NULL) {
                delete getRightChild();
                rightChild = NULL;
            }
            // and delete the current box
            if (theBox != NULL) {
                delete theBox;
                theBox = NULL;
            }

            parent=NULL;

            theBox=new ivector(*rhs.theBox);
            dimension = rhs.dimension;
            label = rhs.label;
            nodeName = rhs.nodeName;

            summary = rhs.summary;

            //recursion on the children
            if (rhs.leftChild) {
                nodeAddLeft(new CollatorSPnode(
                    *(rhs.getLeftChild())));
            }
            else leftChild=NULL;

            if (rhs.rightChild) {
                nodeAddRight(new CollatorSPnode(
                    *(rhs.getRightChild())));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return *this;

    }

    // make a CollatorSPnode which represents an average of the summary
    CollatorSPnode* CollatorSPnode::makeAverageCollation() const
    {
        CollatorSPnode* newnode = NULL;

        try {
            VecDbl newsumm;
            double summ = 0;
            VecDblIt it;
            // should change this to use for_each
            summ = accumulate(summary.begin(), summary.end(), summ);

            newsumm.push_back(
                summ/(1.0*(summary).size()));

            ivector v = getBox();

            // make the new node
            newnode = new CollatorSPnode(v, label, newsumm);
            newnode->setNodeName(nodeName); // set name to this node name

            //recursion on the children
            if (hasLCwithBox()) {
                newnode->nodeAddLeft(getLeftChild()->makeAverageCollation());
            }

            if (hasRCwithBox()) {
                newnode->nodeAddRight(getRightChild()->makeAverageCollation());
            }
        }

        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return newnode;

    }


    // Accessor for the parent of a node.
    //Returns a copy of the pointer to parent node.
    CollatorSPnode* CollatorSPnode::getParent() const
    { return (CollatorSPnode*) parent; }

    // Accessor for the left child of a node.
    // Returns a copy of the pointer to leftChild node, cast to this type
    CollatorSPnode* CollatorSPnode::getLeftChild() const
    { return (CollatorSPnode*) leftChild; }

    // Accessor for the right child of a node.
    // Returns a copy of the pointer to rightChild node, cast this type
    CollatorSPnode* CollatorSPnode::getRightChild() const
    { return (CollatorSPnode*) rightChild; }

    // Accessor for the summary.
    VecDbl CollatorSPnode::getSummary() const
    { return summary; }

    // Get number of subpavings summarised.
    size_t CollatorSPnode::getNumberSummarised() const
    { return summary.size(); }


    // real value of accumulated absolute leaf node areas
    RealVec CollatorSPnode::getLeafNodeAbsAreaAccumulations() const
    {
        size_t n = summary.size();
        dotprecision emptyDP(0.0);
        VecDotPrec areaAcc(n, emptyDP); // n copies of empty dot prec

        areaAcc = getLeafNodeAbsAreaAccumulations(areaAcc);
        RealVec retvalues;

        VecDotPrecIt it;
        for (it = areaAcc.begin(); it < areaAcc.end(); it++) {
            retvalues.push_back(rnd(*it)); // round to nearest
        }

        return retvalues;
    }

    // real value of accumulated leaf node summary
    RealVec CollatorSPnode::getLeafNodeSummaryAccumulations() const
    {
        size_t n = summary.size();
        dotprecision emptyDP(0.0);
        VecDotPrec areaAcc(n, emptyDP); // n copies of empty dot prec

        areaAcc = getLeafNodeSummaryAccumulations(areaAcc);
        RealVec retvalues;

        VecDotPrecIt it;
        for (it = areaAcc.begin(); it < areaAcc.end(); it++) {
            retvalues.push_back(rnd(*it)); // round to nearest
        }

        return retvalues;
    }

    // Print the details of a single leaf node, using tab delimiters
    std::ostream& CollatorSPnode::leafOutputTabs(
                    std::ostream &os) const
    {

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            // output the nodeName, nodeVolume
            os << nodeName;
            double vol = nodeVolume();
            os << "\t" << vol;
            // followed by the summary
            VecDblIt it;
            for (it = summary.begin(); it< summary.end(); it++) {
                os << "\t" << (*it);
            }
            // followed by intervals making up box using Inf & Sup
            // ie unlike cxsc output, there is no [  ] around them
            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }

        }

    }

    //Output for all the  leaf boxes in this, using tab delimiters
    std::ostream& CollatorSPnode::leavesOutputTabs(
                            std::ostream &os) const
    {
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafOutputTabs(os);
            return (os << "\n");

        }

            //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesOutputTabs(os);
        }

        if (getRightChild()!=NULL) {
            getRightChild()->leavesOutputTabs(os);
        }

    }


    //Output for all the  leaf boxes in this, using tab delimiters
    std::ostream& CollatorSPnode::leavesAverageOutputTabs(
                            std::ostream &os) const
    {
        // uses  member function leafAverageOutputTabs for node output
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafAverageOutputTabs(os);
            return (os << "\n");
        }
            //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesAverageOutputTabs(os);
        }
        if (getRightChild()!=NULL) {
            getRightChild()->leavesAverageOutputTabs(os);
        }
    }

    // Print the details of a of a specific node in a subpaving
    std::ostream& CollatorSPnode::nodePrint(std::ostream &os) const
    {
        // output for box in form:
        // box, volume, summary data

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            os << "Box is :";

            for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
                // c-xsc default output for intervals
                os << "  " << thisBox[i];   }

            os << std::endl;
            os << "Box volume is " << nodeVolume() << std::endl;
            os << "Summary data: " ;

            VecDblIt it;
            for (it = summary.begin(); it< summary.end(); it++) {
                os << *it << " ";
            }

            os << std::endl;
        }
        return os;

    }


    // add two sibling nodes to this provided this is a leaf
    // comp argument is passed to Upper() and Lower()
    // split the box in half normal to dimension set by comp
    void CollatorSPnode::nodeExpand(int comp)
    {
        try
        {
            // only do something if this CollatorSPnode is a leaf
            if (isLeaf()) {
                // ivectors to become boxes for new children
                ivector lC, rC;
                // Call Lower() and Upper() to put split boxes
                // into lC and rC respectively
                Lower(getBox(), lC, comp);
                Upper(getBox(), rC, comp);

                // make and add the new children
                nodeAddLeft(new CollatorSPnode(
                    lC, label, summary));

                nodeAddRight(new CollatorSPnode(
                    rC, label, summary));

                //name the new children
                getLeftChild()->setNodeName(nodeName + "L");
                getRightChild()->setNodeName(nodeName + "R");

                // new children have summary from this

            }
        }

        catch (bad_alloc&)
        {
            std::cerr << "Error allocating memory in "
                << "CollatorSPnode::nodeExpand()"
                << std::endl;
            throw;
        }


    }

    // add two sibling nodes to this provided this is a leaf
    // finds its own comp argument then calls nodeExpand(int comp)
    void CollatorSPnode::nodeExpand()
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(maxdiamcomp); // complete nodeExpand

    }

    // nodeReabsorbChildren() can use the base class implementation
    // (the summary for this will be correct so just delete the children ...)


    // computes a minimal subpaving from two sibling subpavings
    // a subpaving is minimal if it has no sibling leaves
    // a minimal subpaving is created by discarding sibling leaves
    // and create summary data for new parent from children
    // warning: nodeReunite would not normally be used with
    // CollatorSPnodes but is in the base class and is
    // reimplemented to try do it appropriately for this
    // derived class should it be needed.
    // This function is untested.
    void CollatorSPnode::nodeReunite(CollatorSPnode *lChild,
                                    CollatorSPnode *rChild)
    // lChild and rChild are the two subpavings to be reunited
    {
        // *this is the node which will become the parent

        // check that the labels match and exit if not
        if ((lChild->label != label ) || (rChild->label != label)) {
            throw SPnodeException("Labels do not match");
        }

        // if both subpavings are leaves and hull of boxes is x
        // discard them: *this is a leaf
        if (lChild->isLeaf() && rChild->isLeaf()) {
            if (*theBox != (*(lChild->theBox) |
                            *(rChild->theBox))) {
                throw SPnodeException("Boxes cannot be combined");

            }
            // how many elements in the left child's summary
            size_t n = (lChild->summary).size();
            // elements in summary for each child should be same
            if ((rChild->summary).size() != n) {
                throw SPnodeException("Summaries cannot be combined");

            }

            size_t i = 0;

            // put into this summary the average of summary of children
            for (i=0; i < n; i++) {
                summary.push_back(((lChild->summary)[i]
                                +(rChild->summary)[i])/2.0);
            }

            //discard the two subpavings given
            delete lChild;
            delete rChild;

        }

        else {  // at least one child is not a leaf
            // this has to adopt them rather than reuniting them
            nodeAdoptLeft(lChild);
            nodeAdoptRight(rChild);
            recursiveRename();
        }
    }

    // graft lChild onto this node
    // lChild could be a leaf or a non-leaf
    // takes care of the summary associated with lChild/its descendents
    // used when we are building a collator statistical subpaving upwards
    void CollatorSPnode::nodeAdoptLeft(CollatorSPnode *lChild)
    {
        // *this is the node which will become the parent

        size_t i = 0;

        size_t n = (lChild->summary).size();

        if (summary.empty()) { // no summary in this box already
            // put into this summary 0.5* the summary of the new child
            for (i=0; i < n; i++) {
                summary.push_back(0.5*(lChild->summary)[i]);
            }
        }
        else { // has summary already

            // we have to make summary for this match
            // that of the children
            // number of elements in this summary should = child
            if (summary.size() != n) {
                throw SPnodeException("Summaries do not match");
            }

            // store current summary temporarily
            VecDbl temp_summary = summary;

            summary.clear();

            // put into this summary the average of
            // the summary of the new child and the old summary
            for (i=0; i < n; i++) {
                summary.push_back(((lChild->summary)[i]
                                + temp_summary[i])/2.0);
            }
        }

        // point parent and child pointers in the right directions
        // nodeAddLeft() checks labels, hull size, present children
        nodeAddLeft(lChild);
    }

    // graft rChild onto this node
    // rChild could be a leaf or a non-leaf
    // takes care of the data associated with rChild/its descendents
    // used when we are building a collator statistical subpaving upwards
    void CollatorSPnode::nodeAdoptRight(CollatorSPnode *rChild)
    {
        // *this is the node which will become the parent
        size_t i = 0;

        size_t n = (rChild->summary).size();

        if (summary.empty()) { // no summary in this box already
            // put into this summary 0.5 * the summary of the new child
            for (i=0; i < n; i++) {
                summary.push_back(0.5*(rChild->summary)[i]);
            }
        }
        else { // has summary already

            // we have to make summary for this match
            // that of the children
            // number of elements in this summary should = child
            if (summary.size() != n) {
                throw SPnodeException("Summaries do not match");
            }
            // store current summary temporarily
            VecDbl temp_summary = summary;

            summary.clear();

            // put into this summary the average of the summary
            // of the new child and the old summary
            for (i=0; i < n; i++) {
                summary.push_back(((rChild->summary)[i]
                                + temp_summary[i])/2.0);
            }
        }

        // point parent and child pointers in the right directions
        // nodeAddRight() checks labels, hull size, present children
        nodeAddRight(rChild);
    }


    CollatorSPnode* CollatorSPnode::addPavings(
                    const CollatorSPnode * const lhs,
                    const CollatorSPnode * const rhs)
    {
        CollatorSPnode* newCollator = NULL;
        CollatorSPnode* temp = NULL;

        bool done = false;

        try {

            if (lhs == NULL && rhs == NULL) done = true; // return null

            // if exactly one is null we can return a copy of the non-null one
            if (!done && lhs == NULL && rhs != NULL) {

                newCollator = new CollatorSPnode(*rhs);
                done = true;

            }
            if (!done && lhs != NULL && rhs == NULL) {

                newCollator = new CollatorSPnode(*lhs);
                done = true;
            }
            // both not null
            if (!done && lhs != NULL && rhs != NULL) {

                if ((lhs->getBox() != NULL) &&
                                (lhs->getBox() != rhs->getBox())) {
                    throw SPnodeException("Boxes do not match");
                }

                newCollator = new CollatorSPnode(*lhs);
                temp = new CollatorSPnode(*rhs);

                newCollator->addPaving(temp);
                delete temp;
                temp = NULL;
            }
        }
        catch (bad_alloc& ba) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            if (temp != NULL) {
                delete temp;
                temp = NULL;
            }
            std::cerr << "Error allocating memory in addPavings" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
        catch (exception& e) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            if (temp != NULL) {
                delete temp;
                temp = NULL;
            }
            throw;
        }

        return newCollator;
    }


    // subtract pavings
    CollatorSPnode* CollatorSPnode::subtractPavings(
                    const CollatorSPnode * const lhs,
                    const CollatorSPnode * const rhs, double c)
    {
        CollatorSPnode* newCollator = NULL;

        bool done = false;

        try {

            if (lhs == NULL && rhs == NULL) done = true; // return null

            if (!done && lhs == NULL && rhs != NULL) {

                newCollator = new CollatorSPnode(*rhs);
                newCollator->nodeNegate(c);

                done = true;

            }
            if (!done && lhs != NULL && rhs == NULL) {

                newCollator = new CollatorSPnode(*lhs);
                done = true;
            }
            // both not null
            if (!done && lhs != NULL && rhs != NULL) {

                if ((lhs->getBox() != NULL) &&
                                (lhs->getBox() != rhs->getBox())) {
                    throw SPnodeException("Boxes do not match");
                }

                newCollator = new CollatorSPnode(*lhs);
                newCollator->addNegatedPaving(rhs, c);  // uses a temp copy of rhs
            }
        }
        catch (bad_alloc& ba) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            std::cerr << "Error allocating memory in addPavings" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
        catch (exception& e) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            throw;
        }

        return newCollator;
    }


    // make a tree which has the same structure as this tree but where the
    // summary values are the difference between the summary values for this
    // and the average summary value over the summary for this
    // ie if this has summary <h1..hn> and the mean is h. then
    // then the new summary becomes <h1-h., .. Hn-h.>
    CollatorSPnode* CollatorSPnode::makeDifferencesToAveragePaving() const
    {
        CollatorSPnode* newCollator = NULL;
        CollatorSPnode* avg = NULL;

        try {

            avg = makeAverageCollation(); // must delete at end!

            newCollator = new CollatorSPnode(*this);

            newCollator->dotDiffPaving(avg);

            delete avg;
            avg = NULL;
        }
        catch (bad_alloc& ba) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            if (avg != NULL) {
                delete avg;
                avg = NULL;
            }
            std::cerr
                << "Error allocating memory in makeDifferencesToAveragePavings"
                << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
        catch (exception& e) {
            if (newCollator != NULL) {
                delete newCollator;
                newCollator = NULL;
            }
            if (avg != NULL) {
                delete avg;
                avg = NULL;
            }
            throw;
        }

        return newCollator;
    }

    // get the sum of the variances of a scalar summary.
    // the variance of the scalar summary of one of the collated subpaving trees 
    // is the square of the sum of the 'areas' of difference between that collated
    // subpaving tree and the average subpaving tree over the collation,
    // where 'area' is summary value * node volume
    real CollatorSPnode::getSumVarsAreaScalar() const
    {
        CollatorSPnode* differences = NULL;
        dotprecision accSquares(0.0);

        // take this collation
        try {
            // make a tree of differences between this and average over this
            differences = makeDifferencesToAveragePaving();

            size_t n = summary.size(); //also the number of collations? 
            dotprecision emptyDP(0.0);
            VecDotPrec areaAcc(n, emptyDP); // n copies of empty dot prec

            areaAcc = differences->getLeafNodeAbsAreaAccumulations(areaAcc);

            VecDotPrecIt it;
            for (it = areaAcc.begin(); it < areaAcc.end(); it++) {
					// add rnd(*it) times rnd(*it) to accSquares
                cxsc::accumulate(accSquares, rnd(*it), rnd(*it)); // round to nearest
            }

            delete differences;
            differences = NULL;
        }
        catch (bad_alloc& ba) {
            if (NULL != differences) {
                delete differences;
                differences = NULL;
            }
            string oldmsg(ba.what());
            string msg = "Error allocating memory summing variances.";
            msg += " Orginal error: " + oldmsg;
            throw SPnodeException(msg);
        }
        catch (exception& e) {
            if (NULL != differences) {
                delete differences;
                differences = NULL;
            }
            throw;
        }

        // return sum of squared sum of areas of difference
        return rnd(accSquares); // CXSC round nearest

    }

    // get the sum of the variances of a scalar summary where the scalar
    // summary for one of the collated histograms is the total value of the
    // summary values for that histogram over all the leaf nodes.
    real CollatorSPnode::getSumVarsTotalSummarisedValueScalar() const
    {
        CollatorSPnode* average = NULL;
        dotprecision accSqrdDiffSumSummary(0.0);

        // take this collation
        try {
            // make a tree of the averages
            average = makeAverageCollation();

            size_t n = summary.size();
            dotprecision emptyDP(0.0);
            VecDotPrec summAcc(n, emptyDP); // n copies of empty dot prec

            // the accumulated summaries for this collation
            summAcc = getLeafNodeSummaryAccumulations(summAcc);

            VecDotPrec avAcc(1, emptyDP); // only 1 value in av summaries
            avAcc = average->getLeafNodeSummaryAccumulations(avAcc);

            VecDotPrecIt it;
            for (it = summAcc.begin(); it < summAcc.end(); it++) {

                real diffsumsummary = rnd((*it)) - rnd(avAcc[0]); // rnd nearest
                cxsc::accumulate(accSqrdDiffSumSummary, diffsumsummary,
                                                diffsumsummary); // acc squares
            }

            delete average;
            average = NULL;
        }
        catch (bad_alloc& ba) {
            if (NULL != average) {
                delete average;
                average = NULL;
            }
            string oldmsg(ba.what());
            string msg = "Error allocating memory summing variances.";
            msg += " Orginal error: " + oldmsg;
            throw SPnodeException(msg);
        }
        catch (exception& e) {
            if (NULL != average) {
                delete average;
                average = NULL;
            }
            throw;
        }

        // return sum of squared sum of areas of difference
        return rnd(accSqrdDiffSumSummary); // CXSC round nearest

    }


    // incorporate a subpaving to this summmary,
    // adjusts this summary for the contents of the subpaving added
    // have not specifed const data for the CollatorSPnode pointer,
    // because if we do that we can't expand it
    // but note that the CollatorSPnode passed in CAN BE ALTERED
    bool CollatorSPnode::addPaving(CollatorSPnode * const spn)
    {

        bool retValue = false;
        bool done = false;  // indicator for done adding

        try {

            if (spn == NULL) {
                done = true;

            }

            // if the boxes are not the same we can't do anything
            if (!done && (theBox != NULL) && (*theBox != spn->getBox())) {
                throw SPnodeException("Boxes do not match");
            }

            // if this has no box yet it has not incorporated anything
            // and so we just use spn to construct this
            if (!done && (theBox == NULL)) {

                ivector v = spn->getBox();

                theBox = new ivector(v);
                dimension = Ub(v) - Lb(v) + 1;
                label = spn->getLabel();

                summary = spn->summary;

                //recursion on the children
                if (spn->leftChild) {
                    nodeAddLeft(new CollatorSPnode(
                        *(spn->getLeftChild())));

                }
                else leftChild=NULL;

                if (spn->rightChild) {
                    nodeAddRight(new CollatorSPnode(
                        *(spn->getRightChild())));
                 }
                else rightChild=NULL;

                done = true;
                retValue = true;

            } // end if theBox==NULL

            // do the rest only if done is not true

            // if this is a leaf and the paving to be added is a leaf,
            // this just sucks in the counter from spn
            if (!done && !done && isLeaf() && spn->isLeaf()) {

                VecDbl temp = spn->getSummary();
                summary.insert(summary.end(), temp.begin(),temp.end());
                done = true;
                retValue = true;

            }

            // else not done and not both leaves,
            // if this is not a leaf or the paving to be added
            // is not a leaf, we may need to split
            // and we will need to recurse further
            else if (!done && (!isLeaf() || !(spn->isLeaf()))) {

                // if this is leaf and spn not we need to split this
                if (isLeaf()) { // so spn can't be a leaf

                    nodeExpand();
                }

                // if spn is leaf and this is not we need to split spn
                // THIS WILL CHANGE the CollatorSPnode pointed to by spn

                if (spn->isLeaf()) { // so this can't be a leaf

                    spn->nodeExpand();

                }

                // put in the data
                VecDbl temp = spn->getSummary();
                    summary.insert(summary.end(), temp.begin(),temp.end());

                // if they are were neither leaves originally
                // we go straight on to recursing with the children
                // otherwise expansions above are followed by recursion

                // recurse with children
                retValue = getLeftChild()->addPaving(spn->getLeftChild());
                retValue = getRightChild()->addPaving(spn->getRightChild());

            } // end of dealing with case where at least one is not a leaf
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory in addPaving" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return retValue;
    }

    // incorporate the negative of a subpaving to this summmary,
    // adjusts this summary for the contents of the subpaving added
    // have not specifed const data for the CollatorSPnode pointer,
    // because if we do that we can't expand it
    // but note that the CollatorSPnode passed in CAN BE ALTERED
    void CollatorSPnode::addNegatedPaving(const CollatorSPnode * const spn, double c)
    {
        SPSnode* temp = NULL;

        try {
            CollatorSPnode* temp = new CollatorSPnode (*spn);

            // negate the node passed in
            temp->nodeNegate(c);

            addPaving(temp);

            delete temp;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory in addNegatedPaving" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

   //---removed then reinserted
	 // Print the details of a single leaf node, using tab delimiters
    // the sum of the summary is printed out
    std::ostream& CollatorSPnode::leafAccumulationOutputTabs(
                    std::ostream &os) const
    {

        if(theBox != NULL) { // do nothing if there is no box
            ivector thisBox = *theBox; // copy theBox
            // output the nodeName, nodeVolume
       		os << nodeName;
       		double vol = nodeVolume();
            os << "\t" << vol;
	         // followed by the sum of the summary        
				os << "\t" << nodeAccumulation();

            // followed by intervals making up box using Inf & Sup
            // ie unlike cxsc output, there is no [  ] around them
            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {
                os << "\t" << Inf(thisBox[i])
                   << "\t" << Sup(thisBox[i]);
            }
        }
    }

    // find the accumulated summary for a node
    real CollatorSPnode::nodeAccumulation() const
    {
     	dotprecision dpSumm;    // use type dotprecision for summation  
      dpSumm=0.0;
		VecDblIt it;
      // should change this to use for_each
            for (it = summary.begin(); it< summary.end(); it++) {
					 accumulate(dpSumm, (*it), 1.0);
            }
				real summ = rnd(dpSumm);
            return summ;
    }
	 	  
     //Output for all the  leaf boxes in this, using tab delimiters
     std::ostream& CollatorSPnode::leavesAccumulationOutputTabs(
                             std::ostream &os) const
     {
         // uses  member function leafAccumulationOutputTabs for nodes
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
             leafAccumulationOutputTabs(os);
             return (os << "\n");
 
         }
 
             //recurse on the children
         if (getLeftChild()!=NULL) {
             getLeftChild()->leavesAccumulationOutputTabs(os);
         }
 
         if (getRightChild()!=NULL) {
             getRightChild()->leavesAccumulationOutputTabs(os);
         }
 
     }
 
//---gloria's additions------------------------//
// Change the summary of the nodes
void CollatorSPnode::leafMakeNewFhat(double wt, std::vector<double> & fhatNew) 
{
   if(theBox != NULL) { // do nothing if there is no box
      double summ = 0;
      VecDblIt it;
      summ = accumulate(summary.begin(), summary.end(), summ);
      double av =  summ/(1.0*summary.size());          
		fhatNew.push_back((1-wt)*av + wt);
	}	
}
   
//Change the summary of the nodes
void CollatorSPnode::leavesMakeNewFhat(double wt, std::vector<double> & fhatNew)
{
   // uses  member function leafMakeNewFhat for nodes
   if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
      leafMakeNewFhat(wt, fhatNew);				
   }
 
   //recurse on the children
   if (getLeftChild()!=NULL) {
      getLeftChild()->leavesMakeNewFhat(wt, fhatNew);
   }
 
  if (getRightChild()!=NULL) {
      getRightChild()->leavesMakeNewFhat(wt, fhatNew);
   }
}
 
    // return a reference to a container of CollatorSnodes
    // contents being the leaves descended from this, or this if this is a leaf
    // left to right order
    std::vector<CollatorSPnode*>& CollatorSPnode::getLeaves(std::vector<CollatorSPnode*>& leaves) const
    {
        //if children, recurse on the children
        if (hasLCwithBox()) {
            getLeftChild()->getLeaves(leaves);
        }

        if (hasRCwithBox()) {
            getRightChild()->getLeaves(leaves);
        }

        if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
            leaves.push_back(const_cast<CollatorSPnode*>(this));
        }
        return leaves;
    }
    
	 
	  // return a reference to a container of SPSnodes
    // contents being all the nodes in left to right order
    std::vector<CollatorSPnode*>& CollatorSPnode::getAllNodes(vector<CollatorSPnode*>& allNodes) const
    {
        if (!isEmpty()) { // this is not empty
		  //if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
				//cout << nodeName << endl;
            allNodes.push_back(const_cast<CollatorSPnode*>(this));
        }
		  
		  //if children, recurse on the children
        if (hasLCwithBox()) {
            getLeftChild()->getAllNodes(allNodes);
        }

        if (hasRCwithBox()) {
            getRightChild()->getAllNodes(allNodes);
        }       
        return allNodes;
   }
	 
	 
	  // totalise the summaries in the collator node
    void CollatorSPnode::totaliseSummaries()
    {
		if (summary.size() > 1) {
			VecDbl tmp;
			tmp.push_back(getTotalSummary());
			summary.swap(tmp);
		}
		
		//recursion on the children
		if (hasLCwithBox()) {
			getLeftChild()->totaliseSummaries();
		}

		if (hasRCwithBox()) {
			getRightChild()->totaliseSummaries();
		}
   }

	// Marginalise, internal version
    // marginalise from given node downwards
    CollatorSPnode* CollatorSPnode::_marginalise(
			const CollatorSPnode * const rhs,
			const std::vector<int>& outDims)
    {
		CollatorSPnode* marginal = NULL;
		CollatorSPnode * newRC = NULL;
		CollatorSPnode * newLC = NULL;
				
		try {
			
			if (rhs != NULL ) { // if NULL we just return NULL
							
				ivector box = rhs->getBox();
				// will throw an exception if there is no box
				int dim = VecLen(box);
				int boxLB = Lb(box);
				
				int splitDim = rhs->getSplitDim();
				
				// deal with children first
				if (!rhs->isLeaf()) {
															
					newRC = _marginalise(rhs->getRightChild(), outDims);
					newLC = _marginalise(rhs->getLeftChild(), outDims);
				
				} // end isLeaf()
				
				// now deal with this node itself
				// if rhs a leaf, or did not split on any of the given dimensions
				// make a node that contracts rhs
				
				// iterator to vector element:
				std::vector<int>::const_iterator found 
						= find (outDims.begin(), outDims.end(), splitDim);
				
				if (found < outDims.end()) { // split on one of the outDims
					// so this will become the result of adding together
					// the two new children
					marginal = addPavings(newRC, newLC);
					// addPavings takes copies of the nodes to be added
					// so we need to destroy these ones
					delete newRC;
					newRC = NULL;
					delete newLC;
					newLC = NULL;
					
					marginal->totaliseSummaries(); 
				}
				else { // did not split on an outdim or is a leaf
					// have to contract rhs and then add the children
					VecDbl temp = rhs->getSummary();
					int l = rhs->getLabel();
					int newDims = dim - outDims.size();
					// for every missing dimension
					ivector newBox = ivector(newDims); 
					int index = Lb(newBox);
					int oldindex = boxLB;
					
					for (; oldindex <= Ub(box); oldindex++) {
						std::vector<int>::const_iterator fit 
						= find (outDims.begin(), outDims.end(), (oldindex - boxLB + 1));
						if (!(fit < outDims.end())) { // keep this one
							newBox[index] = box[oldindex];
							index++;
							
						}
					}
					
					double missingVol = Volume(box)/Volume(newBox);
					
					transform(temp.begin(), temp.end(), temp.begin(), 
					std::bind1st(multiplies<double>(), missingVol) );
					
					marginal = new CollatorSPnode(newBox, rhs->getLabel(), temp);
					marginal->totaliseSummaries(); 
										
					if (!rhs->isLeaf()) {
						marginal->nodeAddRight(newRC);
						marginal->nodeAddLeft(newLC);
					}
					newRC = NULL;
					newLC = NULL;
					
					} // finished else
			} // end if not NULL
			
			// clean up assuming no problems so far
			if (newRC != NULL) {
				try {
					delete newRC;
				}
				catch (exception& ee) {} // catch and swallow
			}
			if (newLC != NULL) {
				try {
					delete newLC;
				}
				catch (exception& ee) {} // catch and swallow
			}
			
			return marginal;
			
        } // end try
        catch (exception& e) {
            
            
            if (newRC != NULL) {
				try {
					delete newRC;
				}
				catch (exception& ee) {} // catch and swallow
			}
			if (newLC != NULL) {
				try {
					delete newLC;
				}
				catch (exception& ee) {} // catch and swallow
			}
			if (marginal != NULL) {
				try {
					delete marginal;
				}
				catch (exception& ee) {} // catch and swallow
			}
            
            const char* msg = e.what();
            throw SPnodeException("Error in _marginalise:\n " + std::string(msg));
        }
      
    }

// Marginalise
    // marginalise from given node downwards
    CollatorSPnode* CollatorSPnode::marginalise(
			const CollatorSPnode * const rhs,
			const std::vector<int>& reqDims)
    {
		CollatorSPnode * marginal = NULL;
		
		try { // throw exception if it is NULL or if dimensions incompatible
			if (rhs == NULL) {
				throw SPnodeException("Cannot marginalise null subpaving");
			}
			if (reqDims.empty()) {
				throw SPnodeException("No dimensions to marginalise on");
			}
			
			ivector box = rhs->getBox();
			int dim = VecLen(box);
			int boxLB = Lb(box);
			int boxUB = Ub(box);
			//each of the required dims must be there
			std::vector<int> sorted = reqDims;
			sort(sorted.begin(), sorted.end());
			
			// remove any duplicates
			vector<int>::iterator it = unique (sorted.begin(), sorted.end());
			sorted.resize( it - sorted.begin() );

			if ( (*(sorted.begin()) < 1)) {
				throw SPnodeException("Dimensions must be >= 1");
			}
			
			if (*(sorted.rbegin()) > boxUB - boxLB + 1)  {
				throw SPnodeException(
					"At least one dimension too large for subpaving box");
			}
			// could use min and max, but we want the not-req dims anyway
			std::vector<int> outDims;
			
			for (int i = 1; i <= dim; i++) {
				if (!(find(reqDims.begin(), reqDims.end(), i) < reqDims.end())) {
					// dim of box was not in reqDims 
					outDims.push_back(i);
				}
			}
			
			marginal = _marginalise(rhs, outDims);
			marginal->recursiveRename();
			return marginal;
		}
		catch (exception& e) {
			
			if (marginal != NULL) {
				try {
					delete marginal;
				}
				catch (exception& ee) {} // catch and swallow
			}
			const char* msg = e.what();
            throw SPnodeException("Error in marginalise:\n " + std::string(msg));
		}	
	}
	
	// find node's split dimension
    int CollatorSPnode::getSplitDim() const
    {
				
		try {
			
			int splitDim = -1;	
			
			if ( !isLeaf() ) {
			
				ivector box = getBox();
				int dim = VecLen(box);
				int boxLB = Lb(box);
					
				ivector boxChild = getRightChild()->getBox();
										
				int index = 0;
				int boxChildLB = Lb(boxChild);
				while (splitDim < 1 && index < dim) {
					if ((Inf(box[boxLB + index]) 
							!= Inf(boxChild[boxChildLB + index]))
						||
						(Sup(box[boxLB + index]) != 
							Sup(boxChild[boxChildLB + index]))) {
								// found splitDim
								splitDim = boxLB + index;
					}
					index ++;
				} // end while
							
				if (splitDim < 0) {
					throw SPnodeException("Cannot find split dimension");
				}
			} // end isLeaf
								
			return splitDim;
			
        } // end try
        catch (exception& e) {
            
            const char* msg = e.what();
            throw SPnodeException("Error in findParentSplitDim:\n " + std::string(msg));
        }
      
    }

    // get the total of the summary in the collator node
    double CollatorSPnode::getTotalSummary() const
    {
		double summ = 0;
		summ = accumulate(summary.begin(), summary.end(), summ);

		return summ;
	}
 
	 // get the average of the summary in the collator node
    double CollatorSPnode::getTotalSummaryAv() const
    {
		double summ = 0;
		summ = accumulate(summary.begin(), summary.end(), summ)/(summary.size()*1.0);

		return summ;
	 }
 
 // Jenny addition for Gloria's convergence work
	// take a container and return the same container, which has been
	// cleared (if necessary) and re-filled with 
	// L1-distances-to-average, one for each histogram in collation
	RealVec& CollatorSPnode::getL1DistancesToAverage(RealVec& container) const
	{
        CollatorSPnode* differences = NULL;
        
        // take this collation
        try {
            // make a tree of differences between this and average over this
            differences = makeDifferencesToAveragePaving();
			
            size_t n = summary.size();
            dotprecision emptyDP(0.0);
            VecDotPrec areaAcc(n, emptyDP); // n copies of empty dot prec

            areaAcc = differences->getLeafNodeAbsAreaAccumulations(areaAcc);
			// one L1 distance-to-average summary, in the form of a dot precision, for 
			// each histogram 

			RealVec temp(n);
			temp.swap(container); // clear container and minimize memory
			
			VecDotPrecIt it;
            for (int i = 0; i < n; ++i) {
				// put rounded L1 diff into the container
				container.at(i) = rnd( areaAcc.at(i) ); // round to nearest
            }

            delete differences;
            differences = NULL;
			
			return container;
        }
        catch (exception& e) {
            if (NULL != differences) {
                delete differences;
                differences = NULL;
            }
            throw;
        }

    }
 
    //gat41
    // Find the nodes that fulfill the Scheffe condition for rows against 
	// columns in the Yatracos growing matrix. True if condition is fulfilled.
   bool CollatorSPnode::getScheffeNode(int theta1, int theta2)
   { 
     //cout.precision(20);
     //cout << "Checking for Scheffe set at node: " << getNodeName() << endl;
     //cout << "Theta1: " << theta1 << "\t" << "Theta2: " << theta2 << endl;     
     //cout << summary[theta1] << "\t" << summary[theta2] << endl;
     
		//check that this is an ordered pair theta1 < theta2
		if (theta1 < theta2) {
			if ((summary[theta1] > summary[theta2])) {
	        //cout << getNodeName() << " is an element of the Scheffe set.****" << endl;
	        return true;          
			} 
			else { return false; }
		}
		else {
			cerr << "theta1 must be less than theta2." << endl;
			exit(0);
		}
   } 	  
    
    
    //gat41
    // get delta for a specific theta
   double CollatorSPnode::getNodeDelta(int thisTheta, size_t sizeColl)
   { 
     //cout << "get delta for " << nodeName << endl;
     // get empirical measure of the training data
     double muTrain = summary[thisTheta] * nodeVolume();
     //cout << "summary: " << summary[thisTheta] << "\t muTrain: " << muTrain << endl;
      
     // get empirical measure of the validation data      
     double muValid = summary[sizeColl-1] * nodeVolume();
    //cout << "summary: " << summary[sizeColl-1] << "\t muValid: " << muValid << endl;

     double delta= muTrain - muValid;
		//cout << "Delta: " << delta << endl; 

     return delta; 

   } // end of function getNodeDelta
    
    //gat41
    // Find the nodes that fulfill the Scheffe condition for rows against 
	// columns in the Yatracos growing matrix. True if condition is fulfilled.
   bool CollatorSPnode::nodeCheckRowSummary(int theta, int k)
   { 
     //cout << "checking for Yat at node: " << getNodeName() << endl;
     //cout << "theta: " << theta << "\t" << "k: " << k << endl;     
     //cout << summary[theta] << "\t" << summary[k] << endl;
      
      if ((summary[theta] > summary[k])) {
		//cout << "height at " << theta << " larger than height at " << k << endl;
            return true;          
      } // end of filling up rows
      else { return false; }
   }   
	  
   // Find the nodes that fulfill the Scheffe condition for columns against 
	// rows in the Yatracos growing matrix. True if condition is fulfilled
   bool CollatorSPnode::nodeCheckColSummary(int theta, int k)
   { 
	 //  cout << "checking for Yat at node: " << getNodeName() << endl;     
    //  cout << "theta: " << theta << "\t" << "k: " << k << endl;
    //  cout << summary[theta] << "\t" << summary[k] << endl;
      if ((summary[k] > summary[theta])) {
	 // cout << "height at " << k << " larger than height at " << theta << endl;
         return true;
      }   
		else { return false; }    
	}
    
    
    
    //gat41
    //get the Yatracos set for a particular pair.
void CollatorSPnode::getYatSet(
			set<CollatorSPnode*, less<CollatorSPnode*> > & YatSetRow, 
			set<CollatorSPnode*, less<CollatorSPnode*> > & YatSetCol, 
			size_t cand1, size_t cand2)
{
	//iterate through the leaves in both candidate histograms to get the 
	//Yatracos set
   if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
      bool rowInd = nodeCheckRowSummary(cand1, cand2);
  		bool colInd = nodeCheckColSummary(cand1, cand2);
     	// insert the node YatSet if return true
		if (rowInd) { 
			//cout << "inserting " << getNodeName() << " into YatSetRow" << endl; 
			YatSetRow.insert(&(*this));
		}
		if (colInd) { 
			//cout << "inserting " << getNodeName() << " into YatSetCol" << endl; 
			YatSetCol.insert(&(*this));
		}
   }
 
   //recurse on the children
   if (getLeftChild()!=NULL) {
         getLeftChild()->getYatSet(YatSetRow, YatSetCol, cand1, cand2);
   }
   if (getRightChild()!=NULL) {
         getRightChild()->getYatSet(YatSetRow, YatSetCol, cand1, cand2);
   }
} // end of getYatset

    
    //gat41
    //get the Scheffe set for a particular pair.	
	void CollatorSPnode::getScheffeSet(
			set<CollatorSPnode*, less<CollatorSPnode*> > & ScheffeSet, 
			size_t cand1, size_t cand2)
	{
		//iterate through the leaves in both candidate histograms to get the 
		//Yatracos set
	   if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
	      bool ind = getScheffeNode(cand1, cand2);
	  			// insert the node YatSet if return true
			if (ind) { 
				ScheffeSet.insert(&(*this));
			}
		}
			
	   //recurse on the children
	   if (getLeftChild()!=NULL) {
	         getLeftChild()->getScheffeSet(ScheffeSet, cand1, cand2);
	   }
	   if (getRightChild()!=NULL) {
	         getRightChild()->getScheffeSet(ScheffeSet, cand1, cand2);
	   }
	}
	
	//src_trunk_0701
	/*
	void CollatorSPnode::swapCollator(CollatorSPnode& spn) //throw() // don't hide base class version
	{
		/* theBox, parent, leftChild,
        rightChild and nodeName are inherited from base class */
		/*SPnode::swap(spn); // use the base version
		
		std::swap(rangeCollection, spn.rangeCollection);     
	}*/

    // ----------------- non member tools functions ----------------------

    //Output all boxes in collator
    std::ostream & operator<<(std::ostream &os,
                        const CollatorSPnode& spn)
    {
        os << spn.nodesAllOutput(os, 1) << std::endl;
        return os;
    }

   //compare total summaries
	bool nodeCompTotalSummary(const CollatorSPnode * const lhs,
                            const CollatorSPnode * const rhs)
	{
		return (lhs->getTotalSummary() < rhs->getTotalSummary());
	}

	//compare average of total summaries 
	bool nodeCompTotalSummaryAv(const CollatorSPnode * const lhs,
                            const CollatorSPnode * const rhs)
	{
		return (lhs->getTotalSummaryAv() < rhs->getTotalSummaryAv());
	}
	
/*	//src_trunk_0701
	// Full specializations of the templates in std namespace can be added in std namespace.
	template <>
	void std::swap(CollatorSPnode & s1, 
				CollatorSPnode & s2) // throw ()
		{
			s1.swapCollator(s2);
		}
*/
} // end namespace subpavings


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

/*!/ \file:     spsnode.cpp
\brief SPSnode (StatsSubPaving) definitions
*/

#include "spsnode.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <exception>

// include fstream so as to be able to output a file
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>

// format manipulation on streams
#include <iomanip>

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"

// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

// to use spsnode splitting classes
#include "splitdecisionobj.hpp"

// to get determinant and inverse of a matrix
#include "gsl/gsl_linalg.h"

// to access gsl_matrix elements
#include "gsl/gsl_matrix.h"

// to perform vector-matrix operations
#include "gsl/gsl_blas.h"

#include "assert.h"

//src_trunk_0701
#include "subpaving_exception.hpp"
#include "sptools.hpp"

using namespace subpavings;
using namespace std;

    // ---------------------- private member functions -------------------

    // recalculate the counter and accumulated sum
    // and accumulated sumproducts
    void SPSnode::recalculateStats(rvector& newdata) const
    {
        
        counter++;  // update the counter

        if (!countsOnly) {

            recalculateSums(newdata); // update the sums

            recalculateSumProducts(newdata); // update the sumproducts
        }
    }


    // recalculate the accumulated sum
    void SPSnode::recalculateSums(rvector& newdata) const
    {
        if (dpSums.empty()) {   //nothing in the sums yet
            // reserve space in dpSums for all elements of the mean
            dpSums.reserve(dimension);

            // for each dimnsn of data, initialise element
            for (size_t i = 0; i< dimension; i++) {
                dotprecision dp;
                dp = 0.0;
                dpSums.push_back(dp);
            }
        }

        // make a dot precision variable out of the ith element
        // of the rvector of new data and store in dpSums
        for (size_t i = 1; i< dimension + 1; i++) {
            // rvectors indexed 1 to n, vectors indexed 0 to n-1
            accumulate(dpSums[i-1], newdata[i], 1.0);
        }

    }

    // recalculate the accumulated sumproducts
    void SPSnode::recalculateSumProducts(rvector& newdata) const
    {
        /* the sumproducts can be thought of as an nxn matrix,
        which is implemented here as a nxn element vector of
        dotprecision variables, using row-major order.
        Ie the m-th element (m = 0, . . . nxn-1) is in row floor(m/n)
        and column m-rowxn in the matrix configuration.
        Or, the sumproduct of elements i and j in an rvector,
        i,j = 0,...,n-1, is element m=(ixn+j) of the sumproducts
        vector. */

        if (dpSumProducts.empty()) {    //nothing there yet
            // reserve space for all elements
            dpSums.reserve(dimension*dimension);

            // for each dimnsn^2 of data, initialise element
            for (size_t i = 0; i< (dimension*dimension); i++) {
                dotprecision dp;
                dp = 0.0;
                dpSumProducts.push_back(dp);
            }
        }

        // make a dot precision variable out of the ith element
        // and jth element of the of the rvector of new data and
        // store in dpSumProducts.
        for (size_t i = 1; i < dimension + 1; i++) {
            // only need to do columns 1 to i because of symmetry
            for (size_t j = 1; j< i + 1; j++) {

                size_t index = (i-1)*dimension + (j-1);
                // rvectors indexed 1 to n
                accumulate(dpSumProducts[index],
                        newdata[i], newdata[j]);

                //if not on the diagonal of the matrix,
                // we can also fill in the symmetric element
                if (i!=j) {
                    size_t sym_index = (j-1)*dimension
                        + (i-1);
                    dpSumProducts[sym_index] =
                        dpSumProducts[index];
                } // end if
            }// end j-loop
        }// end i-loop

        // sumproducts has been updated for new datapoint
    }






    // Only expand the node - no reallocation of data
    // add two sibling nodes to this provided that this is a leaf
    // comp argument is passed to Upper() and Lower()
    // these functions split box in half normal to dimension set by comp
    void SPSnode::nodeExpansionOnly(int comp)
    {
        try
        {
            // only do something if this SPSnode is a leaf
            if (isLeaf()) {
                // ivectors to become boxes for new children
                ivector lC, rC;

                // Call Lower() and Upper() to put the split
                // boxes into lC and rC respectively
                Lower(getBox(), lC, comp);
                Upper(getBox(), rC, comp);

                // when making new children, use constructor
                // that will give space indication (for data)
                // of the size of this node's dataItrs
                size_t space = dataItrs.size();

                nodeAddLeft(new SPSnode(lC,
                                        space, countsOnly, label));

                nodeAddRight(new SPSnode(rC,
                                        space, countsOnly, label));

                //name the new children
                getLeftChild()->setNodeName(nodeName + "L");
                getRightChild()->setNodeName(nodeName + "R");

                // store the split dimension in this
                splitDim = comp;

                // store the split value in this
                // the split value is the infinum of interval
                // of right child box for dimension split on
                splitValue = _double(Inf(
                    ((getRightChild())->getBox())[comp]));
            }
        }

        catch (bad_alloc&)
        {
            std::cerr << "Error allocating memory in "
                << "SPSnode::nodeExpansionOnly()" << std::endl;
            throw;
        }

    }


    // split data between two new children
    // using an SplitDecisionObj to see if the children should be further split
    void SPSnode::splitData(const SplitDecisionObj& boolTest)
    {

        // check that both children exist
        if (!hasLCwithBox() || !hasRCwithBox()) {
            string msg = "Cannot split data when there are not two ";
            msg += " children";
            throw SPnodeException(msg);
        }

        NodeDataItr dataItr; // iterator

        //divvie the data up amongst the children
        for (dataItr = dataItrs.begin();
            dataItr!= dataItrs.end(); dataItr++) {
            BigDataItr newItr = *dataItr;

            //calls insertOneFind on the children of this node
            // so stats are not recalculated for this node itself
            SPSnode* reinsertedInto = NULL;

            if(rightChild!=NULL && !rightChild->isEmpty()) {

                reinsertedInto =
                    (getRightChild())->insertOneFind(
                    newItr, ON_RIGHT, boolTest);
            }

            // only try the left if it's not on the right
            if(reinsertedInto==NULL && leftChild!=NULL
            && !leftChild->isEmpty()) {

                reinsertedInto =
                    (getLeftChild())->insertOneFind(
                    newItr, ON_LEFT, boolTest);
            }

        }

        clearData();         //clear the data in this node
    }


    // Print the data in a node if any
    std::ostream& SPSnode::nodeDataPrint(std::ostream &os) const
    {
        if (!dataItrs.empty()) {

            NodeDataItr dataItr;

            os << "Data is" << std::endl;
            for (dataItr = dataItrs.begin();
                dataItr!= dataItrs.end(); dataItr++) {

                BigDataItr bigIt = *dataItr;
                rvector theData = *bigIt;

                for (size_t i = 1; i < dimension + 1; i++) {
                    os << label; // print the label
                    os << "  " << theData[i]; // print data
                }   // end loop through data elements

                os << std::endl;

            } // end loop through data container
        } // end if counter > 0
        // if no data, ie counter = 0, then just return os unaltered

        return os;
    }

    // Print the mean of the data in a node
    std::ostream& SPSnode::nodeMeanPrint(std::ostream &os) const
    {

        if ((counter > 0) && !countsOnly) {

            os << "Mean is ";

            // loop through the elements in the dpSums vector
            for (size_t i = 0; i< dimension; i++) {
                // default cxsc rounding of dotprecision
                // to rnd_next
                os << "  " << (rnd(dpSums[i])/(1.0*counter));

            }// end loop through the elements in dpSums

            os << std::endl;

        } // end if
        // if no data, ie counter = 0, or if we are only keeping counts
        // then just return os unaltered

        return os;

    }

    // Print the variance covariance matrix of the data in a node
    std::ostream& SPSnode::nodeVarCovarPrint(std::ostream &os) const
    {
        if ((counter > 0) && !countsOnly) {

            RealVec varCovar;
            varCovar = getVarCovar(varCovar);

            /* element k in the vector representing the
            variance-covariance matrix corresponds to
            row k/n, (row 0 to n-1) and column k-row*n (col 0 to n-1)
            in a matrix view variance-covariance */

            os << "Variance Covariance is " << std::endl;

            // loop through the elements and print as matrix
            for (size_t i = 0; i < dimension; i++) {
                for (size_t j = 0; j < dimension; j++) {
                    os << "  " << varCovar[(i*dimension)+j];
                }
                os << std::endl;
            }
        }
        return os;

    }

    // Print the details of a single leaf node, using tab delimiters
    std::ostream& SPSnode::leafOutputTabs(std::ostream &os) const
    {
        int prec = 5; // precision for output

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy of theBox

            // output the node name, nodeVolume, counter
            os << nodeName;
            os << "\t" << nodeVolume();
            os << "\t" << counter;
            // followed by the intervals of box using Inf and Sup
            // ie unlike cxsc output, there is no [  ] around them

            streamsize oldPrec = os.precision();
            os << setprecision(prec);

            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
            os << setprecision(oldPrec);

        }
    }

    // Print the details of a single leaf node, using tab delimiters
    // including EMP contributions and changes if split
    std::ostream& SPSnode::leafOutputTabsWithEMPs(const size_t bigN,
                            std::ostream &os, const int prec) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy of theBox

            // output the name, (not label), nodeVolume, counter
            os << nodeName;
            //os << label;
            os << "\t" << nodeVolume();
            os << "\t" << counter;
            // EMP contributions and changes if split
            os << "\t" << getEMPContributionCOPERR(bigN);
            os << "\t" << rnd(getSplitChangeEMPCOPERR(bigN));
            os << "\t" << getEMPContributionAIC(bigN);
            os << "\t" << rnd(getSplitChangeEMPAIC());

            // followed by the intervals of box using Inf and Sup
            // ie unlike cxsc output, there is no [  ] around them
            streamsize oldPrec = os.precision();
            os << setprecision(prec);

            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
            os << setprecision(oldPrec);

        }
    }

    // Print the details of a single leaf node, using tab delimiters
    // includes the height = n/(N*vol) where n is count in this leaf node,
    // N is count over whole histogram, vol is volume of this leaf node
    std::ostream& SPSnode::leafOutputTabsWithHistHeight(const size_t bigN,
                            std::ostream &os, const int prec) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy of theBox
            double vol = nodeVolume();

            // output the node name, nodeVolume, counter, counter/(bigN * vol)
            os << nodeName;
            os << "\t" << vol;
            os << "\t" << counter;
            os << "\t" << counter/(vol * bigN);
            // followed by the intervals of box using Inf and Sup
            // ie unlike cxsc output, there is no [  ] around them
            streamsize oldPrec = os.precision();
            os << setprecision(prec);

            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
            os << setprecision(oldPrec);

        }
    }

    // Print the details of a single leaf node, using tab delimiters
    // including EMP contributions and changes if split
    std::ostream& SPSnode::leafOutputTabsWithHistHeightAndEMPs(
                                    const size_t bigN, std::ostream &os,
                                    const int prec) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy of theBox
            double vol = nodeVolume();

            // output the node name, nodeVolume, counter, counter/(n * vol)
            os << nodeName;
            //os << label;
            os << "\t" << nodeVolume();
            os << "\t" << counter;
            os << "\t" << counter/(vol * bigN);
            // EMP contributions and changes if split
            os << "\t" << getEMPContributionCOPERR(bigN);
            os << "\t" << rnd(getSplitChangeEMPCOPERR(bigN));
            os << "\t" << getEMPContributionAIC(bigN);
            os << "\t" << rnd(getSplitChangeEMPAIC());

            // followed by the intervals of box using Inf and Sup
            // ie unlike cxsc output, there is no [  ] around them
            streamsize oldPrec = os.precision();
            os << setprecision(prec);

            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
            os << setprecision(oldPrec);

        }
    }


    // gather up all the data in children of a node
    NodeData& SPSnode::gatherData(NodeData& container, SPSnode * spn)
    {
        if (!spn->isLeaf()) {
            if (spn->hasLCwithBox()) {
                container =
                    gatherData(container,
                            spn->getLeftChild());
            }
            if (spn->hasRCwithBox()) {
                container =
                    gatherData(container,
                            spn->getRightChild());
            }
        }
        else { // is a leaf
            // copy data from spn's dataItrs into temp container
            container.insert(container.end(),
                            (spn->dataItrs).begin(),
                            (spn->dataItrs).end());
        }

        return container;
    }

    // set split dimension and split value in parent when children grafted on
    void SPSnode::setSplits()
    {
        // set the split dimension and split value for this box
        // based on the children which have been added
        ivector childBox;
        bool alreadyDone = false;

        if (hasRCwithBox()) {
            childBox = getRightChild()->getBox();
        }
        else if (hasLCwithBox()) {
            childBox = getLeftChild()->getBox();
        }
        int pLb = Lb(*theBox); // parent box lower bound
        int dim = Ub(*theBox) - Lb(*theBox) + 1;
        int cLb = Lb(childBox); // child box lower bound (should be = pLb)

        if ( splitDim != -1 &&
            (splitValue == Inf(childBox[splitDim - pLb + cLb])
            || splitValue == Sup(childBox[splitDim - pLb + cLb])))
                alreadyDone = true;

        if (!alreadyDone) {
            int d = 1;
            splitDim = -1;
            while ((d <= dim) && (splitDim == -1)) {
                if (diam(childBox[d + cLb - 1]) < diam((*theBox)[d])) {
                    splitDim = d + pLb - 1; // the split dimension
                }
                d++;
            }
            // split value is bottom of right child box on dth dim
            if (hasRCwithBox()) {
                splitValue = Inf(childBox[splitDim - pLb + cLb]);
            }
            // else split value is top of left child box on dth dim
            else if (hasLCwithBox()) {
                splitValue = Sup(childBox[splitDim - pLb + cLb]);
            }
        }
    }


    // add two non-minimal pavings in a union operation,
    // return a pointer to a new non-minimal paving
    // but with no data attached to it - up to the manager to add data
    // label will be 0
    SPSnode* SPSnode:: unionNoData(const SPSnode * const lhs,
                            const SPSnode * const rhs)
    {
        SPSnode* newNode = NULL;

        bool done = false;  // indicator for done adding

        try {

            if (lhs == NULL && rhs == NULL) done = true; // we will return NULL

            // if the lhs is null or has no box, return a tree or node based on rhs
            if (!done && (lhs==NULL || ((lhs != NULL) && (lhs->isEmpty())))) {

                newNode = SPSnode::strippedConstructor(rhs);
                done = true;
            }

            // if the rhs is null or has no box, return a tree or node based on lhs
            if (!done && (rhs==NULL || ((rhs != NULL) && (rhs->isEmpty())))) {

                newNode = SPSnode::strippedConstructor(lhs);
                done = true;
            }

            // by now, if we are not done, both pavings are not null and both have boxes
            // we assume that the boxes are the same

            // we have to check who has children

            // if both are leaves we can just return a node based on say lhs
             // if only rhs is leaf, lhs is not a leaf, return a node based on lhs
            if (!done && rhs->isLeaf()) {
                newNode = SPSnode::strippedConstructor(lhs);
                done = true;
            }

            // if only lhs is leaf, rhs is not a leaf, return a node based on rhs
            if (!done && lhs->isLeaf() && !rhs->isLeaf()) {
                newNode = SPSnode::strippedConstructor(rhs);
                done = true;
            }

            // if neither are leaves
            if (!done && !lhs->isLeaf() && !rhs->isLeaf()) {
                // make a node based on one of them, and add on the results of
                // recursing on the children
                ivector* newPermBox = new ivector(lhs->getBox());
                newNode = new SPSnode(*newPermBox);
                newNode->nodeAdoptRight(unionNoData(lhs->getRightChild(),
                                                            rhs->getRightChild()));
                newNode->nodeAdoptLeft(unionNoData(lhs->getLeftChild(),
                                                            rhs->getLeftChild()));
            }
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return newNode;

    }
    
    //src_trunk_0701
	void SPSnode::_reshapeToUnion(const SPnode * const other)
{
	SPnode::_reshapeToUnion(other);
}

/* reshape this tree to have union of this shape and shape of other
 no checks on boxes since this should be redundant if used by unionNoData...
*/
bool SPSnode::_reshapeToUnion(const SPnode * const other,
							size_t minChildPoints,
							const std::string& errorFilename)
{
	// indictator for being able to do union exactly
	bool success = true;
	
	if ( other != NULL && !(other->isEmpty()) ) {

		// this is not a leaf, other is a leaf
		if (!isLeaf() && other->isLeaf()) {

			// no need to do anything
		}

		// this is a leaf, other is not a leaf
		if (isLeaf() && !other->isLeaf()) {

			//we need to expand this
			if (isSplittableNode(minChildPoints)) nodeExpand();
			else {
				success = false;
				// log file
				std::string line = "Could not split " + getNodeName() 
				 + " because of minChildPoints";
				
				outputFile(errorFilename, line); 
				
			}
			
		}

		// now recurse on the children if both have children
		// note - it won't go here is !success because still isLeaf()
		if (!isLeaf() && !other->isLeaf()) {
			success = getLeftChild()->_reshapeToUnion(
					other->getLeftChild(), minChildPoints, errorFilename);
			success = getRightChild()->_reshapeToUnion(
					other->getRightChild(), minChildPoints, errorFilename)
					&& success;
		}
		
		return success;
	}
}
//src_trunk_0701

    // ------------------------ public member functions -----------------

    // Default constructor
    SPSnode::SPSnode() :  counter(0), splitDim(-1), splitValue(0.0),
                                countsOnly(true)
    {
        try {
            //invokes the base class default constructor
            // then does additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            // reserve space
            spaceIndication = static_cast<size_t>(defaultMaxPts);
            // not sure whether to do this or not - leave for the moment
            dataItrs.reserve(spaceIndication);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            std::cerr << msg << std::endl;
            throw;
        }

    }

    // initialised constructor, initialised with one ivector for the box
    // and optionally with lab for label, defaults to 0 (see declaration)
    // countsOnly will default to false
    SPSnode::SPSnode(ivector& v, int lab) : SPnode(v, lab),
        counter(0), splitDim(-1), splitValue(0.0), countsOnly(true)
    {
        try {
            //invokes the base class constructor with ivector & label
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            spaceIndication = static_cast<size_t>(defaultMaxPts);
            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with one ivector for the box
    // and value for countsOnly
    // and optionally with lab for label, defaults to 0 (see declaration)
    SPSnode::SPSnode(ivector& v, bool cntOnly, int lab) : SPnode(v, lab),
        counter(0), splitDim(-1), splitValue(0.0), countsOnly(cntOnly)
    {
        try {
            //invokes the base class constructor with ivector & label
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            spaceIndication = static_cast<size_t>(defaultMaxPts);
            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with one ivector for the box
    // and max for spaceIndication
    // and optionally with lab for label, defaults to 0 (see declaration)
    // countsOnly defaults to false
    SPSnode::SPSnode(ivector& v, size_t max, int lab) :
        SPnode(v, lab),
        spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
        countsOnly(true)
    {
        try {
            //invokes the base class constructor with ivector argument
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication+1);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }


    // initialised constructor, initialised with one ivector for the box
    // and max for spaceIndication and value for countsOnly
    // and optionally with lab for label, defaults to 0 (see declaration)
    SPSnode::SPSnode(ivector& v, size_t max, bool cntOnly, int lab) :
        SPnode(v, lab),
        spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
        countsOnly(cntOnly)
    {
        try {
            //cout << "node constructor" << "\t" << cntOnly << endl;
            //invokes the base class constructor with ivector argument
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication+1);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with a LabBox (labeled box)
    // and a max for spaceIndication
    // and optionally with cntOnly for countsOnly, defaults to false
    SPSnode::SPSnode(LabBox& lb, size_t max, bool cntOnly) : SPnode(lb),
        spaceIndication(max), counter(0), splitDim(-1), splitValue(0.0),
        countsOnly(cntOnly)
    {
        try {

            //invokes the base class constructor with LabBox argument
            //and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication+1);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with a LabBox (labeled box)
    // and optionally with cntOnly for countsOnly, defaults to false
    SPSnode::SPSnode(LabBox& lb, bool cntOnly) : SPnode(lb), counter(0),
        splitDim(-1), splitValue(0.0), countsOnly(cntOnly)
    {
        try {
            //invokes the base class constructor with LabBox argument
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            spaceIndication = static_cast<size_t>(defaultMaxPts);
            //reserve space - not sure if important - leave for moment
            dataItrs.reserve(spaceIndication);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    //Copy constructor
    // copies from given node downwards
    SPSnode::SPSnode(const SPSnode& other) : SPnode(*(other.theBox),
        other.label), spaceIndication(other.spaceIndication),
        counter(other.counter), dpSums(other.dpSums),
            dpSumProducts(other.dpSumProducts), splitDim(other.splitDim),
        splitValue(other.splitValue), countsOnly(other.countsOnly)
    {
        try {
            //reserve space
            dataItrs.reserve((other.dataItrs).size());
            //copy dataItrs from other to this
            dataItrs = other.dataItrs;
            nodeName = other.nodeName;

            //recursion on the children
            if (other.leftChild) {
                nodeAddLeft(new SPSnode(*(other.getLeftChild())));
            }
            else leftChild=NULL;

            if (other.rightChild) {
                nodeAddRight(new SPSnode(*(other.getRightChild())));
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

    //copy assignment operator
    //copies from this node downwards
    SPSnode& SPSnode::operator=(const SPSnode& rhs)
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

            spaceIndication = rhs.spaceIndication;

            counter = rhs.counter;
            dpSums = rhs.dpSums;
            dpSumProducts = rhs.dpSumProducts;
            splitDim = rhs.splitDim;
            splitValue = rhs.splitValue;
            countsOnly = rhs.countsOnly;

            //reserve space
            dataItrs.reserve((rhs.dataItrs).size());
            //copy dataItrs from other to this
            dataItrs = rhs.dataItrs;

            //recursion on the children
            if (rhs.leftChild) {
                nodeAddLeft(new SPSnode(*(rhs.getLeftChild())));
            }
            else leftChild=NULL;

            if (rhs.rightChild) {
                nodeAddRight(new SPSnode(*(rhs.getRightChild())));
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


    // A stripping 'constructor'
    // copies from given node downwards but with no data
    SPSnode* SPSnode::strippedConstructor(const SPSnode * const other)
    {
        SPSnode* newNode = NULL;
        try {
            if (other != NULL) {
                if (other->isEmpty())
                    newNode = new SPSnode;
                else {
                    ivector* newBox = new ivector(other->getBox());
                    newNode = new SPSnode(*newBox);
                    newNode->splitDim = other->splitDim;
                    newNode->splitValue = other->splitValue;
                }

                newNode->nodeName = other->nodeName;
                newNode->label = 0;
                newNode->countsOnly = true;

                if (other->getLeftChild() != NULL)
                    newNode->nodeAddLeft(strippedConstructor(other->getLeftChild()));
                if (other->getRightChild() != NULL)
                    newNode->nodeAddRight(strippedConstructor(other->getRightChild()));
            }
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

        return newNode;
    }



    // Accessor for the counter.
    size_t SPSnode::getCounter() const
    { return counter; }

    // Accessor for the split dimension.
    int SPSnode::getSplitDim() const
    { return splitDim; }

    // Accessor for the split value.
    real SPSnode::getSplitValue() const
    { return splitValue; }

    // Accessor for the countsOnly value.
    real SPSnode::getCountsOnly() const
    { return countsOnly; }

    // Accessor for the node's data collection.
    // Returns a copy of the node's collection of iterators to the big data set.
    NodeData SPSnode::getData() const
    { return dataItrs; }

    // Clears the node's data collection.
    void SPSnode::clearData() const
    { dataItrs.clear(); }

    // Clears the node's data collection and counters
	 void SPSnode::makeEmptyNode()
	 { dataItrs.clear();
	   counter = 0;
	 }
	 
    // Accessor for the parent of a node
    // Returns a copy of the pointer to parent node.
    SPSnode* SPSnode::getParent() const
    { return (SPSnode*) parent; }

    // Accessor for the left child.
    // Returns a copy of the pointer to leftChild node cast to this node type
    SPSnode* SPSnode::getLeftChild() const
    { return (SPSnode*) leftChild; }

    // Accessor for the right child
    //Returns a copy of the pointer to rightChild node cast to this node type
    SPSnode* SPSnode::getRightChild() const
    { return (SPSnode*) rightChild; }
    
    //--src_trunk_0701
    bool SPSnode::checkTreeStateLegal() const 
{
	return SPnode::checkTreeStateLegal();
}

bool SPSnode::checkTreeStateLegal(size_t minChildPoints, 
								double minVol) const
{
	// check current state is legal by looking at everything not a leaf
	bool legal = true;
	if ( !isLeaf() ) {
		legal = isSplittableNode(minChildPoints, minVol);
		if (legal && hasLCwithBox() ) {
				legal = 
					getLeftChild()->checkTreeStateLegal(minChildPoints,
														minVol);
		}
		if (legal && hasRCwithBox() ) {
				legal = 
					getRightChild()->checkTreeStateLegal(minChildPoints,
														minVol);
		}
	}
	
	return legal;
}

bool SPSnode::checkTreeStateLegal(size_t minChildPoints)
{
	// check current state is legal by looking at everything not a leaf
	bool legal = true;
	if ( !isLeaf() ) {
		legal = isSplittableNode(minChildPoints);
		if (legal && hasLCwithBox() ) {
				legal = 
					getLeftChild()->checkTreeStateLegal(minChildPoints);
		}
		if (legal && hasRCwithBox() ) {
				legal = 
					getRightChild()->checkTreeStateLegal(minChildPoints);
		}
	}
	
	return legal;
}
//--src_trunk_0701
    
    
    
    //src_trunk_0701
    bool SPSnode::isSplittableNode() const
{
	return SPnode::isSplittableNode();
} 


//src_trunk_0701
bool SPSnode::isSplittableNode(size_t minChildPoints, 
								double minVol) const 
{
    
	bool retValue = (nodeVolume() >= minVol);
//	cout << nodeVolume() << "\t" << minVol << endl;
	
	//cout << nodeVolume() << "\t" << minVol << endl;
	
	if (retValue) {
		retValue = isSplittableNode(minChildPoints);
	}
	else {
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode: node failed vol test" << endl;
		#endif
	}
	return retValue;
}


//src_trunk_0701
bool SPSnode::isSplittableNode(size_t minChildPoints) const
{
    bool retValue = isSplittableNode();  //basic check
	if (!retValue) {
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode: node failed basic is splittable test" << endl;
		#endif
		#ifdef DEBUG_MCMC_SPLIT_FAIL
			cout << "Failed isSplittableNode: I am " << nodeName << endl;
			{
				ivector box = getBox();
				interval maxD = box[MaxDiamComp(box)];
	
				cout << cxsc::SaveOpt;
				cout << Scientific << SetPrecision(35,30);
				cout << "interval to be split is " << maxD << endl;
				cout << cxsc::RestoreOpt;
			}
		#endif
	}
	#ifdef DEBUG_CHECK_NODE_COUNT
		cout << "isSplittableNode minChildPoints = " << minChildPoints << endl;
	#endif
	if (retValue && minChildPoints > 0) {
		retValue = false; // need to retest
	
		size_t  minChildCount = getMinChildCountIfSplitNEW();
		
				
		if ( (counter >= minChildPoints) &&
			((minChildCount == 0) || (minChildCount >= minChildPoints)) ) {
				retValue = true;
			}
		#ifdef DEBUG_CHECK_NODE_COUNT
			cout << "isSplittableNode minChildCount = " << minChildCount << endl;
			cout << "(minChildCount >= minChildPoints) = " << (minChildCount >= minChildPoints) << endl;
	
			cout << "isSplittable = " << retValue << endl;
		#endif
	}
    return retValue;
}
    
    

  // get the number of datapoints currently associated with this which would
// be associated with the new left child if this node were to split
// remember that the left child's box is open at the split
size_t SPSnode::getLeftCountIfSplit() const
{
	
	// first find what the dimension for the split would be 
	// if the split were made
	// right hand child's box would be if that child
	// were to be created
	cxsc::ivector box = getBox();
	
	int split = MaxDiamComp(box);
	
	cxsc::real midSplit = cxsc::mid(box[split]);

	// left child would have everything up to but not including
	// midSplit, on the split dimension
	size_t leftCount = 0;
	NodeDataItr it;

	for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
		// DataItrs is a container of iterators to a BigDataCollection
		// increment rightCount if the point is in rC
		if(  (**it)[split] < midSplit ) leftCount++;
	}

	return leftCount;
}

	// Smallest number of points in either child if this was split.
	size_t SPSnode::getMinChildCountIfSplit() const
	{
		try {
			size_t min = getLeftCountIfSplit();
			if ((counter - min) < min) min = counter - min;
			return min;
		}
		catch (exception& e) {
			string msg = string(e.what());
			throw SPnodeException("Error in getMinChildCountIfSplit:\n" + msg);
		}
	}

// Smallest number of points in either child if this was split.
size_t SPSnode::getMinChildCountIfSplitNEW() const
{
	size_t min = 0;
	
	if (isLeaf()) {
		min = getLeftCountIfSplit();
	}
	else {
		min = getLeftChild()->getCounter();
	}
	
	if ((counter - min) < min) min = counter - min;
	
	return min;
	
}


	// get the number of datapoints currently associated with this which would
	// be associated with the left and right children of a new right child
	// if this node were split
	// will return a container of the number of points the children
	// of each child of target might have, in order
	// [0] = left child's left child count, [1] = left child's rght child count,
	// [2] = rght child's left child count, [3] = rght child's rght child count,
	Size_tVec& SPSnode::getChildrensLeftAndRightCountsIfSplit
						(Size_tVec& grandchildCounts) const
	{
		try {

			// first find what the children's boxes would be would be
			int splitMe; // variable to hold first longest dimension
			ivector box = getBox();
			double temp1 = MaxDiam(box, splitMe);

			// ivectors to be new boxes for new children
			ivector rCBox;
			ivector lCBox;
			
			// Call Upper() to get what would be the right hand child box
			Upper(box, rCBox, splitMe);
			// Call Lower() to get what would be the left hand child box
			Lower(box, lCBox, splitMe);
			
			// mid point of my box on first longest dimension
			cxsc::real midSplit = cxsc::mid(box[splitMe]);

			// and if those children were split
			// left Child 
			int splitChildren = MaxDiamComp(lCBox);
			
			cxsc::real midSplitLC = cxsc::mid(lCBox[splitChildren]);
			
			// right child 
			// will split on the same dimension as LC
			
			cxsc::real midSplitRC = cxsc::mid(rCBox[splitChildren]);
			
			// now find how many of this node's data points would go right
			// and left children of left and right children
			size_t rightRightCount = 0;
			size_t rightLeftCount = 0;
			size_t leftRightCount = 0;
			size_t leftLeftCount = 0;
			NodeDataItr it;

			for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
				// DataItrs is a container of iterators to a BigDataCollection
				rvector p = **it;
				// increment left child?
				if ( p[splitMe] < midSplit ) {
					if ( p[splitChildren] < midSplitLC ) leftLeftCount++;
					else rightLeftCount++;
				}
				else { // on right of me
					if ( p[splitChildren] < midSplitRC) leftRightCount++;
					else rightRightCount++;
				}
			}

			grandchildCounts.push_back(leftLeftCount);
			grandchildCounts.push_back(rightLeftCount);
			grandchildCounts.push_back(leftRightCount);
			grandchildCounts.push_back(rightRightCount);


			return grandchildCounts;
		}
		catch (exception& e) {
			string msg = string(e.what());
			throw SPnodeException("Error in getChildrensLeftAndRightCountsIfSplit:\n" + msg);
		}
	}


    // fills in container of leaf counts, left to right
    Size_tVec& SPSnode::getLeafNodeCounts(Size_tVec& counts) const
    {

        if (getLeftChild()!=NULL) {
            getLeftChild()->getLeafNodeCounts(counts);
        }
        if (getRightChild()!=NULL) {
            getRightChild()->getLeafNodeCounts(counts);
        }
        if (getLeftChild()==NULL && getRightChild()==NULL) {

            counts.push_back(counter);
        }
        return counts;
    }


    // return a reference to a container of SPSnodes
    // contents being the leaves descended from this, or this if this is a leaf
    // left to right order
    SPSnodePtrs& SPSnode::getLeaves(SPSnodePtrs& leaves) const
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
            leaves.push_back(const_cast<SPSnode*>(this));
        }
        return leaves;
    }

	  // return a reference to a container of SPSnodes
    // contents being all the nodes in left to right order
    SPSnodePtrs& SPSnode::getAllNodes(SPSnodePtrs& allNodes) const
    {
        if (!isEmpty()) { // this is not empty
		  //if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
				//cout << nodeName << endl;
            allNodes.push_back(const_cast<SPSnode*>(this));
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

    // return a reference to a container of SPSnodes
    // contents being the sub-leaf children of the given node
    // sub-leaf nodes are the parents of leaf nodes and only have leaf nodes
    // as children
    // left to right order
    SPSnodePtrs& SPSnode::getSubLeaves(SPSnodePtrs& subleaves) const
    {
        //if children, recurse on the children
        if (hasLCwithBox()) {
            getLeftChild()->getSubLeaves(subleaves);
        }

        if (hasRCwithBox()) {
            getRightChild()->getSubLeaves(subleaves);
        }

        if (isSubLeaf()) { // this is a subleaf
            // arrgh horrible - cast away const if this node is a subleaf
            subleaves.push_back(const_cast<SPSnode*>(this));
        }
        return subleaves;
    }



    //Returns the sum of the count over volume in the leaf nodes
    real SPSnode::getSumLeafCountOverVol() const
    {
        dotprecision sum(0.0);

        if (isLeaf()) {  // this is a leaf
            accumulate(sum, 1.0*counter, (1.0/nodeVolume()));
        }

        else { // this is not a leaf

            SPSnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            SPSnodePtrsItr it;

            for(it = leaves.begin(); it < leaves.end(); it++) {
                accumulate(sum, 1.0*((*it)->getCounter()),
                            (1.0/(*it)->nodeVolume()));
            }
        }
        return rnd(sum);
    }




    //Returns the count in the smallest (by count) leaf node.
    size_t SPSnode::getSmallestLeafCount() const
    {
        size_t smallestCount = 0;

        if (isLeaf()) {  // this is a leaf
            smallestCount = counter;
        }

        else { // this is not a leaf
            // set up a container for the leaf children
            SPSnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            // find the smallest child by count
            SPSnodePtrsItr it;
            SPSnode* smallest = *(leaves.begin());

            smallestCount = smallest->counter;

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->counter < smallestCount) {

                    smallestCount = (*it)->counter;
                }
            }
        } // end else not a leaf

        return smallestCount;
    }

    // Returns the count in the largest (by count) leaf node.
    size_t SPSnode::getLargestLeafCount() const
    {
        size_t largestCount = 0;

        if (isLeaf()) {  // this is a leaf
            largestCount = counter;
        }

        else { // this is not a leaf

            // set up a container for the leaf children
            SPSnodePtrs leaves;
            // fill the container with the leaf children
            // could be just this if no children
            getLeaves(leaves);

            // find the largest child by volume
            SPSnodePtrsItr it;
            largestCount = (*(leaves.begin()))->counter;

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->counter > largestCount) {
                    largestCount = (*it)->counter;
                }
            }
        } // end else not a leaf

        return largestCount;
    }

    // get the count in the ultimate root node ancestor of this node
    size_t SPSnode::getRootCounter() const
    {
        size_t retValue = 0;
        if (parent == NULL) { // this is root
            retValue = counter;
        }
        else {
            // recurse upwards
            retValue = getParent()->getRootCounter();
        }
        return retValue;
    }

    // Get the mean of the data covered by the box of a node
    rvector SPSnode::getMean() const
    {
       // cout << "Get Mean for " << getNodeName() << endl;
		//	cout << getCountsOnly() << "\t" << getCounter() << endl;
			
        // set up an rvector retMean of the correct dimensions
        rvector retMean(dimension);
        // loop through the elements in the dpSums vector
        for (size_t i = 0; i< dimension; i++) {

            // if no data elements each element or if only counts are held,
            // that element of the mean is 0.0
            if (countsOnly || (counter == 0)) {
                // cxsc::rvector is indexed 1 to n
                retMean[i+1] = 0.0;
            }
            // if data elements, find the element-by-element mean
            else {
                // default cxsc rounding dotprecision rnd_next
                retMean[i+1] = rnd(dpSums[i])/(1.0*counter);
            }
        }// end loop through the elements in dpSums

        return retMean;

    }

	//gat41
	// Get the uniform mean of the box of a node.
	rvector SPSnode::getUniformMean() const
	{
		rvector unifMean(dimension);
		
		// loop through the coordinates of this box to get the midpoint at each
		// coordinate
		ivector thisBox = getBox();
		for (size_t i = 1; i <= dimension; i++) {
			unifMean[i] = mid(thisBox[i]);
		}
		return unifMean;
	}

	//gat41
	real SPSnode::getChebDistMean() const
	{
		rvector Mean = getMean();
		rvector MeanUnif = getUniformMean();
		real ChebDist = 0;
		//loop through the means and get the Chebyshev distances
		for (int i = 1; i <= dimension; i++) {
				real temp = abs(Mean[i] - MeanUnif[i]);
				//std::cout << temp << "\t" << ChebDist << std::endl;
				ChebDist  = ( temp > ChebDist) ? temp : ChebDist;
			}
			
			//cout << "getChebDist: " << endl;
			//std::cout << getNodeName() << "\t" << getMean() << "\t" << getUniformMean() << std::endl;
		return ChebDist;
	}
	
	
	//gat41
	real SPSnode::getChebDistCovar() const
	{
		RealVec Covar = getVarCovar();
		RealVec unifCovar = getUniformVarCovar();

		real ChebDist = 0;
		//loop through the real vector and get the Chebyshev distances
		
		for (int i = 0; i < dimension*dimension; i++) {
				real temp = abs(Covar[i] - unifCovar[i]);
				//std::cout << temp << "\t" << ChebDist << std::endl;
				ChebDist  = ( temp > ChebDist) ? temp : ChebDist;
			}
			
			//cout << "getChebDist: " << endl;
		return ChebDist;
	}
	
	//gat41
	// Get the empirical mass.
	double SPSnode::getEmpMass() const
	{
		int n = getRootCounter();
		double empMass = (getCounter()*1.0)/(1.0*n);
		//cout << nodeName << "\t" << empMass << endl;
		return empMass;
	}
	
   //gat41
   // Get the Battharchya distance.
   real SPSnode::getHellingerDist() const
   {
		RealVec Covar = getVarCovar(); //get the covariance matrix/
		real HD = 0.0; //initialize hellinger distance to 0.

		// if there are no points, cov should be undefined. But since we want to push
		// this node to the bottom of the queue, hence let HD = 0.
		// if there is one point, the variance is 0. At the moment, we do not 
		// want to split boxes with only one point and so also let HD = 0.
		if ( getCounter() == 0 || getCounter() == 1 ) { return HD = 0.0; } 

		else {
	//		cout << "===========================" << getNodeName() << "\t" << getCounter() << endl;
			//get the differences of the mean vectors
			rvector diffMean = getMean() - getUniformMean();
			//cout << "mean differences: " << diffMean << endl;

			//make a gsl matrix for the mean difference
			gsl_matrix * diffMeanMat = gsl_matrix_alloc(dimension, 1);
			for (int i = 0; i < dimension; i++) {
				for (int j = 0; j < 1; j++) {
					gsl_matrix_set (diffMeanMat, i, j, _double(diffMean[i+1]));
				}
			}

			// get the variances
			RealVec	unifCovar = getUniformVarCovar();

			// initialize matrix objects
			gsl_matrix * CovarMat = gsl_matrix_alloc(dimension, dimension);
			gsl_matrix * CovarMatMult = gsl_matrix_alloc(dimension, dimension);
			gsl_matrix * UnifCovarMat = gsl_matrix_alloc(dimension, dimension);
			gsl_matrix * PMat = gsl_matrix_alloc(dimension, dimension); //make this same as
																							//CovarMat first
			int n = getRootCounter();
		
			// problem with stably inverting the covariance matrix - if determinant is wrong, will get -DB
			// fill up the matrics for the var-covar
			int k = 0; //counter for RealVec
			for (int i = 0; i < dimension; i++) {
				for (int j=0; j < dimension; j++) {

					if ( i == j ) {
						gsl_matrix_set(CovarMat, i, j, _double(Covar[k]) + 0.00000001); //cast to double
						gsl_matrix_set(CovarMatMult, i, j, 100*(_double(Covar[k]))+0.00000001);
					}
					else {
						gsl_matrix_set(CovarMat, i, j, _double(Covar[k])); //cast to double
						gsl_matrix_set(CovarMatMult, i, j, 100*(_double(Covar[k])));
					}
					
					gsl_matrix_set(PMat, i, j, _double(Covar[k])); //cast to double
					gsl_matrix_set(UnifCovarMat, i, j, _double(unifCovar[k])); //cast to double
					k++;
				}
			}

			// if variance is -ve, atomic data points? treat as only one point (not
			// a very good assumption at the moment) and let HD = 0. 
			for ( int i = 0; i < dimension; i++) {
				for (int j = 0; j < dimension; j++) {
					if ( (i == j) && (gsl_matrix_get(CovarMat, i, j) < 0) ) {
						cerr << "Negative variance!" << endl;
						cout.precision(20);
						cout << getCounter() << endl;
						NodeDataItr dataItr;
						cout << getNodeName() << endl;
						cout.precision(20);
						cout << "Data is" << std::endl;
						for (dataItr = dataItrs.begin();
							dataItr!= dataItrs.end(); dataItr++) {
							BigDataItr bigIt = *dataItr;
							rvector theData = *bigIt;
							cout << theData << endl; 
						} // end loop through data container
				
						//cerr << "Variance cannot be negative." << endl; 
						//exit(1); 
						
						
						//gsl_matrix_free(CovarMat);
						//gsl_matrix_free(UnifCovarMat);
						//gsl_matrix_free(PMat);
						return HD = 0.0;
						 
					}
				}
			}
	
			//else {
				/*cout << "CovarMat: " << endl;
				for (int i = 0; i < dimension; i++) {
					for (int j=0; j < dimension; j++) {
						cout << i << "\t" << j << "\t" << gsl_matrix_get(CovarMat, i, j) << endl; 
					}
				}
				cout << "CovarMatMult: " << endl;
				for (int i = 0; i < dimension; i++) {
					for (int j=0; j < dimension; j++) {
						cout << i << "\t" << j << "\t" << gsl_matrix_get(CovarMatMult, i, j) << endl; 
					}
				}
				cout << "UnifCovarMat: " << endl;
				for (int i = 0; i < dimension; i++) {
					for (int j=0; j < dimension; j++) {
						cout << i << "\t" << j << "\t" << gsl_matrix_get(UnifCovarMat, i, j) << endl; 
					}
				}		
					*/
				//add the two matrices
				gsl_matrix_add(PMat, UnifCovarMat);
				gsl_matrix_scale(PMat, 0.5);
				
				gsl_matrix * PMatForInv = gsl_matrix_alloc(dimension, dimension);
				PMatForInv = PMat; 
				/*cout << "PMat: " << endl;
				for (int i = 0; i < dimension; i++) {
					for (int j=0; j < dimension; j++) {
						cout << i << "\t" << j << "\t" << gsl_matrix_get(PMat, i, j) << endl; 
					}
				}*/

				// get the determinants of CovarMat, UnifCovarMat, PMat
				int s;
				gsl_permutation * p = gsl_permutation_alloc(dimension);
				gsl_linalg_LU_decomp(CovarMatMult, p, &s);
				//cout << "CovarMat LU decomp: " << endl;
				for (int i = 0; i < dimension; i++) {
					for (int j=0; j < dimension; j++) {
						//cout << i << "\t" << j << "\t" << gsl_matrix_get(CovarMat, i, j) << endl; 
					}
				}
				double detCovarMat = gsl_linalg_LU_det(CovarMatMult, s)/(pow(100,dimension));
				//cout << "det covar mat: " << detCovarMat << "\t" << endl;
				gsl_permutation_free(p);
				// it is possible to get negative determinants, and negative DB, temporarily sweeping this under the rug 
				if (detCovarMat < 0) { cout << "Negative determinant." << endl; exit(1); }//return HD = 0.0; }
				else {
					
					int s1;
					gsl_permutation * p1 = gsl_permutation_alloc(dimension);
					gsl_linalg_LU_decomp(UnifCovarMat, p1, &s1);
					//cout << "UnifCovarMat decomp: " << endl;
					for (int i = 0; i < dimension; i++) {
						for (int j=0; j < dimension; j++) {
							//cout << i << "\t" << j << "\t" << gsl_matrix_get(UnifCovarMat, i, j) << endl; 
						}
					}
					double detUnifCovarMat = gsl_linalg_LU_det(UnifCovarMat, s1);
					//cout << "det unif covar mat: " << detUnifCovarMat << endl;
					gsl_permutation_free(p1);

					int s2;
					gsl_permutation * p2 = gsl_permutation_alloc(dimension);
					gsl_linalg_LU_decomp(PMat, p2, &s2);
					double detPMat = gsl_linalg_LU_det(PMat, s2);
					//cout << "det p mat: " << detPMat << endl;
	
					//now get the inverse of P
					gsl_matrix * Pinverse = gsl_matrix_alloc(dimension, dimension); 
					gsl_linalg_LU_invert(PMat, p2, Pinverse);
					gsl_permutation_free(p2);
					//cout << "p inverse " << endl;
					for (int i = 0; i < dimension; i++) {
						for (int j=0; j < dimension; j++) {
							//cout << i << "\t" << j << "\t" << gsl_matrix_get(Pinverse, i, j) << endl; 
						}
					}
			
					//now get the Bhattacharya coefficient
					//DB = 1.0/8.0 * diffMeanMat1 * invP * diffMeanMat2 + 0.5*log(detP/sqrt{detP1}{detP2});
					gsl_matrix * diffMeanMatTranspose = gsl_matrix_alloc(1, dimension);
					gsl_matrix_transpose_memcpy(diffMeanMatTranspose, diffMeanMat);
				
					// now convert all to rmatrix for easier computations
					rmatrix diffMeanR(0, dimension-1, 0, 0);
					for (int i = 0; i < dimension; i++) {
						for (int j=0; j < 1; j++) {
							diffMeanR[i][j] = gsl_matrix_get(diffMeanMat, i, j); 
						}
					}
				
					rmatrix diffMeanTransR(0, 0, 0, dimension-1);
					for (int i = 0; i < 1; i++) {
						for (int j=0; j < dimension; j++) {
							diffMeanTransR[i][j] = gsl_matrix_get(diffMeanMatTranspose, i, j); 
						}
					}
				
					rmatrix PinvR(0, dimension-1, 0, dimension-1);
					for (int i = 0; i < dimension; i++) {
						for (int j=0; j < dimension; j++) {
							PinvR[i][j] = gsl_matrix_get(Pinverse, i, j); 
						}
					}
				
					//free the gsl_matrices
					gsl_matrix_free(CovarMat);
					gsl_matrix_free(UnifCovarMat);
					gsl_matrix_free(PMat);
				
					//cout << diffMeanR << endl;
					//cout << diffMeanTransR << endl;
					//cout << PinvR << endl;
				
					//now get the Bhattacharya coefficient
					//DB = 1.0/8.0 * diffMeanMat1 * invP * diffMeanMat2 + 0.5*log(detP/sqrt{detP1}{detP2});
					//cout << (diffMeanTransR*PinvR)*diffMeanR << endl;
					rmatrix MatOp = (diffMeanTransR*PinvR)*diffMeanR;
					//cout << MatOp << "\t" << MatOp[0][0] << endl;
					assert(MatOp >= 0);
					real DB = 1.0/8.0 * MatOp[0][0] + 0.5*log(detPMat/sqrt(detCovarMat*detUnifCovarMat));
					//cout << "DB: " << DB << endl;
					if (DB < 0) { return HD = 0.0; }
					else {
						real BC = exp(-DB);
						assert(BC >= 0);
					//cout << "BC: " << BC << endl;
						real HD = sqrt(1-BC);
					//cout << "HD: " << HD << endl;
						assert(HD >=0);
						return HD;
					}
				} // end of determinant is not zero
			}
	}   

	//gat41
   // Get the Hellingr distance for 1D data.
   real SPSnode::getHellingerDist1D() const
   {
		real HD = 0.0; 
		RealVec Covar = getVarCovar();

		// if the variance is negative - need to investiage this more
		
		if (Covar[0] < 0) {
			cout << getCounter() << endl;
			cout << Covar[0] << endl;
			NodeDataItr dataItr;
			cout << getNodeName() << endl;
			cout.precision(20);
			cout << "Data is" << std::endl;
          for (dataItr = dataItrs.begin();
                dataItr!= dataItrs.end(); dataItr++) {

                BigDataItr bigIt = *dataItr;
                rvector theData = *bigIt;

                cout << theData << endl; 
			} // end loop through data container
			
			cerr << "Variance cannot be negative." << endl; 
		//	exit(1);
			return HD = 0.0; 
		}

		// can continue if variance is not negative

			// if there are no points, should be undefined. But since we want to push
			// this node to the bottom of the queue, hence let HD = 0.
			// if there is one point, the variance is 0. At the moment, we do not 
			//want to split boxes with only one point and so also let HD = 0.
			// if variance is -ve, atomic data points? treat as only one point (not
			// a very good assumption at the moment) and let HD = 0. 
			if ( getCounter() == 0  || Covar[0] <= 0 ) { return HD = 0.0; } 
			/*else if ( Covar[0] == 0 ) { 
				cout << getCounter() << endl;
				cout << getEmpMass() << endl;
				cout << nodeVolume() << endl;
				cout << getMean() << endl;
				//cerr << "no variance. check!" << endl;
				//exit(1);
			}*/

			else {
				//get the differences of the mean vectors
				rvector diffMean = getMean() - getUniformMean();
				//cout << "mean differences: " << diffMean[1] << endl;
				
				// get the variances
				RealVec	unifCovar = getUniformVarCovar();
				//cout << "Covar: " << Covar[0] <<  endl;
				//cout << "unifCovar: " << unifCovar[0] << endl;
				//if all the elements for CovarMat are all zero, we do not have any points in 
				//this leaf node - so return hellinger distance as 0
	
				// use the sqrt of the squared hellinger distance for two normal distributions
				// 1 - sqrt(2*sigma1*sigma2/(sigma1^2 + sigma2^2))*exp(-0.25*(mu1-mu2)^2/(sigma1^2+sigma2^2))
				//cout << diffMean << "\t";
				interval covarI = interval(Covar[0]);
				interval unifCovarI = interval(unifCovar[0]);
				interval sumVar = covarI + unifCovarI;
				//cout << sumVar << "\t";
				interval insqrt = 2*sqrt(covarI)*sqrt(unifCovarI)/sumVar;
				//cout << insqrt << "\t";
				interval H2 = interval(1,1) - sqrt(insqrt) *exp((-0.25*diffMean[1]*diffMean[1])/sumVar);
				//cout << "H2: " << H2 << endl;
				HD = mid(sqrt(H2));
				//cout << HD << endl;
				if ( HD > 1 || HD < 0) { 
					cerr << "HD should be between 0 and 1." << endl;
					exit(0);
				}
				return HD;
			}

	}

   
    // Get the variance-covariance vector of the data covered
    // by the box of a node
    RealVec& SPSnode::getVarCovar(RealVec& varCovar) const
    {
        varCovar.clear();
        varCovar.reserve(dimension*dimension);

        // loop through the elements in the dpSumProducts vector
        for (size_t k = 0; k < dimension*dimension; k++) {

            // counts only held or if 0 or 1 data points
            // each element of the var-covar is 0.0
            if (countsOnly || (counter <= 1)) {
                varCovar.push_back(0.0);
            }
            // if >1 data points find element-by-element var-covar

            /*the var-covar is the sample var-covar
            which is
            [sumproduct(i,j)-sum(i)sum(j)/counter]/(counter-1)

            element k in the vector of dotprecison sumproducts
            corresponds to row k/n, (row 0 to n-1)
            and column k-row*n (col 0 to n-1)
            in a matrix view of the sumproducts */

            else {
                size_t i = k/dimension; // row  (int/int = int)
                size_t j = k - i*dimension; // column

                // make another dotprecision variable
                dotprecision temp1 = dpSumProducts[k];

                dotprecision temp2(0.0);
                // sum(i) x sum(j)
                // default cxsc rounding dotprecision rnd_next
                accumulate(temp2,  rnd(dpSums[i]),
                        rnd(dpSums[j]));

                real div = -1.0/counter;

                // sumproduct(i,j) - sum(i)(sum(j)/counter
                // default cxsc rounding
                accumulate(temp1, rnd(temp2), div);
                // calculate the variance covariance element
                varCovar.push_back(rnd(temp1)/(1.0*(counter-1)));
            }
        }// end loop through the elements in dpSumProducts

        return varCovar;
		}
	
	    RealVec SPSnode::getVarCovar() const
		{
        RealVec retVarCovar;
        retVarCovar = getVarCovar(retVarCovar);
        return retVarCovar;
		}

	RealVec SPSnode::getUniformVarCovar() const
    {
        RealVec retVarCovar;
        retVarCovar = getUniformVarCovar(retVarCovar);
        return retVarCovar;
    }
    
	 // Get the uniform variance-covar./R	F	iance vector of the data covered
    // by the box of a node
    RealVec& SPSnode::getUniformVarCovar(RealVec& unifVarCovar) const
    {
			unifVarCovar.reserve(dimension*dimension);
			ivector thisBox = getBox();
			
			// fill in the matrix where the diag are (1/12)*(b-a)^2 and off-diag 
			// are 0.
			for (size_t i = 0; i < dimension*dimension; i++) {
					unifVarCovar.push_back(0.0); //first fill up the container with 0
			}
			// then fill up the diags
			for (size_t i = 0; i < dimension; i++) {
				int pos = i*dimension + i;
				unifVarCovar[pos] = 1.0/12.0 * (Sup(thisBox[i+1]) - Inf(thisBox[i+1])) 
											* (Sup(thisBox[i+1]) - Inf(thisBox[i+1]));
			}
			
			//for (size_t i = 0; i < dimension*dimension; i++) {cout << unifVarCovar[i] << endl;}

        return unifVarCovar;
		}


		

    // Print the details of a of a specific node in a subpaving
    std::ostream& SPSnode::nodePrint(std::ostream &os) const
    {
        // output for box in form:
        // box, volume, counter, mean, variance covariance, and data

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            os << "Box is :";

            for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
                // c-xsc default output for intervals
                os << "  " << thisBox[i];
            }

            os << std::endl;
            os << "Box volume is " << nodeVolume() << std::endl;
            os << "Counter = " << counter << std::endl;

            nodeMeanPrint(os);
            nodeVarCovarPrint(os);
            nodeDataPrint(os);


            os << std::endl;
        }
        return os;

    }


    // Get this node's contribution to loglikelihood
    real SPSnode::getLogLik(const size_t n) const
    {
        // contribution to loglikelihood is counter*log(counter/(n * vol))

        dotprecision contribution(0.0);
        if ((n > 0) && (counter > 0)) {

            accumulate(contribution, 1.0*counter, log(1.0*counter));
            accumulate(contribution, -1.0*counter, log(1.0*n));
            accumulate(contribution, -1.0*counter, log(nodeVolume()));
        }
        // default cxsc rounding to nearest
        return rnd(contribution);
    }

    // Get change in scaled contribution to log likelihood split
    dotprecision SPSnode::getSplitChangeLogLik() const
    {
        
       // cout << "caling get split change llk" << endl;
        
        dotprecision change;
        change = 0.0;

        // if counter is 0 there can be no change on splitting
        if (counter > 0) {

            // first find what the left hand child's counter would be if
            // that child were to be created
            size_t leftCount = getLeftCountIfSplit();
            // current number of data points associated to node is counter
            size_t rightCount = counter-leftCount;
		
            // current node volume from nodeVolume; each child will have half
            // change is
            //      (lc_count*ln(lc_count) + rc_count*ln(rc_count) +
            //          count*ln(2)) - (count*ln(count)
            // if we split and lc_count, rc_count were the new counts in
            // left and right children respectively
            // note that the terms involving the total count in the histogram
            // and the volume of this node cancel so this change
            // is effectively scaled and does not need to use n

            dotprecision currentEMP(0.0);
            dotprecision childEMP(0.0);

            if (leftCount > 0) accumulate(childEMP, 1.0*leftCount,
                             log(1.0*leftCount));

            if (rightCount > 0) accumulate(childEMP, 1.0*rightCount,
                             log(1.0*rightCount));

            accumulate(childEMP, 1.0*counter, log(2.0));

            accumulate(currentEMP, 1.0*counter, log(1.0*counter));

            change = childEMP - currentEMP;
        }
        return change;
    }

    // Get change in contribution to log likelihood on merge of leaf children
    dotprecision SPSnode::getMergeChangeLogLik() const
    {
        dotprecision change;
        change = 0.0;

        // if counter is 0 there can be no change on merging
        if (counter > 0) {

            // first find what the left hand child's counter is
            size_t leftCount = getLeftChild()->getCounter();

            // and right child
            size_t rightCount = getRightChild()->getCounter();

            // change is (count*ln(count)
            //      - (lc_count*ln(lc_count) + rc_count*ln(rc_count) +
            //          count*ln(2))
            // note that the terms involving the total count in the histogram
            // and the volume of this node cancel so this change
            // is effectively scaled and does not need to use n

            dotprecision currentEMP(0.0);
            dotprecision childEMP(0.0);
            if (leftCount > 0) accumulate(childEMP, 1.0*leftCount,
                             log(1.0*leftCount));

            if (rightCount > 0) accumulate(childEMP, 1.0*rightCount,
                             log(1.0*rightCount));

            accumulate(childEMP, 1.0*counter, log(2.0));

            accumulate(currentEMP, 1.0*counter, log(1.0*counter));

            change = currentEMP - childEMP;
        }

        return change;
    }



    // Returns the best (smallest positive or most negative) change in EMP
    // from splitting any leaf node under COPERR
    // n is total points, used for scaling
    dotprecision SPSnode::getBestSplitChangeEMPCOPERR(const size_t n) const
    {
        dotprecision bestEMPChange;
        bestEMPChange = 0.0;

        if (isLeaf()) {  // this is a leaf
            bestEMPChange = getSplitChangeEMPCOPERR(n);
        }

        else { // this is not a leaf
            // set up a container for the leaf children
            SPSnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            // find the best child for splitting
            SPSnodePtrsItr it;
            SPSnode* best = *(leaves.begin());

            bestEMPChange = best->getSplitChangeEMPCOPERR(n);

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->getSplitChangeEMPCOPERR(n) < bestEMPChange) {

                    bestEMPChange = (*it)->getSplitChangeEMPCOPERR(n);
                }
            }
        } // end else not a leaf

        return bestEMPChange;
    }

    // Returns the best (smallest positive or most negative) change in EMP
    // from splitting any leaf node under AIC
    // term involving n cancels out of change
    dotprecision SPSnode::getBestSplitChangeEMPAIC() const
    {
        dotprecision bestEMPChange;
        bestEMPChange = 0.0;

        if (isLeaf()) {  // this is a leaf
            bestEMPChange = getSplitChangeEMPAIC();
        }

        else { // this is not a leaf
            // set up a container for the leaf children
            SPSnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            // find the best child for splitting
            SPSnodePtrsItr it;
            SPSnode* best = *(leaves.begin());

            bestEMPChange = best->getSplitChangeEMPAIC();

            for(it = leaves.begin(); it < leaves.end(); it++) {
                if ((*it)->getSplitChangeEMPAIC() < bestEMPChange) {

                    bestEMPChange = (*it)->getSplitChangeEMPAIC();
                }
            }
        } // end else not a leaf

        return bestEMPChange;
    }

// Returns the best (smallest positive or most negative) change in EMP
    // from merging any subleaf node under COPERR
    // n is total points, used for scaling
    dotprecision SPSnode::getBestMergeChangeEMPCOPERR(const size_t n) const
    {
        dotprecision bestEMPChange;
        bestEMPChange = 0.0;

        if (isSubLeaf()) {  // this is a subleaf
            bestEMPChange = getMergeChangeEMPCOPERR(n);
        }

        else { // this is not a subleaf
            // set up a container for the subleaf children
            SPSnodePtrs subleaves;
            // fill the container with the leaf children
            getSubLeaves(subleaves);

            // find the best child for splitting
            SPSnodePtrsItr it;
            SPSnode* best = *(subleaves.begin());

            bestEMPChange = best->getMergeChangeEMPCOPERR(n);

            for(it = subleaves.begin(); it < subleaves.end(); it++) {
                if ((*it)->getMergeChangeEMPCOPERR(n) < bestEMPChange) {

                    bestEMPChange = (*it)->getMergeChangeEMPCOPERR(n);
                }
            }
        } // end else not a subleaf

        return bestEMPChange;
    }

    // Returns the best (smallest positive or most negative) change in EMP
    // from merging any subleaf node under AIC
    // term involving n cancels out of change
    dotprecision SPSnode::getBestMergeChangeEMPAIC() const
    {
        dotprecision bestEMPChange;
        bestEMPChange = 0.0;

        if (isSubLeaf()) {  // this is a subleaf
            bestEMPChange = getMergeChangeEMPAIC();
        }

        else { // this is not a sub leaf
            // set up a container for the subleaf children
            SPSnodePtrs subleaves;
            // fill the container with the leaf children
            getSubLeaves(subleaves);

            // find the best child for splitting
            SPSnodePtrsItr it;
            SPSnode* best = *(subleaves.begin());

            bestEMPChange = best->getMergeChangeEMPAIC();

            for(it = subleaves.begin(); it < subleaves.end(); it++) {
                if ((*it)->getMergeChangeEMPAIC() < bestEMPChange) {

                    bestEMPChange = (*it)->getMergeChangeEMPAIC();
                }
            }
        } // end else not a subleaf

        return bestEMPChange;
    }

    // Get this node's scaled contribution to sum term
    real SPSnode::getEMPContributionCOPERR(const size_t n) const
    {
        // current number of data points associated to node is counter
        // current node volume from nodeVolume, and each child will have half

        dotprecision contribution(0.0);
        if ((n > 0) && (counter > 0)) {
            accumulate(contribution, -(1.0*counter)/(1.0*n),
                    (1.0*counter)/(n*nodeVolume()));
        }

        // contribution is -counter^2/(n^2 * vol)
        // default cxsc rounding to nearest

        //return contribution;
        return rnd(contribution);
    }

    // Get this node's scaled contribution to EMP under AIC.
    real SPSnode::getEMPContributionAIC(const size_t n) const
    {
        return -getLogLik(n);
    }

    // Get change in scaled contribution to EMP under COPERR on split
    dotprecision SPSnode::getSplitChangeEMPCOPERR(const size_t n) const
    {
        // first find what the left hand child's counter would be if that child
        // were to be created
        size_t leftCount = getLeftCountIfSplit();

        // current number of data points associated to node is counter
        // current node volume from nodeVolume, and each child will have half

        // change is 1/(n^2 * vol) x (counter^2 - 2(lc_count^2 + rc_count^2))
        // if we split and lc_count, rc_count were the new counts in
        // left and right children respectively
        // Change is scaled by n, total points in histogram
        dotprecision change;
        change = 0.0;
        if (n > 0) {
            accumulate(change, (1.0*counter)/(n*nodeVolume()),
                                (1.0*counter)/(1.0*n));
            accumulate(change, (1.0*leftCount)/(1.0*n),
                                -(2.0*leftCount)/(n*nodeVolume()));
            accumulate(change, (1.0*(counter - leftCount))/(1.0*n),
                                -(2.0*(counter - leftCount))/(n*nodeVolume()));
        }
        return change;
    }

    // Get change in sum term in EMP under AIC on split.
    dotprecision SPSnode::getSplitChangeEMPAIC() const
    {
        return -getSplitChangeLogLik();

    }

    // Get change in scaled contribution to EMP under COPERR on merge
    dotprecision SPSnode::getMergeChangeEMPCOPERR(const size_t n) const
    {
        // first find what the left hand child's counter is
        size_t leftCount = getLeftChild()->getCounter();

        // and right child
        size_t rightCount = getRightChild()->getCounter();

        // current number of data points associated to node is counter
        // current node volume from nodeVolume, and each child will have half

        // change is 1/(n^2 * vol) x (2(lc_count^2 + rc_count^2) - counter^2)
        // Change is scaled by n, total points in histogram
        dotprecision change;
        change = 0.0;
        if (n > 0) {
            accumulate(change, (1.0*leftCount)/(1.0*n),
                                (2.0*leftCount)/(n*nodeVolume()));
            accumulate(change, (1.0*(counter - leftCount))/(1.0*n),
                                (2.0*(counter - leftCount))/(n*nodeVolume()));
            accumulate(change, -(1.0*counter)/(n*nodeVolume()),
                                    (1.0*counter)/(1.0*n));
        }
        return change;
    }



    // Get change in sum term in EMP under AIC on merge.
    dotprecision SPSnode::getMergeChangeEMPAIC() const
    {
        return -getMergeChangeLogLik();
    }

	//src_trunk_0701
	void SPSnode::reshapeToUnion(const SPnode& other)
{
	SPnode::reshapeToUnion(other);
	
}

void SPSnode::reshapeToUnion(const SPnode& other,
						size_t minChildPoints)
{
	if (isEmpty() || other.isEmpty()) {
		throw NoBox_Error(
		"SPnode::reshapeToUnion(const SPnode&, size_t)");
	}
	if ( getBox() != other.getBox() )  {
		throw IncompatibleDimensions_Error(
		"SPnode::reshapeToUnion(const SPnode&, size_t)");
	}
	if ( !other.checkTreeStateLegal() )
		throw runtime_error(
		"SPnode::reshapeToUnion(const SPnode&, size_t) : other has illegal tree state");
	
	std::string baseErrorFilename("ReshapeErrors");
	std::string errorFilename = getUniqueFilename(baseErrorFilename);
	
	bool success = this->_reshapeToUnion(
					&other, minChildPoints, errorFilename);
	
	// if we returned success there should be no file with that name
	if(!success) {
		std::cerr << "\nCould not exactly reshape this to the union:"
			<< " check " << errorFilename << " for errors\n" << endl;
	}
	
}
//--src_trunk_0701

    //Output for all the  leaf boxes in this, using tab delimiters
    std::ostream& SPSnode::leavesOutputTabs(std::ostream &os) const
    {
        // uses  member function leafOutputTabs to generate node output
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
    // including unscaled EMP contributions and changes if split
    std::ostream& SPSnode::leavesOutputTabsWithEMPs(const size_t bigN,
                            std::ostream &os, const int prec) const
    {
        if (parent == NULL) { // root
            std::string headers = "node \t vol \t count \t EMP COPERR ";
            headers += "\t &change \t EMP AIC \t &change \t dimensions \n";
            os << headers;
        }

        // uses  member function leafOutputTabsWithEMPs to generate node output
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafOutputTabsWithEMPs(bigN, os, prec);
            return (os << "\n");

        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesOutputTabsWithEMPs(bigN, os, prec);
        }

        if (getRightChild()!=NULL) {
            getRightChild()->leavesOutputTabsWithEMPs(bigN, os, prec);
        }

    }


    //Output for all the  leaf boxes in this, using tab delimiters
    //including output for the height of histogram bins for a
    // normalised histogram based on this tree with this as root
    std::ostream& SPSnode::leavesOutputTabsWithHistHeight(
                        std::ostream &os, int prec) const
    {
        leavesOutputTabsWithHistHeight(counter, os, prec);
        return (os);
    }

    //Output for all the  leaf boxes in this, using tab delimiters
    //including output for the height of histogram bins for a
    // normalised histogram based on tree with total number of data points bigN.
    std::ostream& SPSnode::leavesOutputTabsWithHistHeight(const size_t bigN,
                        std::ostream &os, const int prec) const
    {
        // uses  member function leafOutputTabs to generate node output
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafOutputTabsWithHistHeight(bigN, os, prec);
            return (os << "\n");

        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesOutputTabsWithHistHeight(bigN, os, prec);
        }

        if (getRightChild()!=NULL) {
            getRightChild()->leavesOutputTabsWithHistHeight(bigN, os, prec);
        }
    }

    //Output for all the  leaf boxes in this, using tab delimiters
    // including unscaled EMP contributions and changes if split
    std::ostream& SPSnode::leavesOutputTabsWithHistHeightAndEMPs(
                    const size_t bigN, std::ostream &os, const int prec) const
    {
        if (parent == NULL) { // root
            std::string headers = "node \t vol \t count \t height ";
            headers += "\t EMP COPERR \t &change \t EMP AIC \t &change ";
            headers += "\t dimensions \n";
            os << headers;
        }

        // uses  member function leafOutputTabsWithEMPs to generate node output
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            leafOutputTabsWithHistHeightAndEMPs(bigN, os, prec);
            return (os << "\n");

        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            getLeftChild()->leavesOutputTabsWithHistHeightAndEMPs(bigN,
                                                                    os, prec);
        }

        if (getRightChild()!=NULL) {
            getRightChild()->leavesOutputTabsWithHistHeightAndEMPs(bigN,
                                                                    os, prec);
        }

    }

    // Get the scaled EMP sum term under COPERR for the tree rooted at this
    dotprecision SPSnode::getEMPSumCOPERR(const size_t n) const
    {
        dotprecision retValue;
        retValue = 0.0;

        // uses  member function getEMPContributionCOPERR for leaf value
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            retValue = getEMPContributionCOPERR(n);
        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            retValue = retValue + getLeftChild()->getEMPSumCOPERR(n);
        }

        if (getRightChild()!=NULL) {
            retValue = retValue + getRightChild()->getEMPSumCOPERR(n);
        }

        return retValue;

    }

    // Get the unscaled EMP sum term under AIC for the tree rooted at this
    dotprecision SPSnode::getEMPSumAIC(const size_t n) const
    {
        dotprecision retValue;
        retValue = 0.0;

        // uses  member function getEMPContributionAIC for leaf result
        if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
            retValue = getEMPContributionAIC(n);
        }

        //recurse on the children
        if (getLeftChild()!=NULL) {
            retValue+=getLeftChild()->getEMPSumAIC(n);
        }

        if (getRightChild()!=NULL) {
            retValue+=getRightChild()->getEMPSumAIC(n);
        }
        return retValue;
    }


    // Check if a node contains a datapoint
	// it is assumed that the node will have a box
	// childInd is an indicator for which child is being checked
	bool SPSnode::nodeContains(const rvector& p,
							OPERATIONS_ON childInd) const
	{
		bool retValue = false;
		
		// only check for total containment if this is a parent node
		// or to be treated as such
		if ((parent == NULL) || (childInd == ON_PARENT)) {
		
			// cast p to an ivector
			ivector pvec = _ivector(p);
		
			//find if p is in the box
			if (pvec <= getBox()) {
				retValue = true;
			}
		}
		// if not to be treated as a whole box, we assume it was in the parent
		// and only check ourselves with respect to the split dimension
		// and if this is a right child it can be anywhere 
		// but  this is a left child,
		// we need to check the split
		// find what dimension the parent was split on and what
		// the split value was
		// pvector must be strictly less than parentSplitValue
		// on the split dimension
		else if (parent != NULL) { // truly not a parent node
		
			int parentSplitDim = (getParent())->getSplitDim();
			real parentSplitValue = (getParent())->getSplitValue();
			
			bool strictlyLeft = (p[parentSplitDim] < parentSplitValue);

			if (childInd == ON_LEFT) { 
				retValue = strictlyLeft;
				}
			if (childInd == ON_RIGHT) { 
				retValue = !strictlyLeft; 
			}
		}
		
		return retValue;
	}

   // Get the number of counts in any given box
	int SPSnode::spsContains(ivector & z, int countBox, int countInBox) const
   {
      // Query box is z and box to be compared is theBox accessed by getBox()
      // z is assumed not to be empty
      // nb Intersection() gives error if unequal index sets passed
         
      BOOL_INTERVAL retValue = BI_FALSE; // for the return value
        
      // case of no intersection, i.e. the intersection of theBox and z is
      // empty 
        ivector r; // temporary ivector
         if (Intersection(r, z, getBox())==0){
            retValue = BI_FALSE;
            return countBox; 
         }        
  
         // case of isEmpty() being true
         if (isEmpty()){
         BOOL_INTERVAL retValue = BI_FALSE; // for the return value
         }
 
         // case of a non-empty leaf
         if (!isEmpty() && isLeaf()) {           
            ivector r; // temporary, to be passed to Intersection                 
            if (z<getBox()) {
            // z is contained inside theBox but borders are not included
             
               // make z a SPSnode
               SPSnode queryNode(z, 1);
  
               // access data in theBox
               NodeData boxData = getData();
               NodeDataItr boxDataItr;
                 
               // check if queryNode contains data from theBox
               for(boxDataItr = boxData.begin(); boxDataItr != boxData.end(); boxDataItr++){
                  BigDataItr bigIt = *boxDataItr; 
                  rvector theData = *bigIt;  // convert NodeData to rvector
                    
                  // if nodeContains is true, increment countInBox
                  if(queryNode.nodeContains(theData, ON_PARENT)){
                     countInBox += 1;
                  }
               }                             
                 
                // final count;
                countBox = countInBox;
                retValue = BI_TRUE;
            } 
          
            else if (z==getBox()) {
            // z fully covers theBox, including the borders
               countBox += getCounter();
               retValue = BI_TRUE;
            }
       
            else if (Intersection(r, z, getBox())) {
            // result is indeterminate if there is an intersection but z is not             
				// wholly in theBox"
                
               if (Inf(r) == Sup(r)){
                  retValue = BI_INDET;
               }
 
               else{
						// make z a SPSnode
						SPSnode queryNode(z, 1);
  
						// access data in theBox
						NodeData boxData = getData();
						NodeDataItr boxDataItr;
                 
						// check if queryNode contains data from theBox
						for(boxDataItr = boxData.begin(); boxDataItr != boxData.end(); boxDataItr++){
							BigDataItr bigIt = *boxDataItr; 
							rvector theData = *bigIt;  // convert NodeData to rvector
                     
							// if nodeContains is true, increment countInBox
							if(queryNode.nodeContains(theData, ON_PARENT)){
								countInBox += 1;
							}
						}                             
                
						// final count;
						countBox += countInBox;
						retValue = BI_INDET;
               }
            }
 
            else {
					// Case that there is no intersection
					retValue = BI_FALSE;
            }
            
           return countBox;
         } // end (!isEmpty() && isLeaf())
          
         //case of an non-empty non-leaf
         if (!isEmpty() && !isLeaf()) {
            ivector Lz, Rz; // ivectors passed to Intersection()
            // will contain intersection after Intersection() call
 
            // to hold results of tests on left and right children
            BOOL_INTERVAL Ltest = BI_FALSE;
            BOOL_INTERVAL Rtest = BI_FALSE;
 
            // indicators for tested left and right sides
            bool LtestSuccess = false;
            bool RtestSuccess = false;

 //        cout << "    // Find if there is a leftChild with a box" << endl;
           if (hasLCwithBox() &&
               Intersection(Lz, z, getLeftChild()->getBox())) {
               // Lz contains intersctn of z & leftChild box
               // test Lz and left child node
               //Ltest = 
               countBox =  (getLeftChild()->spsContains(Lz, countBox, countInBox));
               LtestSuccess = true;
            }
                  
   //    cout << "   // Find if there is a rightChild with a box" << endl;
           if (hasRCwithBox() &&
               Intersection(Rz, z, getRightChild()->getBox())) {
               // Rz contains intersctn of z & rightChild box 
               // test Rz and right child node
               //Rtest = 
               countBox =  getRightChild()->spsContains(Rz, countBox, countInBox);
               RtestSuccess = true;
            }
        
 //          cout << "  // if both children tested " << endl;
           if (LtestSuccess && RtestSuccess) {
               //return value is the result of both tests
               // if the same or BI_INDET if diff
               Ltest==Rtest ?
               retValue = Ltest : retValue=BI_INDET;
            }
          
 //           cout << " // if has two children but neither was tested" << endl;
            // ie neither Intersection() returned true
            if (hasRCwithBox() && hasLCwithBox()
               && !LtestSuccess && !RtestSuccess) {
               retValue = BI_FALSE;
               // note that the AIA book has BI_TRUE here
               // but this can't be correct
            }
  
   //       cout << "  // if has two children but only right was tested" << endl;
             // ie left Intersection() returned false
             if (hasRCwithBox() && hasLCwithBox()
                && !LtestSuccess && RtestSuccess) {
                retValue = Rtest;
                // return value result of test of right side
             }
            
     //   cout << "    // if has two children but only left was tested" << endl;
             // ie right Intersection() returned false
             if (hasRCwithBox() && hasLCwithBox()
                 && LtestSuccess && !RtestSuccess) {
                 retValue = Ltest;
                 // return value result of test of left side
             }
           
 //    cout << "   // if has right child only and that child was tested" << endl;
             // ie Intersection() returned true
             if (hasRCwithBox()
                 && !hasLCwithBox() && RtestSuccess) {
                 // if all of z contained in right child's box
                 if (Rz==z) {
                     retValue = Rtest;
                 }
                 // return false if Rtest false, else INDET
                 else {
                       Rtest==BI_FALSE ? retValue = BI_FALSE
                       : retValue = BI_INDET;
                 }
            }
            
 //    cout << "  // if has right child only and that child not tested" << endl;
             // ie Intersection() returned false
             if (hasRCwithBox()
                 && !hasLCwithBox() && !RtestSuccess) {
                 retValue = BI_FALSE;
             }
             
 //    cout << "   // if has left child only and that child was tested" << endl;
             // ie Intersection() returned true
             if (!hasRCwithBox() && hasLCwithBox()
                 && LtestSuccess) {
                 // if whole of z contained in left child's box
                 if (Lz==z) {
                     retValue = Ltest;
                 }
                 // return false if Ltest false, otherwise INDET
                 else {
                       Ltest==BI_FALSE ?
                       retValue = BI_FALSE :
                       retValue = BI_INDET;
                  }
             }
           
 //   cout << "   // if has left child only & that child was not tested" << endl;
             // ie Intersection() returned false
             if (!hasRCwithBox() && hasLCwithBox()
                 && !LtestSuccess) {
                 retValue = BI_FALSE;
             }
 
           //  case no children covered by isLeaf() block above
         } // end of (!isEmpty() && !isLeaf())
 
     return countBox;
       
   } // end of spsContains for ivector


    // Expand and split data
    // should be called on a node with no associated data
    void SPSnode::nodeExpand(int comp)
    {
        nodeExpansionOnly(comp);    // expand the node
        SplitNever sn;              // dummy split decision maker
        splitData(sn);            // split the data with no further splitting


    }

    // Expand and split the data with further splitting
    // should be called on a node with no associated data
    void SPSnode::nodeExpand(const SplitDecisionObj& boolTest, int comp)
    {
        nodeExpansionOnly(comp);    // expand the node
        // split the data, allowing for further splitting

        splitData(boolTest);
    }


    // Expand and split
    // finds its own comp argument
    void SPSnode::nodeExpand()
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(maxdiamcomp); // complete nodeExpand

    }

    // Expand and split with further splitting
    // finds its own comp argument
    void SPSnode::nodeExpand(const SplitDecisionObj& boolTest)
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(boolTest, maxdiamcomp); // complete nodeExpand
    }

    // reabsorb both the children of this node
    // the stats in this node will be correct so all we have to do is to
    // associate the child data with this and delete the children
    // This would work even if children are not leaves since gatherData()
    // gets data from descendents of the node it operates on
    void SPSnode::nodeReabsorbChildren()
    {
        // first recursively deal with the children of the children
        if (hasLCwithBox())
            getLeftChild()->nodeReabsorbChildren();
        if (hasRCwithBox())
            getRightChild()->nodeReabsorbChildren();

        if (hasLCwithBox()) {
            gatherData(dataItrs, getLeftChild());
            delete leftChild;
            leftChild = NULL;
        }

        if (hasRCwithBox()) {
            gatherData(dataItrs, getRightChild());
            delete rightChild;
            rightChild = NULL;
        }

        // reset splitDim and splitValue to their defaults
        splitDim = -1;
        splitValue = 0.0;

        leftChild = NULL;
        rightChild = NULL;
    }

    // computes a minimal subpaving from two sibling subpavings
    // a subpaving is minimal if it has no sibling leaves
    // a minimal subpaving is created by discarding sibling leaves
    // and move the the data from the children up to the new leaf
    // warning: nodeReunite would not normally be used with SPSnodes
    // but is in the base class and is reimplmented to try do it
    // appropriately for this derived class should it be needed.
    // This function is untested.

    /* Raaz's comment:
    Perhaps we should do what the base class does with the additional twist
    of passing the new information, if any, from the leaves to the parent:
    The counts and sample mean won't change but the other node features
    may change.  In particular, if we further derive SPSnode for
    functional plug-in estimation, with interval range enclosures of box
    range enclosures for some given target functional (some function
    from root box to R^m), then we want to propagate the hull of the
    intersection of the children to the parent with the an operation
    similar to the base class version that combines two leaf siblings
    into this, otherwise grafts them onto this while still propogating
    any information to the parent. */

    void SPSnode::nodeReunite(SPSnode *lChild, SPSnode *rChild)
        // lChild and rChild are the two subpavings to be reunited
    {
        // redo the box, move the data up,

        // *this is the node which will become the parent

        // check that the labels match and exit if not
        if ((lChild->label != label ) || (rChild->label != label)) {
            throw SPnodeException("Labels do not match");
        }

        // if both subpavings are leaves and hull of boxes is x,
        // discard them: *this is a leaf
        if (lChild->isLeaf() && rChild->isLeaf()) {
            if (*theBox !=
                (*(lChild->theBox) | *(rChild->theBox))) {
                throw SPnodeException("Boxes to be reunited do not fit");
            }

            // we have to collect all the data from the children,
            // and fire it into this to make sure that the
            // stats for this (this node) are correct
            NodeData tempContainer;
            gatherData(tempContainer, lChild);
            gatherData(tempContainer, rChild);

            // reserve capacity in this
            dataItrs.reserve(tempContainer.size());

            NodeDataItr it;

            SPSnode* insertedInto = NULL;

            for (it = tempContainer.begin();
                it < tempContainer.end(); it++) {
                // insert with no splitting
                SplitNever sn;
                insertedInto = insertOneFind(*it, ON_PARENT, sn);
                if (insertedInto == NULL) {
                    std::cerr << "Check "
                        << "SPSnode::nodeReunite: "
                        << " data " << **it
                        << " from nodes to be adopted "
                        << "rejected by new parent"
                        << std::endl;
                }
            }

            // the stats for this should now be right
            // and this will be a leaf so the data
            // should stay associated with it

            //discard the two subpavings given
            delete lChild;
            delete rChild;

        }

        else {  // at least one of the children is not a leaf
            // this has to adopt them rather than reuniting them
            nodeAdoptRight(rChild);
            nodeAdoptLeft(lChild);
            recursiveRename(); // recursively rename child branches
        }
    }



    // graft lChild onto this node
    // lChild could be a leaf or a non-leaf
    // takes care of the data associated with lChild/its descendents
    // used when building a statistical subpaving upwards
    // All stats are assumed to be recalculated (sums and sumproduct as well
    // as counts)
    void SPSnode::nodeAdoptLeft(SPSnode *lChild)
    {
        // *this is the node which will become the parent

        // we have to collect all the data from the child,
        // and fire it into this to make sure that the stats
        // for this (this node) are correct

        NodeData tempContainer;
        gatherData(tempContainer, lChild);

        // reserve capacity in this
        dataItrs.reserve(tempContainer.size());

        NodeDataItr it;

        SPSnode* insertedInto = NULL;

        for (it = tempContainer.begin();
            it < tempContainer.end(); it++) {

            SplitNever sn; // dummy split decision maker
            insertedInto = insertOneFind(*it, ON_PARENT, sn);
            // insert with no splitting
            if (insertedInto == NULL) {
                std::cerr << "Check SPSnode::nodeAdoptLeft: "
                    << "data " << **it << " from node "
                    << "to be adopted rejected by new "
                    << "parent" << std::endl;
            }
        }

        // the stats for this should now be right
        // but the data is associated with its descendent nodes
        // so we need to clear the actual data
        dataItrs.clear();

        // point parent and child pointers in the right directions
        // nodeAddLeft() checks labels, hull size , present children
        nodeAddLeft(lChild);
        setSplits(); // set the split dimension and split value

    }

    // graft rChild onto this node
    // rChild could be a leaf or a non-leaf
    // takes care of the data associated with lChild/its descendents
    // used when building a statistical subpaving upwards
    void SPSnode::nodeAdoptRight(SPSnode *rChild)
    {
        //* this is the node which will become the parent

        // we have to collect all the data from the child,
        // and fire it into this to make sure that the stats
        // for this (this node) are correct
        NodeData tempContainer;
        gatherData(tempContainer, rChild);

        // reserve capacity in this
        dataItrs.reserve(tempContainer.size());

        NodeDataItr it;

        SPSnode* insertedInto = NULL;

        for (it = tempContainer.begin();
            it < tempContainer.end(); it++) {

            SplitNever sn; // dummy split decision maker
            insertedInto = insertOneFind(*it, ON_PARENT, sn);
            // insert with no splitting
            if (insertedInto == NULL) {
                std::cerr << "Check SPSnode::nodeAdoptRight: "
                    << "data " << **it << " from node to "
                    << "be adopted rejected by new parent"
                    << std::endl;
            }
        }

        // the stats for this should now be right
        // but the data is associated with its descendent nodes
        // so we need to clear the actual data
        dataItrs.clear();

        // point parent and child pointers in the right directions
        // nodeAddRight() checks labels, hull size, present children
        nodeAddRight(rChild);
        setSplits(); // set the split dimension and split value
    }


    // Inserts data into this node
    // we are actually inserting an iterator to the data
    // childInd is an indicator for which child is being checked
    SPSnode* SPSnode::insertOneFind(BigDataItr newItr,
                                    OPERATIONS_ON childInd,
                                    const SplitDecisionObj& boolTest)
    {
        rvector newData = *newItr;

        // start at the top

        SPSnode* retObj = NULL;

        if(nodeContains(newData, childInd)) {

            recalculateStats(newData);

            bool wasLeaf = (isLeaf());

            // if it is a leaf, add the data and return this object
            if(wasLeaf) {

                dataItrs.push_back(newItr);

                // give this node as return value
                retObj = this;

                // split if we need to
                if (boolTest(this)) {
                    // expand and split data to children

                    nodeExpand(boolTest);

                } // end if we need to split

            } // end of isLeaf

            // if not a leaf before we had split, and contains data
            // recurse on the children if any
            if (!wasLeaf) {

                if(rightChild!=NULL && !rightChild->isEmpty()){

                    retObj =
                    (getRightChild())->insertOneFind(
                        newItr, ON_RIGHT, boolTest);
                }
                // only try left if we did not find on the right
                if(retObj == NULL && leftChild!=NULL &&
                                    !leftChild->isEmpty()) {

                    retObj =
                    (getLeftChild())->insertOneFind(newItr,
                    ON_LEFT, boolTest);
                }
            }

        } // end if node contains

        // will return null if does not contain the data

        return retObj;
    }


    // add two non-minimal pavings in a union operation,
    // return a pointer to a new non-minimal paving
    // but with no data attached to it - up to the manager to add data
    // label will be 0
    SPSnode* SPSnode:: unionTreeStructure(const SPSnode * const lhs,
                            const SPSnode * const rhs)
    {
        SPSnode* newNode = NULL;

        if ((lhs != NULL) && (rhs != NULL) && (lhs->getBox() != rhs->getBox()))
        {
            throw SPnodeException("Union unequal subpavings");
        }
        else {
            try {

                newNode = unionNoData(lhs, rhs);
                newNode->recursiveRename();
            }
            catch (bad_alloc& a) {
                cerr << a.what() << endl;
                cerr << "Error allocating memory" << endl;
                throw;
            }
            catch (SPnodeException& e) {
                string msg(e.what());
                throw SPnodeException("Error in union: original error " + msg);
            }
        }

        return newNode;

    }
    
    //src_trunk_0701
    void SPSnode::swapSPS(SPSnode& spn) //throw() // don't hide base class version
	{	
		/* theBox, parent, leftChild,
		rightChild and nodeName are inherited from base class */
		SPnode::swap(spn); // use the base version
		
		std::swap(spaceIndication, spn.spaceIndication);
		std::swap(counter, spn.counter);
		std::swap(dpSums, spn.dpSums);
		std::swap(dpSumProducts, spn.dpSumProducts);
		std::swap(dataItrs, spn.dataItrs);
		std::swap(splitDim, spn.splitDim);
		std::swap(splitValue, spn.splitValue);
		std::swap(countsOnly, spn.countsOnly);
		
	}

//src_trunk_0701
// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::SPSnode & s1, 
			subpavings::SPSnode & s2) // throw ()
{
	s1.swapSPS(s2);
}



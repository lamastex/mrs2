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

/*!/ \file:     SPSVnode.cpp
\brief SPSVnode (StatsSubPavingVal) definitions
*/

#include "spsnode.hpp"
#include "spsvnode.hpp"

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

// to use SPSVnode splitting classes
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


using namespace std;

namespace subpavings {

    // ---------------------- private member functions -------------------

    // recalculate stats for hold out data
    // recalculate the counter and accumulated sum
    // and accumulated sumproducts
    void SPSVnode::recalculateStats(rvector& newdata, bool boolVal) const
    {
	      if (boolVal==false) { counter++; } // update the counter  
			
			else { Vcounter++; } // update  the Vcounter
			
			//cout << "incrementing counters for node " << getNodeName() << endl;
    
      if (!countsOnly) {
            //cout << "mean/var calc is on" << endl;
            if (boolVal == false) {
					recalculateSums(newdata); // update the sums
					recalculateSumProducts(newdata); // update the sumproducts
			}
		}
    }
   
    // recalculate the accumulated sum
    void SPSVnode::recalculateSums(rvector& newdata) const
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
    void SPSVnode::recalculateSumProducts(rvector& newdata) const
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
    void SPSVnode::nodeExpansionOnly(int comp)
    {
        
        try
        {
            // only do something if this SPSVnode is a leaf
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
                
                nodeAddLeft(new SPSVnode(lC, space, countsOnly, label));
                nodeAddRight(new SPSVnode(rC,space, countsOnly, label));


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
                << "SPSVnode::nodeExpansionOnly()" << std::endl;
            throw;
        }

    }


   // split both validation and training data between two new children
   // uses a SplitDecisionObj to see if the children should be further split
   void SPSVnode::splitData(const SplitDecisionObj& boolTest, bool boolVal)
   {
        // check that both children exist
        if (!hasLCwithBox() || !hasRCwithBox()) {
            string msg = "Cannot split data when there are not two ";
            msg += " children";
            throw SPnodeException(msg);
        }

        NodeDataItr dataItr; // iterator
     
        boolVal = false;
        for (dataItr = dataItrs.begin();
            dataItr!= dataItrs.end(); dataItr++) {
            BigDataItr newItr = *dataItr;

            //calls insertOneFind on the children of this node
            // so stats are not recalculated for this node itself
            SPSVnode* reinsertedInto = NULL;

            if(rightChild!=NULL && !rightChild->isEmpty()) {

                reinsertedInto =
                    (getRightChild())->insertOneFind(
                    newItr, ON_RIGHT, boolTest, boolVal);
            }

            // only try the left if it's not on the right
            if(reinsertedInto==NULL && leftChild!=NULL
            && !leftChild->isEmpty()) {

                reinsertedInto =
                    (getLeftChild())->insertOneFind(
                    newItr, ON_LEFT, boolTest, boolVal);
            }
        }

        //divide the data up amongst the children
         boolVal = true;
         for (dataItr = VdataItrs.begin();
            dataItr!= VdataItrs.end(); dataItr++) {
            BigDataItr newItr = *dataItr;

            //calls insertOneFind on the children of this node
            // so stats are not recalculated for this node itself
            SPSVnode* reinsertedInto = NULL;

            if(rightChild!=NULL && !rightChild->isEmpty()) {

                reinsertedInto =
                    (getRightChild())->insertOneFind(
                    newItr, ON_RIGHT, boolTest, boolVal);
            }

            // only try the left if it's not on the right
            if(reinsertedInto==NULL && leftChild!=NULL
            && !leftChild->isEmpty()) {

                reinsertedInto =
                    (getLeftChild())->insertOneFind(
                    newItr, ON_LEFT, boolTest, boolVal);
            }
        }
        clearData();         //clear the data in this node
        
    }


    // Print the data in a node if any
    std::ostream& SPSVnode::nodeDataPrint(std::ostream &os) const
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
    std::ostream& SPSVnode::nodeMeanPrint(std::ostream &os) const
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
    std::ostream& SPSVnode::nodeVarCovarPrint(std::ostream &os) const
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
    std::ostream& SPSVnode::leafOutputTabs(std::ostream &os) const
    {
        int prec = 5; // precision for output

        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy of theBox

            // output the node name, nodeVolume, counter
            os << nodeName;
            os << "\t" << nodeVolume();
            os << "\t" << counter;
            os << "\t" << Vcounter;
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
    std::ostream& SPSVnode::leafOutputTabsWithHistHeight(const size_t bigN,
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

    // gather up all the data in children of a node
    NodeData& SPSVnode::gatherData(NodeData& container, SPSVnode * spn)
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
    void SPSVnode::setSplits()
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
    SPSVnode* SPSVnode:: unionNoData(const SPSVnode * const lhs,
                            const SPSVnode * const rhs)
    {
        SPSVnode* newNode = NULL;

        bool done = false;  // indicator for done adding

        try {

            if (lhs == NULL && rhs == NULL) done = true; // we will return NULL

            // if the lhs is null or has no box, return a tree or node based on rhs
            if (!done && (lhs==NULL || ((lhs != NULL) && (lhs->isEmpty())))) {

                newNode = SPSVnode::strippedConstructor(rhs);
                done = true;
            }

            // if the rhs is null or has no box, return a tree or node based on lhs
            if (!done && (rhs==NULL || ((rhs != NULL) && (rhs->isEmpty())))) {

                newNode = SPSVnode::strippedConstructor(lhs);
                done = true;
            }

            // by now, if we are not done, both pavings are not null and both have boxes
            // we assume that the boxes are the same

            // we have to check who has children

            // if both are leaves we can just return a node based on say lhs
             // if only rhs is leaf, lhs is not a leaf, return a node based on lhs
            if (!done && rhs->isLeaf()) {
                newNode = SPSVnode::strippedConstructor(lhs);
                done = true;
            }

            // if only lhs is leaf, rhs is not a leaf, return a node based on rhs
            if (!done && lhs->isLeaf() && !rhs->isLeaf()) {
                newNode = SPSVnode::strippedConstructor(rhs);
                done = true;
            }

            // if neither are leaves
            if (!done && !lhs->isLeaf() && !rhs->isLeaf()) {
                // make a node based on one of them, and add on the results of
                // recursing on the children
                ivector* newPermBox = new ivector(lhs->getBox());
                newNode = new SPSVnode(*newPermBox);
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

    // ------------------------ public member functions -----------------

    // Default constructor
    SPSVnode::SPSVnode() :  Vcounter(0), justSplit(false)
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
				VdataItrs.reserve(spaceIndication);
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }

    }

    // initialised constructor, initialised with one ivector for the box
    // and optionally with lab for label, defaults to 0 (see declaration)
    SPSVnode::SPSVnode(ivector& v, int lab) : SPSnode(v, lab),
														   Vcounter(0), justSplit(false)
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
    SPSVnode::SPSVnode(ivector& v, bool cntOnly, int lab) : 
                                                SPSnode(v, cntOnly, lab),
																Vcounter(0), justSplit(false)
    {
        try {
            //invokes the base class constructor with ivector, cntOnly & label
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            spaceIndication = static_cast<size_t>(defaultMaxPts);
            //reserve space - not sure if important - leave for moment
            VdataItrs.reserve(spaceIndication);
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
    SPSVnode::SPSVnode(ivector& v, size_t max, int lab) :
        SPSnode(v, max, lab), Vcounter(0), justSplit(false)
    {
        try {
            //invokes the base class constructor with ivector, max and label
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
            VdataItrs.reserve(spaceIndication+1);
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
    SPSVnode::SPSVnode(ivector& v, size_t max, bool cntOnly, int lab) :
        SPSnode(v, max, cntOnly, lab), Vcounter(0), justSplit(false)
    {
        try {
            //invokes the base class constructor with ivector, max, cntOnly and    
            //label
            //and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
          //  VdataItrs.reserve(spaceIndication+1); // this will somehow affect the constructor 
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory in SPSVnode constructor" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // initialised constructor, initialised with a LabBox (labeled box)
    // and a max for spaceIndication
    // and optionally with cntOnly for countsOnly, defaults to false
    SPSVnode::SPSVnode(LabBox& lb, size_t max, bool cntOnly) : 
			SPSnode(lb, max, cntOnly), Vcounter(0), justSplit(false)
    {
        try {

            //invokes the base class constructor with LabBox, max and cntOnly
            //and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            //reserve space - not sure if important - leave for moment
            VdataItrs.reserve(spaceIndication+1);
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
    SPSVnode::SPSVnode(LabBox& lb, bool cntOnly) : SPSnode(lb, cntOnly), 
																	Vcounter(0), justSplit(false)
    {
        try {
            //invokes the base class constructor with LabBox and cntOnly 
            //argument
            // and then initialises additional data members

            //dpSums, a vector of dotprecision terms, is not initialised
            //dpSumProducts, similarly not initialised

            spaceIndication = static_cast<size_t>(defaultMaxPts);
            //reserve space - not sure if important - leave for moment
            VdataItrs.reserve(spaceIndication);
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
    SPSVnode::SPSVnode(const SPSVnode& other) : SPSnode(*(other.theBox),
        other.label), Vcounter(other.Vcounter), justSplit(other.justSplit)
    {
        try {
			  // cout << "copy constructor for " << nodeName << " called:" << endl;
            //reserve space
            dataItrs.reserve((other.dataItrs).size());
			   VdataItrs.reserve((other.VdataItrs).size());

				//dataItrs = other.dataItrs;
				counter = other.counter;
				spaceIndication = other.spaceIndication;
            dpSums = other.dpSums;
            dpSumProducts = other.dpSumProducts; 
				splitDim = other.splitDim;
            splitValue = other.splitValue;
				countsOnly = other.countsOnly;
				
				VdataItrs = other.VdataItrs;
            nodeName = other.nodeName;;

            //recursion on the children
            if (other.leftChild) {
                nodeAddLeft(new SPSVnode(*(other.getLeftChild())));
            }
            else leftChild=NULL;

            if (other.rightChild) {
                nodeAddRight(new SPSVnode(*(other.getRightChild())));
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
    SPSVnode& SPSVnode::operator=(const SPSVnode& rhs)
    {
        try {

          //  cout << "copy assignment operator for node " << nodeName << " called:" << endl;
	   
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
            Vcounter = rhs.Vcounter;
            dpSums = rhs.dpSums;
            dpSumProducts = rhs.dpSumProducts;
            splitDim = rhs.splitDim;
            splitValue = rhs.splitValue;
            countsOnly = rhs.countsOnly;

            //reserve space
            dataItrs.reserve((rhs.dataItrs).size());
            VdataItrs.reserve((rhs.VdataItrs).size());
            //copy dataItrs from other to this
            dataItrs = rhs.dataItrs;
				VdataItrs = rhs.VdataItrs;
				
            //recursion on the children
            if (rhs.leftChild) {
                nodeAddLeft(new SPSVnode(*(rhs.getLeftChild())));
            }
            else leftChild=NULL;

            if (rhs.rightChild) {
                nodeAddRight(new SPSVnode(*(rhs.getRightChild())));
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
    SPSVnode* SPSVnode::strippedConstructor(const SPSVnode * const other)
    {
        SPSVnode* newNode = NULL;
        try {
            if (other != NULL) {
                if (other->isEmpty())
                    newNode = new SPSVnode;
                else {
                    ivector* newBox = new ivector(other->getBox());
                    newNode = new SPSVnode(*newBox);
                    newNode->splitDim = other->splitDim;
                    newNode->splitValue = other->splitValue;
                }

                newNode->nodeName = other->nodeName;
                newNode->label = 0;
                newNode->countsOnly = false;

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
    size_t SPSVnode::getVcounter() const
    { return Vcounter; }
    
    // Accessor for the node's justSplit boolean flag.
    bool SPSVnode::getJustSplit() const
    { return justSplit; }

    // Accessor for the node's data collection.
    // Returns a copy of the node's collection of iterators to the big data set.
    NodeData SPSVnode::getVdata() const
    { return VdataItrs; }
    
    // Clears the node's data collection.
    void SPSVnode::clearData() const
    { 
		dataItrs.clear(); 
      VdataItrs.clear();   
    }

    // Accessor for the parent of a node
    // Returns a copy of the pointer to parent node.
    SPSVnode* SPSVnode::getParent() const
    { return (SPSVnode*) parent; }

    // Accessor for the left child.
    // Returns a copy of the pointer to leftChild node cast to this node type
    SPSVnode* SPSVnode::getLeftChild() const
    { return (SPSVnode*) leftChild; }

    // Accessor for the right child
    //Returns a copy of the pointer to rightChild node cast to this node type
    SPSVnode* SPSVnode::getRightChild() const
    { return (SPSVnode*) rightChild; }


    // get the number of datapoints currently associated with this which would
    // be associated with the new right child if this node were to split
    // remember that the right child's box is closed at the split
    size_t SPSVnode::getRightCountIfSplit() const
    {
        // first find what the right hand child's box would be if that child
        // were to be created
            int maxdiamcomp; // variable to hold first longest dimension
            double temp = ::MaxDiam(getBox(), maxdiamcomp);

            // ivectors to be new boxes for new children
            ivector rCBox;

            // Call Upper() to get what would be the right hand child box
            Upper(getBox(), rCBox, maxdiamcomp);

        // now find how many of this node's data points would go right
        // all the rest of them would go left
        size_t rightCount = 0;
        NodeDataItr it;

        for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
            // DataItrs is a container of iterators to a BigDataCollection
            ivector pvec = _ivector((*(*it)));
            // increment rightCount if the point is in rC
            if (pvec <= rCBox) rightCount++;
        }

        return rightCount;
    }

    // The count the left child would have if this node was split.
    // Does not split the nodes, just calculates how many of the data points
    // currently associated with this node would go to the left child
    // if the node were to be split.
    size_t SPSVnode::getLeftCountIfSplit() const
    {
        return counter - getRightCountIfSplit();
    }

    // Smallest number of points in either child if this was split.
    size_t SPSVnode::getMinChildCountIfSplit() const
    {
        size_t min = getRightCountIfSplit();
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
    Size_tVec& SPSVnode::getChildrensLeftAndRightCountsIfSplit
                        (Size_tVec& grandchildCounts) const
    {
        // first find what the right hand child's box would be if that child
        // were to be created
            int maxdiamcomp1; // variable to hold first longest dimension
            double temp1 = ::MaxDiam(getBox(), maxdiamcomp1);

            // ivectors to be new boxes for new children
            ivector rCBox;
            ivector lCBox;
            // Call Upper() to get what would be the right hand child box
            Upper(getBox(), rCBox, maxdiamcomp1);
            // Call Lower() to get what would be the left hand child box
            Lower(getBox(), lCBox, maxdiamcomp1);

            int maxdiamcomp2; // variable to hold first longest dimension
            double temp2 = ::MaxDiam(rCBox, maxdiamcomp2);

            // ivectors to be new boxes for new children
            ivector rCrCBox;
            ivector rClCBox;

            // Call Upper() to get what would be the right hand child box
            Upper(getBox(), rCrCBox, maxdiamcomp2);

            // Call Lower() to get what would be the left hand child box
            Lower(getBox(), rClCBox, maxdiamcomp2);

            int maxdiamcomp3; // variable to hold first longest dimension
            double temp3 = ::MaxDiam(lCBox, maxdiamcomp3);

            // ivectors to be new boxes for new children
            ivector lCrCBox;
            ivector lClCBox;

            // Call Upper() to get what would be the right hand child box
            Upper(getBox(), lCrCBox, maxdiamcomp3);

            // Call Lower() to get what would be the left hand child box
            Lower(getBox(), lClCBox, maxdiamcomp3);


        // now find how many of this node's data points would go right
        // and left children of left and right children
        size_t rightRightCount = 0;
        size_t rightLeftCount = 0;
        size_t leftRightCount = 0;
        size_t leftLeftCount = 0;
        NodeDataItr it;

        for (it = dataItrs.begin(); it < dataItrs.end(); it++) {
            // DataItrs is a container of iterators to a BigDataCollection
            ivector pvec = _ivector((*(*it)));
            // increment rightCount if the point is in rC
            if (pvec <= rCBox) {
                if (pvec <= rCrCBox) rightRightCount++;
                else leftRightCount++;
            }
            else {
                if (pvec <= lCrCBox) rightLeftCount++;
                else leftLeftCount++;
            }
        }

        grandchildCounts.push_back(leftLeftCount);
        grandchildCounts.push_back(rightLeftCount);
        grandchildCounts.push_back(leftRightCount);
        grandchildCounts.push_back(rightRightCount);

        return grandchildCounts;
    }


    // fills in container of leaf counts, left to right
    Size_tVec& SPSVnode::getLeafNodeCounts(Size_tVec& counts) const
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


    // return a reference to a container of SPSVnodes
    // contents being the leaves descended from this, or this if this is a leaf
    // left to right order
    SPSVnodePtrs& SPSVnode::getLeaves(SPSVnodePtrs& leaves) const
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
            leaves.push_back(const_cast<SPSVnode*>(this));
        }
        return leaves;
    }

    // return a reference to a container of SPSVnodes
    // contents being the sub-leaf children of the given node
    // sub-leaf nodes are the parents of leaf nodes and only have leaf nodes
    // as children
    // left to right order
    SPSVnodePtrs& SPSVnode::getSubLeaves(SPSVnodePtrs& subleaves) const
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
            subleaves.push_back(const_cast<SPSVnode*>(this));
        }
        return subleaves;
    }

    // return a reference to a container of SPSnodes
    // contents being all the nodes in left to right order
    SPSVnodePtrs& SPSVnode::getAllNodes(SPSVnodePtrs& allNodes) const
    {
        if (!isEmpty()) { // this is not empty
		  //if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
				//cout << nodeName << endl;
            allNodes.push_back(const_cast<SPSVnode*>(this));
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




    //Returns the sum of the count over volume in the leaf nodes
    real SPSVnode::getSumLeafCountOverVol() const
    {
        dotprecision sum(0.0);

        if (isLeaf()) {  // this is a leaf
            accumulate(sum, 1.0*counter, (1.0/nodeVolume()));
        }

        else { // this is not a leaf

            SPSVnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            SPSVnodePtrsItr it;

            for(it = leaves.begin(); it < leaves.end(); it++) {
                accumulate(sum, 1.0*((*it)->getCounter()),
                            (1.0/(*it)->nodeVolume()));
            }
        }
        return rnd(sum);
    }




    //Returns the count in the smallest (by count) leaf node.
    size_t SPSVnode::getSmallestLeafCount() const
    {
        size_t smallestCount = 0;

        if (isLeaf()) {  // this is a leaf
            smallestCount = counter;
        }

        else { // this is not a leaf
            // set up a container for the leaf children
            SPSVnodePtrs leaves;
            // fill the container with the leaf children
            getLeaves(leaves);

            // find the smallest child by count
            SPSVnodePtrsItr it;
            SPSVnode* smallest = *(leaves.begin());

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
    size_t SPSVnode::getLargestLeafCount() const
    {
        size_t largestCount = 0;

        if (isLeaf()) {  // this is a leaf
            largestCount = counter;
        }

        else { // this is not a leaf

            // set up a container for the leaf children
            SPSVnodePtrs leaves;
            // fill the container with the leaf children
            // could be just this if no children
            getLeaves(leaves);

            // find the largest child by volume
            SPSVnodePtrsItr it;
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
    size_t SPSVnode::getRootCounter() const
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

    // get the count in the ultimate root node ancestor of this node
    size_t SPSVnode::getRootVcounter() const
    {
        size_t retValue = 0;
        if (parent == NULL) { // this is root
            retValue = Vcounter;
        }
        else {
            // recurse upwards
            retValue = getParent()->getRootVcounter();
        }
        return retValue;
    }


    // Get the mean of the data covered by the box of a node
    rvector SPSVnode::getMean() const
    {
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

    // Get the variance-covariance vector of the data covered
    // by the box of a node
    RealVec& SPSVnode::getVarCovar(RealVec& varCovar) const
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


    RealVec SPSVnode::getVarCovar() const
    {
        RealVec retVarCovar;
        retVarCovar = getVarCovar(retVarCovar);
        return retVarCovar;
    }

	//gat41
	// Get the uniform mean of the box of a node.
	rvector SPSVnode::getUniformMean() const
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
	real SPSVnode::getChebDistMean() const
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
	real SPSVnode::getChebDistCovar() const
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
	double SPSVnode::getEmpMass() const
	{
		int n = getRootCounter();
		double empMass = getCounter()*1.0/(n*1.0);
		//cout << nodeName << "\t" << empMass << endl;
		return empMass;
	}
	
   //gat41
   // Get the Battharchya distance.
   real SPSVnode::getHellingerDist() const
    {
		RealVec Covar = getVarCovar(); //get the covariance matrix/
		real HD = 0.0; //initialize hellinger distance to 0.

		// if there are no points, cov should be undefined. But since we want to push
		// this node to the bottom of the queue, hence let HD = 0.
		// if there is one point, the variance is 0. At the moment, we do not 
		// want to split boxes with only one point and so also let HD = 0.
		if ( getCounter() == 0 || getCounter() == 1 ) { return HD = 0.0; } 

		else {
			//cout << "===========================" << getNodeName() << "\t" << getCounter() << endl;
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
				
						cerr << "Variance cannot be negative." << endl; 
						exit(1); 
						
						/*
						gsl_matrix_free(CovarMat);
						gsl_matrix_free(UnifCovarMat);
						gsl_matrix_free(PMat);
						return HD = 0.0;
						*/ 
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
				}*/

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
					if ( DB <  0 ) { return HD = 0.0; }
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
   real SPSVnode::getHellingerDist1D() const
   {
		real HD = 0.0; 
		RealVec Covar = getVarCovar();

		// if the variance is negative - need to investiage this more
		/*
		if (Covar[0] < 0) {
			cout << getCounter() << endl;
			cout << Covar[0] << endl;
			NodeDataItr dataItr;
			cout << getNodeName() << endl;
			cout << "Data is" << std::endl;
          for (dataItr = dataItrs.begin();
                dataItr!= dataItrs.end(); dataItr++) {

                BigDataItr bigIt = *dataItr;
                rvector theData = *bigIt;

                cout << theData << endl; 
			} // end loop through data container
			
			cerr << "Variance cannot be negative." << endl; 
			exit(1); 
		}
	*/
	
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
				cerr << "no variance. check!" << endl;
				exit(1);
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
	
	RealVec SPSVnode::getUniformVarCovar() const
    {
        RealVec retVarCovar;
        retVarCovar = getUniformVarCovar(retVarCovar);
        return retVarCovar;
    }
    
	 // Get the uniform variance-covar./R	F	iance vector of the data covered
    // by the box of a node
    RealVec& SPSVnode::getUniformVarCovar(RealVec& unifVarCovar) const
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
    std::ostream& SPSVnode::nodePrint(std::ostream &os) const
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
				os << "VCounter = " << Vcounter << std::endl;
				
            nodeMeanPrint(os);
            nodeVarCovarPrint(os);
            nodeDataPrint(os);
            os << std::endl;
        }
        return os;
    }

    //Output for all the  leaf boxes in this, using tab delimiters
    std::ostream& SPSVnode::leavesOutputTabs(std::ostream &os) const
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
    //including output for the height of histogram bins for a
    // normalised histogram based on this tree with this as root
    std::ostream& SPSVnode::leavesOutputTabsWithHistHeight(
                        std::ostream &os, int prec) const
    {
        leavesOutputTabsWithHistHeight(counter, os, prec);
        return (os);
    }

    //Output for all the  leaf boxes in this, using tab delimiters
    //including output for the height of histogram bins for a
    // normalised histogram based on tree with total number of data points bigN.
    std::ostream& SPSVnode::leavesOutputTabsWithHistHeight(const size_t bigN,
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
  
    // Expand and split data
    // should be called on a node with no associated data
    void SPSVnode::nodeExpand(int comp)
    {
        nodeExpansionOnly(comp);    // expand the node
        SplitNever sn;      
        bool boolVal = true;        // dummy split decision maker
        splitData(sn, boolVal);            // split the data with no further splitting


    }

    // Expand and split the data with further splitting
    // should be called on a node with no associated data
    void SPSVnode::nodeExpand(const SplitDecisionObj& boolTest, int comp)
    {
        nodeExpansionOnly(comp);    // expand the node
        // split the data, allowing for further splitting
        bool boolVal = true; 
        splitData(boolTest, boolVal);
    }


    // Expand and split
    // finds its own comp argument
    void SPSVnode::nodeExpand()
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
        nodeExpand(maxdiamcomp); // complete nodeExpand

    }

    // Expand and split with further splitting for both training and validation   data
    void SPSVnode::nodeExpand(const SplitDecisionObj& boolTest, bool boolVal)
    {
        int maxdiamcomp; // variable to hold first longest dimension
        double temp = ::MaxDiam(getBox(), maxdiamcomp);
      
        nodeExpansionOnly(maxdiamcomp);
        justSplit = true;
        splitData(boolTest, boolVal);
    }
    
    // Expand and split with further splitting for both training and validation  data
    void SPSVnode::nodeExpand(bool boolVal)
    {
              
        int maxdiamcomp; // variable to hold first longest dimension      
        double temp = ::MaxDiam(getBox(), maxdiamcomp);        
        nodeExpansionOnly(maxdiamcomp);        
        justSplit = true;                
        SplitNever sn;
        splitData(sn, boolVal);
    }

    // reabsorb both the children of this node
    // the stats in this node will be correct so all we have to do is to
    // associate the child data with this and delete the children
    // This would work even if children are not leaves since gatherData()
    // gets data from descendents of the node it operates on
    void SPSVnode::nodeReabsorbChildren()
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
    // warning: nodeReunite would not normally be used with SPSVnodes
    // but is in the base class and is reimplmented to try do it
    // appropriately for this derived class should it be needed.
    // This function is untested.

    /* Raaz's comment:
    Perhaps we should do what the base class does with the additional twist
    of passing the new information, if any, from the leaves to the parent:
    The counts and sample mean won't change but the other node features
    may change.  In particular, if we further derive SPSVnode for
    functional plug-in estimation, with interval range enclosures of box
    range enclosures for some given target functional (some function
    from root box to R^m), then we want to propagate the hull of the
    intersection of the children to the parent with the an operation
    similar to the base class version that combines two leaf siblings
    into this, otherwise grafts them onto this while still propogating
    any information to the parent. */

  
    void SPSVnode::nodeReunite(SPSVnode *lChild, SPSVnode *rChild)
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

            SPSVnode* insertedInto = NULL;

            for (it = tempContainer.begin();
                it < tempContainer.end(); it++) {
                // insert with no splitting
                SplitNever sn;
                bool boolVal = false;
                insertedInto = insertOneFind(*it, ON_PARENT, sn, boolVal);
                if (insertedInto == NULL) {
                    std::cerr << "Check "
                        << "SPSVnode::nodeReunite: "
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
    void SPSVnode::nodeAdoptLeft(SPSVnode *lChild)
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

        SPSVnode* insertedInto = NULL;

        for (it = tempContainer.begin();
            it < tempContainer.end(); it++) {

            SplitNever sn; // dummy split decision maker
            bool boolVal = true;
            insertedInto = insertOneFind(*it, ON_PARENT, sn, boolVal);
            // insert with no splitting
            if (insertedInto == NULL) {
                std::cerr << "Check SPSVnode::nodeAdoptLeft: "
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
    void SPSVnode::nodeAdoptRight(SPSVnode *rChild)
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

        SPSVnode* insertedInto = NULL;

        for (it = tempContainer.begin();
            it < tempContainer.end(); it++) {

            SplitNever sn; // dummy split decision maker
            bool boolVal = false;
            insertedInto = insertOneFind(*it, ON_PARENT, sn, boolVal);
            // insert with no splitting
            if (insertedInto == NULL) {
                std::cerr << "Check SPSVnode::nodeAdoptRight: "
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


   // Inserts data into this node - boolVal to check if it's test or train data
   // we are actually inserting an iterator to the data
   // childInd is an indicator for which child is being checked
   SPSVnode* SPSVnode::insertOneFind(BigDataItr newItr,
                                    OPERATIONS_ON childInd,
                                    const SplitDecisionObj& boolTest, 
                                    bool boolVal)
   {
      rvector newData = *newItr;
      // start at the top
      SPSVnode* retObj = NULL;
      if(nodeContains(newData, childInd)) {
            recalculateStats(newData, boolVal);
            bool wasLeaf = (isLeaf());
            // if it is a leaf, add the data and return this object
            if(wasLeaf) { 
                if (boolVal==true) {
                   VdataItrs.push_back(newItr);
                }
                else {
                     dataItrs.push_back(newItr);
                }            
                // give this node as return value
                retObj = this;
                // split if we need to
                if (boolTest(this)) {
                    // expand and split data to children
                    nodeExpand(boolTest, boolVal);
                } // end if we need to split
            } // end of isLeaf
            // if not a leaf before we had split, and contains data
            // recurse on the children if any
            if (!wasLeaf) {
               if(rightChild!=NULL && !rightChild->isEmpty()){
                    retObj =
                    (getRightChild())->insertOneFind(
                        newItr, ON_RIGHT, boolTest, boolVal);
               }
               // only try left if we did not find on the right
               if(retObj == NULL && leftChild!=NULL &&
                                    !leftChild->isEmpty()) {
                    retObj =
                    (getLeftChild())->insertOneFind(newItr,
                    ON_LEFT, boolTest, boolVal);
               }
            }
         } // end if node contains
        // will return null if does not contain the data
        return retObj;
      }
 

      // Check if a node contains a datapoint
    // it is assumed that the node will have a box
    // childInd is an indicator for which child is being checked
    bool SPSVnode::nodeContains(const rvector& p,
                            OPERATIONS_ON childInd) const
    {
        bool retValue = false; // for the return value

        // cast p to an ivector
        ivector pvec = _ivector(p);

        //find if p is in the box
        if (pvec <= getBox()) {
            retValue = true;
        }

        // if true and there is a parent and this is a left child,
        // we need to check the split
        // find what dimension the parent was split on and what
        // the split value was
        // pvector must be strictly less than parentSplitValue
        // on the split dimension
        if (retValue && parent!=NULL && childInd == ON_LEFT) {

            int parentSplitDim = (getParent())->getSplitDim();
            real parentSplitValue = (getParent())->getSplitValue();

            if (!(p[parentSplitDim] < parentSplitValue)) {
            retValue = false;
            }
        }

        return retValue;
    }

    // add two non-minimal pavings in a union operation,
    // return a pointer to a new non-minimal paving
    // but with no data attached to it - up to the manager to add data
    // label will be 0
    SPSVnode* SPSVnode:: unionTreeStructure(const SPSVnode * const lhs,
                            const SPSVnode * const rhs)
    {
        SPSVnode* newNode = NULL;

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

} // end namespace subpavings


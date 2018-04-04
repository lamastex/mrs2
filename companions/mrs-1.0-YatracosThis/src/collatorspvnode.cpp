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

/*!/ \file:     CollatorSPVnode.cpp
\brief CollatorSPVnode definitions
*/

#include "collatorspvnode.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <exception>

// to use vectors
#include <vector>

//to use algorithms
#include <algorithm>

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"

#include "spsnode.hpp"
#include "spsvnode.hpp"

using namespace std;

namespace subpavings {


    //============================private member functions===================//
    // private initialised constructor, initialised with a pointer to an SPSnode
    // and a normalising constant, eg sum of counts in each node for a histogram
    // the summary becomes count /(normalising constant * vol) for the SPSnode
	 // Vemp becomes count/total number of validation data points
    CollatorSPVnode::CollatorSPVnode(const SPSVnode * const spn, 
	                                 size_t bigN, size_t bigM, int whatSum)
    {
        try {
         //cout << "private initialised constructor with validation summ" << endl;
            theBox = new ivector(spn->getBox());
            dimension = spn->getDimension();
            label = spn->getLabel();
            nodeName = spn->getNodeName();
       //     cout << nodeName << endl;

            if (whatSum==1){
            // add the summary to the vector summary
        // cout << "counter from private constructor: " << (spn->getCounter()) << "\t";
        // cout << "bigN: "<< bigN << "\t";
        // cout << "spn->nodeVolume(): " << spn->nodeVolume() << endl;
        // cout << "height " << spn->getCounter()/(1.0*bigN*spn->nodeVolume()) << endl;
           summary.push_back(spn->getCounter()/(1.0*bigN*spn->nodeVolume()));
            // for validation data, add the summary to the vector summary
        // cout << "Vemp from private constructor: " << 
          //         (spn->getVcounter())/(1.0 * bigM)<< endl;
				Vemp = ((spn->getVcounter())/(1.0 * bigM));        

				}
				
				else if (whatSum==2) {
					summary.push_back(spn->getCounter());
					Vemp = 0;
				}        

            //recursion on the children
            if (spn->hasLCwithBox()) {
                nodeAddLeft(new CollatorSPVnode(spn->getLeftChild(), bigN, 
                            bigM, whatSum));
            }
            else leftChild=NULL;

            if (spn->hasRCwithBox()) {
                nodeAddRight(new CollatorSPVnode(spn->getRightChild(), bigN, 
                             bigM, whatSum));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cout << msg << std::endl;
            throw;
        }
     }

    // Print the average of the summary of a single leaf node, using tab delimiters
    std::ostream& CollatorSPVnode::leafAverageOutputTabs(
                                        std::ostream &os) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            double summ = 0;
            VecDblIt it;
            // should change this to use for_each
            for (it = summary.begin(); it< summary.end(); it++) {
                summ+=(*it);
            }
            double vol = nodeVolume();
            double av =  summ/(1.0*summary.size());

            // output the nodeName, nodeVolume
            os << nodeName;
            os << "\t" << vol;
				
				// followed by the empirical measure
				os << "\t" << getVemp();
            
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

    // Print the details of a single leaf node, using tab delimiters
    // the sum of the summary is printed out
    std::ostream& CollatorSPVnode::leafAccumulationOutputTabs(
                    std::ostream &os) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            double summ = 0;
            VecDblIt it;
            // should change this to use for_each
            for (it = summary.begin(); it< summary.end(); it++) {
                summ+=(*it);
            }

            // output the nodeName, nodeVolume
            os << nodeName;
            double vol = nodeVolume();
            os << "\t" << vol;
				// followed by the empirical measure
			//	os << "\t" << getVemp();
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

 // Print the details of a single leaf node, using tab delimiters
    // the sum of the summary is printed out
    std::ostream& CollatorSPVnode::leafDifferenceOutputTabs(
                    std::ostream &os) const
    {
        if(theBox != NULL) { // do nothing if there is no box

            ivector thisBox = *theBox; // copy theBox

            double diff = summary[1] + summary[2];
           
			   // output the nodeName, nodeVolume
            os << nodeName;
            double vol = nodeVolume();
            os << "\t" << vol;
				// followed by the empirical measure
			//	os << "\t" << getVemp();
            // followed by the sum of the summary
            os << "\t" << diff;

            // followed by intervals making up box using Inf & Sup
            // ie unlike cxsc output, there is no [  ] around them
            for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {

                os << "\t" << Inf(thisBox[i])
                    << "\t" << Sup(thisBox[i]);
            }
        }
    }


    // find the accumulated summary for a node
    double CollatorSPVnode::nodeAccumulation() const
    {
            double summ = 0;
            VecDblIt it;
            // should change this to use for_each
            for (it = summary.begin(); it< summary.end(); it++) {
                summ+=(*it);
            }
            return summ;
    }

    // find the difference in summary for a node
    double CollatorSPVnode::nodeDifference() const
    {
            double diff = summary[1] + summary[2]; 
            return diff;
    }

	
    // find the absolute accumulated summary for a node
    double CollatorSPVnode::nodeAbsAccumulation() const
    {
            double absSumm = 0;
            VecDblIt it;
            // should change this to use for_each
            // add together the absolute values of each summary
            for (it = summary.begin(); it< summary.end(); it++) {
                *it < 0.0 ? absSumm-=(*it) : absSumm+=(*it);
            }
            return absSumm;
    }

    // find the accumulated summary multiplied by volume for a node
    double CollatorSPVnode::nodeAccumulationMultVol() const
    {
            return (nodeAccumulation() * nodeVolume());
    }

    // find the absolute accumulated summary multiplied by volume for a node
    double CollatorSPVnode::nodeAbsAccumulationMultVol() const
    {
            return (nodeAbsAccumulation() * nodeVolume());
    }

    // ------------------------ public member functions ----------------

    //default constructor,
     CollatorSPVnode::CollatorSPVnode() {}

    // initialised constructor, initialised with a pointer to an SPSnode
    // the summary becomes the k/(N*vol) for the SPSnode
    CollatorSPVnode::CollatorSPVnode(const SPSVnode * const spn, int whatSum)
    {
        try {
         //   cout << "default constructor for validation data" << endl;
            theBox = new ivector(spn->getBox());
            dimension = spn->getDimension();
            label = spn->getLabel();
            nodeName = spn->getNodeName();
            size_t rootCounter = spn->getRootCounter();
				size_t rootVCounter = spn->getRootVcounter(); //get the total number
				                                              // of validation data
				
				if (whatSum==1) {
					// get the empirical measure of the validation data in this node 
					Vemp = (spn->getVcounter()/(1.0*rootVCounter)); 
					// add the normalised count/volume to the vector summary					
					summary.push_back(spn->getCounter()/
                                (1.0*rootCounter * spn->nodeVolume()));       
				}
				else if (whatSum==2) {
					summary.push_back(spn->getCounter());
					Vemp = 0;
				}
				
            //recursion on the children using constructor which normalises 
				// counts
            if (spn->hasLCwithBox()) {
                nodeAddLeft(new CollatorSPVnode(spn->getLeftChild(),
                                               rootCounter, rootVCounter, 
                                               whatSum));
            }
            else leftChild=NULL;

            if (spn->hasRCwithBox()) {
                nodeAddRight(new CollatorSPVnode(spn->getRightChild(),
                                                rootCounter, rootVCounter,
                                                whatSum));
            }
            else rightChild=NULL;
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << endl;
            const char* msg = ba.what();
            std::cout << msg << std::endl;
            throw;
        }

    }

     // initialised constructor
     // initialised with a box, a label, and a vector summary
     CollatorSPVnode::CollatorSPVnode(ivector& v, int lab, VecDbl summ, 
													double Vempp) : SPnode(v, lab)
     {
         try {
           //   cout << "initialized constructor with box, label, summary, vsumm" << endl;
             // copy the vector summary
             summary = summ;
				 // copy the empirical measure
             Vemp = Vempp;
         }
         catch (bad_alloc& ba) {
             std::cout << "Error allocating memory" << endl;
             const char* msg = ba.what();
             std::cout << msg << std::endl;
             throw;
         }
     }
    
     //Copy constructor for hold out estimation
     // copies from given node downwards
     CollatorSPVnode::CollatorSPVnode(const CollatorSPVnode& other)
         : SPnode(*(other.theBox), other.label)
     {
        try {
             summary = other.summary;
             nodeName = other.nodeName;
             Vemp = other.Vemp;
 
            //recursion on the children
             if (other.leftChild) {
                nodeAddLeft(new CollatorSPVnode(
                     *(other.getLeftChild())));
           }
             else leftChild=NULL;

           if (other.rightChild) {
                 nodeAddRight(new CollatorSPVnode(
                     *(other.getRightChild())));
            }
             else rightChild=NULL;
         }
         catch (bad_alloc& ba) {
             std::cout << "Error allocating memory" << endl;
             const char* msg = ba.what();
             std::cout << msg << std::endl;
             throw;
        }
 
     }
	  
     //Copy constructor for subtracted nodes
     // copies from given node downwards
     CollatorSPVnode::CollatorSPVnode(const CollatorSPVnode& other, int toSubtract)
         : SPnode(*(other.theBox), other.label)
     {
        try {
				 
				 summary.push_back(other.nodeAccumulation()); 
				 nodeName = other.nodeName;
				 Vemp = other.Vemp;
 
            //recursion on the children
             if (other.leftChild) {
                nodeAddLeft(new CollatorSPVnode(
                     *(other.getLeftChild()), toSubtract));
           }
             else leftChild=NULL;

           if (other.rightChild) {
                 nodeAddRight(new CollatorSPVnode(
                     *(other.getRightChild()), toSubtract));
            }
             else rightChild=NULL;
         }
         catch (bad_alloc& ba) {
             std::cout << "Error allocating memory" << endl;
             const char* msg = ba.what();
             std::cout << msg << std::endl;
             throw;
        }
 
     }

     //Copy assignment operator
     // copies from given node downwards
     CollatorSPVnode& CollatorSPVnode::operator=(const CollatorSPVnode& rhs)
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
                nodeAddLeft(new CollatorSPVnode(
                     *(rhs.getLeftChild())));
             }
             else leftChild=NULL;
 
            if (rhs.rightChild) {
                 nodeAddRight(new CollatorSPVnode(
                     *(rhs.getRightChild())));
             }
             else rightChild=NULL;
         }
         catch (bad_alloc& ba) {
            std::cout << "Error allocating memory" << endl;
             const char* msg = ba.what();
             std::cout << msg << std::endl;
            throw;
         }
 
         return *this;

     }

     // make a CollatorSPVnode which represents an average of the summary
     CollatorSPVnode* CollatorSPVnode::makeAverageCollation() const
     {
         CollatorSPVnode* newnode = NULL;

         try {
             VecDbl newsumm;
             double summ = 0;
             VecDblIt it;
            // should change this to use for_each
             for (it = (summary.begin()); it < (summary.end()); it++) {
                 summ+=(*it);
             }
             newsumm.push_back(
                 summ/(1.0*(summary).size()));

             ivector v = getBox();
 
             // make the new node
             newnode = new CollatorSPVnode(v, label, newsumm, getVemp());
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
             std::cout << "Error allocating memory" << endl;
             const char* msg = ba.what();
            std::cout << msg << std::endl;
             throw;
         }
 
       return newnode;

    }

 
    // Accessor for the parent of a node.
    //Returns a copy of the pointer to parent node.
     CollatorSPVnode* CollatorSPVnode::getParent() const
     { return (CollatorSPVnode*) parent; }
 
     // Accessor for the left child of a node.
     // Returns a copy of the pointer to leftChild node, cast to this type
    CollatorSPVnode* CollatorSPVnode::getLeftChild() const
     { return (CollatorSPVnode*) leftChild; }
 
     // Accessor for the right child of a node.
     // Returns a copy of the pointer to rightChild node, cast this type
    CollatorSPVnode* CollatorSPVnode::getRightChild() const
     { return (CollatorSPVnode*) rightChild; }
 
     // Accessor for the summary.
     VecDbl CollatorSPVnode::getSummary() const
     { return summary; }
 
     // Accessor for the validation summary.
     double CollatorSPVnode::getVemp() const
     { return Vemp; }
     
     // Get number of subpavings summarised.
     size_t CollatorSPVnode::getNumberSummarised() const
     { return summary.size(); }

     // find the sum over leaf nodes of the absolute accumulated summary			     
	  // multiplied by volume for each leaf node
     double CollatorSPVnode::leavesAbsAccumulationMultVol() const
     {
         double retValue = 0;
 
         // uses  member function leafSummaryMultVolOutputTabs for node output
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
             retValue = nodeAbsAccumulationMultVol();
         }
         //recurse on the children
         if (getLeftChild()!=NULL) {
             getLeftChild()->leavesAbsAccumulationMultVol();
         }
 
         if (getRightChild()!=NULL) {
             getRightChild()->leavesAbsAccumulationMultVol();
         }
 
         return retValue;
     }
 
     // Print the details of a single leaf node, using tab delimiters
     // the summary is printed out remultiplied by the node volume
     std::ostream& CollatorSPVnode::leafOutputTabs(
                     std::ostream &os) const
     {
 
         if(theBox != NULL) { // do nothing if there is no box
 
             ivector thisBox = *theBox; // copy theBox
 
             // output the nodeName, nodeVolume
             os << nodeName;
             double vol = nodeVolume();
             os << "\t" << vol;
				 // followed by the empirical measure
				 os << "\t" << getVemp();
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
	  
	   std::ostream& CollatorSPVnode::leafOutputTabs(
                     std::ostream &os, int whichColl) const
     {
 
         if(theBox != NULL) { // do nothing if there is no box
 
             ivector thisBox = *theBox; // copy theBox
 
             // output the nodeName, nodeVolume
             os << nodeName;
             double vol = nodeVolume();
             os << "\t" << vol;
				 // followed by the empirical measure
				// os << "\t" << getVemp();
             // followed by the summary
             os << "\t" << summary[whichColl];
             
             // followed by intervals making up box using Inf & Sup
             // ie unlike cxsc output, there is no [  ] around them
             for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {
 
                 os << "\t" << Inf(thisBox[i])
                     << "\t" << Sup(thisBox[i]);
             }
 
         }
 
     }
	  
 
     //Output for all the  leaf boxes in this, using tab delimiters
     std::ostream& CollatorSPVnode::leavesOutputTabs(
                             std::ostream &os) const
     {
         // uses  member function leafSummaryMultVolOutputTabs for node output
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
     //    if (!(isEmpty())) { // this is a non-empty leaf
     //        cout << *theBox << endl; 
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
     std::ostream& CollatorSPVnode::leavesOutputTabs(
                             std::ostream &os, int whichColl) const
     {
         // uses  member function leafSummaryMultVolOutputTabs for node output
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
     //    if (!(isEmpty())) { // this is a non-empty leaf
     //        cout << *theBox << endl; 
             leafOutputTabs(os, whichColl);
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
     std::ostream& CollatorSPVnode::leavesAverageOutputTabs(
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
 
 
     //Output for all the  leaf boxes in this, using tab delimiters
     std::ostream& CollatorSPVnode::leavesAccumulationOutputTabs(
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
	  
	   //Output for all the  leaf boxes in this, using tab delimiters
     std::ostream& CollatorSPVnode::leavesDifferenceOutputTabs(
                             std::ostream &os) const
     {
         // uses  member function leafAccumulationOutputTabs for nodes
         if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
             leafDifferenceOutputTabs(os);
             return (os << "\n");
 
         }
 
             //recurse on the children
         if (getLeftChild()!=NULL) {
             getLeftChild()->leavesDifferenceOutputTabs(os);
         }
 
         if (getRightChild()!=NULL) {
             getRightChild()->leavesDifferenceOutputTabs(os);
         }
 
     }
	  
 
     // Print the details of a of a specific node in a subpaving
     std::ostream& CollatorSPVnode::nodePrint(std::ostream &os) const
     {
         // output for box in form:
         // box, volume, summary data
 
         if(theBox != NULL) { // do nothing if there is no box
 
            ivector thisBox = *theBox; // copy theBox
 
            NodeDataItr dataItr;
 
            os << "Box is :";
            for (int i = Lb(thisBox); i <= Ub(thisBox) ; i++) {
                 // c-xsc default output for intervals
                 os << "  " << thisBox[i];   
				}
 
				os << std::endl;
				os << "Box volume is " << nodeVolume() << std::endl;
				os << "Empirical measure is " << getVemp() << std::endl;
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
     // also bring down Vemp
     void CollatorSPVnode::nodeExpand(int comp)
     {
         try
         {
             // only do something if this CollatorSPVnode is a leaf
				if (isLeaf()) {
                 // ivectors to become boxes for new children
                 ivector lC, rC;
                 // Call Lower() and Upper() to put split boxes
                 // into lC and rC respectively
                 Lower(getBox(), lC, comp);
                 Upper(getBox(), rC, comp);
 
                 // make and add the new children
                 nodeAddLeft(new CollatorSPVnode(lC, label, summary, Vemp));
                 nodeAddRight(new CollatorSPVnode(rC, label, summary, Vemp));
 
                 //name the new children
                 getLeftChild()->setNodeName(nodeName + "L");
                 getRightChild()->setNodeName(nodeName + "R");

                // new children have summary from this
				}
			}

		catch (bad_alloc&)
		{
         std::cout << "Error allocating memory in "
            << "CollatorSPVnode::nodeExpand()"
            << std::endl;
         throw;
		}
	}

	// add two sibling nodes to this provided this is a leaf
	// finds its own comp argument then calls nodeExpand(int comp)
   void CollatorSPVnode::nodeExpand()
   {
      int maxdiamcomp; // variable to hold first longest dimension
		double temp = ::MaxDiam(getBox(), maxdiamcomp);
      nodeExpand(maxdiamcomp); // complete nodeExpand
 
   }
 
     // computes a minimal subpaving from two sibling subpavings
     // a subpaving is minimal if it has no sibling leaves
     // a minimal subpaving is created by discarding sibling leaves
     // and create summary data for new parent from children
     // warning: nodeReunite would not normally be used with
     // CollatorSPVnodes but is in the base class and is
     // reimplemented to try do it appropriately for this
     // derived class should it be needed.
     // This function is untested.
     //void CollatorSPVnode::nodeReunite(CollatorSPVnode *lChild,
     //                               CollatorSPVnode *rChild)
      void CollatorSPVnode::nodesReunite()
		// lChild and rChild are the two subpavings to be reunited
		{
      //  cout << "calling nodesReunite" << endl;
		 
		 /*
		  // *this is the node which will become the parent
        
		  // check that the labels match and exit if not
        if ((lChild->label != label ) || (rChild->label != label)) {
            throw SPnodeException("Labels do not match");
        }
		  
        // if both subpavings are leaves and hull of boxes is x
        // discard them: *this is a leaf
        if (lChild->isLeaf() && rChild->isLeaf()) {
			  cout << "both are leaves" << endl;
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
            // put into this summary the average of summary of children
            size_t i = 0;
				for (i=0; i < n; i++) {
					 //to take note
                summary.push_back(((lChild->summary)[i]
                                +(rChild->summary)[i])/2.0);
            }		
            //discard the two subpavings given
				delete lChild;
            delete rChild;
        }

        else {  // at least one child is not a leaf
            // this has to adopt them rather than reuniting them
				cout << "at least one child is not a leaf" << endl;
            nodeAdoptLeft(lChild);
            nodeAdoptRight(rChild);
            recursiveRename();
        }
		  */
		  
		  //only propagate if nodeAccumulation for both left and right child are
		  // the same
		  
		 // cout << getNodeName() << endl;
	
		  // first recursively deal with the children of the children		  		  
	     if (hasLCwithBox())
            getLeftChild()->nodesReunite();
        if (hasRCwithBox())
            getRightChild()->nodesReunite();

			// now deal with this
			if (hasLCwithBox() && hasRCwithBox()) {
            if (getLeftChild()->isLeaf() && getRightChild()->isLeaf()) {
					//cout << "both are leaves" << endl;
					//cout << getLeftChild()->getNodeName() << "\t" << getRightChild()->getNodeName() << endl;
				//	cout << getLeftChild()->nodeAccumulation() << "\t" << getRightChild()->nodeAccumulation() << endl;
					
					if ( ((getLeftChild()->nodeAccumulation())==0) && 
							((getRightChild()->nodeAccumulation())==0) ) {
				//		cout << "----reunite to get " << getNodeName() << endl;
						
						size_t i = 0;
						size_t n = summary.size();
						for (i=0; i < n; i++) {
							//to take note
							summary[i]=(((getLeftChild()->summary)[i]
											  +(getRightChild()->summary)[i]));
						//	cout << "reunited summary: " << ((getLeftChild()->summary)[i]
						//					  +(getRightChild()->summary)[i]) << endl;
						}		
						//discard the two children
						delete leftChild;
						delete rightChild;
						leftChild = NULL;
						rightChild = NULL;
					} 
				} 
			}
   } // end of function
	 
 // return a reference to a container of CollatorSPVnodes
    // contents being the leaves descended from this, or this if this is a leaf
    // left to right order
    std::vector<CollatorSPVnode*>& CollatorSPVnode::getLeaves(std::vector<CollatorSPVnode*>& leaves) const
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
            leaves.push_back(const_cast<CollatorSPVnode*>(this));
        }
        return leaves;
    }
    
	 
	  // return a reference to a container of SPSnodes
    // contents being all the nodes in left to right order
    std::vector<CollatorSPVnode*>& CollatorSPVnode::getAllNodes(std::vector<CollatorSPVnode*>& allNodes) const
    {
        if (!isEmpty()) { // this is not empty
		  //if (!hasLCwithBox() && !hasRCwithBox()) { // this is a leaf
            // arrgh horrible - cast away const if this node is a leaf
				//cout << nodeName << endl;
            allNodes.push_back(const_cast<CollatorSPVnode*>(this));
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
	 
	 

    // graft lChild onto this node
    // lChild could be a leaf or a non-leaf
    // takes care of the summary associated with lChild/its descendents
    // used when we are building a collator statistical subpaving upwards
    void CollatorSPVnode::nodeAdoptLeft(CollatorSPVnode *lChild)
    {
        // *this is the node which will become the parent

        size_t i = 0;

        size_t n = (lChild->summary).size();

        if (summary.empty()) { // no summary in this box already
            // put into this summary the summary of the new child
            for (i=0; i < n; i++) {
                summary.push_back((lChild->summary)[i]);
            }
        }
        else { // has summary already

            // we have to make summary for this match
            // that of the children
            // number of elements in this summary should = child
            if (summary.size() != n) {
                throw SPnodeException("Summaries do not match");;
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
    void CollatorSPVnode::nodeAdoptRight(CollatorSPVnode *rChild)
    {
        // *this is the node which will become the parent
        size_t i = 0;

        size_t n = (rChild->summary).size();

        if (summary.empty()) { // no summary in this box already
            // put into this summary the summary of the new child
            for (i=0; i < n; i++) {
                summary.push_back((rChild->summary)[i]);
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
                     
    // incorporate a subpaving to this summmary,
    // adjusts this summary for the contents of the subpaving added
    // have not specifed const data for the CollatorSPVnode pointer,
    // because if we do that we can't expand it
    // but note that the CollatorSPVnode passed in CAN BE ALTERED
    // this is for ALL nodes, i.e. including internal nodes
    bool CollatorSPVnode::addPavingWithVal(CollatorSPVnode * const spn)
    {
       
        bool retValue = false;
        bool done = false;  // indicator for done adding

        VecDblIt summaryIt;

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
                //cout << getNodeName() << "\t" << spn->getNodeName() << endl;
                ivector v = spn->getBox();
                theBox = new ivector(v);
                dimension = Ub(v) - Lb(v) + 1;
                label = spn->getLabel();
                summary = spn->summary;               
                Vemp = spn->Vemp;
                //cout << Vemp << endl;
                                
                //recursion on the children
                if (spn->leftChild) {
                    nodeAddLeft(new CollatorSPVnode(*(spn->getLeftChild())));
                }
                else leftChild=NULL;

                if (spn->rightChild) {
                    nodeAddRight(new CollatorSPVnode(*(spn->getRightChild())));
                 }
                else rightChild=NULL;

                done = true;
                retValue = true;
            } // end if theBox==NULL
            
            // do the rest only if done is not true

            // if this is a leaf and the paving to be added is a leaf we don't
            // this just sucks in the counter from spn
            if (!done && !done && isLeaf() && spn->isLeaf()) {
					 //cout << getNodeName() << "\t" << spn->getNodeName() << endl;
                //cout << "both are leaves" << endl;
                VecDbl temp = spn->getSummary();
                summary.insert(summary.end(), temp.begin(),temp.end());
                
                //danger!
                // if this is split and became a leaf because of addition into 
                // collator, then we need to use spn's Vemp
                // assuming that spn's Vemp is smaller than Vemp
					if ( Vemp > (spn->Vemp) ) {
						//cout << Vemp << endl;
						Vemp = spn->Vemp;
						//cout << Vemp << endl;
					}
               // if spn became a leaf because of addition, then we need to use 
               // Vemp of this
               // if Vemp = spn->Vemp, doesn't matter which Vemp is used
               else if ( Vemp <= (spn->Vemp) ) {
						//cout << Vemp << endl;
						Vemp = Vemp;
						//cout << Vemp << endl;
					}
               
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
                     //cout << getNodeName() << "\t" << spn->getNodeName() << endl;
                     //cout << "this is a leaf " << endl;
                     nodeExpand();
                     //cout << Vemp << endl;
							Vemp = spn->Vemp;    
							//cout << Vemp << endl;  
                }

                // if spn is leaf and this is not we need to split spn
                // THIS WILL CHANGE the CollatorSPVnode pointed to by spn

                if (spn->isLeaf()) { // so this can't be a leaf
							//cout << getNodeName() << "\t" << spn->getNodeName() << endl;
                     //cout << "spn is a leaf " << endl;
							
							spn->nodeExpand();
                     //cout << Vemp << endl;
							Vemp = Vemp;    
							//cout << Vemp << endl;   
                }

                // put in the data
                VecDbl temp = spn->getSummary();
                summary.insert(summary.end(), temp.begin(),temp.end());
                done = true;
                               
                // if they are were neither leaves originally
                // we go straight on to recursing with the children
                // otherwise expansions above are followed by recursion

                // recurse with children
                //cout << "recursing with children" << endl;
                retValue=getLeftChild()->addPavingWithVal(spn->getLeftChild());
                retValue=getRightChild()->addPavingWithVal(spn->getRightChild());
            } // end of dealing with case where at least one is not a leaf
        }
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory in addPaving" << endl;
            const char* msg = ba.what();
            std::cout << msg << std::endl;
            throw;
        }
        return retValue;
    }
    
   
   //adding CollatorSPVnodes
   CollatorSPVnode* CollatorSPVnode::addPavings(
															const CollatorSPVnode * const lhs,
															const CollatorSPVnode * const rhs)
    {
        CollatorSPVnode* newCollator = NULL;
        CollatorSPVnode* temp = NULL;
        bool done = false;
        try {
            if (lhs == NULL && rhs == NULL) done = true; // return null

            // if exactly one is null we can return a copy of the non-null one
            if (!done && lhs == NULL && rhs != NULL) {
                newCollator = new CollatorSPVnode(*rhs);
                done = true;
            }
            if (!done && lhs != NULL && rhs == NULL) {
                newCollator = new CollatorSPVnode(*lhs);
                done = true;
            }
            // both not null
            if (!done && lhs != NULL && rhs != NULL) {
                if ((lhs->getBox() != NULL) &&
                                (lhs->getBox() != rhs->getBox())) {
                    throw SPnodeException("Boxes do not match");
                }
                newCollator = new CollatorSPVnode(*lhs);
                temp = new CollatorSPVnode(*rhs);
                newCollator->addPavingWithVal(temp);
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
    CollatorSPVnode* CollatorSPVnode::subtractPavings(
                    const CollatorSPVnode * const lhs,
                    const CollatorSPVnode * const rhs, double c)
    {
        cout << "subtract pavings called" << endl;
        CollatorSPVnode* newCollator = NULL;

        bool done = false;

        try {

            if (lhs == NULL && rhs == NULL) done = true; // return null

            if (!done && lhs == NULL && rhs != NULL) {

                newCollator = new CollatorSPVnode(*rhs);
                newCollator->nodeNegate(c);

                done = true;

            }
            if (!done && lhs != NULL && rhs == NULL) {

                newCollator = new CollatorSPVnode(*lhs);
                done = true;
            }
            // both not null
            if (!done && lhs != NULL && rhs != NULL) {

                if ((lhs->getBox() != NULL) &&
                                (lhs->getBox() != rhs->getBox())) {
                    throw SPnodeException("Boxes do not match");
                }

                newCollator = new CollatorSPVnode(*lhs);
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

    // negates the summary for every node in tree rooted at this
    void CollatorSPVnode::nodeNegate(double c) 
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
    
     // incorporate the negative of a subpaving to this summmary,
    // adjusts this summary for the contents of the subpaving added
    // have not specifed const data for the CollatorSPnode pointer,
    // because if we do that we can't expand it
    // but note that the CollatorSPnode passed in CAN BE ALTERED
    void CollatorSPVnode::addNegatedPaving(const CollatorSPVnode * const spn, double c)
    {
        SPSnode* temp = NULL;

        try {
            CollatorSPVnode* temp = new CollatorSPVnode (*spn);

            // negate the node passed in
            temp->nodeNegate(c);

            addPavingWithVal(temp);

            delete temp;
        }
        catch (bad_alloc& ba) {
            std::cerr << "Error allocating memory in addNegatedPaving" << endl;
            const char* msg = ba.what();
            std::cerr << msg << std::endl;
            throw;
        }
    }

    // Get a CollatorSPVnode pointer to the corresponding SPSVnode that was split 
    bool CollatorSPVnode::getSplitNodePtrCSPV(CollatorSPVnode * &splitCollNode, SPSVnode * const spn)
    {
	//	cout << "getSplitNodePtrCSPV for: " <<spn->getNodeName() << endl;
		bool retvalue=false;
		bool done=false;
		//cout << "done on top: " << done << endl;
				
		try {
			
	         if (theBox == NULL) { done = true;
	         retvalue = true; } // if there is no memory allocated to this pointer
	        
	         else {
		//			   cout << "now looking at node: " << getNodeName() << endl;
				//		cout << "done: " << done << endl;			
					//get the pointer to the CollatorSPVnode that was split
					if ( !done && spn->getBox() == getBox() ){ 
			//			cout << "splitCollNode: " << spn->getNodeName() << "\t" << spn->getJustSplit() << endl;
						splitCollNode = (&(*this));
						done = true; 
						retvalue = true;
					}
             
					// if not done, recurse with children
					if ( done == false ) { 
				//		cout << "recursing with children" << endl;	
					//	cout << "still not done " << !done << endl;
						if (getLeftChild() != NULL&& retvalue == false) {
						retvalue = getLeftChild()->getSplitNodePtrCSPV(splitCollNode, 																spn);
					//	cout << "retvalue for left: " << retvalue << "done: " << done << endl;
						}
						if (getRightChild() != NULL && retvalue == false) {
						retvalue = getRightChild()->getSplitNodePtrCSPV(splitCollNode, 															spn);
						//cout << "retvalue for right: " << retvalue << "done: " << done << endl;
						} 
					}
				}
			}
        catch (bad_alloc& ba) {
            std::cout << "Error allocating memory in getSplitNodePtrCSPV" << endl;
            const char* msg = ba.what();
            std::cout << msg << std::endl;
            throw;
        }
        
        return retvalue;
    }


    // Find the nodes that fulfill the Scheffe condition for rows against 
	// columns in the Yatracos growing matrix. True if condition is fulfilled.
   bool CollatorSPVnode::getScheffeNode(int theta1, int theta2)
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

   // Find the nodes that fulfill the Scheffe condition for rows against 
	// columns in the Yatracos growing matrix. True if condition is fulfilled.
   bool CollatorSPVnode::nodeCheckRowSummary(int theta, int k)
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
   bool CollatorSPVnode::nodeCheckColSummary(int theta, int k)
   { 
	   //cout << "checking for Yat at node: " << getNodeName() << endl;     
      //cout << "theta: " << theta << "\t" << "k: " << k << endl;
      //cout << summary[theta] << "\t" << summary[k] << endl;
      if ((summary[k] > summary[theta])) {
		//cout << "height at " << k << " larger than height at " << theta << endl;
         return true;
      }   
		else { return false; }    
	}
     
   // get delta for a specific theta
   double CollatorSPVnode::getNodeDelta(int thisTheta)
   { 
    // cout << "get delta for " << nodeName << " at theta = " << thisTheta << endl;
     // get empirical measure of the training data
     double muTrain = summary[thisTheta] * nodeVolume();
     //cout << "summary: " << summary[thisTheta] << "\t muTrain: " << muTrain << endl;
      
     // get empirical measure of the validation data      
     //cout << "muValid: " << getVemp() << endl;

     double delta= muTrain - getVemp();
		//cout << "Delta: " << delta << endl; 

     return delta; 

   } // end of function getNodeDelta

//get the Yatracos set for a particular pair.
void CollatorSPVnode::getYatSet(
			set<CollatorSPVnode*, less<CollatorSPVnode*> > & YatSetRow, 
			set<CollatorSPVnode*, less<CollatorSPVnode*> > & YatSetCol, 
			size_t cand1, size_t cand2)
{
	//cout << "iterating through the leaves in candidate histograms" << endl;
	//iterate through the leaves in both candidate histograms to get the 
	//Yatracos set
   if (!(isEmpty()) && isLeaf()) { // this is a non-empty leaf
      bool rowInd = nodeCheckRowSummary(cand1, cand2);
     	// insert the node YatSet if return true
		if (rowInd) { 
			//cout << "inserting " << getNodeName() << " into YatSetRow" << endl; 
			YatSetRow.insert(&(*this));
		}
		else {
			bool colInd = nodeCheckColSummary(cand1, cand2);
			if (colInd) { 
				//cout << "inserting " << getNodeName() << " into YatSetCol" << endl; 
				YatSetCol.insert(&(*this));
			}
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

//get the Scheffe set for a particular pair.	
void CollatorSPVnode::getScheffeSet(
			set<CollatorSPVnode*, less<CollatorSPVnode*> > & ScheffeSet, 
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

   // ----------------- non member tools functions ----------------------

    //Output all boxes in collator
    std::ostream & operator<<(std::ostream &os,
                        const CollatorSPVnode& spn)
    {
        os << spn.nodesAllOutput(os, 1) << std::endl;
        return os;
    }
 
    // function for transform algorithm
    double opNegate(double d)
    {
        return -d;
    }

} // end namespace subpavings


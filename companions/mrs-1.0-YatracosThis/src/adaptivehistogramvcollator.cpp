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

/*! \file AdaptiveHistogramVCollator.cpp
\brief AdaptiveHistogramVCollator definitions
*/

#include "adaptivehistogramvcollator.hpp"

#include <string>   // to use the C++ string class
#include <vector>   // to use stl::vector container
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <exception> // use exceptions
#include <gsl/gsl_rng.h>        // to use the gsl random number generator
#include <gsl/gsl_randist.h>
#include <math.h> // math library
#include <limits> // to use negative infinity
// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
// to use LabBox and RSSample objects
#include "SmallClasses.hpp"
//to use subpavings
#include "sptools.hpp"
#include "spalgorithms.hpp"
// to use stats subpavings
#include "spsvnode.hpp" // includes spnode.hpp includes sptypes.hpp
#include "collatorspvnode.hpp"
#include "adaptivehistogram.hpp"
#include "adaptivehistogramvalidation.hpp"

using namespace subpavings;
using namespace std;


// ---------- implementation of AdaptiveHistogramVCollator class -------------

// ----------------private methods

// initialised constructor, initialised with a subpaving pointer
AdaptiveHistogramVCollator::AdaptiveHistogramVCollator(CollatorSPVnode * spn)
{
    if (NULL == spn) {
        throw HistException("Null CollatorSPVnode pointer in constructor");
    }
    rootVCollator = spn;
}    

// --------------- end private methods

// ---------------- public methods

// default constructor
AdaptiveHistogramVCollator::AdaptiveHistogramVCollator()
{
    try {
        rootVCollator = new CollatorSPVnode();
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cout << "Error allocating memory in constructor: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }

}

// initialised constructor, initialised with an AdaptiveHistogramVal object
AdaptiveHistogramVCollator::AdaptiveHistogramVCollator(
                                        const AdaptiveHistogramValidation& adh,
                                        int whatSum)
{
    try {
        rootVCollator = new CollatorSPVnode(adh.getSubPaving(), whatSum);
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cout << "Error allocating memory in constructor: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }
}

// copy constructor
AdaptiveHistogramVCollator::AdaptiveHistogramVCollator(
                const AdaptiveHistogramVCollator& other)
{
    try {
        rootVCollator = new CollatorSPVnode(*(other.rootVCollator));
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cout << "Error allocating memory in constructor: original error "
                                            << msg << std:: endl;
        throw HistException("Memory allocation error in constructor: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cout << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in constructor: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cout << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor: " + msg);
    }
}

// copy constructor
AdaptiveHistogramVCollator::AdaptiveHistogramVCollator(
                const AdaptiveHistogramVCollator& other, int toSubtract)
{
    try {
        rootVCollator = new CollatorSPVnode(*(other.rootVCollator), toSubtract);
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cout << "Error allocating memory in constructor: original error "
                                            << msg << std:: endl;
        throw HistException("Memory allocation error in constructor: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cout << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in constructor: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cout << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor: " + msg);
    }
}

// assignment operator
AdaptiveHistogramVCollator& AdaptiveHistogramVCollator::operator=(const
    AdaptiveHistogramVCollator& rhs)
{
    try {

        // we have to make sure we delete the current paving
        if (NULL != rootVCollator) {
            delete rootVCollator;
            rootVCollator = NULL;
        }

        if (NULL != rhs.rootVCollator)
            rootVCollator = new CollatorSPVnode(*(rhs.rootVCollator));
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cout << "Error allocating memory in assignment: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cout << "SPnodeExcepton in assignment: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in assignment: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cout << "Error in assignment: original error "
                                            << msg << std::endl;
        throw HistException("Error in assignment: " + msg);
    }
}

// increment addition operator
AdaptiveHistogramVCollator& AdaptiveHistogramVCollator::operator+=(const
    AdaptiveHistogramVCollator& rhs)
{
    try {
        // get the subpaving out of rhs to form a new CollatorSPnode
        CollatorSPVnode toAdd(*rhs.getSubPaving());
        // add the new CollatorSPnode into the collation
        // note that addPaving will alter toAdd, but that is okay because
        // toAdd is a temporary object created and deleted in this procedure
        rootVCollator->addPavingWithVal(&toAdd);
    }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory adding to collation.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    return *this;
}

// subtraction operator
AdaptiveHistogramVCollator AdaptiveHistogramVCollator::operator-(const
    AdaptiveHistogramVCollator& rhs) const
{
    if ((NULL != rootVCollator) && (NULL != rhs.rootVCollator) &&
    ((Ub(rootVCollator->getBox()) != Ub(rhs.rootVCollator->getBox()))
    || (Lb(rootVCollator->getBox()) != Lb(rhs.rootVCollator->getBox()))))
        throw HistException("Histogram collators have unequal dimensions");

       CollatorSPVnode* newnode = NULL;

    try {
		  double c = -1.0;
         newnode =
            CollatorSPVnode::subtractPavings(rootVCollator, rhs.rootVCollator, c);
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in subtraction: original error "
                                    << msg << std::endl;
        throw HistException("Memory allocation error in subtraction: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std::cerr << "SPnodeExcepton in subtraction: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in subtraction: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std::cerr << "Error in subtraction: original error "
                                            << msg << std::endl;
        throw HistException("Error in subtraction: " + msg);
    }

    AdaptiveHistogramVCollator newCollatorHist(newnode);

    return newCollatorHist;
}




// Destructor.
AdaptiveHistogramVCollator::~AdaptiveHistogramVCollator()
{ delete rootVCollator; }

// averaging method
AdaptiveHistogramVCollator AdaptiveHistogramVCollator::makeAverage() const
{
    if (NULL == rootVCollator) {
            string msg = "Cannot average this: rootVCollator is NULL";
            throw HistException(msg);
    }

    //average only makes sense if all values in the summary are positive
    VecDbl mySummary = rootVCollator->getSummary();

    VecDblIt it = find_if(mySummary.begin(), mySummary.end(), isNegative);
    if (it < mySummary.end()) {
            string msg = "Cannot average this: the collation contains negatives";
            throw HistException(msg);
    }

    AdaptiveHistogramVCollator newCollator;
    try {
        newCollator.rootVCollator = (getSubPaving())->makeAverageCollation();
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory making average.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error making average.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException making average.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error making average.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return newCollator;
}

// Return a pointer to the CollatorPSnode this manages
CollatorSPVnode* AdaptiveHistogramVCollator::getSubPaving() const
{return rootVCollator;} // boost::shared_ptr might be better

// make a .dot file for the histogram
bool AdaptiveHistogramVCollator::outputGraphDot() const
{
    if (NULL == rootVCollator) {

        throw HistException("No root paving for graph output");
    }
    bool success = rootVCollator->outputGraphDot();
    return success;
}

// Get the number of Adaptive Histogram objects collated.
size_t AdaptiveHistogramVCollator::getNumberCollated() const
{ return rootVCollator->getNumberSummarised(); }

// Output the collated normalised histogram heights and bins data to a txt file
void AdaptiveHistogramVCollator::outputToTxtTabs(const std::string& s) const
{
    // To generate a file output of the AdaptiveHistogramVCollator object
    ofstream os(s.c_str());         // Filename, c-string version

    rootVCollator->leavesOutputTabs(os);
    std::cout << "The output of the AdaptiveHistogramVCollator has been "
        << "written to " << s << std::endl << std::endl;
}

void AdaptiveHistogramVCollator::outputToTxtTabs(const std::string& s, 
																			int whichColl) const
{
    // To generate a file output of the AdaptiveHistogramVCollator object
    ofstream os(s.c_str());         // Filename, c-string version

    rootVCollator->leavesOutputTabs(os, whichColl);
    std::cout << "The output of the AdaptiveHistogramVCollator has been "
        << "written to " << s << std::endl << std::endl;
}


// Output the average data over the collation to a txt file
// this outputs the normalised average histogram heights and bins
void AdaptiveHistogramVCollator::outputAverageToTxtTabs(const
    std::string& s) const
{
    try {

        if (NULL != rootVCollator) {

            //average only makes sense if all values in the summary are positive
            VecDbl mySummary = rootVCollator->getSummary();

            VecDblIt it = find_if(mySummary.begin(), mySummary.end(), isNegative);
            if (it < mySummary.end()) {
                    string msg = "Cannot average this: the collation contains negatives";
                    std::cout << "\n" << msg << "\n" << std::endl;
            }
            else {
                // To generate a file output of the AdaptiveHistogramVCollator object
                ofstream os(s.c_str());         // Filename, c-string version

                if (NULL != rootVCollator) {
                    rootVCollator->leavesAverageOutputTabs(os);
                    std::cout << "The output of the average AdaptiveHistogram has been "
                        << "written to " << s << std::endl << std::endl;
                }
                else {
                    std::cout << "Sorry, nothing is in collation to average"
                        << std::endl;
                }
            }
        }
    }

    catch (exception& e) {
        std::cout << "Problem averaging: " << e.what() << std::endl;
    }

}

// Output the accumulated data over the collation to a txt file
// this outputs the sum over the collation summary
void AdaptiveHistogramVCollator::outputAccumulationToTxtTabs(const
    std::string& s) const
{
    try {

        if (NULL != rootVCollator) {

            // To generate a file output of the AdaptiveHistogramVCollator object
            ofstream os(s.c_str());         // Filename, c-string version

            if (NULL != rootVCollator) {
                rootVCollator->leavesAccumulationOutputTabs(os);
                std::cout << "The output of the accumulated AdaptiveHistograms "
                    << "has been written to " << s << std::endl << std::endl;
            }
            else {
                std::cout << "Sorry, nothing is in collation to accumulate"
                    << std::endl;
            }
        }
    }
    catch (SPnodeException& e) {
        std::cout << "Problem acccumulating: " << e.what() << std::endl;
    }

}

// Output the accumulated data over the collation to a txt file
// this outputs the sum over the collation summary
void AdaptiveHistogramVCollator::outputDifferenceToTxtTabs(const
    std::string& s) const
{
    try {

        if (NULL != rootVCollator) {

            // To generate a file output of the AdaptiveHistogramVCollator object
            ofstream os(s.c_str());         // Filename, c-string version

            if (NULL != rootVCollator) {
                rootVCollator->leavesDifferenceOutputTabs(os);
                std::cout << "The output of the accumulated AdaptiveHistograms "
                    << "has been written to " << s << std::endl << std::endl;
            }
            else {
                std::cout << "Sorry, nothing is in collation to accumulate"
                    << std::endl;
            }
        }
    }
    catch (SPnodeException& e) {
        std::cout << "Problem acccumulating: " << e.what() << std::endl;
    }

}

// Add an AdaptiveHistogram into the collation with Vcounter only
void AdaptiveHistogramVCollator::addToCollationWithVal(
                                        const AdaptiveHistogramValidation& adh,
                                        int whatSum,
													 size_t & agg)
{
     try {
        // make the AdaptiveHistogram into a new CollatorSPVnode
        CollatorSPVnode toAdd(adh.getSubPaving(), whatSum);
        // add the new CollatorSPVnode into the collation
        // note that addPaving will alter toAdd, but that is okay because
        // toAdd is a temporary object created and deleted in this procedure
        bool successfullyAdded = rootVCollator->addPavingWithVal(&toAdd);
		  
		  //for space complexity purposes
		  agg = spTotalNodes(rootVCollator);
		  
        if (!successfullyAdded) { // addition returned false
            std::cout << "Nothing added - check console output "
                << "for error messages" << std::endl;
        }
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory adding to collation.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
}

// Get a CollatorSPVnode pointer to the corresponding SPSVnode that was split 
void AdaptiveHistogramVCollator::getSplitNodePtr(CollatorSPVnode* &splitCollNode, SPSVnode * spn)
{
	bool success = rootVCollator->getSplitNodePtrCSPV(splitCollNode, spn);
   if (!success) { // no pointers obtained
            cerr << "No pointers obtained." << std::endl;
            exit(0);
	}
}

// Obtain the Yatracos set 
void AdaptiveHistogramVCollator::getYatracosClassAll(
CollatorSPVnode * const splitNode,
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecRowYatSet,
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecColYatSet,
list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & listYatSet )
{
	//================setting up containers====================================//
   //iterator for vector of sets
   vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator vecIt;
   set<CollatorSPVnode*, less < CollatorSPVnode* > > emptySetCSP; //empty set
   set<CollatorSPVnode*, less < CollatorSPVnode* > >::iterator setIt; //iterator
   
	//=============initializations============================================//
	int numAdd = getNumberCollated(); 
	// the number of histograms collated including the first histogram
	size_t theta = numAdd-1; // the current number of splits

  	//cout << "initialize with empty set " << endl;
    //initialize vecRowYatSet and vecColYatSet with the empty set (by definition)
      vecRowYatSet.push_back(emptySetCSP);
      vecColYatSet.push_back(emptySetCSP); 
	
	//cout << "get the left and right children of: " << splitNode->getNodeName() << endl;
    // get the left and right children of the split node
	CollatorSPVnode * lChild = splitNode->getLeftChild(); 
	CollatorSPVnode * rChild = splitNode->getRightChild();
	//cout << lChild->getNodeName() << "\t" << rChild->getNodeName() << endl;

	//============begin pairwise comparisons===================================//
	for (size_t k=0; k < theta; k++) {
	//	cout << "k= " << k << "\t theta = " << theta << endl;
		// insert the previous vecRowYatSet and vecColYatSet into listYatSet
		// remove the node that was split from vecRowYatSet/vecColYatSet 
      if (!(vecRowYatSet[k]).empty()) {
			listYatSet.push_back(vecRowYatSet[k]);  
			 //cout << "YatRowSet:" << endl;
			//for (setIt = vecRowYatSet[k].begin(); setIt != vecRowYatSet[k].end(); setIt++) { 
			//cout << (*setIt)->getNodeName() << endl; }
			vecRowYatSet[k].erase(splitNode); //if splitNode is not inside, how can it be removed? HMM!
		}  
      if (!(vecColYatSet[k]).empty()) {
			listYatSet.push_back(vecColYatSet[k]);
			//cout << "YatColSet: " << endl;
			//for (setIt = vecColYatSet[k].begin(); setIt != vecColYatSet[k].end(); setIt++) { 
			//cout << (*setIt)->getNodeName() << endl; }
			vecColYatSet[k].erase(splitNode);
		}

      // check summaries of lChild and rChild for row
		//cout << "check row: " << endl;
		bool leftRowInd = lChild->nodeCheckRowSummary(theta, k);
  		bool rightRowInd = rChild->nodeCheckRowSummary(theta, k);
     	// insert the node into vecRowYat Set if return true
     	if (leftRowInd) { 
			//cout << "inserting " << lChild->getNodeName() << " into vecRowYatSet" << endl; 
			vecRowYatSet[k].insert(lChild);
		}
		else {
			//cout << "inserting " << lChild->getNodeName() << " into vecColYatSet" << endl; 
			vecColYatSet[k].insert(lChild);
		}
		if (rightRowInd) { 
			//cout << "inserting " << rChild->getNodeName() << " into vecRowYatSet" << endl; 
			vecRowYatSet[k].insert(rChild);
		}
		else {
		//	cout << "inserting " << rChild->getNodeName() << " into vecColYatSet" << endl; 
			vecColYatSet[k].insert(rChild);
		}
		
		/*
		// check summaries of lChild and rChild for columns
		//cout << "check col: " << endl;
		bool leftColInd = lChild->nodeCheckColSummary(theta, k);
		bool rightColInd = rChild->nodeCheckColSummary(theta, k);
      // insert the node into vecColYat Set if return true
		  	if (leftColInd) { 
			cout << "inserting " << lChild->getNodeName() << " into vecColYatSet" << endl; 
			vecColYatSet[k].insert(lChild);
		}
		if (rightColInd) { 
			cout << "inserting " << rChild->getNodeName() << " into vecColYatSet" << endl; 
			vecColYatSet[k].insert(rChild);
		}*/	
			 
	} // end of pairwise comparisons

   // checking the length of vecRowYatSet and vecColYatSet
  /* if ( vecColYatSet.size() != (theta) || vecRowYatSet.size() != theta) {
	   cerr << "Length of vecRowYatSet and/or vecColYatSet is incorrect. Should have length " << numAdd << "." << endl;
		exit(0); 
   }*/

	//cout << "push back empty sets in index theta." << endl;
	//push back empty sets at index theta (because that corresponds to theta-theta comparison)
	vecColYatSet.push_back(emptySetCSP); 
	vecRowYatSet.push_back(emptySetCSP);
	//cout << "size of vecColYatSet " << vecColYatSet.size() << endl;

	//see if i can get rid of this step
	//cout << "Getting rid of repetitions" << endl;
   // get rid of repeated nodes/unions of nodes by sorting the list and 
   // checking for uniqueness
	listYatSet.sort();
	listYatSet.unique();
	
	/*list< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;     
	cout << "**Current Yatracos set has " << listYatSet.size() << " nodes." << endl;
	for (listIt = (listYatSet).begin(); listIt != (listYatSet).end(); listIt++) {
		for (setIt = (*listIt).begin(); setIt != (*listIt).end(); setIt++) {
			cout << (*setIt)->getNodeName() << endl;
		}
	}*/
} // end of function getYatracos

// get delta for each node(or union of nodes). 
double AdaptiveHistogramVCollator::getNodesDelta(
set<CollatorSPVnode*, less<CollatorSPVnode*> > & YatSet, int thisTheta)
{
  // iterator for Yatracos set
  set<CollatorSPVnode*, less<CollatorSPVnode*> >::iterator YatSetIt;  

  //gloria - think about dotprecision summation
  //initialization
  double delta = 0;
  //dotprecision deltaDP = 0;
  
  //go through each node in this set to get delta
  for (YatSetIt = YatSet.begin(); YatSetIt != YatSet.end(); YatSetIt++) {
		//cout << (*YatSetIt)->getNodeName() << endl;
		//cout << "union " << endl;
		delta += (*YatSetIt)->getNodeDelta(thisTheta);
		//accumulate(deltaDP, (*YatSetIt)->getNodeDelta(k, thisTheta), 1.0);
   }
 //  cout << "-------end of union-------" << endl;

  // take the absolute value of the sums
 //cout << "Delta: " << fabs(delta) << endl;
  return fabs(delta);
}

// get delta for each node(or union of nodes). 
//think about dotprecision summation
double AdaptiveHistogramVCollator::getNodesMaxDelta(
			vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecYatSet, 
			int thisTheta)
{
  // iterators  
  vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator YatSetIt;  

	double DeltaMax = 0;
	//dotprecision deltaDP = 0;
	//set<CollatorSPVnode*, less<CollatorSPVnode*> > YatSet;
	//cout << "Size of vecYatSet: " << vecYatSet.size() << endl;
	
	//go through each node in this set to get delta
	for (YatSetIt = vecYatSet.begin(); YatSetIt < vecYatSet.end(); YatSetIt++){
		double delta = getNodesDelta((*YatSetIt), thisTheta);
		//accumulate(deltaDP, (*YatSetIt)->getNodeDelta(k, thisTheta), 1.0);
		delta = fabs(delta);
		DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
	//	cout << "DeltaMax: " << DeltaMax << endl;
	}

  // take the absolute value of the sums
  return fabs(DeltaMax);
	
}

// get delta_theta for all theta
void AdaptiveHistogramVCollator::getYatracosDelta(
list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & listYatSet, 
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecRowYatSet, 
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecColYatSet, 
vector< vector<double> > & vecMaxDeltaVec)
{
  // iterators
  list< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator listYatSetIt;
  vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator vecIt;
  
  //set for CollatorSPVnode*
  set <CollatorSPVnode*, less<CollatorSPVnode*> > YatSet; 
  
  double delta = 0; // the delta value
  double DeltaMax = 0; // the maximum delta
  size_t theta = getNumberCollated() - 1;  // number of splits up to now
   
  /*	// check that the size of vecMaxDeltaVec is theta   
	if (vecMaxDeltaVec.size() != theta) {
		cerr << "Size of vecMaxDeltaVec before getting delta should be " << theta 
		     << "." << endl;
		exit(0);
	}*/

	//initialize vecMaxDeltaVec at the largest index
	//cout << "initialize vecMaxDeltaVec at the largest index: " << endl;
	vector<double> dummyTheta;
	for (size_t i = 0; i <=theta; i++) {
		dummyTheta.push_back(-1*(numeric_limits<double>::infinity())); // supremum of an empty set is -infinity
	}
	vecMaxDeltaVec.push_back(dummyTheta);

		// go through each histogram i to update DeltaMax for the new Scheffe 
		// sets in vecRowYatSet and vecColYatSet 
		for (size_t i = 0; i <= theta; i++) {
		//	cout << "\n =====Checking histogram  " << i << " for the " << theta << "-th split======" << endl;
			if ( i < theta) { // only can compare with deltas in indices from 0:theta-1
				DeltaMax = vecMaxDeltaVec[theta-1][i];
				//cout << "The previous Delta is: " << DeltaMax << endl;
			}
			else {
				DeltaMax = vecMaxDeltaVec[theta][theta];
			//	cout << "The previous Delta is: " << DeltaMax << endl;	
			}

			//go through each node/union of nodes in the Yatracos set
			for (vecIt = vecRowYatSet.begin(); vecIt < vecRowYatSet.end(); vecIt++){
				YatSet = *(vecIt); //dereference the pointer to set of CollatorSPVnode
				   if (!YatSet.empty()) {
					//	cout << "this set is not empty" <<endl;
						//get delta
				      delta = getNodesDelta(YatSet, i);
						//cout << "delta for this node/union of nodes is: " << delta << endl;
						//DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
						//cout << "Delta now: " << delta << "\tPrev Delta: " << DeltaMax << endl;
						DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
						//cout << "Checking the comparisons of maximums: " << endl;
						//cout << delta << "\t" << DeltaMax << endl;
					}
			//	 else { cout << " i am empty!" << endl; } 
			} // end of going through each YatSet in vecRowYatSet for histogram i

		   //cout << "DeltaMax after going through VecRowYatSet: " << DeltaMax << endl;	
			
			//go through each node/union of nodes in the Yatracos set
			for (vecIt = vecColYatSet.begin(); vecIt < vecColYatSet.end(); vecIt++){
				YatSet = *(vecIt); //dereference the pointer to set of CollatorSPVnode
				
				if (!YatSet.empty()) {
				//	cout << "this set is not empty" <<endl;
					//get delta
					delta = getNodesDelta(YatSet, i);
				//	cout << "delta for this node/union of nodes is: " << delta << endl;
            
					//DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
					DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
					//cout << "Checking the comparisons of maximums: " << endl;
					//cout << delta << "\t" << DeltaMax << endl;           
				}
				
			//	else { cout << " i am empty!" << endl; }
			} // end of going through each YatSet in vecColYatSet for histogram i
			
			// keep the updated DeltaMax in vecMaxDeltaVec at the [theta][i] position 
		//   cout << "the updated DeltaMax for the edges is: " << DeltaMax <<endl;  
			vecMaxDeltaVec[theta][i] = DeltaMax;
	} // end of going through each histogram
			
	// now go through listYatSet for the histogram at the largest index theta
	if (!(listYatSet.empty())) { //listYatSet is not empty
       DeltaMax = vecMaxDeltaVec[theta][theta];
     //  cout << "\n =====Checking histogram  " << theta << "for listYatSet=====" << endl;
       // go through each node/unions of nodes in listYatSet
       for (listYatSetIt = listYatSet.begin(); listYatSetIt != listYatSet.end(); listYatSetIt++) {
		 	YatSet = *(listYatSetIt);
				//get delta
			   delta = getNodesDelta(YatSet, theta);
       //     cout << "delta for this node/union of nodes is: " << delta << endl;
            
            //DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
            DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
            //cout << "Checking the comparisons of maximums: " << endl;
          //  cout << delta << "\t" << DeltaMax << endl;   
		} // end of going through listYatSet
	//	cout << "the updated DeltaMax for the edges is: " << DeltaMax <<endl; 
	   vecMaxDeltaVec[theta][theta] = DeltaMax; 
	} // end of if listYatSet is not empty     
} // end of method getYatracosDelta


// get delta_theta for all theta at the end of all splits
void AdaptiveHistogramVCollator::getYatracosDeltaEnd(
list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & listYatSet, 
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecRowYatSet, 
vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecColYatSet, 
vector<double> & vecMaxDelta)
{
  // iterators
  list< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator listYatSetIt;
  vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator vecIt;
  
  //set for CollatorSPVnode*
  set <CollatorSPVnode*, less<CollatorSPVnode*> > YatSet; 
  double delta = 0;
  size_t theta = getNumberCollated() - 1;  // number of splits up to now
  
  // go through each histogram i
  for (size_t i = 0; i <= theta; i++) {
	double DeltaMax = 0; // the maximum delta
    //vector<double>* myPtr = new vector<double>;

	cout << "\n =====Checking histogram  " << i << " for the " << theta << "-th split======" << endl;
	cout << vecRowYatSet.size() << "\t" << vecColYatSet.size() << "\t" << listYatSet.size() << endl;
			
			//go through each node/union of nodes in the Yatracos set
			for (vecIt = vecRowYatSet.begin(); vecIt < vecRowYatSet.end(); vecIt++){
				YatSet = *(vecIt); //dereference the pointer to set of CollatorSPVnode
				   if (!YatSet.empty()) {
					//trying: 
					//	cout << "this set is not empty" <<endl;
						//get delta
				        delta = getNodesDelta(YatSet, i);
				        //vecRowYatSet.erase(vecRowYatSet.begin());
				        //(*myPtr).push_back(delta);
						//cout << "delta for this node/union of nodes is: " << delta << endl;
						//DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
						//cout << "Delta now: " << delta << "\tPrev Delta: " << DeltaMax << endl;
						DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
						//cout << "Checking the comparisons of maximums: " << endl;
						//cout << delta << "\t" << DeltaMax << endl;
					}
			//	 else { cout << " i am empty!" << endl; } 
			} // end of going through each YatSet in vecRowYatSet for histogram i

		   //cout << "DeltaMax after going through VecRowYatSet: " << DeltaMax << endl;	
			
			//go through each node/union of nodes in the Yatracos set
			for (vecIt = vecColYatSet.begin(); vecIt < vecColYatSet.end(); vecIt++){
				YatSet = *(vecIt); //dereference the pointer to set of CollatorSPVnode
				
				if (!YatSet.empty()) {
				//	cout << "this set is not empty" <<endl;
					//get delta
					delta = getNodesDelta(YatSet, i);
			        //vecColYatSet.erase(vecColYatSet.begin());
				        
			        //(*myPtr).push_back(delta);
					//cout << "delta for this node/union of nodes is: " << delta << endl;
            
					//DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
					DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
					//cout << "Checking the comparisons of maximums: " << endl;
					//cout << delta << "\t" << DeltaMax << endl;           
				}
				
			//	else { cout << " i am empty!" << endl; }
			} // end of going through each YatSet in vecColYatSet for histogram i
			
			// keep the updated DeltaMax in vecMaxDeltaVec at the [theta][i] position 
		//   cout << "the updated DeltaMax for the edges is: " << DeltaMax <<endl;  

     //  cout << "\n =====Checking histogram  " << theta << "for listYatSet=====" << endl;
       // go through each node/unions of nodes in listYatSet
       for (listYatSetIt = listYatSet.begin(); listYatSetIt != listYatSet.end(); listYatSetIt++) {
		 	YatSet = *(listYatSetIt);
				//get delta
			   delta = getNodesDelta(YatSet, i);
     	        //listYatSet.erase(listYatSet.begin());
     	        //(*myPtr).push_back(delta);

       //     cout << "delta for this node/union of nodes is: " << delta << endl;
            
            //DeltaMax is delta if delta > DeltaPrev; else DeltaMax is still DeltaPrev
            DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
            //cout << "Checking the comparisons of maximums: " << endl;
          //  cout << delta << "\t" << DeltaMax << endl;   

	//	cout << "the updated DeltaMax for the edges is: " << DeltaMax <<endl; 

	} // end of going through listYatSet
	

	//store the DeltaMax
	//DeltaMax = *std::max_element((*myPtr).begin(), (*myPtr).end());
	vecMaxDelta.push_back(DeltaMax);   
	//delete myPtr;
	
	
	
	} // end of going through each histogram
	listYatSet.erase(listYatSet.begin(), listYatSet.end());

} // end of method getYatracosDeltaEnd


// get delta_theta for all theta
void AdaptiveHistogramVCollator::getScheffeWinner(
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecScheffeSet, 
		vector< std::vector<int> > & vecWinnerVec,
		vector< std::vector<double> > & vecDeltaWinnerVec)
{	
	//	cout << "Get Scheffe Winner" << endl;
		size_t theta = getNumberCollated() - 1;  // number of splits up to now

		// check that the size of vecMaxDeltaVec is theta   
		if (vecDeltaWinnerVec.size() != theta) {
			cerr << "Size of vecDeltaWinnerVec before getting delta should be " 
			     << theta << "." << endl;
			exit(0);
		}
		if (vecWinnerVec.size() != theta) {
			cerr << "Size of vecWinnerVec before running tournament should be "
			     << theta << "." << endl;
			exit(0);
		}	
	
   //check if the ScheffeSet for pair {theta-1, theta} is empty
	if ( vecScheffeSet[theta-1].empty() ) {
		//initialize vecDeltaWinnerVec at the largest index
		vector<double> dummyTheta;
		for (size_t i = 0; i <=theta; i++) {
			dummyTheta.push_back(-1*(numeric_limits<double>::infinity())); 
			// supremum of an empty set is -infinity
		}
		vecDeltaWinnerVec.push_back(dummyTheta);
		//initialize at the previous indices
		for (size_t i = 0; i < theta; i++){
			vecDeltaWinnerVec[i].push_back(-1*(numeric_limits<double>::infinity()));
		}		
		//push back 0 in vecWinnerVec at the largest index
		vector<int> vecWinner;
		for (size_t i = 0; i <=theta; i++) {
			vecWinner.push_back(0); 
		}
		vecWinnerVec.push_back(vecWinner);
		//push back 0 in vecWinnerVec at the previous indices
		for (size_t i = 0; i < theta; i++){
			vecWinnerVec[i].push_back(0);
		}
	} // end of vecScheffeSet[theta-1] is empty	
	
	else {
		// iterators
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >::iterator vecIt;
		set <CollatorSPVnode*, less<CollatorSPVnode*> >::iterator setIt; 

		//initialize vecDeltaWinnerVec at the largest index
		vector<double> dummyTheta;
		for (size_t i = 0; i <=theta; i++) {
			dummyTheta.push_back(-1*(numeric_limits<double>::infinity())); 
			// supremum of an empty set is -infinity
		}
		vecDeltaWinnerVec.push_back(dummyTheta);
		//initialize at the previous indices
		for (size_t i = 0; i < theta; i++){
			vecDeltaWinnerVec[i].push_back(-1*(numeric_limits<double>::infinity()));
		}
		
		//initialize vecWinnerVec at the largest index
		vector<int> vecWinner;
		for (size_t i = 0; i <=theta; i++) {
			vecWinner.push_back(0); 
		}
		vecWinnerVec.push_back(vecWinner);		
		//initialize vecWinnerVec at the previous indices
		for (size_t i = 0; i < theta; i++){
			vecWinnerVec[i].push_back(0);
		}
				
		//get the current theta		  
		double deltaTheta = getNodesDelta(vecScheffeSet[theta-1], theta);    
		
		// go through each histogram i to get the winner
		for (size_t i = 0; i < theta; i++) {
				if (!(vecScheffeSet[i]).empty()) {				
					//get delta
					double delta = getNodesDelta(vecScheffeSet[theta-1], i);
					//perform competition
					bool i_win = delta < deltaTheta;
					if ( i_win == true ) {
						cout << i << "\t" << delta << "\n" << theta << "\t" << deltaTheta << endl;
						cout << "Winner is: " << i << endl;
						// winner is i
						vecWinnerVec[theta][i]=(0);		
						vecWinnerVec[i][theta]=(1);								
						vecDeltaWinnerVec[theta][i]=(deltaTheta);		
						vecDeltaWinnerVec[i][theta]=(delta);			
					}
					else if ( i_win == false ) { 
						cout << i << "\t" << delta << "\n" << theta << "\t" << deltaTheta << endl;
						cout << "Winner is: " << theta << endl;
						// winner is theta
						vecWinnerVec[i][theta]=(0);
						vecWinnerVec[theta][i]=(0);						
						vecDeltaWinnerVec[i][theta]=(delta);		
						vecDeltaWinnerVec[theta][i]=(deltaTheta);		
					}
					else { // deltaTheta = delta 
					   cout << i << "\t" << delta << "\n" << theta << "\t" << deltaTheta << endl;
						cout << "tie" << endl;
						vecWinnerVec[theta][i]=(0);
						vecWinnerVec[i][theta]=(0);						
						vecDeltaWinnerVec[theta][i] = (delta);		
						vecDeltaWinnerVec[i][theta]= (delta);								
					}  // end of competition
				} // end of ScheffeSet not empty				           
		} // end of going through each histogram				
	} // end of vecScheffeSet[theta-1] not empty   
} // end of getScheffeDelta

//get the theta of the minimum of delta_theta
void AdaptiveHistogramVCollator::getMinDistTheta(									
						std::vector< std::vector<int> > & vecMinDistTheta, 
   std::vector< std::vector<double> > & vecMaxDeltaVec, int n)
{  
//   cout << "Get the minimum:" << "\t" << endl;
	//want the vector of DeltaMax at the last entry of vecMaxDeltaVec 
   vector<double> vecDeltaMax = vecMaxDeltaVec.back();
   double DeltaInf = *(min_element(vecDeltaMax.begin(), vecDeltaMax.end()));
   //cout << DeltaInf << endl;
	//double upperLimit = DeltaInf + 1.0/n*1.0;

   // keep the associated theta in a vector and then push back in the bigger vector
   vector<int> minDistTheta;
   
     for (size_t i=0; i < vecDeltaMax.size(); i++) {
      //   cout << "i: " << i << "\t" << DeltaInf << " vs " << vecDeltaMax[i] << endl;
         if (vecDeltaMax[i] == DeltaInf) { 
			//if (vecDeltaMax[i] < upperLimit) 
      //       cout << "this is at split " << i << endl;
             minDistTheta.push_back(i); 
			}  
		}
    vecMinDistTheta.push_back(minDistTheta);
 } 

// get delta_theta for all theta
void AdaptiveHistogramVCollator::getHistScheffeWinner(
		vector< vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > > & vecScheffeSetVec, 
		vector< std::vector<int> > & vecWinnerVec,
		vector< std::vector<double> > & vecDeltaWinnerVec)
{	
	// go through each ordered pair to get the winner
	for (size_t i = 0; i < vecScheffeSetVec.size(); i++) {
		
		vector<int> WinnerVec(vecScheffeSetVec[i].size());
		vector<double> DeltaVec;

		for (size_t j = 0; j < vecScheffeSetVec[i].size(); j++) {

			if ( vecScheffeSetVec[i][j].empty() ) {
				//size_t cand1 = i;
				//size_t cand2 = j+i+1;
				WinnerVec[j] = (-1);
				DeltaVec.push_back(-1*(numeric_limits<double>::infinity()));
				//cout << "no scheffe set at position " << i << "\t" << j << endl;
				//cout << "========" << endl;
			} 
		
			else {
				//cout << "scheffe at position " << i << "\t" << j << endl;
				size_t cand1 = i;
				size_t cand2 = j+i+1;
				cout << cand1 << "\t" << cand2 << endl;

				//cout << "----------get delta for " << cand1 << endl;
				double deltaI = getNodesDelta(vecScheffeSetVec[i][j], cand1);
				//cout << "---------get delta for " << cand2 << endl;
				double deltaJ = getNodesDelta(vecScheffeSetVec[i][j], cand2);

				// perform competition
				if ( deltaI < deltaJ ) {
					//cout << cand1 << "\t" << deltaI << "\n" << cand2 << "\t" << deltaJ << endl;
					//cout << "Winner is: " << cand1 << endl;
					// winner is i
					WinnerVec[j] = (1);
					DeltaVec.push_back(deltaI);
				}
				else { // deltaTheta >= delta 
					//cout << cand1 << "\t" << deltaI << "\n" << cand2 << "\t" << deltaJ << endl;
					//cout << "Winner is: " << cand2 << endl;
					WinnerVec[j] = (0);
					DeltaVec.push_back(deltaJ);
				}  // end
			} // end of set not empty
		} // end of going through j
		vecWinnerVec.push_back(WinnerVec);
		vecDeltaWinnerVec.push_back(DeltaVec);
	} // end of going through vecScheffeSet

} // end of getHistScheffeWinner




//get the minimum delta_theta
bool AdaptiveHistogramVCollator::getMinDelta(int maxCheck,
												vector< vector<double> > & vecMaxDeltaVec)
{
	//want the vector of DeltaMax at the last entry of vecMaxDeltaVec 
	vector<double> vecDeltaMax = vecMaxDeltaVec.back();
	double minDelta = *(min_element(vecDeltaMax.begin(), vecDeltaMax.end()));
	size_t minPos = -1;
	
	vector<bool> signs;
	//subtract each element from vecDeltaMax
	for (size_t i = 0; i < vecDeltaMax.size(); i++) {
//		cout << i << "\t" << vecDeltaMax[i] << "\t" << minDelta << endl;
		if ( vecDeltaMax[i]  > minDelta ) { signs.push_back(1); }
		else if ( vecDeltaMax[i] == minDelta ) {  
			signs.push_back(0); 
			// get the position at which the first time this happens
			if ( minPos == -1 ) { minPos = i; }
		}

		if ( vecDeltaMax[i] < minDelta ) {
			cerr << "vecDeltaMax[i] < minDelta!" << endl;
			exit(0);
		}
	}
	
//	cout << "Min Delta at position " << minPos << endl;

	bool toStop = false;

	//if the infimum is at the last theta - continue splitting
//	cout << "======check if minDelta is at the last row====" << endl;
	if ( minPos == (getNumberCollated()-1) ) { 
//		cout << "continue splitting" << endl;
		return toStop;
	}
	
	else if ( (getNumberCollated()-1) != minPos ) {
//	cout << "======minDelta is not in the last row, check all the following deltas====" << endl;
//	cout << "=======now check if at least k of the following deltas are equal or greater =====" << endl;
		int flag = 0;
		for (size_t i = (minPos+1); i < signs.size(); i++) {
//			cout << "i: " << i << " sign: " << signs[i] << endl;
			flag++;
			cout << "flag:"  << flag << endl;
		}
	
		if ( flag >= maxCheck ) { return toStop = true;}
		else { return toStop = false; }
	} // end of else if
} // end of function

 // Get infimum delta_theta
 vector<double> AdaptiveHistogramVCollator::getInfDelta(
												vector<double> & vecInfDelta, 
												vector< vector<double> > & vecMaxDeltaVec, 
												int n)
 {
   //want the vector of DeltaMax at the last entry of vecMaxDeltaVec 
   vector<double> vecDeltaMax = vecMaxDeltaVec.back();

  //  cout << "getting the upper limit" << endl;
  //  double DeltaInf = *(min_element(vecDeltaMax.begin(), vecDeltaMax.end()));
    double DeltaInf = *(min_element(vecDeltaMax.begin(), vecDeltaMax.end()));
    
    vecInfDelta.push_back(DeltaInf);      
    
    return vecInfDelta;
}

// get the root box of the collator
ivector AdaptiveHistogramVCollator::getRootBox() 
{
	ivector rootBox = rootVCollator->getBox();	
	return rootBox;
}

// Get Scheffe Tournament All Winner
void AdaptiveHistogramVCollator::getScheffeSetAll(
												CollatorSPVnode * const splitNode,
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecScheffeSet,
		list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & listScheffeSet)
{
	//================setting up containers====================================//
  // cout << "set up containers" << endl;
	//iterator for vector
   vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator vecIt;
   	
   //set for CollatorSPVnode* and iterator for sets
   set<CollatorSPVnode*, less < CollatorSPVnode* > > setCSP;
   set<CollatorSPVnode*, less < CollatorSPVnode* > >::iterator setIt;
   	
   //empty set
   set<CollatorSPVnode*, less < CollatorSPVnode* > > emptySetCSP;
   //=============end of setting up containers================================//
	//cout << "getNumberCollated: " << endl; 	
	int numAdd = getNumberCollated(); // the number of histograms collated including the 0-th histogram
	size_t theta = numAdd-1; // the current number of splits
	
   // initialize the vecScheffeSet
   if ( splitNode->getBox() == getRootBox() ) { //initialize vecScheffeSet with the empty set (by definition)
      setCSP.insert(splitNode); // urm what is this for?
      vecScheffeSet.push_back(emptySetCSP); 
	}
	//cout << "get the left and right children of: " << splitNode->getNodeName() << endl;		
   // get the left and right children of the split node
  	CollatorSPVnode * lChild = splitNode->getLeftChild(); 
  	CollatorSPVnode * rChild = splitNode->getRightChild();
  	//cout << lChild->getNodeName() << "\t" << rChild->getNodeName() << endl;
 	
  	//============begin pairwise comparisons===================================//
  	for (size_t k=0; k < theta; k++) {	
   //   cout << "k= " << k << endl; 	   
	   // insert the previous vecScheffeSet into listScheffeSet	and
		// remove the node that was split from vecScheffeSet 
      if (!(vecScheffeSet[k]).empty()) {
         listScheffeSet.push_back(vecScheffeSet[k]);  
      //  for (setIt = vecScheffeSet[k].begin(); setIt != vecScheffeSet[k].end(); setIt++) { cout << (*setIt)->getNodeName() << endl; }
			vecScheffeSet[k].erase(splitNode);
		}	 
   
	   // check summaries of lChild and rChild for row
		bool leftRowInd = lChild->getScheffeNode(k, theta);
  		bool rightRowInd = rChild->getScheffeNode(k, theta);
     	// insert the node into vecScheffeSet if return true
     	if (leftRowInd) { 
	//		cout << "inserting " << lChild->getNodeName() << " into vecScheffeSet" << endl; 
			vecScheffeSet[k].insert(lChild);
		}
		if (rightRowInd) { 
	//		cout << "inserting " << rChild->getNodeName() << " into vecScheffeSet" << endl; 
			vecScheffeSet[k].insert(rChild);
		}
	} // end of pairwise comparisons

   // checking the length of vecRowYatSet and vecColYatSet
   if ( vecScheffeSet.size() != theta) {
	   cerr << "Length of vecRowYatSet is incorrect. Should have length " << numAdd << "." << endl;
		exit(0); 
   }
		
	//push back empty sets at index theta (because that corresponds to theta-theta comparison)
	vecScheffeSet.push_back(emptySetCSP); 
			   
   // get rid of repeated nodes/unions of nodes
  	// sort the list
  	listScheffeSet.sort();
	
	// check for uniqueness
	listScheffeSet.unique();
	
} // end of function getScheffeSetAll


// Get Scheffe set from sub-pavings.
void AdaptiveHistogramVCollator::getHistScheffeSet(
		vector < vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > > & vecScheffeSetVec)
{
   //=============end of setting up containers================================// 
	int numAdd = getNumberCollated(); // the number of histograms collated including the 0-th histogram
	//cout << "getNumberCollated: " << numAdd << endl;
	size_t theta = numAdd-1; // the current number of splits

  	//============begin pairwise comparisons===================================//
  	for (size_t k=0; k < numAdd; k++) {
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecScheffeSet;
      for (size_t j = 1; j < numAdd; j++) {
			if ( (k != j) && (k < j) ) {
				set<CollatorSPVnode*, less < CollatorSPVnode* > > currentScheffeSet;
				cout << "k= " << k << "\t" << "theta = " << j << endl;
				getSubPaving()->getScheffeSet(currentScheffeSet, k, j);
				//if (currentScheffeSet.empty()) { cout << "nothing here" << endl; }
				vecScheffeSet.push_back(currentScheffeSet);
			}
		}
		vecScheffeSetVec.push_back(vecScheffeSet);
	} // end of pairwise comparisons
} // end of function getHistScheffeSet

void AdaptiveHistogramVCollator::getHistYatSet(
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecYatSet)
{
   //=============end of setting up containers================================// 
	int numAdd = getNumberCollated(); // the number of histograms collated including the 0-th histogram
	//============begin pairwise comparisons===================================//
	for (size_t k= 0; k < numAdd; k++) {
		// get A_ij
      for (size_t j = 0; j < numAdd; j++) {
			if ( (k != j) && (k < j) ) {
				set<CollatorSPVnode*, less < CollatorSPVnode* > > RowSet;
				set<CollatorSPVnode*, less < CollatorSPVnode* > > ColSet;
				//cout << "k= " << k << "\t" << "theta = " << j << endl;	
				getSubPaving()->getYatSet(RowSet, ColSet, k, j);
				vecYatSet.push_back(RowSet);
				vecYatSet.push_back(ColSet);
			}
		}
	} // end of pairwise comparisons
	
	sort( vecYatSet.begin(), vecYatSet.end() );
	vecYatSet.erase( unique( vecYatSet.begin(), vecYatSet.end() ), vecYatSet.end() );
} // end of function getHistYatSet

// Get winner of the knock out.
/*
void AdaptiveHistogramVCollator::getBisectionSearchEstimate(
														int method,
														std::vector<size_t> & winner, 
														std::vector<size_t> & loser,
														std::vector<double> & deltaWinner,
														std::vector<double> & deltaLoser)
{
	//growing Yatracos list 
	list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > listYatSet;	
	//iterators
	set<CollatorSPVnode*, less<CollatorSPVnode*> >::iterator setIt;
     
   //initialize candidates
   //always let cand1 be winner and cand2 be updated
   size_t cand1 = 0, cand2 = getNumberCollated() - 1;
   int finalWinner;
   
	//start the match   
   while ( (cand1 >= 1 && cand2 <= (getNumberCollated()-2)) || 
			  (cand2 >= 1 && cand1 <= (getNumberCollated()-2)) ) {	
		if (cand1 != cand2) {
			cout << "----------------------------------------------------" << endl;		
			cout << "Cand1: " << cand1 << "\t" << "Cand2: " << cand2 << endl;
			//set up a vector for the Yatracos set	
			set<CollatorSPVnode*, less<CollatorSPVnode*> > ScheffeSet;
			set<CollatorSPVnode*, less<CollatorSPVnode*> > YatSetRow;  
			set<CollatorSPVnode*, less<CollatorSPVnode*> > YatSetCol;
						
			if (method==1) {
				//get the Scheffe set
				if (cand1 > cand2) { //for the sake of definition
					(*this).getSubPaving()->getScheffeSet(ScheffeSet, cand2, cand1);
				}
				else { 
					(*this).getSubPaving()->getScheffeSet(ScheffeSet, cand1, cand2); 
				}
			}
			else if (method==2) {
				//get the Yatracos set
				(*this).getSubPaving()->getYatSet(YatSetRow, YatSetCol, cand1, cand2);
				listYatSet.push_back(YatSetRow);
				listYatSet.push_back(YatSetCol);
				listYatSet.sort();
				listYatSet.unique();
			}
			
			if ( YatSetRow.empty() && YatSetCol.empty() && ScheffeSet.empty() ) { // check if it's empty
				cout << "No Scheffe/Yatracos set for this round. Bisect again." << endl;
				size_t temp = ceil((cand2+cand1)/2.0);
			//	cout << "Potential candidate is: " << temp << endl;
				if ( cand1 == temp || cand2 == temp ) {  
						cout << "cand1 == temp || cand2 == temp. Cannot compete " <<
									"against ownself. Tournament ended." << endl;
						break; 
				} // to prevent repeating candidates
				else {
					cand2 = temp;
				//	cout << "New candidate is: " << cand2 << endl;
				} 
			}
			else {
				double deltaCand1, deltaCand2;
				if (method==1) {
					//get delta for cand1
					deltaCand1 = getNodesDelta(ScheffeSet, cand1);
					//get delta for cand2
					deltaCand2 = getNodesDelta(ScheffeSet, cand2);
				}
			   else if (method==2) {					
					//get max delta for cand1				
				//	cout << "\nMaxDelta for " << cand1 << endl; 	
					deltaCand1 = getNodesMaxDelta(listYatSet, cand1);					
					//get max delta for cand2
				//	cout << "\nMaxDelta for " << cand2 << endl;
					deltaCand2 = getNodesMaxDelta(listYatSet, cand2);
				}
										
				//get the minimum of deltaCand1 and deltaCand2 along with the 
				//corresponding theta that gave the minimum
				if ( deltaCand1 < deltaCand2 ) { 
					cout << cand1 << "\t" << cand2 << endl;
					cout << deltaCand1 << "\t" << deltaCand2 << endl;
					cout << "Winner is " << cand1 << endl;
					
					//put into containers
					winner.push_back(cand1);
					loser.push_back(cand2);
					deltaWinner.push_back(deltaCand1);
					deltaLoser.push_back(deltaCand2);
													
					//winner = cand1 so change cand2
					size_t temp = ceil((cand2+cand1)/2.0);
				//	cout << "Potential candidate is: " << temp << endl;
					if ( cand1 == temp || cand2 == temp ) { break; } // to prevent repeating candidates
					else {
						finalWinner = cand1;
						cand2 = temp;
				//		cout << "New candidate, cand2 is: " << cand2 << endl;
					} 
				}
				else if ( deltaCand2 < deltaCand1 ) { 
					cout << cand1 << "\t" << cand2 << endl;
					cout << deltaCand1 << "\t" << deltaCand2 << endl;
					cout << "Winner is " << cand2 << endl;
					
					//put into containers
					winner.push_back(cand2);
					loser.push_back(cand1);
					deltaWinner.push_back(deltaCand2);
					deltaLoser.push_back(deltaCand1);										
					
					//winner = cand2 so change cand1
					size_t temp = ceil((cand2+cand1)/2.0);
				//   cout << "Potential candidate is: " << temp << endl;
					if ( cand1 == temp || cand2 == temp ) { break; } // to prevent repeating candidates
					else {
						cand1 = cand2; //let cand1 be the winner
						cand2 = temp;  //let cand2 be the updated candidate
						finalWinner = cand1;						
				//		cout << "New candidate is: " << cand2 << endl;
					} 
				}
				//maybe can use this as stopping criteria
				else {
					cout << "Both deltas are the same." << endl;
					//bisect again
					size_t temp = ceil((cand2+cand1)/2.0);
				//	cout << "Potential candidate is: " << temp << endl;
					if ( cand1 == temp || cand2 == temp ) {  
						cout << "cand1 == temp || cand2 == temp. Cannot compete " <<
									"against ownself. Tournament ended." << endl;
						break;
					} // to prevent repeating candidates
					else {					
						//cand1 should be the one with the smallest index (break ties by taking the smallest index)					
						if ( cand1 > cand2) { cand1 = cand2; }					
						cand2 = temp;
				//      cout << "New candidate is: " << cand2 << endl;				
					} 
				}	   
			} // end of if Yatracos set not empty
		} // end of cand1 != cand2
	} // end of while loop
}
*/

//make minimal sub-pavings
void AdaptiveHistogramVCollator::makeMinimal()
{
  try {
		  rootVCollator->nodesReunite();	  
   }
   catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory adding to collation.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error adding to collation.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
}
		 
	 // Get a string of the leaf node levels.
    std::string AdaptiveHistogramVCollator::getLeafLevelsString() const
   {
		string retValue = "";
		if (NULL != rootVCollator)
        retValue = rootVCollator->getLeafNodeLevelsString();

		return retValue;
	}

	//get total number of nodes
	size_t AdaptiveHistogramVCollator::getTotalNodes()
	{
		size_t numNodes = 0;
		if (NULL != rootVCollator)
        numNodes = spTotalNodes(rootVCollator);
		
		return numNodes;
	}


void AdaptiveHistogramVCollator::getMinDistEst(vector<double> & maxDelta, 	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecYatSet)
{
	//get the yatracos class for ALL the candidates
	getHistYatSet(vecYatSet); 
	
	//get the maximum delta at each "theta" - here theta refers to the position of the
	//candidate in the collator
	//cout << "Number collated: " << getNumberCollated() << endl;
	for (int i = 0; i < getNumberCollated(); i++) {
		//get the maximum delta at this candidate
		double deltaMax = getNodesMaxDelta(vecYatSet, i);
		maxDelta.push_back(deltaMax);
	}

}
  
// ---------- end implementation of AdaptiveHistogramVCollators -----------

//Output all boxes in AdaptiveHistogramVCollator adhc
std::ostream & operator<<(std::ostream &os,
                    const AdaptiveHistogramVCollator& adhc)
{
    if (NULL != adhc.getSubPaving()) {
        os << (adhc.getSubPaving())->nodesAllOutput(os, 1) << std::endl;
    }

    return os;
}

    // function for for_each algorithm
    double isNegative(double d)
    {
        return (d<0.0);
    }


   //output all boxes in collator to text file
   void outputAllNodesToTxtTabs(const std::string& s, const AdaptiveHistogramVCollator& adhc) 
   {
    // To generate a file output of the AdaptiveHistogramVCollator object
    ofstream os(s.c_str());         // Filename, c-string version

    adhc.getSubPaving()->nodesAllOutput(os, 1);
    std::cout << "The output of ALL the nodes of the AdaptiveHistogramVCollator has been "
        << "written to " << s << std::endl << std::endl;
    }

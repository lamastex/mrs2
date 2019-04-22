/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
* Copyright (C) 2011 Gloria Teng
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

/*! \file AdaptiveHistogramValidation.cppas
\brief AdaptiveHistogramValidation definitions
*/

#include "adaptivehistogramvalidation.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <exception> // use exceptions

#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h> // to use the constant M_PI 
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
#include "spsnode.hpp"
#include "collatorspvnode.hpp"
#include "adaptivehistogram.hpp"
#include "adaptivehistogramvcollator.hpp"

// to use mcmc function objects
#include "histmcmcobjs.hpp"

// to use histogram evaluation objects
#include "histevalobjval.hpp"

// to use error functions 
#include "errorfunc.hpp"

// to use 2D integration using taylor methods
#include "../examples/StatsSubPav/ExactInt/Int.h"
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"

// to use assert
#include <assert.h>

using namespace subpavings;
using namespace std;

// a class for comparison between spsvnodes
class MyCompare
{
    const NodeCompObjVal& myNC;

    public:
    MyCompare(const NodeCompObjVal& nc) : myNC(nc) {}

    bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
    { return myNC(lhs, rhs); }

};

//=====implementation of AdaptiveHistogramValidationValidation class ==========

// --------------------------- private ---------------------------------------

// a constant for padding a box if it is tailor-made for data
const real AdaptiveHistogramValidation::padding = 0.000005;

// initialised constructor, initialised with a subpaving pointer
AdaptiveHistogramValidation::AdaptiveHistogramValidation(SPSVnode * spn, bool as)
        : holdAllStats(as)
{
    if (NULL == spn) {
        throw HistException("Cannot use null SPSVnode pointer in constructor");
    }
    rootVpaving = spn;
    creationString = rootVpaving->getNodeName();
    creationString += rootVpaving->getChildNodeNames();

    rootBox = spn->getBox();

    // nothing happens to dataCollection when object is constructed
}


// complete insertion of data from a vector of data for hold out
// given a container of rvectors of the data  to insert
bool AdaptiveHistogramValidation::completeDataInsertionFromVec(const RVecData& theData,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging, size_t holdOutCount, 
										  std::vector<size_t> & numNodes)
{

    bool retValue = false;

    //find the data dimensions from the first datapoint
    int dataDim = Ub(*theData.begin()) - Lb(*theData.begin()) + 1;

    // ensure the paving exists
    bool hadToMakePaving = haveMadePaving(theData, dataDim);

    // if we did not make the paving we have to check data dimensions
    if (!hadToMakePaving) {
        if(dataDim != (Ub(rootBox) - Lb(rootBox)) + 1) {

            throw HistException("Dimensions of data do not match paving");
        }
    }

    // insert the data
  //  std::cout << "calling insertData..." << endl;
    size_t dataCountInserted
            = insertDataFromContainer(theData, boolTest, logging, holdOutCount,
				                          numNodes);

    if (dataCountInserted > 0) {
        retValue = true;
        // switch on for more output during histogram creation "
        /*
        std::cout << "End of inserting data: " << dataCountInserted
            << " data points inserted to dataCollection "
            << std:: endl;
        std::cout << "and associated with the tree if "
            << "they fit into the root box" << std::endl;
        std::cout << "(check console output for possible records "
            << "of datapoints which did not fit)" << std::endl;
        */
        }

    if (dataCountInserted == 0) {
        throw HistException("No data inserted");

    }

    return retValue;
}

// check if we need to make a paving for the histogram object
// make it if we need to, matching the dimensions of the data
// as given in function argument
// return true if needed to make the paving
bool AdaptiveHistogramValidation::haveMadePaving(const RVecData& theData,
                                    const size_t dim)
{

    bool retValue = false;

    try {

        // check if we need to make the paving on the basis of the data
        if (isEmpty(rootVpaving)) {

            rootBox = makeBox(theData, dim);

            // point rootVpaving to a new SPSVnode with box myBox
            // and also pass in the not value of holdAllStats which controls
            // whether all available statistics are maintained in the
            // rootVpaving (true) or just counts (false)
            rootVpaving = new SPSVnode(rootBox, !holdAllStats);
            creationString = rootVpaving->getNodeName();

            retValue = true;
        }
    }

    catch (bad_alloc& e)
    {
        const char* msg = e.what();
        std::cerr << msg << std::endl;
        std::cerr << "Error allocating memory in "
            << "AdaptiveHistogramValidation::haveMadePaving()"
            << std::endl;
        throw;
    }

    return retValue;
    // end of making the subpaving if there was not one
}

// make a box to fit all the data
ivector AdaptiveHistogramValidation::makeBox(const RVecData& theData, const size_t dim)
{
    // set up a vector of maxes
    vector<real> maxs;

    // give maxs starting values from the first element in the rvectors
    rvector first = *theData.begin();

    for (size_t i = 1; i <=dim; i++) {
        maxs.push_back(first[i]);
    }

    // make mins the same as maxes to start with
    vector<real> mins = maxs;

    RVecDataCItr cit;

    // go over the rest of the container
    for(cit = theData.begin()+1; cit < theData.end(); cit++) {
        for (size_t i = 1; i <= dim; i++) {
            real r = (*cit)[i];
            // vectors indexed 0 - n-1, rvectors ndexed 1 - n
            if(r < mins[i-1]) {
                mins[i-1] = r;
            }
            if(r > maxs[i-1]) {
                maxs[i-1] = r;
            }
        } // end going through rvector elements
    } // end going through rvectors

    ivector retVal(dim);    // set up an ivector to become the return value

    // and make each interval the (min, max) of the corresponding elements
    // of the rvectors -/+ some padding

    std::cout << "A box is being made for the data.  "
        << "The box is " << std::endl;  // standard output message

    // make intervals and make them elements of the ivector
    for (size_t i = 1; i <=dim; i++) {
        interval myInterval(mins[i-1]-padding, maxs[i-1]+padding);
        std::cout << myInterval << "  ";    // output
        retVal[i]=myInterval;
    }
    std::cout << std::endl;

    return retVal;

}

// insert data from a container
// return number of data points inserted into dataCollection
// (and for which insertion into subpaving was attempted)
// used by all the other bulk-insert methods
// creates log file of process if logging is true.
//size_t AdaptiveHistogramValidation::insertDataFromContainer(const RVecData& theData,
//                                    const SplitDecisionObj& boolTest,
//                                   LOGGING_LEVEL logging, size_t holdOutCount)
//temporarily for air traffic trajectory analysis
size_t AdaptiveHistogramValidation::insertDataFromContainer(const RVecData& theData,
                                    const SplitDecisionObj& boolTest,
                                    LOGGING_LEVEL logging, size_t holdOutCount,
												vector<size_t> & numNodes)												
{
    size_t counter = 0;    // to count the input
    bool boolVal; // to  indicate if this data is a training or validation point
    // for logging output to keep track of splits if necessary
    int i = 0;
    std::string baseFileName = "";
    std::string s = "";
    // if we are splitting as we go and logging, set up a log file
    if ((logging != NOLOG) && (boolTest() == true)) {
        baseFileName = "splitOutput";
        s = getUniqueFilename(baseFileName);
       // outputLogStart(s);
        // log the current state of the histogram
      //  outputLog(s, i);
        i++;
    }
	 
	 //int hist = 0;
	 RVecDataCItr cit;
    // feed the data to myHist
    for(cit = theData.begin(); cit < theData.end(); cit++) {
        // put it into dataCollection
        BigDataItr it = dataCollection.end();
        it = dataCollection.insert(it, *cit);
        SPSVnode* insertedInto = NULL;
        if (counter < holdOutCount) {
          // try inserting 
          boolVal = true;
          insertedInto =
                rootVpaving->insertOneFind(it,ON_PARENT, boolTest, boolVal);
        }
        else {
         boolVal = false;
			// try inserting
         insertedInto =
                rootVpaving->insertOneFind(it,ON_PARENT, boolTest, boolVal);
					 numNodes.push_back(spTotalNodes(rootVpaving));
					 
					 /* 
					 string histFileName;
					 std::ostringstream stm1;
					 stm1 << hist;
					 histFileName = "Hist";
					 histFileName += stm1.str();
					 histFileName += ".txt";
					 cout << "get histogram:" << histFileName << endl;
					 outputToTxtTabs(histFileName);			
					 hist++;	  
					 */ 		
			}

        //insertOneFind returns either NULL if no insert possible
        // or a pointer to the node the data goes to before that node
        // is split (it could be split more than once)
        if (NULL == insertedInto) { // failed to insert
            std::cout << "Failed to insert point "
                << *cit << std::endl;
            std::cout << "Root node of subpaving has box "
                << rootVpaving << std::endl;
        }
        // successful insertion, and we are splitting as we go
        else if (boolTest() == true && boolVal == false ) {
            std::string newNames = insertedInto->getChildNodeNames();
     
            if(newNames.length() > 0) { // there are new nodes
                //add the new child names if any
                creationString += newNames;

                if (logging) { // log the current state of the histogram
//                    outputLog(s, i);
                    i++;
                }				 
					 
            }
       }
        counter++;
    }
    if (counter > 0) { // data inserted
         if ((logging != NOLOG) && (boolTest() == true))  {
            // add leaf node levels string to log
            outputFile(s, getLeafLevelsString());
        }
    }
    return counter;
}

// Method to put opening line into a log file
void AdaptiveHistogramValidation::outputLogStart(const std::string& s) const
{
    // Make a string with filename and timestamp to start log file
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    ofstream os(s.c_str());         // replace data
    if (os.is_open()) {
        os << "File " << s << " created " <<  asctime (timeinfo) << std::endl;
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

bool AdaptiveHistogramValidation::checkNodeCountForSplit(const SPSVnode * const spn,
                bool volChecking, double minVol, size_t minChildPoints)
{
    bool retValue = false;

    if ((!volChecking || (volChecking && (spn->nodeVolume() >= minVol)))
        && ((minChildPoints == 0)
            || (minChildPoints > 0
                &&
                ((spn->getCounter() >= minChildPoints) &&
                    ((spn->getMinChildCountIfSplit() == 0)
                    ||
                    (spn->getMinChildCountIfSplit() >= minChildPoints))
                ))
            )
        ) { retValue = true; }

    return retValue;
}

// ----------- histogram public methods

// default constructor
// holdAllStats defaults to false.
AdaptiveHistogramValidation::AdaptiveHistogramValidation()
        : holdAllStats(false), creationString("")
{
    rootVpaving = NULL;
    rootBox = ivector();    // ivector with length 1 and undefined elements


    // nothing happens to dataCollection when object is constructed
}

// initialised constructor with bool to control whether all stats maintained
// in root paving
AdaptiveHistogramValidation::AdaptiveHistogramValidation(bool as)
        : holdAllStats(as), creationString("")
{
    rootVpaving = NULL;
    rootBox = ivector();    // ivector with length 1 and undefined elements


    // nothing happens to dataCollection when object is constructed
}

// initialised constructor, initialised with ivector for box
// and with bool to control whether all stats are maintained in root paving.
// (defaults to false which means that only counts are maintained in rootVpaving)
AdaptiveHistogramValidation::AdaptiveHistogramValidation(ivector& v, bool as)
        : holdAllStats(as)
{
    try {
        rootVpaving = new SPSVnode(v, !as);
        creationString = rootVpaving->getNodeName();

        rootBox = v;
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor" << std::endl;
        throw HistException("Memory allocation error in constructor: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cerr << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Memory allocation error in constructor:" + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cerr << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Memory allocation error in constructor:" + msg);
    }

    // nothing happens to dataCollection when object is constructed
}

// copy constructor
AdaptiveHistogramValidation::AdaptiveHistogramValidation(const AdaptiveHistogramValidation& other)
        : rootBox(other.rootBox), holdAllStats(other.holdAllStats) 
{
    try {
       // cout << "calling copy constructor" << endl;
		  rootVpaving = new SPSVnode(*(other.rootVpaving));
        creationString = rootVpaving->getNodeName();
        creationString += rootVpaving->getChildNodeNames();

        //copy dataCollection from other to this
        dataCollection = other.dataCollection;

    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor: original error "
                                            << msg << std:: endl;
        throw HistException("Memory allocation error in constructor: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cerr << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in constructor: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cerr << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor: " + msg);
    }
}


//copy assignment operator
//deep copy of the whole histogram
AdaptiveHistogramValidation&
            AdaptiveHistogramValidation::operator=(const AdaptiveHistogramValidation& rhs)
{
    try {
        //cout << "copy assignment operator" << endl;
        // we have to make sure we delete the current paving
        if (NULL != rootVpaving) {
            delete rootVpaving;
            rootVpaving = NULL;
        }

        if (NULL != rhs.rootVpaving) {
            rootVpaving = new SPSVnode(*(rhs.rootVpaving));
            creationString = rootVpaving->getNodeName();
            creationString += rootVpaving->getChildNodeNames();

            //copy dataCollection from other to this
            dataCollection = rhs.dataCollection;
            holdAllStats = rhs.holdAllStats;
        }
        rootBox = rhs.rootBox;
        return *this;
    }
    catch (bad_alloc& ba) {
        std::cerr << "Error allocating memory in constructor" << std::endl;
        throw HistException("Error in constructor");
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cerr << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor");
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cerr << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor");
    }
}

//Destructor
AdaptiveHistogramValidation::~AdaptiveHistogramValidation()
{
    delete rootVpaving;
}

//src_trunk_0701
// get whether this has a subpaving.
bool AdaptiveHistogramValidation::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

//src_trunk_0701
int AdaptiveHistogramValidation::getLabel() const
{
	return 0; //this is temporarily for gat41 src
	//return label;
	}

// Return a pointer to the SPSVnode this manages.
SPSVnode* AdaptiveHistogramValidation::getSubPaving() const
{return rootVpaving;}


// Gets count in the rootVpaving in the root paving.
size_t AdaptiveHistogramValidation::getRootVcounter() const
{ return rootVpaving->getVcounter(); }

// Gets number of leaf nodes in the root paving.
size_t AdaptiveHistogramValidation::getRootLeaves() const
{ return spLeaves(rootVpaving); }

// Gets the sum of leaf count over volume in root paving.
real AdaptiveHistogramValidation::getRootSumLeafCountOverVol() const
{ return rootVpaving->getSumLeafCountOverVol(); }

// get the value of the minimum volume for a splittable node.
double AdaptiveHistogramValidation::getMinVol(double minVolB) const
{
    double retValue = 0.0;

    if (NULL != rootVpaving) {

        size_t counter = rootVpaving->getCounter();
        retValue =  minVolB * log(1.0*counter)*log(1.0*counter)/counter;
    }
    return retValue;
}

// get the value of holdAllStats field.
bool AdaptiveHistogramValidation::getHoldAllStats() const
{
    return holdAllStats;
}

// Get a string of the leaf node levels.
std::string AdaptiveHistogramValidation::getLeafLevelsString() const
{
    string retValue = "";
    if (NULL != rootVpaving)
        retValue = rootVpaving->getLeafNodeLevelsString();

    return retValue;
}

// method to insert rvectors from a txt file
bool AdaptiveHistogramValidation::insertRvectorsFromTxt(const std::string& s,
                                std::vector<size_t> & numNodes,
										  const SplitDecisionObj& boolTest,										  
								        const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
    bool retValue = false;

    try {
        RVecData myDataRvectors; // container for the rvectors we take in

        // try to read in the file
        retValue = readRvectorsFromTxt(myDataRvectors, s, headerlines);

        if (retValue) {
            retValue = completeDataInsertionFromVec(myDataRvectors,
                                                    boolTest, logging, 0,
																	 numNodes);
        }
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        cout << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}


// method to insert all rvectors from an RVecData object
bool AdaptiveHistogramValidation::insertFromRVec(const RVecData& rvec,
                            const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging)
{
    bool retValue = false;

    try {

        RVecData myDataRvectors; // container for the rvectors we take in

        size_t numberFound = 0;

        if (rvec.empty()) { // no data points to get
            throw HistException("No data to insert");
        }

        else { // there is data to get

            // get data from the container and check how many data points found
            size_t numberFound = getRvectorsFromRVec(myDataRvectors, rvec);


            if (numberFound > 0) {
                /*
                // confirm the amount of data taken from the container
                std::cout << "End of taking data from container of rvectors: "
                    << numberFound << " data points found" << std::endl;
                */
                // complete the data insertion
					 vector<size_t> temp;
                retValue = completeDataInsertionFromVec(myDataRvectors,
                                                        boolTest, logging, 0,
																		  temp);
            }
        }
    }
    catch (bad_alloc& ba) {
         string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;

}

// method to insert rvectors from an RVecData object for hold-out estimation
bool AdaptiveHistogramValidation::insertFromRVecForHoldOut(const RVecData& rvec, 
const SplitDecisionObj& boolTest, int holdOutCount, LOGGING_LEVEL logging)
{
    bool retValue = false;

  //  bool boolVal=true;
    
    try {

        RVecData myDataRvectors; // container for the rvectors we take in

//        size_t numberFound = 0;

        if (rvec.empty()) { // no data points to get
            throw HistException("No data to insert");
        }

        else { // there is data to get

            // get data from the container and check how many data points found
            size_t numberFound = getRvectorsFromRVec(myDataRvectors, rvec);


            if (numberFound > 0) {
                /*
                // confirm the amount of data taken from the container
                std::cout << "End of taking data from container of rvectors: "
                    << numberFound << " data points found" << std::endl;
                */
                // complete the data insertion

               // cout << "calling complateDataInsertion..." << endl;
					vector<size_t> temp;
                retValue = completeDataInsertionFromVec(myDataRvectors,
                                     boolTest, logging, holdOutCount, temp);
            }
        }
    }
    catch (bad_alloc& ba) {
         string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
} // end of insertRVecForHoldOut 

//All RSSample are associated with the root paving for hold out estimation, no splitting. */
bool AdaptiveHistogramValidation::insertFromRSSampleForHoldOut(
                                 const RSSample& rss, int label,  
											const SplitDecisionObj& boolTest,
											int holdOutCount,
											LOGGING_LEVEL logging)
{
    bool retValue = false;

  //  bool boolVal=true;
    
    try {

        RVecData myDataRvectors; // container for the rvectors we take in

            // get data from the container and check how many data points found
            size_t numberFound = getRvectorsFromRSSample(myDataRvectors, rss, label);


            if (numberFound > 0) {
                /*
                // confirm the amount of data taken from the container
                std::cout << "End of taking data from container of rvectors: "
                    << numberFound << " data points found" << std::endl;
                */
                // complete the data insertion

               // cout << "calling complateDataInsertion..." << endl;
					vector<size_t> temp;
                retValue = completeDataInsertionFromVec(myDataRvectors,
                                     boolTest, logging, holdOutCount, temp);
            }
      }
    
    catch (bad_alloc& ba) {
         string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
} // end of insertFromRSSampleForHoldOut

// prioritySplit
// method for data splitting 
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a stopping criteria is fulfilled
// outputs to a log file if logging is true
// makes its own random number generator
bool AdaptiveHistogramValidation::prioritySplit(
                const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
				LOGGING_LEVEL logging, size_t minChildPoints, 
				double minVolB, size_t maxLeafNodes)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplit(compTest, he, logging, minChildPoints, 
								    minVolB, rgsl, maxLeafNodes);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}


// prioritySplit
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a  splitting criteria is satisfied
// outputs to a log file if logging required
bool AdaptiveHistogramValidation::prioritySplit(
              const NodeCompObjVal& compTest, const HistEvalObjVal& he,  
              LOGGING_LEVEL logging, size_t minChildPoints, 
			  double minVolB, gsl_rng * rgsl,
			  size_t maxLeafNodes)
{
     //cout << "calling prioritySplit:" << endl;
     //cout << "---split 0 " << endl;
     int n = getSubPaving()->getCounter();
	 bool cancontinue = false;
	 bool TooManyLeaves = false;
	 
    //boolean for validation data
    bool boolVal = true;
    
   // check if the root box is empty
    if (NULL == rootVpaving) {
            throw HistException("No root paving for prioritySplit");
    }
    try {       

		//============checks  for splittable nodes=============================//
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
      // a multiset for the queue (key values are not necessarily unique)
      multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
      int i=0;
      if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
      }
      // put nodes into the starting set IF they meet minVol test AND IF either
      // there are enough points in the whole node
      // and minChildCountIfSplit is 0 (ie all points go to one child)
      // or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //==================end of checks=====================================//
  		
        //=========start priority queue====================================//
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	  
		  while (bigEnough && !he(this) && !TooManyLeaves) {          
            SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSVnode* chosenLargest;
            // find if there are any more equal to largest around
            multiset<SPSVnode*, MyCompare>::iterator mit;
            pair<multiset<SPSVnode*, MyCompare>::iterator,
                multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
            equalLargest = pq.equal_range(largest); // everything that = largest
            size_t numberLargest = pq.count(largest); // number of =largest

            if (numberLargest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);
                real sum = 0.0;
                // random selection of the =largest node to chose
                for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
                    sum += 1.0/(1.0*numberLargest);
                    if (rand < sum) {
                        break;
                    }
                }
                chosenLargest = *(mit); // the chosen largest in the set
                pq.erase(mit);// take the iterator to chosen largest out of the set
            }
            else {
                chosenLargest = *(pq.rbegin ()); // the only largest
                multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }
            // split the biggest one and divide up its training and validation 
            // data
            
           ExpandWithValid(chosenLargest, boolVal);
                          
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                pq.insert(chosenLargest->getLeftChild());
            }
            if (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
                pq.insert(chosenLargest->getRightChild());
            }
            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

			//==========checks to see if need to split again=========//
            //checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
					std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
			} // end of while loop
			cout << "===========End of splitting=============" << endl;
         	
   } // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }

   return (cancontinue); 
}

bool AdaptiveHistogramValidation::prioritySplitAndEstimateWithSwitch(
                   const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						 LOGGING_LEVEL logging, size_t minChildPoints, 
						 double minVolB, bool stopCrit, 
						 FinMix& mixt, int method, size_t hist, 
						 size_t maxLeafNodes, int maxCheck, double tol, int deg,
						 AdaptiveHistogramValidation& optHist)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplitAndEstimateWithSwitch(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, stopCrit, mixt, method, hist,
											  maxLeafNodes, maxCheck, tol, deg, optHist);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}

// prioritySplitAndEstimate for uniform mixtures
// method for data splitting and hold out estimation
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a stopping criteria is fulfilled
// outputs to a log file if logging is true
// makes its own random number generator
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
                   const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						 LOGGING_LEVEL logging, size_t minChildPoints, 
						 double minVolB, bool stopCrit, 
						 AdaptiveHistogram& myPart, double weight, vector<int> holesLoc,
						 int method, size_t hist, 
						 size_t maxLeafNodes, int maxCheck,
						 int& minTheta)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplitAndEstimate(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, stopCrit, myPart, weight, holesLoc,
											  method, hist,
											  maxLeafNodes, maxCheck, minTheta);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}

// prioritySplitAndEstimate for uniform mixtures
// hold out estimation based on Devroye and Lugosi 2006
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a  splitting criteria is satisfied
// outputs to a log file if logging required
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
                     const NodeCompObjVal& compTest, const HistEvalObjVal& he,  
                     LOGGING_LEVEL logging, size_t minChildPoints, 
							double minVolB, gsl_rng * rgsl, bool stopCrit, 
							AdaptiveHistogram & myPart, double weight, 
							vector<int> holesLoc,
							int method, size_t hist,
							size_t maxLeafNodes, int maxCheck,
							int& minTheta)
{
     cout << "calling prioritySplitAndEstimate:" << endl;
     cout << "---split 0 " << endl;
     int n = getSubPaving()->getCounter();
	 
	 bool cancontinue = false;
	 bool TooManyLeaves = false;
	 
    //boolean for validation data
    bool boolVal = true;
    
    // for stopping criteria
    size_t flagStop = 0;
    int currentSmallest = 0;
    
    //set up collator to keep the histograms as splits happen
    AdaptiveHistogramVCollator coll;

	 //=======initializing containers======================================
	//set up a list for the Yatracos set 
	list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > listYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (row)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecRowYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (col)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecColYatSet;    
	//set up a vector for maximum Delta_theta vectors
	vector< vector<double> > vecMaxDeltaVec;
	//initializing the vector - to allow the delta vector to be in 
	// right order  since the first histogram does not have a 
	// Yatracos set
	//the first element in this vector will not be plotted since 
	// the first histogram is an empty set
	vector<double> theta0;
	theta0.push_back(-1*(numeric_limits<double>::infinity())); 
	//the supremum of an empty set is -Infimum 
	vecMaxDeltaVec.push_back(theta0);
	//set up a vector of the corresponding theta with the minimum 
	// distance estimates
	vector< vector<int> > vecMinDistTheta;
	// set up a vector for the infimum 
	vector<double> vecInfDelta;
	// set up a vector for the integrated absolute error for each histogram
   vector<real> vecIAE; 
   vector<real> vecIAEFull;
   real minIAE = 1000.00;
   
   vector<real> TrueDelta;
   TrueDelta.push_back(-1); 
   
   real trueDeltaCurrent = 0;
   
   // to keep the histograms
   //to remove temphist
   //vector<AdaptiveHistogramValidation> tempHist;
   //==============end of initializing containers=============================//   
   // check if the root box is empty
    if (NULL == rootVpaving) {
            throw HistException("No root paving for prioritySplit");
    }
    try {       
        // add the histogram before any split happens into the collator
        size_t agg = 0;
		  coll.addToCollationWithVal(*this, 1, agg);
		  //tempHist.push_back(*this);
		  // calculate the IAE 
		  real IAE = getUnifIAE(myPart, weight, holesLoc, 0);
		  // push back into vecIAE 
		  vecIAE.push_back(IAE);
			minIAE = (IAE < minIAE) ? IAE : minIAE;
			
			//get the IAE for the full data set
			real IAEF = getUnifIAE(myPart, weight, holesLoc, 1);
		  // push back into vecIAE 
		  vecIAEFull.push_back(IAEF);

		//============checks  for splittable nodes=============================//
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
      // a multiset for the queue (key values are not necessarily unique)
      multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
      int i=0;
      if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
      }
      // put nodes into the starting set IF they meet minVol test AND IF either
      // there are enough points in the whole node
      // and minChildCountIfSplit is 0 (ie all points go to one child)
      // or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //==================end of checks=====================================//
			
        //=========start priority queue====================================//
        size_t ch = 0; //to keep track of split number
		
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	 
		  while (bigEnough && !he(this) && !TooManyLeaves) {          
            SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSVnode* chosenLargest;
            // find if there are any more equal to largest around
            multiset<SPSVnode*, MyCompare>::iterator mit;
            pair<multiset<SPSVnode*, MyCompare>::iterator,
                multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
            equalLargest = pq.equal_range(largest); // everything that = largest
            size_t numberLargest = pq.count(largest); // number of =largest

            if (numberLargest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);
                real sum = 0.0;
                // random selection of the =largest node to chose
                for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
                    sum += 1.0/(1.0*numberLargest);
                    if (rand < sum) {
                        break;
                    }
                }
                chosenLargest = *(mit); // the chosen largest in the set
                pq.erase(mit);// take the iterator to chosen largest out of the set
            }
            else {
                chosenLargest = *(pq.rbegin ()); // the only largest
                multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }
            // split the biggest one and divide up its training and validation 
            // data
	
            cout << "--------------Split " << coll.getNumberCollated() << endl;
            ExpandWithValid(chosenLargest, boolVal);       
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                pq.insert(chosenLargest->getLeftChild());
            }
            if (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
                pq.insert(chosenLargest->getRightChild());
            }
            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

		      //==========get IAE for this histogram======================//
				//cout << "get IAE for histogram " << coll.getNumberCollated() << endl;
				real IAE = getUnifIAE(myPart, weight, holesLoc, 0);
				minIAE = (IAE < minIAE) ? IAE : minIAE;
				vecIAE.push_back(IAE); 
				
				real IAEF =getUnifIAE(myPart, weight, holesLoc, 1);
				vecIAEFull.push_back(IAEF); 
				
				//cout << "add into collator" << endl;
				// add current histogram to collation
				size_t agg = 0;
				coll.addToCollationWithVal(*this, 1, agg);
				//outputToTxtTabs("SplitHist.txt");
					
			//	cout << "get the split node" << endl;
				// first we need a pointer to the corresponding CollatorSPVnode 
				// of the SPSVnode* chosenLargest     
				CollatorSPVnode * splitCollNode;
				coll.getSplitNodePtr(splitCollNode, chosenLargest);
			//	cout << chosenLargest->getNodeName() << "\t" 
			//	<< splitCollNode->getNodeName() << endl;
				
 				//cout << "get the yat class" << endl;
				// get the Yatracos class for this collation
				coll.getYatracosClassAll(splitCollNode, vecRowYatSet,
										 vecColYatSet, listYatSet);
				
				//if want to keep track of delta max, can call
				coll.getYatracosDelta(listYatSet, vecRowYatSet, 
										vecColYatSet, vecMaxDeltaVec);
															
				// keep this histogram in a container 
				//tempHist.push_back(*this);
				
				// /*optional: output this histogram
				/*string fileName = "QueueHist";
				ostringstream stm;
				stm << ch;
				fileName += stm.str();
				fileName += ".txt";
				outputToTxtTabs(fileName);
				*/
				ch++;
				
				
				/*// get the true delta
				real trueDelta = 0;
				vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;   
				//cout << "Current Yatracos set has " << (*tempList).size() << " nodes." << endl;
				for (listIt = (vecRowYatSet).begin(); listIt < vecRowYatSet.end(); listIt++) {
					if ( !(*listIt).empty() ) {
						real trueDeltaR = getUnifTrueDelta(myPart, weight, holesLoc, (*listIt));
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						//TrueDelta.push_back(trueDelta);

					}
				}

				for (listIt = (vecColYatSet).begin(); listIt < vecColYatSet.end(); listIt++) {
					if ( !(*listIt).empty() ) {
						real trueDeltaR = getUnifTrueDelta(myPart, weight, holesLoc, (*listIt));
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						//TrueDelta.push_back(trueDelta);

					}
				}
				
				//store in container for output purposes
				TrueDelta.push_back(trueDelta);

				if ( vecRowYatSet.empty() && vecColYatSet.empty() ) 
				{ trueDelta = -1; TrueDelta.push_back(trueDelta); } 

				//check theorem 10.1
				//cout << "check theorem: " << endl;
				//cout << IAE << "\t" << minIAE << "\t" << trueDelta << "\t" << TrueDelta.size() << endl;
				//if (trueDelta >= 0) {	assert(IAE <= (3*minIAE + 4*trueDelta)); }
*/
				//stopping criteria
				if (stopCrit == true) {
					//cout << "checking stopping criteria: " << endl;
					bool toStop = coll.getMinDelta(maxCheck, vecMaxDeltaVec);
					if (toStop == true) {
						cout << "Stopping criteria met." << endl;
						break;
					} 
				}

				//==========checks to see if need to split again=========//
            //checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
					std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
			} // end of while loop
			//cout << "===========End of splitting=============" << endl;
        
			//do the merging here based on vecMinDistTheta
         
			//================Outputs to .txt files=================== 
			ofstream os;         // ofstream object
			os << scientific;  // set formatting for input to oss
			os.precision(5);
			ostringstream stm1, stm2;
			int Theta=0;
			
			cout << "get delta theta using YatEnd" << endl;
			// get delta_theta
			vector<double> vecMaxDelta;
			coll.getYatracosDeltaEnd(listYatSet, vecRowYatSet, vecColYatSet, 
												vecMaxDelta);
			for (size_t i =  0; i < vecMaxDelta.size(); i++){
				cout << "Theta: " << Theta << "\t" << vecMaxDelta[i] << endl;
				Theta++;
			}
             
             
             cout << "using the old YatDelta" << endl;
			 // get the minimum delta to get the MDE histogram
			vector< vector<double> >::iterator it1; 
			vector<double>::iterator it2;
			Theta = 0;
			size_t F = vecMaxDeltaVec.size(); 
			//cout << "size of vecMaxDeltaVec: " << F << endl;
			double minDelta = 1000;
			for (size_t i = 0; i < F; i++){
				cout << "Theta: " << Theta << "\t" << vecMaxDeltaVec[F-1][i] << endl;
				if ( vecMaxDeltaVec[F-1][i] < minDelta ) { 
					minDelta = vecMaxDeltaVec[F-1][i]; 
					minTheta = Theta; 
				} 
				Theta++;
			}

			cout << "MDE at " << minTheta << " with IAE " << vecIAE[minTheta] << endl; 
            
            //optHist = tempHist[minTheta];
			
			// output vecDeltaMaxVec into .txt 
			stm1 << hist;
			stm2 << method;
			string fileNameDelta = "UnifMethod";
			fileNameDelta += stm2.str();
			fileNameDelta += "DeltaMax";
			fileNameDelta += stm1.str();
			fileNameDelta += ".txt";  
			os.open(fileNameDelta.c_str());
			
			//if want DeltaMax at each split
			/*for (it1 = vecMaxDeltaVec.begin(); it1 < vecMaxDeltaVec.end(); it1++){ 
				for (it2 = (*it1).begin(); it2 < (*it1).end(); it2++){
					os << (*it2) << "\t";
				}
				os << "\n";
			} 
			         
			for (size_t i = 0; i < F; i++){
				os << vecMaxDeltaVec[F-1][i] << endl;
			}
			os << flush;
			os.close();
			std::cout << "DeltaMax for each theta output to " << fileNameDelta << "." << endl;
			*/
			//----------------end of output for vecDeltaMaxVec-------------
 
         //output vecIAE to .txt file
			string outputFileName;// for output file
			outputFileName = "UnifMethod";
			outputFileName += stm2.str();
			outputFileName += "IAEandTrueDelta";
			outputFileName += stm1.str();
			outputFileName += ".txt";
			os.open(outputFileName.c_str());
			for (size_t i = 0; i < vecIAE.size(); i++){
				os << vecIAE[i] << "\t" << vecIAEFull[i] << endl;
			}
			os << flush;
			os.close();
			std::cout << "IAE output to " << outputFileName << endl;
			//=================end of output for vecIAE---------------------------			
   } // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority stage split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    
   return (cancontinue);
}

// hold out estimation based on Devroye and Lugosi 2006
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a  splitting criteria is satisfied
// outputs to a log file if logging required
bool AdaptiveHistogramValidation::prioritySplitAndEstimateWithSwitch(
                     const NodeCompObjVal& compTest, const HistEvalObjVal& he,  
                     LOGGING_LEVEL logging, size_t minChildPoints, 
							double minVolB, gsl_rng * rgsl, bool stopCrit, 
							FinMix& mixt, int method, size_t hist,
							size_t maxLeafNodes, int maxCheck, double tol, int deg,
							AdaptiveHistogramValidation& optHist)
{
    //cout << "calling prioritySplitAndEstimate:" << endl;
	 bool cancontinue = false;
	 bool TooManyLeaves = false;
	 
    //boolean for validation data
    bool boolVal = true;
    
    // for stopping criteria
    size_t flagStop = 0;
    int currentSmallest = 0;
    
    int n = getSubPaving()->getCounter();
    
    //set up collator to keep the histograms as splits happen
    AdaptiveHistogramVCollator coll;
    
	 //=======initializing containers======================================
	//set up a list for the Yatracos set 
	list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > listYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (row)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecRowYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (col)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecColYatSet;    
	//set up a vector for maximum Delta_theta vectors
	vector< vector<double> > vecMaxDeltaVec;
	//initializing the vector - to allow the delta vector to be in 
	// right order  since the first histogram does not have a 
	// Yatracos set
	//the first element in this vector will not be plotted since 
	// the first histogram is an empty set
	vector<double> theta0;
	theta0.push_back(-1*(numeric_limits<double>::infinity())); 
	//the supremum of an empty set is -Infimum 
	vecMaxDeltaVec.push_back(theta0);
	//set up a vector of the corresponding theta with the minimum 
	// distance estimates
	vector< vector<int> > vecMinDistTheta;
	// set up a vector for the infimum 
	vector<double> vecInfDelta;
	// set up a vector for the integrated absolute error for each histogram
   vector<real> vecIAE; 
   vector<real> vecIAEFull;
   real minIAE = 1000.00;
   
   vector<real> TrueDelta;
   TrueDelta.push_back(-1); 

	real trueDeltaCurrent = 0;
	
   // to keep the histograms
   vector<AdaptiveHistogramValidation> tempHist;
   //==============end of initializing containers=============================//   
   // check if the root box is empty
    if (NULL == rootVpaving) {
            throw HistException("No root paving for prioritySplit");
    }
    try {       
        // add the histogram before any split happens into the collator
        size_t agg = 0;
		  coll.addToCollationWithVal(*this, 1, agg);
		  tempHist.push_back(*this);
		  // calculate the IAE 
		  real IAE = mid(getFinMixIntervalIAE(mixt, tol, deg, 0));
		  // push back into vecIAE 
		  vecIAE.push_back(IAE);
			minIAE = (IAE < minIAE) ? IAE : minIAE;
			
			//get the IAE for the full data set
			real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
		  // push back into vecIAE 
		  vecIAEFull.push_back(IAEF);

		//============checks  for splittable nodes=============================//
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
      // a multiset for the queue (key values are not necessarily unique)
      multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
      int i=0;
      if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
      }
      // put nodes into the starting set IF they meet minVol test AND IF either
      // there are enough points in the whole node
      // and minChildCountIfSplit is 0 (ie all points go to one child)
      // or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //==================end of checks=====================================//
  
        //=========start priority queue====================================//
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	  
		  while (bigEnough && !he(this) && !TooManyLeaves) {          
            SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSVnode* chosenLargest;
            // find if there are any more equal to largest around
            multiset<SPSVnode*, MyCompare>::iterator mit;
            pair<multiset<SPSVnode*, MyCompare>::iterator,
                multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
            equalLargest = pq.equal_range(largest); // everything that = largest
            size_t numberLargest = pq.count(largest); // number of =largest

            if (numberLargest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);
                real sum = 0.0;
                // random selection of the =largest node to chose
                for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
                    sum += 1.0/(1.0*numberLargest);
                    if (rand < sum) {
                        break;
                    }
                }
                chosenLargest = *(mit); // the chosen largest in the set
                pq.erase(mit);// take the iterator to chosen largest out of the set
            }
            else {
                chosenLargest = *(pq.rbegin ()); // the only largest
                multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }
            // split the biggest one and divide up its training and validation 
            // data
            ExpandWithValid(chosenLargest, boolVal);
           
            cout << "---------split " << coll.getNumberCollated() << endl;
                           
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

				// remove empty boxes AND
				// but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if ( ((chosenLargest->getLeftChild())->getCounter() > 0) 
						&& (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new left child into the multiset
                //cout << (chosenLargest->getLeftChild())->getCounter()  << endl;
                pq.insert(chosenLargest->getLeftChild());
            }

            if ( ((chosenLargest->getRightChild())->getCounter() > 0) 
						&& (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new right child into the multiset
               //cout << (chosenLargest->getRightChild())->getCounter()  << endl;
                pq.insert(chosenLargest->getRightChild());
            }

            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

		      //==========get IAE for this histogram======================//
		     
				//cout << "get IAE for histogram " << coll.getNumberCollated() << endl;
				real IAE = mid(getFinMixIntervalIAE(mixt, tol, deg, 0));
				minIAE = (IAE < minIAE) ? IAE : minIAE;
				vecIAE.push_back(IAE); 
				
				real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
				vecIAEFull.push_back(IAEF); 

				// keep this histogram in a container 
				tempHist.push_back(*this);
				
				//cout << "add into collator" << endl;
				// add current histogram to collation
				size_t agg = 0;
				coll.addToCollationWithVal(*this, 1, agg);
					
				//cout << "get the split node" << endl;
				// first we need a pointer to the corresponding CollatorSPVnode 
				// of the SPSVnode* chosenLargest     
				CollatorSPVnode * splitCollNode;
				coll.getSplitNodePtr(splitCollNode, chosenLargest);
				//cout << chosenLargest->getNodeName() << "\t" << splitCollNode->getNodeName() << endl;
				
 				//cout << "get the yat class" << endl;
				// get the Yatracos class for this collation
				coll.getYatracosClassAll(splitCollNode, vecRowYatSet,
														vecColYatSet, listYatSet);

				//cout << "get delta theta" << endl;
				// get delta_theta for each theta
				coll.getYatracosDelta(listYatSet, vecRowYatSet, vecColYatSet, 
												vecMaxDeltaVec);

				// get the true delta
				real trueDelta = 0.0;
				vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;   
				//cout << "Current Yatracos set has " << (*tempList).size() << " nodes." << endl;
				for (listIt = (vecRowYatSet).begin(); listIt < vecRowYatSet.end(); listIt++) {
					if ( !((*listIt).empty()) ) {
						interval trueDeltaI = getFinMixIntervalTrueDelta(mixt, tol, deg, (*listIt));
						real trueDeltaR = mid(trueDeltaI);
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						TrueDelta.push_back(trueDelta);
					}
				}
				
				for (listIt = (vecColYatSet).begin(); listIt < vecColYatSet.end(); listIt++) {
					if ( !((*listIt).empty()) ) {
						interval trueDeltaI = getFinMixIntervalTrueDelta(mixt, tol, deg, (*listIt));
						real trueDeltaR = mid(trueDeltaI);
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						TrueDelta.push_back(trueDelta);
					}
				}
				
				
				if ( vecRowYatSet.empty() && vecColYatSet.empty() ) 
				{ trueDelta = -1; TrueDelta.push_back(trueDelta); } 
				
			//check theorem 10.1
			//if ( trueDelta >= 0) {	assert(IAE <= (3*minIAE + 4*trueDelta)); }

				//stopping criteria
				if (stopCrit == true) {
					//cout << "checking stopping criteria: " << endl;
					bool toStop = coll.getMinDelta(maxCheck, vecMaxDeltaVec);
					if (toStop == true) {
						cout << "Stopping criteria met." << endl;
						break;
					} 
				}

				//==========checks to see if need to split again=========//
            //checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
					std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
			} // end of while loop
			//cout << "===========End of splitting=============" << endl;
        
			//do the merging here based on vecMinDistTheta
         
			//================Outputs to .txt files=================== 
			ofstream os;         // ofstream object
			os << scientific;  // set formatting for input to oss
			os.precision(5);

			 // get the minimum delta to get the MDE histogram
			vector< vector<double> >::iterator it1; 
			vector<double>::iterator it2;
			int Theta=0;
			//cout << "MaxDelta" << endl;
			size_t F = vecMaxDeltaVec.size(); 
			double minDelta = 1000;
			int minTheta = 0;
			for (size_t i = 0; i < F; i++){
				//cout << "Theta: " << Theta << "\t" << vecMaxDeltaVec[F-1][i] << endl;
				if ( vecMaxDeltaVec[F-1][i] < minDelta ) { 
					minDelta = vecMaxDeltaVec[F-1][i]; 
					minTheta = Theta; 
				} 
				Theta++;
			}

			cout << "MDE at " << minTheta << " with IAE " << vecIAE[minTheta] << endl; 
         optHist = tempHist[minTheta];

			// output vecDeltaMaxVec into .txt 
			ostringstream stm1, stm2;
			stm1 << hist;
			stm2 << method;
			string fileNameDelta = "FinMixMethod";
			fileNameDelta += stm2.str();
			fileNameDelta += "DeltaMax";
			fileNameDelta += stm1.str();
			fileNameDelta += ".txt";  
			os.open(fileNameDelta.c_str());
			for (it1 = vecMaxDeltaVec.begin(); it1 < vecMaxDeltaVec.end(); it1++){ 
				for (it2 = (*it1).begin(); it2 < (*it1).end(); it2++){
					os << (*it2) << "\t";
				}
				os << "\n";
			}          
			os << flush;
			os.close();
			std::cout << "DeltaMax for each theta output to " << fileNameDelta << "." << endl;
			//----------------end of output for vecDeltaMaxVec-------------
 
         //output vecIAE to .txt file
			string outputFileName;// for output file
			outputFileName = "FinMixMethod";
			outputFileName += stm2.str();
			outputFileName += "IAEandTrueDelta";
			outputFileName += stm1.str();
			outputFileName += ".txt";
			os.open(outputFileName.c_str());
			for (size_t i = 0; i < vecIAE.size(); i++){
				os << vecIAE[i] << "\t" << vecIAEFull[i] << "\t" << TrueDelta[i] << endl;
			}
			os << flush;
			os.close();
			std::cout << "IAE output to " << outputFileName << endl;
			//=================end of output for vecIAE---------------------------			
   } // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority stage split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    
   return (cancontinue);
}


bool AdaptiveHistogramValidation::prioritySplitAndEstimateWithSwitch(
                   const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						 LOGGING_LEVEL logging, size_t minChildPoints, 
						 double minVolB, bool stopCrit, 
						AdaptiveHistogram& myPart, double weight,	std::vector<int> holesLoc,
						 int method, size_t hist, 
						 size_t maxLeafNodes, int maxCheck,
						 AdaptiveHistogramValidation& optHist)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplitAndEstimateWithSwitch(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, stopCrit, myPart, weight, holesLoc,
											  method, hist,
											  maxLeafNodes, maxCheck, optHist);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}


// hold out estimation based on Devroye and Lugosi 2006
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a  splitting criteria is satisfied
// outputs to a log file if logging required
bool AdaptiveHistogramValidation::prioritySplitAndEstimateWithSwitch(
                     const NodeCompObjVal& compTest, const HistEvalObjVal& he,  
                     LOGGING_LEVEL logging, size_t minChildPoints, 
							double minVolB, gsl_rng * rgsl, bool stopCrit, 
							AdaptiveHistogram& myPart, double weight,
							std::vector<int> holesLoc, int method, size_t hist,
							size_t maxLeafNodes, int maxCheck, 
							AdaptiveHistogramValidation& optHist)
{
    //cout << "calling prioritySplitAndEstimate:" << endl;
	 bool cancontinue = false;
	 bool TooManyLeaves = false;
	 
    //boolean for validation data
    bool boolVal = true;
    
    // for stopping criteria
    size_t flagStop = 0;
    int currentSmallest = 0;
    
    int n = getSubPaving()->getCounter();
    
    //set up collator to keep the histograms as splits happen
    AdaptiveHistogramVCollator coll;
    
	 //=======initializing containers======================================
	//set up a list for the Yatracos set 
	list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > listYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (row)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecRowYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (col)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecColYatSet;    
	//set up a vector for maximum Delta_theta vectors
	vector< vector<double> > vecMaxDeltaVec;
	//initializing the vector - to allow the delta vector to be in 
	// right order  since the first histogram does not have a 
	// Yatracos set
	//the first element in this vector will not be plotted since 
	// the first histogram is an empty set
	vector<double> theta0;
	theta0.push_back(-1*(numeric_limits<double>::infinity())); 
	//the supremum of an empty set is -Infimum 
	vecMaxDeltaVec.push_back(theta0);
	//set up a vector of the corresponding theta with the minimum 
	// distance estimates
	vector< vector<int> > vecMinDistTheta;
	// set up a vector for the infimum 
	vector<double> vecInfDelta;
	// set up a vector for the integrated absolute error for each histogram
   vector<real> vecIAE; 
   vector<real> vecIAEFull;
   real minIAE = 1000.00;
   
   vector<real> TrueDelta;
   TrueDelta.push_back(-1); 

	real trueDeltaCurrent = 0;
	
   // to keep the histograms
   vector<AdaptiveHistogramValidation> tempHist;
   //==============end of initializing containers=============================//   
   // check if the root box is empty
    if (NULL == rootVpaving) {
            throw HistException("No root paving for prioritySplit");
    }
    try {       
        // add the histogram before any split happens into the collator
        size_t agg = 0;
		  coll.addToCollationWithVal(*this, 1, agg);
		  tempHist.push_back(*this);
		  // calculate the IAE 
		  real IAE = getUnifIAE(myPart, weight, holesLoc, 0);
		  // push back into vecIAE 
		  vecIAE.push_back(IAE);
			minIAE = (IAE < minIAE) ? IAE : minIAE;
			
			//get the IAE for the full data set
			real IAEF = getUnifIAE(myPart, weight, holesLoc, 1);
		  // push back into vecIAE 
		  vecIAEFull.push_back(IAEF);

		//============checks  for splittable nodes=============================//
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
      // a multiset for the queue (key values are not necessarily unique)
      multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
      int i=0;
      if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
      }
      // put nodes into the starting set IF they meet minVol test AND IF either
      // there are enough points in the whole node
      // and minChildCountIfSplit is 0 (ie all points go to one child)
      // or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //==================end of checks=====================================//
  
        //=========start priority queue====================================//
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	  
		  while (bigEnough && !he(this) && !TooManyLeaves) {          
            SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSVnode* chosenLargest;
            // find if there are any more equal to largest around
            multiset<SPSVnode*, MyCompare>::iterator mit;
            pair<multiset<SPSVnode*, MyCompare>::iterator,
                multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
            equalLargest = pq.equal_range(largest); // everything that = largest
            size_t numberLargest = pq.count(largest); // number of =largest

            if (numberLargest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);
                real sum = 0.0;
                // random selection of the =largest node to chose
                for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
                    sum += 1.0/(1.0*numberLargest);
                    if (rand < sum) {
                        break;
                    }
                }
                chosenLargest = *(mit); // the chosen largest in the set
                pq.erase(mit);// take the iterator to chosen largest out of the set
            }
            else {
                chosenLargest = *(pq.rbegin ()); // the only largest
                multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }
            // split the biggest one and divide up its training and validation 
            // data
            ExpandWithValid(chosenLargest, boolVal);
           
            cout << "---------split " << coll.getNumberCollated() << endl;
                           
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

				// remove empty boxes AND
				// but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if ( ((chosenLargest->getLeftChild())->getCounter() > 0) 
						&& (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new left child into the multiset
                //cout << (chosenLargest->getLeftChild())->getCounter()  << endl;
                pq.insert(chosenLargest->getLeftChild());
            }

            if ( ((chosenLargest->getRightChild())->getCounter() > 0) 
						&& (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new right child into the multiset
               //cout << (chosenLargest->getRightChild())->getCounter()  << endl;
                pq.insert(chosenLargest->getRightChild());
            }

            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

		      //==========get IAE for this histogram======================//
		     
				//cout << "get IAE for histogram " << coll.getNumberCollated() << endl;
				real IAE = getUnifIAE(myPart, weight, holesLoc, 0);
				minIAE = (IAE < minIAE) ? IAE : minIAE;
				vecIAE.push_back(IAE); 
				
				real IAEF = getUnifIAE(myPart, weight, holesLoc, 1);
				vecIAEFull.push_back(IAEF); 

				// keep this histogram in a container 
				tempHist.push_back(*this);
				
				//cout << "add into collator" << endl;
				// add current histogram to collation
				size_t agg = 0;
				coll.addToCollationWithVal(*this, 1, agg);
					
				//cout << "get the split node" << endl;
				// first we need a pointer to the corresponding CollatorSPVnode 
				// of the SPSVnode* chosenLargest     
				CollatorSPVnode * splitCollNode;
				coll.getSplitNodePtr(splitCollNode, chosenLargest);
				//cout << chosenLargest->getNodeName() << "\t" << splitCollNode->getNodeName() << endl;
				
 				//cout << "get the yat class" << endl;
				// get the Yatracos class for this collation
				coll.getYatracosClassAll(splitCollNode, vecRowYatSet,
														vecColYatSet, listYatSet);

				//cout << "get delta theta" << endl;
				// get delta_theta for each theta
				//coll.getYatracosDelta(listYatSet, vecRowYatSet, vecColYatSet, 
				//								vecMaxDeltaVec);

				// get the true delta
				/*real trueDelta = 0.0;
				vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;   
				//cout << "Current Yatracos set has " << (*tempList).size() << " nodes." << endl;
				for (listIt = (vecRowYatSet).begin(); listIt < vecRowYatSet.end(); listIt++) {
					if ( !(*listIt).empty() ) {
						real trueDeltaR = getUnifTrueDelta(myPart, weight, holesLoc, (*listIt));
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						TrueDelta.push_back(trueDelta);
					}
				}

				for (listIt = (vecColYatSet).begin(); listIt < vecColYatSet.end(); listIt++) {
					if ( !(*listIt).empty() ) {
						real trueDeltaR = getUnifTrueDelta(myPart, weight, holesLoc, (*listIt));
						trueDeltaR = abs(trueDeltaR);
						trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
						//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
						trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
						//cout << "delta after comparison: " << trueDelta << endl;
						trueDeltaCurrent = trueDelta;
						TrueDelta.push_back(trueDelta);
					}
				}
				
				if ( vecRowYatSet.empty() && vecColYatSet.empty() ) 
				{ trueDelta = -1; TrueDelta.push_back(trueDelta); } 
				
				//check theorem 10.1
				//if ( trueDelta >= 0 ) {	assert(IAE <= (3*minIAE + 4*trueDelta)); }

				//stopping criteria
				if (stopCrit == true) {
					//cout << "checking stopping criteria: " << endl;
					bool toStop = coll.getMinDelta(maxCheck, vecMaxDeltaVec);
					if (toStop == true) {
						cout << "Stopping criteria met." << endl;
						break;
					} 
				}*/
				


				//==========checks to see if need to split again=========//
            //checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
					std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
			} // end of while loop
			//cout << "===========End of splitting=============" << endl;
        
			//do the merging here based on vecMinDistTheta
         
			//================Outputs to .txt files=================== 
			ofstream os;         // ofstream object
			os << scientific;  // set formatting for input to oss
			os.precision(5);

		/*	 // get the minimum delta to get the MDE histogram
			vector< vector<double> >::iterator it1; 
			vector<double>::iterator it2;
			int Theta=0;
			//cout << "MaxDelta" << endl;
			size_t F = vecMaxDeltaVec.size(); 
			double minDelta = 1000;
			int minTheta = 0;
			for (size_t i = 0; i < F; i++){
				//cout << "Theta: " << Theta << "\t" << vecMaxDeltaVec[F-1][i] << endl;
				if ( vecMaxDeltaVec[F-1][i] < minDelta ) { 
					minDelta = vecMaxDeltaVec[F-1][i]; 
					minTheta = Theta; 
				} 
				Theta++;
			}

			cout << "MDE at " << minTheta << " with IAE " << vecIAE[minTheta] << endl; 
            optHist = tempHist[minTheta];

			// output vecDeltaMaxVec into .txt 
			ostringstream stm1, stm2;
			stm1 << hist;
			stm2 << method;
			string fileNameDelta = "UnifMethod";
			fileNameDelta += stm2.str();
			fileNameDelta += "DeltaMax";
			fileNameDelta += stm1.str();
			fileNameDelta += ".txt";  
			os.open(fileNameDelta.c_str());
			for (it1 = vecMaxDeltaVec.begin(); it1 < vecMaxDeltaVec.end(); it1++){ 
				for (it2 = (*it1).begin(); it2 < (*it1).end(); it2++){
					os << (*it2) << "\t";
				}
				os << "\n";
			}          
			os << flush;
			os.close();
			std::cout << "DeltaMax for each theta output to " << fileNameDelta << "." << endl;
			//----------------end of output for vecDeltaMaxVec-------------
 
         //output vecIAE to .txt file
			string outputFileName;// for output file
			outputFileName = "UnifMethod";
			outputFileName += stm2.str();
			outputFileName += "IAEandTrueDelta";
			outputFileName += stm1.str();
			outputFileName += ".txt";
			os.open(outputFileName.c_str());
			for (size_t i = 0; i < vecIAE.size(); i++){
				os << vecIAE[i] << "\t" << vecIAEFull[i] << "\t" << TrueDelta[i] << endl;
			}
			os << flush;
			os.close();
			std::cout << "IAE output to " << outputFileName << endl;
			//=================end of output for vecIAE---------------------------			
            */
   } // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority stage split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    
   return (cancontinue);
}

// prioritySplitAndEstimate for mapped functions
// method for data splitting and hold out estimation
// outputs to a log file if logging is true
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
             const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						 LOGGING_LEVEL logging, size_t minChildPoints, 
						 double minVolB, 
						 PiecewiseConstantFunction& nodeEst, 
						 size_t maxLeafNodes, bool computeIAE,
						 vector<int> sequence,
						 vector<double> & vecMaxDelta, vector<real> & vecIAE)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplitAndEstimate(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, nodeEst, 
											  maxLeafNodes, computeIAE, sequence,
											  vecMaxDelta, vecIAE);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}

//what i need to clean
//prioritySplitAndEstimate for mapped functions
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
            const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						LOGGING_LEVEL logging, size_t minChildPoints, 
						double minVolB, gsl_rng * rgsl, 
						PiecewiseConstantFunction& nodeEst, 
						size_t maxLeafNodes, bool computeIAE,
						vector<int> sequence,
						vector<double> & vecMaxDelta, vector<real> & vecIAE)
 {
		cout << "Calling prioritySplitAndEstimate..." << endl;
		int n = getSubPaving()->getCounter();
		bool cancontinue = false;
		bool TooManyLeaves = false;
		bool boolVal = true;     //boolean for validation data
	  size_t numHist = 0; //a counter to track the number of histograms
	    
		//set up collator to keep the histograms as splits happen
		AdaptiveHistogramVCollator coll;
    
		//initializing containers
		//container for getMinDistEst()
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecYatSet;

		//set up a list for the Yatracos set for ...
		list< set<CollatorSPVnode*, less<CollatorSPVnode*> > >* listYatSet 
		= new list< set<CollatorSPVnode*, less<CollatorSPVnode*> > >;
	
		//set up a vector for sets of pointers to CollatorSPVnode (row)
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >* vecRowYatSet
		 = new vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >;
	
		//set up a vector for sets of pointers to CollatorSPVnode (col)
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >* vecColYatSet
		 = new vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > >;    
	
		//set up a vector for maximum Delta_theta vectors
		//vector< vector<double> > vecMaxDeltaVec;
		//initializing the vector - to allow the delta vector to be in 
		// right order  since the first histogram does not have a 
		// Yatracos set
		//the first element in this vector will not be plotted since 
		// the first histogram is an empty set
		vector<double> theta0;
		theta0.push_back(-1*(numeric_limits<double>::infinity())); 
		//the supremum of an empty set is -Infimum 
		//vecMaxDeltaVec.push_back(theta0);
		//set up a vector of the corresponding theta with the minimum 
		// distance estimates
		//vector< vector<int> > vecMinDistTheta;
		// set up a vector for the infimum 
		vector<double> vecInfDelta;
		// set up a vector for the integrated absolute error for each histogram
		//vector<real>* vecIAE = new vector<real>; 
		vector<real> vecIAEFull;
		real minIAE = 1000.00;

   
		vector<real> TrueDelta;
		TrueDelta.push_back(-1); 
		real trueDeltaCurrent = 0;
   	//end of initializing containers//
   	   
		// check if the root box is empty
		if (NULL == rootVpaving) {
				throw HistException("No root paving for prioritySplit");
    }
    try {       
        // add the histogram before any split happens into the collator
        size_t agg = 0;
				coll.addToCollationWithVal(*this, 1, agg);
				numHist += 1;
		
				if (computeIAE == TRUE) {
					PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(*this);
					real IAE = nodeEst.getIAE(*tempPCF);
					delete tempPCF;
					(vecIAE).push_back(IAE);
				}
				
				//get the IAE for the full data set
				//	real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
				// push back into vecIAE 
				// vecIAEFull.push_back(IAEF);

				//checks for splittable nodes//
				bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
				// a multiset for the queue (key values are not necessarily unique)
				multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
				int i=0;
				if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
				}
      
				// put nodes into the starting set IF they meet minVol test AND IF either
				// there are enough points in the whole node
				// and minChildCountIfSplit is 0 (ie all points go to one child)
				// or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //end of checks//
		
				//start priority queue//
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	  
				while (bigEnough && !he(this) && !TooManyLeaves) {          
					SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
					SPSVnode* chosenLargest;
					// find if there are any more equal to largest around
					multiset<SPSVnode*, MyCompare>::iterator mit;
					pair<multiset<SPSVnode*, MyCompare>::iterator,
							multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
					equalLargest = pq.equal_range(largest); // everything that = largest
					size_t numberLargest = pq.count(largest); // number of =largest

					if (numberLargest > 1) {
							// draw a random number in [0,1)
							double rand = gsl_rng_uniform(rgsl);
							real sum = 0.0;
							// random selection of the =largest node to chose
							for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
									sum += 1.0/(1.0*numberLargest);
									if (rand < sum) {
											break;
									}
							}
							chosenLargest = *(mit); // the chosen largest in the set
							pq.erase(mit);// take the iterator to chosen largest out of the set
           }
           else {
							chosenLargest = *(pq.rbegin ()); // the only largest
							multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
							it--;
							pq.erase(it);// take this largest out of the set
            }
            
            // split the biggest one and divide up its training and validation 
            // data
            ExpandWithValid(chosenLargest, boolVal);
                          
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
            // the child's minChildCountIfSplit is 0 (ie all points go to
            // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                pq.insert(chosenLargest->getLeftChild());
            }
            if (checkNodeCountForSplit(chosenLargest->getRightChild(),
                volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
                pq.insert(chosenLargest->getRightChild());
            }
            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

						//get IAE for this histogram if computeIAE == TRUE
						if (computeIAE == TRUE) {
							PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(*this);
							// gotta double check this (that it is only using the height 
							// from counter and not from Vcounter
							real IAE = nodeEst.getIAE(*tempPCF);
							(vecIAE).push_back(IAE); 
							delete tempPCF;
	          	//real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
							//vecIAEFull.push_back(IAEF); 
						}
						numHist += 1;
				
						// only collate the k-th histogram and obtain the delta values
						if (find(sequence.begin(), sequence.end(), numHist) != sequence.end()) {
							coll.addToCollationWithVal(*this, 1, agg);
							cout << "---- Hist " << numHist << "-----" << endl;
						}
					
						//cout << "get the split node" << endl;
						// first we need a pointer to the corresponding CollatorSPVnode 
						// of the SPSVnode* chosenLargest     
						//CollatorSPVnode * splitCollNode;
						//coll.getSplitNodePtr(splitCollNode, chosenLargest);
						//cout << chosenLargest->getNodeName() << "\t" << splitCollNode->getNodeName() << endl;
						
		 				//cout << "get the yat class" << endl;
						// get the Yatracos class for this collation
						//coll.getYatracosClassAll(splitCollNode, *vecRowYatSet,
						//										*vecColYatSet, 
						//										*listYatSet);
						//} //end of every k-th histogram
						
						//cout << "get delta theta" << endl;
						// get delta_theta for each theta
						//coll.getYatracosDelta(listYatSet, vecRowYatSet, vecColYatSet, 
							//							vecMaxDeltaVec);
		
						/*later...
						// get the true delta
						real trueDelta = 0.0;
						vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;   
						//cout << "Current Yatracos set has " << (*tempList).size() << " nodes." << endl;
						for (listIt = (vecRowYatSet).begin(); listIt < vecRowYatSet.end(); listIt++) {
							if ( !(*listIt).empty() ) {
									real trueDeltaR = getMappedFunctionTrueDelta(nodeEst, (*listIt));
									trueDeltaR = abs(trueDeltaR);
									trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
									//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
									trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
									//cout << "delta after comparison: " << trueDelta << endl;
									trueDeltaCurrent = trueDelta;
									//TrueDelta.push_back(trueDelta);
								}
							}
		
						for (listIt = (vecColYatSet).begin(); listIt < vecColYatSet.end(); listIt++) {
							if ( !(*listIt).empty() ) {
								real trueDeltaR = getMappedFunctionTrueDelta(nodeEst, (*listIt));
								trueDeltaR = abs(trueDeltaR);
								trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
								//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
								trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
								//cout << "delta after comparison: " << trueDelta << endl;
								trueDeltaCurrent = trueDelta;
								//TrueDelta.push_back(trueDelta);
							}
						}
										
						if ( vecRowYatSet.empty() && vecColYatSet.empty() ) 
						{ trueDelta = -1; TrueDelta.push_back(trueDelta); } 
						
						TrueDelta.push_back(trueDelta);
						
						//check theorem 10.1
						//cout << "check theorem: " << endl;
						//cout << IAE << "\t" << minIAE << "\t" << trueDelta << endl;
						//if ( trueDelta >= 0) {	assert(IAE <= (3*minIAE + 4*trueDelta)); }
						*/
						
						//checks to see if need to split again
						//checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
							std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
						// check if number of leaf nodes in subpaving > maxLeafNodes
						// maximum number of leaf nodes allowed
						//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
						TooManyLeaves = (getRootLeaves() > maxLeafNodes);
						if ( TooManyLeaves) {
							std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
						}
					} // end of while loop
				
					//get the Delta values
					coll.getMinDistEst(vecMaxDelta, vecYatSet);				
					
					//Outputs to .txt files
					//ofstream os;         // ofstream object
					//os << scientific;  // set formatting for input to oss
					//os.precision(5);
          // int Theta=0;
          //cout << "get delta theta using YatEnd" << endl;
					// get delta_theta
					//vector<double>* vecMaxDelta = new vector<double>;
			
					//coll.getYatracosDeltaEnd(*listYatSet, *vecRowYatSet, 
					//						 *vecColYatSet, *vecMaxDelta);
			
					//int F = (vecMaxDelta).size();
					//cout << F << endl;
					//double minDelta = 1000;
					//int minTheta = 0;
					/*
					Theta = 0;
					for (size_t i = 0; i < F; i++){
					//	cout << "Theta: " << i << "\t" << vecMaxDelta[i] << endl;
						if ( vecMaxDelta[i] < minDelta ) { 
							minDelta = vecMaxDelta[i]; 
							minTheta = i; 
						} 
					}
					 cout << "minDelta = " << minDelta << " giving MDE at " 
						<< minTheta << endl; 
					*/
				
					/* // get the minimum delta to get the MDE histogram
					vector< vector<double> >::iterator it1; 
					vector<double>::iterator it2;
		
					//cout << "MaxDelta" << endl;
					size_t F = vecMaxDeltaVec.size(); 
					double minDelta = 1000;
					//int minTheta = 0;
					Theta = 0;
					for (size_t i = 0; i < F; i++){
						cout << "Theta: " << Theta << "\t" << vecMaxDeltaVec[F-1][i] << endl;
						if ( vecMaxDeltaVec[F-1][i] < minDelta ) { 
							minDelta = vecMaxDeltaVec[F-1][i]; 
							minTheta = Theta; 
						} 
						Theta++;
					}
	
				  cout << "minDelta = " << minDelta << " giving MDE at " 
						<< minTheta << endl; 
		          *//*
					// output vecDeltaMaxVec into .txt 
					ostringstream stm1, stm2;
					//stm1 << hist;
					stm2 << method;
					string fileNameDelta = "Mapped";
					//fileNameDelta += stm2.str();
					fileNameDelta += "DeltaMax";
					fileNameDelta += stm2.str();
					fileNameDelta += ".txt";  
					os.open(fileNameDelta.c_str());
					for (it1 = vecMaxDeltaVec.begin(); it1 < vecMaxDeltaVec.end(); it1++){ 
						for (it2 = (*it1).begin(); it2 < (*it1).end(); it2++){
							os << (*it2) << "\t";
						}
						os << "\n";
					}*/         
				
					//comment this out later
					/*cout << "output to txt" << endl;
					for (size_t i = 0; i < F; i++){
						os << (vecMaxDelta)[i] << endl;
					}
					os << flush;
					os.close();
					
					fileNameDelta = "MappedIAE";
					fileNameDelta += stm2.str();
					fileNameDelta += ".txt";  
					os.open(fileNameDelta.c_str());
					for (size_t i = 0; i < (vecIAE).size(); i++){
						os << (vecIAE)[i] << endl;
					}			 
					os << flush;
					os.close();
					std::cout << "Files written." << endl;
					*/
					delete listYatSet, vecRowYatSet, vecColYatSet;
							
		      /* //output vecIAE to .txt file
					string outputFileName;// for output file
					outputFileName = "Mapped";
					//outputFileName += stm2.str();
					outputFileName += "IAEandTrueDelta";
					outputFileName += stm2.str();
					outputFileName += ".txt";
					os.open(outputFileName.c_str());
					for (size_t i = 0; i < vecIAE.size(); i++){
						//os << vecIAE[i] << "\t" << vecIAEFull[i] << TrueDelta[i] << endl;
						os << vecIAE[i] << endl;
					}
					os << flush;
					os.close();
					std::cout << "IAE output to " << outputFileName << endl;*/

		} // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority stage split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    
   return (cancontinue);
}



// prioritySplitAndEstimate for finite mixtures
// method for data splitting and hold out estimation
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a stopping criteria is fulfilled
// outputs to a log file if logging is true
// makes its own random number generator
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
                   const NodeCompObjVal& compTest, const HistEvalObjVal& he, 
						 LOGGING_LEVEL logging, size_t minChildPoints, 
						 double minVolB, bool stopCrit, 
						 FinMix& mixt, int method, size_t hist, 
						 size_t maxLeafNodes, int maxCheck, double tol, int deg,
						 int& minTheta)
						 //AdaptiveHistogramValidation& optHist)
{
    gsl_rng * rgsl = NULL;
    bool cancontinue;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        cancontinue = prioritySplitAndEstimate(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, stopCrit, mixt, method, hist,
											  maxLeafNodes, maxCheck, tol, deg, minTheta);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return cancontinue;
}

// prioritySplitAndEstimate for finite mixtures
// hold out estimation based on Devroye and Lugosi 2006
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a  splitting criteria is satisfied
// outputs to a log file if logging required
bool AdaptiveHistogramValidation::prioritySplitAndEstimate(
                     const NodeCompObjVal& compTest, const HistEvalObjVal& he,  
                     LOGGING_LEVEL logging, size_t minChildPoints, 
							double minVolB, gsl_rng * rgsl, bool stopCrit, 
							FinMix& mixt, int method, size_t hist,
							size_t maxLeafNodes, int maxCheck, double tol, int deg,
							int& minTheta)
							//AdaptiveHistogramValidation& optHist)
{
    //cout << "calling prioritySplitAndEstimate:" << endl;
	 bool cancontinue = false;
	 bool TooManyLeaves = false;
	 
    //boolean for validation data
    bool boolVal = true;
    
    // for stopping criteria
    size_t flagStop = 0;
    int currentSmallest = 0;
    
    int n = getSubPaving()->getCounter();
    
    //set up collator to keep the histograms as splits happen
    AdaptiveHistogramVCollator coll;
    
	 //=======initializing containers======================================
	//set up a list for the Yatracos set 
	list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > listYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (row)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecRowYatSet;
	//set up a vector for sets of pointers to CollatorSPVnode (col)
	vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > vecColYatSet;    
	//set up a vector for maximum Delta_theta vectors
	vector< vector<double> > vecMaxDeltaVec;
	//initializing the vector - to allow the delta vector to be in 
	// right order  since the first histogram does not have a 
	// Yatracos set
	//the first element in this vector will not be plotted since 
	// the first histogram is an empty set
	vector<double> theta0;
	theta0.push_back(-1*(numeric_limits<double>::infinity())); 
	//the supremum of an empty set is -Infimum 
	vecMaxDeltaVec.push_back(theta0);
	//set up a vector of the corresponding theta with the minimum 
	// distance estimates
	vector< vector<int> > vecMinDistTheta;
	// set up a vector for the infimum 
	vector<double> vecInfDelta;
	// set up a vector for the integrated absolute error for each histogram
   vector<real> vecIAE; 
   vector<real> vecIAEFull;
   real minIAE = 1000.00;
   
   vector<real> TrueDelta;
   TrueDelta.push_back(-1); 
   
   real trueDeltaCurrent = 0;
   
   // to keep the histograms
   vector<AdaptiveHistogramValidation> tempHist;
   //==============end of initializing containers=============================//   
   // check if the root box is empty
    if (NULL == rootVpaving) {
            throw HistException("No root paving for prioritySplit");
    }
    try {       
        // add the histogram before any split happens into the collator
        size_t agg = 0;
		  coll.addToCollationWithVal(*this, 1, agg);
		  //tempHist.push_back(*this);
		  // calculate the IAE 
		  real IAE = mid(getFinMixIntervalIAE(mixt, tol, deg, 0));
		  // push back into vecIAE 
		  vecIAE.push_back(IAE);
			minIAE = (IAE < minIAE) ? IAE : minIAE;
			
			//get the IAE for the full data set
			real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
		  // push back into vecIAE 
		  vecIAEFull.push_back(IAEF);

		//============checks  for splittable nodes=============================//
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        //logging
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqOutput";
            s = getUniqueFilename(baseFileName);
        }
        // make volChecking true if minVolB is > 0.0
        if (minVolB > 0.0) {
            // minimum volume of a splittable node is minVolB(log n)^2/n
            minVol = getMinVol(minVolB);
            volChecking = true;
        }
      // a multiset for the queue (key values are not necessarily unique)
      multiset<SPSVnode*, MyCompare> pq((MyCompare(compTest)));
      int i=0;
      if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);    
            i++;
      }
      // put nodes into the starting set IF they meet minVol test AND IF either
      // there are enough points in the whole node
      // and minChildCountIfSplit is 0 (ie all points go to one child)
      // or the minChildCountIfSplit test passed
        if (rootVpaving->isLeaf()) {
            // check to insert a copy of the rootVpaving pointer into the set
           if (checkNodeCountForSplit(rootVpaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootVpaving);
            }
        }
        else { // root is not a leaf
            SPSVnodePtrs leaves;
            rootVpaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSVnodePtrsItr sit;            
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
						   pq.insert(*sit);
                }
            }
        }
        cancontinue = (!pq.empty());
        bool bigEnough = cancontinue;
        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }        
        //==================end of checks=====================================//
  
			size_t ch = 0;
			
        //=========start priority queue====================================//
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out	  
		  while (bigEnough && !he(this) && !TooManyLeaves) {          
            SPSVnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSVnode* chosenLargest;
            // find if there are any more equal to largest around
            multiset<SPSVnode*, MyCompare>::iterator mit;
            pair<multiset<SPSVnode*, MyCompare>::iterator,
                multiset<SPSVnode*, MyCompare>::iterator> equalLargest;
            equalLargest = pq.equal_range(largest); // everything that = largest
            size_t numberLargest = pq.count(largest); // number of =largest

            if (numberLargest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);
                real sum = 0.0;
                // random selection of the =largest node to chose
                for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
                    sum += 1.0/(1.0*numberLargest);
                    if (rand < sum) {
                        break;
                    }
                }
                chosenLargest = *(mit); // the chosen largest in the set
                pq.erase(mit);// take the iterator to chosen largest out of the set
            }
            else {
                chosenLargest = *(pq.rbegin ()); // the only largest
                multiset<SPSVnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }
            // split the biggest one and divide up its training and validation 
            // data
            
            cout << "--------------Split " << coll.getNumberCollated() << endl;
            ExpandWithValid(chosenLargest, boolVal);
                          
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            if (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                pq.insert(chosenLargest->getLeftChild());
            }
            if (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
                pq.insert(chosenLargest->getRightChild());
            }
            if (logging != NOLOG) {
                // To add current state of histogram to log file                   
                i++;
            }

		      //==========get IAE for this histogram======================//
				//cout << "get IAE for histogram " << coll.getNumberCollated() << endl;
				real IAE = mid(getFinMixIntervalIAE(mixt, tol, deg, 0));
				minIAE = (IAE < minIAE) ? IAE : minIAE;
				vecIAE.push_back(IAE); 
				
				real IAEF = mid(getFinMixIntervalIAE(mixt, tol, deg, 1));
				vecIAEFull.push_back(IAEF); 
							
				// keep this histogram in a container 
				//tempHist.push_back(*this);
				
				/*
				string fileName = "QueueHist";
				ostringstream stm;
				stm << ch;
				fileName += stm.str();
				fileName += ".txt";
				outputToTxtTabs(fileName);
				
				ch++;
				*/
				//cout << "add into collator" << endl;
				// add current histogram to collation
				size_t agg = 0;
				coll.addToCollationWithVal(*this, 1, agg);
					
				//cout << "get the split node" << endl;
				// first we need a pointer to the corresponding CollatorSPVnode 
				// of the SPSVnode* chosenLargest     
				CollatorSPVnode * splitCollNode;
				coll.getSplitNodePtr(splitCollNode, chosenLargest);
				//cout << chosenLargest->getNodeName() << "\t" << splitCollNode->getNodeName() << endl;
				
 				//cout << "get the yat class" << endl;
				// get the Yatracos class for this collation
				coll.getYatracosClassAll(splitCollNode, vecRowYatSet,
														vecColYatSet, listYatSet);

				//cout << "get delta theta" << endl;
				// get delta_theta for each theta
				coll.getYatracosDelta(listYatSet, vecRowYatSet, vecColYatSet, 
												vecMaxDeltaVec);

				// get the true delta
				real trueDelta = 0.0;
				TrueDelta.push_back(trueDelta);
			/*	vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator listIt;   
				//cout << "Current Yatracos set has " << (*tempList).size() << " nodes." << endl;
				for (listIt = (vecRowYatSet).begin(); listIt < vecRowYatSet.end(); listIt++) {
					if ( !(*listIt).empty() ) {
							interval TrueDeltaI = getFinMixIntervalTrueDelta(mixt, tol, deg, (*listIt));
							real trueDeltaR = mid(TrueDeltaI);
							trueDeltaR = abs(trueDeltaR);
							trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
							//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
							trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
							//cout << "delta after comparison: " << trueDelta << endl;
							trueDeltaCurrent = trueDelta;
							//TrueDelta.push_back(trueDelta);
						}
					}

					for (listIt = (vecColYatSet).begin(); listIt < vecColYatSet.end(); listIt++) {
						if ( !(*listIt).empty() ) {
							interval TrueDeltaI = getFinMixIntervalTrueDelta(mixt, tol, deg, (*listIt));
							real trueDeltaR = mid(TrueDeltaI);
							trueDeltaR = abs(trueDeltaR);
							trueDelta = (trueDeltaR > trueDelta) ? trueDeltaR : trueDelta;
							//cout << "previous: " << trueDeltaCurrent << "\t current: " << trueDelta << endl;
							trueDelta = (trueDeltaCurrent > trueDelta) ? trueDeltaCurrent : trueDelta;
							//cout << "delta after comparison: " << trueDelta << endl;
							trueDeltaCurrent = trueDelta;
							//TrueDelta.push_back(trueDelta);
						}
					}

				if ( vecRowYatSet.empty() && vecColYatSet.empty() ) 
				{ trueDelta = -1; TrueDelta.push_back(trueDelta); } 
				
				TrueDelta.push_back(trueDelta);
				*/
				//check theorem 10.1
				//cout << "check theorem: " << endl;
				//cout << IAE << "\t" << minIAE << "\t" << trueDelta << endl;
				//if ( trueDelta >= 0) {	assert(IAE <= (3*minIAE + 4*trueDelta)); }

				//stopping criteria
				if (stopCrit == true) {
					//cout << "checking stopping criteria: " << endl;
					bool toStop = coll.getMinDelta(maxCheck, vecMaxDeltaVec);
					if (toStop == true) {
						cout << "Stopping criteria met." << endl;
						break;
					} 
				}
				//==========checks to see if need to split again=========//
            //checking if there are any more 'largest' nodes in the priority queue
            bigEnough = (!pq.empty());
            if (!bigEnough){    
					std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
            }
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
			} // end of while loop
			//cout << "===========End of splitting=============" << endl;
        
			//do the merging here based on vecMinDistTheta
         
			//================Outputs to .txt files=================== 
			ofstream os;         // ofstream object
			os << scientific;  // set formatting for input to oss
			os.precision(5);

			 // get the minimum delta to get the MDE histogram
			vector< vector<double> >::iterator it1; 
			vector<double>::iterator it2;
			int Theta=0;
			//cout << "MaxDelta" << endl;
			size_t F = vecMaxDeltaVec.size(); 
			double minDelta = 1000;
			for (size_t i = 0; i < F; i++){
				//cout << "Theta: " << Theta << "\t" << vecMaxDeltaVec[F-1][i] << endl;
				if ( vecMaxDeltaVec[F-1][i] < minDelta ) { 
					minDelta = vecMaxDeltaVec[F-1][i]; 
					minTheta = Theta; 
				} 
				Theta++;
			}

		  cout << "minDelta = " << minDelta << " giving MDE at " 
				<< minTheta << " with IAE " << vecIAE[minTheta] << endl; 
         //optHist = tempHist[minTheta];

			// output vecDeltaMaxVec into .txt 
			ostringstream stm1, stm2;
			stm1 << hist;
			stm2 << method;
			string fileNameDelta = "FinMixMethod";
			fileNameDelta += stm2.str();
			fileNameDelta += "DeltaMax";
			fileNameDelta += stm1.str();
			fileNameDelta += ".txt";  
			os.open(fileNameDelta.c_str());
			/*for (it1 = vecMaxDeltaVec.begin(); it1 < vecMaxDeltaVec.end(); it1++){ 
				for (it2 = (*it1).begin(); it2 < (*it1).end(); it2++){
					os << (*it2) << "\t";
				}
				os << "\n";
			}*/         
			
			for (size_t i = 0; i < F; i++){
				os << vecMaxDeltaVec[F-1][i] << endl;
			}
			 
			os << flush;
			os.close();
			std::cout << "DeltaMax for each theta output to " << fileNameDelta << "." << endl;
			//----------------end of output for vecDeltaMaxVec-------------
 
         //output vecIAE to .txt file
			string outputFileName;// for output file
			outputFileName = "FinMixMethod";
			outputFileName += stm2.str();
			outputFileName += "IAEandTrueDelta";
			outputFileName += stm1.str();
			outputFileName += ".txt";
			os.open(outputFileName.c_str());
			for (size_t i = 0; i < vecIAE.size(); i++){
				os << vecIAE[i] << "\t" << vecIAEFull[i] << TrueDelta[i] << endl;
			}
			os << flush;
			os.close();
			std::cout << "IAE output to " << outputFileName << endl;
			//=================end of output for vecIAE---------------------------			
   } // end of try
    
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority stage split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    
   return (cancontinue);
}

//splits histogram according to string instruction
//returns true if some splitting was achieved
bool AdaptiveHistogramValidation::splitToShape(std::string instruction)
{
    bool success = false;

    // checks:  is there a root paving, is the string properly formed?
    if (NULL == rootVpaving) {
        throw HistException("No root paving for splitToShape");
    }

    if (instruction.length() == 0) {
		  throw HistException("No instruction");
    }

    std::string legal(", 0123456789");
    if (instruction.find_first_not_of(legal) != std::string::npos) {
        throw HistException("Illegal character in instruction");
    }

    try { // all seems to be okay, we can start spliting the root paving
        // specify what to look for as numbers or decimal point or + or -

       success = rootVpaving->splitRootToShape(instruction);

        if (success) {
            // update the creation string
            creationString = rootVpaving->getNodeName()
                + rootVpaving->getChildNodeNames();
        }
        else {
            std::cerr << std::endl;
            std::cerr << "Your instruction does not describe a proper tree.";
            std::cerr << "  Please check your instruction and try again."
            << std::endl;
       }
    }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory in splitToShape.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in splitToShape.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in splitToShape.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in splitToShape.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return success;
}

// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec AdaptiveHistogramValidation::getLeafLevels() const
{
    IntVec levels; // empty container

    if (NULL != rootVpaving) {
        rootVpaving->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// returns a vector of leaf counts
// left to right
Size_tVec AdaptiveHistogramValidation::getLeafCounts() const
{
    Size_tVec counts; // empty container
    if (NULL != rootVpaving) {
        rootVpaving->getLeafNodeCounts(counts);
        //levels has now been filled in
    }
    return counts;
}


// make a .dot file for the histogram
bool AdaptiveHistogramValidation::outputGraphDot() const
{
    bool success = false;

    if (NULL != rootVpaving) {
        success = rootVpaving->outputGraphDot();

    }
    else {
        std::cerr << "Sorry, you can't make a graph without a root paving"
                << std::endl;
    }
    return success;
}


// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void AdaptiveHistogramValidation::outputToTxtTabs(const std::string& s,
                            bool confirm) const
{
    if (NULL != rootVpaving) {

        // To generate a file output of the AdaptiveHistogramValidation object
        ofstream os(s.c_str());         // Filename, c-string version
        if (os.is_open()) {

            getSubPaving()->leavesOutputTabs(os); // the output
            if (confirm)
                std::cout << "The output of the AdaptiveHistogramValidation "
                    << "has been written to " << s << std::endl << std::endl;
        }
        else {
            std::cerr << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }
}


// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void AdaptiveHistogramValidation::outputRootToTxt(const std::string& s,
                                            bool confirm) const
{
    if (NULL != rootVpaving) {

        // To generate a file output of root node of the AdaptiveHistogramValidation
        ofstream os(s.c_str());         // Filename, c-string version
        getSubPaving()->nodePrint(os); // the output
        if (confirm)
            std::cout << "Details of the root paving of the AdaptiveHistogramValidation "
                << "has been written to " << s << std::endl << std::endl;
    }

}

/*! Get the IAE of the corresponding distribution based on function arguments.
*/
real AdaptiveHistogramValidation::getIAE(int distr)
{
      real IAE = 0;
		taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval);
		switch(distr)
		{ 
         
			case 1: //bivariate gaussian mixtures
         testpnt = BiGOP;
			IAE = get2DIAE(testpnt);
			break;
			
			case 2: // Levy 2D
			testpnt = LevyOP;
			IAE = get2DIAE(testpnt);
			break;
			
			case 3: //Rosenbrock 2D
			testpnt = RosenOP;
			IAE = get2DIAE(testpnt);
			break;
		}		
		return IAE;
}

/*! Get the IAE for the unform distribution
*/
real AdaptiveHistogramValidation::getUnifIAE(AdaptiveHistogram & myPart, double weight,
															vector<int> holesLoc, bool full)
{
   // get the true height, f of the corresponding box in myPart
	SPSnodePtrs trueLeaves;
	SPSnodePtrsItr trueIt;
	//AdaptiveHistogram * adhPtr;
	//adhPtr = &myPart;
	(myPart).getSubPaving()->getLeaves(trueLeaves);

	// setting up containers for the leaves
	SPSVnodePtrs leaves;
	SPSVnodePtrsItr it;
	getSubPaving()->getLeaves(leaves); // fill the container

	double trueF; //true density
	ivector temp;
	
	dotprecision dpIAE;    // use type dotprecision for summation  
	dpIAE=0.0;

	int n = getSubPaving()->getRootCounter();
	int allN = n + rootVpaving->getVcounter();
	//go through all the leaves in this
	for(it = leaves.begin(); it < leaves.end(); it++) {
		ivector thisBox = (*it)->getBox();
		//cout << "====checking " << (*it)->getBox() << endl;
      
		// get the height of this leaf
		double fhat;
		if ( full == 0 ) { 
			fhat = (*it)->getCounter()/(*it)->nodeVolume()/n;
			//cout << (*it)->getCounter() << "\t" << (*it)->nodeVolume() <<"\t" << fhat << endl;
		}
		else {
			size_t totalCount = (*it)->getCounter() + (*it)->getVcounter();
			fhat = totalCount/(*it)->nodeVolume()/allN;
		}

		//cout << full << "\tfhat for box " << ":" << fhat << endl;

		size_t L = 0;
		for (trueIt = trueLeaves.begin(); trueIt < trueLeaves.end(); trueIt++) {
			//cout << "----True leaf: " << (*trueIt)->getBox() << "\t" << endl;
			ivector trueBox = (*trueIt)->getBox();

			if (  holesLoc[L] == 0 ) { trueF = 0; }
			else { trueF = weight/((*trueIt)->nodeVolume()); }
			//cout << "pdf: " << trueF  << endl;
			
			// if this is contained in trueBox
			if ( (*it)->getBox() <= (*trueIt)->getBox() || (*it)->getBox() == (*trueIt)->getBox() ) {
				//use the volume of this
				real r = ((*it)->nodeVolume())*(fhat - trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
				//can move on to next leaf rather than iterating thru all trueBoxes
				//think about this later
			} //end of if this box is in trueBox
			
			// if this contains trueBox
			else if ((*trueIt)->getBox() <= (*it)->getBox()) {
				//use the volume of trueBox
				real r = ((*trueIt)->nodeVolume())*(fhat - trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
			} //end of if trueBox is in this box
			
			// if this is partially contained in trueBox 
			else if 	(Intersection(temp, thisBox, trueBox)) {
				if (Inf(temp) != Sup(temp)){
					double volume = Volume(temp);
					real r = volume*(fhat - trueF);
					//cout << "r: " << r << "\t" << abs(r) << endl;
					accumulate(dpIAE, abs(r), 1.0);
				}
			}
			L++;
		} // end of going through trueBoxes
	} // end of going through thisBoxes
	
   //cast dotprecision to real
   real unifIAE = rnd(dpIAE);
   //cout << "IAE: " << unifIAE << endl;
	return unifIAE;

}
  
  
  
  
  

/*! Get the IAE for a finite mixture distribution
*/
real AdaptiveHistogramValidation::getFinMixIAE(FinMix& mixt)
{
	//---------fill in containers for names, volumes and heights of leaves
	SPSVnodePtrs leaves; // set up empty container for leaf node pointers
	SPSVnodePtrsItr it; // and an iterator over the container
	(*this).getSubPaving()->getLeaves(leaves); // fill the container
	//a container for the counts
	IntVec counts;  // IntVec is a typedef for vector<int>
	//a container for the boxes
	vector<ivector> boxes; vector<ivector>::iterator itBoxes;
	//a container for the volumes
	vector<double> volumes;
	// a container for fhat
	vector<double> fhat; vector<double>::iterator itFhat;
	//number of points
	int n = rootVpaving->getCounter();
	          
	for(it = leaves.begin(); it < leaves.end(); it++) {
		// remember that it points to a pointer, so *it is still a ptr
		// get the counts in all the leaves
	   counts.push_back((*it)->getCounter());
	   // get the boxes from all the leaves
	   boxes.push_back((*it)->getBox());
	   // get the volumes of all the leaves
	   volumes.push_back((*it)->nodeVolume());
	   // get fhat for all leaves
	   fhat.push_back(((*it)->getCounter())/((*it)->nodeVolume())/n);
	} // end of iterating through leaves 
	
	//----------------get the IAE-----------------------------------------------
	dotprecision dpIAE, dpIAEBoun;
	dpIAE = 0.0;
	int Nbin=counts.size();
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	double result = 0.0;
	double error;
	gsl_function F;
	
	F.function = &FinMixAbs;
	F.params =  &mixt;
	
	for (int j=0; j< Nbin; j++){
	  mixt.fhat = fhat[j];
	  rvector xuppVec = Sup(boxes[j]);
	  double xupp = _double(xuppVec[1]);
	  rvector xlowVec = Inf(boxes[j]);
	  double xlow = _double(xlowVec[1]);
	  gsl_integration_qags(&F, xlow, xupp, 0, 1e-7, 1000, w, &result, &error);
	  accumulate(dpIAE, result, 1.0);
	}
	
/*	// Accounting for the boundaries
	rvector xuppVec1 = Sup(boxes[Nbin-1]);
	double xupp1 = _double(xuppVec1[1]);
	rvector xlowVec1 = Inf(boxes[0]);
	double xlow1 = _double(xlowVec1[1]);
	dpIAEBoun = dpFinMixIAEBoun(xlow1, xupp1, Weight, Mean, Sigma);
	dpIAE += dpIAEBoun;
*/	
	// cast dot precision to real
	real FinMixIAE = rnd(dpIAE);
	
	// free the workspace
	gsl_integration_workspace_free (w);
	
	return FinMixIAE;
}
	
/*! Get the IAE for a finite mixture distribution
*/
real AdaptiveHistogramValidation::get2DIAE(taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval))
{
  //number of points
  int n = rootVpaving->getCounter();
	 
  SPSVnodePtrs leaves; // set up empty container for leaf node pointers
  SPSVnodePtrsItr it; // and an iterator over the container
  (*this).getSubPaving()->getLeaves(leaves); // fill the container

  // set up for taylor integration
 // taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval);
  real tol=1e-6;
  int o=16;
 // testpnt=BiGOP;

  real result = 0; 
  for (it=leaves.begin(); it < leaves.end(); it++)
  {	
	  //get fhat	
     interval fhat = interval(real((*it)->getCounter()/
	                           (((*it)->nodeVolume())*1.0*n))); 
	   
	  //get domain
	  ivector domain = (*it)->getBox();
     // get the integrated absolute error at this box
	  interval resultInt = integrateWithSplitting(testpnt, fhat, domain, o, tol);    
     //add the errors    
     result += Sup(resultInt);
   }
	//accounting for boundaries - will have to think about this later perhaps
	//accumulate(dpIAE,gsl_cdf_ugaussian_P(xlow[0]),1.0);
	//accumulate(dpIAE,gsl_cdf_ugaussian_Q(xupp[nLeaves-1]),1.0);
	return result;	
}

// For samples drawn from MappedSPnodes
// method for data splitting and hold out estimation
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// or until a stopping criteria is fulfilled
// outputs to a log file if logging is true
// makes its own random number generator
/*
AdaptiveHistogramVCollator AdaptiveHistogramValidation::prioritySplitAndEstimate(
									const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, int maxCheck, size_t hist,
									size_t maxLeafNodes, RealMappedSPnode& nodeEst)
{
    gsl_rng * rgsl = NULL;

    AdaptiveHistogramVCollator coll;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // call the function with a random number generator
        coll = prioritySplitAndEstimate(compTest, he, logging, minChildPoints, 
											  minVolB, rgsl, stopCrit,maxCheck, hist,
											  maxLeafNodes, nodeEst);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority stage split.  Orginal error: "
                                     + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority stage split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
   
   return coll;
}
*/

// Get the IAE for a finite gaussian mixture distribution using interval 
// techniques.
cxsc::interval AdaptiveHistogramValidation::getFinMixIntervalIAE(FinMix& mixt, double tol, int deg, bool full)
{
	//cout << "get finmix interval IAE" << endl;
	
	interval totalArea(0.0); //initialize
	int n = rootVpaving->getCounter();
	int N = n + rootVpaving->getVcounter();
	
	// need to iterate through the leaves
	SPSVnodePtrs leaves; // set up empty container for leaf node pointers
	SPSVnodePtrsItr it; // and an iterator over the container
	getSubPaving()->getLeaves(leaves); // fill the container
	
	// container is filled by reading leaves off tree from left to right
	for(it = leaves.begin(); it < leaves.end(); it++) {
		//cout << "-----IAE for " << (*it)->getNodeName() << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;
		
		//get the height in this leaf node
		double fhat;
		if ( full == 0 ) { 
			fhat = (*it)->getCounter()/(*it)->nodeVolume()/n;
			//cout << (*it)->getCounter() << "\t" << (*it)->nodeVolume() <<"\t" << fhat << endl;
		}
		else {
			size_t totalCount = (*it)->getCounter() + (*it)->getVcounter();
			fhat = totalCount/(*it)->nodeVolume()/N;
		}

		//get the box of this leaf node
		ivector thisBox = (*it)->getBox();
		//cout << (*it)->getBox() << endl;
		
		//---------find the root at this domain
		// make an intervalw object using thisBox
		rvector lb = Inf(thisBox);
		rvector ub = Sup(thisBox);
		intervalw thisIntW(_double(lb[1]), _double(ub[1]));
		interval thisInt(_double(lb[1]), _double(ub[1]));
		
		// find the root
		//cout << "finding roots at this node " << thisInt << endl;
		bisect(thisIntW, tol, fhat, rootVec, mixt.W, mixt.M, mixt.S); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				//cout << "root " << i << ": " << rootVec[i] << endl;
				//first insert into uniqueRootVec
				uniqueRootVec.push_back(rootVec[i]);
				//cout << i-1 << "\t" << i << "\t" << i+1 << endl;
				//now check for uniqueness
				if (((i-1) >= 0) && (i < rootVec.size())) {
					//cout << rootVec[i] << "\t" << rootVec[i-1] << endl;
					bool uniq = (subset(abs(rootVec[i] - rootVec[i-1]), intervalw(0, 1e-10)));
					if ( uniq ) { 
						//cout << "this root has a duplicate" << endl;
						uniqueRootVec.pop_back(); }
				}
			}
			//cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;
			
			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				//cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
					//cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = getL1error(fhat, thisSubIntLeft, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = getL1error(fhat, thisSubIntRight, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				//cout << "let's have a look at all the roots:" << endl;
				//for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					//cout << uniqueRootVec[i] << endl;
				//}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast, deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty
	} // end of iterating through the leaf nodes
	
	//cout << "IAE: " << totalArea << endl;
	return totalArea;
}

/*
// get the IAE for mapped functions
real AdaptiveHistogramValidation::getMappedFunctionIAE(RealMappedSPnode& nodeEst, bool full)
{
	ivector thisBox = getSubPaving()->getBox();
	RealMappedSPnode histMap(thisBox);
	int n = rootVpaving->getCounter();
	int allN = n + rootVpaving->getVcounter();
	
	// split the root box into the shape of myHist
	string leafLevelString = getLeafLevelsString();
	int depth = atoi(leafLevelString.c_str());
	if (depth != 0) {
		histMap.splitToShape(leafLevelString); 
	}

	//container to store heights for histNodes 
	vector< RangeCollectionClass<real> > heightHist;
	//get all the nodes in the histogram 
	SPSVnodePtrs histNodes;
	SPSVnodePtrsItr histNodeIt;
	getSubPaving()->getAllNodes(histNodes); 

	//traverse the tree and get the heights 
	//cout << "get the height at each node" << endl;
	for (histNodeIt = histNodes.begin(); histNodeIt < histNodes.end(); 
			histNodeIt++) {
				
		//get the height in this leaf node
		real fhat;
		if ( full == 0 ) { 
			fhat = (*histNodeIt)->getCounter()/((*histNodeIt)->nodeVolume()*n);
			//cout << (*histNodeIt)->getCounter() << "\t" << ((*histNodeIt)->nodeVolume())<< "\t" << n << endl;
		}
		else {
			size_t totalCount = (*histNodeIt)->getCounter() + (*histNodeIt)->getVcounter();
			fhat = totalCount/((*histNodeIt)->nodeVolume()*allN);
			//cout << totalCount << "\t" << fhat << endl;
		}
		
		//cout << (*histNodeIt)->getNodeName() << "\tfhat: " << fhat << endl;
		
		//get the height at each node
		RangeCollectionClass<real> height(fhat);
		heightHist.push_back(height);
	} // end of traversing all nodes in histogram
	
	//allocate ranges for histNode
	histMap.allocateRanges(heightHist, 0);
	return nodeEst.getMappedSPIAE(histMap);
}
*/


// ----------------------------- non member functions
//Output all boxes in AdaptiveHistogramValidation adh
std::ostream & operator<<(std::ostream &os, const AdaptiveHistogramValidation& adh)
{
    os << (adh.getSubPaving())->nodesAllOutput(os, 1) << std::endl;

    return os;
}

// check whether we can stop splitting using some stopping criteria
size_t checkNumValley(vector<double> vecMaxDelta, vector<int>& valleyHistPos, 
							bool& plateau, int& smallestDeltaInd)
{
	size_t flagValley = 0;
	size_t flagSame = 0;
	size_t flagSameLarger = 0;
	int Prev = 1;
	double currentSmallest = vecMaxDelta[0];
	
	for (size_t i = 1; i < (vecMaxDelta.size()); i++) {
		//cout << "====current: " << vecMaxDelta[i] << "\t previous: " << vecMaxDelta[i-1] << "=====" << endl;
		//cout << "Prev before checks: " << Prev << endl;
		double stopCritCurrent = vecMaxDelta[i];
		double stopCritPrevious = vecMaxDelta[i-1];
	
		// check if it is a local minimum
		if ( ((stopCritCurrent > stopCritPrevious) && (Prev == 1)) ) { //|| ((flagSame >= 20) && (Prev == 1)) ) {
			//cout << "larger: " << (stopCritCurrent > stopCritPrevious) << "\t Prev: " << Prev << endl;
			
			if ( stopCritPrevious <= currentSmallest) {
				currentSmallest = stopCritPrevious;
				//valleyHistPos.pop_back(); //only keep the smallest delta
				valleyHistPos.push_back(i-1);
				smallestDeltaInd = i-1;
				Prev = 0;
				flagSame = 0;
				flagValley++;
				//cout << "!!! Local minima previously at split !!! " << i-1 << endl;
			}
		}
		// check if stopCrit is decreasing
		else if (stopCritCurrent < stopCritPrevious) {
			flagSameLarger = 0;
			//cout << "current < previous " << endl;
			Prev = 1;
		}
		// if stopCritCurrent = stopCritPrevious
		else if ( (stopCritCurrent == stopCritPrevious) && (Prev == 1) ) {
			flagSame++;
			//cout << "same smaller: " << flagSame << endl;
		}
		//  increasing
		else if ( (stopCritCurrent == stopCritPrevious) && (Prev == 0) ) {
			flagSameLarger++;
			//cout << "same larger: " << flagSameLarger << endl;
		}
		//cout << "Prev after checks: " << Prev << endl;
		
		if (flagSameLarger >= 20) { 
			plateau = 0;
			//cout << "Perhaps have arrived at a plateau." << endl; 
		}
	}
	
	//cout << "There are " << flagValley << " valleys." << endl;
	return flagValley;
}

void getCurrentYatClass(vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >& vecRowYatSet,
vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >& vecColYatSet, 
list < set<CollatorSPVnode*, less < CollatorSPVnode* > > >& listYatSet)
{
	vector< set<CollatorSPVnode*, less < CollatorSPVnode* > > >::iterator vecIt;
	
	//insert vecRowYatSet and vecColYatSet into listYatSet
	
	if (!vecRowYatSet.empty()) {
		for (vecIt =  vecRowYatSet.begin(); vecIt < vecRowYatSet.end(); vecIt++){
				listYatSet.push_back(*vecIt);
		}
	}
	
	if (!vecColYatSet.empty()) {
		for (vecIt =  vecColYatSet.begin(); vecIt < vecColYatSet.end(); vecIt++){
			listYatSet.push_back(*vecIt);
		}
	}	
	
	if (!listYatSet.empty()) {
		//sort and get unique
		listYatSet.sort();
		listYatSet.unique();
	}
}

// Get the true delta for a finite gaussian mixture distribution using interval 
// techniques.
cxsc::interval getFinMixIntervalTrueDelta(
		FinMix& mixt, double tol, int deg, 
		std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet)
{
	interval totalArea(0.0); //initialize
	interval muValids(0.0);

	// need to iterate through the nodes
	set<CollatorSPVnode*, less < CollatorSPVnode* > >::iterator it;

	// container is filled by reading leaves off tree from left to right
	for(it = YatSet.begin(); it != YatSet.end(); it++) {
		//cout << "-----------------" << endl;
		//cout << (*it)->getNodeName() << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;

		//get the box of this leaf node
		ivector thisBox = (*it)->getBox();
		//cout << (*it)->getBox() << endl;
		
		//---------find the root at this domain
		// make an intervalw object using thisBox
		rvector lb = Inf(thisBox);
		rvector ub = Sup(thisBox);
		intervalw thisIntW(_double(lb[1]), _double(ub[1]));
		interval thisInt(_double(lb[1]), _double(ub[1]));
		
		// find the root
		//cout << "finding roots at this node " << thisInt << endl;
		double fhat = 0;
		bisect(thisIntW, tol, fhat, rootVec, mixt.W, mixt.M, mixt.S); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				//cout << "root " << i << ": " << rootVec[i] << endl;
				//first insert into uniqueRootVec
				uniqueRootVec.push_back(rootVec[i]);
				//cout << i-1 << "\t" << i << "\t" << i+1 << endl;
				//now check for uniqueness
				if (((i-1) >= 0) && (i < rootVec.size())) {
					//cout << rootVec[i] << "\t" << rootVec[i-1] << endl;
					bool uniq = (subset(abs(rootVec[i] - rootVec[i-1]), intervalw(0, 1e-10)));
					if ( uniq ) { 
						//cout << "this root has a duplicate" << endl;
						uniqueRootVec.pop_back(); }
				}
			}
			//cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;
			
			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				//cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
					//cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = getL1error(fhat, thisSubIntLeft, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = getL1error(fhat, thisSubIntRight, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				//cout << "let's have a look at all the roots:" << endl;
				//for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					//cout << uniqueRootVec[i] << endl;
				//}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast, deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty
		
		// get the empirical measure of the validation data at this node
		interval muValid((*it)->getVemp(), (*it)->getVemp());
		muValids = muValids + muValid;
	} // end of iterating through the leaf nodes

	// get the difference
	interval trueDelta = totalArea - muValids;
	
	return trueDelta;
}


// Get the true delta for a uniform mixture distribution 
cxsc::real getUnifTrueDelta(
		AdaptiveHistogram& myPart, double weight, vector<int> holesLoc, 
		std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet)
{
	
	//cout << "Calling dunction: " << endl;
	real totalArea = 0; //initialize
	real muValids = 0;

	// get the true height, f of the corresponding box in myPart
	SPSnodePtrs trueLeaves;
	SPSnodePtrsItr trueIt;
	(myPart).getSubPaving()->getLeaves(trueLeaves);
	double trueF; //true density
	ivector temp;

	dotprecision dpIAE;    // use type dotprecision for summation  
	dpIAE=0.0;
	
	// need to iterate through the nodes
	set<CollatorSPVnode*, less < CollatorSPVnode* > >::iterator it;

	// container is filled by reading leaves off tree from left to right
	for(it = YatSet.begin(); it != YatSet.end(); it++) {
		ivector thisBox = (*it)->getBox();
		//cout  << "=====checking " << thisBox << "======" << endl;
		size_t L = 0;
		for (trueIt = trueLeaves.begin(); trueIt < trueLeaves.end(); trueIt++) {
			//cout << "----True leaf: " << (*trueIt)->getBox() << "\t" << endl;
			ivector trueBox = (*trueIt)->getBox();

			if (  holesLoc[L] == 0 ) { trueF = 0; }
			else { trueF = weight/((*trueIt)->nodeVolume()); }
			//cout << "pdf: " << trueF << endl;
			
			// if this is contained in trueBox
			if ( (*it)->getBox() <= (*trueIt)->getBox() || (*it)->getBox() == (*trueIt)->getBox() ) {
				//use the volume of this
				real r = ((*it)->nodeVolume())*(trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
				//can move on to next leaf rather than iterating thru all trueBoxes
				//think about this later
			} //end of if this box is in trueBox
			
			// if this contains trueBox
			else if ((*trueIt)->getBox() <= (*it)->getBox()) {
				//use the volume of trueBox
				real r = ((*trueIt)->nodeVolume())*(trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
			} //end of if trueBox is in this box
			
			// if this is partially contained in trueBox 
			else if 	(Intersection(temp, thisBox, trueBox)) {
				if (Inf(temp) != Sup(temp)){
					double volume = Volume(temp);
					real r = volume*(trueF);
					//cout << "r: " << r << "\t" << abs(r) << endl;
					accumulate(dpIAE, abs(r), 1.0);
				}
			}
			L++;
		} // end of going through trueBoxes

		// get the empirical measure of the validation data at this node
		real muValid = (*it)->getVemp();
		muValids = muValids + muValid;
		//cout << "Area \t Emp Mass" << endl;
		//cout << rnd(dpIAE) << "\t" << muValids << endl;
	} // end of iterating through the leaf nodes

	// get the difference
	totalArea = rnd(dpIAE);
	real trueDelta = (totalArea - muValids);
	
	return trueDelta;
}

// get the IAE for mapped functions
real getMappedFunctionTrueDelta(PiecewiseConstantFunction& nodeEstHist,
					std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet)
{
	//cout << "-------function called---------" << endl;
	
	//iterator for YatSet
	std::set<CollatorSPVnode*, less < CollatorSPVnode* > >::iterator histNodeIt;

	//traverse the tree and get the heights 
	real trueArea = 0.0;
	real muValid = 0.0;

	for (histNodeIt = YatSet.begin(); histNodeIt != YatSet.end(); 
			histNodeIt++) {
		
		real thisArea = 0.0;
		
		cout << (*histNodeIt)->getNodeName() << endl;
		ivector thisBox = (*histNodeIt)->getBox();

		// need to get the area of the nodes of nodeEst in thisBox
		//cout << nodeEstHist.getDimensions() << "\t" << nodeEstHist.getRootBox() << endl;
		//thisArea = nodeEstHist.getArea(thisArea, thisBox);

		//cout << thisArea << "\t" << (*histNodeIt)->getVemp() << endl;
		 
		trueArea += thisArea;
		muValid += (*histNodeIt)->getVemp();

	} // end of traversing iterating through YatSet

	//cout << "Final: " << endl;
	//cout << trueArea << "\t" << muValid << endl;
	
	real trueDelta = trueArea - muValid;
	return abs(trueDelta);
}

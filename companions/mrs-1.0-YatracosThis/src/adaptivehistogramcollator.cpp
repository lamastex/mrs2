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

/*! \file adaptivehistogramcollator.cpp
\brief AdaptiveHistogramCollator definitions
*/

#include "adaptivehistogramcollator.hpp"

#include <string>   // to use the C++ string class
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <exception> // use exceptions
#include <functional> // use functionals

#include <gsl/gsl_rng.h>        // to use the gsl random number generator
#include <gsl/gsl_randist.h>

#include <math.h> // math library

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

//to use subpavings
#include "sptools.hpp"
#include "spalgorithms.hpp"

// to use stats subpavings
#include "spsnode.hpp" // includes spnode.hpp includes sptypes.hpp
#include "collatorspnode.hpp"

#include "adaptivehistogram.hpp"

// to use histogram penalty function objects
#include "histpenalty.hpp"

#include <limits> // to use negative infinity

using namespace subpavings;
using namespace std;


// ---------- implementation of AdaptiveHistogramCollator class -------------

// ----------------private methods

// initialised constructor, initialised with a subpaving pointer
AdaptiveHistogramCollator::AdaptiveHistogramCollator(CollatorSPnode * spn)
{
    if (NULL == spn) {
        throw HistException("Null CollatorSPnode pointer in constructor");
    }
    rootCollator = spn;
}

// Make a collated histogram from a container of rvectors
bool AdaptiveHistogramCollator::collateFromRVec(size_t samplesize,
                size_t numberSamples, const RVecData rv, ivector pavingBox,
                int indImmedSplit, const SplitDecisionObj& boolTest,
                const NodeCompObj& compTest, const HistEvalObj& he,
                size_t minChildPoints, double minVolB)
{

    bool retValue = false;
    gsl_rng * rgsl = NULL;

    try {

        size_t countIn = 0; // track the number of histograms made and added

        // set up a random number sampler
        const gsl_rng_type * tgsl;

        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();

        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        std::string fileName; // a name for the files to use


        // for loop to generate histograms and add to collation
        for (size_t j=1; j<=numberSamples; j++) {

            bool successfulInsertion; // recognise successes

            // make an Adaptive Histogram object with a specified box
            AdaptiveHistogram myHist(pavingBox);

            if (indImmedSplit == 1) { // doing immediate splitting

                successfulInsertion = myHist.insertSampleFromRVec(samplesize,
                    rgsl, rv, boolTest);
            }

            if (indImmedSplit == 0) { // doing priority queue splitting

                successfulInsertion = myHist.insertSampleFromRVec(samplesize,
                    rgsl, rv);

                bool successfulPQSplit;

                if (successfulInsertion) {

                    successfulPQSplit = myHist.prioritySplit(compTest, he,
                                                    NOLOG, minChildPoints, minVolB);
                }

                successfulInsertion = successfulInsertion && successfulPQSplit;

            }

            // only do more if some data was fed in
            if(successfulInsertion) {

                // create a name for the file to output
                fileName = "Hist";
                //convert j to a string
                std::ostringstream stm2;
                stm2 << j;
                // add the stringed j to the filename
                fileName += stm2.str();
                fileName += ".txt"; // and finish the filename


                // To realize a file output
                myHist.outputToTxtTabs(fileName);

                // add the histogram to the collection represented by this
                addToCollation(myHist);

                countIn++; // increment the counter
            }
        } // end of for loop creating histograms

        if (countIn == numberSamples) {

            retValue = true;
        }
        else { // did not add required number of histograms
            std::cerr << "Problem in collateFromRVec(): check "
                    << "console for error reports " << std::endl;
        }
        // free the random number generator
        gsl_rng_free (rgsl);

    }
    catch (exception&) {
        if (NULL != rgsl) // free the random number generator
            gsl_rng_free (rgsl);
        throw;
    }

    return retValue;
}


// Make a collated histogram from an RSSample
bool AdaptiveHistogramCollator::collateFromRSSample(size_t samplesize,
                size_t numberSamples, const RSSample rss, ivector pavingBox,
                int indImmedSplit, const SplitDecisionObj& boolTest,
                const NodeCompObj& compTest, const HistEvalObj& he,
                size_t minChildPoints, double minVolB, int label)
{
    // container to put the rvectors into
    RVecData allData;

    //get the container of rvectors
    //use getRvectorsFromRSSample to put rvectors from labeled points in
    // rss.Samples into allData where the labeled point label matches label
    size_t numberFound = getRvectorsFromRSSample(allData, rss, label);

    bool cancontinue = (numberFound > 0);
    // cancontinue will be false if there was a problem getting data points
    // if cancontinue is true data should contain at least some data points

    bool retValue = false;

    if (cancontinue) {

        // use the RVec method to complete the process of histogram
        // creation and averaging
        retValue = collateFromRVec(samplesize, numberSamples, allData,
                                pavingBox, indImmedSplit, boolTest,
                                compTest, he, minChildPoints, minVolB);
    }

    return retValue;
}

// Method to add current state of the histogram collator to a log file
// Output goes to file named according to argument s
void AdaptiveHistogramCollator::outputLog(const std::string& s, const int i) const
{
    // To add output of the AdaptiveHistogramCollator object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        getSubPaving()->leavesOutputTabs(os); // the output
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}


// --------------- end private methods

// ---------------- public methods

// default constructor
AdaptiveHistogramCollator::AdaptiveHistogramCollator()
{
    try {
        rootCollator = new CollatorSPnode();
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }

}


// initialised constructor, initialised with an AdaptiveHistogram object
AdaptiveHistogramCollator::AdaptiveHistogramCollator(const
    AdaptiveHistogram& adh)
{
    try {
        rootCollator = new CollatorSPnode(adh.getSubPaving());
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }
}


// copy constructor
AdaptiveHistogramCollator::AdaptiveHistogramCollator(
                const AdaptiveHistogramCollator& other)
{
    try {
        rootCollator = new CollatorSPnode(*(other.rootCollator));
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor: original error "
                                            << msg << std:: endl;
        throw HistException("Memory allocation error in constructor: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std::cerr << "SPnodeExcepton in constructor: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in constructor: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std::cerr << "Error in constructor: original error "
                                            << msg << std::endl;
        throw HistException("Error in constructor: " + msg);
    }
}

// assignment operator
AdaptiveHistogramCollator& AdaptiveHistogramCollator::operator=(const
    AdaptiveHistogramCollator& rhs)
{
    try {

        // we have to make sure we delete the current paving
        if (NULL != rootCollator) {
            delete rootCollator;
            rootCollator = NULL;
        }

        if (NULL != rhs.rootCollator)
            rootCollator = new CollatorSPnode(*(rhs.rootCollator));
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in assignment: original error "
                            << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std::cerr << "SPnodeExcepton in assignment: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in assignment: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std::cerr << "Error in assignment: original error "
                                            << msg << std::endl;
        throw HistException("Error in assignment: " + msg);
    }

}


// Destructor.
AdaptiveHistogramCollator::~AdaptiveHistogramCollator()
{ delete rootCollator; }

// addition operator
AdaptiveHistogramCollator AdaptiveHistogramCollator::operator+(const
    AdaptiveHistogramCollator& rhs) const
{
    if ((NULL != rootCollator) && (NULL != rhs.rootCollator) &&
    ((Ub(rootCollator->getBox()) != Ub(rhs.rootCollator->getBox()))
    || (Lb(rootCollator->getBox()) != Lb(rhs.rootCollator->getBox()))))
        throw HistException("Added histogram collators have unequal dimensions");

    CollatorSPnode* newnode = NULL;

    try {

        newnode =
            CollatorSPnode::addPavings(rootCollator, rhs.rootCollator);
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in addition: original error "
                                    << msg << std::endl;
        throw HistException("Memory allocation error in addition: " + msg);
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std::cerr << "SPnodeExcepton in addition: original error "
                                            << msg << std::endl;
        throw HistException("SPnodeException in addition: " + msg);
    }
    catch (exception& e) {
        string msg(e.what());
        std::cerr << "Error in addition: original error "
                                            << msg << std::endl;
        throw HistException("Error in addition: " + msg);
    }


    AdaptiveHistogramCollator newCollatorHist(newnode);

    return newCollatorHist;
}

// increment addition operator
AdaptiveHistogramCollator& AdaptiveHistogramCollator::operator+=(const
    AdaptiveHistogramCollator& rhs)
{
    try {
        // get the subpaving out of rhs to form a new CollatorSPnode
        CollatorSPnode toAdd(*rhs.getSubPaving());
        // add the new CollatorSPnode into the collation
        // note that addPaving will alter toAdd, but that is okay because
        // toAdd is a temporary object created and deleted in this procedure
        rootCollator->addPaving(&toAdd);
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
AdaptiveHistogramCollator AdaptiveHistogramCollator::operator-(const
    AdaptiveHistogramCollator& rhs) const
{
    if ((NULL != rootCollator) && (NULL != rhs.rootCollator) &&
    ((Ub(rootCollator->getBox()) != Ub(rhs.rootCollator->getBox()))
    || (Lb(rootCollator->getBox()) != Lb(rhs.rootCollator->getBox()))))
        throw HistException("Histogram collators have unequal dimensions");

    CollatorSPnode* newnode = NULL;

    try {

        double c = -1.0;
        newnode =
            CollatorSPnode::subtractPavings(rootCollator, rhs.rootCollator, c);
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


    AdaptiveHistogramCollator newCollatorHist(newnode);

    return newCollatorHist;
}


// Return a pointer to the CollatorPSnode this manages
CollatorSPnode* AdaptiveHistogramCollator::getSubPaving() const
{return rootCollator;} // boost::shared_ptr might be better


// averaging method
AdaptiveHistogramCollator AdaptiveHistogramCollator::makeAverage() const
{
    if (NULL == rootCollator) {
            string msg = "Cannot average this: rootCollator is NULL";
            throw HistException(msg);
    }

    //average only makes sense if all values in the summary are positive
    VecDbl mySummary = rootCollator->getSummary();

    VecDblIt it = find_if(mySummary.begin(), mySummary.end(),
                    bind2nd(less<double>(), 0.0));
    if (it < mySummary.end()) {
            string msg = "Cannot average this: the collation contains negatives";
            throw HistException(msg);
    }

    AdaptiveHistogramCollator newCollator;
    try {
        newCollator.rootCollator = (getSubPaving())->makeAverageCollation();
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




// Add an AdaptiveHistogram into the collation
void AdaptiveHistogramCollator::addToCollation(const AdaptiveHistogram& adh)
{
    try {
        // get the subpaving out of adh to form a new CollatorSPnode
        CollatorSPnode toAdd(adh.getSubPaving());
        // add the new CollatorSPnode into the collation
        // note that addPaving will alter toAdd, but that is okay because
        // toAdd is a temporary object created and deleted in this procedure

        //rootCollator->addPaving(&toAdd);

        bool successfullyAdded = rootCollator->addPaving(&toAdd);
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



// Add the negation of an AdaptiveHistogram into the collation
void AdaptiveHistogramCollator::addNegationToCollation(
                            const AdaptiveHistogram& adh, double c)
{
    try {
        // make the AdaptiveHistogram into a new CollatorSPnode
        CollatorSPnode toNeg(adh.getSubPaving());
       // put the negation of the new CollatorSPnode into the collation
       // note that addNegatedPaving will alter toNeg, but that is okay because
       // toNeg is a temporary object created and deleted in this procedure
        rootCollator->addNegatedPaving(&toNeg, c);
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory adding negation to collation. ";
        msg += "Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error adding negation to collation. ";
        msg += "Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException adding negation to collation. ";
        msg += "Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error adding negation to collation. ";
        msg += "Orginal error: " + oldmsg;
        throw HistException(msg);
    }
}


// Make a collated histogram from a container of rvectors
// immediate splitting
bool AdaptiveHistogramCollator::collateFromRVecSplitNow(size_t samplesize,
                    size_t numberSamples, const RVecData rv, ivector pavingBox,
                    const SplitDecisionObj& boolTest)
{
    bool retValue = false;

    try {
        int indImmedSplit = 1; // immediate splitting

        CompNothing compTest; // dummy comparison test
        CritStopAll he; // dummy stopping function object
        size_t minPoints = 0; // dummy minPoints
        double minVolB = 0.0; // dummy minVolB


        retValue = collateFromRVec(samplesize, numberSamples, rv,
                                    pavingBox, indImmedSplit, boolTest,
                                    compTest, he, minPoints, minVolB);
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory averaging.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return retValue;
}

// Make a collated histogram from an RSSample
// priority queue splitting
bool AdaptiveHistogramCollator::collateFromRVecSplitPQ(size_t samplesize,
                    size_t numberSamples, const RVecData rv, ivector pavingBox,
                    const NodeCompObj& compTest, const HistEvalObj& he,
                    size_t minChildPoints, double minVolB)
{
    bool retValue = false;

    try {

        int indImmedSplit = 0; // pq rather than immediate splitting

        SplitNever sn; // dummy

        retValue = collateFromRVec(samplesize, numberSamples, rv,
                                    pavingBox, indImmedSplit, sn,
                                    compTest, he, minChildPoints, minVolB);
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory averaging.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return retValue;
}

// Make an collated histogram from an RSSample
// immediate splitting
bool AdaptiveHistogramCollator::collateFromRSSampleSplitNow(size_t samplesize,
                    size_t numberSamples, const RSSample rss, ivector pavingBox,
                    const SplitDecisionObj& boolTest, int label)
{
    bool retValue = false;

    try {

        int indImmedSplit = 1; // immediate splitting

        CompNothing compTest; // dummy comparison test
        CritStopAll he; // dummy stopping function object
        size_t minPoints = 0; // dummy minPoints
        double minVolB = 0.0; // dummy minVolB

        retValue = collateFromRSSample(samplesize, numberSamples, rss,
                                    pavingBox, indImmedSplit, boolTest,
                                    compTest, he, minPoints, minVolB, label);
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory averaging.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return retValue;
}

// Make a collated histogram from an RSSample
// priority queue splitting
bool AdaptiveHistogramCollator::collateFromRSSampleSplitPQ(size_t samplesize,
                    size_t numberSamples, const RSSample rss,
                    ivector pavingBox, const NodeCompObj& compTest,
                    const HistEvalObj& he,
                    size_t minChildPoints, double minVolB, int label)
{
    bool retValue = false;

    try {
        int indImmedSplit = 0; // pq rather than immediate splitting

        SplitNever sn; // dummy

        retValue = collateFromRSSample(samplesize, numberSamples, rss,
                                    pavingBox, indImmedSplit, sn,
                                    compTest, he, minChildPoints, minVolB, label);
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory averaging.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error averaging.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return retValue;
}

// make a .dot file for the histogram
bool AdaptiveHistogramCollator::outputGraphDot() const
{


    if (NULL == rootCollator) {

        throw HistException("No root paving for graph output");
    }

    bool success = rootCollator->outputGraphDot();

    return success;
}

// Get the number of Adaptive Histogram objects collated.
size_t AdaptiveHistogramCollator::getNumberCollated() const
{ return rootCollator->getNumberSummarised(); }

// Output the collated normalised histogram heights and bins data to a txt file
void AdaptiveHistogramCollator::outputToTxtTabs(const std::string& s,
                                        bool confirm) const
{
    // To generate a file output of the AdaptiveHistogramCollator object
    ofstream os(s.c_str());         // Filename, c-string version

    rootCollator->leavesOutputTabs(os);
    if (confirm)
        std::cout << "The output of the AdaptiveHistogramCollator has been "
            << "written to " << s << std::endl << std::endl;

}


// get number of leaf nodes
size_t AdaptiveHistogramCollator::getNumLeaves()
{
	vector<CollatorSPnode*> leaves;
	getSubPaving()->getLeaves(leaves); 

	return leaves.size();
}

// Output the average data over the collation to a txt file
// this outputs the normalised average histogram heights and bins
void AdaptiveHistogramCollator::outputAverageToTxtTabs(const
    std::string& s, bool confirm) const
{
    try {

        if (NULL != rootCollator) {

            //average only makes sense if all values in the summary are positive
            VecDbl mySummary = rootCollator->getSummary();

            VecDblIt it = find_if(mySummary.begin(), mySummary.end(),
                                                bind2nd(less<double>(), 0.0));
            if (it < mySummary.end()) {
                    string msg = "Cannot average this: the collation contains negatives";
                    std::cerr << "\n" << msg << "\n" << std::endl;
            }
            else {

                // To generate a file output of the AdaptiveHistogramCollator object
                ofstream os(s.c_str());         // Filename, c-string version

                if (NULL != rootCollator) {
                    rootCollator->leavesAverageOutputTabs(os);
                    if (confirm)
                        std::cout << "The output of the average AdaptiveHistogram has been "
                            << "written to " << s << std::endl << std::endl;
                }
                else {
                    std::cerr << "Sorry, nothing is in collation to average"
                        << std::endl;
                }
            }
        }
    }

    catch (exception& e) {
        std::cerr << "Problem averaging: " << e.what() << std::endl;
    }

}

void AdaptiveHistogramCollator::publicOutputLog(const std::string& s, const int i) const
{
    // use the private version
    outputLog(s, i);
}

//--this was removed then re-inserted-----
// Output the accumulated data over the collation to a txt file
// this outputs the sum over the collation summary
void AdaptiveHistogramCollator::outputAccumulationToTxtTabs(const
    std::string& s) const
{
    try {
         if (NULL != rootCollator) {

            // To generate a file output of the AdaptiveHistogramVCollator object
            ofstream os(s.c_str());         // Filename, c-string version

            if (NULL != rootCollator) {
                rootCollator->leavesAccumulationOutputTabs(os);
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

/*! Distribution-free Likelihood Estimation
*/
real AdaptiveHistogramCollator::getEstLogLikelihoodFromRSSample(
								RSSample& labSampledData, double dx, double wt,
								double WeightHist,
						std::map<rvector, double, std::less<rvector> >& WeightsPM)
{
	if (NULL != rootCollator) {			
		// container to store heights
		vector<double> fhatNew;         
      // make new fhat by adding some mass
      // maybe can update the summary instead of outputting as vector of
      // doubles
		rootCollator->leavesMakeNewFhat(wt, fhatNew);

		// get pointers to the leaf nodes
		vector<CollatorSPnode*> leaves;
		vector<CollatorSPnode*>::iterator it;
		getSubPaving()->getLeaves(leaves); 
		vector<LabPnt>::iterator dataIt;

		dotprecision dpEstLogLik;
		dpEstLogLik = 0.0;
		//iterated through labSampledData
		for (dataIt = labSampledData.Samples.begin(); 
			dataIt < labSampledData.Samples.end(); dataIt++) {
			bool done = false;
			while ( !done ) {
				// start the counts for the leaves so that it corresponds to the 
				// right fhat
				size_t pos = 0;
				//iterate through the leaves
				for (it = leaves.begin(); it < leaves.end(); it++) {
					//get the box of this leaf node
					ivector theBox = (*it)->getBox();
					//make a pointer to an SPSnode with the box
					SPSnode* newNode = NULL;
					newNode = new SPSnode(theBox, 0);
				
					//now check if this data is inside the box
					if ( newNode->nodeContains((*dataIt).Pnt, ON_PARENT) ) {  
						
						//if data is a point mass in current simulated dataset from theta
						if ( (*dataIt).L == 0 && WeightsPM[(*dataIt).Pnt] != 0 ) {			 
							accumulate(dpEstLogLik, 
							log(WeightsPM[(*dataIt).Pnt] + dx*WeightHist*fhatNew[pos]), 1);				
						}
						
						//if data is a point mass but not in current simulated dataset from theta
						else if ( (*dataIt).L == 0 && WeightsPM[(*dataIt).Pnt] == 0 ) {			 
							accumulate(dpEstLogLik, log(dx*WeightHist*fhatNew[pos]), 1);				
						}	

						//if data is unique
						else {
							accumulate(dpEstLogLik, log(dx*WeightHist*fhatNew[pos]), 1);
						}
						
						done = true; //so that don't have to iterate through ALL the leaves
										  //just for this data
						break;
					}
					pos++; // increment position
					delete newNode; //free memory
				} // end of iterating though leaves
			} // end of while loop 
			
		} // end of going through the data
			
		real estLogLik = rnd(dpEstLogLik);
		return estLogLik;
	} // end of if NULL != rootCollator
	else { cerr << "Empty collator" << endl; exit(1); }
}

//marginalisation
AdaptiveHistogramCollator AdaptiveHistogramCollator::marginalise(
								const std::vector<int>& reqDims) const
{
	// take the root of the other and marginalise it
	return AdaptiveHistogramCollator(CollatorSPnode::marginalise(getSubPaving(),
	reqDims));
}

//find density region with cov coverage
void AdaptiveHistogramCollator::findDensityRegion(double cov, double weightPM,
														vector<CollatorSPnode*> & covNodes,
														string covFileName)
{
	try {
			
			if ( (cov - weightPM) <= 0) {
				cout << cov << " percent of the mass are already covered by the point masses" << endl;
			}
			else {
			
				// put the leaves into a vector and sort it, smallest to largest
				vector<CollatorSPnode*> leaves;
				getSubPaving()->getLeaves(leaves);
				CompHeight compheight;	
				//sort according to average height
				sort(leaves.begin(), leaves.end(), nodeCompTotalSummaryAv);
			
				//start iterating from the largest
				vector<CollatorSPnode*>::reverse_iterator rit = leaves.rbegin();
				bool found = FALSE; //found the boxes that gives cov density region
				
				dotprecision totalCov;
				totalCov = 0.0;
				
				while (!found && rit < leaves.rend()) {
					// double check this:
					// height is box counts/box vol
					// ie count is unnormalised vol of an individual element of histogram
					// box vol * height == box vol * (count / box vol) == count
				
					//accumulate the summary * box vol
					accumulate(totalCov, (*rit)->getTotalSummaryAv()*(*rit)->nodeVolume(), 1); 

					//push back the node that fulfill the condition totalCov <= cov 
					//into the container covNodes
					if (totalCov <= (cov - weightPM) ) { 
						covNodes.push_back((*rit)); 
					} 
					
					// check that totalCov is at most cov
					if (totalCov >= (cov - weightPM) ) { found = TRUE; } // break out of loop
					++rit;				
				} // end while 
			
			//output covNodes to .txt	
				ofstream os;
				os.open(covFileName.c_str());		
				vector<CollatorSPnode*>::iterator vit;
				for (vit = covNodes.begin(); vit < covNodes.end(); vit++) {
					ivector thisBox = (*vit)->getBox(); // copy theBox         
					double vol = (*vit)->nodeVolume();
					// output the nodeName, nodeVolume
					os << (*vit)->getNodeName();
					os << "\t" << vol;
					// followed by the average
					os << "\t" << (*vit)->getTotalSummaryAv();
					// followed by intervals making up box using Inf & Sup
					// ie unlike cxsc output, there is no [  ] around them
					for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {
						 os << "\t" << Inf(thisBox[i]) << "\t" << Sup(thisBox[i]);
					}
					os << endl;
				}
				os << flush;
				os.close();
			} 
		} // end of try			
	catch (exception& e) {
		throw HistException(
		"Error in AdaptiveHistogramCollation::coverage :\n"
		+ string( e.what() ) );
	}
}

// Get a string of the leaf node levels.
std::string AdaptiveHistogramCollator::getLeafLevelsString() const
{
    string retValue = "";
    if (NULL != rootCollator)
        retValue = rootCollator->getLeafNodeLevelsString();

    return retValue;
}

// get the sum of the variances over samples contained in a collation
// for a scalar summary of a histogram where the variance of the scalar
// for a histogram sample is the square of the sum of the areas of difference
// between a histogram sample and the average histogram over the samples.
real AdaptiveHistogramCollator::getSumVarianceAreaScalar() const
{
    if (NULL == getSubPaving()) {
            string msg = "Cannot get variances for this: rootCollator is NULL";
            throw HistException(msg);
    }

    //variances only makes sense if all values in the summary are positive
    VecDbl mySummary = getSubPaving()->getSummary();
    VecDblIt it = find_if(mySummary.begin(), mySummary.end(),
                                        bind2nd(less<double>(), 0.0));
    if (it < mySummary.end()) {
            string msg = "Cannot get variances: collation contains negatives";
            throw HistException(msg);
    }

    real sumVars = 0.0;
    // take this collation
    try {
        sumVars = getSubPaving()->getSumVarsAreaScalar();
        //cout << "sumVars: " << sumVars << endl;
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error summing variances.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    return sumVars;
}


// get the sample variance over sample contained in a collation for a scalar
// summary value defined so that the variance of a sample histogram in the
// collation is the sum of the squares of the areas of difference between the
// histogram and the average histogram over the sample
// the sample variance is (sum of the sample variances)/(number samples - 1)
real AdaptiveHistogramCollator::getSampleVarianceAreaScalar() const
{
    real sumVars = getSumVarianceAreaScalar();
    size_t numberColl = getNumberCollated();

    if (numberColl < 2) {
        throw HistException("Cannot do sample variance for sample size < 2");
    }

    return sumVars/(1.0*(numberColl - 1));
}


// get the sum of the variances over samples contained in a collation
// for a scalar summary of a histogram which is the total heights of all
// the bins in the histogram.
real AdaptiveHistogramCollator::getSumVarianceTotalHeightScalar()
                                                                        const
{
    if (NULL == getSubPaving()) {
            string msg = "Cannot get variances for this: rootCollator is NULL";
            throw HistException(msg);
    }

    //variances only make sense if all values in the summary are positive
    VecDbl mySummary = getSubPaving()->getSummary();

    VecDblIt it = find_if(mySummary.begin(), mySummary.end(),
                                            bind2nd(less<double>(), 0.0));
    if (it < mySummary.end()) {
            string msg = "Cannot get variances: collation contains negatives";
            throw HistException(msg);
    }

    real sumVars = 0.0;

    // take this collation
    try {
        sumVars = getSubPaving()->getSumVarsTotalSummarisedValueScalar();
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException summing variances.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error summing variances.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }

    return sumVars;

}

// get the sample variance over sample contained in a collation
// for a scalar summary of a histogram which is the total heights of all
// the bins in the histogram.
// the sample variance is (sum of the sample variances)/(number samples - 1)
real AdaptiveHistogramCollator::getSampleVarianceTotalHeightScalar()
                                                                        const
{
    real sumVars = getSumVarianceTotalHeightScalar();
    size_t numberColl = getNumberCollated();

    if (numberColl < 2) {
        throw HistException("Cannot do sample variance for sample size < 2");
    }

    return sumVars/(1.0*(numberColl - 1));
}

// Jenny addition for Gloria's convergence work
// make a collator that is the differences of each element in a sample to the sample average
AdaptiveHistogramCollator AdaptiveHistogramCollator::makeDifferencesToAverage() const
{
	try {
		if (NULL == getSubPaving()) {
				string msg = "Cannot make differences to average for this: rootCollator is NULL";
				throw HistException(msg);
		}

		AdaptiveHistogramCollator newCollator;
    
        newCollator.rootCollator = (getSubPaving())->makeDifferencesToAveragePaving();
    
		return newCollator;

    }
    
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = string("Error in AdaptiveHistogramCollator::getL1DistancesToAverage.")
				+  string("Orginal error: ") + oldmsg;
        throw HistException(msg);
    }
} //end of function makeDifferencesToAverage()

// Jenny addition for Gloria's convergence work
// take a container and return the same container, which has been
// cleared (if necessary) and re-filled with 
// L1-distances-to-average, one for each histogram in collation
RealVec& AdaptiveHistogramCollator::getL1DistancesToAverage(RealVec& container) const
{
	try {
		if (NULL == getSubPaving()) {
				string msg = "Cannot get L1 distances for this: rootCollator is NULL";
				throw HistException(msg);
		}

		container = getSubPaving()->getL1DistancesToAverage(container);
	   
		return container;

    }
    
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = string("Error in AdaptiveHistogramCollator::getL1DistancesToAverage.")
				+  string("Orginal error: ") + oldmsg;
        throw HistException(msg);
    }

}

//gloria's addition to add to trunk later
/*! Get the IAE for a uniform mixture.
*/
real AdaptiveHistogramCollator::getUnifIAE(AdaptiveHistogram & myPart, 
									vector<int> holesLoc, double weight)
{                                               
   // get the true height, f of the corresponding box in myPart
	SPSnodePtrs trueLeaves;
	SPSnodePtrsItr trueIt;
	//AdaptiveHistogram * adhPtr;
	//adhPtr = &myPart;
	(myPart).getSubPaving()->getLeaves(trueLeaves);

	// setting up containers for the leaves
	vector<CollatorSPnode*> leaves; // set up empty container for leaf node pointers
	vector<CollatorSPnode*>::iterator it; // and an iterator over the container
	(*this).getSubPaving()->getLeaves(leaves); // fill the container

	double trueF; //true density
	ivector temp;
	
	 dotprecision dpIAE;    // use type dotprecision for summation  
   dpIAE=0.0;

	
	//go through all the leaves in this
	for(it = leaves.begin(); it < leaves.end(); it++) {
		ivector thisBox = (*it)->getBox();
		//cout << "====checking " << (*it)->getBox() << endl;
      
		// get the height of this leaf
		double fhat = (*it)->getSummary()[0]; //hopefully this is a robust way to
														  //get the needed value
		
		//cout << "fhat for box " << ":" << fhat << endl;

		size_t L = 0;
		for (trueIt = trueLeaves.begin(); trueIt < trueLeaves.end(); trueIt++) {
			//cout << "----True leaf: " << (*trueIt)->getBox() << "\t" << endl;
			ivector trueBox = (*trueIt)->getBox();

			if (  holesLoc[L] == 0 ) { trueF = 0; }
			else { trueF = weight/((*trueIt)->nodeVolume()); }
			//cout << "pdf: " << trueF << "------" << endl;
			
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
	return unifIAE;														  

} //end of function getUnifIAE()

/*! Get the IAE for a standard uniform distribution.
*/
real AdaptiveHistogramCollator::getUnifIAE()
{
	// setting up containers for the leaves
	vector<CollatorSPnode*> leaves; // set up empty container for leaf node pointers
	vector<CollatorSPnode*>::iterator it; // and an iterator over the container
	(*this).getSubPaving()->getLeaves(leaves); // fill the container

   dotprecision dpIAE;    // use type dotprecision for summation  
   dpIAE=0.0;
   //go through all the leaves in this
   for(it = leaves.begin(); it < leaves.end(); it++) {
      // get the height of this leaf
      // the height is contained in the first position of the summary vector
      // so we need to access the vector
      double fhat = (*it)->getSummary()[0]; //hopefully this is a robust way to
														  //get the needed value
		// cout << "fhat for box " << (*it)->getNodeName() << ":" << fhat << endl;

		//now calculate the IAE
      if ((1 - fhat) < 0.0){
			real r = ((*it)->nodeVolume())*(fhat - 1);
			//cout << "r: " << r << endl;
			accumulate(dpIAE, r, 1.0);
		}

		else if ((1 - fhat) > 0.0){
			real r = ((*it)->nodeVolume())*(1 - fhat);
			//cout << "r: " << r << endl;
			accumulate(dpIAE, r, 1.0);
		}
	} // end of going through all the leaves in this

   //cast dotprecision to real
   real unifIAE = rnd(dpIAE);
	return unifIAE;
} //end of function getUnifIAE()

// Get the IAE for a finite gaussian mixture distribution using interval 
// techniques.
cxsc::interval AdaptiveHistogramCollator::getFinMixIntervalIAE(FinMix& mixt, double tol, int deg)
{
	interval totalArea(0.0); //initialize
	
	// need to iterate through the leaves
	vector<CollatorSPnode*> leaves; // set up empty container for leaf node pointers
	vector<CollatorSPnode*>::iterator it; // and an iterator over the container
	getSubPaving()->getLeaves(leaves); // fill the container
	
	// container is filled by reading leaves off tree from left to right
	for(it = leaves.begin(); it < leaves.end(); it++) {
		//cout << "-----------------" << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;
		
		//get the height in this leaf node
		double fhat = ((*it)->getSummary())[0];
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








// gloria addition 
/*! \brief Returns the IAE between an AdaptiveHistogram object and a RealMappedSPnode.
*/
/*cxsc::real AdaptiveHistogramCollator::getMappedIAE(RealMappedSPnode& nodeEst, 
												ivector pavingBox) const
{
	RealMappedSPnode histMap(pavingBox);
	
	// split the root box into the shape of this 
	string leafLevelString = getLeafLevelsString();
	int depth = atoi(leafLevelString.c_str());
			if (depth != 0) {
				histMap.splitToShape(leafLevelString); 
	}

	//container to store heights for histNodes 
	vector< RangeCollectionClass<real> > heightHist;

	//get all the nodes in the histogram 
	vector<CollatorSPnode*> histNodes;
	vector<CollatorSPnode*>::iterator histNodeIt;
	(*this).getSubPaving()->getAllNodes(histNodes); 

	vector<double>::iterator vecIt;
	
	//traverse the tree and get the heights 
	for (histNodeIt = histNodes.begin(); histNodeIt < histNodes.end(); 
			histNodeIt++) {
		//get the height at each node
		RangeCollectionClass<real> height((*histNodeIt)->getSummary()[0]);
		heightHist.push_back(height);
	}

	//allocate ranges for histNode
	histMap.allocateRanges(heightHist, 0);

	return nodeEst.getMappedSPIAE(histMap);
}
*/

//gat41
// get delta for each node(or union of nodes). 
double AdaptiveHistogramCollator::getNodesDelta(
set<CollatorSPnode*, less<CollatorSPnode*> > & YatSet, int thisTheta, size_t sizeColl)
{
  // iterator for Yatracos set
  set<CollatorSPnode*, less<CollatorSPnode*> >::iterator YatSetIt;  

  //gloria - think about dotprecision summation
  
  //initialization
  double delta = 0;
  //dotprecision deltaDP = 0;
  
  //go through each node in this set to get delta
  for (YatSetIt = YatSet.begin(); YatSetIt != YatSet.end(); YatSetIt++) {
		//cout << (*YatSetIt)->getNodeName() << endl;
		//cout << "union " << endl;
		delta += (*YatSetIt)->getNodeDelta(thisTheta, sizeColl);
		//accumulate(deltaDP, (*YatSetIt)->getNodeDelta(k, thisTheta), 1.0);
   }
   //cout << "end of union" << endl;

  // take the absolute value of the sums
  //cout << "Delta: " << fabs(delta) << endl;
  return fabs(delta);
}

// Get Scheffe set from sub-pavings.
void AdaptiveHistogramCollator::getHistScheffeSet(
		vector < vector< set<CollatorSPnode*, less<CollatorSPnode*> > > > & vecScheffeSetVec)
{
   //=============end of setting up containers================================// 
	int numAdd = getNumberCollated()-1; // the number of histograms collated including the 0-th histogram
	//cout << "getNumberCollated: " << numAdd << endl;
	size_t theta = numAdd-1; // the current number of splits

  	//============begin pairwise comparisons===================================//
  	for (size_t k=0; k < numAdd; k++) {
		vector< set<CollatorSPnode*, less<CollatorSPnode*> > > vecScheffeSet;
      for (size_t j = 1; j < numAdd; j++) {
			if ( (k != j) && (k < j) ) {
				set<CollatorSPnode*, less < CollatorSPnode* > > currentScheffeSet;
				cout << "k= " << k << "\t" << "theta = " << j << endl;
				getSubPaving()->getScheffeSet(currentScheffeSet, k, j);
				//if (currentScheffeSet.empty()) { cout << "nothing here" << endl; }
				vecScheffeSet.push_back(currentScheffeSet);
			}
		}
		vecScheffeSetVec.push_back(vecScheffeSet);
	} // end of pairwise comparisons
} // end of function getHistScheffeSet

// get delta_theta for all theta
void AdaptiveHistogramCollator::getHistScheffeWinner(
		vector< vector< set<CollatorSPnode*, less<CollatorSPnode*> > > > & vecScheffeSetVec, 
		vector< std::vector<int> > & vecWinnerVec,
		vector< std::vector<double> > & vecDeltaWinnerVec)
{	
	size_t sizeColl = getNumberCollated();
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
				double deltaI = getNodesDelta(vecScheffeSetVec[i][j], cand1, sizeColl);
				//cout << "---------get delta for " << cand2 << endl;
				double deltaJ = getNodesDelta(vecScheffeSetVec[i][j], cand2, sizeColl);

				// perform competition
				if ( deltaI < deltaJ ) {
					//cout << cand1 << "\t" << deltaI << "\n" << cand2 << "\t" << deltaJ << endl;
					cout << "Winner is: " << cand1 << endl;
					// winner is i
					WinnerVec[j] = (1);
					DeltaVec.push_back(deltaI);
				}
				else { // deltaTheta >= delta 
					//cout << cand1 << "\t" << deltaI << "\n" << cand2 << "\t" << deltaJ << endl;
					cout << "Winner is: " << cand2 << endl;
					WinnerVec[j] = (0);
					DeltaVec.push_back(deltaJ);
				}  // end
			} // end of set not empty
		} // end of going through j
		vecWinnerVec.push_back(WinnerVec);
		vecDeltaWinnerVec.push_back(DeltaVec);
	} // end of going through vecScheffeSet

} // end of getHistScheffeWinner


void AdaptiveHistogramCollator::getHistYatSet(
		vector< set<CollatorSPnode*, less<CollatorSPnode*> > > & vecYatSet)
{
	int numAdd = getNumberCollated()-1; 
	// the number of histograms collated including the 0-th histogram
	// but not the last histogram as it is the validation histogram
	//============begin pairwise comparisons===================================//
	for (size_t k= 0; k < numAdd; k++) {
		// get A_ij
      for (size_t j = 0; j < numAdd; j++) {
			if ( (k != j) && (k<j) ) {
				set<CollatorSPnode*, less < CollatorSPnode* > > RowSet;
				set<CollatorSPnode*, less < CollatorSPnode* > > ColSet;
				//cout << "k= " << k << "\t" << "theta = " << j << endl;
				getSubPaving()->getYatSet(RowSet, ColSet, k, j);
				vecYatSet.push_back(RowSet);
				vecYatSet.push_back(ColSet);
			}
		}
	} // end of pairwise comparisons
} // end of function getHistYatSet

void AdaptiveHistogramCollator::getMinDistEst(
																				vector<double> & maxDelta, 	
vector< set<CollatorSPnode*, less<CollatorSPnode*> > > & vecYatSet)
{
	//cout << "calling getmindistest" << endl;
	//get the yatracos class for ALL the candidates
	getHistYatSet(vecYatSet); 
	
	//get the maximum delta at each "theta" - here theta refers to the position of the
	//candidate in the collator
	for (int i = 0; i < getNumberCollated(); i++) {
		//get the maximum delta at this candidate
		double deltaMax = getNodesMaxDelta(vecYatSet, i);
		maxDelta.push_back(deltaMax);
	}
}

// get delta for each node(or union of nodes). 
//think about dotprecision summation
double AdaptiveHistogramCollator::getNodesMaxDelta(
			vector< set<CollatorSPnode*, less<CollatorSPnode*> > > & vecYatSet, 
			int thisTheta)
{
  // iterators  
  vector< set<CollatorSPnode*, less<CollatorSPnode*> > >::iterator YatSetIt;  

	double DeltaMax = 0;
	//dotprecision deltaDP = 0;
	set<CollatorSPnode*, less<CollatorSPnode*> > YatSet;
	
	//go through each node in this set to get delta
	for (YatSetIt = vecYatSet.begin(); YatSetIt < vecYatSet.end(); YatSetIt++){
		double delta = getNodesDelta((*YatSetIt), thisTheta, getNumberCollated());
		//accumulate(deltaDP, (*YatSetIt)->getNodeDelta(k, thisTheta), 1.0);
		delta = fabs(delta);
		DeltaMax = (delta > DeltaMax) ? delta : DeltaMax; 
		//cout << "DeltaMax: " << DeltaMax << endl;
	}

  // take the absolute value of the sums
  return fabs(DeltaMax);
	
}


// ---------- end implementation of AdaptiveHistogramCollators -----------

//Output all boxes in AdaptiveHistogramCollator adhc
std::ostream & operator<<(std::ostream &os,
                    const AdaptiveHistogramCollator& adhc)
{
    if (NULL != adhc.getSubPaving()) {
        os << (adhc.getSubPaving())->nodesAllOutput(os, 1) << std::endl;
    }

    return os;
}

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

/*! \file adaptivehistogram.cpp
\brief AdaptiveHistogram definitions
*/

#include "adaptivehistogram.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <exception> // use exceptions

#include <gsl/gsl_math.h> // to use the constant M_PI 
#include <math.h> // math library
#include <gsl/gsl_rng.h>        // to know about the gsl random number generator
#include <gsl/gsl_randist.h>    // to get the IAE 
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

//to use subpavings
#include "sptools.hpp"
#include "spalgorithms.hpp"

// to use stats subpavings
#include "collatorspnode.hpp"

// to use mcmc function objects
#include "histmcmcobjs.hpp"

// to use histogram penalty function objects
#include "histpenalty.hpp"

// to use histogram evaluation objects
#include "histevalobj.hpp"

// to use collator objects
#include "adaptivehistogramcollator.hpp"

// to use error functions
#include "errorfunc.hpp"

// to use 2D integration using taylor methods
#include "../examples/StatsSubPav/ExactInt/Int.h"
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"

using namespace subpavings;
using namespace std;

// a class for comparison between spsnodes
class MyCompare
{
    const NodeCompObj& myNC;

    public:
    MyCompare(const NodeCompObj& nc) : myNC(nc) {}

    bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
    { return myNC(lhs, rhs); }

};

// -------------------implementation of AdaptiveHistogram class --------------

// --------------------------- private ---------------------------------------

// a constant for padding a box if it is tailor-made for data
const real AdaptiveHistogram::padding = 0.000005;

// initialised constructor, initialised with a subpaving pointer
AdaptiveHistogram::AdaptiveHistogram(SPSnode * spn, bool as)
        : holdAllStats(as),
          scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    if (NULL == spn) {
        throw HistException("Cannot use null SPSnode pointer in constructor");
    }
    rootPaving = spn;
    creationString = rootPaving->getNodeName();
    creationString += rootPaving->getChildNodeNames();

    rootBox = spn->getBox();

    // nothing happens to dataCollection when object is constructed
}


// complete insertion of data from a vector of data
// given a container of rvectors of the data  to insert
bool AdaptiveHistogram::completeDataInsertionFromVec(const RVecData& theData,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging)
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
    size_t dataCountInserted
            = insertDataFromContainer(theData, boolTest, logging);

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
bool AdaptiveHistogram::haveMadePaving(const RVecData& theData,
                                    const size_t dim)
{

    bool retValue = false;

    try {

        // check if we need to make the paving on the basis of the data
        if (isEmpty(rootPaving)) {

            rootBox = makeBox(theData, dim);

            // point rootPaving to a new SPSnode with box myBox
            // and also pass in the not value of holdAllStats which controls
            // whether all available statistics are maintained in the
            // rootPaving (true) or just counts (false)
            rootPaving = new SPSnode(rootBox, !holdAllStats);
            creationString = rootPaving->getNodeName();

            retValue = true;
        }
    }

    catch (bad_alloc& e)
    {
        const char* msg = e.what();
        std::cerr << msg << std::endl;
        std::cerr << "Error allocating memory in "
            << "AdaptiveHistogram::haveMadePaving()"
            << std::endl;
        throw;
    }

    return retValue;
    // end of making the subpaving if there was not one
}




// make a box to fit all the data
ivector AdaptiveHistogram::makeBox(const RVecData& theData, const size_t dim)
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
// recalculates the scaled EMP Sums for COPERR and AIC
// creates log file of process if logging is true.
size_t AdaptiveHistogram::insertDataFromContainer(const RVecData& theData,
                                    const SplitDecisionObj& boolTest,
                                    LOGGING_LEVEL logging)
{
    size_t counter = 0;    // to count the input

    // for logging output to keep track of splits if necessary
    int i = 0;
    std::string baseFileName = "";
    std::string s = "";

    // if we are splitting as we go and logging, set up a log file
    if ((logging != NOLOG) && (boolTest() == true)) {
        baseFileName = "splitOutput";
        s = getUniqueFilename(baseFileName);
        outputLogStart(s);
        // log the current state of the histogram
        outputLog(s, i);
        i++;
    }

    RVecDataCItr cit;

    // feed the data to myHist
    for(cit = theData.begin(); cit < theData.end(); cit++) {

        // put it into dataCollection
        BigDataItr it = dataCollection.end();
        it = dataCollection.insert(it, *cit);

        SPSnode* insertedInto = NULL;

        // try inserting
        insertedInto =
                rootPaving->insertOneFind(it,ON_PARENT, boolTest);

        //insertOneFind returns either NULL if no insert possible
        // or a pointer to the node the data goes to before that node
        // is split (it could be split more than once)
        if (NULL == insertedInto) { // failed to insert
            std::cout << "Failed to insert point "
                << *cit << std::endl;
            std::cout << "Root node of subpaving has box "
                << rootPaving << std::endl;
        }
        // successful insertion, and we are splitting as we go
        else if (boolTest() == true) {
            std::string newNames = insertedInto->getChildNodeNames();

            if(newNames.length() > 0) { // there are new nodes
                //add the new child names if any
                creationString += newNames;

                if (logging) { // log the current state of the histogram
                    outputLog(s, i);
                    i++;
                }
            }
       }

        counter++;
    }

    if (counter > 0) { // data inserted
        //recalculate the scaled EMP sum values;
        recalcScaledEMPSumCOPERR();
        recalcScaledEMPSumAIC();

        if ((logging != NOLOG) && (boolTest() == true))  {
            // add leaf node levels string to log
            outputFile(s, getLeafLevelsString());
        }
    }

    return counter;
}


// Recalculate the scaled EMP part of COPERR score.
void AdaptiveHistogram::recalcScaledEMPSumCOPERR() const
{
    // use the scaled EMP Sum from the root node's getEMPSumCOPERR()
    scaledEMPSumCOPERR = rootPaving->getEMPSumCOPERR(
                        rootPaving->getCounter());
}

// Recalculate the unscaled EMP part of AIC score.
void AdaptiveHistogram::recalcScaledEMPSumAIC() const
{
    // use the scaled EMP Sum from the root node's getEMPSumCOPERR()
    scaledEMPSumAIC = rootPaving->getEMPSumAIC(
                        rootPaving->getCounter());
}

// Update the scaled EMP part COPERR score given change.
void AdaptiveHistogram::updateScaledEMPSumCOPERR(dotprecision change) const
{
    scaledEMPSumCOPERR = scaledEMPSumCOPERR + change;
}

// Update the the scaled EMP part AIC score given change.
void AdaptiveHistogram::updateScaledEMPSumAIC(dotprecision change) const
{
    scaledEMPSumAIC = scaledEMPSumAIC + change;
}



// Method to add current state of the histogram during splitting to a log file
// Output goes to file named according to argument s
// Output includes scaled EMP contributions under COPERR and AIC
// and changes if split
void AdaptiveHistogram::outputLog(const std::string& s, const int i) const
{
    // To add output of the AdaptiveHistogram object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        size_t n = rootPaving->getCounter();

        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        os << creationString << std::endl; // creation string so far
        getSubPaving()->leavesOutputTabsWithEMPs(n, os); // the output
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}


// Method to put opening line into a log file
void AdaptiveHistogram::outputLogStart(const std::string& s) const
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

// log changes in log of MCMC proposal probabilty compared to random (0,1)
void AdaptiveHistogram::logMCMCDeltas(std::string s, int i,
                            real deltaL, real deltaP, real deltaQ, real deltaPi,
                            double randChange)
{
    RealVec vals;
    vals.push_back(deltaL);
    vals.push_back(deltaP);
    vals.push_back(deltaQ);
    vals.push_back(deltaPi+deltaQ);
    vals.push_back(log(randChange));
    vals.push_back(tryExp(deltaL));
    vals.push_back(tryExp(deltaP));
    vals.push_back(tryExp(deltaQ));
    vals.push_back(tryExp(deltaPi+deltaQ));
    vals.push_back(randChange);
    outputFile(s, vals, i);
}


// start a log file for MCMC
// s is the file name
void AdaptiveHistogram::MCMCStartLogFile(std::string s, int i,
                const MCMCProposal& proposal, const LogMCMCPrior& logPrior)
{
    // Start log file with filename and timestamp
    outputLogStart(s);
    // put in the name of the proposal and prior
    std::string line = "Prior is " + logPrior.getName();
    line += ", proposal is " + proposal.getName();
    outputFile(s, line);
    // log the current state of the histogram
    outputLog(s, i);
    // output AIC score information
    outputLogEMPAIC(s);
    // output COPERR score information
    outputLogEMPCOPERR(s);
    std::string headers = "deltaL \t deltaP \t deltaQ \t deltaPi&Q \t ln(rand)";
    headers += "\t ratioL \t ratioP \t ratioQ \t ratioPi&Q \t rand";
    outputFile(s, headers);

}

// output the state of this histogram as an MCMC sample
void AdaptiveHistogram::outputMCMCStateSample(int i)
{
    // create a name for the file to output
    std::string sampleFileName = "MCMCSample";
    //convert i to a string
    std::ostringstream stm;
    stm << i;

    // add the stringed i to the filename
    sampleFileName += stm.str();
    sampleFileName += ".txt"; // and finish the filename

    // To realize a file output
    outputToTxtTabs(sampleFileName);

}

// To add final state of histogram to log file
void AdaptiveHistogram::MCMCLogFinalState(std::string s, int i)
{
    outputLog(s, i);
    // output AIC score information
    outputLogEMPAIC(s);
    // output COPERR score information
    outputLogEMPCOPERR(s);
    // log the leaf levels line
    outputFile(s, getLeafLevelsString());

}



// Method to append COPERR EMP Score values to output log file
// Output goes to file named according to argument s
void AdaptiveHistogram::outputLogEMPCOPERR(const std::string& s) const
{
    // To add output of the AdaptiveHistogram object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {

        real emp = rnd(scaledEMPSumCOPERR);
        os << std::endl << "COPERR EMP is \t" << emp;
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

// Method to append AIC EMP Score values to output log file
// Output goes to file named according to argument s
void AdaptiveHistogram::outputLogEMPAIC(const std::string& s) const
{
    // To add output of the AdaptiveHistogram object to file
    ofstream os(s.c_str(), ios::app);         // append

    if (os.is_open()) {
        real emp = rnd(scaledEMPSumAIC);
        os << std::endl << "AIC EMP is \t" << emp;
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

// Returns an iterator to node to propose for changes
// alters haveNode to true if a proposal node has been found (otherwise the
// iterator just points to the beginning of nodes and should not be used).
SPSnodeListItr AdaptiveHistogram::proposeChangeMCMCState
                        (const MCMCProposal& proposal, SPSnodeList & nodes,
                        size_t numLeaves, size_t numCherries,
                        gsl_rng* rgsl, bool& haveNode)
{
    //cout << "---> calling propose change" << endl;
    
    SPSnodeListItr it;

    RealVec probs;
    RealVecItr pit;

    // fillNodeProposalProbs returns the sum of the probabilities and
    // also fills in probs
    // the sum of the probabilities may be < 1 if for instance we fix the
    // probability of a split or merge in advance
    real psum = proposal.fillNodeProposalProbs(numLeaves, numCherries, probs);

    // check we got back the right number of probabilites
    if (!nodes.empty() && nodes.size() == probs.size()) {

        // pick a node at random  by drawing a random number in [0,1)
        double rand = gsl_rng_uniform(rgsl);
       
        if ((numLeaves > 0) && (numCherries > 0)) {
            
          //  cout << "----> cherry and leaf " << endl;
            
            pit = probs.begin();
            it = nodes.begin();
            real sum = 0.0;

            for (pit = probs.begin(); pit < probs.end(); pit++) {

                    sum += *pit;
                    if (rand < sum) {
                        haveNode = true;
                        break;
                    }
                    it++; // using this to point to nodes relies on having got
                            // the right number of probabilities to match nodes
             }
        }
        // if we only have leaves (which should mean one leaf, no cherries)
        // then we'll only pick a leaf if rand < psum
        if ((numLeaves > 0) && (numCherries == 0) && (rand < psum)) {
            
            
            pit = probs.begin();
            it = nodes.begin();
            real sum = 0.0;

            for (pit = probs.begin(); pit < probs.end(); pit++) {
                    
                    sum += *pit;
                    if (rand < sum) {
                        haveNode = true;
                       break;
                    }
                    it++; // using this to point to nodes relies on having got
                            // the right number of probabilities to match nodes
            }
        }
        // if we only have cherries
        // then we'll only pick a cherry if rand >= 1-psum
        if ((numLeaves == 0) && (numCherries > 0) && (rand >= 1.0-psum)) {
            
            pit = probs.begin();
            it = nodes.begin();
            real sum = 1.0-psum; // note that sum starts at 1-psum

            for (pit = probs.begin(); pit < probs.end(); pit++) {
         
                   
                    sum += *pit;
                    if (rand < sum) {
                        haveNode = true;
                        break;
                    }
                    it++; // using this to point to nodes relies on having got
                            // the right number of probabilities to match nodes
                   
            }
        }
        // it iterator should now point to the node we want to target
        //else we've not taken a node

    }

    return it;
}

// returns true/false decision on whether to split or not
bool AdaptiveHistogram::decisionMCMCSplit(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        size_t numLeaves, size_t numCherries, size_t minPoints,
                        LOGGING_LEVEL logging, const std::string& s, int i) const
{
    bool willSplit = false;

	//cout << "split or not split with minpoints = " << minPoints << endl;
	
    if ((logging == TXT) || (logging == TXTANDGRAPH))
        outputFile(s, "grabbing leaf " + target->getNodeName());

    // change in contribution to log likelihood for this node on split
    real deltaL = rnd(target->getSplitChangeLogLik());

    size_t realNumLeaves = getRootLeaves(); // realNumLeaves is the number of 
														//leaf nodes of the 'old' state

    // use the prior distribution object to find the change in prior
    // here realNumLeaves = #leaves in new state - 1 
    real deltaP = logPrior(realNumLeaves) - logPrior(realNumLeaves - 1);

    // posterior is proportional to likelihood * prior
    real deltaPi = deltaL + deltaP;

    // new numbers of leaves and cherries under proposal depends on minPoints
    // because this determines whether the new leaf children will go into the
    // nodes container
    size_t newNumLeaves = numLeaves - 1; // current number of leaves less this
    // increase the number of new leaves for each new child that can
    // go into the nodes container

    //childrensSpread will be a container of the number of points the children
    // of each child of target might have, in order
    // [0] = left child's left child count, [1] = left child's rght child count,
    // [2] = rght child's left child count, [3] = rght child's rght child count,
    Size_tVec childrensSpread;
    childrensSpread =
                target->getChildrensLeftAndRightCountsIfSplit(childrensSpread);

    if ((childrensSpread[2] >= minPoints) &&
    (childrensSpread[3] >= minPoints)) {
     //   cout << childrensSpread[2] << "\t" << childrensSpread[3] << endl;        
        newNumLeaves++;
    } // this will add one to the leaf numbers if minPoints == 0

    // we will also be prepared to put the right child into the container if
    // there is a minPoints > 0 but one of its
    // children would take all the points, the other getting none
    size_t rightChildCount = childrensSpread[2] + childrensSpread[3];
    if ((minPoints > 0) && (rightChildCount >= minPoints) &&
        ((childrensSpread[2] == 0) || (childrensSpread[3] == 0))) {
        //cout << rightChildCount << endl;
        newNumLeaves++;
    }

    if ((childrensSpread[0] >= minPoints) &&
    (childrensSpread[1] >= minPoints)) {
		 
			//cout << childrensSpread[0] << "\t" << childrensSpread[1] << endl;
        newNumLeaves++;
    }  // this will add one to the leaf numbers if minPoints == 0

    // we would also be prepared to put the left child into the container if
    // there is a minPoints > 0 but one of its
    // children would take all the points, the other getting none
    size_t leftChildCount = childrensSpread[0] + childrensSpread[1];
    if ((minPoints > 0) && (leftChildCount >= minPoints) &&
        ((childrensSpread[0] == 0) || (childrensSpread[1] == 0))) {
        newNumLeaves++;
    }

    size_t newNumCherries = numCherries;
    if (!(target->hasLeafSibling())) newNumCherries = numCherries + 1;

    // Using proposal distribution object
    real deltaQ = proposal.getLogQRatioSplitProposal(numLeaves, numCherries,
                                                newNumLeaves, newNumCherries);
    //get another random number
    double randChange = gsl_rng_uniform(rgsl);
//	cout << log(randChange) << "\t" << (deltaPi + deltaQ) << endl;

    if (log(randChange) < deltaPi + deltaQ) { willSplit = true; }

    if ((logging == TXT) || (logging == TXTANDGRAPH)) { // log these values
        logMCMCDeltas(s, i, deltaL, deltaP, deltaQ, deltaPi, randChange);
        if (willSplit) outputFile(s, "Splitting");
        else outputFile(s, "Not splitting");
    }

    return willSplit;
}

// returns true false decision on whether to merge or not
bool AdaptiveHistogram::decisionMCMCMerge(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        size_t numLeaves, size_t numCherries, size_t minPoints,
                        LOGGING_LEVEL logging, const std::string& s, int i) const
{
    bool willMerge = false;

	//cout << "merge or not merge?" << endl;
	//cout << target->getLeftChild()->getCounter() << "\t" << 
	//target->getRightChild()->getCounter() << endl;

    // cherry so we are merging
    if ((logging == TXT) || (logging == TXTANDGRAPH))
        outputFile(s,"grabbing cherry " + target->getNodeName());

    // change in log likelihood on merge is getMergeChangeLogLik
    // for this node
    real deltaL = rnd(target->getMergeChangeLogLik());

    size_t realNumLeaves = getRootLeaves();

    real deltaP = logPrior(realNumLeaves - 2) - logPrior(realNumLeaves - 1);

    // posterior is proportional to likelihood * prior
    real deltaPi = deltaL + deltaP;

    // calculate the number of leaves and cherries after proposed merge
    // we have to take into account minPoints and the effect that this will have
    // had on whether the target's children are in the nodes container
    size_t newNumLeaves = numLeaves + 1; // current number of leaves plus target
    // but decrement newNumLeaves for each of the target's children that comes
    // out of the container
    if (target->getLeftChild()->getMinChildCountIfSplit() >= minPoints) {
        newNumLeaves--;
    }
    // the left child would also have been in the container if it had enough
    // points and all of them went to one child, the other getting nothing
    if ((minPoints > 0)
        && (target->getLeftChild()->getCounter() >= minPoints)
        && (target->getLeftChild()->getMinChildCountIfSplit() == 0)) {
            newNumLeaves--;
    }

    if (target->getRightChild()->getMinChildCountIfSplit() >= minPoints) {
        newNumLeaves--;
    }
    // the right child would also have been in the container if it had enough
    // points and all of them went to one child, the other getting nothing
    if ((minPoints > 0)
        && (target->getRightChild()->getCounter() >= minPoints)
        && (target->getRightChild()->getMinChildCountIfSplit() == 0)) {
            newNumLeaves--;
    }

    size_t newNumCherries = numCherries;
    if (!(target->hasLeafSibling())) { 
		 
		 //cout << "only remove cherry if it doesn't have a sibling leaf node" << endl;
		 newNumCherries = numCherries - 1; }

    // Using proposal distribution object
    real deltaQ = proposal.getLogQRatioMergeProposal(numLeaves, numCherries,
                                                newNumLeaves, newNumCherries);
    //get another random number
    double randChange = gsl_rng_uniform(rgsl);
    if (log(randChange) < deltaPi + deltaQ) willMerge = true;

    if ((logging == TXT) || (logging == TXTANDGRAPH)) {
        logMCMCDeltas(s, i, deltaL, deltaP, deltaQ, deltaPi, randChange);
        if (willMerge) outputFile(s, "Merging");
        else outputFile(s, "Not merging");
    }

    return willMerge;
}


// change histogram state by splitting the target leaf node
bool AdaptiveHistogram::changeStateForSplit(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries, size_t minPoints)
{
    bool success = true;

    size_t points = rootPaving->getCounter(); // need to recalculate COPERR EMP

    // accumulate the changes in scaled EMP sums that will result
    // from this expansion
   // updateScaledEMPSumCOPERR(target->getSplitChangeEMPCOPERR(points));
   // updateScaledEMPSumAIC(target->getSplitChangeEMPAIC());

    // split the target and divvie up its data
    Expand(target);
    numLeaves--;

    // add the new child names to the creation string
    creationString += target->getChildNodeNames();

    // but only put the children into the container if they can be split, which
    // means if their children would have more than minPoints points in them
    // or if minPoints > 0 and right child has enough points and splitting it would
    // give all its points to one child, none to the other
    if ((target->getRightChild()->getMinChildCountIfSplit() >= minPoints)
        ||
        ((minPoints > 0) && (target->getRightChild()->getCounter() >= minPoints)
        && (target->getRightChild()->getMinChildCountIfSplit() == 0))
        ) {

        // insert the new children ptrs into the list at the beginning
        nodes.push_front(target->getRightChild());
        numLeaves++;
    }

    if ((target->getLeftChild()->getMinChildCountIfSplit() >= minPoints)
        ||
        ((minPoints > 0) && (target->getLeftChild()->getCounter() >= minPoints)
        && (target->getLeftChild()->getMinChildCountIfSplit() == 0))
        ) {
        // insert the new children ptrs into the list at the beginning
        nodes.push_front(target->getLeftChild());// left goes first
        numLeaves++;
    }

    // if sibling was a leaf, take parent out of cherries
    if(target->hasLeafSibling()) {

        // how to find parent? - search the cherries?
        bool foundParent = false;
        std::string nodeParent = target->getParent()->getNodeName();

        SPSnodeListItr git = nodes.begin();
        advance(git, numLeaves); // advance to the cherries
        // break out of loop if we find parent
        SPSnodeListItr it;
        for (it = git ; it != nodes.end(); it++ ) {
            if ((*it)->getNodeName() == nodeParent) {
                nodes.erase(it);
                numCherries--;
                foundParent = true;
                break;
            }
        }
        success = foundParent;

    }

    // put this node ptr into the cherries, ie at end of list
    nodes.push_back(target);

    numCherries++;

    return success;
}

// change histogram state by merging the target cherry node
bool AdaptiveHistogram::changeStateForMerge(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries, size_t minPoints)
{
    bool success = true;

    size_t points = rootPaving->getCounter(); // need to recalculate COPERR EMP

    // accumulate the changes in scaled EMP sums that will result
    // from this expansion
    //updateScaledEMPSumCOPERR(target->getMergeChangeEMPCOPERR(points));
    //updateScaledEMPSumAIC(target->getMergeChangeEMPAIC());

    // subtract the child names from the creation string
    creationString += (" -(" + target->getChildNodeNames() + ")");

    // take the children out of the list of leaves if they are there
    // each child will only be in the list of leaves if splitting that child
    // would give children with at least the minimum number of data points
    // associated with them, or if the child itself has enough points but
    // and splitting would give one child with 0 points were minPoints > 0
    std::string lcName = target->getLeftChild()->getNodeName();
    std::string rcName = target->getRightChild()->getNodeName();
    int lcPos = 0;
    bool foundLeft = false;
    bool foundRight = false;

    SPSnodeListItr it;

    if ((target->getLeftChild()->getMinChildCountIfSplit() >= minPoints)
        || (minPoints > 0 && target->getLeftChild()->getCounter() >= minPoints
            && target->getLeftChild()->getMinChildCountIfSplit() == 0)) {
        std::string lcName = target->getLeftChild()->getNodeName();
        SPSnodeListItr it;
        // break out of loop if we find left child or get to cherries
        for (it=nodes.begin() ; it != nodes.end(); it++ ) {

            if ((*it)->isSubLeaf()) break;
            if ((*it)->getNodeName() == lcName) {
                nodes.erase(it); // can't keep using iterator now
                numLeaves--;
                foundLeft = true;
                break;
            }
            lcPos++;  // gives position at which lc was found
        }
    }
    else foundLeft = true;

    // now try to find right child - could be immediately after left
    if ((target->getRightChild()->getMinChildCountIfSplit() >= minPoints)
        || (minPoints > 0 && target->getRightChild()->getCounter() >= minPoints
            && target->getRightChild()->getMinChildCountIfSplit() == 0)) {

        if (foundLeft) {
            SPSnodeListItr git = nodes.begin();
            advance(git, lcPos);
            // break out of loop if we find right child or get to cherries
            for (it=git ; it != nodes.end(); it++ ) {

                if ((*it)->isSubLeaf()) break;
                if ((*it)->getNodeName() == rcName) {
                    nodes.erase(it);
                    numLeaves--;
                    foundRight = true;
                    break;
                }
            } // just in case right child was before left
            if (!foundRight) {
                // break out of loop if we find right child
                // or get to cherries
                for (it=nodes.begin() ; it != nodes.end(); it++ ) {
                    if ((*it)->isSubLeaf()) break;
                    if ((*it)->getNodeName() == rcName) {
                        nodes.erase(it);
                        numLeaves--;
                        foundRight = true;
                        break;
                    }
                }
            }
        }
    }
    else foundRight = true;

    success = foundRight;

    // merge the target
    target->nodeReabsorbChildren();
    numCherries--;

    // insert the new leaf ptr into the list at the beginning
    nodes.push_front(target);
    numLeaves++;

    // if sibling was a leaf, add parent to cherries, at end
    if(target->hasLeafSibling()) { // returns false if no parent
	   
        nodes.push_back(target->getParent());
        numCherries++;
    }

    return success;
}

bool AdaptiveHistogram::checkNodeCountForSplit(const SPSnode * const spn,
                bool volChecking, double minVol, size_t minChildPoints)
{
    bool retValue = false;
	
	size_t minChildCount = 0;
	size_t counter = spn->getCounter();
	if (spn->isLeaf()) minChildCount = spn->getMinChildCountIfSplit();
	else {
		size_t lcCount = spn->getLeftChild()->getCounter();
		size_t rcCount = counter - lcCount;
		minChildCount = (lcCount < rcCount ? lcCount : rcCount);
	}
	
	if ((!volChecking || (volChecking && (spn->nodeVolume() >= minVol)))
        && ((minChildPoints == 0)
            || (minChildPoints > 0
                &&
                ((counter >= minChildPoints) &&
                    ((minChildCount == 0)
                    ||
                    (minChildCount >= minChildPoints))
                ))
            )
        ) { retValue = true; }

    return retValue;
}

// ----------- histogram public methods

// default constructor
// holdAllStats defaults to false.
AdaptiveHistogram::AdaptiveHistogram()
        : holdAllStats(false), creationString(""),
          scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    rootPaving = NULL;
    rootBox = ivector();    // ivector with length 1 and undefined elements


    // nothing happens to dataCollection when object is constructed
}

// initialised constructor with bool to control whether all stats maintained
// in root paving
AdaptiveHistogram::AdaptiveHistogram(bool as)
        : holdAllStats(as), creationString(""),
          scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    rootPaving = NULL;
    rootBox = ivector();    // ivector with length 1 and undefined elements


    // nothing happens to dataCollection when object is constructed
}

// initialised constructor, initialised with ivector for box
// and with bool to control whether all stats are maintained in root paving.
// (defaults to false which means that only counts are maintained in rootpaving)
AdaptiveHistogram::AdaptiveHistogram(ivector& v, bool as)
        : holdAllStats(as),
          scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    try {
        rootPaving = new SPSnode(v, !as);
        creationString = rootPaving->getNodeName();

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
AdaptiveHistogram::AdaptiveHistogram(const AdaptiveHistogram& other)
        : holdAllStats(other.holdAllStats), rootBox(other.rootBox)
{
    try {
        rootPaving = new SPSnode(*(other.rootPaving));
        creationString = rootPaving->getNodeName();
        creationString += rootPaving->getChildNodeNames();

        //copy dataCollection from other to this
        dataCollection = other.dataCollection;

        other.recalcScaledEMPSumAIC();
        other.recalcScaledEMPSumCOPERR();
        scaledEMPSumCOPERR = other.scaledEMPSumCOPERR;
        scaledEMPSumAIC = other.scaledEMPSumAIC;
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
AdaptiveHistogram&
            AdaptiveHistogram::operator=(const AdaptiveHistogram& rhs)
{
    try {

        // we have to make sure we delete the current paving
        if (NULL != rootPaving) {
            delete rootPaving;
            rootPaving = NULL;
        }

        if (NULL != rhs.rootPaving) {
            rootPaving = new SPSnode(*(rhs.rootPaving));
            creationString = rootPaving->getNodeName();
            creationString += rootPaving->getChildNodeNames();

            //copy dataCollection from other to this
            dataCollection = rhs.dataCollection;

            rhs.recalcScaledEMPSumAIC();
            rhs.recalcScaledEMPSumCOPERR();
            scaledEMPSumCOPERR = rhs.scaledEMPSumCOPERR;
            scaledEMPSumAIC = rhs.scaledEMPSumAIC;
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

// overloading of + operator
AdaptiveHistogram
            AdaptiveHistogram::operator+(const AdaptiveHistogram& rhs)
{
    if (((NULL != rootPaving) && (NULL != rhs.rootPaving)) &&
    ((Ub(rootPaving->getBox()) != Ub(rhs.rootPaving->getBox()))
    || (Lb(rootPaving->getBox()) != Lb(rhs.rootPaving->getBox()))))
        throw HistException("Added histograms have unequal dimensions");

    SPSnode* newRoot = NULL;

    try {


        newRoot = SPSnode::unionTreeStructure(rootPaving,
                                            rhs.rootPaving);
    }
    catch (bad_alloc& ba) {
        string msg(ba.what());
        std::cerr << "Error allocating memory in constructor: original error "
                                    << msg << std::endl;
        throw HistException("Memory allocation rrror in constructor: " + msg);
    }

    bool hold = (getHoldAllStats() && rhs.getHoldAllStats());

    AdaptiveHistogram newHist(newRoot, hold);

    // put all the data from the two histograms into this one.
    RVecData allData;

    allData.reserve( dataCollection.size() + rhs.dataCollection.size() );
    // copy from this dataCollection into allData;
    allData.assign(dataCollection.begin(), dataCollection.end());
    allData.insert(allData.end(), rhs.dataCollection.begin(),
            rhs.dataCollection.end());

    // and put the data into the histogram
    newHist.insertFromRVec(allData, NOLOG);

    return newHist;
}


//Destructor
AdaptiveHistogram::~AdaptiveHistogram()
{
    delete rootPaving;
}


// Return a pointer to the SPSnode this manages.
SPSnode* AdaptiveHistogram::getSubPaving() const
{return rootPaving;}

//src_trunk_0701
int AdaptiveHistogram::getLabel() const
{
	return 0; //this is temporarily for gat41 src
	//return label;
	}

// Gets the mean from the root box of the paving this manages.
rvector AdaptiveHistogram::getRootPavingMean() const
{
    if (!holdAllStats) std::cout << "Note, holdAllStats is false."
                << std::endl;
    return rootPaving->getMean();
}

// Gets variance covariance vector from root box of rootpaving.
RealVec AdaptiveHistogram::getRootPavingVarCovar() const
{
    if (!holdAllStats) std::cout << "Note, holdAllStats is false."
            << std::endl;
    return rootPaving->getVarCovar();
}

// Gets count in the rootpaving in the root paving.
size_t AdaptiveHistogram::getRootCounter() const
{ return rootPaving->getCounter(); }

// Gets number of leaf nodes in the root paving.
size_t AdaptiveHistogram::getRootLeaves() const
{ return spLeaves(rootPaving); }

// Gets the sum of leaf count over volume in root paving.
real AdaptiveHistogram::getRootSumLeafCountOverVol() const
{ return rootPaving->getSumLeafCountOverVol(); }


// get the penalty value
real AdaptiveHistogram::getPENValue(const PenObj& pen,
                                            int deltaLeaf) const
{
    return pen(this, deltaLeaf);
}


// get the EMP part of COPERR score
real AdaptiveHistogram::getEMPScoreCOPERR() const
{
    //recalcScaledEMPSumCOPERR();
    // default cxsc rounding dotprecision rnd_next
    return rnd(scaledEMPSumCOPERR);
}

// Get the EMP part of AIC score.
real AdaptiveHistogram::getEMPScoreAIC() const
{
    //recalcScaledEMPSumAIC();
    // default cxsc rounding dotprecision rnd_next
    return rnd(scaledEMPSumAIC);
}


// get the COPERR score
// verbose = true gives additional output.
real AdaptiveHistogram::getScoreCOPERR(const PenObj& pen,
                                        bool verbose) const
{
    dotprecision temptotal = scaledEMPSumCOPERR;

    if (NULL != rootPaving) {

        size_t counter = rootPaving->getCounter();

        if (verbose) {
            std::cout << "COPERR EMP is " <<
                        rnd(scaledEMPSumCOPERR) << std::endl;
            std::cout << "COPERR penalty is " << getPENValue(pen, 0) << std::endl;
            std::cout << "Total COPERR score is " <<
                        rnd(scaledEMPSumCOPERR) + getPENValue(pen, 0) << std::endl;

        }

       // getPENValue(pen, 0) gives value of PEN under pen
        accumulate(temptotal, getPENValue(pen), 1.0);
        // temptotal now holds scaledEMPSumCOPERR + PEN)
    }

    // default cxsc rounding dotprecision rnd_next
    return rnd(temptotal);
}

// Get the AIC score.
// verbose = true gives additional output.
real AdaptiveHistogram::getScoreAIC(const PenObj& pen,
                                        bool verbose) const
{
    dotprecision temptotal = scaledEMPSumAIC;

    if (NULL != rootPaving) {

        size_t counter = rootPaving->getCounter();

        if (verbose) {
            std::cout << "AIC EMP is " <<
                        rnd(scaledEMPSumAIC) << std::endl;
            std::cout << "AIC penalty is " << getPENValue(pen, 0) << std::endl;
            std::cout << "Total AIC score is " <<
                        rnd(scaledEMPSumAIC) + getPENValue(pen, 0) << std::endl;

        }

        // getPENValue(pen, 0) gives value of PEN under pen
        accumulate(temptotal, getPENValue(pen, 0), 1.0);
        // temptotal now holds scaledEMPSumAIC + PEN)
    }

    // default cxsc rounding dotprecision rnd_next
    return rnd(temptotal);
}

// get the value of the minimum volume for a splittable node.
double AdaptiveHistogram::getMinVol(double minVolB) const
{
    double retValue = 0.0;

    if (NULL != rootPaving) {

        size_t counter = rootPaving->getCounter();
        retValue = minVolB;
       // retValue =  minVolB * log(1.0*counter)*log(1.0*counter)/counter;
    }
    return retValue;
}

// get the value of holdAllStats field.
bool AdaptiveHistogram::getHoldAllStats() const
{
    return holdAllStats;
}

//src_trunk_0701
// get whether this has a subpaving.
bool AdaptiveHistogram::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

//src_trunk_0701
cxsc::ivector AdaptiveHistogram::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
							"AdaptiveHistogram::getRootBox()");
	}
	return getSubPaving()->getBox();
}

// Get a string of the leaf node levels.
std::string AdaptiveHistogram::getLeafLevelsString() const
{
    string retValue = "";
    if (NULL != rootPaving)
        retValue = rootPaving->getLeafNodeLevelsString();

    return retValue;
}


//insert a single data point into the AdaptiveHistogram object
void AdaptiveHistogram::insertOne(rvector newData,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging)
{
    // make sure we have a paving and then try inserting
    if (isEmpty(rootPaving)) {
        throw HistException("Trying to insert to empty or NULL node");
    }

    // check the dimensions
    if((Ub(newData)-Lb(newData)) != (Ub(rootBox) - Lb(rootBox))) {
        throw HistException("Dimensions of data do not match paving");

    }

    // for logging output to keep track of splits
    int i = 0;
    std::string baseFileName = "";
    std::string s = "";

    if (logging != NOLOG) {
        baseFileName = "splitOutput";
        s = getUniqueFilename(baseFileName);
        outputLogStart(s);
        // log the current state of the histogram
        outputLog(s, i);
        i++;
    }


    BigDataItr it = dataCollection.end();
    it = dataCollection.insert(it, newData);

    SPSnode* insertedInto = NULL;

    // try inserting
    insertedInto = rootPaving->insertOneFind(it, ON_PARENT,
                                            boolTest);

    if (insertedInto==NULL) { // failed to insert
        std::cout << "Failed to insert point " << newData << std::endl;
        std::cout << "Root node of subpaving has box " << rootPaving
            << std::endl;
    }
    else { // insertion succeeded

        std::string newNames = insertedInto->getChildNodeNames();

        if(newNames.length() > 0) { // there are new nodes
            //add the new child names if any
            creationString += newNames;
            if (logging != NOLOG) {
                // log the current state of the histogram
                outputLog(s, i);
                i++;
            }
        }

        //recalculate the scaled EMP sum values;
        recalcScaledEMPSumCOPERR();
        recalcScaledEMPSumAIC();
    }
    if (logging != NOLOG) {
        // add leaf node levels string to log
        outputFile(s, getLeafLevelsString());
    }
}


// method to insert one dimensional data from a txt file
// should be able to deal with integers or doubles
bool AdaptiveHistogram::insertOneDimDataFromTxt(const std::string& s,
                                const SplitDecisionObj& boolTest,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
    bool retValue = false;

    try {
        RVecData myDataRvectors; // container for the rvectors we take in

        // try to read in the file
        retValue = readOneDimDataFromTxt(myDataRvectors, s, headerlines);

        if (retValue) {
            // complete the data insertion
            retValue = completeDataInsertionFromVec(myDataRvectors,
                                                    boolTest, logging);
        }
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        throw HistException(msg);
    }
    // recalc of EMP sums takes place further down

    return retValue;

}


/*
// method to insert doubles from a txt file
// use the insert One Dim data method
bool AdaptiveHistogram::insertDoublesFromTxt(const std::string& s,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging)
{
    return insertOneDimDataFromTxt(s, boolTest, logging)
}
*/

// method to insert rvectors from a txt file
bool AdaptiveHistogram::insertRvectorsFromTxt(const std::string& s,
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
                                                    boolTest, logging);
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
bool AdaptiveHistogram::insertFromRVec(const RVecData& rvec,
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
                retValue = completeDataInsertionFromVec(myDataRvectors,
                                                        boolTest, logging);
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

// method to insert a sample of rvectors from a container of rvectors
// this version takes a random number generator
bool AdaptiveHistogram::insertSampleFromRVec(size_t samplesize,
            gsl_rng * rgsl, const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
    bool retValue = false;

    try {

        RVecData myDataRvectors; // container for the rvectors we take in

        // try to sample data from the container, check how many data points found
        size_t numberTaken = getSampleRvectorsFromRVec(myDataRvectors,
                                rgsl, samplesize, rvec);

        // complete the data insertion
        if (numberTaken > 0) {
            /*  Switch on for more output
            // confirm the amount of data taken from the RSSample
            std::cout << "End of taking sample from data from the container: "
            << numberTaken << " data points used for sample" << std::endl;
            */
            retValue = completeDataInsertionFromVec(myDataRvectors,
                                                    boolTest, logging);
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

// method to insert a sample of rvectors from a container of rvectors
// this version takes seed for a random number generator
bool AdaptiveHistogram::insertSampleFromRVec(size_t samplesize,
            int seed, const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {

        const gsl_rng_type * tgsl;

        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();

        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed
        gsl_rng_set (rgsl, seed); // change the seed

        retValue = insertSampleFromRVec(samplesize, rgsl, rvec, boolTest,
            logging);

        gsl_rng_free(rgsl); // free the random number generator

    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }


    return retValue;
}

// method to insert a sample of rvectors from a container of rvectors
// this version will set up a random number generator with default seed
bool AdaptiveHistogram::insertSampleFromRVec(size_t samplesize,
            const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {

        const gsl_rng_type * tgsl;

        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();

        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = insertSampleFromRVec(samplesize, rgsl, rvec, boolTest,
                logging);

        gsl_rng_free(rgsl); // free the random number generator
    }
    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}

// method to insert all rvectors from an RSSample object
bool AdaptiveHistogram::insertFromRSSample(const RSSample& rss,
                                        const SplitDecisionObj& boolTest,
                                        LOGGING_LEVEL logging, int label)
{
    bool retValue = false;

    try {
        RVecData myDataRvectors; // container for the rvectors we take in

        // try to get data from rss.Samples and check how many data points found
        size_t numberFound = getRvectorsFromRSSample(myDataRvectors, rss, label);

        if (numberFound > 0) {
            /*
            // confirm the amount of data taken from the RSSample
            std::cout << "End of taking data from RSSample: "
                << numberFound << " data points with label "
                << label << " found" << std::endl;
            */
            // complete the data insertion
            retValue = completeDataInsertionFromVec(myDataRvectors,
                                                    boolTest, logging);

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

// method to insert a sample of rvectors from an RSSample object
// this version takes a random number generator
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
            gsl_rng * rgsl, const RSSample& rss,
            const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging, int label)
{
    bool retValue = false;

    try {

        RVecData myDataRvectors; // container for the rvectors we take in

        // try to sample data from rss.Samples and check how many data points found
        size_t numberTaken = getSampleRvectorsFromRSSample(myDataRvectors,
                                rgsl, samplesize, rss, label);

        if (numberTaken > 0) {
            /* switch on for more output during histogram creation
            // confirm the amount of data taken from the RSSample
            std::cout << "End of taking sample from data from RSSample: "
                << numberTaken << " data points used for sample" << std::endl;
            */

            retValue = completeDataInsertionFromVec(myDataRvectors,
                                                    boolTest, logging);
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

// method to insert a sample of rvectors from an RSSample object
// this version takes seed for a random number generator
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
            int seed, const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging, int label)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {

        const gsl_rng_type * tgsl;

        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();

        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed
        gsl_rng_set (rgsl, seed); // change the seed

        retValue = insertSampleFromRSSample(samplesize, rgsl, rss, boolTest,
                logging, label);

        gsl_rng_free(rgsl); // free the random number generator
    }
    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}

// method to insert a sample of rvectors from an RSSample object
// this version will set up a random number generator with default seed
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
            const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging, int label)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {

        const gsl_rng_type * tgsl;

        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();

        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = insertSampleFromRSSample(samplesize, rgsl, rss, boolTest,
                logging, label);

        gsl_rng_free(rgsl); // free the random number generator
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory inserting data.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error inserting data.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error inserting data.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}



// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// outputs to a log file if logging is true
// makes its own random number generator
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB, 
										  size_t maxLeafNodes)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplit(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}


// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// outputs to a log file if logging required
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes)
{    
    //cout << minChildPoints << endl;
    //cout << minVolB << endl;
    //cout << maxLeafNodes << endl;
    
    bool cancontinue = false;
    bool TooManyLeaves = false;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {

        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

        cancontinue = (!pq.empty());
		  
        bool bigEnough = cancontinue;
	     TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

			
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // accumulate the changes in scaled EMP sums that will result
            // from this expansion
            //updateScaledEMPSumCOPERR(chosenLargest->getSplitChangeEMPCOPERR(n));
            //updateScaledEMPSumAIC(chosenLargest->getSplitChangeEMPAIC());

            // split the biggest one and divvie up its data

          // cout << "chosenLargest: " << chosenLargest->getNodeName() << "\t" << chosenLargest->getCounter() << endl;
           
           Expand(chosenLargest);

           
           //cout << getLeafLevelsString() << endl;

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
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
				
				
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
		  
	}
         
			   
        if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());

        }

        // EMPSums will have been adjusted during the splitting process
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}

//gat41
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// outputs to a log file if logging is true
// makes its own random number generator
bool AdaptiveHistogram::prioritySplitWithSwitches(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB, 
										  size_t maxLeafNodes, int removeBox)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplitWithSwitches(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes,
												removeBox);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}

//gat41
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// outputs to a log file if logging required
bool AdaptiveHistogram::prioritySplitWithSwitches(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes, int removeBox)
{    
    //cout << minChildPoints << endl;
    //cout << minVolB << endl;
    //cout << maxLeafNodes << endl;

    
    bool cancontinue = false;
    bool TooManyLeaves = false;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {

        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

        cancontinue = (!pq.empty());
		  
        bool bigEnough = cancontinue;
	     TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

			
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        
        size_t temp = 0;
        
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // accumulate the changes in scaled EMP sums that will result
            // from this expansion
            //updateScaledEMPSumCOPERR(chosenLargest->getSplitChangeEMPCOPERR(n));
            //updateScaledEMPSumAIC(chosenLargest->getSplitChangeEMPAIC());

            // split the biggest one and divvie up its data
				//cout << "===============" << endl;
				
				//cout << "chosenLargest: " << chosenLargest->getNodeName() << 
				//"\t" << chosenLargest->getBox() << "\t"
				//<< chosenLargest->getHellingerDist1D() << endl;
					
				//rvector diffMean = chosenLargest->getMean() - chosenLargest->getUniformMean();
				//cout << "mean differences: " << diffMean[1] << endl;

				//RealVec	unifCovar = chosenLargest->getUniformVarCovar();
				//RealVec Covar = chosenLargest->getVarCovar();
				//cout << "variance difference: " << Covar[0] - unifCovar[0] << endl;
				
				Expand(chosenLargest);

				/*
				string fileName = "QueueHist";
				ostringstream stm;
				stm << temp;
				fileName += stm.str();
				fileName += ".txt";
				outputToTxtTabs(fileName);
				cout << "===============" << endl;
				*/
				temp++;
				 
            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();
            
            // only insert the children into the queue if they have more than or
            // equal to removeBox number of points. If true, then
               // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed
            
            if ( ((chosenLargest->getLeftChild())->getCounter() > removeBox) 
						&& (checkNodeCountForSplit(chosenLargest->getLeftChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new left child into the multiset
                pq.insert(chosenLargest->getLeftChild());
            }

            if ( ((chosenLargest->getRightChild())->getCounter() > removeBox) 
						&& (checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) ) {
                // insert the new right child into the multiset
                pq.insert(chosenLargest->getRightChild());
            }

            if (logging != NOLOG) {
                // To add current state of histogram to log file
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
				
				
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
		  
	}
         
			   
        if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());

        }

        // EMPSums will have been adjusted during the splitting process
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}

// method to make a multi-node histogram into one with possibly fewer nodes
// by prioritising which subleaf node to merge first
// keeps merging until the stopTest is satisfied or runs out of subleaves
// outputs to a log file if logging is true
bool AdaptiveHistogram::priorityMerge(const NodeCompObj& compTest,
                        const HistEvalObj& he, LOGGING_LEVEL logging)
{
    bool cancontinue = false;

    if (NULL == rootPaving) {
            throw HistException("No root paving for priorityMerge");
    }

    gsl_rng * rgsl = NULL;

    try {
        size_t n = 0; // number of points in histogram

        int i = 0;
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqMergeOutput";
            s = getUniqueFilename(baseFileName);
        }

        if (rootPaving->isLeaf()) {
            cancontinue = false;
            std::cerr << "Error in priorityMerge: trying to do "
                << "priority merge where the rootPaving "
                << "has no children" << std::endl;
        }
        else cancontinue = true;

        // a multiset for the queue (key values are not necessarily unique)
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        if (cancontinue) {
            n = rootPaving->getCounter(); // number of points in histogram

            if (logging != NOLOG) {
                 // Start log file with filename and timestamp
                outputLogStart(s);
                // log the current state of the histogram
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            SPSnodePtrs subleaves;
            rootPaving->getSubLeaves(subleaves);

            // insert a copy of the current set of subleaves into the multiset
            pq.insert(subleaves.begin(), subleaves.end());

            cancontinue = (pq.size()>0);
        }

        if(!cancontinue) {
            throw HistException("Error in priority merging, aborting merge");
        }

        bool canmerge = cancontinue;

        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // merge until the HistEvalObj he () operator returns true
        while (canmerge && !he(this)) {

            SPSnode* smallest = *(pq.begin ()); // the first smallest in the set

            // find if there are any more equal to smallest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalSmallest;

            equalSmallest = pq.equal_range(smallest); // everything that = smallest
            size_t numberSmallest = pq.count(smallest); // number of =smallest

            if (numberSmallest > 1) {
                // draw a random number in [0,1)
                double rand = gsl_rng_uniform(rgsl);

                real sum = 0.0;

                // random selection of the =largest node to chose
                for (mit=equalSmallest.first; mit!=equalSmallest.second; ++mit) {
                    sum += 1.0/(1.0*numberSmallest);
                    if (rand < sum) {
                        break;
                    }
                }
            }
            else mit = pq.begin ();

            SPSnode* chosenSmallest = *(mit); // the chosen smallest in the set

            updateScaledEMPSumCOPERR(chosenSmallest->getMergeChangeEMPCOPERR(n));
            updateScaledEMPSumAIC(chosenSmallest->getMergeChangeEMPAIC());

            // subtract the child names from the creation string
            creationString += (" -(" + chosenSmallest->getChildNodeNames() + ")");

            // merge the biggest one
            chosenSmallest->nodeReabsorbChildren();
            pq.erase(mit);// take the iterator to chosen smallest out of the set

            // if smallest had a leaf sibling, smallest's parent is now a cherry
            // and should be inserted into the multiset
            if (chosenSmallest->hasLeafSibling()) {

                pq.insert(chosenSmallest->getParent());
           }

            canmerge = (pq.size()>0);

            if (logging != NOLOG) {
                // To add current state of histogram to log file
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }
        }

        if(!canmerge) {
            std::cout << "No more subleaves left to merge" << std::endl;
        }

        if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
        }

        // EMPSums will have been adjusted during the merging process
        if (NULL != rgsl) gsl_rng_free(rgsl);
    }
    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority merge.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error  in priority merge.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException  in priority merge.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error  in priority merge.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }


    return (cancontinue);
}

// method to make a multi-node histogram into a single node histogram
// by merging up to the root box
bool AdaptiveHistogram::mergeUp()
{
    bool cancontinue = false;

    if (NULL == rootPaving) {
        throw HistException("No root paving for mergeUp");

    }

    try {

        if (rootPaving->isLeaf()) {
            cancontinue = false;
            std::cerr << "Nothing to be done - root paving is already a leaf "
                    << std::endl;
        }
        else cancontinue = true;

        SPSnodePtrs subleaves;

        if (cancontinue) {
            rootPaving->getSubLeaves(subleaves); // subleaves contains the subleaves
            cancontinue = (subleaves.size()>0);
        }

        if(!cancontinue) {
            throw HistException("Error in mergeUp getting subleaves, aborting merge");
        }

        bool canmerge = cancontinue;

        // merge until there is only one leaf
        while (canmerge) {

            SPSnode* target = *(subleaves.rbegin ()); // the last in the vector

            // subtract the child names from the creation string
                creationString += (" -(" + target->getChildNodeNames() + ")");

            // merge the biggest one
            target->nodeReabsorbChildren();
            SPSnodePtrsItr it = subleaves.end();
            it--;
            subleaves.erase(it, subleaves.end());// take the last out of the vector

            // if target had a leaf sibling, target's parent is now a cherry
            // and should be inserted into the multiset
            if (target->hasLeafSibling()) {

                subleaves.push_back(target->getParent());
           }

            canmerge = (subleaves.size()>0);
        }

        if(!canmerge) {
            std::cerr << "Merged to root" << std::endl;
        }

        // EMPSums are not adjusted during the merging process
        recalcScaledEMPSumAIC();
        recalcScaledEMPSumCOPERR();
    }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory in mergeUp.  Orginal error: "
                                            + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in mergeUp.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in mergeUp.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in mergeUp.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}

//splits histogram according to string instruction
//returns true if some splitting was achieved
bool AdaptiveHistogram::splitToShape(std::string instruction)
{
    bool success = false;

    // checks:  is there a root paving, is the string properly formed?
    if (NULL == rootPaving) {
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

       success = rootPaving->splitRootToShape(instruction);

        if (success) {
            // update the creation string
            creationString = rootPaving->getNodeName()
                + rootPaving->getChildNodeNames();
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

// Generating MCMC samples from histogram state space.

/* Note that when minPoints > 0, proposals are effectively drawn from set of
leaf and cherry nodes which does not include any leaf which, if split, would
have a child whose number of points is < minPoints.  Thus the implementation
needs to distinguish between the overall state of the tree and the
set of splittable leaf nodes.  */
bool AdaptiveHistogram::MCMC(unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging)
{
	try {

		bool thinning = (thinout > 0);

		std::vector < AdaptiveHistogram > samples;
		
		samples = MCMCsamples(samples, 
						loops, burnin,
						thinout,
						proposal, logPrior,
						minPoints, logging);
		
		bool good = true;
		
		if (thinning && ( loops > 0)
				&& ( samples.size() < (loops - burnin)/thinout + 1) ) {
			good = false;
		}
		else if (thinning) {
			
			// make a collation object, empty at present
			AdaptiveHistogramCollator coll;

			for (size_t i = 0; i < samples.size(); ++i) {
				// output and collate the sample state;
				(samples[i]).outputMCMCStateSample(i);
				coll.addToCollation(samples[i]);
			}
		
		    std::string collFileName = "CollatorMCMC.txt";
            coll.outputToTxtTabs(collFileName); // output the collation to file

            //  Average the sampled histograms
            std::string avFileName = "AverageMCMC.txt";     // provide a filename

            coll.outputAverageToTxtTabs(avFileName);  // output the average to file
        }

		return good;
        
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error  in MCMC:\n " + oldmsg;
        throw HistException(msg);
    }
}

/* Note that when minPoints > 0, proposals are effectively drawn from set of
leaf and cherry nodes which does not include any leaf which, if split, would
have a child whose number of points is < minPoints.  Thus the implementation
needs to distinguish between the overall state of the tree and the
set of splittable leaf nodes.  */
std::vector < AdaptiveHistogram >& AdaptiveHistogram::MCMCsamples(
						std::vector < AdaptiveHistogram >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging)
{
	
    gsl_rng * rgsl = NULL;

    try {
		
		bool cancontinue = false;

		if (NULL == getSubPaving()) {
			throw HistException("No root paving for MCMC histogram");
		}
		
		if (burnin > loops) {
			throw HistException("burnin > loops");
		}
		
		bool thinning = (thinout > 0);
	
		if (thinning) {
			samples.reserve(samples.size() + (loops - burnin)/thinout + 1);
		}

        // for logging
        int i = 0;
        std::string s = "";
        std::string dot = "";
        if ((logging == TXT) || (logging == TXTANDGRAPH)
                || (logging == LOGSAMPLES) || (logging == LOGANDGRAPHSAMPLES)
				|| thinning) {

            // pass to log output to keep track of splits
            std::string baseFileName = "MCMCOutput";
            s = getUniqueFilename(baseFileName);

            if (logging == TXTANDGRAPH) {

                // for dot graph
                baseFileName = "graph";
                std::string suffix = ".dot";
                dot = getUniqueFilename(baseFileName, suffix);
                outputFile(dot, "digraph G {"); // opening line
            }
        }

        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        // check input parameters
        if (loops - burnin < 0) {
            std::cerr << "Try again: number of loops is less than burnin"
                    <<std::endl;
        }
        else if (thinout > 0 && ((loops-burnin)/thinout < 1)) {
            std::cerr << "Try again: thinout value means no samples will be taken"
                    <<std::endl;
        }
        else cancontinue = true;

        // set up a container for the leaf children
        SPSnodePtrs leafVec;
        // set up a container for the subleaf children
        SPSnodePtrs cherryVec;
        SPSnodeList nodes;
         //lists better than vectors for random access removal
        size_t numLeaves = 0;
        size_t numCherries = 0;

        if (cancontinue) {

            // fill the container with the leaf children
            getSubPaving()->getLeaves(leafVec);

            // fill the container with the subleaf children
            getSubPaving()->getSubLeaves(cherryVec);

            numCherries = cherryVec.size();

            if (!leafVec.empty()) {
                // but only put into the container the leaves which, if split,
                // would have at least minPoints data points associated with them

                SPSnodePtrsItr lit;
                for (lit = leafVec.begin(); lit < leafVec.end(); lit++) {
                    if (checkNodeCountForSplit((*lit), false, 0.0, minPoints))
                    {
                        // leaf can go into container
                        nodes.push_back(*lit);
                        numLeaves++;
                    }
                }
            }

            // no need to check on cherries - they can all go in
            if (numCherries > 0)
                nodes.insert(nodes.end(), cherryVec.begin(),cherryVec.end());

            if (nodes.size() == 0) {
                cancontinue = false;
                std::cerr << "No changeable nodes given minPoints = "
                            << minPoints << ". Sorry, aborting MCMC." << std::endl;
            }
        }

        bool goodLoop = cancontinue;

        if (cancontinue && ((logging == TXT) || (logging == TXTANDGRAPH)
                                            || thinning)) {

            MCMCStartLogFile(s, i, proposal, logPrior);
        }

        i++;

        std::string stateNow = "";
        std::string stateAfter = "";

        // make a collation object, empty at present
        //AdaptiveHistogramCollator coll;

        // loop from here conditional on good loop and cancontinue
        while (goodLoop && (loops > 0) ) {

            // capture state now
            stateNow = "\"" + getLeafLevelsString() + "\"";
            loops--;

            // changeMCMCState updates nodes, numLeaves, numCherries, i
            goodLoop = changeMCMCState(nodes, numLeaves, numCherries, proposal,
                        logPrior, minPoints, rgsl, logging, s, i);

            if (goodLoop) {
                if ((logging == TXT) || (logging == TXTANDGRAPH)) {
                    // log the current state of the histogram
                    outputLog(s, i);
                    outputLogEMPAIC(s);
                }

                if (i >= burnin && (logging == TXTANDGRAPH)) {
                    // capture state after split or merge to graph file

                    stateAfter = "\"" + getLeafLevelsString() + "\"";
                    std::string line = "\t " + stateNow + " -> " + stateAfter + ";";
                    outputFile(dot, line);

                    outputGraphDot(); // and make a graph of current state
                }

                if ((numLeaves == 0 && numCherries == 0)) {
                    throw HistException("No more leaves or cherries in MCMC");

                }

                // if we are taking samples take the sample here
                if (goodLoop && thinning && (i >= burnin) &&
                        ((i-burnin)%thinout == 0)) {

                    // output and collate the sample state;
                    //outputMCMCStateSample(i);
                    
                    //coll.addToCollation(*this);
					samples.push_back(*this);
					
                    if (logging == LOGSAMPLES || logging == LOGANDGRAPHSAMPLES) {
                        // log this sample to log file
                        outputLog(s, i);
                        outputLogEMPAIC(s); // add AIC scores
                    }
                    if (logging == LOGANDGRAPHSAMPLES) outputGraphDot();

                }
            }

            i++;
            // back into loop
        }
        // finished loop
        cancontinue = goodLoop;

        if (cancontinue && ((logging == TXT) || (logging == TXTANDGRAPH))) {
            i--; // this histogram is full state from last pass
            // To add final state of histogram to log file
            MCMCLogFinalState(s, i);
            if (logging == TXTANDGRAPH) {
                // close the dot graph of the changes process
                outputFile(dot, "}"); // closing line

                // make the graph image
                makeDotImage(dot);
            }
        }

		/*
        if (cancontinue && thinning) {
            std::string collFileName = "CollatorMCMC.txt";
            coll.outputToTxtTabs(collFileName); // output the collation to file

            //  Average the sampled histograms
            std::string avFileName = "AverageMCMC.txt";     // provide a filename

            coll.outputAverageToTxtTabs(avFileName);  // output the average to file
        }
		*/
        // free the random number generator
        
        gsl_rng_free (rgsl);
        
        if (!cancontinue) { //empty out the samples if loop failed
			std::vector< AdaptiveHistogram > tmp;
			tmp.swap(samples);
		}
        
        return samples;
    }
    
    catch (exception& e) {
		try {
			if (NULL != rgsl) gsl_rng_free(rgsl); 
			// free the random number generator
		}
		catch (exception& ee) {} // catch and swallow
		
        throw HistException("Error  in MCMCsamples:\n" 
							+ string(e.what()));
    }

}

void AdaptiveHistogram::publicOutputMCMCStateSample(int ci, int i,
                                bool confirm)
{
    // create a name for the file to output
    std::ostringstream stm;
    stm << "MCMC_" << ci << "_Sample" << i << ".txt";

    string sampleFileName = stm.str();

    // To realize a file output
    outputToTxtTabs(sampleFileName, confirm);
}

void AdaptiveHistogram::publicLogMCMCSample(std::string s, int i)
{
    outputLog(s, i);
    outputLogEMPAIC(s); // add AIC scores
}

// Changes the state of this Adapative Histogram using MCMC
// nodes, numLeaves, numCherries are passed by reference

/* Note that the container of splittable leaf nodes and cherries is
maintained separately from the overall state of the tree to save having to
repeatedly assess whether a leaf can be split (given minPoints).  The
number of splittable leaves and cherries is maintained for convenience to
avoid repeatedly counting nodes of different types in the container. */
bool AdaptiveHistogram::changeMCMCState (SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
                        gsl_rng* rgsl, LOGGING_LEVEL logging,
                        std::string s, int i)
{
    bool success = true; // start by assuming that all will be well...
    //bool success = false;
	
	//cout << "change mcmc state called for " << endl;
	//cout << "node list size: " << nodes.size() << endl;
	//cout << "#leaves: " << numLeaves << "\t" << "#cherries: " << numCherries << endl;

    try {

			double rand = gsl_rng_uniform(rgsl);
        
        bool haveNode = false;

        // use proposal to fill the proposal probabilities

        // this changes haveNode as well
        //cout << "proposing" << endl;
        SPSnodeListItr it = proposeChangeMCMCState (proposal, nodes,
                                numLeaves, numCherries,
                                rgsl, haveNode);
			
			/*if (haveNode) 
			{ 
				cout << "a node is proposed"  << (*it)->getNodeName() << endl;
			}
			else {cout << "nothing proposed. stay." << endl;}
			*/
			
        // only do more if haveNode is true, which means that it points to something
        if (haveNode && (*it)->isLeaf()) {

            //grab the leaf
            SPSnode* target = *it; // don't change where it points until erased


            // leaf so we are splitting
            bool willSplit = decisionMCMCSplit(target, proposal, logPrior, rgsl,
                            numLeaves, numCherries, minPoints,
                            logging, s, i);            

            if (willSplit) {
                // take the target out of the list
					//cout << "split and remove from list " << (*it)->getNodeName() << endl;
               nodes.erase(it);

                // try to change the state according to the proposed split
                // nodes, numLeaves and numCherries are passed by reference
                success = changeStateForSplit(target, nodes,
                            numLeaves, numCherries, minPoints);
            } // end of willSplit

        } // end of isLeaf

        else if (haveNode && (*it)->isSubLeaf()) {

            // grab the cherry
            SPSnode* target = *it; // don't change where it points until deleted

            bool willMerge = decisionMCMCMerge(target, proposal, logPrior, rgsl,
                            numLeaves, numCherries, minPoints, logging, s, i);

            if (willMerge) {
                // take the target out of the list of cherries
					//cout << "merge and remove from list " << (*it)->getNodeName() << endl;
               nodes.erase(it);

               success = changeStateForMerge(target,
                            nodes, numLeaves, numCherries, minPoints);

            } // end willMerge
        } // end if cherry


        if (((logging == TXT) || (logging == TXTANDGRAPH)) && !haveNode) {

            std::string line = "No node grabbed (possible if proposal has fixed ";
            line += "probability of split): state stays the same";

            outputFile(s, line, true); // append to log file
        }

        if (success) {
            // is it worth doing this?  on balance I think yes
            success = (nodes.size() == (numLeaves+numCherries));
            //cout << "successful state change/stay." << endl;
				//cout << "---------------------------------" << endl;
        }
        
        else {
            throw HistException("Nodes muddled in changeMCMCstate");

        }
    }
    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory in MCMC.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error  in MCMC.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException  in MCMC.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error  in MCMC.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return success;
}


// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec AdaptiveHistogram::getLeafLevels() const
{
    IntVec levels; // empty container

    if (NULL != rootPaving) {
        rootPaving->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// returns a vector of leaf counts
// left to right
Size_tVec AdaptiveHistogram::getLeafCounts() const
{
    Size_tVec counts; // empty container
    if (NULL != rootPaving) {
        rootPaving->getLeafNodeCounts(counts);
        //levels has now been filled in
    }
    return counts;
}


// make a .dot file for the histogram
bool AdaptiveHistogram::outputGraphDot() const
{
    bool success = false;

    if (NULL != rootPaving) {
        success = rootPaving->outputGraphDot();

    }
    else {
        std::cerr << "Sorry, you can't make a graph without a root paving"
                << std::endl;
    }
    return success;
}


// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void AdaptiveHistogram::outputToTxtTabs(const std::string& s,
                            bool confirm) const
{
    if (NULL != rootPaving) {

        // To generate a file output of the AdaptiveHistogram object
        ofstream os(s.c_str());         // Filename, c-string version
        if (os.is_open()) {

            getSubPaving()->leavesOutputTabs(os); // the output
            if (confirm)
                std::cout << "The output of the AdaptiveHistogram "
                    << "has been written to " << s << std::endl << std::endl;
        }
        else {
            std::cerr << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
// Output includes scaled EMP contributions under COPERR and AIC
// and changes if split
void AdaptiveHistogram::outputToTxtTabsWithEMPs(const std::string& s,
                                                    bool confirm) const
{
    if (NULL != rootPaving) {

        // To generate a file output of the AdaptiveHistogram object
        ofstream os(s.c_str());         // Filename, c-string version
        if (os.is_open()) {
            size_t n = rootPaving->getCounter();
            getSubPaving()->leavesOutputTabsWithEMPs(n, os); // the output
            if (confirm)
                std::cout << "The output of the AdaptiveHistogram with scaled EMPs "
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
void AdaptiveHistogram::outputRootToTxt(const std::string& s,
                                            bool confirm) const
{
    if (NULL != rootPaving) {

        // To generate a file output of root node of the AdaptiveHistogram
        ofstream os(s.c_str());         // Filename, c-string version
        getSubPaving()->nodePrint(os); // the output
        if (confirm)
            std::cout << "Details of the root paving of the AdaptiveHistogram "
                << "has been written to " << s << std::endl << std::endl;
    }

}

//src_trunk_0701
void AdaptiveHistogram::reshapeToUnion(
					const PiecewiseConstantFunction& other)
{
	if ( !hasSubPaving() || !other.hasSubPaving() ) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::reshapeToUnion(const PiecewiseConstantFunction&)");
	}
	
	getSubPaving()->reshapeToUnion(other.getCopySubPaving());
}

void AdaptiveHistogram::reshapeToUnion(
					const PiecewiseConstantFunction& other,
					size_t minChildPoints)
{
	if ( !hasSubPaving() || !other.hasSubPaving() ) {
		throw NullSubpavingPointer_Error(
		"AdaptiveHistogram::reshapeToUnion(const PiecewiseConstantFunction&, size_t)");
	}
	
	getSubPaving()->reshapeToUnion(other.getCopySubPaving(), minChildPoints);
}
//--src_trunk_0701


/*! Clear the histogram's data and counter
 */
void AdaptiveHistogram::makeEmpty()
{
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
	SPSnodePtrsItr it; // and an iterator over the container
   getSubPaving()->getLeaves(leaves); // fill the container
   for(it = leaves.begin(); it < leaves.end(); it++) {
		(*it)->makeEmptyNode();
	}
} 

/*! Distribution-free Likelihood Estimation
*/
real AdaptiveHistogram::getEstLogLikelihoodFromRSSample(
								RSSample& labSampledData, double dx, double wt,
								double WeightHist,
						std::map<rvector, double, std::less<rvector> >& WeightsPM)
{
		dotprecision dpEstLogLik;
		dpEstLogLik = 0.0;
			
		//Get log-likelihoods for points in model 1, 
		//and points both in model 0 and 1
		
		//first multiply the box heights with the weight corresponding to the
		// histogram model
		vector<double> fhat; //container to store heights
		SPSnodePtrs leaves; // set up empty container for leaf node pointers
      SPSnodePtrsItr it; // and an iterator over the container
      getSubPaving()->getLeaves(leaves); // fill the container
      size_t N = getRootCounter();
		// add wt to the heights at each box
		for(it = leaves.begin(); it < leaves.end(); it++) {
			fhat.push_back(WeightHist*
			(((1-wt)*(*it)->getCounter()/(N*1.0)/(*it)->nodeVolume()) + wt));
		}
		
		
		// clear the current histogram's data and counters and 
		// insert new data into the histogram 
		AdaptiveHistogram tempHist = (*this);
		tempHist.makeEmpty();	
		
		//cout << "insert from model 1" << endl;		
		// insert data from model 1 into the empty hist	
	   bool hasData = false;
		hasData = tempHist.insertFromRSSample(labSampledData, NOLOG, 1);
		if (hasData) {
			//now get the estimated likelihood
			size_t pos = 0;
			SPSnodePtrs leavesTemp; // set up empty container for leaf node pointers
			tempHist.getSubPaving()->getLeaves(leavesTemp); // fill the container
			for(it = leavesTemp.begin(); it < leavesTemp.end(); it++) {
				if ((*it)->getCounter() != 0) {       
					//cout << fhat[pos] << "\t" << (*it)->getCounter() <<"\t" << dx*fhat[pos] << "\t" << ((*it)->getCounter())*log(dx*fhat[pos]) << endl;
					accumulate(dpEstLogLik, ((*it)->getCounter())*log(dx*fhat[pos]), 1);
					//cout << dpEstLogLik << endl;
				}
				pos++; 
			}   
		} // end of if hasData for label 1
		
		//make tempHist empty again
		tempHist.makeEmpty();
		//insert data from model 0
		//cout << "inserting from model 0" << endl;
		hasData = false;
		hasData =tempHist.insertFromRSSample(labSampledData, NOLOG, 0);
	   if (hasData) {
		   size_t pos = 0;
			SPSnodePtrs leavesTemp; // set up empty container for leaf node pointers
			tempHist.getSubPaving()->getLeaves(leavesTemp); // fill the container
			for(it = leavesTemp.begin(); it < leavesTemp.end(); it++) {
				if ((*it)->getCounter() != 0) {       
					//get the node's data
					NodeData nodeData = (*it)->getData();
					//go through each data in node
					NodeDataItr dit;
					for (dit = nodeData.begin(); dit != nodeData.end(); dit++){
						BigDataItr bigIt = *dit; 
                  rvector theData = *bigIt;  // convert NodeData to rvector
						
						//get the EMF
					   if ( WeightsPM[theData] != 0 ) {
						//the check for WeightsPM[theData] is needed because there 
						//may be point mass in the sampled data but not in the 
						//histogram that you are observing the data with
							//cout << (*it)->getCounter() << "\t" << (*bigIt) << endl;
							//cout << dx*fhat[pos] << "\t" << WeightsPM[theData] << endl; 
							accumulate(dpEstLogLik, log(dx*fhat[pos] + WeightsPM[theData]), 1);
							//cout << dpEstLogLik << endl;
						} // end of weight check
						else {
							//cout << (*it)->getCounter() << "\t" << (*bigIt) << endl;
							//cout << dx*fhat[pos] << endl; 
							accumulate(dpEstLogLik, log(dx*fhat[pos]), 1);
							//cout << dpEstLogLik << endl;
						}
					} // end of going through node's data			
				} // if counter < 0
				pos++; 
			} // end of iterating over leaves			
		} // end of if hasData for label 0   
				
		real estLogLik = rnd(dpEstLogLik);
		return estLogLik;
}


/*! Get the IAE of the corresponding distribution based on function arguments.
*/
real AdaptiveHistogram::getIAE(int distr)
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

/*! Get the IAE for a unform mixture 
*/
real AdaptiveHistogram::getUnifIAE(AdaptiveHistogram & myPart, double weight, vector<int> holesLoc)
{
	 // get the true height, f of the corresponding box in myPart
	SPSnodePtrs trueLeaves;
	SPSnodePtrsItr trueIt;
	//AdaptiveHistogram * adhPtr;
	//adhPtr = &myPart;
	(myPart).getSubPaving()->getLeaves(trueLeaves);

	// setting up containers for the leaves
	SPSnodePtrs leaves;
	SPSnodePtrsItr it;
	(*this).getSubPaving()->getLeaves(leaves); // fill the container

	double trueF; //true density
	ivector temp;
	
	dotprecision dpIAE;    // use type dotprecision for summation  
	dpIAE=0.0;

	int n = getSubPaving()->getRootCounter();
	//go through all the leaves in this
	for(it = leaves.begin(); it < leaves.end(); it++) {
		ivector thisBox = (*it)->getBox();
		//cout << "====checking " << (*it)->getBox() << endl;
      
		// get the height of this leaf
		double fhat = (*it)->getCounter()/(n*1.0)/(*it)->nodeVolume(); 
		
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
   //cout << "IAE: " << unifIAE << endl;
	return unifIAE;														  

}

/*! Get the IAE for the unform distribution
*/
real AdaptiveHistogram::getUnifIAE()
{
  	//-------setting up containers-------------------------------
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
   SPSnodePtrsItr it; // and an iterator over the container
   (*this).getSubPaving()->getLeaves(leaves); // fill the container		
	int nSample = rootPaving->getCounter(); // get total number of training 
	                                         // data
   dotprecision dpIAE;    // use type dotprecision for summation  
   dpIAE=0.0;
	double f = 1;
   //go through all the leaves in this
   for(it = leaves.begin(); it < leaves.end(); it++) {
      // get the height of this leaf
      //cout << (*it)->getCounter() << "\t" << (*it)->nodeVolume() << endl;
      double fhat = (((*it)->getCounter())*1.0)/((*it)->nodeVolume())/
		               (nSample*1.0);

      //now calculate the IAE
		if ((f - fhat) < 0.0){
			real r = ((*it)->nodeVolume())*(fhat - f);
			//cout << "r: " << r << endl;
			accumulate(dpIAE, r, 1.0);
		}

		else if ((f - fhat) > 0.0){
			real r = ((*it)->nodeVolume())*(f - fhat);
			//cout << "r: " << r << endl;
			accumulate(dpIAE, r, 1.0);
		}

	} // end of going through all the leaves in this

   //cast dotprecision to real
   real unifIAE = rnd(dpIAE);
	return unifIAE;
}


/*! Get the IAE for a finite mixture distribution
*/
real AdaptiveHistogram::getFinMixIAE(FinMix& mixt)
{
	cout << "GEtting IAE for Finite Mixture: " << endl; 

	//---------fill in containers for names, volumes and heights of leaves
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
	SPSnodePtrsItr it; // and an iterator over the container
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
	int n = rootPaving->getCounter();
	          
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
	
	/*
	// Accounting for the boundaries
	rvector xuppVec1 = Sup(boxes[Nbin-1]);
	double xupp1 = _double(xuppVec1[1]);
	rvector xlowVec1 = Inf(boxes[0]);
	double xlow1 = _double(xlowVec1[1]);
	dpIAEBoun = dpFinMixIAEBoun(xlow1, xupp1, mixt);
	dpIAE += dpIAEBoun;
	*/
	
	// cast dot precision to real
	real FinMixIAE = rnd(dpIAE);
	
	// free the workspace
	gsl_integration_workspace_free (w);
	
	return FinMixIAE;
}

// Get the IAE for a finite gaussian mixture distribution using interval 
// techniques.
cxsc::interval AdaptiveHistogram::getFinMixIntervalIAE(FinMix& mixt, double tol, int deg)
{
	interval totalArea(0.0); //initialize
	int n = getRootCounter();

	// need to iterate through the leaves
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
	SPSnodePtrsItr it; // and an iterator over the container
	getSubPaving()->getLeaves(leaves); // fill the container
	
	// container is filled by reading leaves off tree from left to right
	for(it = leaves.begin(); it < leaves.end(); it++) {
		//cout << "-----------------" << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;
		
		//get the height in this leaf node
		double fhat = (*it)->getCounter()/(*it)->nodeVolume()/n;
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
real AdaptiveHistogram::getMappedFunctionIAE(RealMappedSPnode nodeEst)
{
	ivector thisBox = getSubPaving()->getBox();
	RealMappedSPnode histMap(thisBox);
	int n = getRootCounter();

	// split the root box into the shape of myHist
	string leafLevelString = getLeafLevelsString();
	int depth = atoi(leafLevelString.c_str());
	if (depth != 0) {
		histMap.splitToShape(leafLevelString); 
	}

	//container to store heights for histNodes 
	vector< RangeCollectionClass<real> > heightHist;
	//get all the nodes in the histogram 
	SPSnodePtrs histNodes;
	SPSnodePtrsItr histNodeIt;
	getSubPaving()->getAllNodes(histNodes); 

	//traverse the tree and get the heights 
	//cout << "get the height at each node" << endl;
	for (histNodeIt = histNodes.begin(); histNodeIt < histNodes.end(); 
			histNodeIt++) {
		//cout << "IAE at node " << (*histNodeIt)->getNodeName() << endl;
		RangeCollectionClass<real> height((*histNodeIt)->getCounter()/
							((*histNodeIt)->nodeVolume()*n));
		
		//cout << (*histNodeIt)->getCounter()/
		//					((*histNodeIt)->nodeVolume()*n) << endl;
		
		heightHist.push_back(height);
	} // end of traversing all nodes in histogram
	
	//allocate ranges for histNode
	histMap.allocateRanges(heightHist, 0);
	return nodeEst.getMappedSPIAE(histMap);
}
*/

// Get the IAE for a laplace distribution using interval 
// techniques.
cxsc::interval AdaptiveHistogram::getLaplaceIntervalIAE(double tol, int deg)
{
	interval totalArea(0.0); //initialize
	int n = getRootCounter();

	// need to iterate through the leaves
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
	SPSnodePtrsItr it; // and an iterator over the container
	getSubPaving()->getLeaves(leaves); // fill the container
	
	// container is filled by reading leaves off tree from left to right
	for(it = leaves.begin(); it < leaves.end(); it++) {
		cout << "-----------------" << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;
		
		//get the height in this leaf node
		double fhat = (*it)->getCounter()/(*it)->nodeVolume()/n;
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
		cout << "finding roots at this node " << thisInt << endl;
		cout << "bisection" << endl;
		LaplaceBisect(thisIntW, tol, fhat, rootVec); 

		cout << "get L1 error" << endl;
		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = LaplaceGetL1error(fhat, thisInt, deg, tol);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				cout << "root " << i << ": " << rootVec[i] << endl;
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
			cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;

			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
				cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = LaplaceGetL1error(fhat, thisInt, deg, tol);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntLeft, deg, tol);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = LaplaceGetL1error(fhat, thisSubIntRight, deg, tol);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				cout << "let's have a look at all the roots:" << endl;
				for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					cout << uniqueRootVec[i] << endl;
				}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty

	} // end of iterating through the leaf nodes
	
	cout << "IAE: " << totalArea << endl;
	return totalArea;
}

// Get the IAE for a lognormal distribution using interval 
// techniques.
cxsc::interval AdaptiveHistogram::getLognormalIntervalIAE(double tol, int deg)
{
	interval totalArea(0.0); //initialize
	int n = getRootCounter();

	// need to iterate through the leaves
	SPSnodePtrs leaves; // set up empty container for leaf node pointers
	SPSnodePtrsItr it; // and an iterator over the container
	getSubPaving()->getLeaves(leaves); // fill the container
	
	// container is filled by reading leaves off tree from left to right
	for(it = leaves.begin(); it < leaves.end(); it++) {
		//cout << "-----------------" << endl;
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;
		
		//get the height in this leaf node
		double fhat = (*it)->getCounter()/(*it)->nodeVolume()/n;
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
		LognormalBisect(thisIntW, tol, fhat, rootVec); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = LognormalGetL1error(fhat, thisInt, deg, tol);
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
					interval diffArea = LognormalGetL1error(fhat, thisInt, deg, tol);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = LognormalGetL1error(fhat, thisSubIntLeft, deg, tol);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = LognormalGetL1error(fhat, thisSubIntRight, deg, tol);
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
					interval diffArea = LognormalGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea =LognormalGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea =LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty

	} // end of iterating through the leaf nodes
	
	//cout << "IAE: " << totalArea << endl;
	return totalArea;
}


/*! Get the IAE for 2D distribution
*/
real AdaptiveHistogram::get2DIAE(taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval))
{
  //number of points
  int n = rootPaving->getCounter();
	 
  SPSnodePtrs leaves; // set up empty container for leaf node pointers
  SPSnodePtrsItr it; // and an iterator over the container
  (*this).getSubPaving()->getLeaves(leaves); // fill the container
  real tol=1e-6;
  int o=16;
  real result = 0; 
  for (it=leaves.begin(); it < leaves.end(); it++)
  {	
      //get domain
	  ivector domain = (*it)->getBox();
	  cout << "Get IAE for node " << domain << endl;
	  //get fhat	
	  cout << "fhat:" << endl;
     interval fhat = interval(real((*it)->getCounter()/
	                           (((*it)->nodeVolume())*1.0*n))); 
	   
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

void AdaptiveHistogram::findDensityRegion(double cov, double weightPM,
														vector<SPSnode*> & covNodes,
														string covFileName)
{
	try {
			if ( (cov - weightPM) <= 0) {
				cout << cov << " percent of the mass are already covered by the point masses" << endl;
			}
			else {
				// put the leaves into a vector and sort it, smallest to largest
				SPSnodePtrs leaves;
				getSubPaving()->getLeaves(leaves);
				CompHeight compheight;
				//sort according to average height
				sort(leaves.begin(), leaves.end(), MyCompare(compheight));
				
				//start iterating from the largest
				SPSnodePtrs::reverse_iterator rit = leaves.rbegin();
				bool found = FALSE; //found the boxes that gives cov density region
				
				dotprecision totalCov;
				totalCov = 0.0;
				
				size_t totalN = getRootCounter();
				
				while (!found && rit < leaves.rend()) {
					//accumulate the height * box vo\l 
					accumulate(totalCov, (1.0*(*rit)->getCounter())/(1.0*totalN), 1); 

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
				SPSnodePtrsItr vit;
				for (vit = covNodes.begin(); vit < covNodes.end(); vit++) {
					ivector thisBox = (*vit)->getBox(); // copy theBox         
					double vol = (*vit)->nodeVolume();
					// output the nodeName, nodeVolume
					os << (*vit)->getNodeName();
					os << "\t" << vol;
					// followed by the height
					os << "\t" << (1.0*(*vit)->getCounter())/(1.0*totalN)/vol;
					// followed by intervals making up box using Inf & Sup
					// ie unlike cxsc output, there is no [  ] around them
					for (int i= Lb(thisBox); i <= Ub(thisBox) ; i++) {
						 os << "\t" << Inf(thisBox[i]) << "\t" << Sup(thisBox[i]);
					}
					os << endl;
				}
				os << flush;
				os.close();
			}	//end of going through the continuous part
		} // end of try			
	catch (exception& e) {
		throw HistException(
		"Error in AdaptiveHistogram::findDensityRegion :\n"
		+ string( e.what() ) );
	}
}

// Jenny addition for Gloria's convergence work
// Method to add current state of the histogram during splitting to a log file
// Output goes to file named according to argument s
// Output is plain, just textToTabs
void AdaptiveHistogram::outputLogPlain(const std::string& s, const int i) const
{
    // To add output of the AdaptiveHistogram object to file
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
 
// gloria addition 
/*! \brief Returns the IAE between an AdaptiveHistogram object and a RealMappedSPnode.
*/
/*
cxsc::real AdaptiveHistogram::getMappedIAE(RealMappedSPnode& nodeEst, ivector pavingBox) const
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
	SPSnodePtrs histNodes;
	SPSnodePtrsItr histNodeIt;
	getSubPaving()->getAllNodes(histNodes); 

	size_t n = getSubPaving()->getRootCounter(); 
	//traverse the tree and get the heights 
	for (histNodeIt = histNodes.begin(); histNodeIt < histNodes.end(); 
			histNodeIt++) {
		//get the height at each node
		RangeCollectionClass<real> height((*histNodeIt)->getCounter()/
										((*histNodeIt)->nodeVolume()*n));
		heightHist.push_back(height);
	}

	//allocate ranges for histNode
	histMap.allocateRanges(heightHist, 0);

	return nodeEst.getMappedSPIAE(histMap);
}
*/

// gloria's addition
bool AdaptiveHistogram::prioritySplitWithTotalVar(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB, 
										  size_t maxLeafNodes, int StopVal,
										  vector<AdaptiveHistogram>& HistAtValley, int simNum)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplitWithTotalVar(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes, 
                                    StopVal, HistAtValley, simNum);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}


// gloria addition
// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until the function object he returns true
// or until there are no more splittable nodes
// outputs to a log file if logging required
bool AdaptiveHistogram::prioritySplitWithTotalVar(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes, int StopVal,
                                vector<AdaptiveHistogram>& HistAtValley, int simNum)
{    
    bool cancontinue = false;
    bool TooManyLeaves = false;
    bool shouldStop = false;
    int flagStop = 0;
    vector<double> TotalVarDist;
    int Prev = 1;
    AdaptiveHistogram probValley;
    size_t split = 0;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {
        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

        cancontinue = (!pq.empty());
		  
        bool bigEnough = cancontinue;
	     TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // accumulate the changes in scaled EMP sums that will result
            // from this expansion
            //updateScaledEMPSumCOPERR(chosenLargest->getSplitChangeEMPCOPERR(n));
            //updateScaledEMPSumAIC(chosenLargest->getSplitChangeEMPAIC());

            // split the biggest one and divide up its data
           cout << "===============" << endl;
           cout << "chosenLArgest: " << chosenLargest->getNodeName() << "\t" << chosenLargest->getCounter() << endl;
           Expand(chosenLargest);
           cout << "===============" << endl;
           
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
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }


				// get the total variation distance
				SPSnodePtrs leaves;
				SPSnodePtrsItr leavesIt;
				getSubPaving()->getLeaves(leaves);

				double totalVarDist = 0;
				double FUnif = 1.0/leaves.size()*1.0;
				for (leavesIt = leaves.begin(); leavesIt < leaves.end(); leavesIt++) {
					double leafVol = (*leavesIt)->nodeVolume();
					//get the total variation distance
					//calculate \mu_n - \mu
               double fhat = ((*leavesIt)->getCounter())/leafVol/n;
               double diffMu = fabs(fhat*leafVol - FUnif);
               //cout << "previous: " << totalVarDist << "\t current: " << diffMu << endl;
					totalVarDist += diffMu;
					//totalVarDist = (diffMu > totalVarDist) ? diffMu : totalVarDist; 
				}
				split++;
				TotalVarDist.push_back(totalVarDist);
				//cout << "---------Split " << split << ": " << totalVarDist << "--------------" << endl;

				if (split == 1) { probValley = *this; } // keep the second state

				// start the checks after 1 split
				if ( split > 1 ) {
					// use the total var distance as stopping criteria
					size_t vecSize = TotalVarDist.size();
					shouldStop = checkStopCrit(TotalVarDist[vecSize-1], TotalVarDist[vecSize-2], Prev);
					if (shouldStop) { 
						flagStop++; 
						HistAtValley.push_back(probValley); //keep this histogram
					}
					if (flagStop == StopVal) { 
						cout << "Stopping criteria met. There are " << StopVal << " valleys." << endl;
						break; 
					}
					if (Prev == 1) // keep this histogram if prev = 1 
						{ probValley = *this; }
				}
				
            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;

				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}  
			} // end of while loop

			//output the total variation distance for plotting purposes
			ofstream oss;         // ofstream object
			oss << scientific;  // set formatting for input to oss
			oss.precision(5);
			std::ostringstream stm;
			stm << simNum;
			string FileName = "PQTotalVariationOutput";
			FileName += stm.str();
			FileName += ".txt";
			oss.open(FileName.c_str());
			for (size_t j=0; j < TotalVarDist.size(); j++) {
				oss << TotalVarDist[j] << endl;
			}
			oss << flush;
			oss.close();
			cout << "Total variation distance output to " << FileName << endl;

			if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
			}

        // EMPSums will have been adjusted during the splitting process
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}

//new
//gloria addition
bool AdaptiveHistogram::prioritySplitMCMC(const NodeCompObj& compTest,
														const HistEvalObj& he,
														LOGGING_LEVEL logging,
														size_t minChildPoints, 
														double minVolB, 
														size_t maxLeafNodes, 
														std::vector<real>& Posterior,
														LogMCMCPrior& logPrior)
{
    bool retValue = false;
    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplitMCMC(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes,
                                    Posterior, logPrior);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}


// gat41
bool AdaptiveHistogram::prioritySplitMCMC(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes,
                                std::vector<real>& Posterior,
                                LogMCMCPrior& logPrior)
{   
    bool cancontinue = false;
    bool TooManyLeaves = false;
    
    real deltaL = 0;
    real deltaP = 0;
    int removeBox = 0;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {

        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

		  // get the log-likelihood
		  SPSnodePtrs leaves; // set up empty container for leaf node pointers
		  SPSnodePtrsItr it; // and an iterator over the container
		  getSubPaving()->getLeaves(leaves); // fill the container
        for(it = leaves.begin(); it < leaves.end(); it++) {
				//cout << (*it)->getNodeName() << "\t" << (*it)->getCounter() << endl; 
				deltaL += (*it)->getLogLik(n);
		  }

		  // use the prior distribution object to find the prior 
		  deltaP = logPrior(getRootLeaves()-1);
		  // posterior is proportional to likelihood * prior
		  real deltaPi = deltaL + deltaP;
		  // push back into vector
		  Posterior.push_back(deltaPi);

		  cancontinue = (!pq.empty());
        
        bool bigEnough = cancontinue;
        TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // accumulate the changes in scaled EMP sums that will result
            // from this expansion
            //updateScaledEMPSumCOPERR(chosenLargest->getSplitChangeEMPCOPERR(n));
            //updateScaledEMPSumAIC(chosenLargest->getSplitChangeEMPAIC());

            // split the biggest one and divide up its data
             Expand(chosenLargest);
             
				// get the log-likelihood for the children node 
				real changeL = chosenLargest->getLeftChild()->getLogLik(n);
				changeL += chosenLargest->getRightChild()->getLogLik(n);

				// compute the likelihood * prior to get posterior
				size_t realNumLeaves = getRootLeaves(); 
				
				//cout << "---------------------------------------" << endl;
				//cout << "Number of leaves " << realNumLeaves << endl;

				// log likelihood
				deltaL = deltaL - chosenLargest->getLogLik(n) + changeL; 

				// use the prior distribution object to find the prior 
				deltaP = logPrior(realNumLeaves-1);
				
				//cout << deltaL << "\t" << deltaP << endl;
				
				// posterior is proportional to likelihood * prior
				real deltaPi = deltaL + deltaP;
				
				// push back into vector
				Posterior.push_back(deltaPi);

            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed

            if (((chosenLargest->getLeftChild())->getCounter() > removeBox) &&
					   checkNodeCountForSplit(chosenLargest->getLeftChild(),
                  volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                
                pq.insert(chosenLargest->getLeftChild());
            }

            if ( ((chosenLargest->getRightChild())->getCounter() > removeBox) &&
                 checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
               
                pq.insert(chosenLargest->getRightChild());
            }

            if (logging != NOLOG) {
                // To add current state of histogram to log file
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
				
				
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
		  
			}
          

          
          if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
          }

          // EMPSums will have been adjusted during the splitting process
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}

//new
//gloria addition
bool AdaptiveHistogram::prioritySplitGet(const NodeCompObj& compTest,
														const HistEvalObj& he,
														LOGGING_LEVEL logging,
														size_t minChildPoints, 
														double minVolB, 
														size_t maxLeafNodes, 
													std::vector<AdaptiveHistogram>& States, 
													std::vector<size_t>& Sampled)
{
    bool retValue = false;
    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplitGet(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes,
                                    States, Sampled);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}


// gat41
bool AdaptiveHistogram::prioritySplitGet(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes,
										  std::vector<AdaptiveHistogram>& States, 
                                std::vector<size_t>& Sampled)
{   
    bool cancontinue = false;
    bool TooManyLeaves = false;
    
    int removeBox = 0;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {

        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

		  cancontinue = (!pq.empty());
        
        bool bigEnough = cancontinue;
        TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // accumulate the changes in scaled EMP sums that will result
            // from this expansion
            //updateScaledEMPSumCOPERR(chosenLargest->getSplitChangeEMPCOPERR(n));
            //updateScaledEMPSumAIC(chosenLargest->getSplitChangeEMPAIC());

				// collect the needed states here
				if ( getRootLeaves() == (Sampled[0]+1) ) {
					States.push_back(*this);
				}
				else if ( getRootLeaves() == (Sampled[1]+1) ) {
					States.push_back(*this);
					break;
				}

            // split the biggest one and divide up its data
             Expand(chosenLargest);

            // add the new child names to the creation string
            creationString += chosenLargest->getChildNodeNames();

            // but only put the children into the container if they can be
            // split, which means IF the child meets the min vol test AND IF
            // either there are enough points in the whole child and
                // the child's minChildCountIfSplit is 0 (ie all points go to
                // one child of the child)
            // or the child's minChildCountIfSplit test is passed

            if (((chosenLargest->getLeftChild())->getCounter() > removeBox) &&
					   checkNodeCountForSplit(chosenLargest->getLeftChild(),
                  volChecking, minVol, minChildPoints)) {
                // insert the new left child into the multiset
                
                pq.insert(chosenLargest->getLeftChild());
            }

            if ( ((chosenLargest->getRightChild())->getCounter() > removeBox) &&
                 checkNodeCountForSplit(chosenLargest->getRightChild(),
                    volChecking, minVol, minChildPoints)) {
                // insert the new right child into the multiset
               
                pq.insert(chosenLargest->getRightChild());
            }

            if (logging != NOLOG) {
                // To add current state of histogram to log file
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
				
				
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
		  
			}
          

          
          if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
          }

          // EMPSums will have been adjusted during the splitting process
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}




// check whether we can stop splitting using some stopping criteria
bool AdaptiveHistogram::checkStopCrit(double stopCritCurrent, double stopCritPrevious, int& Prev)
{
	bool valley = false;
	//cout << "current: " << stopCritCurrent << "\t previous: " << stopCritPrevious << endl;
	//cout << "Prev before checks: " << Prev << endl;
	
	// check if it is a local minimum
	if ( (stopCritCurrent > stopCritPrevious) && (Prev == 1) ) {
		//cout << "larger: " << (stopCritCurrent > stopCritPrevious) << "\t Prev: " << Prev << endl;
		//cout << "!!! Local minima previously !!! " << endl;
		Prev = 0;
		valley = true;
	}
	
	// check if stopCrit is decreasing
	else if (stopCritCurrent < stopCritPrevious) {
		Prev = 1;
		valley = false;
	}

	else { 
		valley = false; }

	//cout << "Prev after checks: " << Prev << endl;
	return valley;
}

//src_trunk_0701
void AdaptiveHistogram::swap(AdaptiveHistogram& adh) // throw()
{
	//std::swap(label, adh.label);
	std::swap(dataCollection, adh.dataCollection); // use stl specialisation of swap
    std::swap(holdAllStats, adh.holdAllStats);
	std::swap(creationString, adh.creationString);
    
	// cxsc don't seem to have a swap for dot precisions
    dotprecision tempCOPERR(adh.scaledEMPSumCOPERR);
	adh.scaledEMPSumCOPERR = scaledEMPSumCOPERR;
	scaledEMPSumCOPERR = tempCOPERR;
	
	dotprecision tempAIC(adh.scaledEMPSumAIC);
	adh.scaledEMPSumAIC = scaledEMPSumAIC;
	scaledEMPSumAIC = tempAIC;
	
	std::swap(rootPaving, adh.rootPaving); // just swap the pointers
}

//----------------
bool AdaptiveHistogram::prioritySplitMappedIAE(
														const NodeCompObj& compTest,
                            const HistEvalObj& he,
                            LOGGING_LEVEL logging,
                            size_t minChildPoints, double minVolB, 
														size_t maxLeafNodes, 
														PiecewiseConstantFunction& nodeEst,
														std::vector<real>& vecIAE)
{
    bool retValue = false;

    gsl_rng * rgsl = NULL;

    try {
        // set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed

        retValue = prioritySplitMappedIAE(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl, maxLeafNodes,
                                    nodeEst, vecIAE);
        gsl_rng_free (rgsl);
    }

    catch (bad_alloc& ba) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(ba.what());
        string msg = "Error allocating memory in priority split.  Orginal error: "
                                     + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        if (NULL != rgsl) gsl_rng_free(rgsl); // free the random number generator
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return retValue;
}

bool AdaptiveHistogram::prioritySplitMappedIAE(
																const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl, size_t maxLeafNodes,
                                PiecewiseConstantFunction& nodeEst,
																std::vector<real>& vecIAE)
{    
		size_t numHist = 0; //a counter to track the number of histograms
	    
    bool cancontinue = false;
    bool TooManyLeaves = false;
    
    if (NULL == rootPaving) {
            throw HistException("No root paving for prioritySplit");
    }

    try {
				numHist += 1;
				cout << "---- Hist " << numHist << "-----" << endl;
				// get the IAE of the first histogram
				PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(*this);
				real IAE = nodeEst.getIAE(*tempPCF);
				delete tempPCF;
				(vecIAE).push_back(IAE);

        bool volChecking = false; // record if we need to check volume before split
        double minVol = -1.0; // minimum volume (used only if checking)
        size_t n; // for number of points in histogram

        int i = 0;
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
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

        n = rootPaving->getCounter(); // number of points in histogram

        if (logging != NOLOG) {
             // Start log file with filename and timestamp
            outputLogStart(s);
            // log the current state of the histogram
            outputLog(s, i);
            outputLogEMPAIC(s); // add AIC scores
            i++;
        }

        // put nodes into the starting set IF they meet minVol test AND IF either
        // there are enough points in the whole node
                // and minChildCountIfSplit is 0 (ie all points go to one child)
        // or the minChildCountIfSplit test passed

        if (rootPaving->isLeaf()) {
            // check to insert a copy of the rootPaving pointer into the set
            if (checkNodeCountForSplit(rootPaving, volChecking, minVol,
                minChildPoints)) {
                    pq.insert(rootPaving);
            }
        }
        else { // root is not a leaf
            SPSnodePtrs leaves;
            rootPaving->getLeaves(leaves);
            // check to insert each of the leaves into the set
            SPSnodePtrsItr sit;
            for (sit = leaves.begin(); sit < leaves.end(); sit++) {
                if (checkNodeCountForSplit((*sit), volChecking, minVol,
                minChildPoints)) {
                    pq.insert(*sit);
                }
            }
        }

        cancontinue = (!pq.empty());
		  
        bool bigEnough = cancontinue;
	     TooManyLeaves = (getRootLeaves() > maxLeafNodes);

        if(!cancontinue) {
            std::cout << "No splittable leaves to split - aborting" << std::endl;
        }

			
        // split until the HistEvalObj he () operator returns true
        // we only put splittable nodes into the set, so we don't have to check
        // that they are splittable when we take them out
        while (bigEnough && !he(this) && !TooManyLeaves) {
            
            SPSnode* largest = *(pq.rbegin ()); // the last largest in the set
            SPSnode* chosenLargest;
            
            // find if there are any more equal to largest around
            multiset<SPSnode*, MyCompare>::iterator mit;
            pair<multiset<SPSnode*, MyCompare>::iterator,
                multiset<SPSnode*, MyCompare>::iterator> equalLargest;

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
                multiset<SPSnode*, MyCompare>::iterator it = pq.end();
                it--;
                pq.erase(it);// take this largest out of the set
            }

            // split the biggest one and divide up its data
						// cout << "chosenLargest: " << chosenLargest->getNodeName() << "\t" << chosenLargest->getCounter() << endl;
            Expand(chosenLargest);
            
          	numHist += 1;
          	cout << "---- Hist " << numHist << "-----" << endl;
						// get the IAE of the first histogram
						PiecewiseConstantFunction* tempPCF = new PiecewiseConstantFunction(*this);
						real IAE = nodeEst.getIAE(*tempPCF);
						delete tempPCF;
						(vecIAE).push_back(IAE);
            
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
                outputLog(s, i);
                outputLogEMPAIC(s); // add AIC scores
                i++;
            }

            bigEnough = (!pq.empty());
            if (!bigEnough)
                std::cout << "Terminated splitting: no splittable nodes left"
                    << std::endl;
				
				
				// check if number of leaf nodes in subpaving > maxLeafNodes
				// maximum number of leaf nodes allowed
				//n^B, A+B > 1, 0  < A < 1, 0 < B < 1 - refer Prop. 1 in PQ paper
				TooManyLeaves = (getRootLeaves() > maxLeafNodes);
				if ( TooManyLeaves) {
					std::cout << "Terminated splitting: maximum number of leaf nodes = "<< maxLeafNodes << " reached"
                          << std::endl;
				}
		  
	}

        if (cancontinue && (logging != NOLOG)) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());

        }
   }

    catch (bad_alloc& ba) {
        string oldmsg(ba.what());
        string msg = "Error allocating memory iin priority split.  Orginal error: "
                                    + oldmsg;
        std::cout << msg << std::endl;
        throw HistException(msg);
    }
    catch (HistException& e) {
        string oldmsg(e.what());
        string msg = "HistException error in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (SPnodeException& spe) {
        string oldmsg(spe.what());
        string msg = "SPnodeException in priority split.  Orginal error: "
                                    + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }
    catch (exception& e) {
        string oldmsg(e.what());
        string msg = "Error in priority split.  Orginal error: " + oldmsg;
        std::cerr << msg << std::endl;
        throw HistException(msg);
    }

    return (cancontinue);
}
//---------




// ------------    implementation of HistDescription class -----------------
// ----------  private methods

// ----------  public methods:

// Copy assignment operator.
HistDescription& HistDescription::operator=(const HistDescription& rhs)
{
    depthString = rhs.depthString;
    goodString = rhs.goodString;
    return *this;
}

//  Get the first level in the description.
int HistDescription::peekFirst() const
{
    // parse string to get first level out
    std::string sep = ", ";
    size_t startpos = depthString.find_first_not_of(sep);
    size_t endpos = std::string::npos;
    size_t newstartpos = std::string::npos;
    int depth = 0;
    std::string str = "";
    if (startpos != std::string::npos) {
        endpos = depthString.find_first_of(sep, startpos);
        //not the last digit
        if (endpos != std::string::npos) {
            str = depthString.substr(startpos, endpos-startpos);
        }
        //last digit
        else {
            str = depthString.substr(startpos);
        }
    }
    depth = atoi(str.c_str()); // 0 if not valid integer
    if (depth == 0) {
        goodString = false;
    }

    return depth; // 0 if not valid integer
}

//  Strip off the first level in the description.
bool HistDescription::popFirst()
{
    bool stripped = false;
    // parse string to get first level out
    std::string sep = ", ";
    size_t startpos = depthString.find_first_not_of(sep);
    size_t endpos = std::string::npos;
    size_t newstartpos = std::string::npos;
    if (startpos != std::string::npos) {
        endpos = depthString.find_first_of(sep, startpos);
        //not the last digit
        if (endpos != std::string::npos) {
            newstartpos = depthString.find_first_not_of(sep,
                                                    endpos + 1);
        }

        if (newstartpos != std::string::npos) {
            depthString = depthString.substr(newstartpos);
        }
        else depthString = "";
        stripped = true;
    }

    return stripped;
}

// Output the HistDescription object.
std::ostream & operator<<(std::ostream &os,
                                const HistDescription& hd)
{
    os << hd.getDepthString();
    return os;
}

// Comparison operator for the histogram description.
bool operator<(const HistDescription& lhs,
                        const HistDescription& rhs)
{
    return lhs.getDepthString() < rhs.getDepthString();
}

// -----------  end of implementation of HistDescription class ---------------


// ----------------------------- non member functions

//Output all boxes in AdaptiveHistogram adh
std::ostream & operator<<(std::ostream &os, const AdaptiveHistogram& adh)
{
    os << (adh.getSubPaving())->nodesAllOutput(os, 1) << std::endl;

    return os;
}

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap (AdaptiveHistogram & a1, 
		AdaptiveHistogram & a2) // throw ()
{
	a1.swap(a2);
}


// ----------------------------- Histogram exceptions definitions

HistException::HistException(std::string ss) : s(ss) {}
HistException::~HistException () throw () {}
const char* HistException::what() const throw() { return s.c_str(); }

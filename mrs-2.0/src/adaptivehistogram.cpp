/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
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

/*! \file
\brief AdaptiveHistogram definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "adaptivehistogram.hpp"
#include "MCMCPartitionGenerator.hpp"
#include "collatorspnode.hpp"
#include "adaptivehistogramcollator.hpp"

#include "adaptivehistogram_changeofstateinfo_basic.hpp"

// to use toolz includes (including std::vector) and toolz methods
#include "toolz.hpp"
#include "sptools.hpp"

#include "subpaving_exception.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <numeric>	// to use accumulate etc
#include <list>     // to use stl:: lists
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <cassert> // must be after definitn of NDEBUG if there is one
#include <climits>

#include <gsl/gsl_randist.h>

#include <math.h> // math library
#include <cfloat>

//#define MYDEBUG

//#define OLDLOGSTYLE // to get the old logging style - effectively obselete now
//#define DEBUG_NEWMCMC
//#define DEBUG_IMHMCMC

//#define LOGEMPS // if you want lots of extra computations for emps

#define MCMCTRACK //tracking number of loops completed in new MCMC

//#define MYDEBUG_MCMCPQ
#define MYDEBUG_MCMCPQ_MIN
//#define MYDEBUG_MCMCPQ_EXTRA
//#define MYDEBUG_MCMCPQ_CHECK_QUEUE
//#define MYDEBUG_MCMCPQ_LOOP
//#define MYDEBUG_MCMCPQ_OUTPUT
//#define DISABLE_MCMC_PCF_LOGGING

/* to look at what happens if we continue the carving queue
 *  from the local maxes*/
//#define MYDEBUG_MCMCPQ_CONTINUE

//#define DEBUGLOGPOSTERIOR
//#define DEBUGLOGPOSTERIORAVG


#ifdef NDEBUG
	#undef MYDEBUG
	#undef DEBUG_NEWMCMC
	#undef DEBUG_IMHMCMC
	#undef MCMCTRACK
	#undef MYDEBUG_MCMCPQ
	#undef MYDEBUG_MCMCPQ_MIN
	#undef MYDEBUG_MCMCPQ_EXTRA
	#undef MYDEBUG_MCMCPQ_CHECK_QUEUE
	#undef MYDEBUG_MCMCPQ_LOOP
	#undef MYDEBUG_MCMCPQ_OUTPUT
	#undef MYDEBUG_MCMCPQ_CONTINUE
	#undef DEBUGLOGPOSTERIOR
	#undef DEBUGLOGPOSTERIORAVG
#endif


using namespace subpavings;
using namespace std;

// a class for comparison between spsnodes
MyCompare::MyCompare(const NodeCompObj& nc) : myNC(nc) {}

bool MyCompare::operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
    { return myNC(lhs, rhs); }

// --- inner class for evaluating queues

AdaptiveHistogram::PrioritySplitQueueEvaluator::PrioritySplitQueueEvaluator(
		SPSNodeMeasure& m, real cs, size_t maxL)
		: measurer(m), critStop(cs), maxLeaves(maxL), usingCritStop(true) {}
			
AdaptiveHistogram::PrioritySplitQueueEvaluator::PrioritySplitQueueEvaluator(
		SPSNodeMeasure& m, size_t cs, size_t maxL)
		: measurer(m), critStop(1.0*cs), maxLeaves(maxL), usingCritStop(true) {}

AdaptiveHistogram::PrioritySplitQueueEvaluator::PrioritySplitQueueEvaluator(
		SPSNodeMeasure& m, int cs, size_t maxL)
		: measurer(m), critStop(1.0*cs), maxLeaves(maxL), usingCritStop(true) {}
		
AdaptiveHistogram::PrioritySplitQueueEvaluator::PrioritySplitQueueEvaluator(
		SPSNodeMeasure& m, size_t maxL)
		: measurer(m), critStop(0.0), maxLeaves(maxL), usingCritStop(false) {}

SPSNodeMeasure& AdaptiveHistogram::PrioritySplitQueueEvaluator::getMeasurer() const
{ return measurer; }
			
real AdaptiveHistogram::PrioritySplitQueueEvaluator::getCritStop() const
{ return critStop; }

size_t AdaptiveHistogram::PrioritySplitQueueEvaluator::getMaxLeaves() const
{ return maxLeaves; }

bool AdaptiveHistogram::PrioritySplitQueueEvaluator::getUsingCritStop() const
{ return usingCritStop; }

void AdaptiveHistogram::PrioritySplitQueueEvaluator::setMaxLeaves(size_t ml)
{ maxLeaves = ml; }

void AdaptiveHistogram::PrioritySplitQueueEvaluator::setUsingCritStop(bool b)
{ usingCritStop = b; }

// --- end the inner class for evaluating queues


// -------------------implementation of AdaptiveHistogram class --------------


// ----------- histogram public methods

// default constructor
// holdAllStats defaults to false.
AdaptiveHistogram::AdaptiveHistogram()
        : label(0), rootPaving(NULL), 
			holdAllStats(false),
			scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    // nothing happens to dataCollection when object is constructed
}

// initialised constructor with bool to control whether all stats maintained
// in root paving
AdaptiveHistogram::AdaptiveHistogram(bool as, int lab)
        : label(lab), rootPaving(NULL), 
			holdAllStats(as),
			scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    // nothing happens to dataCollection when object is constructed
}

// initialised constructor, initialised with ivector for box
// and with bool to control whether all stats are maintained in root paving.
// (defaults to false which means that only counts are maintained in rootpaving)
AdaptiveHistogram::AdaptiveHistogram(const ivector& v, bool as)
        : label(0), rootPaving(NULL), 
			holdAllStats(as),
			scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"AdaptiveHistogram::AdaptiveHistogram(const ivector&, bool)");
		}
        
		rootPaving = new SPSnode(v, !as);

    }
    catch (exception const& e) {
		constructor_error_handler();
    }
}

// initialised constructor, initialised with ivector for box and label,
// and with bool to control whether all stats are maintained in root paving.
// (defaults to false which means that only counts are maintained in rootpaving)
AdaptiveHistogram::AdaptiveHistogram(const ivector& v, bool as, int lab)
        : label(lab), rootPaving(NULL), 
			holdAllStats(as),
			scaledEMPSumCOPERR(0.0), scaledEMPSumAIC(0.0)
{
    try {
        // check the box here
        if (!checkBox(v)) {
			throw subpavings::MalconstructedBox_Error(
			"AdaptiveHistogram::AdaptiveHistogram(const ivector&, bool)");
		}
        
		rootPaving = new SPSnode(v, !as);
		
    }
    catch (exception const& e) {
		constructor_error_handler();
    }
}



// copy constructor
/*
 * Copy constructor needs to make sure that the dataItrs held by
 * the subpaving correctly refer to the data collection held by
 * the new histogram.  The best way to do this seems to be to
 * save the data in another storage container, 
 * copy the subpaving, strip it of data (ie only shape is preserved)
 * make sure copy his has empty data collection, and then fire
 * the data back in from the storage container ...
 * */
AdaptiveHistogram::AdaptiveHistogram(const AdaptiveHistogram& other)
        : label(other.label), 
			rootPaving(NULL),
			holdAllStats(other.holdAllStats)
{
    try {
		
		if (other.hasSubPaving()) {
			rootPaving = new SPSnode(*(other.getSubPaving()));
			
			// and then strip the data
			rootPaving->clearAllDataHeld();
		
			// preserve the data collection
			RVecData allData;
			BigDataCollection tmp = other.getDataCollection();
			allData.reserve( tmp.size() );
			// copy from other dataCollection into allData;
			allData.assign(tmp.begin(), tmp.end());
			
			// and put the data into this subpaving
			insertFromRVec(allData, NOLOG);
			// insert from RVec returns false if no insertion, which
			// includes if there is no data in the vector
			
			//recalc the emps
			recalcScaledEMPSumAIC();
			recalcScaledEMPSumCOPERR();
			
		} // else subpaving is NULL
		
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}



//Destructor
AdaptiveHistogram::~AdaptiveHistogram()
{
	try {
		delete rootPaving;
		rootPaving = NULL;

	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in AdaptiveHistogram destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

//copy assignment operator
//rhs passed by value
AdaptiveHistogram&
            AdaptiveHistogram::operator=(AdaptiveHistogram rhs)
{
	rhs.swap(*this);
	return *this;
	
}


// overloading of += operator
// no safety copy made before operation - user could do this before op
AdaptiveHistogram&
            AdaptiveHistogram::operator+=(const AdaptiveHistogram& rhs)
{
	//nothing to add, just return this
	if ( !rhs.hasSubPaving() || rhs.getSubPaving()->isEmpty() ) {
		return *this;
	}
	
	// we now know that rhs subpaving exists and is not empty
	
	/* to avoid expense of multiple copies for large hists	
	 * I no longer use a temp as a work space, and at the end copy it into this
	 * AdaptiveHistogram temp, but that means that failure cannot be reversed */
	
	// if this has no subpaving or an empty one, it can just copy the other one
	if ( !hasSubPaving() || getSubPaving()->isEmpty() ) {
		
		bool thisHoldAllStats =  getHoldAllStats();
		int lab = getLabel();
		
		(*this) = rhs;
		// reset label to be the one we had before
		setLabel(lab);
		// but only finished if the hold all stats are the same
		// otherwise we have to reverse the new holdAllStats.
		if ( thisHoldAllStats != rhs.getHoldAllStats() ) {
			
			setHoldAllStats(thisHoldAllStats);
		}
	}
	else {
		
		getSubPaving()->unionTreeStructure(rhs.getSubPaving());
		
		// put all the data from the two histograms into this one.
		RVecData allData;

		BigDataCollection tmp1;
		tmp1.swap(dataCollection);
		// this clears my data collection and puts the contents
		// into tmp1
		assert(dataCollection.empty());
		
		BigDataCollection tmp2 = rhs.getDataCollection();
		
		allData.reserve( tmp1.size() + tmp2.size() );
		// copy from the dataCollections into allData;
		allData.assign(tmp1.begin(), tmp1.end());
		allData.insert(allData.end(), tmp2.begin(),
				tmp2.end());
				
		// and put the data into the histogram
		insertFromRVec(allData, NOLOG);
	}
		
	return *this;
}

// overloading of + operator
const AdaptiveHistogram
            AdaptiveHistogram::operator+(const AdaptiveHistogram& rhs) const
{
	/* result of += will have the holdAllStats of the lhs operand
	 * so if just one of them holds all stats, we want this
	 * one to be on the left, and if both or neither hold all stats, either
	 * can be on the left*/
		
	AdaptiveHistogram temp;
	
	if (!getHoldAllStats() && rhs.getHoldAllStats() ) {
		
		temp = rhs;
		temp += (*this);
	}
	
	else {
		
		temp = (*this);
		temp += rhs;
	}
	
	/* temp will have label of rhs or this: if these are different
	 * we want to reset temp's label to 0*/
	if ( getLabel() != rhs.getLabel() ) temp.setLabel(0);
	
	return temp;
	
}


// Return a pointer to the SPSnode this manages.
// this is bad bad bad ...
SPSnode* AdaptiveHistogram::getSubPaving() const
{return rootPaving;}

int AdaptiveHistogram::getLabel() const
{return label;}

BigDataCollection AdaptiveHistogram::getDataCollection() const
{return dataCollection;}

// get the value of holdAllStats field.
bool AdaptiveHistogram::getHoldAllStats() const
{
    return holdAllStats;
}

// get whether this has a subpaving.
bool AdaptiveHistogram::hasSubPaving() const
{
    return ( getSubPaving() != NULL );
}

cxsc::ivector AdaptiveHistogram::getRootBox() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
							"AdaptiveHistogram::getRootBox()");
	}
	return getSubPaving()->getBox();
}

cxsc::real AdaptiveHistogram::getRootVolume() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
							"AdaptiveHistogram::getRootVolume()");
	}
	return getSubPaving()->nodeRealVolume();
}

int AdaptiveHistogram::getDimensions() const
{
	int retValue = 0;
	if (hasSubPaving()) {
		retValue = getSubPaving()->getDimension();
	}
	return retValue;
}

// Gets the mean from the root box of the paving this manages.
rvector AdaptiveHistogram::getRootPavingMean() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::getRootPavingMean()");
	}
	return getSubPaving()->getMean();
}

// Gets variance covariance vector from root box of rootpaving.
RealVec AdaptiveHistogram::getRootPavingVarCovar() const
{
    if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::getRootPavingVarCovar()");
		
	}
	return getSubPaving()->getVarCovar();
}

// Gets the mean from the data collection.
rvector AdaptiveHistogram::getDataCollectionMean() const
{
	return calculateMeanFromBigData();
}

// Gets variance covariance vector the data collection.
RealVec AdaptiveHistogram::getDataCollectionVarCovar() const
{
    
	return calculateVarCovarFromBigData();
}


// Gets count in the rootpaving in the root paving.
size_t AdaptiveHistogram::getRootCounter() const
{ 
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::getRootCounter()");
		
	}
	return getSubPaving()->getCounter();
}

// Gets number of leaf nodes in the root paving.
size_t AdaptiveHistogram::getRootLeaves() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::getRootLeaves()");
		
	}
	return getSubPaving()->getNumberLeaves();
}

// Gets number of cherry nodes in the root paving.
size_t AdaptiveHistogram::getRootCherries() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
							"AdaptiveHistogram::getRootCherries()");
		
	}
	return getSubPaving()->getNumberCherries();
}

// Gets accumulated depth over all leaf nodes of paving
unsigned long int AdaptiveHistogram::getRootTotalLeafDepth() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
							"AdaptiveHistogram::getRootTotalLeafDepth()");
		
	}
	return getSubPaving()->getTotalLeafDepth();
}

// Gets the sum of leaf count over volume in root paving.
real AdaptiveHistogram::getRootSumLeafCountOverVol() const
{ 
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::getRootSumLeafCountOverVol()");
		
	}
	return getSubPaving()->getSumLeafCountOverVol();
}


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

    if (hasSubPaving()) {

        
        
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

    if (hasSubPaving()) {

        
        
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

    if (hasSubPaving()) {

        size_t counter = getRootCounter();
        retValue =  minVolB * log(1.0*counter)*log(1.0*counter)/counter;
    }
    return retValue;
}


// Get a string of the leaf node levels.
std::string AdaptiveHistogram::getLeafLevelsString() const
{
    string retValue = "";
    if (hasSubPaving())
        retValue = getSubPaving()->getLeafNodeLevelsString();

    return retValue;
}


//insert a single data point into the AdaptiveHistogram object
bool AdaptiveHistogram::insertOne(rvector newData,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging)
{
    // make sure we have a paving and then try inserting
    if (!hasSubPaving()) {
        throw NullSubpavingPointer_Error(
        "insertOne(rvector, const SplitDecisionObj&, LOGGING_LEVEL)");
    }
    // make sure we have a paving and then try inserting
    if (getSubPaving()->isEmpty()) {
        throw NoBox_Error(
        "insertOne(rvector, const SplitDecisionObj&, LOGGING_LEVEL)");
    }

    // check the dimensions
    if( (Ub(newData)-Lb(newData) + 1) != getDimensions() ) {
        throw IncompatibleDimensions_Error(
        "insertOne(rvector, const SplitDecisionObj&, LOGGING_LEVEL)");
   }

    // for logging output to keep track of splits
    int i = 0;
    std::string baseFileName = "";
    std::string s = "";

    if (logging != NOLOG) {
        baseFileName = "splitOutput";
        s = getUniqueFilename(baseFileName);
		
		#ifdef INSERTDEBUG
			cout << "logging to " << s << endl;
		#endif
        outputLogStart(s);
        // log the current state of the histogram
        outputLogPlain(s, i);
		#ifdef INSERTDEBUG
			cout << "logged state with i = " << i << endl;
		#endif
        i++;
    }

	bool retValue = false;

    BigDataItr it = dataCollection.end();
    it = dataCollection.insert(it, newData);

    SPSnode* insertedInto = NULL;

    // try inserting
    insertedInto = getSubPaving()->insertOneFind(it, ON_PARENT,
                                            boolTest);

    if (insertedInto==NULL) { // failed to insert
	
		dataCollection.erase(it);
		
        std::cerr << "Failed to insert point ";
		prettyPrint(std::cerr, newData);
		std::cerr << std::endl;
    }
    else { // insertion succeeded
	
		retValue = true;

        // if we split on insert we may want to log
        if(!(insertedInto->isLeaf()) && logging != NOLOG) {
			// log the current state of the histogram
			#ifdef INSERTDEBUG
				cout << "insertion successful, about to log current state" << endl;
			#endif
			outputLogPlain(s, i);
			#ifdef INSERTDEBUG
				cout << "logged state with i = " << i << endl;
			#endif
			i++;
        }

        //recalculate the scaled EMP sum values;
        recalcScaledEMPSumCOPERR();
        recalcScaledEMPSumAIC();
    }
    if (logging != NOLOG) {
        // add leaf node levels string to log
        outputFile(s, getLeafLevelsString());
		#ifdef INSERTDEBUG
			cout << "added leaf node levels string to log s= " << s << endl;
		#endif
    }
	
	return retValue;
}



// deprecated
bool AdaptiveHistogram::insertOneDimDataFromTxt(const std::string& s,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	int dim = 1;
		
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines, dim);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// deprecated
bool AdaptiveHistogram::insertOneDimDataFromTxt(const std::string& s,
				const SplitDecisionObj& boolTest,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	int dim = 1;
		
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}


// No splitting, just insertion into root box
// deprecated
bool AdaptiveHistogram::insertRvectorsFromTxt(const std::string& s,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	return insertRvectorsFromTxtOrd(s, headerlines, logging);
}


// Adaptive splitting with each data point inserted. */
// deprecated
bool AdaptiveHistogram::insertRvectorsFromTxt(const std::string& s,
				const SplitDecisionObj& boolTest,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	return insertRvectorsFromTxtOrd(s, boolTest, headerlines, logging);
}

// No splitting, just insertion into root box
// deprecated
bool AdaptiveHistogram::insertRvectorsFromTxt(const std::string& s,
								int dim, 
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	return insertRvectorsFromTxtOrd(s, dim, headerlines, logging);
}

// Adaptive splitting with each data point inserted. */
// deprecated
bool AdaptiveHistogram::insertRvectorsFromTxt(const std::string& s,
				const SplitDecisionObj& boolTest,
				int dim,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	return insertRvectorsFromTxtOrd(s, boolTest, dim, headerlines, logging);
}

// Ordinary level checking
// No splitting, just insertion into root box
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}


// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
				const SplitDecisionObj& boolTest,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// No splitting, just insertion into root box
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
								int dim, 
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines, dim);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}


// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
				const SplitDecisionObj& boolTest,
				int dim,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, headerlines, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

//new
// no splitting, reqDims
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
								const vector < int>& reqDims,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, reqDims, headerlines);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

//new
// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtOrd(const std::string& s,
								const SplitDecisionObj& boolTest,
								const vector < int>& reqDims,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtOrd(theData, s, reqDims, headerlines);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// paranoid
// No splitting, just insertion into root box
bool AdaptiveHistogram::insertRvectorsFromTxtParanoid(const std::string& s,
								int dim, 
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtParanoid(theData, s, headerlines, dim);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtParanoid(const std::string& s,
				const SplitDecisionObj& boolTest,
				int dim,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtParanoid(theData, s, headerlines, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

//new
// no splitting, reqDims
bool AdaptiveHistogram::insertRvectorsFromTxtParanoid(const std::string& s,
								const vector < int>& reqDims,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtParanoid(theData, s, reqDims, headerlines);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}
//new
// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtParanoid(const std::string& s,
								const SplitDecisionObj& boolTest,
								const vector < int>& reqDims,
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtParanoid(theData, s, reqDims, headerlines);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// fast
// No splitting, just insertion into root box
bool AdaptiveHistogram::insertRvectorsFromTxtFast(const std::string& s,
								int dim, 
								const std::size_t headerlines,
                                LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtFast(theData, s, headerlines, dim);

	if (retValue) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// Adaptive splitting with each data point inserted. */
bool AdaptiveHistogram::insertRvectorsFromTxtFast(const std::string& s,
				const SplitDecisionObj& boolTest,
				int dim,
				const std::size_t headerlines,
				LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	bool retValue = readRvectorsFromTxtFast(theData, s, headerlines, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}

	return retValue;
}

// no splitting
bool AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls(
								const std::vector < VecDbl >& inputData,
								LOGGING_LEVEL logging)
{
	SplitNever boolTest; // a dummy split decision object
        
	RVecData theData; // container for the rvectors we take in

	// try to get data from the vector
	bool retValue = getRvectorsFromVectorOfVecDbl(theData, inputData);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}


// adaptive splitting
bool AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls(
					const std::vector < VecDbl >& inputData,
                    const SplitDecisionObj& boolTest,
					LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	// try to get data from the vector
	bool retValue = getRvectorsFromVectorOfVecDbl(theData, inputData);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}

// no splitting
bool AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls(
								const std::vector < VecDbl >& inputData,
								int dim,
								LOGGING_LEVEL logging)
{
	SplitNever boolTest; // a dummy split decision object
        
	RVecData theData; // container for the rvectors we take in

	// try to get data from the vector
	bool retValue = getRvectorsFromVectorOfVecDbl(theData, inputData, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}


// adaptive splitting
bool AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls(
					const std::vector < VecDbl >& inputData,
                    const SplitDecisionObj& boolTest,
					int dim,
					LOGGING_LEVEL logging)
{
	RVecData theData; // container for the rvectors we take in

	// try to get data from the vector
	bool retValue = getRvectorsFromVectorOfVecDbl(theData, inputData, dim);

	if (retValue) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}

// no splitting
bool AdaptiveHistogram::insertFromRVec(const RVecData& inputData, 
							LOGGING_LEVEL logging)
{
	bool retValue = false;

	RVecData theData; // container for the rvectors we take in
	
	size_t numberFound = getRvectorsFromRVec(theData, inputData);

	if (numberFound) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}

// method to insert all rvectors from an RVecData object
bool AdaptiveHistogram::insertFromRVec(const RVecData& inputData,
                            const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging)
{
	bool retValue = false;

	RVecData theData; // container for the rvectors we take in
	
	size_t numberFound = getRvectorsFromRVec(theData, inputData);

	if (numberFound) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}

// no splitting
bool AdaptiveHistogram::insertFromRVec(const RVecData& inputData,
							bool checkDims,
							LOGGING_LEVEL logging)
{
	bool retValue = false;

	RVecData theData; // container for the rvectors we take in
	
	size_t numberFound = getRvectorsFromRVec(theData, inputData, checkDims);

	if (numberFound) {
		SplitNever boolTest; // a dummy split decision object
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
	}
	
	return retValue;
}

// method to insert all rvectors from an RVecData object
bool AdaptiveHistogram::insertFromRVec(const RVecData& inputData,
                            const SplitDecisionObj& boolTest,
							bool checkDims,
                            LOGGING_LEVEL logging)
{
	bool retValue = false;

	RVecData theData; // container for the rvectors we take in
	
	size_t numberFound = getRvectorsFromRVec(theData, inputData, checkDims);

	if (numberFound) {
		retValue = completeDataInsertionFromVec(theData,
												boolTest, logging);
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

	return retValue;
}

// method to insert a sample of rvectors from a container of rvectors
// this version takes seed for a random number generator
bool AdaptiveHistogram::insertSampleFromRVec(size_t samplesize,
            int seed, const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
	gsl_rng * rgsl = NULL;

    try {
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        gsl_rng_set (rgsl, seed); // change the seed

        bool retValue = insertSampleFromRVec(samplesize, rgsl, rvec, boolTest,
            logging);

        gsl_rng_free(rgsl); // free the random number generator
		rgsl = NULL;

		return retValue;

    }
    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}

// method to insert a sample of rvectors from a container of rvectors
// this version will set up a random number generator with default seed
bool AdaptiveHistogram::insertSampleFromRVec(size_t samplesize,
            const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
    gsl_rng * rgsl = NULL;

    try {
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
        bool retValue = insertSampleFromRVec(samplesize, rgsl, rvec, boolTest,
                logging);

        gsl_rng_free(rgsl); // free the random number generator
		rgsl = NULL;
		
		return retValue;
    }
    catch (exception const& e) {
        
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }

}

// method to insert rvectors from an RSSample object
bool AdaptiveHistogram::insertFromRSSample(const RSSample& rss,
										int lab,
                                        const SplitDecisionObj& boolTest,
                                        LOGGING_LEVEL logging)
{
	bool retValue = false;
	
	RVecData myDataRvectors; // container for the rvectors we take in

	// try to get data from rss.Samples and check how many data points found
	// uses the label from this
	size_t numberFound = getRvectorsFromRSSample(myDataRvectors, rss, lab);

	if (numberFound > 0) {
		
		// complete the data insertion
		retValue = completeDataInsertionFromVec(myDataRvectors,
												boolTest, logging);

	}

	return retValue;
}

// method to insert a sample of rvectors from an RSSample object
// this version takes a random number generator
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
			gsl_rng * rgsl, const RSSample& rss, int lab,
            const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
	bool retValue = false;

	RVecData myDataRvectors; // container for the rvectors we take in

	// try to sample data from rss.Samples and check how many data points found
	size_t numberTaken = getSampleRvectorsFromRSSample(myDataRvectors,
							rgsl, samplesize, rss, lab);

	if (numberTaken > 0) {
		/* switch on for more output during histogram creation
		// confirm the amount of data taken from the RSSample
		std::cout << "End of taking sample from data from RSSample: "
			<< numberTaken << " data points used for sample" << std::endl;
		*/

		retValue = completeDataInsertionFromVec(myDataRvectors,
												boolTest, logging);
	}
	return retValue;
}

// method to insert a sample of rvectors from an RSSample object
// this version takes seed for a random number generator
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
			int seed, const RSSample& rss, int lab,
            const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
    gsl_rng * rgsl = NULL;

    try {
		
        // make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        gsl_rng_set (rgsl, seed); // change the seed

        bool retValue = insertSampleFromRSSample(samplesize, 
		rgsl, rss, lab, boolTest, logging);

        gsl_rng_free(rgsl); // free the random number generator
		rgsl = NULL;
		
		return retValue;
    }
    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}

// method to insert a sample of rvectors from an RSSample object
// this version will set up a random number generator with default seed
bool AdaptiveHistogram::insertSampleFromRSSample(size_t samplesize,
			const RSSample& rss, int lab,
            const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
{
	gsl_rng * rgsl = NULL;

    try {

       // make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
        bool retValue = insertSampleFromRSSample(samplesize, 
				rgsl, rss, lab, boolTest, logging);

        gsl_rng_free(rgsl); // free the random number generator
		rgsl = NULL;
		
		return retValue;
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}



/* method to make a leaf node histogram into a multi-node histogram
 * by prioritising which node to split first
 * keeps splitting until the function object he returns true
 * or until there are no more splittable nodes
 * outputs to a log file if logging is true
 * makes its own random number generator */
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB)
{
    gsl_rng * rgsl = NULL;

    try {
       // make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
        bool retValue = prioritySplit(compTest, he, logging,
                                    minChildPoints, minVolB, rgsl);
        gsl_rng_free (rgsl);
		rgsl = NULL;
		
		return retValue;
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}


/* method to make a leaf node histogram into a multi-node histogram
 * by prioritising which node to split first
 * keeps splitting until the function object he returns true
 * or until there are no more splittable nodes
 * outputs to a log file if logging required */
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                const HistEvalObj& he,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl)
{
	string errorMsg("AdaptiveHistogram::prioritySplit(const NodeCompObj&, const HistEvalObj&, LOGGING_LEVEL, size_t, double, gsl_rng *)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	bool volChecking = false; // record if we need to check volume before split
	double minVol = -1.0; // minimum volume (used only if checking)
	size_t n = 0; // for number of points in histogram

	// make volChecking true if minVolB is > 0.0
	if (minVolB > 0.0) {
		// minimum volume of a splittable node is minVolB(log n)^2/n
		minVol = getMinVol(minVolB);
		volChecking = true;
		
	}
	
	#ifdef MYDEBUG
		cout << "logging = " << logging << endl;
		cout << "minChildPoints = " << minChildPoints << endl;
		cout << "volChecking = " << volChecking << endl;
		cout << "minVol = " << minVol << endl;
	#endif

	// a multiset for the queue (key values are not necessarily unique)
	multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

	_setupPrioritySplit(pq, minChildPoints, minVol);
	// throws an exception if current state not legal

	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	n = getRootCounter(); // number of points in histogram
	
	#ifdef MYDEBUG
		cout << "rootCounter = " << n << endl;
		
	#endif

	int i = 0;
	
	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		i++;
	}

	size_t numLeaves = getRootLeaves(); // actual number of leaves

	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	
	/* split until the HistEvalObj he () operator returns true
	 we only put splittable nodes into the set, so we don't have to check
	 that they are splittable when we take them out */
	while (canContinue && !he(this)) {
		
		canContinue = _prioritySplitLoop(pq, n, minChildPoints, minVol, rgsl);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			i++;
		}
		
		/*any successful split will have increased actual number of 
		leaves */
		if (canContinue) numLeaves++;

	}

	if (canContinue && (logging != NOLOG)) {
		// log the leaf levels line
		outputFile(s, getLeafLevelsString());

	}

	// EMPSums will have been adjusted during the splitting process
	
	return canContinue; // true if fe satisfied
}


/* method to make a leaf node histogram into a multi-node histogram
 * by prioritising which node to split first
 * keeps splitting until the function object he returns true
 * or until there are no more splittable nodes
 * outputs to a log file if logging is true
 * makes its own random number generator */
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                size_t maxLeaves,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB)
{
    gsl_rng * rgsl = NULL;

    try {
       // make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
        bool retValue = prioritySplit(compTest, maxLeaves, logging,
                                    minChildPoints, minVolB, rgsl);
        gsl_rng_free (rgsl);
		rgsl = NULL;
		
		return retValue;
    }

    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}


// method to make a leaf node histogram into a multi-node histogram
// by prioritising which node to split first
// keeps splitting until we have maxLeaves returns true
// or until there are no more splittable nodes
// outputs to a log file if logging required
bool AdaptiveHistogram::prioritySplit(const NodeCompObj& compTest,
                                size_t maxLeaves,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, double minVolB,
                                gsl_rng * rgsl)
{
	string errorMsg("AdaptiveHistogram::prioritySplit(const NodeCompObj&, size_t, LOGGING_LEVEL, size_t, double, gsl_rng *)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	bool volChecking = false; // record if we need to check volume before split
	double minVol = -1.0; // minimum volume (used only if checking)
	size_t n = 0; // for number of points in histogram

	// make volChecking true if minVolB is > 0.0
	if (minVolB > 0.0) {
		// minimum volume of a splittable node is minVolB(log n)^2/n
		minVol = getMinVol(minVolB);
		volChecking = true;
		
	}
	
	#ifdef MYDEBUG
		cout << "logging = " << logging << endl;
		cout << "minChildPoints = " << minChildPoints << endl;
		cout << "volChecking = " << volChecking << endl;
		cout << "minVol = " << minVol << endl;
	#endif

	// a multiset for the queue (key values are not necessarily unique)
	multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));

	_setupPrioritySplit(pq, minChildPoints, minVol);
	// throws an exception if current state not legal

	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	n = getRootCounter(); // number of points in histogram
	
	#ifdef MYDEBUG
		cout << "rootCounter = " << n << endl;
		
	#endif

	int i = 0;
	
	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		i++;
	}

	size_t numLeaves = getRootLeaves(); // actual number of leaves

	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	
	/* split until we have maxLeaves leaves
	 we only put splittable nodes into the set, so we don't have to check
	 that they are splittable when we take them out */
	while (canContinue && numLeaves < maxLeaves) {
		
		canContinue = _prioritySplitLoop(pq, n, minChildPoints, minVol, rgsl);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			i++;
		}
		
		/*any successful split will have increased actual number of 
		leaves */
		if (canContinue) numLeaves++;

	}

	if (canContinue && (logging != NOLOG)) {
		// log the leaf levels line
		outputFile(s, getLeafLevelsString());

	}

	// EMPSums will have been adjusted during the splitting process
	
	return canContinue; // true if fe satisfied
}

	
/* method to make a multi-node histogram into one with possibly fewer nodes
 * by prioritising which subleaf node to merge first
 * keeps merging until the stopTest is satisfied or runs out of subleaves
 * outputs to a log file if logging is true */
bool AdaptiveHistogram::priorityMerge(const NodeCompObj& compTest,
                        const HistEvalObj& he, LOGGING_LEVEL logging)
{

    gsl_rng * rgsl = NULL;

    try {
		
		if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"AdaptiveHistogram::priorityMerge(const NodeCompObj&, const HistEvalObj&, LOGGING_LEVEL)");
		}

        if (getSubPaving()->isLeaf()) {
			throw UnfulfillableRequest_Error(
			"AdaptiveHistogram::priorityMerge(const NodeCompObj&, const HistEvalObj&, LOGGING_LEVEL)");
        }
 		
	    size_t n = getRootCounter(); // number of points in histogram

        int i = 0;
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqMergeOutput";
            s = getUniqueFilename(baseFileName);
        }

		if (logging != NOLOG) {
			 // Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the histogram
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			i++;
		}
		
		// a multiset for the queue (key values are not necessarily unique)
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));


		_setupPriorityMerge(pq);
		
		size_t numLeaves = getRootLeaves(); // actual number of leaves
		
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
		bool canContinue = !pq.empty();
		
		if (canContinue) {
			
			bool bigEnough = true;
		
			// merge until the HistEvalObj he () operator returns true
			while (bigEnough && !he(this)) {
				
				bigEnough = AdaptiveHistogram::_priorityMergeLoop(
															pq, n,rgsl);

				if (logging != NOLOG) {
					// To add current state of histogram to log file
					outputLogPlain(s, i);
					#ifdef LOGEMPS
						outputLogEMPAIC(s); // add AIC scores
					#endif
					i++;
				}
				
				if (bigEnough) numLeaves--;
			}
		}
		else {
			std::cerr << "No mergable cherries - aborting" << std::endl;
		}

        if (canContinue && logging != NOLOG) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
        }

        // EMPSums will have been adjusted during the merging process
        if (NULL != rgsl) {
			gsl_rng_free(rgsl);
			rgsl = NULL;
		}
		
		return canContinue;  // true unless we could not even start merge
    }
    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}

/* method to make a multi-node histogram into one with possibly fewer nodes
 * by prioritising which subleaf node to merge first
 * keeps merging until the stopTest is satisfied or runs out of subleaves
 * outputs to a log file if logging is true */
bool AdaptiveHistogram::priorityMerge(const NodeCompObj& compTest,
                        size_t minLeaves, LOGGING_LEVEL logging)
{

    gsl_rng * rgsl = NULL;

    try {
		
		if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
		"AdaptiveHistogram::priorityMerge(const NodeCompObj&, const HistEvalObj&, LOGGING_LEVEL)");
		}

        if (getSubPaving()->isLeaf()) {
			throw UnfulfillableRequest_Error(
			"AdaptiveHistogram::priorityMerge(const NodeCompObj&, size_t, LOGGING_LEVEL)");
        }
 		
	    size_t n = getRootCounter(); // number of points in histogram

        int i = 0;
        std::string baseFileName = "";
        std::string s = "";
        if (logging != NOLOG) {
            // pass to log output to keep track of splits
            baseFileName = "pqMergeOutput";
            s = getUniqueFilename(baseFileName);
        }

		if (logging != NOLOG) {
			 // Start log file with filename and timestamp
			outputLogStart(s);
			// log the current state of the histogram
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			i++;
		}
		
		// a multiset for the queue (key values are not necessarily unique)
        multiset<SPSnode*, MyCompare> pq((MyCompare(compTest)));


		_setupPriorityMerge(pq);
		
		size_t numLeaves = getRootLeaves(); // actual number of leaves
		
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        long unsigned int seed = 1234;
		gsl_rng_set (rgsl, seed); // change the seed
		
		bool canContinue = !pq.empty();
		
		if (canContinue) {
			
			bool bigEnough = true;
		
			// merge until the HistEvalObj he () operator returns true
			while (bigEnough && numLeaves > minLeaves) {
				
				bigEnough = AdaptiveHistogram::_priorityMergeLoop(
															pq, n,rgsl);

				if (logging != NOLOG) {
					// To add current state of histogram to log file
					outputLogPlain(s, i);
					#ifdef LOGEMPS
						outputLogEMPAIC(s); // add AIC scores
					#endif
					i++;
				}
				
				if (bigEnough) numLeaves--;
			}
		}
		else {
			std::cerr << "No mergable cherries - aborting" << std::endl;
		}

        if (canContinue && logging != NOLOG) {
            // log the leaf levels line
            outputFile(s, getLeafLevelsString());
        }

        // EMPSums will have been adjusted during the merging process
        if (NULL != rgsl) {
			gsl_rng_free(rgsl);
			rgsl = NULL;
		}
		
		return canContinue;  // true unless we could not even start merge
    }
    catch (exception const& e) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}

/* method to make a multi-node histogram into a single node histogram
 by merging up to the root box */
bool AdaptiveHistogram::mergeUp()
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error("AdaptiveHistogram::mergeUp()");
	}
			
	bool merged = false;

	if (!getSubPaving()->isLeaf()) {
		
		merged = true;

		SPSnodePtrs subleaves;

		getSubPaving()->getSubLeaves(subleaves); // subleaves contains the subleaves
		
		if(subleaves.empty()) {
            throw std::logic_error(
			"AdaptiveHistogram::mergeUp()) : Error in queue");
		}

		bool canmerge = true;

		// merge until there is only one leaf
		while (canmerge) {

			SPSnode* target = *(subleaves.rbegin ()); // the last in the vector

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
			std::cout << "Merged to root" << std::endl;
		}

		// EMPSums are not adjusted during the merging process
		recalcScaledEMPSumAIC();
		recalcScaledEMPSumCOPERR();
	}
	else {
		std::cerr << "Nothing to be done - root paving is already a leaf "
				<< std::endl;
	}
	
	return merged;
}

/* get a summary of the non-empty boxes in the histogram  in the form of a pair, 
 * (number of empty boxes, accumulated vol of empty boxes) */
std::pair<size_t, cxsc::real> AdaptiveHistogram::getNonEmptyBoxSummary() const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::getNonEmptyBoxSummary()");
	}
	
	
	return getSubPaving()->getNonEmptyBoxSummary();
}

// do checks and use private coverage method to find coverage
double AdaptiveHistogram::findCoverage(const rvector& pt) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::findCoverage(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
				"AdaptiveHistogram::findCoverage(const rvector&)");
	}
	
	return _coverage(pt);
}

// do checks and use private density method to find density
double AdaptiveHistogram::findEmpiricalDensity(const rvector& pt) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::findEmpiricalDensity(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
				"AdaptiveHistogram::findEmpiricalDensity(const rvector&)");
	}
	
	return _empiricalDensity(pt);
}

// Find if box contains point pt
bool AdaptiveHistogram::histBoxContains(const rvector& pt) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::histBoxContains(const rvector&)");
	}
	if (getDimensions() != (Ub(pt) - Lb(pt) + 1)) {
		throw IncompatibleDimensions_Error(
				"AdaptiveHistogram::histBoxContains(const rvector&)");
	}
	
	return ( getSubPaving()->findContainingNode(pt) != NULL ); 
}


//splits histogram according to string instruction
//returns true if some splitting was achieved
bool AdaptiveHistogram::splitToShape(std::string instruction)
{
	
	// checks:  is there a root paving, is the string properly formed?
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::splitToShape()");
	}
	bool success = false;
	SPSnode temp(*getSubPaving()); // copy to temp
	try {
		if (instruction.length() == 0) {
			throw std::invalid_argument(
				"AdaptiveHistogram::splitToShape() : No instruction");
		}

		std::string legal(", 0123456789");
		if (instruction.find_first_not_of(legal) != std::string::npos) {
			throw std::invalid_argument(
				"AdaptiveHistogram::splitToShape() : Illegal character");
		}

		// all seems to be okay, we can start splitting the root paving
		
		success = getSubPaving()->splitRootToShape(instruction);

		if (!success) {
			handleSplitToShapeError(temp);
	   }
	   
	}
	catch (std::invalid_argument const& ia) {
		cerr << ia.what() << endl;
		handleSplitToShapeError(temp);
		success = false;
	}
	catch (std::logic_error const& le) {
		cerr << le.what() << endl;
		handleSplitToShapeError(temp);
		success = false;
	}
	return success;
	// any other exceptions are unhandled
}


/* Use MCMC routing to make an average histogram as a PCF */
PiecewiseConstantFunction AdaptiveHistogram::MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging)
{
	// gsl_rng_default_seed to default environmental vars
    gsl_rng_env_setup();
    long unsigned int seed = gsl_rng_default_seed;
	
	return MCMC(loops, 
					burnin,
					thinout,
					proposal, logPrior,
					minPoints, logging,
					seed);
}

PiecewiseConstantFunction AdaptiveHistogram::MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
		
		
		
		if (thinout == 0)
			throw std::invalid_argument(
				"AdaptiveHistogram::MCMC(...) : thinout == 0");

		// set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed
		gsl_rng_set (rgsl, seed); // change seed
	
		std::vector < PiecewiseConstantFunction > samples;
		bool average = true;
		
		bool volChecking = false;
		real minVol(0.0);
		
		// use the internal method -> non normalised samples
		_MCMCsamples(samples,
						average,
						loops, 
						burnin,
						thinout,
						proposal, logPrior,
						minPoints, 
						volChecking,
						minVol,
						logging,
						rgsl);
		
		gsl_rng_free (rgsl);
		rgsl = NULL;
        
		if (samples.empty())
			throw std::runtime_error(
				"AdaptiveHistogram::MCMC(...) : no samples collected");
		
		samples.back().normalise();
		return samples.back();
		
    }
    
    catch (exception const& e) {
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
		
}

PiecewiseConstantFunction AdaptiveHistogram::MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					gsl_rng * rgsl)
{
	std::vector < PiecewiseConstantFunction > samples;
	bool average = true;
	
	bool volChecking = false;
	real minVol(0.0);
		
	// use the internal method -> non normalised samples
	_MCMCsamples(samples,
					average,
					loops, 
					burnin,
					thinout,
					proposal, logPrior,
					minPoints, 
					volChecking,
					minVol,
					logging,
					rgsl);
					
	samples.back().normalise();
	return samples.back();
		
}


PiecewiseConstantFunction AdaptiveHistogram::MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging)
{
	// gsl_rng_default_seed to default environmental vars
    gsl_rng_env_setup();
    long unsigned int seed = gsl_rng_default_seed;
	
	return MCMC_IMH(maxLeaves,
					loops, 
					burnin,
					thinout,
					logPrior,
					minPoints, logging,
					seed);
}



PiecewiseConstantFunction AdaptiveHistogram::MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
		
		
		
		if (thinout == 0)
			throw std::invalid_argument(
				"AdaptiveHistogram::MCMC_IMH(...) : thinout == 0");

		// set up a random number generator for uniform rvs
        const gsl_rng_type * tgsl;
        // set the library variables *gsl_rng_default and
        // gsl_rng_default_seed to default environmental vars
        gsl_rng_env_setup();
        tgsl = gsl_rng_default; // make tgsl the default type
        rgsl = gsl_rng_alloc (tgsl); // set up with default seed
		gsl_rng_set (rgsl, seed); // change seed
	
		std::vector < PiecewiseConstantFunction > samples;
		bool average = true;
		
		// use the internal method -> non normalised samples
		_MCMCsamplesIMH(samples,
						average,
						maxLeaves,
						loops, 
						burnin,
						thinout,
						logPrior,
						minPoints, logging,
						rgsl);
		
		gsl_rng_free (rgsl);
		rgsl = NULL;
        
		if (samples.empty())
			throw std::runtime_error(
				"AdaptiveHistogram::MCMC_IMH(...) : no samples collected");
		
		samples.back().normalise();
		return samples.back();
		
    }
    
    catch (exception const& e) {
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
		
}

PiecewiseConstantFunction AdaptiveHistogram::MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					gsl_rng * rgsl)
{
	std::vector < PiecewiseConstantFunction > samples;
	bool average = true;
	
	// use the internal method -> non normalised samples
	_MCMCsamplesIMH(samples,
					average,
					maxLeaves,
					loops, 
					burnin,
					thinout,
					logPrior,
					minPoints, logging,
					rgsl);
					
	samples.back().normalise();
	return samples.back();
		
}

/* The MCMC methods as implemented here run a shadow RMSP alongside
 * the SPSnode tree managed by this to try to make the process more
 * efficient.  The shadow is a representation of a virtual state 
 * of this that is what will be used in samples and averages.  The 
 * actual state of this is effectively the most split of any state
 * reached by the shadow.  That is, if a merge proposal is accepted
 * the RMSP is updated for the merge but the actual state of the SPSnode
 * tree managed by this is not updated (although the container
 * holding the leaves and cherries of this on which split and
 * merge operations can take place is updated, so that this container
 * reflects the virtual state of this - the leaves only being the 
 * splittable virtual leaves - and the
 * proposal mechanism just has to consider what is in this containers).  
 * This saves having to repeatedly
 * move data pointer up and down the tree if merges and splits are
 * repeated.   */

/* Note that when minPoints > 0, proposals are effectively drawn from set of
leaf and cherry nodes which does not include any leaf which, if split, would
have a child whose number of points is < minPoints.  Thus the implementation
needs to distinguish between the overall state of the tree and the
set of splittable leaf nodes.  */
std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging)
{
	// gsl_rng_default_seed to default environmental vars
    gsl_rng_env_setup();
    long unsigned int seed = gsl_rng_default_seed;
	
	MCMCsamples(		samples, 
						loops, 
						burnin,
						thinout,
						proposal, logPrior,
						minPoints, logging,
						seed);
	
	return samples;
}

std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
		
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        gsl_rng_set (rgsl, seed); // change seed
		
		
		bool average = false;
		bool volChecking = false;
		real minVol(0.0);
				
		// use the internal method
		_MCMCsamples(samples,
						average,
						loops, 
						burnin,
						thinout,
						proposal, logPrior,
						minPoints, 
						volChecking,
						minVol,
						logging,
						rgsl);
		
		for (std::vector < PiecewiseConstantFunction >::iterator it = samples.begin();
				it < samples.end();
				++it) {
			it->normalise();
		}
		
		// free the random number generator
        
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (...) {} // catch and swallow
		
		return samples;
		
    }
    
    catch (...) { // catch anything
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (...) {} // catch and swallow
		throw; // rethrow original exception
    }
}

std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	bool average = false;
	bool volChecking = false;
	real minVol(0.0);
		
	_MCMCsamples(		samples, 
						average,
						loops, 
						burnin,
						thinout,
						proposal, logPrior,
						minPoints, 
						volChecking,
						minVol,
						logging,
						rgsl);
	
	for (std::vector < PiecewiseConstantFunction >::iterator it = samples.begin();
			it < samples.end();
			++it) {
		it->normalise();
		}
			
	return samples;
}

std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, 
						LOGGING_LEVEL logging)
{
	// gsl_rng_default_seed to default environmental vars
    gsl_rng_env_setup();
    long unsigned int seed = gsl_rng_default_seed;
	
	MCMCsamplesIMH(		maxLeaves,
						samples, 
						loops, 
						burnin,
						thinout,
						logPrior,
						minPoints,
						logging,
						seed);
	
	return samples;
}

std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints,
						LOGGING_LEVEL logging,
						long unsigned int seed)
{
	gsl_rng * rgsl = NULL;

    try {
		
		// make use mt19937 for generator
        rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
        gsl_rng_set (rgsl, seed); // change seed
		
		
		bool average = false;
		
		// use the internal method
		_MCMCsamplesIMH(
						samples,
						average,
						maxLeaves,
						loops, 
						burnin,
						thinout,
						logPrior,
						minPoints,
						logging,
						rgsl);
		
		for (std::vector < PiecewiseConstantFunction >::iterator it = samples.begin();
				it < samples.end();
				++it) {
			it->normalise();
		}
		
		// free the random number generator
        
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (...) {} // catch and swallow
		
		return samples;
		
    }
    
    catch (...) { // catch anything
		try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
		}
		catch (...) {} // catch and swallow
		throw; // rethrow original exception
    }
}

std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	bool average = false;
	
	_MCMCsamplesIMH(	samples, 
						average,
						maxLeaves,
						loops, 
						burnin,
						thinout,
						logPrior,
						minPoints,
						logging,
						rgsl);
	
	for (std::vector < PiecewiseConstantFunction >::iterator it = samples.begin();
			it < samples.end();
			++it) {
		it->normalise();
		}
			
	return samples;
}



/* Returns success or not for change in state.
 * Updates numLeaves, unscaledLoglik, logPosterior*/
bool AdaptiveHistogram::changeMCMCStateIMH(
						size_t maxLeaves,
                        const MCMCPartitionGenerator& partitioner,
						const LogMCMCIMHPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode& rmsp,
						size_t& numLeaves,
						real& unscaledLoglik,
						real& logPosterior,
						gsl_rng* rgsl, 
						LOGGING_LEVEL logging,
                        const std::string& sDec,
						int i,
						const string& failureLogFilename)
{
	#ifdef DEBUG_IMHMCMC
		cout << "\n------ In changeMCMCStateIMH --------- i = " << i << endl;
		cout << "numLeaves " << numLeaves 
				<< ", unscaledLoglik " << unscaledLoglik << endl;
		cout << "logPosterior = " << logPosterior << endl;
	#endif
			
	/* the proposal has two parts, the number of leaves and the
	 * actual state with that number of leaves.  
	 * 
	 * The partitioner takes the number of leaves required and 
	 * generates a random state with that number of leaves.
	 * 
	 * Doing this on the actual subpaving will return true or false
	 * depending on whether
	 * the proposal can actually be carried out (because of min points
	 * or other splittability issues).
	 * 
	 * I will treat a proposal that cannot be made properly (false)
	 * (because of min points or other splittability issues) as
	 * a proposal that is not accepted, ie the change of state is 
	 * is status quo */

	/* get a number of leaves to change to , betweeen [1, maxLeaves]*/
	unsigned long int newNumLeaves 
		= 1 + gsl_rng_uniform_int (rgsl, maxLeaves);
	
	#ifdef DEBUG_IMHMCMC
		cout << "number of leaves proposed is " << newNumLeaves << endl;
	#endif	
	
	/* ask partitioner to get a state */
	dotprecision nlogn(0.0);
	int ndepth = 0.0;
	bool saveInstructions = true;

	bool success = getSubPaving()->randomKnuthMCMCSplitRootAtLeastSPS(
			newNumLeaves,
			partitioner,
			&rmsp,
			nlogn, ndepth,
			minPoints,
			saveInstructions,
			failureLogFilename);
	
	#ifdef DEBUG_IMHMCMC
		cout << "initial success = " << success << endl;
	#endif

	/* If success we have to choose whether to accept */
	if (success) {
		
		/* posterior = likelihood * prior
		 * log posterior = loglikelihood + log prior
		 * posterior ratio = posterior after / posterior before
		 * log posterior ratio = loglikelihood_after + log prior_after
		 * - (loglikelihood_before + log prior_before)
		 * log posterior ratio = unscaled loglikelihood_after + log prior_after
		 * - (unscaled loglikelihood_before + log prior_before)
		 * log posterior ratio = unscaled loglikelihood_after 
		 * 						- unscaled loglikelihood_before
		 * 						+ log prior_after- log prior_before*/
		
		dotprecision ndlog2(0.0);
		accumulate (ndlog2, cxsc::Ln2_real, 1.0*ndepth);
		cxsc::real unscaledLoglikAfter = (cxsc::rnd(nlogn + ndlog2));
		
		/* change in prior adjusted for change in proposal ratio
		 * ie log prior_after- log prior_before + (ln (C_kPrime) - ln(C_k)) */
		
		std::pair < cxsc::real, cxsc::real> changeLogPriorPair = 
						logPrior.logChange(numLeaves-1, newNumLeaves-1);
		
		cxsc::real changeLogPrior = changeLogPriorPair.first;
		cxsc::real adjChangeLogPrior = changeLogPriorPair.second;
		
		cxsc::real logAcceptanceRatio = unscaledLoglikAfter - unscaledLoglik
										+ adjChangeLogPrior;
		
		#ifdef DEBUG_IMHMCMC
			cout << "unscaledLoglikAfter = " << unscaledLoglikAfter << endl;
			cout << "changeLogLik = " << (unscaledLoglikAfter - unscaledLoglik) << endl;
			cout << "changeLogPrior = " << changeLogPrior  << endl;
			cout << "adjChangeLogPrior = " << adjChangeLogPrior  << endl;
			cout << "logAcceptanceRatio = " << logAcceptanceRatio << endl;
		#endif
		
		//get another random number
		double randChange = gsl_rng_uniform(rgsl);
		if (randChange > 0) success = (cxsc::ln(randChange) < logAcceptanceRatio);

		if ((logging == TXT) || (logging == TXTANDGRAPH)) { // log these values
			
			logMCMCIMH(sDec, i,
						numLeaves,
						newNumLeaves,
						unscaledLoglikAfter - unscaledLoglik,
						changeLogPrior,
						adjChangeLogPrior - changeLogPrior,
						logAcceptanceRatio,
						randChange,
						success);
			
		}
		#ifdef DEBUG_IMHMCMC
			cout << "randChange = " << randChange
					<< " cxsc::ln(randChange) = "  << cxsc::ln(randChange) << endl;
			cout << "success = " << success << endl;
								
		#endif
		
		if (success) {
			
			numLeaves = newNumLeaves;
			logPosterior += (unscaledLoglikAfter - unscaledLoglik
										+ changeLogPrior);
			unscaledLoglik = unscaledLoglikAfter;
		}
		#ifdef DEBUG_IMHMCMC
			cout << "at end of change state, logPosterior = " 
						<< logPosterior << endl;
											
		#endif

	}

	
	else if (((logging == TXT) || (logging == TXTANDGRAPH))) {
			
			logMCMCIMH(sDec, i,
				numLeaves,
				newNumLeaves,
				-cxsc::Infinity,
				0.0,
				0.0,
				0.0,
				0.0,
				success);
			
	}

	return success;
	
}

/* public just so that auto-mcmc methods can use it */
/* Change the virtual state of this */
AdaptiveHistogram::ChangeOfStateInformation&
			AdaptiveHistogram::changeMCMCState(
						SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						AdaptiveHistogram::ChangeOfStateInformation& info, //updated
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i)
{
	real minVol(0.0);
		
	return changeMCMCState (nodes,
                        numLeaves, numCherries,
                        proposal, logPrior,
                        minPoints,
						minVol,
						leaves,
						cherries,
						info,
                        rgsl, logging,
                        sHist, i);
	
}

AdaptiveHistogram::ChangeOfStateInformation&
			AdaptiveHistogram::changeMCMCState(
				SPSnodeList& nodes,
				size_t& numLeaves, size_t& numCherries,
				MCMCProposal& proposal, LogMCMCPrior& logPrior,
				size_t minPoints,
				real minVol,
				RealMappedSPnode::ListPtrs& leaves,
				RealMappedSPnode::ListPtrs& cherries,
				AdaptiveHistogram::ChangeOfStateInformation& info, //updated
				gsl_rng* rgsl, LOGGING_LEVEL logging,
				const std::string& sHist, int i)
{
	#ifdef DEBUG_NEWMCMC
		cout << "\n------ In changeMCMCState --------- i = " << i << endl;
		cout << "numLeaves " << numLeaves << endl;
		cout << "shadowLeaves.size() " << leaves.size() << endl;
		cout << "numCherries " << numCherries << endl;
		cout << "shadowCherries.size() " << cherries.size() << endl;
	#endif
	
	bool volChecking = (minVol > 0.0);
			
	bool haveNode = false;

	// use proposal to fill the proposal probabilities

	// this changes haveNode as well
	SPSnodeListItr it = proposeChangeMCMCState (proposal, nodes,
							numLeaves, numCherries,
							rgsl, haveNode);

	// for change in log-posterior resulting from change in state IF ANY
	cxsc::real deltaPi(0.0);
	
	// only do more if haveNode is true, which means that it points to something
	if (haveNode) {
		
		cxsc::real deltaPiTentative(0.0);
		
		/* should be able to tell from distance it is into nodes
		 * whether it is a virtual leaf or a virtual cherry 
		 * (we can't use isLeaf or hasLeafSibling because we want virtual
		 * state not real state */
		SPSnodeListItr beg_it = nodes.begin();
		if (std::distance(beg_it, it) < static_cast<int>(numLeaves)) { // it pts to a leaf
		
			//grab the leaf
			SPSnode* target = *it; // don't change where it points until erased

			int addLeftChild = 0;
			int addRightChild = 0;
			int deductParent = 0;
			size_t leftChildCount = 0;
			size_t rightChildCount = 0;

			// leaf so we are splitting
			/* addLeftChild, addRightChild, deductParent, 
			 * leftChildCount, rightChildCount and deltaPiTentative
			 * all passed by reference,
			 * real number of leaves is leaves.size()*/			
			bool willSplit = decisionMCMCSplitNEW(target,
                        proposal, logPrior, rgsl,
                        numLeaves, numCherries,
						nodes, leaves.size(), 
						minPoints, 
						volChecking,
						minVol,
						addLeftChild, addRightChild,
						deductParent, leftChildCount, rightChildCount,
						deltaPiTentative,
						logging, sHist, i);

			if (willSplit) {
				// take the target out of the list
				nodes.erase(it);
				numLeaves--;
				
				/* update the info before we actually change the states
				 * the tentative change in log posterior 
				 * will be the actual change */
				deltaPi = deltaPiTentative;
				
				// try to change the state according to the proposed split
				// nodes, numLeaves and numCherries are passed by reference
				changeStateForSplitNEW(target, nodes,
							numLeaves, numCherries,
							addLeftChild, addRightChild, deductParent);
							
				/* and also have to change the state of the rmsp,
				 * which means we have to find equivalent node in the leaves
				 * and implement the change*/
				changeStateForSplitRMSP(target,
							leaves, cherries, 
							leftChildCount, rightChildCount,
							info);
				
			} // end of willSplit

		} // end of isLeaf

		else { // it points to a cherry

			// grab the cherry
			SPSnode* target = *it; // don't change where it points until deleted
			
			int deductLeftChild = 0;
			int deductRightChild = 0;
			int addParent = 0;

			/* deductLeftChild, deductRightChild, addParent, deltaPiTentative
			 * all passed by reference,
			 * real number of leaves is leaves.size() */
			bool willMerge = decisionMCMCMergeNEW(target, proposal, logPrior, rgsl,
							numLeaves, numCherries, 
							leaves,
							minPoints, 
							volChecking,
							minVol,
							deductLeftChild, deductRightChild,
							addParent,
							deltaPiTentative,
							logging, sHist, i);

			if (willMerge) {
				// take the target out of the list of cherries
				nodes.erase(it);
				numCherries--;

				changeStateForMergeNEW(target,
							nodes, numLeaves, numCherries,
							deductLeftChild, deductRightChild,
							addParent);
				
				/* and also have to change the state of the rmsp,
				 * which means we have to find equivalent node in the cherries
				 * and implement the change*/
				changeStateForMergeRMSP(target,
							leaves, cherries,
							info);
				
				// the tentative change in log posterior will be the actual change
				deltaPi = deltaPiTentative;

			} // end willMerge
		} // end if cherry
	} // end haveNode

	// value in deltaPi should have been updated IF we changed

	else if (((logging == TXT) || (logging == TXTANDGRAPH)) && !haveNode) {

		#if defined OLDLOGSTYLE
			std::string line = "No node grabbed (possible if proposal has fixed ";
            line += "probability of split): state stays the same";
            outputFile(sHist, line, true); // append to log file
		#else
			outputFile(sHist, i, 0);
		#endif
	}

	bool success = (nodes.size() == (numLeaves+numCherries));
	
	assert( numCherries == cherries.size() );
	/*It is not necessarily true that leaves.size() == numLeaves
	 * because if leaves are considered unsplittable because of minPoints
	 * they will be in leaves but not counted as numLeaves. */ 

	if (!success) {
		std::string errorMsg("AdaptiveHistogram::changeMCMCState(");
		errorMsg += "SPSnodeList&, size_t&, size_t&, ";
		errorMsg += "MCMCProposal&, LogMCMCPrior&, size_t, ";
		errorMsg += "RealMappedSPnode&, ";
		errorMsg += "RealMappedSPnode::Ptrs&, ";
		errorMsg += "RealMappedSPnode::Ptrs&, ";
		errorMsg +=  "gsl_rng*, LOGGING_LEVEL, std::string, int)";
		errorMsg += "\n: Nodes muddled in changeMCMCstateNEW";
		throw std::logic_error(errorMsg);
	}
	
	info.notifyDeltaPi(deltaPi);
	
	return info;
	
}


/* support posteriors and posteriors
 * used in routines like seb-carving queue */
/* return true if there are still nodes to split when routine terminates, 
	false if run out of splittable nodes and this forced end 
	of pq splitting */
bool AdaptiveHistogram::prioritySplitWithSupportPosterior(
								const PrioritySplitQueueEvaluator& psqe,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints,
								LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								std::vector<real>& posteriorVec,
								gsl_rng * rgsl,
                                bool shiftCatalan)
{
	double minVol = 0.0;
	return prioritySplitWithSupportPosterior(
								psqe,
                                logging,
                                minChildPoints,
								minVol, 
                                logPrior,
                                emptyBoxVolsVec,
                                posteriorSupportVec,
								posteriorVec,
								rgsl,
                                shiftCatalan);
}

bool AdaptiveHistogram::prioritySplitWithSupportPosterior(
								const PrioritySplitQueueEvaluator& psqe,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								std::vector<real>& posteriorVec,
								gsl_rng * rgsl,
                                bool shiftCatalan)
{
	string errorMsg("AdaptiveHistogram::prioritySplitWithSupportPosterior(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	size_t n = getRootCounter(); // for number of points in histogram
	
	real volRoot = getRootVolume();

	assert(volRoot > 0.0);

	// a multiset for the queue (key values are not necessarily unique)
	PriorityQueueT pq;
	
	PrivatePrioritySplitQueueEvaluator ppsqe(psqe);

	/* make sure PQ has non-empty cells and go through all
	leaves to accumulate vol of empty cells
	to accumulate the volume of empty boxes */
	
	// set up the priority queue
	real emptyBoxVol = _setupPrioritySplitWithEmptyVolMeasure(pq, 
					ppsqe.measurer, minChildPoints, minVol);
	// throws an exception if current state not legal
	
	size_t numLeaves = getRootLeaves(); // actual number of leaves
	
	// prior:  if we shifted we'd have Catalan for 0 = 1
	real lnPrior(1.0); 
	if ( !shiftCatalan ) lnPrior = logPrior(numLeaves);

	// log-posterior computations
	real lnLik = getLogLikelihood();
			
		// accumulate the ratios of the empty volumes
	real cumEmptyBoxVol= (emptyBoxVol/volRoot);
	emptyBoxVolsVec.push_back(cumEmptyBoxVol);
	
	// log-posterior support computations
	// log-prior
	// use the prior distribution object to find the prior 
	
	// likelihood 
	real invSupportVol(1.0/(volRoot - emptyBoxVol));
	real lnLikSupport = (1.0*n) * cxsc::ln(invSupportVol);
	
	// posterior is proportional to likelihood * prior
	real posteriorSupport = lnLikSupport + lnPrior;
	// store
	posteriorSupportVec.push_back(posteriorSupport);
	
	// actual posterior
	real posterior = lnLik + lnPrior;
	posteriorVec.push_back(posterior);

	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	// a counter for the loop
	int i = 0;
	
	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		
	}
	i++;
	
	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	
	/* split until the ppsqe says to stop */ 
	while (canContinue && !(ppsqe.stopQueueQuery(pq, numLeaves)) ) {

		real loopEmptyBoxVol(0.0);
		real deltaL(0.0);
		
		canContinue = _prioritySplitLoopWithEmptyVolMeasure(pq, 
						ppsqe.measurer,
						n, minChildPoints, minVol,
						loopEmptyBoxVol, deltaL, rgsl);

		emptyBoxVol += loopEmptyBoxVol;
		cumEmptyBoxVol = emptyBoxVol/volRoot;
		emptyBoxVolsVec.push_back(cumEmptyBoxVol);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			
		}
		i++;

		/* if loop was successful number of leaves will increased by 1*/
		if (canContinue) {
			real deltaP = logPrior.changeOnSplitOne(numLeaves-1);
			numLeaves++;
			
			// update the logprior
			lnPrior += deltaP;
			
			// update the support likelihood
			invSupportVol = 1.0/(volRoot - emptyBoxVol);
			lnLikSupport = (1.0*n) * cxsc::ln(invSupportVol);
			
			//update the support posterior
			posteriorSupport = lnLikSupport + lnPrior;
			// store
			posteriorSupportVec.push_back(posteriorSupport);
		
			// update the 'real' likelihood
			lnLik += deltaL;
			// and posterior
			posterior = lnLik + lnPrior;
			posteriorVec.push_back(posterior);
	
		}

		if (canContinue && (logging != NOLOG)) {
			// log the leaf levels line
			outputFile(s, getLeafLevelsString());
		}
	}
	// EMPSums will have been adjusted during the splitting process
	
	return canContinue; 
	/* true if there are still nodes to split, 
	false if run out of splittable nodes and this forced end 
	of pq splitting */
}

/* Find maximum log-posterior point(s) and keep record.
 * Used in routines like seb-carving queue. */
/* Return true if routines launched from various step points got
	to the end of the required queue without running out of leaves*/
/* stopOnMaxPosterior = true does not give immediate stop,
 * still makes sure
* that the posterior has fallen below the starting point of the
* seb queue and that a specfied minimum number of states have  elapsed*/
bool AdaptiveHistogram::prioritySplitWithSupportPosteriorMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
								LOGGING_LEVEL logging,
                                size_t minChildPoints, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								int checkMaxStep,
								bool stopOnMaxPosterior,
								std::vector<size_t>& carvedLaunchPoints,
								std::vector<size_t>& maxPosteriorPoints,
								std::vector<real>& maxPosteriors,
                                gsl_rng * rgsl,
                                bool shiftCatalan)
{
	double minVol = 0.0;
	return prioritySplitWithSupportPosteriorMaxLik(
								psqe,
								psqePosterior,
								logging,
                                minChildPoints,
								minVol, 
                                logPrior,
                                emptyBoxVolsVec,
                                posteriorSupportVec,
								checkMaxStep,
								stopOnMaxPosterior,
								carvedLaunchPoints,
								maxPosteriorPoints,
								maxPosteriors,
                                rgsl,
                                shiftCatalan);
	
}

bool AdaptiveHistogram::prioritySplitWithSupportPosteriorMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
								LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								int checkMaxStep, // if 0, no carving
								bool stopOnMaxPosterior,
								std::vector<size_t>& carvedLaunchPoints,
								std::vector<size_t>& maxPosteriorPoints,
								std::vector<real>& maxPosteriors,
                                gsl_rng * rgsl,
                                bool shiftCatalan)
{
	string errorMsg("AdaptiveHistogram::prioritySplitWithSupportPosteriorMaxLik(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	#ifdef MYDEBUG_MCMCPQ_CHECK_QUEUE
		cout << "\n\nEntering prioritySplitWithSupportPosteriorMaxLik " << endl;
		gsl_rng * rgsl_tmp = gsl_rng_clone (rgsl);
		
		cout << "The next rand would be " << gsl_rng_uniform(rgsl_tmp) << endl;
	
		gsl_rng_free(rgsl_tmp);
	
	#endif
	
	//make sure empty
	std::vector<size_t>().swap(carvedLaunchPoints);
	std::vector<size_t>().swap(maxPosteriorPoints);
	std::vector<real>().swap(maxPosteriors);
	
	size_t n = getRootCounter(); // for number of points in histogram
	
	real volRoot = getRootVolume();

	assert(volRoot > 0.0);

	// a multiset for the queue (key values are not necessarily unique)
	PriorityQueueT pq;

	PrivatePrioritySplitQueueEvaluator ppsqe(psqe);
	
	// set up the priority queue
	real emptyBoxVol = _setupPrioritySplitWithEmptyVolMeasure(pq, 
							ppsqe.measurer, minChildPoints, minVol);
	// throws an exception if current state not legal
	
	size_t numLeaves = getRootLeaves(); // actual number of leaves
	
	// prior:  if we shifted we'd have Catalan for 0 = 1
	real lnPrior(1.0); 
	if ( !shiftCatalan ) lnPrior = logPrior(numLeaves);

	// log-posterior computations
	real lnLik = getLogLikelihood();
			
	// accumulate the ratios of the empty volumes
	real cumEmptyBoxVol= (emptyBoxVol/volRoot);
	emptyBoxVolsVec.push_back(cumEmptyBoxVol);
	
	// log-posterior support computations
	// log-prior
	// use the prior distribution object to find the prior 
	
	// likelihood 
	real invSupportVol(1.0/(volRoot - emptyBoxVol));
	real lnLikSupport = (1.0*n) * cxsc::ln(invSupportVol);
	// posterior is proportional to likelihood * prior
	real posteriorSupport = lnLikSupport + lnPrior;
	// store
	posteriorSupportVec.push_back(posteriorSupport);

	//make copies of all our starting points
	PriorityQueueT stepPQ(pq);
	real stepPosteriorSupport = posteriorSupport;
	size_t stepNumLeaves = numLeaves;
	real stepVolRoot = volRoot;
	real stepEmptyBoxVol = emptyBoxVol;
	real stepLnPrior = lnPrior;
	real stepLnLik = lnLik;
	gsl_rng * rgsl_step = gsl_rng_clone (rgsl);
					
	
	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	// a counter for the loop
	int i = 0;
	// and for steps
	int stepCounter = 0;
	int withinStepCounter = 1;

	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		
	}
	i++;
	
	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	else { // have some leaves
	
		carvedLaunchPoints.push_back(stepNumLeaves);
				
		// want to try the posterior right from root
		cout << "\nlaunching full posterior check from root :" << endl;
		cout << "stepNumLeaves " << stepNumLeaves 
			<< " and posterior support " << stepPosteriorSupport
			<< " and posterior " << (stepLnLik + stepLnPrior) << endl;
		cout << endl;
		
		// always do the one from the root all the way through
		bool rootStopOnMaxPosterior = false;
		vector<real> tmpPosteriorVec;
		vector<real> tmpLoglikVec;
		
		PriorityQueueT rootPQ(pq);
		gsl_rng * rgsl_root = gsl_rng_clone (rgsl);
		_launchPrioritySplitWithPosterior( // ignore return value
					rootPQ,
					psqePosterior,
					logging,
					minChildPoints,
					minVol, 
					logPrior,
					tmpPosteriorVec,
					tmpLoglikVec,
					n,
					stepNumLeaves,
					stepLnPrior,
					stepLnLik,
					rootStopOnMaxPosterior,
					maxPosteriorPoints,
					maxPosteriors,
					rgsl_root);
		
		gsl_rng_free(rgsl_root);		
	}
	
	/* if checkMaxStep > 0, 
	 * split until the ppsqe says we can stop */ 
	while ((checkMaxStep > 0) && canContinue && !(ppsqe.stopQueueQuery(pq,numLeaves)) ) {

		real loopEmptyBoxVol(0.0);
		real deltaL(0.0);
		
		canContinue = _prioritySplitLoopWithEmptyVolMeasure(pq, 
								ppsqe.measurer,
								n, minChildPoints, minVol,
								loopEmptyBoxVol, deltaL, rgsl);

		emptyBoxVol += loopEmptyBoxVol;
		cumEmptyBoxVol = emptyBoxVol/volRoot;
		emptyBoxVolsVec.push_back(cumEmptyBoxVol);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			
		}
		i++;
		
		// is there anything left 
		canContinue = (volRoot - emptyBoxVol > cxsc::MinReal);
		
		/* if loop was successful number of leaves will increased by 1*/
		if (canContinue) {
			real deltaP = logPrior.changeOnSplitOne(numLeaves-1);
			numLeaves++;
			withinStepCounter++;
			
			// update the logprior
			lnPrior += deltaP;
			
			// update the support likelihood
			invSupportVol = 1.0/(volRoot - emptyBoxVol);
			lnLikSupport = (1.0*n) * cxsc::ln(invSupportVol);
			
			//update the support posterior
			real newPosteriorSupport = lnLikSupport + lnPrior;
			// store
			posteriorSupportVec.push_back(newPosteriorSupport);
		
			// update the 'real' likelihood
			lnLik += deltaL;
			
			#ifdef MYDEBUG_MCMCPQ
				cout << "Cqueue log"
					<< "\tqueue\t" << pq.size() << "\tleaves\t" << numLeaves 
					<< "\treal post\t" << (lnLik + lnPrior) << endl; 
			#endif
	
			//need to update if we are a new max in the step or at the very start of the step
			if ( ( (withinStepCounter <= checkMaxStep) &&
					(newPosteriorSupport > posteriorSupport)) || (withinStepCounter == 1)) {
			
				#ifdef MYDEBUG_MCMCPQ_EXTRA
					cout << "\nresetting at step " << stepCounter << " : " << withinStepCounter << endl;
					cout << "size of set is " << pq.size() << endl;
				#endif
				
				stepPQ = pq;
								
				assert(stepPQ.size() == pq.size());
				#ifdef MYDEBUG_MCMCPQ_CHECK_QUEUE
					PriorityQueueItrT sit = stepPQ.begin();
					int j = 1;
					for (PriorityQueueItrT it = pq.begin();
						it != pq.end(); ++it, ++sit, ++j) {
							if (it->nodePtr->getNodeName() != sit->nodePtr->getNodeName()) {
								cout << "\tsets differ" << endl;
								
								PriorityQueueItrT ssit = stepPQ.begin();
								for (PriorityQueueItrT iit = pq.begin();
								iit != pq.end(); ++iit, ++ssit) {
									cout << iit->toString() 
									<<"\t" << ssit->toString()
									<<endl;
								}
								cout << endl;
							}
							break;	
					}
					
				#endif
									
				stepPosteriorSupport = newPosteriorSupport;
				stepNumLeaves = numLeaves;
				stepVolRoot = volRoot;
				stepEmptyBoxVol = emptyBoxVol;
				stepLnPrior = lnPrior;
				stepLnLik = lnLik;
				gsl_rng_memcpy (rgsl_step, rgsl);
				
			}
			
			posteriorSupport = newPosteriorSupport;

			//if are at the end of a step we find who was the max
			if (withinStepCounter == checkMaxStep) {
				
				carvedLaunchPoints.push_back(stepNumLeaves);
				
				#ifdef MYDEBUG_MCMCPQ_MIN
					cout << "\nlocal max in step " << stepCounter << ":" << endl;
					cout << "stepNumLeaves " << stepNumLeaves 
						<< " and posterior support " << stepPosteriorSupport
						<< " and posterior " << (stepLnLik + stepLnPrior) << endl;
					cout << endl;
				#endif
				
				// can't do this if we are debugging looking at the continuation of the original pq process
				#ifndef MYDEBUG_MCMCPQ_CONTINUE 
				{
					vector<real> tmpPosteriorVec;
					vector<real> tmpLoglikVec;
					
					_launchPrioritySplitWithPosterior( // ignore return value
							stepPQ, // okay, reset immediately after
							psqePosterior,
							logging,
							minChildPoints,
							minVol, 
							logPrior,
							tmpPosteriorVec,
							tmpLoglikVec,
							n,
							stepNumLeaves,
							stepLnPrior,
							stepLnLik,
							stopOnMaxPosterior,
							maxPosteriorPoints,
							maxPosteriors,
							rgsl_step);
				}	
				#endif
				
				// Debugging looking at the continuation of the original pq process
				#ifdef MYDEBUG_MCMCPQ_CONTINUE
				{
					vector<real> tmpPosteriorVec;
					vector<real> tmpLoglikVec;
					
					cout << "size of set is " << stepPQ.size() << endl;
					
					bool stopOnMaxPosterior = false;
					size_t leavesForMaxPosterior = 0;
					real maxPosterior(0.0);
														
					bool success = _prioritySplitWithPosterior(
								stepPQ, ppsqe, logging,
								minChildPoints, minVol, 
								logPrior,
								tmpPosteriorVec,
								tmpLoglikVec,
								n, stepNumLeaves,
								stepLnPrior, stepLnLik,
								stopOnMaxPosterior,
								leavesForMaxPosterior,
								maxPosterior,
								rgsl_step);
					
					cout << "\n\nBack in original\n" << endl;
					
					string postFileName = "PosteriorsJennyClone.txt";
					std::vector < const subpavings::RealVec* > dataPtrs;
					dataPtrs.push_back(&tmpPosteriorVec, &tmpLoglikVec);
					std::vector < std::string > colNames;
					colNames.push_back("Posterior");
					colNames.push_back("LnLik");
					
					outputToFileVertical(dataPtrs, colNames, postFileName);
					
				}
				#endif		
				
				#ifdef MYDEBUG_MCMCPQ
					cout << "\nBack in original, resetting for start of new step" << endl;
					cout << "size of this queue is\t" << pq.size() << "\tnumLeaves\t" << numLeaves << endl; 
				#endif
				
				withinStepCounter = 0; // reset the within step counter
								
				++stepCounter;
			}
		}
		
		#ifdef MYDEBUG_MCMCPQ_MIN
			if (!canContinue) 
				cout << "\nCan't continue: leaves " << numLeaves << endl;
		#endif
		
		if (canContinue && (logging != NOLOG)) {
			// log the leaf levels line
			outputFile(s, getLeafLevelsString());
		}
	}
	// EMPSums will have been adjusted during the splitting process
	try {
		gsl_rng_free(rgsl_step);
	}
	catch(...) {} // catch and swallow
	
	return canContinue; 
	/* true if routines launched from various carving points got
	to the end of the required queue without running out of leaves*/
}

//external
// Refind a max point
/* Return true if there were splittable leaves at end of process, false
 * if process terminated because we run out of splittable leaves */
bool AdaptiveHistogram::prioritySplitMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
					            LOGGING_LEVEL logging,
                                size_t minChildPoints, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& posteriorVec,
								std::vector<real>& loglikVec,
                                gsl_rng * rgsl,
                                bool shiftCatalan)
{
	double minVol = 0.0;
	return AdaptiveHistogram::prioritySplitMaxLik(
								psqe,
								psqePosterior,
					            logging,
                                minChildPoints,
								minVol, 
                                logPrior,
                                posteriorVec,
								loglikVec,
                                rgsl,
                                shiftCatalan);
}

bool AdaptiveHistogram::prioritySplitMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
					            LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& posteriorVec,
								std::vector<real>& loglikVec,
                                gsl_rng * rgsl,
                                bool shiftCatalan)
{
	#ifdef MYDEBUG_MCMCPQ_CHECK_QUEUE
		cout << "\n\nEntering prioritySplitMaxLik " << endl;
		gsl_rng * rgsl_tmp = gsl_rng_clone (rgsl);
		
		cout << "The next rand would be " << gsl_rng_uniform(rgsl_tmp) << endl;
	
		gsl_rng_free(rgsl_tmp);
	
	#endif
	
	string errorMsg("AdaptiveHistogram::prioritySplitMaxLik()");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	size_t numLeaves = getRootLeaves(); // actual number of leaves
	
	if (psqe.getMaxLeaves() < numLeaves) {
		
		throw std::invalid_argument(errorMsg 
					+ ": psqe.getMaxLeaves() < numLeaves");
	}
	if (psqePosterior.getMaxLeaves() < psqe.getMaxLeaves()) {
		
		throw std::invalid_argument(errorMsg 
					+ ": psqePosterior.getMaxLeaves() < psqe.getMaxLeaves()");
	}
	
	size_t n = getRootCounter(); // for number of points in histogram
	
	// a multiset for the queue (key values are not necessarily unique)
	PriorityQueueT pq;

	PrivatePrioritySplitQueueEvaluator ppsqe(psqe);

	real emptyBoxVol = _setupPrioritySplitWithEmptyVolMeasure(pq, 
									ppsqe.measurer, minChildPoints, minVol);
	// throws an exception if current state not legal
	
	// prior:  if we shifted we'd have Catalan for 0 = 1
	real lnPrior(1.0); 
	if ( !shiftCatalan ) lnPrior = logPrior(numLeaves);

	// log-posterior computations
	real lnLik = getLogLikelihood();
			
	// posterior is proportional to likelihood * prior
	real posterior = lnLik + lnPrior;
	// store
	posteriorVec.push_back(posterior);
	loglikVec.push_back(lnLik);

	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	if (logging == LOGSAMPLES) {
		ostringstream oss;
		oss << "CarverCarvingQueueState_" << numLeaves << ".txt";
		outputToTxtTabs(oss.str());
		
	}
	
	// a counter for the loop
	int i = 0;
	
	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		
	}
	i++;
	
	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
		
	/* split until ppsqe says we can stop */ 
	while (canContinue && !( ppsqe.stopQueueQuery(pq, numLeaves) )) {

		real loopEmptyBoxVol(0.0);
		real deltaL(0.0);
		
		canContinue = _prioritySplitLoopWithEmptyVolMeasure(
								pq, ppsqe.measurer,
								n, minChildPoints, minVol,
								loopEmptyBoxVol, deltaL, rgsl);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			
		}
		i++;

		/* if loop was successful number of leaves will increased by 1*/
		if (canContinue) {
			real deltaP = logPrior.changeOnSplitOne(numLeaves-1);
			numLeaves++;
			
			// update the logprior
			lnPrior += deltaP;
			
			// update the 'real' likelihood and posterior
			lnLik += deltaL;
			posterior = lnLik + lnPrior;
			posteriorVec.push_back(posterior);
			loglikVec.push_back(lnLik);
			
			#ifdef MYDEBUG_MCMCPQ
				cout << "MLC queue log"
					<< "\tqueue\t" << pq.size() << "\tleaves\t" << numLeaves 
					<< "\treal post\t" << posterior << endl; 
			#endif
			
			if (logging == LOGSAMPLES) {
				ostringstream oss;
				oss << "CarverCarvingQueueState_" << numLeaves << ".txt";
				outputToTxtTabs(oss.str());
				
			}

		}

		if (canContinue && (logging != NOLOG)) {
			// log the leaf levels line
			outputFile(s, getLeafLevelsString());
		}
	}
	
	#ifdef MYDEBUG_MCMCPQ_EXTRA
		cout << "\nend of carving part, set size is " << pq.size() << endl;
		cout << "numLeaves is  " << numLeaves << endl;
	#endif
	#ifdef MYDEBUG_MCMCPQ_CHECK_QUEUE
		cout << "queue at this stage is " << endl;
		for (PriorityQueueItrT it = pq.begin();
			it != pq.end(); ++it) {
				cout << it->toString() << endl;
			}
	#endif
				
	
	if(pq.empty()) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	if (canContinue && !pq.empty()) {
		// now do the other bit
		
		// dummy variables - results not used
		bool stopOnMaxPosterior = false;
		size_t leavesForMaxPosterior = 0; 
		real maxPosterior(0.0); 
		
		/* remake the queue with the new measurer */
		PrivatePrioritySplitQueueEvaluator ppsqePosterior(psqePosterior);
		
		PriorityQueueT pqNew;
		
		_resetPrioritySplitWithEmptyVolMeasure(pq, pqNew, 
						ppsqePosterior.measurer);
		
		/*need to take the last recorded value out of the posteriors
		and log liks because it will get recorded again */
		posteriorVec.resize(posteriorVec.size()-1);
		loglikVec.resize(loglikVec.size()-1);
			
		canContinue = _prioritySplitWithPosterior(
						pqNew, ppsqePosterior,
						logging, minChildPoints, minVol,
						logPrior,
						posteriorVec, 
						loglikVec,
						n, numLeaves,
						lnPrior, lnLik, stopOnMaxPosterior,
						leavesForMaxPosterior, maxPosterior,
						rgsl);
		
	}
	
	return canContinue; 
}

/* does the full posterior using just one type of queue.
 * Return true if there were splittable leaves at end of process,
 * false if process terminated due to lack of splittable leaves */
bool AdaptiveHistogram::prioritySplitWithPosterior(
						const PrioritySplitQueueEvaluator& psqe,
						LOGGING_LEVEL logging,
						size_t minChildPoints,
						LogMCMCPrior& logPrior,
						std::vector<real>& posteriorVec,
						std::vector<real>& loglikVec,
						gsl_rng * rgsl,
						bool shiftCatalan)
{
	double minVol = 0.0;
	return prioritySplitWithPosterior(
						psqe,
						logging,
						minChildPoints,
						minVol, 
						logPrior,
						posteriorVec,
						loglikVec,
						rgsl,
						shiftCatalan);
}


bool AdaptiveHistogram::prioritySplitWithPosterior(
						const PrioritySplitQueueEvaluator& psqe,
						LOGGING_LEVEL logging,
						size_t minChildPoints,
						double minVol, 
						LogMCMCPrior& logPrior,
						std::vector<real>& posteriorVec,
						std::vector<real>& loglikVec,
						gsl_rng * rgsl,
						bool shiftCatalan)
{
	string errorMsg("AdaptiveHistogram::prioritySplitWithPosterior(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}

	size_t n = getRootCounter(); // for number of points in histogram
	
	// a multiset for the queue (key values are not necessarily unique)
	PriorityQueueT pq;

	PrivatePrioritySplitQueueEvaluator ppsqe(psqe);
	
	// set up the priority queue
	real emptyBoxVol = _setupPrioritySplitWithEmptyVolMeasure(pq, 
													ppsqe.measurer,
													minChildPoints,
													minVol);
	// throws an exception if current state not legal

	size_t numLeaves = getRootLeaves(); // actual number of leaves
	
	// log-posterior support computations
	// log-prior
	// use the prior distribution object to find the prior 
	
	// prior:  if we shifted we'd have Catalan for 0 = 1
	real lnPrior(1.0); 
	if ( !shiftCatalan ) lnPrior = logPrior(numLeaves);

	// log-posterior computations
	real lnLik = getLogLikelihood();
	
	size_t leavesForMaxPosterior = 0;
	bool stopOnMaxPosterior = false;
	real maxPosterior(0.0);
	
	bool success = _prioritySplitWithPosterior(
					pq,
					ppsqe,
					logging,
					minChildPoints, 
					minVol,
					logPrior,
					posteriorVec,
					loglikVec,
					n,
					numLeaves,
					lnPrior,
					lnLik,
					stopOnMaxPosterior,
					leavesForMaxPosterior,
					maxPosterior,
					rgsl);

	
	return success; 
}



// uses getUnscaledTreeLogLik() method in SPSnode class
cxsc::real AdaptiveHistogram::getLogLikelihood() const
{
	if ( !hasSubPaving()) {
		throw NullSubpavingPointer_Error(
				"AdaptiveHistogram::getLogLikelihood()");
	}
	cxsc::real loglik(0.0);
	size_t totalPoints = getRootCounter();
	if (totalPoints > 0) {
		
		// unscaled 
		loglik = getSubPaving()->getUnscaledTreeLogLik();
		
		/* scaling is -n x ln (n x vol) where n is 
		the total points and vol is the volume of the root box */
		loglik -= (1.0*totalPoints * cxsc::ln(
				getSubPaving()->nodeRealVolume() * (1.0*totalPoints)));
	}
	
	return loglik;
	
}


// returns a vector of leaf levels as ints
// left to right, 0 is root
IntVec AdaptiveHistogram::getLeafLevels() const
{
    IntVec levels; // empty container

    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeLevels(levels, 0);
        //levels has now been filled in
    }
    return levels;
}


// returns a vector of leaf counts
// left to right
Size_tVec AdaptiveHistogram::getLeafCounts() const
{
    Size_tVec counts; // empty container
    if (hasSubPaving()) {
        getSubPaving()->getLeafNodeCounts(counts);
        //levels has now been filled in
    }
    return counts;
}


// make a .dot file for the histogram
bool AdaptiveHistogram::outputGraphDot() const
{
    bool success = false;

    if (hasSubPaving()) {
        success = getSubPaving()->outputGraphDot();

    }
    else {
        std::cerr << "Sorry, you can't make a graph without a root paving"
                << std::endl;
    }
    return success;
}

// Method to output the subpaving to a txt file
std::ostream & AdaptiveHistogram::outputToStreamTabs(std::ostream & os,
								int prec) const
{
	if (hasSubPaving()) {

		 return getSubPaving()->leavesOutputTabs(os, prec); // the output
	}
	
    else return os;
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void AdaptiveHistogram::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void AdaptiveHistogram::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the AdaptiveHistogram object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {

			getSubPaving()->leavesOutputTabs(os, prec); // the output
		}
		if (confirm)
			std::cout << "The output of the AdaptiveHistogram "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
// Output includes scaled EMP contributions under COPERR and AIC
// and changes if split

void AdaptiveHistogram::outputToTxtTabsWithEMPs(const std::string& s,
                                                    int prec) const
{
	outputToTxtTabsWithEMPs(s, prec, false);
}

void AdaptiveHistogram::outputToTxtTabsWithEMPs(const std::string& s,
										int prec, bool confirm) const
{
	// To generate a file output of the AdaptiveHistogram object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		if (hasSubPaving()) {
			size_t n = getRootCounter();
			getSubPaving()->leavesOutputTabsWithEMPs(n, os, prec); // the output
		}
		if (confirm)
			std::cout << "The output of the AdaptiveHistogram with scaled EMPs "
                    << "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
	
}

void AdaptiveHistogram::outputRootToTxt(const std::string& s,
										int prec) const
{
	outputRootToTxt(s, prec, false);
}

// Method to output details and stats on the root paving to a txt file
// Output goes to file named according to arguement s
void AdaptiveHistogram::outputRootToTxt(const std::string& s,
										int prec, bool confirm) const
{
 	// To generate a file output of root node of the AdaptiveHistogram
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
	
		if (hasSubPaving()) {
			
			os << cxsc::SaveOpt;
			os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
		
			getSubPaving()->nodePrint(os); // the output
			
			os << cxsc::RestoreOpt;
		}
		if (confirm)
			std::cout << "Details of the root paving of the AdaptiveHistogram "
				<< "has been written to " << s << std::endl << std::endl;
				
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

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

// get distance between two histograms
// throws exception if this has null paving
// spsnode routines checks for null pointer in other, and empty pavings, and unequal boxes
// labels do not have to match
cxsc::real AdaptiveHistogram::getL1Distance(const AdaptiveHistogram& other) const
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
			"AdaptiveHistogram::getL1Distance(const AdaptiveHistogram&)");
	}
	
	return getSubPaving()->getL1Distance(other.getSubPaving());
		
}

void AdaptiveHistogram::setHoldAllStats(bool setTo)
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
			"AdaptiveHistogram::setHoldAllStats(bool)");
	}

	if (getHoldAllStats() != setTo) {
		
		//if setTo is false, paving's countsOnly needs to be true
		//if setTo is true, paving's countsOnly needs to be false
		getSubPaving()->setCountsOnly(!setTo);
		holdAllStats = setTo;
		
	}
}

void AdaptiveHistogram::setLabel(int lab)
{
	label = lab;
}

// clear all data from the histogram
void AdaptiveHistogram::clearAllHistData()
{
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(
			"AdaptiveHistogram::clearAllHistData()");
	}

	getSubPaving()->clearAllDataHeld();
	
	BigDataCollection tmp;
	tmp.swap(dataCollection); // swap and clear
	
	scaledEMPSumCOPERR = cxsc::dotprecision(0.0);
	scaledEMPSumAIC = scaledEMPSumCOPERR;
	
}

void AdaptiveHistogram::swap(AdaptiveHistogram& adh) // throw()
{
	std::swap(label, adh.label);
	std::swap(dataCollection, adh.dataCollection); // use stl specialisation of swap
    std::swap(holdAllStats, adh.holdAllStats);
	
	// cxsc don't seem to have a swap for dot precisions
    dotprecision tempCOPERR(adh.scaledEMPSumCOPERR);
	adh.scaledEMPSumCOPERR = scaledEMPSumCOPERR;
	scaledEMPSumCOPERR = tempCOPERR;
	
	dotprecision tempAIC(adh.scaledEMPSumAIC);
	adh.scaledEMPSumAIC = scaledEMPSumAIC;
	scaledEMPSumAIC = tempAIC;
	
	std::swap(rootPaving, adh.rootPaving); // just swap the pointers
}

std::string AdaptiveHistogram::stringSummary() const 
{
	std::ostringstream oss;
	
	oss << "This address = " << (this) << endl;
	oss << "Label = " << label << endl;
	oss << "holdAllStats = " << getHoldAllStats() << endl;
	oss << "AIC score = " << getEMPScoreAIC() << ", and COPERR score = " << getEMPScoreCOPERR() << endl;
	oss << "addresss of dataCollection is " << &dataCollection << " and dataCollection is:" << endl;
	
	std::string delim = "\t";
	for (BigDataConstItr it = dataCollection.begin(); it != dataCollection.end(); ++it) {
		prettyPrint(oss, (*it));
		oss << delim;
	}
	oss << endl;
	
	if (hasSubPaving()) oss << "Address of subpaving is " << getSubPaving() << endl;
	else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
}

//for debugging
std::string AdaptiveHistogram::doubleCheckStringSummary() const
{
	std::ostringstream oss;
	
	oss << "This address = " << this << endl;
	oss << "Label = " << label << endl;
	oss << "holdAllStats = " << getHoldAllStats() << endl;
	oss << "address of data collection is " << &dataCollection << " and data collection addresses are: " << endl;
		std::string delim = "\t";
		for (BigDataConstItr it = dataCollection.begin(); it != dataCollection.end(); ++it) {
			oss << &(*it) << delim;
		}
		oss << endl;
		if (hasSubPaving()) {
			oss << "Subpaving is: " << endl;
			oss << getSubPaving()->doubleCheckStringSummary();
		}
		else oss << "Subpaving is NULL" << endl;
	
	return oss.str();
	
}

// --------------------------- private ---------------------------------------

// a constant for padding a box if it is tailor-made for data
const real AdaptiveHistogram::padding = 0.000005;



// complete insertion of data from a vector of data
// given a container of rvectors of the data  to insert
// return true if at least some data went into the actual subpaving itself
// else false
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
        if(dataDim != getDimensions()) {

            throw IncompatibleDimensions_Error(
            "AdaptiveHistogram::completeDataInsertionFromVec(const RVecData&, const SplitDecisionObj&, LOGGING_LEVEL)");
        }
    }

    // insert the data
    size_t dataCountInserted
            = insertDataFromContainer(theData, boolTest, logging);

    retValue = (dataCountInserted); // will be true if > 1
	
	    // switch on for more output during histogram creation "
        /*
        std::cout << "End of inserting data: " << dataCountInserted
            << " data points inserted to dataCollection "
            << std:: endl;
        std::cout << "and associated with the tree " << std::endl;
        std::cout << "(check console output for possible records "
            << "of datapoints which did not fit)" << std::endl;
        */
	
    return retValue;
}


// check if we need to make a paving for the histogram object
// make it if we need to, matching the dimensions of the data
// as given in function argument
// return true if needed to make the paving
bool AdaptiveHistogram::haveMadePaving(const RVecData& theData,
                                    const size_t dim)
{

    SPSnode* newroot = NULL;
	
	try {
		bool retValue = false;
	
		if (!hasSubPaving() || getSubPaving()->isEmpty()) { // NULL or no box

			ivector rootBox = makeBox(theData, dim);
			
			// if this subpaving is NULL, we want to make the new one
			
			if (!hasSubPaving()) {

				newroot = new SPSnode(rootBox, !holdAllStats);
				rootPaving = newroot;
				// two pointers to the SPSnode now, but newroot is a temp
			}
			/* this case should not happen with the class as it is at 
			present, but best to cover it anyway I guess ...*/
			else if (getSubPaving()->isEmpty() ) // not null but no box
			{
				/* make a temp, then swap - temp goes out of scope at
				 * the end of this method and is deleted normally */
				SPSnode temp(rootBox, !holdAllStats);
				std::swap(*rootPaving, temp); 
			}

			retValue = true;
		}
		return retValue;
		// end of making the subpaving if there was not one
	}
	catch (std::exception const& e) {
		
		try {
			delete newroot;
			newroot = NULL;
		}
		catch (std::exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
	}
}




// make a box to fit all the data
ivector AdaptiveHistogram::makeBox(const RVecData& theData, const size_t dim)
{
    ivector retVal = subpavings::makeBox(theData, static_cast<int>(dim),
					padding);    

    // and make each interval the (min, max) of the corresponding elements
    // of the rvectors -/+ some padding

    std::cout << "A box is being made for the data.  "
        << "The box is " << std::endl;  // standard output message

    // make intervals and make them elements of the ivector
    for (size_t i = 1; i <=dim; i++) {
        std::cout << retVal[i] << "  ";    // output
    }
    std::cout << std::endl;
    
    // check the box here - should not be necessary but do for completeness
	if (!checkBox(retVal)) {
		throw subpavings::MalconstructedBox_Error(
		"AdaptiveHistogram::makeBox(const RVecData&, const size_t)");
	}
	
    return retVal;

}


// insert data from a container
// return number of data points inserted into the actual subpaving
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
        outputLogPlain(s, i);
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
                getSubPaving()->insertOneFind(it,ON_PARENT, boolTest);

        //insertOneFind returns either NULL if no insert possible
        // or a pointer to the node the data goes to before that node
        // is split (it could be split more than once)
        if (NULL == insertedInto) { // failed to insert
		
			// erase it from the list
			dataCollection.erase(it);
			
            std::cerr << "Failed to insert point ";
			prettyPrint(std::cerr, *cit);
			std::cerr << std::endl;
        }
        // successful insertion, and we are splitting as we go
        else {
			counter++;
			
			if (boolTest() == true) {
				// if we split on insertion we may want to log
				if (!(insertedInto->isLeaf()) && logging) { // log the current state of the histogram
						outputLogPlain(s, i);
						i++;
				}
			}
		}
	} // end iteration through data

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
    scaledEMPSumCOPERR = getSubPaving()->getEMPSumCOPERR(
                        getRootCounter());
}

// Recalculate the unscaled EMP part of AIC score.
void AdaptiveHistogram::recalcScaledEMPSumAIC() const
{
    // use the scaled EMP Sum from the root node's getEMPSumCOPERR()
    scaledEMPSumAIC = getSubPaving()->getEMPSumAIC(
                        getRootCounter());
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
void AdaptiveHistogram::outputLog(const std::string& s, 
											int i, int prec) const
{
	// To add output of the AdaptiveHistogram object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        size_t n = getRootCounter();

        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
       
        getSubPaving()->leavesOutputTabsWithEMPs(n, os, prec); // the output
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

// Method to add current state of the histogram during splitting to a log file
// Output goes to file named according to argument s
// Output is plain, just textToTabs
void AdaptiveHistogram::outputLogPlain(const std::string& s, 
											int i, int prec) const
{
    // To add output of the AdaptiveHistogram object to file
    ofstream os(s.c_str(), ios::app);         // append
    if (os.is_open()) {
        os << std::endl;
        os << "Pass " << i << std::endl; // numbering
        getSubPaving()->leavesOutputTabs(os, prec); // the output
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
void AdaptiveHistogram::logMCMCDeltas(const std::string& s, int i,
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

//NEWAPRIL2012
// log changes in log of MCMC proposal probability compared to random (0,1)
// second version - logs more info about node
void AdaptiveHistogram::logMCMCDeltasAugmented(std::string s, int i,
                            int nodeType, int accepted,
							size_t nCount,
							real deltaL,
							real deltaP,
							real deltaQ,
							real deltaPi,
                            double randChange)
{
    RealVec vals;
	vals.push_back(deltaL);
    vals.push_back(deltaP);
    vals.push_back(deltaQ);
    vals.push_back(deltaPi+deltaQ);
    vals.push_back(log(randChange));
    outputFile(s, vals, i, nodeType, accepted, nCount);
}

// log for IMH MCMC proposal probability compared to random [0,1)
void AdaptiveHistogram::logMCMCIMH(std::string s, int i,
                            size_t currentLeaves,
							size_t proposedLeaves,
							real deltaL,
							real deltaP,
							real deltaQ,
							real logAcceptance,
                            double randChange,
							bool accepted)
{
	try {
		ostringstream oss;
		oss << "\n" << i << "\t" << currentLeaves << "\t" << proposedLeaves << "\t" 
				<< deltaL << "\t" << deltaP << "\t"
				<< deltaQ << "\t" << logAcceptance << "\t" 
				<< log(randChange) << "\t" << accepted;
		outputFileString(s, oss.str()); 
	}
	catch (...) {
		outputFileString(s, "(logging failure)"); 
		//swallow exception
	}
}

// start a log file for MCMC
// s is the file name
void AdaptiveHistogram::MCMCStartLogFile(const std::string& s, int i,
                const MCMCProposal& proposal, const LogMCMCPrior& logPrior)
{
    // Start log file with filename and timestamp
    outputLogStart(s);
    // put in the name of the proposal and prior
    std::string line = "Prior is " + logPrior.getName();
    line += ", proposal is " + proposal.getName();
    outputFile(s, line);
    // log the current state of the histogram
    outputLogPlain(s, i);
    
	#ifdef LOGEMPS
		outputLogEMPAIC(s); // output AIC score information
		 // output COPERR score information
		outputLogEMPCOPERR(s);
	#endif
   
    std::string headers 
			= "deltaL \t deltaP \t deltaQ \t deltaPi&Q \t ln(rand)";
    #if defined OLDLOGSTYLE
		headers += "\t ratioL \t ratioP \t ratioQ \t ratioPi&Q \t rand";
	#else
		headers = "\t\t\t\t" + headers;
	#endif
    outputFile(s, headers);

}

// start a log file for MCMC
// s is the file name
void AdaptiveHistogram::MCMCStartLogFile(const std::string& s,
                const LogMCMCIMHPrior& logPrior)
{
    // Start log file with filename and timestamp
    outputLogStart(s);
    // put in the name of the prior
    std::string line = "Prior is " + logPrior.getName();
    outputFile(s, line);
    
	std::string headers 
		= "l\tl'\tdeltaL\tdeltaP\tdeltaQ\tdeltaPi&Q\tln(rand)\taccept";
    headers = "\t" + headers;
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
void AdaptiveHistogram::MCMCLogFinalState(const std::string& s, int i)
{
    outputLogPlain(s, i);
	#ifdef LOGEMPS
		// output AIC score information
		outputLogEMPAIC(s);
		// output COPERR score information
		outputLogEMPCOPERR(s);
	#endif	
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


//internal method to get samples for Independent Metropolis Hastings MCMC
// does NOT normalise the samples
// average == true means we return just one pcf, the average
std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::_MCMCsamplesIMH(
						std::vector < PiecewiseConstantFunction >& samples, 
						bool average,
						size_t maxLeaves,
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, 
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	#ifdef DEBUG_IMHMCMC
		cout << "\nIn _MCMCsamplesIMH" << endl;
		cout << "maxLeaves " << maxLeaves << endl;
		
	#endif
	
    std::string errorMsg("AdaptiveHistogram::_MCMCsamplesIMH(...)");
	if (!hasSubPaving()) {

		throw NullSubpavingPointer_Error(errorMsg);
	}
	size_t bigN = getRootCounter();
	if (!bigN) {
		errorMsg += "\n: no data";
		throw std::runtime_error(errorMsg);
	}
	if (burnin > loops) {
		errorMsg += "\n: burnin > loops";
		throw std::invalid_argument(errorMsg);
	}
	if (thinout > 0 && ((loops-burnin+1)/thinout < 1)) {
		errorMsg += "\n: thinout value means no samples will be taken";
		throw std::invalid_argument(errorMsg);
    }
	
	//check the starting state
	if (!(getSubPaving()->checkTreeStateLegal(minPoints))) {
		errorMsg += "\n: illegal starting state";
		throw std::runtime_error(errorMsg);
	}
   
	if (logging == TXTANDGRAPH) {
		cerr << "Sorry, TXTANDGRAPH not supported, changing to TXT logging" << endl;
		logging = TXT;
	}
	if (logging == LOGANDGRAPHSAMPLES) {
		cerr << "Sorry, LOGANDGRAPHSAMPLES not supported, changing to LOGSAMPLES logging" << endl;
		logging = LOGSAMPLES;
	}
	// set up the partitioner 
    unsigned int seed = gsl_rng_get(rgsl);
	MCMCPartitionGenerator partitioner(seed);
	
	// set up the shadow paving
	RealMappedSPnode shadow(*getSubPaving());
	
	bool thinning = (thinout > 0);

	// for samples
	std::vector< PiecewiseConstantFunction > tmp;
	
	// an average
	PiecewiseConstantFunction avg(getRootBox(), getLabel() );

	if (thinning) {
		if (average) { // averaging
			tmp.push_back( avg );
		}
		else { // sampling 
		 
			tmp.reserve((loops - burnin)/thinout + 1);
		}
	}

	//Failure logging
	string failureLogFilename("FailureLogIMH.txt");
	outputFileStart(failureLogFilename);
	
	// for logging
	unsigned int i = 0;
	std::string sDec = ""; // for logging decisions
	std::string sPCF = ""; // for logging 'real' pcf state
	std::string dotPCF = "";
	std::string tracePCF = ""; // logging leaf and log posterior traces for pcf shadow
	std::string traceAvg = ""; // logging leaf and log posterior traces for pcf avg
	if ( (logging == TXT) ||(logging == LOGSAMPLES)) {

		std::string baseFileNamePCF = "PcfIMHmcmcOutput";
		
		cout << "Log files for decisions and pcf are " << endl;	
		if ( (logging == LOGSAMPLES)) {

			std::string dummySuffix("");
			sPCF = getUniqueFilename(baseFileNamePCF,dummySuffix);
		}
		else {
			sPCF = getUniqueFilename(baseFileNamePCF);
			// pass to log output to keep track of splits
			std::string baseFileNameDec = "DecIMHmcmcOutput";
			sDec = getUniqueFilename(baseFileNameDec);
			cout << "\t" << sDec << endl; 
		}
		
		
		cout << "\t" << sPCF << endl;
	}
	
	if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
		std::string baseFileNameStateTrace = "StateTraceIMH";
		tracePCF = getUniqueFilename(baseFileNameStateTrace);
		cout << "Log file for state leaf and log-posterior trace is " << endl;
		cout << "\t" << tracePCF << endl; 
		
		if (logging == LOGMCMCTRACE) {
			std::string baseFileNameAvgTrace = "AvgTraceIMH";
			traceAvg = getUniqueFilename(baseFileNameAvgTrace);
			cout << "Log file for average leaf and log-posterior trace is " << endl;
			cout << "\t" << traceAvg << endl; 
		}
	}
	
	/* get the starting point */
	size_t numLeaves = getRootLeaves();
	//unscaled log likelihood
	real unscaledLoglik = getSubPaving()->getUnscaledTreeLogLik();
	
	/* Logging for log posteriors:
	 * Get the log likelihood for our starting histogram */
	real vol = getSubPaving()->nodeRealVolume();
	
	/* scaling is -n x ln (n x vol) where n is 
	the total points (bigN) and vol is the volume of the root box */
	real loglik = unscaledLoglik - 1.0*bigN * cxsc::ln(vol * (1.0*bigN));

	/* and use prior to get posterior:
	We will update the log posterior as we go */
	cxsc::real logPosterior = loglik + logPrior(numLeaves - 1);
	
	// stuff we need to track the log posterior of the sample average
	
	cxsc::real logPosteriorSampleAverage(0.0);
	size_t numberOfSamples = 0;
		
	
	#ifdef DEBUG_IMHMCMC
		cout << "\nAt start" << endl;
		cout << "numLeaves " << numLeaves << endl;
		cout << "unscaledLoglik " << unscaledLoglik << endl;
		cout << "Loglik " << loglik << endl;
		cout << "logPrior(numLeaves - 1) " << (logPrior(numLeaves - 1)) << endl;
		cout << "logPosterior " << logPosterior << endl;
		
		
	#endif
	
	
	bool goodLoop = true;

	if (logging == TXT) {
		
		MCMCStartLogFile(sDec, logPrior);
		// and record the shadow state
		outputLogStart(sPCF);
		shadow.outputLog(sPCF, i);
	}
	
	//for logposteriors
	if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
		outputLogStart(tracePCF);
		ostringstream os;
		os << i << "\t" << numLeaves << "\t" << logPosterior << endl;
		outputFileString(tracePCF, os.str());

		if (logging == LOGMCMCTRACE) {
			outputLogStart(traceAvg);
		}
	}

	++i;

	std::string stateNow = "";
	std::string stateAfter = "";
	
	//holder for the partition insructions from the partitioner
	std:: vector < unsigned long int > instructions;
				
	// loop from here conditional on good loop and cancontinue
	while (goodLoop && (loops > 0) ) {

		loops--;

		try {
			/* changeMCMCStateIMH updates numLeaves, unscaledLoglik, logPosterior.
			 * The shadow is changed exactly in line with this. */
			bool success = changeMCMCStateIMH(
						maxLeaves,
                        partitioner,
						logPrior,
                        minPoints,
						shadow,
						numLeaves,
						unscaledLoglik,
						logPosterior,
						rgsl, 
						logging,
                        sDec,
						i,
						failureLogFilename);
			
			// logPosterior etc updated (if success)
			
			#ifdef DEBUG_IMHMCMC
			{
				cout << "After changeMCMCStateIMH," << endl;
				cout << "numLeaves = " << numLeaves << endl;
				cout << "unscaledLoglik = " << unscaledLoglik << endl;
				cout << "logPosterior = " << logPosterior << endl;
			}
			#endif
			
			//get the latest instructions only if change was successful 
			if (success) instructions = partitioner.getInstructions();
		}
		catch (UnfulfillableRequest_Error&  ure) {
			goodLoop = false;
			std::cerr << "Change in mcmc state failed: aborting process" << endl;
		}

		if (goodLoop) {
			
			#ifdef MCMCTRACK
				if (i%10000 == 0) cout << i << endl;
			#endif
			
			if (logging == TXT) {
				/* log the current state of the pcf
				 * 
				 * Note that this is relatively expensive because
				 * we have to make a new temporary object each time */
				
				RealMappedSPnode shadowCopy(shadow);
				shadowCopy.reshapeToMirrorValidSplitPartitions(
							instructions, numLeaves);
				
				shadowCopy.outputLog(sPCF, i);
			}
			if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
				ostringstream os;
				os << i << "\t" << numLeaves << "\t" << logPosterior  << endl;
				outputFileString(tracePCF, os.str());
			}

			// if we are taking samples take the sample here
			if (thinning && (i >= burnin) &&
					((i-burnin)%thinout == 0)) {
						
				/* If we want an actual sample we have to 
				 * take a copy of the shadow, then reshape it according
				 * to the instructions*/
				
				RealMappedSPnode shadowCopy(shadow);
				shadowCopy.reshapeToMirrorValidSplitPartitions(
							instructions, numLeaves);
				
				#ifdef DEBUG_IMHMCMC
				{
					cout << "\nSampling, shadowCopy after reshape has "
						<< shadowCopy.getNumberLeaves() << " leaves" << endl;
				}
				#endif

				/* If we are after the average
				 * or if we want average logs, 
				 * add a PCF with subpaving from shadow and label
				 * from this to avg- it is NOT normalised here
				 * and we have not divided by the number of samples*/
				if (average || (logging == LOGMCMCTRACE)) {
					avg += PiecewiseConstantFunction(shadowCopy, getLabel() );
					#ifdef DEBUG_IMHMCMC
						{
							cout << "Averaging, average now has  "
							<< avg.getRootLeaves() << " leaves" << endl;
						}
					#endif
				}
				++numberOfSamples;
				
				/* If we want samples rather than just the average, we also
				 * make a PCF with subpaving from shadow and label
				 * from this - it is NOT normalised here - to store as the sample*/
				if (!average) {
					
					tmp.push_back( 
						PiecewiseConstantFunction(shadowCopy, getLabel()) );
				}
				
				/* Because calculating the log posterior for the average is 
				 * relatively expensive, I'm only doing it if we want to log it */
				if (logging == LOGMCMCTRACE) {
					
					/* update the logPosteriorSampleAverage
					 * For this we need the likelihood and the prior for
					 * the average.  The log likelihood we get from the
					 * pieces in the average and the counts in this,
					 * but the 'average' is in fact not normalised yet
					 * so we need to deduct bigN(log(#samples) + log(bigN)).*/
					cxsc::real loglikSampleAvgUnnormalised = avg.getLogLikelihood(*this);
					cxsc::real loglikSampleAvg = loglikSampleAvgUnnormalised - 
						(1.0*bigN * (cxsc::ln(1.0*numberOfSamples) + cxsc::ln(1.0*bigN)) );
					/*For the logPrior for the average we need the number of
					leaves in the average.  This is a real pain because there
					is no easy way to keep track of it:  it is not necessarily the 
					most-leaved sample we have ever taken because two different 
					samples split in different ways will add to give a result that has
					more leaves than either of them.  Adding up the number of leaves
					in a big tree takes time ...*/
					size_t sampleAvgLeaves = avg.getRootLeaves();
					
					cxsc::real logPriorSampleAvg = logPrior(sampleAvgLeaves-1);
					
					logPosteriorSampleAverage = loglikSampleAvg + logPriorSampleAvg;
					
					#ifdef DEBUGLOGPOSTERIORAVG
					{
						cout << "\nsampleAvgLeaves = " << sampleAvgLeaves << endl;
						cout << "bigN = " << bigN << endl;
						cout << "numberOfSamples = " << numberOfSamples << endl;
						cout << "loglikSampleAvgUnnormalised = " 
									<< loglikSampleAvgUnnormalised << endl;
						cout << "loglikSampleAvg = " << loglikSampleAvg << endl;
						cout << "logPriorSampleAvg = " << logPriorSampleAvg << endl;
						cout << "logPosteriorSampleAverage = " << logPosteriorSampleAverage << endl;
					}
					#endif
					// output loglik and prior and logposterior separately
					ostringstream os;
					os << numberOfSamples << "\t" << sampleAvgLeaves 
							<< "\t" << loglikSampleAvg
							<< "\t" << logPriorSampleAvg
							<< "\t" << logPosteriorSampleAverage  << endl;
					outputFileString(traceAvg, os.str());
				}
				
				
			
				if (logging == LOGSAMPLES) {
					
					/* want a new log file for each shadow sampled */
					{
						ostringstream oss;
						oss << sPCF << "_" << i << ".txt";
						shadowCopy.outputLog(oss.str(), i);
					}
				}
				
			}
		}

		++i;
		// back into loop
	}
	
	//UPDATED JUNE 2012 for logposteriors
	if (goodLoop) { //no change in samples if loop failed
		//if we are averaging, swap the average into the tmp container
		if (average) {
			
			avg.swap(tmp.back());			
			
		}
		tmp.swap(samples);
	}
	
	return samples;
}


//UPDATED JUNE 2012 for logposteriors
//internal method to get MCMCsamples
// does NOT normalise the samples
// average == true means we return just one pcf, the average
std::vector < PiecewiseConstantFunction >& AdaptiveHistogram::_MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						bool average,
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, 
						bool volChecking,
						real minVol,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl)
{
	
    std::string errorMsg("AdaptiveHistogram::_MCMCsamples(");
			errorMsg += "std::vector < PiecewiseConstantFunction >&, bool, ";
			errorMsg += "unsigned int, ";
			errorMsg += "unsigned int, unsigned int, MCMCProposal&, ";
			errorMsg+= "LogMCMCPrior&, size_t, LOGGING_LEVEL, gsl_rng *)";
	
	if (!hasSubPaving()) {

		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	if (burnin > loops) {
		errorMsg += "\n: burnin > loops";
		throw std::invalid_argument(errorMsg);
	}
	
	
	
	bool thinning = (thinout > 0);

	std::vector< PiecewiseConstantFunction > tmp;
	
	//UPDATED JUNE 2012 for logposteriors
	PiecewiseConstantFunction avg(getRootBox(), getLabel() );
	
	//UPDATED JUNE 2012 for logposteriors
	if (thinning) {
		if (average) { // averaging
			tmp.push_back( avg );
		}
		else { // sampling 
		 
			tmp.reserve((loops - burnin)/thinout + 1);
		}
	}

	// for logging
	unsigned int i = 0;
	std::string sHist = "";
	std::string sPCF = "";
	std::string dotPCF = "";
	//UPDATED JUNE 2012 for logposteriors
	std::string tracePCF = ""; // logging leaf and log posterior traces for pcf shadow
	std::string traceAvg = ""; // logging leaf and log posterior traces for pcf avg
	if ( (logging == TXT) || (logging == TXTANDGRAPH)
				||(logging == LOGSAMPLES) 
				|| (logging == LOGANDGRAPHSAMPLES)) {

		// pass to log output to keep track of splits
		std::string baseFileNameHist = "histMCMCOutput";
		sHist = getUniqueFilename(baseFileNameHist);
		
		std::string baseFileNamePCF = "pcfMCMCOutput";
			
		if ( (logging == LOGSAMPLES) 
			|| (logging == LOGANDGRAPHSAMPLES) ) {

			std::string dummySuffix("");
			sPCF = getUniqueFilename(baseFileNamePCF,dummySuffix);
		}
		else {
			sPCF = getUniqueFilename(baseFileNamePCF);
		}
		cout << "Log files for hist and pcf are " << endl;
		cout << "\t" << sHist << endl; 
		if ( (logging == LOGSAMPLES) 
			|| (logging == LOGANDGRAPHSAMPLES) ) {
			cout << "\t" << sPCF << endl;
		}
		else {
			#ifndef DISABLE_MCMC_PCF_LOGGING
				cout << "\t" << sPCF << endl;
			#endif
		}

		if (logging == TXTANDGRAPH) {

			// for dot graph
			std::string baseFileNameGraph = "graph";
			std::string suffix = ".dot";
			dotPCF = getUniqueFilename(baseFileNameGraph, suffix);
			outputFile(dotPCF, "digraph G {"); // opening line
		}
	}
	

	//UPDATED JUNE 2012 for logposteriors
	if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
		std::string baseFileNameStateTrace = "stateTrace";
		tracePCF = getUniqueFilename(baseFileNameStateTrace);
		cout << "Log file for state leaf and log-posterior trace is " << endl;
		cout << "\t" << tracePCF << endl; 
		
		if (logging == LOGMCMCTRACE) {
			std::string baseFileNameAvgTrace = "avgTrace";
			traceAvg = getUniqueFilename(baseFileNameAvgTrace);
			cout << "Log file for average leaf and log-posterior trace is " << endl;
			cout << "\t" << traceAvg << endl; 
		}
	}

	// set up a container for the leaf children
	SPSnodePtrs leafVec;
	// set up a container for the subleaf children
	SPSnodePtrs cherryVec;
	SPSnodeList nodes;
	 //lists better than vectors for random access removal
	size_t numLeaves = 0;
	size_t numCherries = 0;

	// fill the container with the leaf children
	getSubPaving()->getLeaves(leafVec);

	// fill the container with the subleaf children
	getSubPaving()->getSubLeaves(cherryVec);
	
	//check the cherries are all "legal", ie pass isSplittableNode and minVol
	if (!cherryVec.empty()) {
		SPSnodePtrsItr cit;
		for (cit = cherryVec.begin(); cit < cherryVec.end(); ++cit) {
			if ((volChecking && ((*cit)->nodeRealVolume() < 2*minVol))
				||
				!((*cit)->isSplittableNode(minPoints)))
			{
				throw std::logic_error(errorMsg 
				+ string(
				"\nIllegal state - cherries do not satisfy minPoints and minVol for split"));
			}
		}
	}

	numCherries = cherryVec.size();

	if (!leafVec.empty()) {
		/* but only put into the container the leaves which, if split,
		 * would have at least minPoints data points associated with them
		 * and which can be split to give children with vol >= minVol */

		SPSnodePtrsItr lit;
		for (lit = leafVec.begin(); lit < leafVec.end(); ++lit) {
			if  ((!volChecking || 
			
				!((((*lit)->getParent() == NULL) && ((*lit)->nodeRealVolume() < 2*minVol))
					|| 
				(*lit)->getParent()->nodeRealVolume() < 4*minVol))
				&&
				((*lit)->isSplittableNode(minPoints)) )
			{
				// leaf can go into container
				nodes.push_back(*lit);
				++numLeaves;
			}
		}
	}

	// no need to check on cherries now - they can all go in
	if (numCherries > 0)
		nodes.insert(nodes.end(), cherryVec.begin(),cherryVec.end());

	if (nodes.size() == 0) {
		throw UnfulfillableRequest_Error(errorMsg);
	}
	
	/* Make a RMP out of the current subpaving and run that
	 * alongside changes to this; when sampling this is 
	 * copied into a PCF ...*/
	
	RealMappedSPnode shadow(*getSubPaving());
	RealMappedSPnode::Ptrs shadowLeavesVec;
	RealMappedSPnode::Ptrs shadowCherriesVec;
	shadow.getLeaves(shadowLeavesVec);
	shadow.getSubLeaves(shadowCherriesVec);
	RealMappedSPnode::ListPtrs 
			shadowLeaves(shadowLeavesVec.begin(),shadowLeavesVec.end());
	RealMappedSPnode::ListPtrs 
			shadowCherries(shadowCherriesVec.begin(),shadowCherriesVec.end());
	
	//UPDATED JUNE 2012 for logposteriors		
	/* Get the log likelihood for our starting histogram */
	cxsc::real loglik = getLogLikelihood();
	/* and use prior to get posterior:  leafVec has all our actual leaves in
	so we can use that to tell us the number of leaves 
	We will update the log posterior as we go */
	cxsc::real logPosterior = loglik + logPrior(leafVec.size() - 1);
	
	//UPDATED JUNE 2012 for logposteriors
	// stuff we need to track the log posterior of the sample average
	size_t bigN = getRootCounter();
	cxsc::real logPosteriorSampleAverage(0.0);
	size_t numberOfSamples = 0;
		
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nAt start" << endl;
		cout << "shadowCherries.size() " << (shadowCherries.size()) << endl;
		cout << "shadowLeaves.size() " << (shadowLeaves.size()) << endl;
		
	#endif
	
	//UPDATED JUNE 2012 for logposteriors
	#ifdef DEBUGLOGPOSTERIOR
		cout << "leafVec.size() = " << leafVec.size() << endl;
		cout << "log-posterior = " << logPosterior << endl;
	#endif
	
	bool goodLoop = true;

	if ((logging == TXT) || (logging == TXTANDGRAPH) ) {

		MCMCStartLogFile(sHist, i, proposal, logPrior);
		// and record the shadow state
		#ifndef DISABLE_MCMC_PCF_LOGGING
			outputLogStart(sPCF);
			shadow.outputLog(sPCF, i);
		#endif
		
	}
	
	//UPDATED JUNE 2012 for logposteriors
	if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
		outputLogStart(tracePCF);
		ostringstream os;
		os << i << "\t" << shadowLeaves.size() << "\t" << logPosterior << endl;
		outputFileString(tracePCF, os.str());

		if (logging == LOGMCMCTRACE) {
			outputLogStart(traceAvg);
		}
	}

	++i;

	std::string stateNow = "";
	std::string stateAfter = "";
	
	// info will keep track of changes in state for us
	ChangeOfStateInformationBasic info(shadowLeaves.size(),
										shadowCherries.size());

	// loop from here conditional on good loop and cancontinue
	while (goodLoop && (loops > 0) ) {

		if (logging == TXTANDGRAPH) {
			// capture state now
			stateNow = "\"" + (shadow.getLeafNodeLevelsString()) + "\"";
		}
		
		loops--;

		try {
			//UPDATED JUNE 2012 for logposteriors - change now returns deltaPI
			/* changeMCMCState updates nodes, numLeaves, numCherries, i
			shadow, shadowLeaves, shadowCherries, and --
			UPDATED OCT 2102 -- crucial 
			change in log posterior */
			changeMCMCState(nodes, 
						numLeaves, numCherries,
						proposal, logPrior, minPoints, 
						minVol,
						shadowLeaves, shadowCherries,
						info, // updated
						rgsl, logging, sHist, i);
			
			//UPDATED JUNE 2012 for logposteriors
			// update the log posterior
			logPosterior += info.getDeltaPi();
			
			assert(info.getCurrentLeaves() == shadowLeaves.size());
			assert(info.getCurrentCherries() == shadowCherries.size());
			
			//UPDATED JUNE 2012 for logposteriors
			#ifdef DEBUGLOGPOSTERIOR
			{
				cout << "deltaPi = " << info.getDeltaPi() << endl;
				cout << "logPosterior = " << logPosterior << endl;
			}
			#endif
			
		}
		catch (UnfulfillableRequest_Error&  ure) {
			goodLoop = false;
			std::cerr << "Change in mcmc state failed: aborting process" << endl;
		}

		if (goodLoop) {
			
			#ifdef MCMCTRACK
				if (i%10000 == 0) cout << i << endl;
			#endif
			
			if ((logging == TXT) || (logging == TXTANDGRAPH)) {
				// log the current state of the pcf
				#ifndef DISABLE_MCMC_PCF_LOGGING
					shadow.outputLog(sPCF, i);
				#endif
				#ifdef LOGEMPS
					outputLogEMPAIC(s);
				#endif
			}
			//UPDATED JUNE 2012 for logposteriors
			// UPDATED OCT 2012 to use leaves from info
			if (logging == LOGSTATETRACE || logging == LOGMCMCTRACE) {
				ostringstream os;
				os << i << "\t" << info.getCurrentLeaves() << "\t" << logPosterior  << endl;
				outputFileString(tracePCF, os.str());
			}

			if (i >= burnin && (logging == TXTANDGRAPH)) {
				// capture state after split or merge to graph file

				stateAfter = "\"" + shadow.getLeafNodeLevelsString() + "\"";
				std::string line = "\t " + stateNow + " -> " + stateAfter + ";";
				outputFile(dotPCF, line);

				outputGraphDot(); // and make a graph of current state
				
			}

			if ((numLeaves == 0 && numCherries == 0)) {
				errorMsg += "\n: No more leaves or cherries in MCMC";
				throw std::logic_error(errorMsg);

			}

			// if we are taking samples take the sample here
			if (thinning && (i >= burnin) &&
					((i-burnin)%thinout == 0)) {

				//UPDATED JUNE 2012 for logposteriors
				/* If we are after the average
				 * or if we want average logs, 
				 * add a PCF with subpaving from shadow and label
				 * from this to avg- it is NOT normalised here
				 * and we have not divided by the number of samples*/
				if (average || (logging == LOGMCMCTRACE)) {
					avg += PiecewiseConstantFunction(shadow, getLabel() );
				}
				++numberOfSamples;
				
				//UPDATED JUNE 2012 for logposteriors
				/* If we want samples rather than just the average, we also
				 * make a PCF with subpaving from shadow and label
				 * from this - it is NOT normalised here - to store as the sample*/
				if (!average) {
					
					tmp.push_back( 
						PiecewiseConstantFunction(shadow, getLabel()) );
				}
				
				//UPDATED JUNE 2012 for logposteriors - right down to end of degugging output section
				/* Because calculating the log posterior for the average is 
				 * relatively expensive, I'm only doing it if we want to log it */
				if (logging == LOGMCMCTRACE) {
					
					/* update the logPosteriorSampleAverage
					 * For this we need the likelihood and the prior for
					 * the average.  The log likelihood we get from the
					 * pieces in the average and the counts in this,
					 * but the 'average' is in fact not normalised yet
					 * so we need to deduct bigN(log(#samples) + log(bigN)).*/
					cxsc::real loglikSampleAvg = avg.getLogLikelihood(*this);
					loglikSampleAvg -= 
						(1.0*bigN * (cxsc::ln(1.0*numberOfSamples) + cxsc::ln(1.0*bigN)) );
					/*For the logPrior for the average we need the number of
					leaves in the average.  This is a real pain because there
					is no easy way to keep track of it:  it is not necessarily the 
					most-leaved sample we have ever taken because two different 
					samples split in different ways will add to give a result that has
					more leaves than either of them.  Adding up the number of leaves
					in a big tree takes time ...*/
					size_t sampleAvgLeaves = avg.getRootLeaves();
					
					cxsc::real logPriorSampleAvg = logPrior(sampleAvgLeaves-1);
					
					logPosteriorSampleAverage = loglikSampleAvg + logPriorSampleAvg;
					
					#ifdef DEBUGLOGPOSTERIORAVG
					{
						cout << "\n\nsampleAvgLeaves = " << sampleAvgLeaves << endl;
						cout << "bigN = " << bigN << endl;
						cout << "numberOfSamples = " << numberOfSamples << endl;
						cout << "loglikSampleAvg = " << loglikSampleAvg << endl;
						cout << "logPriorSampleAvg = " << logPriorSampleAvg << endl;
						cout << "logPosteriorSampleAverage = " << logPosteriorSampleAverage << endl;
					}
					#endif
					//UPDATED 4JULY2012 to output loglik and prior and logposterior separately
					ostringstream os;
					os << numberOfSamples << "\t" << sampleAvgLeaves 
							<< "\t" << loglikSampleAvg
							<< "\t" << logPriorSampleAvg
							<< "\t" << logPosteriorSampleAverage  << endl;
					outputFileString(traceAvg, os.str());
				}
				
				
			
				if (logging == LOGSAMPLES || logging == LOGANDGRAPHSAMPLES) {
					// log this sample to log file
					outputLogPlain(sHist, i);
					
					/* want a new log file for each shadow sampled */
					{
						ostringstream oss;
						oss << sPCF << "_" << i << ".txt";
						shadow.outputLog(oss.str(), i);
					}
					
					#ifdef LOGEMPS
						outputLogEMPAIC(sHist); // add AIC scores
					#endif
				}
				if (logging == LOGANDGRAPHSAMPLES) outputGraphDot();
			}
		}

		++i;
		// back into loop
	}
	
	if (logging == TXTANDGRAPH) {
		// close the dot graph of the changes process
		outputFile(dotPCF, "}"); // closing line

		// make the graph image
		makeDotImage(dotPCF);
	}

	//UPDATED JUNE 2012 for logposteriors
	if (goodLoop) { //no change in samples if loop failed
		//if we are averaging, swap the average into the tmp container
		if (average) {
			
			avg.swap(tmp.back());			
			
		}
		tmp.swap(samples);
	}
	
	return samples;
}

cxsc::real AdaptiveHistogram::changeMCMCState(SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i)
{
	/*horrible - give it a dummy number of leaves and cherries
	because we don't actually use that! */
	ChangeOfStateInformationBasic info(1, 0);
	real minVol(0.0);
						
	
	changeMCMCState (nodes,
                        numLeaves, numCherries,
                        proposal, logPrior,
                        minPoints,
						minVol,
						leaves,
						cherries,
						info,
                        rgsl, logging,
                        sHist, i);
	
	return info.getDeltaPi();
}

cxsc::real AdaptiveHistogram::changeMCMCState(SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						real minVol,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i)
{
	/*horrible - give it a dummy number of leaves and cherries
	because we don't actually use that! */
	ChangeOfStateInformationBasic info(1, 0);
	
	changeMCMCState (nodes,
                        numLeaves, numCherries,
                        proposal, logPrior,
                        minPoints,
						minVol,
						leaves,
						cherries,
						info,
                        rgsl, logging,
                        sHist, i);
	
	return info.getDeltaPi();
}


// Returns an interator to node to propose for changes
// alters haveNode to true if a proposal node has been found (otherwise the
// iterator just points to the beginning of nodes and should not be used).
SPSnodeListItr AdaptiveHistogram::proposeChangeMCMCState
                        (const MCMCProposal& proposal,
                        SPSnodeList& nodes,
                        size_t numLeaves, size_t numCherries,
                        gsl_rng* rgsl, bool& haveNode)
{
    SPSnodeListItr it;

    RealVec probs;
    RealVecItr pit;

	// fillNodeProposalProbs returns the sum of the probabilities and
    // also fills in probs
    // the sum of the probabilities may be < 1
    real psum = proposal.fillNodeProposalProbs(numLeaves, numCherries, probs);

	if (probs.empty()) throw std::runtime_error(
		"AdaptiveHistogram::proposeChangeMCMCState(...) : probs empty");
    else {

        // pick a node at random  by drawing a random number in [0,1)
        double rand = gsl_rng_uniform(rgsl);
		
		// we'll only pick a node if rand < psum
		// if we have at least some leaves
		if ((numLeaves > 0) && rand < psum) {
			pit = probs.begin();
			it = nodes.begin();
			real sum = 0.0;

			for (pit = probs.begin(); pit < probs.end(); ++pit) {
					sum += *pit;
					if (rand < sum) {
						haveNode = true;
						break;
					}
					++it; // using this to point to nodes relies on having got
							// the right number of probabilities to match nodes
			}
		}
		
		// if we only have cherries (possible if no leaves are splittable)
		/* then we'll only pick a cherry if rand >= 1-psum
		 * I've now realised that this is unnecessarily complex but 
		 * it still works even with the stay-split-merge proposal because
		 * the probability of rand >= 1-psum = is still psum
		 * so keep it so that we still get the same results as before*/
		if ((numLeaves == 0) && (numCherries > 0) && (rand >= 1.0-psum)) {
			pit = probs.begin();
			it = nodes.begin();
			real sum = 1.0-psum; // note that sum starts at 1-psum

			for (pit = probs.begin(); pit < probs.end(); ++pit) {
					sum += *pit;
					if (rand < sum) {
						haveNode = true;
						break;
					}
					++it; // using this to point to nodes relies on having got
							// the right number of probabilities to match nodes
			}
		}
		// it iterator should now point to the node we want to target
		//else we've not taken a node

    }
	
	
	return it;
}



//UPDATED JUNE 2012 for logposteriors - whole function 
// returns true/false decision on whether to split or not
bool AdaptiveHistogram::decisionMCMCSplitNEW(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        const size_t numLeaves, const size_t numCherries,
						SPSnodeList& nodes, 
						size_t realNumLeaves,
						size_t minPoints,
						bool volChecking,
						real minVol,
                        int& addLeftChild, int& addRightChild,
						int& deductParent,
						size_t& leftChildCount, size_t& rightChildCount,
						cxsc::real& deltaPi,
						LOGGING_LEVEL logging, 
						const std::string& s, int i) const
{
    bool willSplit = false;
	
	addLeftChild = 0;
	addRightChild = 0;
	deductParent = 0;

    #if defined OLDLOGSTYLE
    if ((logging == TXT) || (logging == TXTANDGRAPH)) {
        if (target->isLeaf()) 
			outputFile(s, "grabbing leaf " + target->getNodeName());
		else 
			outputFile(s, "grabbing virtual leaf " + target->getNodeName());
	}
	#endif
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn decisionMCMCSplitNEW, numLeaves = " << numLeaves << " numCherries = " << numCherries << endl;
		cout << "realNumLeaves = " << realNumLeaves << endl;
		if (target->isLeaf()) 
			cout << "grabbing leaf " <<( target->getNodeName() ) << endl;
		else 
			cout << "grabbing virtual leaf " <<( target->getNodeName() ) << endl;
			
		cout << "VolChecking = " << volChecking
			<< " and minVol = " << _double(minVol) 
			<< " (2 * minVol = " << _double(2*minVol)
			<< " and 4 * minVol = " << _double(4*minVol) << ")" << endl;
			
		cout << "This target node's volume is = " << _double(target->nodeRealVolume()) << endl;
	
	#endif
	
    //childrensSpread will be a container of the number of points the children
    // of each child of target might have, in order
    // [0] = left child's left child count, [1] = left child's rght child count,
    // [2] = rght child's left child count, [3] = rght child's rght child count,
    Size_tVec childrensSpread;
    childrensSpread =
                target->getChildrensLeftAndRightCountsIfSplitNEW(childrensSpread);

	// update references passed in
	leftChildCount = childrensSpread[0] + childrensSpread[1];
	rightChildCount = childrensSpread[2] + childrensSpread[3];

	#ifdef DEBUG_NEWMCMC
		cout << "childrensSpread[0] = " << childrensSpread[0] << " childrensSpread[1] = " << childrensSpread[1] << endl;
		cout << "childrensSpread[2] = " << childrensSpread[2] << " childrensSpread[3] = " << childrensSpread[3] << endl;
		cout << "using grandchild counts, leftChildCount = " << leftChildCount;
		cout << " rightChildCount = " << rightChildCount << endl;
	
	#endif

	// change in contribution to log likelihood for this node on split
	real deltaL = getSplitChangeLogLik(leftChildCount, rightChildCount);

    // realNumLeaves has to be supplied

    // use the prior distribution object to find the change in prior
    real deltaP = logPrior.changeOnSplitOne(realNumLeaves - 1);

    // posterior is proportional to likelihood * prior
    // deltaPi is passed in by reference; here we update it
    deltaPi = deltaL + deltaP;

    /* new numbers of leaves and cherries under proposal depends on minPoints
		and minVol
     * because these determine whether the new leaf children will go into the
     * nodes container */
    size_t newNumLeaves = numLeaves - 1; // current number of leaves less this
	
    /* increase the number of new leaves for each new child that can
     * go into the nodes container
	 * -- if we are vol checking then neither can go in if they don't 
	 * have vol >= 2*minVol, ie if this does not have vol >= 4*minVol */
    if ( !volChecking || !(target->nodeRealVolume() < 4*minVol) ) {
		#ifdef DEBUG_NEWMCMC
			cout << "Node volume is okay children to be splittable" << endl;
		#endif
		if (((childrensSpread[0] >= minPoints) && (childrensSpread[1] >= minPoints))
			// this will add one to the leaf numbers if minPoints == 0
			||
			// we will also be prepared to put the right child into the container if
			// there is a minPoints > 0 but one of its
			// children would take all the points, the other getting none
			((minPoints > 0) && (leftChildCount >= minPoints) &&
			((childrensSpread[0] == 0) || (childrensSpread[1] == 0))) ) {
		 
			++newNumLeaves;
			addLeftChild = 1; // change value referenced
		} 

		if (((childrensSpread[2] >= minPoints) && (childrensSpread[3] >= minPoints))
			// this will add one to the leaf numbers if minPoints == 0
			||
			// we will also be prepared to put the right child into the container if
			// there is a minPoints > 0 but one of its
			// children would take all the points, the other getting none
			((minPoints > 0) && (rightChildCount >= minPoints) &&
			((childrensSpread[2] == 0) || (childrensSpread[3] == 0))) ) {
		 
			++newNumLeaves;
			addRightChild = 1; // change value referenced
		}
	}
	#ifdef DEBUG_NEWMCMC
		cout << "deducting 1 from leaves for split,";
		if (addLeftChild) cout << " addLeftChild ";
		if (addRightChild) cout << " addRightChild ";
		cout << " newNumLeaves = " << newNumLeaves << endl;
		cout << "adding 1 to cherries for new cherry" << endl;
		
	#endif 
	
    size_t newNumCherries = numCherries + 1;
	
	/*But if the parent is in the cherries it would have to come out.*/
	SPSnode* nodeParent = target->getParent();
	
	if (nodeParent != NULL) {
		
		SPSnodeListItr git = nodes.begin();
		advance(git, numLeaves); // advance to the cherries
		// break out of loop if we find parent
		SPSnodeListItr it;
		for (it = git ; it != nodes.end(); ++it ) {
			if ((*it) == nodeParent) {
				newNumCherries--;
				deductParent = 1;
				break;
			}
		}
	
		#ifdef DEBUG_NEWMCMC
			cout << "node's parent is = " << (nodeParent->getNodeName()) << endl;
			if (deductParent) cout << "And this is in current cherries, so newNumCherries = " << newNumCherries << endl;
		#endif 
	}
	
    
    // Using proposal distribution object
    real deltaQ = proposal.getLogQRatioSplitProposal(numLeaves, numCherries,
                                                newNumLeaves, newNumCherries);
    //get another random number
    double randChange = gsl_rng_uniform(rgsl);
    if (log(randChange) < deltaPi + deltaQ) willSplit = true;

    if ((logging == TXT) || (logging == TXTANDGRAPH)) { // log these values
	
		#if defined OLDLOGSTYLE
			logMCMCDeltas(s, i, deltaL, deltaP, deltaQ, deltaPi, randChange);
			if (willSplit) outputFile(s, "Splitting");
			else outputFile(s, "Not splitting");
		#else // log with 10 for leaf as nodeType, or 11 if a virtual leaf
			int nodeType = 10;
			if (!target->isLeaf() ) nodeType = 11;
			logMCMCDeltasAugmented(s, i, nodeType, (willSplit ? 1: 0), 
					leftChildCount + rightChildCount,
					deltaL, deltaP, deltaQ, deltaPi, randChange);
		#endif
		
    }
	#ifdef DEBUG_NEWMCMC
		cout << "willSplit = " << willSplit << endl;
	#endif 
	
    return willSplit;
}


//UPDATED JUNE 2012 for logposteriors
// returns true false decision on whether to merge or not
bool AdaptiveHistogram::decisionMCMCMergeNEW(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        size_t numLeaves, size_t numCherries, 
						RealMappedSPnode::ListPtrs& leaves, 
						size_t minPoints,
                        bool volChecking,
						real minVol,
                        int& deductLeftChild, int& deductRightChild,
						int& addParent,
						cxsc::real& deltaPi,
						LOGGING_LEVEL logging, 
						const std::string& s, int i) const
{
    bool willMerge = false;

	deductLeftChild = 0;
	deductRightChild = 0;
	addParent = 0;

    size_t realNumLeaves = leaves.size();
	
	// cherry so we are merging
    #ifdef OLDLOGSTYLE
    if ((logging == TXT) || (logging == TXTANDGRAPH)) {
        if (target->isSubLeaf() )
			outputFile(s,"grabbing cherry " + target->getNodeName());
		
		else
			outputFile(s,"grabbing virtual cherry " + target->getNodeName());
	}
	#endif
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn decisionMCMCMergeNEW, numLeaves = " << numLeaves << " numCherries = " << numCherries << endl;
		cout << "realNumLeaves = " << realNumLeaves << endl;
		if (target->isSubLeaf()) 
			cout << "grabbing cherry " <<( target->getNodeName() ) << endl;
		else 
			cout << "grabbing virtual cherry " <<( target->getNodeName() ) << endl;
	
		cout << "VolChecking = " << volChecking
			<< " and minVol = " << _double(minVol) 
			<< " (2 * minVol = " << _double(2*minVol)
			<< " and 4 * minVol = " << _double(4*minVol) << ")" << endl;
			
		cout << "This target node's volume is = " << _double(target->nodeRealVolume()) << endl;
	
	
	#endif
	
   
    // change in log likelihood on merge is getMergeChangeLogLik
    // for this node
	real deltaL = rnd(target->getMergeChangeLogLik());

    // realNumLeaves must be given
	real deltaP = logPrior.changeOnMergeOne(realNumLeaves - 1);

    // posterior is proportional to likelihood * prior
	// deltaPi is passed in by reference; here we update it
    deltaPi = deltaL + deltaP;

    /* calculate the number of leaves and cherries after proposed merge
     * we have to take into account minPoints and the effect that this will have
     * had on whether the target's children are in the nodes container */
    
	size_t newNumLeaves = numLeaves + 1; // current number of leaves plus target
    
	/* but decrement newNumLeaves for each of the target's children that comes
	 * out of the container
     * -- if we are vol checking then children can only be in if they
	 * have vol >= 2*minVol, ie if this has vol >= 4*minVol */
    if ( !volChecking || !(target->nodeRealVolume() < 4*minVol) ) {
		#ifdef DEBUG_NEWMCMC
			cout << "Node volume is okay children to be splittable" << endl;
		#endif
		if ( (target->getLeftChild()->getMinChildCountIfSplitNEW() >= minPoints)
			||
			// the left child would also have been in the container if it had enough
			// points and all of them went to one child, the other getting nothing
			((minPoints > 0)
			&& (target->getLeftChild()->getCounter() >= minPoints)
			&& (target->getLeftChild()->getMinChildCountIfSplitNEW() == 0)) ){
				
				newNumLeaves--;
				deductLeftChild = 1; // changes value by reference 
		}
		
		if ( (target->getRightChild()->getMinChildCountIfSplitNEW() >= minPoints) 
			||
			// the right child would also have been in the container if it had enough
			// points and all of them went to one child, the other getting nothing
			((minPoints > 0)
			&& (target->getRightChild()->getCounter() >= minPoints)
			&& (target->getRightChild()->getMinChildCountIfSplitNEW() == 0)) ) {
				newNumLeaves--;
				deductRightChild = 1; // changes value by reference 
		}
	}
	else {
		 
		#ifdef DEBUG_NEWMCMC
			cout << "(target->nodeRealVolume() < 4*minVol) = " 
				<< ((target->nodeRealVolume() < 4*minVol)) << endl;
		#endif
	}
	#ifdef DEBUG_NEWMCMC
		cout << "adding 1 to leaves for merge,";
		if (deductLeftChild) cout << " deductLeftChild ";
		if (deductRightChild) cout << " deductRightChild ";
		cout << " newNumLeaves = " << newNumLeaves << endl;
		cout << "deducting 1 from cherries" << endl;
		
	#endif 

    size_t newNumCherries = numCherries - 1;
	
	/* but if target's sibling is a virtual leaf then parent would get added
	to the cherries - irrespective of whether sibling is a splittable leaf or not*/
   if (target->getParent() != NULL) { // no sibling if has no parent
	
		// want to find the target's sibling
		SPSnode* siblingNode = target->getParent()->getLeftChild();
		if (target == siblingNode) {
			siblingNode = target->getParent()->getRightChild();
		}
		std::string siblingNodeName = siblingNode->getNodeName();
		
		RealMappedSPnode::ListPtrsItr siblingIt;
		// break out of loop if we find sibling's name in 'real' leaves
        for (siblingIt=leaves.begin() ; siblingIt != leaves.end(); ++siblingIt ) {

            if ((*siblingIt)->getNodeName() == siblingNodeName) {
				// sibling is a virtual leaf, parent would be a new cherry
                ++newNumCherries;
				addParent = 1;
                break;
            }
        }
		#ifdef DEBUG_NEWMCMC
			cout << "node's sibling is = " << siblingNodeName << endl;
			if (addParent) cout << "And this is in current virtual leaves, so newNumCherries = " << newNumCherries << endl;
		#endif 
		
	}

    // Using proposal distribution object
    real deltaQ = proposal.getLogQRatioMergeProposal(numLeaves, numCherries,
                                                newNumLeaves, newNumCherries);
    //get another random number
    double randChange = gsl_rng_uniform(rgsl);
    if (log(randChange) < deltaPi + deltaQ) willMerge = true;

    if ((logging == TXT) || (logging == TXTANDGRAPH)) { // log these values
	
		#if defined OLDLOGSTYLE
			logMCMCDeltas(s, i, deltaL, deltaP, deltaQ, deltaPi, randChange);
			if (willSplit) outputFile(s, "Merging");
			else outputFile(s, "Not merging");
		#else // log with 2 for leaf as nodeType
			int nodeType = 20;
			if (!target->isSubLeaf()) nodeType = 21;
			logMCMCDeltasAugmented(s, i, nodeType, (willMerge ? 1: 0), 
					target->getCounter(),
					deltaL, deltaP, deltaQ, deltaPi, randChange);
		#endif
		
    }
	#ifdef DEBUG_NEWMCMC
		cout << "willMerge = " << willMerge << endl;
	#endif 

    return willMerge;
}



// change histogram state by splitting the target leaf node
void AdaptiveHistogram::changeStateForSplitNEW(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries,
						int addLeftChild, int addRightChild,
						int deductParent)
{
	#ifdef LOGEMPS
		size_t points = getRootCounter(); // need to recalculate COPERR EMP

		// accumulate the changes in scaled EMP sums that will result
		// from this expansion
		updateScaledEMPSumCOPERR(target->getSplitChangeEMPCOPERR(points));
		updateScaledEMPSumAIC(target->getSplitChangeEMPAIC());
	#endif
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn changeStateForSplitNEW, numLeaves = " << (numLeaves+1) << " numCherries = " << numCherries << endl;
		cout << " addLeftChild = " << addLeftChild;
		cout << " addRightChild = " << addRightChild;
		cout << " deductParent = " << deductParent << endl;
		if (target->isLeaf()) 
			cout << "target is a leaf : " <<( target->getNodeName() ) << endl;
		else 
			cout << "target is only a virtual leaf : " <<( target->getNodeName() ) << endl;

	#endif

    /* if the node is not already actually split
	 * then split the target and divvie up its data */
	try {
		if (target->isLeaf()) target->nodeExpand();
		
	}
	catch (UnfulfillableRequest_Error& ure) {
		std::cerr << "Something has gone wrong in the splitting.  Error reported is:"
			<< std::endl;
		std::cerr << ure.what()	<< std::endl;
		std::cerr << "(It may be that a split is requested on a box so narrow that child == parent)\n" << std::endl;
		
		throw;								// rethrow the exception
	}
	
	// but only put the children into the container if they can be split
    if (addRightChild) {

        // insert the new children ptrs into the list at the beginning
        nodes.push_front(target->getRightChild());
        ++numLeaves;
    }

    if (addLeftChild) {
        // insert the new children ptrs into the list at the beginning
        nodes.push_front(target->getLeftChild());// left goes first
        ++numLeaves;
    }

	/*If the parent is in the cherries it must come out.  */
	if (deductParent) {
		
		bool foundParent = false;
	
		SPSnode* nodeParent = target->getParent();
	
		SPSnodeListItr it = nodes.begin();
		advance(it, numLeaves); // advance to the cherries
		// break out of loop if we find parent
		for ( ; it != nodes.end(); ++it ) {
			if ((*it) == nodeParent) {
				nodes.erase(it);
				numCherries--;
				foundParent = true;
				break;
			}
		}
		if (!foundParent) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForSplitNEW(...) : parent lost");
		}
	}
    // put this node ptr into the cherries, ie at end of list
    nodes.push_back(target);
    ++numCherries;

	#ifdef DEBUG_NEWMCMC
		cout << "End changeStateForSplitNEW, numLeaves = " << numLeaves << " numCherries = " << numCherries << endl;
		
	#endif
}

// change histogram state by merging the target cherry node
void AdaptiveHistogram::changeStateForMergeNEW(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries,
						int deductLeftChild, int deductRightChild,
						int addParent)
{
	#ifdef LOGEMPS
		size_t points = getRootCounter(); // need to recalculate COPERR EMP

		// accumulate the changes in scaled EMP sums that will result
		// from this expansion
		updateScaledEMPSumCOPERR(target->getMergeChangeEMPCOPERR(points));
		updateScaledEMPSumAIC(target->getMergeChangeEMPAIC());
	#endif
	
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn changeStateForMergeNEW, numLeaves = " << numLeaves << " numCherries = " << (numCherries+1) << endl;
		cout << " deductLeftChild = " << deductLeftChild;
		cout << " deductRightChild = " << deductRightChild;
		cout << " addParent = " << addParent << endl;
		if (target->isSubLeaf()) 
			cout << "target is a sub leaf : " <<( target->getNodeName() ) << endl;
		else 
			cout << "target is only a virtual subleaf : " <<( target->getNodeName() ) << endl;

	#endif
	
    // take the children out of the list of virtual leaves if they are there
	int lcPos = 0;
   
    if (deductLeftChild) {
		
		bool foundLeft = false;
        
		SPSnode* lcNode = target->getLeftChild();
        SPSnodeListItr it;
		
		SPSnodeListItr end_it = nodes.end();
		advance(end_it, -numCherries);
        // break out of loop if we find left child or get to cherries
        for ( SPSnodeListItr it = nodes.begin() ; it != end_it; ++it ) {

            if ((*it) == lcNode) {
                nodes.erase(it); // can't keep using iterator now
                numLeaves--;
                foundLeft = true;
                break;
            }
            ++lcPos;  // gives position at which lc was found
        }
		if (!foundLeft) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForMergeNEW(...) : left child lost");
		}
    }
    
    // now try to find right child - could be immediately after left
    if (deductRightChild) {

		bool foundRight = false;

    	SPSnodeListItr it = nodes.begin();
		advance(it, lcPos);
		// break out of loop if we find right child or get to cherries
		SPSnodeListItr end_it = nodes.end();
		advance(end_it, -numCherries);
       
		SPSnode* rcNode = target->getRightChild();
        
		for ( ; it != end_it; ++it ) {

			if ((*it) == rcNode) {
				nodes.erase(it);
				numLeaves--;
				foundRight = true;
				break;
			}
		} // just in case right child was before left
		if (!foundRight) {
			// break out of loop if we find right child
			// or get to cherries
			for (it=nodes.begin() ; it != end_it; ++it ) {
				if ((*it) == rcNode) {
					nodes.erase(it);
					numLeaves--;
					foundRight = true;
					break;
				}
			}
		}

		if (!foundRight) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForMergeNEW(...) : right child lost");
		}
    }
    
    // BUT we do NOT actually merge the target
   
    // insert the new "virtual leaf" ptr into the list at the beginning
    nodes.push_front(target);
    ++numLeaves;

	// if sibling is a virtual leaf, add parent to cherries
    if (addParent) {
	
		nodes.push_back(target->getParent());
		++numCherries;
        
	}
	#ifdef DEBUG_NEWMCMC
		cout << "End changeStateForMergeNEW, numLeaves = " << numLeaves << " numCherries = " << numCherries << endl;
		
	#endif


}


// change realmapped paving state by splitting the target leaf node
void AdaptiveHistogram::changeStateForSplitRMSP(
						const SPSnode * const spsTarget,
                        RealMappedSPnode::ListPtrs& leaves, 
						RealMappedSPnode::ListPtrs& cherries, 
						const size_t leftChildCount,
						const size_t rightChildCount,
						AdaptiveHistogram::ChangeOfStateInformation& info)
{
	#ifdef DEBUG_NEWMCMC
		cout << "\nIn changeStateForSplitRMSP, ";
		cout << " target is " << spsTarget->getNodeName() << endl;
		cout << " leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		
	#endif
	
	// find the equivalent to the sps target in the rmsp
	RealMappedSPnode* target = NULL;
	{
		std::string spsName = spsTarget->getNodeName();
		bool foundTarget = false;
		for (RealMappedSPnode::ListPtrsItr it = leaves.begin();
				it != leaves.end();
				++it) {
			
			if ( (*it)->getNodeName() == spsName) {
				target = (*it);		
				leaves.erase(it);
				foundTarget = true;
				break;
			}
		}
		if (!foundTarget || (target == NULL)) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForSplitRMSP(...) : target rmsp lost");
		}
	}
	
	#ifdef DEBUG_NEWMCMC
		cout << "After finding and erasing MSP target = " << (target->getNodeName()) << endl;
		cout << "leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		cout << "leftChildCount = " << leftChildCount;
		cout << " rightChildCount = " << rightChildCount << endl;
		if (target->isLeaf()) 
			cout << "target is a leaf" << endl;
		else 
			cout << "target is not a leaf -- OHHHH NOOOO" << endl;

	#endif
	
    //notify the info object that this target is about to be split
	info.notifySplit(target);                    
	
    // split the target rmsp node
	try {
		target->nodeExpand();
	}
	catch (UnfulfillableRequest_Error& ure) {
		std::cerr << "Something has gone wrong in the splitting.  Error reported is:"
			<< std::endl;
		std::cerr << ure.what()	<< std::endl;
		std::cerr << "(It may be that a split is requested on a box so narrow that child == parent)\n" << std::endl;
		
		throw;								// rethrow the exception
	}
	real vol = target->nodeRealVolume(); 
	
	//range (unnormalised) will be new count / (0.5vol) = 2*count/vol
	target->getLeftChild()->setRange(2.0*leftChildCount/vol);
	target->getRightChild()->setRange(2.0*rightChildCount/vol);
	
	leaves.push_back(target->getLeftChild());
	leaves.push_back(target->getRightChild());
	
	// if sibling was a leaf, take parent out of cherries
    if(target->hasLeafSibling()) {

        // how to find parent? - search the cherries?
        bool foundParent = false;
        RealMappedSPnode* nodeParent = target->getParent();

        RealMappedSPnode::ListPtrsItr it = cherries.begin();
        for (RealMappedSPnode::ListPtrsItr it = cherries.begin();
				it != cherries.end();
				++it ) {
			if ((*it) == nodeParent) {
                leaves.erase(it);
                foundParent = true;
                break;
            }
        }
        if (!foundParent) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForSplitRMSP(...) : parent cherry lost");
		}
		#ifdef DEBUG_NEWMCMC
			cout << "Target had a leaf sibling, so taking parent " << nodeParent->getNodeName() << " out of cherries " << endl;

		#endif

    }

	#ifdef DEBUG_NEWMCMC
		cout << "Adding new left child " << target->getLeftChild()->getNodeName() << " to leaves " << endl;
		cout << "Adding new right child " << target->getRightChild()->getNodeName() << " to leaves " << endl;
		cout << "Put this node into the cherries " << endl;
	#endif
    // put this node ptr into the cherries
    cherries.push_back(target);
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nAt end of changeStateForSplitRMSP, ";
		cout << " leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		
	#endif
}

// change realmapped paving state by splitting the target leaf node
void AdaptiveHistogram::changeStateForMergeRMSP(
						const SPSnode * const spsTarget,
                        RealMappedSPnode::ListPtrs& leaves, 
						RealMappedSPnode::ListPtrs& cherries,
						AdaptiveHistogram::ChangeOfStateInformation& info)
{

	#ifdef DEBUG_NEWMCMC
		cout << "\nIn changeStateForMergeRMSP";
		cout << " leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		
	#endif
	
	// find the equivalent to the sps target in the rmsp
	RealMappedSPnode* target = NULL;
	{
		std::string spsName = spsTarget->getNodeName();
		bool foundTarget = false;
		for (RealMappedSPnode::ListPtrsItr it = cherries.begin();
				it != cherries.end();
				++it) {
			
			if ( (*it)->getNodeName() == spsName) {
				target = (*it);		
				cherries.erase(it);
				foundTarget = true;
				break;
			}
		}
		if (!foundTarget || (target == NULL)) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForMergeRMSP(...) : target rmsp lost");
		}
	}
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nAfter finding and erasing RMSP target = " << (target->getNodeName()) << endl;
		cout << "leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		if (target->isSubLeaf()) 
			cout << "target is a subleaf" << endl;
		else 
			cout << "target is not a subleaf -- OHHHH NOOOO" << endl;

	#endif
    
	// take the children out of the leaves 
	int lcPos = 0;
   
    { // left child
		
		bool foundLeft = false;
        
		RealMappedSPnode* lcNode = target->getLeftChild();
        // break out of loop if we find left child
        for (RealMappedSPnode::ListPtrsItr it =leaves.begin();
					it != leaves.end();
					++it ) {

           if ((*it) == lcNode) {
                leaves.erase(it); // can't keep using iterator now
                foundLeft = true;
                break;
            }
            ++lcPos;  // gives position at which lc was found
        }
		if (!foundLeft) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForMergeRMSP(...) : left child lost");
		}
    }
    
    // now try to find right child - should be immediately after left
    {

		bool foundRight = false;

    	RealMappedSPnode::ListPtrsItr it = leaves.begin();
		advance(it, lcPos);
		// break out of loop if we find right child or get to cherries
		
		RealMappedSPnode* rcNode = target->getRightChild();
        
		for ( ; it != leaves.end(); ++it) {

			if ((*it) == rcNode) {
				leaves.erase(it);
				foundRight = true;
				break;
			}
		} // just in case right child was before left
		if (!foundRight) {
			// break out of loop if we find right child
			// or get to cherries
			for (it=leaves.begin() ; it != leaves.end(); ++it ) {
				if ((*it) == rcNode) {
					leaves.erase(it);
					foundRight = true;
					break;
				}
			}
		}

		if (!foundRight) {
			throw std::logic_error(
				"AdaptiveHistogram::changeStateForMergeRMSP(...) : right child lost");
		}
    }
    
    //notify the info oject that this target is about to be merged
	info.notifyMerge(target);                    

	// merge the target
    target->nodeReabsorbChildren();
    
    // insert the new leaf ptr into the leaves
    leaves.push_back(target);
	#ifdef DEBUG_NEWMCMC
		cout << "Taking old left child and right child out of leaves " << endl;
		cout << "Put this node into the leaves " << endl;
	#endif
    
    // if sibling is a leaf, add parent to cherries
    if(target->hasLeafSibling()) { // returns false if no parent
        cherries.push_back(target->getParent());
		
		#ifdef DEBUG_NEWMCMC
			cout << "Adding the parent " << target->getParent()->getNodeName() << " to the cherries" << endl;
		#endif
    }
	
	#ifdef DEBUG_NEWMCMC
		cout << "\nAt end of changeStateForSplitRMSP, ";
		cout << " leaves.size() = " << (leaves.size());
		cout << " cherries.size() = " << (cherries.size()) << endl;
		
	#endif
}



//NEWAPRIL2012
// Get change in scaled contribution to log likelihood of a node on a split
real AdaptiveHistogram::getSplitChangeLogLik(
				size_t leftCount, size_t rightCount)
{
	
	dotprecision change(0.0);
	
	size_t totalCount = leftCount + rightCount;
	
	// if total count is 0 there can be no change on splitting
	if (totalCount > 0) {

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

		accumulate(childEMP, 1.0*totalCount, log(2.0));

		accumulate(currentEMP, 1.0*totalCount, log(1.0*totalCount));

		change = childEMP - currentEMP;
	}
	return rnd(change);
}

// internal method to set up priority queue
std::multiset< SPSnode*, subpavings::MyCompare>&
	AdaptiveHistogram::_setupPrioritySplit(
		std::multiset< SPSnode*, subpavings::MyCompare>& pq,
        size_t minChildPoints, double minVol)
{

	bool legalState = true; // use to check current state

	if (minVol > 0.0) {
		legalState = getSubPaving()->checkTreeStateLegal(
												minChildPoints, minVol);
	}
	else { // not vol checking
		legalState = getSubPaving()->checkTreeStateLegal(minChildPoints);
	}
	
	// check current state is legal
	if (!legalState) {
		throw std::logic_error(
			"_setupPrioritySplit(...) : Illegal starting state");
	}
	
	// put nodes into the starting set IF they meet minVol test AND IF either
	// there are enough points in the whole node
			// and minChildCountIfSplit is 0 (ie all points go to one child)
	// or the minChildCountIfSplit test passed

	if (getSubPaving()->isLeaf()) {
		#ifdef MYDEBUG
			cout << "root is a leaf" << endl;
			
		#endif
		// check to insert a copy of the rootPaving pointer into the set
		if (getSubPaving()->isSplittableNode(minChildPoints, minVol)) {
				pq.insert(getSubPaving());
				#ifdef MYDEBUG
					cout << "inserted to pq, pq.size = " << pq.size() << endl;
					
				#endif
		}
	}
	else { // root is not a leaf
	
		#ifdef MYDEBUG
			cout << "root is not a leaf" << endl;
			
		#endif
		SPSnodePtrs leaves;
		getSubPaving()->getLeaves(leaves);
		// check to insert each of the leaves into the set
		#ifdef MYDEBUG
			cout << "root has " << leaves.size() << " leaves" << endl;
			
		#endif
		SPSnodePtrsItr sit;
		for (sit = leaves.begin(); sit < leaves.end(); ++sit) {
			if ((*sit)->isSplittableNode(minChildPoints, minVol)) {
				pq.insert(*sit);
			}
		}
		
	}
	#ifdef MYDEBUG
		cout << "pq now has size " << pq.size() << endl;
		
	#endif
		
	return pq;
}
	

// internal method do do priority split loop
bool AdaptiveHistogram::_prioritySplitLoop(
								std::multiset< SPSnode*, subpavings::MyCompare>& pq,
								size_t n,
                                size_t minChildPoints, double minVol,
                                gsl_rng * rgsl)
{
	bool canContinue = true;
	
	if (!pq.empty()) {
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
		dotprecision changeCOPERR = chosenLargest->getSplitChangeEMPCOPERR(n);
		dotprecision changeAIC = chosenLargest->getSplitChangeEMPAIC();

		// split the biggest one and divvie up its data

		chosenLargest->nodeExpand();
			
		updateScaledEMPSumCOPERR(changeCOPERR);
		updateScaledEMPSumAIC(changeAIC);

		
		/*but only put the children into the container if they can be
		 split */
		if (chosenLargest->getLeftChild()->isSplittableNode(
					minChildPoints, minVol)) {
			// insert the new left child into the multiset
			pq.insert(chosenLargest->getLeftChild());
		}

		if (chosenLargest->getRightChild()->isSplittableNode(
					minChildPoints, minVol)) {
			// insert the new right child into the multiset
			pq.insert(chosenLargest->getRightChild());
		}
	}
	else {
		canContinue = false;
		std::cerr << "Terminated splitting: no splittable nodes left"
				<< std::endl;
	}
	
	return canContinue;
}

// internal method to set up priority merge
std::multiset< SPSnode*, MyCompare>&
 AdaptiveHistogram::_setupPriorityMerge(
						std::multiset< SPSnode*, MyCompare>& pq)
{
	SPSnodePtrs subleaves;
	getSubPaving()->getSubLeaves(subleaves);

	// insert a copy of the current set of subleaves into the multiset
	pq.insert(subleaves.begin(), subleaves.end());

	return pq;
}


// internal method for priority merge loop
bool AdaptiveHistogram::_priorityMergeLoop(
							std::multiset< SPSnode*, MyCompare>& pq,
							size_t n,
							gsl_rng * rgsl)
{
	bool canContinue = true;
	
	if (!pq.empty()) {
		
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

		
		// merge the biggest one
		chosenSmallest->nodeReabsorbChildren();
		pq.erase(mit);// take the iterator to chosen smallest out of the set

		// if smallest had a leaf sibling, smallest's parent is now a cherry
		// and should be inserted into the multiset
		if (chosenSmallest->hasLeafSibling()) {

			pq.insert(chosenSmallest->getParent());
	    }

    }
	
	else { // pq empty
		std::cerr << "No subleaves left to merge" << std::endl;
		
		canContinue = false;
	}
	
	return canContinue;  
}

	

// internal method to find coverage
double AdaptiveHistogram::_coverage(const rvector& pt) const
{
	double cov = 0.0;
	
	// if there is nothing in the histogram, coverage will always be 0.0
	
	if ( getRootCounter() ) {
	
		const SPSnode * container = 
				getSubPaving()->findContainingNode(pt);
		if (container != NULL) {
			
			size_t totCount = getRootCounter();
			size_t culmCounts = totCount;
		
			// put the leaves into a vector and sort it, smallest to largest
			SPSnodePtrs leaves;
			getSubPaving()->getLeaves(leaves);
			CompHeight compheight;
			
			sort(leaves.begin(), leaves.end(), MyCompare(compheight));
			
			bool found = FALSE;
			
			SPSnodePtrs::reverse_iterator rit = leaves.rbegin();
			
			double containerHeightNonNorm = 
						container->getCountOverVolume();
			// height is box counts/box vol
			// ie count is unnormalised vol of an individual element of histogram
			// box vol * height == box vol * (count / box vol) == count
		
			while (!found && rit < leaves.rend()) {
				
				double thisHeightNonNorm = (*rit)->getCountOverVolume();
				
				// check the boxes
				// stop at the first box not taller than this one
				// this makes sure answer is not influenced by chance
				// ordering of nodes of same height
				
				if (thisHeightNonNorm > containerHeightNonNorm) {
					culmCounts -= (*rit)->getCounter();
				}
				/*
				if ( fabs( thisHeightNonNorm 
						- containerHeightNonNorm ) > 
						DBL_EPSILON*max(thisHeightNonNorm,containerHeightNonNorm) ) {
					// decrement cumulative counts (unnormalised histgram volume)
					culmCounts -= (*rit)->getCounter();
				}
				*/
				else { 
					found = TRUE;	// break out of loop
				}
				++rit;
				
			} // end while
			
			// if we have not found we have a problem,
			// since findContainingNode said it was here somewhere
			
			if (!found) {
				throw std::logic_error("_coverage(const rvector&) : lost containing node");
			}
			cov = static_cast< double > (culmCounts)/getRootCounter();
		}
	}
	return cov;
}

// internal method to find empiricalDensity
double AdaptiveHistogram::_empiricalDensity(const rvector& pt) const
{
	double den = 0.0;
	
	// if there is no data in the histogram, empirical density will always be 0
	
	if ( getRootCounter() ) {
	
		const SPSnode * container = getSubPaving()->findContainingNode(pt);
		if (container != NULL) {
			
			den = container->getCountOverVolume()/getRootCounter();
			
		}
	}
	return den;
}


RealVec AdaptiveHistogram::calculateVarCovarFromBigData() const
{
	int dimension = getDimensions();
	VecDotPrec dpSums(dimension, cxsc::dotprecision(0.0));
	VecDotPrec dpSumProducts(dimension*dimension, cxsc::dotprecision(0.0));
	
	/* the sumproducts can be thought of as an nxn matrix,
	which is implemented here as a nxn element vector of
	dotprecision variables, using row-major order.
	Ie the m-th element (m = 0, . . . nxn-1) is in row floor(m/n)
	and column m-rowxn in the matrix configuration.
	Or, the sumproduct of elements i and j in an rvector,
	i,j = 0,...,n-1, is element m=(ixn+j) of the sumproducts
	vector. */
	
	for (BigDataConstItr it = dataCollection.begin();
			it != dataCollection.end();
			++it) {
		
		rvector newdata = *it;
		// make a dot precision variable out of the ith element
		// and jth element of the of the rvector of new data and
		// store in dpSumProducts.
		for (int i = 1; i < dimension + 1; ++i) {
			
			accumulate(dpSums[i-1], newdata[i], 1.0);
			
			// only need to do columns 1 to i because of symmetry
			for (int j = 1; j< i + 1; ++j) {

				int index = (i-1)*dimension + (j-1);
				// rvectors indexed 1 to n
				accumulate(dpSumProducts[index],
						newdata[i], newdata[j]);

				//if not on the diagonal of the matrix,
				// we can also fill in the symmetric element
				if (i!=j) {
					int sym_index = (j-1)*dimension
						+ (i-1);
					dpSumProducts[sym_index] =
						dpSumProducts[index];
				} // end if
			}// end j-loop
		}// end i-loop
	} // end iteration through data
	
	RealVec varCovar;
	varCovar.reserve(dimension*dimension);
	
	size_t n = dataCollection.size();

	// loop through the elements in the dpSumProducts vector
	for (int k = 0; k < dimension*dimension; k++) {

		/*the var-covar is the sample var-covar
		which is
		[sumproduct(i,j)-sum(i)sum(j)/counter]/(counter-1)

		element k in the vector of dotprecison sumproducts
		corresponds to row k/n, (row 0 to n-1)
		and column k-row*n (col 0 to n-1)
		in a matrix view of the sumproducts */
		if (n < 2) varCovar.push_back(cxsc::SignalingNaN);
		else {
			int i = k/dimension; // row  (int/int = int)
			int j = k - i*dimension; // column

			// make another dotprecision variable
			dotprecision temp1 = dpSumProducts[k];

			dotprecision temp2(0.0);
			// sum(i) x sum(j)
			// default cxsc rounding dotprecision rnd_next
			accumulate(temp2,  rnd(dpSums[i]),
					rnd(dpSums[j]));

			real div = -1.0/n;

			// sumproduct(i,j) - sum(i)(sum(j)/counter
			// default cxsc rounding
			accumulate(temp1, rnd(temp2), div);
			// calculate the variance covariance element
			varCovar.push_back(rnd(temp1)/(1.0*(n-1)));

		}
	}// end loop through the elements in dpSumProducts

	return varCovar;
}

cxsc::rvector AdaptiveHistogram::calculateMeanFromBigData() const
{
	int dimension = getDimensions();
	RealVec sums(dimension, cxsc::real(0.0));
	
	for (BigDataConstItr it = dataCollection.begin();
			it != dataCollection.end();
			++it) {
		rvector newdata = *it;
		
		for (int i = 1; i < dimension + 1; ++i) {
			
			sums[i-1]+=newdata[i];
			
		}// end i-loop
	} // end iteration through data
	
	cxsc::rvector mean(dimension);
	size_t n = dataCollection.size();
	for (int i = 0; i < dimension; ++i) {
			
			if (n == 0) mean[i+1] = cxsc::SignalingNaN;
			else mean[i+1] = sums[i]/(1.0*n);
			
	}// end i-loop
	
	return mean;
}


//check that the box is okay for the histogram
bool AdaptiveHistogram::checkBox(const cxsc::ivector& box)
{
	return subpavings::checkBox(box);
}

void AdaptiveHistogram::handleSplitToShapeError(SPSnode& spn)
{
	// restore our spn to the supplied copy
	std::swap(*(getSubPaving()), spn);
	
	std::cerr << std::endl;
			std::cerr << "Your instruction does not describe a proper tree.";
			std::cerr << "  Please check your instruction and try again."
			<< std::endl;
}

// ensure rootPaving is deleted if constructed in failed constructor
void AdaptiveHistogram::constructor_error_handler() 
{
	try {
		
			delete rootPaving;
			rootPaving = NULL;
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}




/* true if there were splittable leaves at end of process, false
 * if process terminated because we run out of splittable leaves */
bool AdaptiveHistogram::_prioritySplitWithPosterior(
					PriorityQueueT& pq,
					const PrivatePrioritySplitQueueEvaluator& ppsqe,
					LOGGING_LEVEL logging,
					size_t minChildPoints,
					double minVol, 
					LogMCMCPrior& logPrior,
					std::vector<real>& posteriorVec,
					std::vector<real>& loglikVec,
					size_t n,
					size_t numLeaves,
					real lnPrior,
					real lnLik,
					bool stopOnMaxPosterior,
					size_t& leavesForMaxPosterior,
					real& maxPosterior,
					gsl_rng * rgsl)
{
	string errorMsg("AdaptiveHistogram::_prioritySplitWithPosterior(...)");
	if (!hasSubPaving()) {
		throw NullSubpavingPointer_Error(errorMsg);
	}
	
	#ifdef MYDEBUG_MCMCPQ_CHECK_QUEUE
		cout << "entering _prioritySplitWithPosterior" << endl;
		cout << "size of set is " << pq.size() << " and set is " << endl;
		for (PriorityQueueItrT it = pq.begin();
			it != pq.end(); ++it) {
				cout << it->toString() << endl;
			}
		
		gsl_rng * rgsl_tmp = gsl_rng_clone (rgsl);
		
		cout << "The next rand would be " << gsl_rng_uniform(rgsl_tmp) << endl;
	
		gsl_rng_free(rgsl_tmp);
	
	#endif
	
	/* a parameter to control how many states we will look at, at least,
	 * if stopOnMaxPosterior, before we decide to quit if the posterior
	 * has fallen below the starting posterior.  Should not really be hard
	 * coded here */
	const size_t minStatesToStop = 500;
	
	real posterior = lnLik + lnPrior;
	// store
	loglikVec.push_back(lnLik);
	posteriorVec.push_back(posterior);

	// maxs at this point
	size_t lastLeavesForMaxPosterior = numLeaves;
	real lastMaxPosterior = posterior;
	real firstMaxPosterior = posterior;
	size_t firstNumLeaves = numLeaves;
	
	#ifdef MYDEBUG_MCMCPQ
		cout << "entering _prioritySplitWithPosterior" << endl;
		cout << "numleaves is " <<numLeaves << " and posterior is " << posterior << endl;
		
	
	#endif

	std::string baseFileName = "";
	std::string s = "";
	if (logging != NOLOG) {
		// pass to log output to keep track of splits
		baseFileName = "pqOutput";
		s = getUniqueFilename(baseFileName);
	}
	
	// a counter for the loop
	int i = 0;

	if (logging != NOLOG) {
		 // Start log file with filename and timestamp
		outputLogStart(s);
		// log the current state of the histogram
		outputLogPlain(s, i);
		#ifdef LOGEMPS
			outputLogEMPAIC(s); // add AIC scores
		#endif
		
	}
	if (logging == LOGSAMPLES) {
		ostringstream oss;
		oss << "CarverSEBQueueState_" << numLeaves << ".txt";
		outputToTxtTabs(oss.str());
		
	}
	i++;

	
	bool canContinue = !pq.empty(); 
	if(!canContinue) {
		std::cerr << "No splittable leaves to split - aborting" << std::endl;
	}
	
	/* split until the ppsqe says to stop */ 
	while (canContinue && (!ppsqe.stopQueueQuery(pq, numLeaves))) {

		real loopEmptyBoxVol(0.0);
		real deltaL(0.0);
		
		canContinue = _prioritySplitLoopWithEmptyVolMeasure(pq, ppsqe.measurer,
								n, minChildPoints, minVol,
								loopEmptyBoxVol, deltaL, rgsl);

		if (logging != NOLOG) {
			// To add current state of histogram to log file
			outputLogPlain(s, i);
			#ifdef LOGEMPS
				outputLogEMPAIC(s); // add AIC scores
			#endif
			
		}
		i++;

		/* if loop was successful number of leaves will increased by 1*/
		if (canContinue) {
			real deltaP = logPrior.changeOnSplitOne(numLeaves-1);
			numLeaves++;
		
			
			// update the logprior
			lnPrior += deltaP;
			
			// update the 'real' likelihood
			lnLik += deltaL;

			posterior = lnLik + lnPrior;
			// store
			loglikVec.push_back(lnLik);
			posteriorVec.push_back(posterior);
			
			// check for new max
			if (posterior > lastMaxPosterior) {
				lastLeavesForMaxPosterior = numLeaves;
				lastMaxPosterior = posterior;
			}
			// if we've started to go down we might want to stop now
			else if (stopOnMaxPosterior && (posterior < firstMaxPosterior) 
							&& (numLeaves-firstNumLeaves > minStatesToStop)) {
				// break out of the loop: canContinue will still be true
				#ifdef MYDEBUG_MCMCPQ
					cout << "breaking out :" << endl;
					cout << "numleaves is " <<numLeaves << " and posterior is " << posterior << endl;
					cout << "last posterior is " << lastMaxPosterior << endl;
					
				
				#endif
				break; 
			}
			

			if (logging == LOGSAMPLES) {
				ostringstream oss;
				oss << "CarverSEBQueueState_" << numLeaves << ".txt";
				outputToTxtTabs(oss.str());
				
			}
			#ifdef MYDEBUG_MCMCPQ_EXTRA
				cout << "i = " << i << "\tnumLeaves\t" << numLeaves << "\tlnLik\t" <<lnLik << endl;
			#endif
		}

		
		if (canContinue && (logging != NOLOG)) {
			// log the leaf levels line
			outputFile(s, getLeafLevelsString());
		}
	}
	// EMPSums will have been adjusted during the splitting process
	
	leavesForMaxPosterior = lastLeavesForMaxPosterior;
	maxPosterior = lastMaxPosterior;
	
	return canContinue; 
}

//new
bool AdaptiveHistogram::_launchPrioritySplitWithPosterior(
					PriorityQueueT& stepPQ,
					const PrioritySplitQueueEvaluator& psqePosterior,
					LOGGING_LEVEL logging,
					size_t minChildPoints,
					double minVol, 
					LogMCMCPrior& logPrior,
					std::vector<real>& posteriorVec,
					std::vector<real>& loglikVec,
					size_t n,
					size_t stepNumLeaves,
					real stepLnPrior,
					real stepLnLik,
					bool stopOnMaxPosterior,
					std::vector<size_t>& maxPosteriorPoints,
					std::vector<real>& maxPosteriors,
                    gsl_rng * rgsl_step)
{
	size_t leavesForMaxPosterior = 0;
	real maxPosterior(0.0);
	
	#ifdef MYDEBUG_MCMCPQ
	cout << "\nLaunching a posterior queue:  size of queue coming in is\t" 
		<< stepPQ.size() << "\tstepNumLeaves is " << stepNumLeaves << endl; 
	#endif
	
	/* remake the queue with the new measure */
	PrivatePrioritySplitQueueEvaluator ppsqePosterior(psqePosterior);
	
	PriorityQueueT pqNew;
	
	_resetPrioritySplitWithEmptyVolMeasure(stepPQ, pqNew, 
										ppsqePosterior.measurer);
	
	bool success = _prioritySplitWithPosterior(
					pqNew, // new queue
					ppsqePosterior, 
					logging,
					minChildPoints, 
					minVol,
					logPrior,
					posteriorVec,
					loglikVec,
					n, stepNumLeaves,
					stepLnPrior, stepLnLik,
					stopOnMaxPosterior,
					leavesForMaxPosterior,
					maxPosterior,
					rgsl_step);
	
	if (!success) {
		cerr << "Note: posterior search started at leaves "
		<< stepNumLeaves << " ended early did not return success" << endl;
	}
	/* If leavesForMaxPosterior == stepNumLeaves
	 * that means posterior went down right from the start
	 * If leavesForMaxPosterior - stepNumLeaves + 1
	 * == posteriorVec.size() that indicates 
	 * that we were still going up when we ended the 
	 * queue. */
	if (leavesForMaxPosterior == stepNumLeaves) {
		cerr << "Note: posterior search started at leaves "
		<< stepNumLeaves 
		<< " never found a higher max than its starting point" << endl;
	} 
	if ((leavesForMaxPosterior == stepNumLeaves) ||
		(leavesForMaxPosterior - stepNumLeaves + 1 < 
		posteriorVec.size())) {
		
		maxPosteriorPoints.push_back(leavesForMaxPosterior);
		maxPosteriors.push_back(maxPosterior);
		#ifdef MYDEBUG_MCMCPQ_MIN
			cout << "Logging max posterior" << maxPosterior 
				<< " at leaves " << leavesForMaxPosterior << endl;;
		#endif	
	}
	else { 
	
		#ifdef MYDEBUG_MCMCPQ_MIN
			cout << "Failed to find a proper max posterior" << endl;;
		#endif	
		
		#if(0)
			// put in the largest we can
			maxPosteriorPoints.push_back(ULONG_MAX);
			maxPosteriors.push_back(cxsc::Infinity);
		
		#endif
		#if(1)
			/* log the values we did get: this means that if we keep reaching the 
			 * max we will at least get some thing 'valid' */
			maxPosteriorPoints.push_back(leavesForMaxPosterior);
			maxPosteriors.push_back(maxPosterior);
		#endif
		
	}
	
	
	#ifdef MYDEBUG_MCMCPQ
		ostringstream oss;
		oss << "PosteriorsFrom_" << stepNumLeaves << ".txt";
		string postFileName = oss.str();
		std::vector < const subpavings::RealVec* > dataPtrs;
		dataPtrs.push_back(&posteriorVec);
		dataPtrs.push_back(&loglikVec);
		std::vector < std::string > colNames;
		colNames.push_back("Posterior");
		colNames.push_back("lnLik");
		
		outputToFileVertical(dataPtrs, colNames, postFileName);
		
		cout << "\n\nBack in original, posteriors are in " 
			<< postFileName << "\n" << endl;
	#endif
	
	return success;
}	







//based on Gloria
/* internal method to set up priority queue when we want to know a
about empty volumes */
real AdaptiveHistogram::_setupPrioritySplitWithEmptyVolMeasure(
		PriorityQueueT& pq,
		SPSNodeMeasure& measurer,
        size_t minChildPoints,
		double minVol 
		) const
{
	//return value is total volume of empty boxes when queue is set up
	real emptyBoxVol = 0.0;
	
	bool legalState = getSubPaving()->checkTreeStateLegal(minChildPoints,
															minVol);
	
	// check current state is legal
	if (!legalState) {
		throw std::logic_error(
			"_setupPrioritySplitVol(...) : Illegal starting state");
	}
	
	/* put nodes into the starting set IF they are not empty 
	   AND IF either
	   there are enough points in the whole node
			and minChildCountIfSplit is 0 (ie all points go to one child)
	   or the minChildCountIfSplit test passed */

	
	SPSnodePtrs leaves;
	getSubPaving()->getLeaves(leaves);
	// check to insert each of the leaves into the set
	#ifdef MYDEBUG
		cout << "root has " << leaves.size() << " leaves" << endl;
		
	#endif
	SPSnodePtrsItr sit;
	for (sit = leaves.begin(); sit < leaves.end(); sit++) {
		//do not put empty nodes into the starting set
		if ( (*sit)->getCounter() > 0 ) {
		
			/* isSplittableNode chekcs that node has >=2*minVol, ie
			 * child will have vol >= minVol */
			if ((*sit)->isSplittableNode(minChildPoints,
										minVol)) { 
			
				real m = measurer( (*sit) );
				pq.insert( NodePtrMeasurePair((*sit), m) );
			}
		}
		else {
			emptyBoxVol += (*sit)->nodeVolume();
		}
	}
		

	#ifdef MYDEBUG
		cout << "pq now has size " << pq.size() << endl;
		
	#endif
		
	return emptyBoxVol;
}



/* internal method to reset a priority queue using a different measure */
void AdaptiveHistogram::_resetPrioritySplitWithEmptyVolMeasure(
		const PriorityQueueT& pqOld,
		PriorityQueueT& pqNew,
		SPSNodeMeasure& measurer)
{
	// comes from a current queue - we are just changing the measure
	for (PriorityQueueConstItrT it = pqOld.begin(); it != pqOld.end(); ++it) {
		real m = measurer(it->nodePtr);
		pqNew.insert( NodePtrMeasurePair(it->nodePtr, m) );
	}	
	assert(pqNew.size() == pqOld.size());
	
	#ifdef MYDEBUG_MCMCPQ
		cout << "\n reset the queue to new measure: queue size " << pqOld.size() << endl;
	#endif
	
	
		
}



/* return false if no splittable nodes left, true otherwise */
bool AdaptiveHistogram::_prioritySplitLoopWithEmptyVolMeasure(
					PriorityQueueT& pq,
					SPSNodeMeasure& measurer,
					size_t n, // only needed for emps
					size_t minChildPoints, 
					double minVol, 
					real& loopEmptyBoxVolume, real& deltaL, 
					gsl_rng * rgsl)
{
	bool canContinue = true;
	
	//initialize empty box vol and deltaL
	loopEmptyBoxVolume = 0.0;
	deltaL = 0.0;

	if (!pq.empty()) {
		
		NodePtrMeasurePair largest = *(pq.rbegin ()); // the last largest in the set
		
		// find if there are any more equal to largest around
		pair< PriorityQueueItrT, PriorityQueueItrT > equalLargest;

		equalLargest = pq.equal_range(largest); // everything that = largest
		size_t numberLargest = pq.count(largest); // number of =largest
		assert(numberLargest > 0);

		PriorityQueueItrT mit = equalLargest.first;
		
		if (numberLargest > 1) {

			// draw a random number in [0,1)
			double rand = gsl_rng_uniform(rgsl);
			
			real sum = 0.0;
			bool found = false;
			// random selection of the =largest node to chose
			for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
				
				sum += 1.0/(1.0*numberLargest);
				#ifdef MYDEBUG_MCMCPQ_LOOP
					cout << "considering " << mit->toString();
					cout << ":\t sum = " << sum << endl;
				#endif
				if (rand < sum) {
					#ifdef MYDEBUG_MCMCPQ_LOOP
						cout << "choosing this one: " << mit->toString() << endl;
					#endif
					found = true;
					break;
				}
			}
			assert(found);
			
			largest = *(mit); // the chosen largest in the set
			#ifdef MYDEBUG_MCMCPQ_LOOP
				cout << "chosen one is " << largest.toString() << endl;;
			#endif
			
			#ifdef MYDEBUG_MCMCPQ_EXTRA
				cout << "rand\t" << rand << endl;
			#endif
		}
		#ifdef MYDEBUG_MCMCPQ_EXTRA
				cout << "chosen\t" << largest.toString() << endl;
		#endif
		pq.erase(mit);// take the iterator to chosen largest out of the set

		/* accumulate the changes in scaled EMP sums that will result
		 * from this expansion */
		#ifdef LOGEMPS
			dotprecision changeCOPERR = largest.nodePtr->getSplitChangeEMPCOPERR(n);
			dotprecision changeAIC = largest.nodePtr->getSplitChangeEMPAIC();
			updateScaledEMPSumCOPERR(changeCOPERR);
			updateScaledEMPSumAIC(changeAIC);
		#endif
		
		// split the biggest one and divide up its data
		try {
			largest.nodePtr->nodeExpand();
			#ifdef MYDEBUG_MCMCPQ_LOOP
				cout << "should have now split it" << endl;
			#endif
		}
		catch (UnfulfillableRequest_Error& ure) {
			std::cerr << "Something has gone wrong in the splitting.  Error reported is:"
				<< std::endl;
			std::cerr << ure.what()	<< std::endl;
			std::cerr << "Process aborted with " << std::endl;
			canContinue = false;
		}
		assert(!largest.nodePtr->isLeaf());
		
		if (canContinue && !largest.nodePtr->isLeaf()) {
			size_t lcount = 0;
			size_t rcount = 0;
			
			
			bool volOkay = !( (minVol > 0.0) && (largest.nodePtr->nodeRealVolume() < 4*minVol) ); 
			
			{
				SPSnode* child 
							= largest.nodePtr->getLeftChild();
				lcount = child->getCounter();
				
				if ( lcount > 0) {
					/* isSplittableNode checks that node has at least 2 x min vol
					* ie that child will have >= minVol.*/
					if (volOkay && (child->isSplittableNode(minChildPoints))) {
						real m = measurer( child );
						#ifdef MYDEBUG_MCMCPQ_LOOP
							cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
						#endif
						pq.insert(NodePtrMeasurePair(child, m));
					}
				}
				else {
					loopEmptyBoxVolume += child->nodeRealVolume();
					
				}
			}
			{
				SPSnode* child 
							= largest.nodePtr->getRightChild();
				rcount = child->getCounter();
				
				if ( rcount > 0) {
					if (volOkay && (child->isSplittableNode(minChildPoints))) {
						real m = measurer( child );
						#ifdef MYDEBUG_MCMCPQ_LOOP
							cout << "inserting " << child->getNodeName() << " with measure " << m << " into the set" << endl;
						#endif
						pq.insert(NodePtrMeasurePair(child, m));
					}
				}
				else {
					loopEmptyBoxVolume += child->nodeRealVolume();
					
				}
			}
			/* we can easily get the change in likelihood using the
			 * left and right child counts*/
			deltaL = getSplitChangeLogLik(lcount, rcount);
		
		} // end can continue and largest not a leaf
		
	}
	else {
		canContinue = false;
		std::cerr << "Terminated splitting: no splittable nodes left"
				<< std::endl;
	}

	return canContinue;
}


// ----------- private inner classes

AdaptiveHistogram::PrivatePrioritySplitQueueEvaluator::
PrivatePrioritySplitQueueEvaluator(const 
					PrioritySplitQueueEvaluator& psqe)
	: measurer(psqe.getMeasurer()), critStop(psqe.getCritStop()), 
		maxLeaves(psqe.getMaxLeaves()),
		usingCritStop(psqe.getUsingCritStop()) {}
							
cxsc::real AdaptiveHistogram::PrivatePrioritySplitQueueEvaluator::
		measure(const SPSnode * const spn) const
{ return measurer(spn); }

bool AdaptiveHistogram::PrivatePrioritySplitQueueEvaluator::
		stopQueueQuery(const PriorityQueueT& pq,
						size_t numLeaves) const
{
	bool retValue = (numLeaves >= maxLeaves);
	
	if (!retValue && usingCritStop && !pq.empty()) {
		/* get the largest in the queue and compare its measure to critical value
		 * stop the queue if the measure on the largest is not larger
		 * than the critical value */
		retValue = (!(pq.rbegin()->measure > critStop)); 
		
	}	
	
	return retValue;
}

bool AdaptiveHistogram::PrivatePrioritySplitQueueEvaluator::
		stopQueueQuery(size_t numLeaves) const
{ return numLeaves >= maxLeaves; }				
		

// end inner class for evaluating queues

// implementations for private inner class for nodes and measurements pairs

AdaptiveHistogram::NodePtrMeasurePair::NodePtrMeasurePair(
	SPSnode *p, real m) : nodePtr(p), measure(m) {}

bool AdaptiveHistogram::NodePtrMeasurePair::operator<(
		const NodePtrMeasurePair& rhs) const
{
	return measure < rhs.measure;
}


std::string 
	AdaptiveHistogram::NodePtrMeasurePair::toString() const
{
	ostringstream oss;
	oss << (nodePtr->getNodeName()) << ":\t" << measure;
	return oss.str();
}



// ----------------------------- non member functions

//Output all boxes in AdaptiveHistogram adh
std::ostream & subpavings::operator<<(std::ostream &os, const subpavings::AdaptiveHistogram& adh)
{
    if (adh.hasSubPaving()) {
		adh.getSubPaving()->nodesAllOutput(os, 1);
		os << std::endl;
	}
    return os;
}

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap (subpavings::AdaptiveHistogram & a1, 
		subpavings::AdaptiveHistogram & a2) // throw ()
{
	a1.swap(a2);
}


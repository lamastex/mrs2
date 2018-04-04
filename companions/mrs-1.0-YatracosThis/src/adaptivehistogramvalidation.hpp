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

/*! \file      adaptivehistogramvalidation.hpp
\brief AdaptiveHistogramValidation declarations.
*/

#ifndef ___ADAPTIVEHISTVAL_HPP__
#define ___ADAPTIVEHISTVAL_HPP__

#include "spsvnode.hpp"
#include "splitdecisionobj.hpp"
#include "nodecompobjval.hpp"
#include "adaptivehistogram.hpp"
#include "realmappedspnode.hpp"
#include "piecewise_constant_function.hpp"
#include "pcfval.hpp"

#include <gsl/gsl_rng.h>        // to know about the gsl random number generator
#include <gsl/gsl_randist.h>    // we need these libraries to get the IAE for
										  // finite mixtures
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include "errorfunc.hpp"
#include "../examples/StatsSubPav/ExactInt/Int.h"
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"

struct RegHist;
struct FinMix;
class RSSample;

namespace subpavings {
	
//! Forward class declarations
class AdaptiveHistogram;
class AdaptiveHistogramValidation;
class AdaptiveHistogramVCollator;
class HistEvalObjVal;

using namespace subpavings;

/*! \brief A wrapper or manager for an SPSVnode aka StatsSubPavingVal in
conjunction with massive amounts of sample data.

Here sample data is multi-dimensional point-valued data in a cxsc::rvector
container.  The AdaptiveHistogramValidation class
manages \link subpavings::SPSVnode SPSVnode \endlink
objects (\link subpavings::StatsSubPaving StatsSubPavings \endlink)
for the purpose of creating adaptive histograms from sample data and also
manages the container for the sample data itself.

An SPSVnode (a pointer to an SPSVnode is aliased as StatsSubPaving) is a
binary tree representation of a regular subpaving which can be used for
processing statistical sample data.  SPSVnodes do not actually hold data,
they only need to know where the data they are associated with is stored.
The leaf nodes in the SPSVnode tree controlled by an AdaptiveHistogramValidation object
have a vector of iterators into the
\link subpavings::BigDataCollection BigDataCollection\endlink, a
dataCollection managed by that AdaptiveHistogramValidation object.

The AdaptiveHistogramValidation class uses the C-XSC library class rvector
for sample data points.  rvectors can have 1 or many dimensions.
*/

class AdaptiveHistogramValidation {
private:

    /*! \brief a constant for padding a box if it is tailor-made for data.

    The padding if the size of the root box is obtained from the min
    and max of the data that is fully fed in.
    */
    static const real padding;

    /*! \brief Pointer to the root node of the subpaving tree.

    An SPSVnode is a binary tree representation of a subpaving, designed for
    processing statistical data.
    */
    SPSVnode* rootVpaving;

    /*! \brief The root box used to form the subpaving tree.

    We may not need this, at present it gets passed to the SPSVnode
    constructor.
    */
    ivector rootBox;

    /*! \brief A container for all sample data passed to this.

    The sample that has come in thus far.
    */
    BigDataCollection dataCollection;

    /*! \brief Controls whether all available statistics are maintained in
    the rootPaving.  If set to false (default) only counts are maintained.
    */
    bool holdAllStats;

    /*! \brief A string showing the order of creation of the rootPaving.
    */
    std::string creationString;

    /*! \brief Private initialised constructor.
    Initialised  with pointer to subpaving and value for holdAllStats.
    */
    explicit AdaptiveHistogramValidation(SPSVnode * spn, bool as);

    /*! \brief Complete insertion of training and validation data from a vector container.   
   
    First checks if the box exists and makes it otherwise, then
    checks box dimensions against data dimensions if box already exists,
    and finally inserts the data. A boolean boolVal indicates if the data should be in the training or validation set.

    \param theData a reference to a container of rvector data.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether a log file is created for
    histogram creation process; set to TXT for logging to a txt file.
    \param boolVal a boolean to indicate if the data point should be in the validation set (true) or not (false).
    \pre A container theData containing at least one rvector.
    \pre this AdaptiveHistogram object may have an initialised SPSnode
    pointed to by rootPaving, but rootPaving can also be NULL.
    \post If the rootPaving was NULL when the function was entered, then
    rootPaving is pointed to a new SPSnode object whose root node has a box
    tailored to contain all the data read in.
    \post The data in theData has been put into the AdaptiveHistogram's
    dataCollection and also associated with  rootPaving's leaves via
    iterators to dataCollection.
    \return true if data successfully put into a dataCollection and
    associated with the rootPaving's leaves, false otherwise.
    */
    bool completeDataInsertionFromVec(const RVecData& theData,
                    const SplitDecisionObj& boolTest, LOGGING_LEVEL logging,
                    size_t holdOutCount, std::vector<size_t> & numNodes);

    /*! \brief Checks if we need to make root paving for the histogram object.

    Points rootPAving to a new SPSVnode if rootPaving is NULL, with the box
    of the root node matching the dimensions of the data as given in
    function argument and tailored to fit the data in theData.

    \param theData a reference to a container of rvector data.
    \param dim the dimensions of the data.
    \post if return value is true, rootPaving points to a new SPSVnode and a
    new ivector has been assigned to rootBox, and the rootBox ivector is
    the box of the root node of the SPSVnode to which rootPaving points.
    \return true if function needed to make a new SPSVnode, false if
    rootPaving already pointed to an SPS node.
    */
    bool haveMadePaving(const RVecData& theData, const size_t dim);

    /*! \brief Make a box to contain all the data.

    Used if a box has not already been provided. Makes a box tailored to
    contain all of the data.  So all the data has to be available for input.

    \param theData a reference to a container of rvector data
    \param dim the dimensions of the data
    \return an ivector same dimensions as the data and to fit all the data
    including an allowance for padding.
    */
    static ivector makeBox(const RVecData& theData, const size_t dim);

    /*! \brief Insert training and validation data from a container.

    Attempts to insert data from a container theData into this
    AdaptiveHistogram object's dataCollection and to associate the data
    with the leaves of the subpaving tree pointed to by this's rootPaving.
    Data in theData which falls outside the boundaries of the rootBox will
    not be inserted and a message will be printed to standard output. A boolean is used to indicate if the data point should be in the training or validation set.

    \param theData a reference to a container of rvector data.
    \param logging an enum controlling whether a log file is created for
    histogram creation process; set to TXT for logging to a txt file.
    \param boolTest is a reference to an object providing a function
    \param boolVal a boolean to indicate if the data point should be in the validation set (true) or not (false). 
    operator determining whether to split a node when a data point arrives.
    \pre A container theData containing at least one rvector.
    \pre this AdaptiveHistogramValidation object must have a rootPaving 
    pointing to an initialised SPSnode; this SPSnode may already have data associated with it.
    \post The data in theData which is within the boundaries of the rootBox
    has been put into the AdaptiveHistogram's dataCollection and also
    associated with  rootPaving's leaves via iterators to dataCollection.
    \return number of datapoints for which insertion has been attempted.
    */
    size_t insertDataFromContainer(const RVecData& theData,
                    const SplitDecisionObj& boolTest, LOGGING_LEVEL logging,
                    size_t holdOutCount, std::vector<size_t> & numNodes);
      
    /*! \brief Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.
    */
    void outputLogStart(const std::string& s) const;   
   
    /*! \brief Method to do checking for whether to split a node.

    Used in prioritySplit.

    Decides whether to split node based on checking volume and number of points
    that would result in child nodes.

    Node volume must be >=minVol to split
        and
    If there is a minChildPoints>0 specified, then
        either the node must have at least minChildPoints and all the points go
            to one of the children (the other getting none)
        or the smallest number of points which would go to the either of the
        prospective new children must be >= minChildPoints

    Thus in general the method will only return true if the given node satisfies
    both the minVol test and, if it were to be split, both children would have
    at least minChildPoints data points, but if all the data points would go
    to one child (none) to the other, this is considered to also satisfy the
    minChildPoints test.

    \param spn is a pointer to the target node.
    \param volChecking indicates whether volume is being checked
    \param minVol is the minimum volume allowed to be tested for.
    \param minChildPoints is the minimum number of points that there would be
    in the children if the node were to be split.
    \return true if has been a test conditions satisfied, false otherwise.
    */
    static bool checkNodeCountForSplit(const SPSVnode * const spn,
                bool volChecking, double minVol, size_t minChildPoints);

    public:

    /*! \brief Default constructor

    By default, only counts are maintained in subpaving this manages, rather
    than all available stats.
    */
    explicit AdaptiveHistogramValidation();

    /*! \brief Initialised constructor.

    Initialised with parameter controlling whether all available
    statistics be maintained in the SPSVnode tree managed by this
    AdaptiveHistogramValidation (true for all stats, false for counts only).
    */
    explicit AdaptiveHistogramValidation  (bool as);

    /*! \brief Initialised constructor.

    Initialised  with domain box.
    By default, only counts are maintained as stats in the in the SPSVnode tree
    managed by this AdaptiveHistogramValidation.

    Ideal constructor when the support domain of data is known a priori
    or has been transformed to a known domain but splitting criteria have
    not been determined a priori.
    */
    explicit AdaptiveHistogramValidation(ivector& v, bool as = false);

    /*! \brief  Copy constructor.
    */
    AdaptiveHistogramValidation(const AdaptiveHistogramValidation& other);

    /*! \brief Copy assignment operator.
    */
    AdaptiveHistogramValidation& operator=(const AdaptiveHistogramValidation& rhs);

    /*! \brief Overloaded addition operator.

    Makes a new histogram by adding this and rhs together.
    The subpaving will be the union of the subpavings of this and rhs.
    Data is reinserted into the new histogram so that the counts in each box
    of the subpaving are exactly right for that subpaving.  This is in contrast
    to the way that the AdaptiveHistogramValidationCollator works (divides a count evenly
    in two when apportioning it between bisected boxes).
    holdAllStats will be set to the logical and of the values for this and rhs.
    */
    AdaptiveHistogramValidation operator+(const AdaptiveHistogramValidation& rhs);


    //! Destructor
    ~AdaptiveHistogramValidation();

    //src_trunk_0701
    /*! \brief Get whether this has a subpaving to manage.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return true if this has a subpaving to manage.
	false otherwise.*/
    bool hasSubPaving() const;
    
      //src_trunk_0701
    /*! \brief Return the label for this.
	*/
    int getLabel() const;

    /*! \brief Return a pointer to the SPSVnode this manages.
    */
    SPSVnode* getSubPaving() const;

    /*! \brief Gets count in the rootpaving in the root paving.
    */
    size_t getRootVcounter() const;

    /*! \brief Gets number of leaf nodes in the root paving.
    */
    size_t getRootLeaves() const;

    /*! \brief Gets the sum of leaf count over volume in root paving.
    */
    real getRootSumLeafCountOverVol() const;
       
    /*! \brief get the value of the minimum volume for a splittable node.

    Minimum volume = minVolB * (log n) ^2/n where n is points in histogram.
    Minimum volume is used in COPERR or AIC priority queue splitting to
    limit which nodes can be split.
    \param minVolB the multiplier applied to log n) ^2/n to find the
    minimum allowed node volume at which a node can be split (children will
    half the volume of the parent node).
    */
    double getMinVol(double minVolB) const;

    /*! \brief get the value of holdAllStats field.

    This determines whether the histrogram's rootPaving will maintain all
    available stats (true) or just the counts (false).
    */
    bool getHoldAllStats() const;

    /*! Get a vector of the leaf node levels.

    Root is level 0, next level down is 1, etc.

    \return a vector of leaf levels, left to right order,
    */
    IntVec getLeafLevels() const;


    /*! Get a string of the leaf node levels.

    Root is level 0, next level down is 1, etc.
    Example return string "3,3,2,1"

    \return a comma separated string of leaf levels, left to right order
    */
    std::string getLeafLevelsString() const;


    /*! Get a vector of the leaf node counts.

    \return a vector of leaf counts, left to right order,
    */
    Size_tVec getLeafCounts() const;

    /*! \brief Append current state of histogram to a txt log file.

    Format is a tab-delimited file of numeric data.
    Output includes node contributions to unscaled EMP under COPERR and AIC
    and the changes in EMP that would result from splitting the node.

    \param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process
    */
    void outputLog(const std::string& s, const int i) const;
    
    /** @name Insert rvectors in a txt file into AdaptiveHistogram object.

    A group of overloaded functions which read in lines of data
    representing rvectors from a txt file.  The dimensions of the
    rvector are deduced from the input format and all the data then
    expected to have the same dimension.  Any data not matching the
    expected dimensions, based on assessing the first valid line found,
    will be rejected.  Expects one line per rvector with the elements
    separated by white space (space or tabs), with no non-numeric
    characters.  Carries out some basic data checking through
    checkString().  Input lines which do not pass are printed to standard
    output with an error message but the entire file will continue to be
    processed and valid lines converted to rvectors which are stored in
    theData.  Conversion to rvectors is via the cxsc::operator<< which
    allows an rvector to be constructed from a stream.

    Can read 1-d rvector data from doubles or floats but insists on a decimal
    point in each number (otherwise the number is rejected).

    For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    3-dimensional rvector, a string "-30.0006" will be read as a
    1-dimensional rvector and a string "30" will be rejected.

    \param s the name of the txt file to read data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
	\param headerlines is number of headerlines to skip before reading 
	data.  Defaults to 0.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)
    \pre A file with filename s in the same directory as main() or with the
    filename incorporating directory location.
    \return true if at least some data inserted, false if txt file cannot
    be read, is empty, or contains no valid data, or if something failed
    after data is read and when it is being inserted.
    */
    //@{
    /** All rvectors are associated with the root paving, no spliting. */
    bool insertRvectorsFromTxt(const std::string& s,
										   std::vector<size_t>& numNodes,
										  const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertRvectorsFromTxt(s, numNodes, sn, headerlines, logging);
    }

    /** Adaptive splitting with each data point inserted. */
    bool insertRvectorsFromTxt(const std::string& s,
	   								   std::vector<size_t>& numNodes,
										  const SplitDecisionObj& boolTest,
										  const std::size_t headerlines = 0,
										  LOGGING_LEVEL logging=NOLOG);
    //@}

 
    /** @name Insert all rvectors from a container of rvectors.
    */
    //@{
   /*! Adaptive splitting with each data point inserted. */
    bool insertFromRVec(const RVecData& rvec, const SplitDecisionObj& boolTest,
                             LOGGING_LEVEL logging=NOLOG); 
                                                             
    /*! All rvectors are associated with the root paving for hold out estimation, no splitting. */
  	bool insertFromRVecForHoldOut(const RVecData& rvec, 
											const SplitDecisionObj& boolTest,
											int holdOutCount,
											LOGGING_LEVEL logging=NOLOG);

  /*! All RSSample are associated with the root paving for hold out estimation, no splitting. */
  	bool insertFromRSSampleForHoldOut(const RSSample& rss, int label,   
											const SplitDecisionObj& boolTest,
											int holdOutCount,
											LOGGING_LEVEL logging);

    //@}

   /** @name Data splitting method to obtain the "best" estimate.

    Implementation of the data splitting method based on Devroye and Lugosi, 2001, p. 98.
    
    This takes a histogram and progressively split using a priority
    queue to determine which node to split first. The current histogram is 
	 added into an \link AdaptiveHistogramCollator AdaptiveHistogramCollator 
	 \endlink object. The Yatracos set as defined in Devroye and Lugosi, 2001 is 
	 obtained and stored in a list of pointers to sets of 
	 \link subpavings::CollatorSPnode \endlink. The corresponding delta value is 
	 obtained by taking the difference of the empirical measure of the training 
	 data and the empirical measure of the validation data.  
    
    Splitting continues until some criteria applying either to individual nodes 
	 or to the histogram as a whole is satisfied, or there are no more 
	 splittable nodes, or if some criteria in relation to the delta value is 
	 satisfied.

    Nodes are not considered to be splittable if they satisfy two criteria:
    First, their volume is greater than the minimum volume specified for the
    histogram as a whole, minVolB.
    Second if both prospective children would have at least the parameter
    minChildPoints data points associated with them.

    If more than one node is equally 'large', on the basis of the node
    comparison compTest used, then a random choice is made between all equally
    large nodes to find the node which will be split.

    The random number generator used for random selection between equally
    'large' nodes uses a default seed to ensure that results can be
    replicated.  If you are looking at distributions of results across
    mulitple histograms, supply the random number generator to the priority
    queue to ensure that each histogram will make different random choices.

    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to determine
    when to stop splitting.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \param rgsl is a pointer to a gls random number generator.
    \param tol is the tolerance for the stopping criteria.
    \param distr is an integer that indicates what distribution is being used to calculate the integrated absolute error (IAE) for simulation purposes. If 0, there will be no calculations for the IAE.
     \return coll an \link AdaptiveHistogramCollator AdaptiveHistogramCollator \endlink object.
    */
    //@{
    /** 2D shapes: minVolB and minChildPoints supplied but no random number generator.*/
    AdaptiveHistogramVCollator prioritySplitAndEstimate(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, int distr, int method, size_t hist,
									size_t maxLeafNodes);							
									
    /** 2D shapes: With random number generator. All other parameters supplied.*/
    AdaptiveHistogramVCollator prioritySplitAndEstimate(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, gsl_rng * rgsl, 
									bool stopCrit, int distr, int method, size_t hist,
									size_t maxLeafNodes);
	
	//@}
	
	/** @name Data splitting method to obtain the "best" estimate for uniform 
	          mixtures.
	 
   The implementation is the same as above except that the IAE calculations are
	for uniform mixtures.
	 
   \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
   \param he is an instance of a class which provides a function to determine
    when to stop splitting.
   \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
   \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
   \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
   \param rgsl is a pointer to a gls random number generator.
   \param tol is the tolerance for the stopping criteria.
   \param myPart is a reference to to StatsSubPaving that is split into a
	uniform mixture. 
   \return coll an \link AdaptiveHistogramCollator AdaptiveHistogramCollator 
	  \endlink object.
    */
    //@{
    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									AdaptiveHistogram& myPart, double weight,
									std::vector<int> holesLoc,
									int method, size_t hist,
									size_t maxLeafNodes, int maxCheck,
									int& minTheta);
									//AdaptiveHistogramValidation& optHist);

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									AdaptiveHistogram& myPart, double weight,
									std::vector<int> holesLoc,
									int method, size_t hist,
									size_t maxLeafNodes, int maxCheck,
									int& minTheta);
									//AdaptiveHistogramValidation& optHist);
									
	bool prioritySplitAndEstimateWithSwitch(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									AdaptiveHistogram& myPart, double weight,
									std::vector<int> holesLoc, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck,
									AdaptiveHistogramValidation& opthist);
    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimateWithSwitch(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									AdaptiveHistogram& myPart, double weight,
									std::vector<int> holesLoc, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck, 
									AdaptiveHistogramValidation& opthist);
									
   //@}
   
   /** @name Data splitting method - no MDE estimation being done
	  
   \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
   \param he is an instance of a class which provides a function to determine
    when to stop splitting.
   \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
   \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
   \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
   \param rgsl is a pointer to a gls random number generator.
   \param tol is the tolerance for the stopping criteria.
    */
    //@{
    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplit(const NodeCompObjVal& compTest,
					   const HistEvalObjVal& he, LOGGING_LEVEL logging,
                       size_t minChildPoints, double minVolB, 
					   size_t maxLeafNodes);

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplit(const NodeCompObjVal& compTest, 
                       const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                       size_t minChildPoints, double minVolB, 
						gsl_rng * rgsl, size_t maxLeafNodes);
									
   //@}
   
   /** @name Data splitting method to obtain the "best" estimate for finite
	          mixtures.
	 
   The implementation is the same as above except that the IAE calculations are
	for finite mixtures.
	 
   \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
   \param he is an instance of a class which provides a function to determine
    when to stop splitting.
   \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
   \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
   \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
   \param rgsl is a pointer to a gls random number generator.
   \param tol is the tolerance for the stopping criteria.
   \param myPart is a reference to to StatsSubPaving that is split into a
	uniform mixture. 
   \return coll an \link AdaptiveHistogramCollator AdaptiveHistogramCollator 
	  \endlink object.
    */
    //@{
    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									FinMix& mixt, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck, double tol, int deg,
									int& minTheta);
									//AdaptiveHistogramValidation& opthist);

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									FinMix& mixt, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck, double tol, int deg,
									int& minTheta);
									//AdaptiveHistogramValidation& opthist);
									
	 bool prioritySplitAndEstimateWithSwitch(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
                           size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									FinMix& mixt, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck, double tol, int deg,
									AdaptiveHistogramValidation& opthist);
    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimateWithSwitch(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									FinMix& mixt, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck, double tol, int deg,
									AdaptiveHistogramValidation& opthist);
   //@}
	
	/** @name Data splitting method to obtain the "best" estimate for uniform 
	          mixtures.
	 
   The implementation is the same as above except that the IAE calculations are
	for uniform mixtures.
	 
   \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
   \param he is an instance of a class which provides a function to determine
    when to stop splitting.
   \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
   \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
   \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
   \param rgsl is a pointer to a gls random number generator.
   \param tol is the tolerance for the stopping criteria.
   \param myPart is a reference to to StatsSubPaving that is split into a
	uniform mixture. 
   \return coll an \link AdaptiveHistogramCollator AdaptiveHistogramCollator 
	  \endlink object.
    */
    //@{
    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
									size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									RealMappedSPnode& nodeEst, int method, 
									size_t maxLeafNodes, int maxCheck,
									AdaptiveHistogramValidation& optHist);
	//20160830								
    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimate(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									PiecewiseConstantFunction& nodeEst, 
									int method, 
									size_t maxLeafNodes, int maxCheck,
									int& minTheta, vector<int> sequence,
									vector<double> & vecMaxDelta, vector<real> & vecIAE);
									
	 bool prioritySplitAndEstimate(const NodeCompObjVal& compTest,
									const HistEvalObjVal& he, LOGGING_LEVEL logging,
									size_t minChildPoints, double minVolB, 
									bool stopCrit, 
									PiecewiseConstantFunction& nodeEst, 
									int method, 
									size_t maxLeafNodes, int maxCheck,
									int& minTheta, vector<int> sequence,
									vector<double> & vecMaxDelta, vector<real> & vecIAE);
    
    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitAndEstimateWithSwitch(const NodeCompObjVal& compTest, 
                           const HistEvalObjVal& he, LOGGING_LEVEL logging, 
                           size_t minChildPoints, double minVolB, 
									gsl_rng * rgsl, bool stopCrit,
									RealMappedSPnode& nodeEst, int method, size_t hist,
									size_t maxLeafNodes, int maxCheck,
									AdaptiveHistogramValidation& optHist);
   
   //@}
	

			
   /*! \brief Merge a multileaf histogram up to just root box.

    No prioritistion, just brute force
    */
    bool mergeUp();

    /*! \brief Split a histogram  to a specified shape.

    \param instruction specifies the required shape, eg "3, 3. 2, 1"
    */
    bool splitToShape(std::string instruction);

    /*! \brief Make a .dot graph file from histogram structure.

    Makes a simple .dot graph from the histogram using node names and the
    .png image for this graph.

    \pre a constructed histogram
    \post a .dot file and a .png in the same directory as the program creating
    it was run in.
    */
    bool outputGraphDot() const;


    /*! \brief Output the subpaving managed by this to a txt file.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then the description of the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param s the name of the txt file to send output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputToTxtTabs(const std::string& s, bool confirm = false) const;

    /*! \brief Output details of full sample (from root) to txt tile.

    Format is a mixture of alpha and  numeric data.

    \param s the name of the txt file to send output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputRootToTxt(const std::string& s, bool confirm = false) const;
    
   /** @name Get the IAE of a distribution
	 
	Get the integrated absolute error of the specified distribution.
	\param distr is an integer that indicates which distribution is used.
					1: bivariate gaussian distribution (not implemented yet)
					2: Levy 2D distribution (not implemented yet)
					3: Rosenbrock 2D distribution (not implemented yet) 
	\return the integrated absolute error for this realization
	*/
	//@{
	/** Get the IAE of a distribution. */
	cxsc::real getIAE(int distr);

	/** Get the IAE for a finite gaussian mixture distribution. */ 
	cxsc::real getFinMixIAE(FinMix& mixt); 
	
	/** Get the IAE for a finite gaussian mixture distribution. */ 
	cxsc::interval getFinMixIntervalIAE(FinMix& mixt, double tol, int deg, bool full);
	
	
	
	/** Get the IAE of a bivariate gaussian/Levy 2D/Rosen 2D distribution. */
	cxsc::real get2DIAE(taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval));
	
	/** Get the IAE for a uniform (mixture) distribution. */ 
	cxsc::real getUnifIAE(AdaptiveHistogram& myPart, double weight, 
									vector<int> holesLoc, bool full);
	
	/** Get the IAE for mapped functions. */
	//cxsc::real getMappedFunctionIAE(RealMappedSPnode& nodeEst, bool full);

	//@}
	

}; // end of AdaptiveHistogramValidation class declarations

} // end of namespace subpavings

// ----------  declarations of non-member functions ----------------------

/*! \brief Output the contents of an AdaptiveHistogramValidation object.

Verbose output for an AdaptiveHistogramValidation object, including all boxes
(not just leaves), data, and summary statistics.
*/
std::ostream & operator<<(std::ostream &os, const subpavings::AdaptiveHistogramValidation& adh);

size_t checkNumValley(std::vector<double> vecMaxDelta, 
							std::vector<int>& valleyHistPos, bool& plateau, int& currentSmallest);
							
void getCurrentYatClass(std::vector< std::set<CollatorSPVnode*, std::less < CollatorSPVnode* > > >& vecRowYatSet,
std::vector< std::set<CollatorSPVnode*, std::less < CollatorSPVnode* > > >& vecColYatSet, 
std::list < std::set<CollatorSPVnode*, std::less < CollatorSPVnode* > > >& listYatSet);

/** Get the true delta for a finite gaussian mixture distribution. */ 
	cxsc::interval getFinMixIntervalTrueDelta(FinMix& mixt, double tol, int deg, 
		std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet);


/** Get the true delta for a finite gaussian mixture distribution. */ 
	cxsc::real getUnifTrueDelta(AdaptiveHistogram& myPart, double weight, vector<int> holesLoc, 
		std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet);

/** Get the true delta for mapped functions. */
	cxsc::real getMappedFunctionTrueDelta(PiecewiseConstantFunction& nodeEstHist, 
				std::set<CollatorSPVnode*, less < CollatorSPVnode* > >& YatSet);



#endif



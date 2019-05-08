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

/*! \file      adaptivehistogram.hpp
\brief AdaptiveHistogram declarations.
*/

#ifndef ___ADAPTIVEHIST_HPP__
#define ___ADAPTIVEHIST_HPP__

#include "spsnode.hpp"
#include "splitdecisionobj.hpp"
#include "nodecompobj.hpp"
#include <gsl/gsl_rng.h>        // to know about the gsl random number generator
#include "errorfunc.hpp"
#include "../examples/StatsSubPav/ExactInt/Int.h"
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"
#include "realmappedspnode.hpp"
#include "piecewise_constant_function.hpp"

//using namespace subpavings;
struct RegHist;
struct FinMix;

namespace subpavings {

//! Forward class declarations
class AdaptiveHistogramCollator;
class LogMCMCPrior;
class MCMCProposal;
class HistEvalObj;
class PenObj;

/*! \brief A wrapper or manager for an SPSnode aka StatsSubPaving in
conjunction with massive amounts of sample data.

Here sample data is multi-dimensional point-valued data in a cxsc::rvector
container.  The AdaptiveHistogram class
manages \link subpavings::SPSnode SPSnode \endlink
objects (\link subpavings::StatsSubPaving StatsSubPavings \endlink)
for the purpose of creating adaptive histograms from sample data and also
manages the container for the sample data itself.

An SPSnode (a pointer to an SPSnode is aliased as StatsSubPaving) is a
binary tree representation of a regular subpaving which can be used for
processing statistical sample data.  SPSnodes do not actually hold data,
they only need to know where the data they are associated with is stored.
The leaf nodes in the SPSnode tree controlled by an AdaptiveHistogram object
have a vector of iterators into the
\link subpavings::BigDataCollection BigDataCollection\endlink, a
dataCollection managed by that AdaptiveHistogram object.

The AdaptiveHistogram class uses the C-XSC library class rvector
for sample data points.  rvectors can have 1 or many dimensions.
*/

class AdaptiveHistogram {
private:

    /*! \brief a constant for padding a box if it is tailor-made for data.

    The padding if the size of the root box is obtained from the min
    and max of the data that is fully fed in.
    */
    static const real padding;

    /*! \brief Pointer to the root node of the subpaving tree.

    An SPSnode is a binary tree representation of a subpaving, designed for
    processing statistical data.
    */
    SPSnode* rootPaving;

    /*! \brief The root box used to form the subpaving tree.

    We may not need this, at present it gets passed to the SPSnode
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

    /*! \brief A value for the unscaled EMP part of COPERR score.

    pi for a histogram is c * e(-1/t * energy)
    energy is EMP + PEN

    Under COPERR scoring, EMP is
    -1/n^2 x sum over leaves of (counts in leaf squared / volume of leaf)
    where n is the total number of data points in the histogram

    Scaling means multiplying this unscaled Emp sum term by 1/n^2

    scaledEMPSumCOPERR will default to zero until set by some method of
    the object.

    Dotprecision accumulator used to mitigate rounding errors in calculations.
    */
    mutable dotprecision scaledEMPSumCOPERR;

    /*! \brief A value for the unscaled EMP part of AIC score.

    pi for a histogram is c * e(-1/t * energy)
    energy is EMP + PEN

    Under AIC, EMP is
    - 1 x sum over leaves of
    (counts in leaf x -ln(count in leaf / (n x vol of leaf)))
    where n is the total number of data points in the histogram

    scaledEMPSumAIC will default to zero until set by some method of
    the object.

    Dotprecision accumulator used to mitigate rounding errors in calculations.
    */
    mutable dotprecision scaledEMPSumAIC;



    /*! \brief Private initialised constructor.

    Initialised  with pointer to subpaving and value for holdAllStats.
    */
    explicit AdaptiveHistogram(SPSnode * spn, bool as);

    /*! \brief Complete insertion of data from a vector container.

    First checks if the box exists and makes it otherwise, then
    checks box dimensions against data dimensions if box already exists,
    and finally inserts the data.

    \param theData a reference to a container of rvector data.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether a log file is created for
    histogram creation process; set to TXT for logging to a txt file.
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
                    const SplitDecisionObj& boolTest, LOGGING_LEVEL logging);

    /*! \brief Checks if we need to make root paving for the histogram object.

    Points rootPAving to a new SPSnode if rootPaving is NULL, with the box
    of the root node matching the dimensions of the data as given in
    function argument and tailored to fit the data in theData.

    \param theData a reference to a container of rvector data.
    \param dim the dimensions of the data.
    \post if return value is true, rootPaving points to a new SPSnode and a
    new ivector has been assigned to rootBox, and the rootBox ivector is
    the box of the root node of the SPSnode to which rootPaving points.
    \return true if function needed to make a new SPSnode, false if
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

    /*! \brief Insert data from a container.

    Attempts to insert data from a container theData into this
    AdaptiveHistogram object's dataCollection and to associate the data
    with the leaves of the subpaving tree pointed to by this's rootPaving.
    Data in theData which falls outside the boundaries of the rootBox will
    not be inserted and a message will be printed to standard output.

    \param theData a reference to a container of rvector data.
    \param logging an enum controlling whether a log file is created for
    histogram creation process; set to TXT for logging to a txt file.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \pre A container theData containing at least one rvector.
    \pre this AdaptiveHistogram object must have a rootPaving pointing to
    an initialised SPSnode; this SPSnode may already have data associated
    with it.
    \post The data in theData which is within the boundaries of the rootBox
    has been put into the AdaptiveHistogram's dataCollection and also
    associated with  rootPaving's leaves via iterators to dataCollection.
    \return number of datapoints for which insertion has been attempted.
    */
    size_t insertDataFromContainer(const RVecData& theData,
                    const SplitDecisionObj& boolTest, LOGGING_LEVEL logging);

    /*! \brief Recalculate the scaled EMP part of COPERR score.

    Uses whole histogram
    */
    void recalcScaledEMPSumCOPERR() const;

    /*! \brief Recalculate the unscaled EMP part of AIC score.

    Uses whole histogram.
    */
    void recalcScaledEMPSumAIC() const;

    /*! \brief Update the scaled EMP part COPERR score given change.
    */
    void updateScaledEMPSumCOPERR(dotprecision change) const;

    /*! \brief Update the the scaled EMP part AIC score given change.
    */
    void updateScaledEMPSumAIC(dotprecision change) const;


    /*! \brief Add COPERR EMP score to log file.

    Tab delimited
    \param s the name of the txt file to send output to.
    */
    void outputLogEMPCOPERR(const std::string& s) const;

    /*! \brief Add AIC EMP score to log file.

    Tab delimited
    \param s the name of the txt file to send output to.
    */
    void outputLogEMPAIC(const std::string& s) const;

    /*! \brief Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.
    */
    void outputLogStart(const std::string& s) const;

    /*! \brief Send a collection of changes in MCMC probabilities to log file.

    The probabilities being logged are the changes in the components of the
    probability of acceptance of a proposed state under MCMC

    \param s the name of the txt file to send output to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param deltaL is the log change in likelihood.
    \param deltaP is the log change in prior probability.
    \param deltaQ is the log change in transition probability.
    \param deltaPi is the log change in posterior probability.
    \param randChange is the draw from the Uniform (0,1) distribution which is
    compared to the changes in combined posterior and transition probabilities
    to determine whether the proposed change takes place
    */
    static void logMCMCDeltas(std::string s, int i,
                        real deltaL, real deltaP, real deltaQ, real deltaPi,
                        double randChange);


    /*! \brief Put header in a log file for MCMC.

    \param s is the name of the file to log the final state to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param proposal a reference to a proposal function object.
    \param logPrior a reference to a log prior function object.
    */
    void MCMCStartLogFile(std::string s, int i, const MCMCProposal& proposal,
                                    const LogMCMCPrior& logPrior);

    /*! \brief Output the state of this histogram as an MCMC sample.

    \param i is a number indicating the state number in a chain or series
    of changes in state.
    */
    void outputMCMCStateSample(int i);

    /*! \brief Capture the final state of this histogram after MCMC.

    \param s is the name of the file to log the final state to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    */
    void MCMCLogFinalState(std::string s, int i);


    /*! \brief Finds the node to target for change in MCMC on SPSnode trees.

    This method changes the value of the bool haveNode and returns an iterator
    into a container of SPSnode pointers.  If, after the method has completed,
    haveNode is positive, then the iterator points to a pointer to a node
    proposed for change (split or merge) under MCMC on the SPSnode tree.

    Proposals are made by selecting a splittable leaf or cherry at random.
    If a cherry is chosen the proposed change in state is to merge the two
    leaf children of that cherry; if a splittable leaf is chosen the proposal
    is to split (bisect) the leaf.

    \pre haveNode is false.
    \param proposal is a reference to a proposal distribution object.
    \param nodes is a reference to a container of pointers to the splittable
    leaf and cherry nodes of the SPSnode tree managed by this AdaptiveHistogram,
    which will be updated if the method results in change in tree state.
    \param numLeaves is the number of splittable leaves in the SPSnode tree.
    in tree state.
    \param numCherries is the number of cherries in the SPSnode tree.
    in tree state.
    \param rgsl is a pointer to a uniform random number generator.
    \param haveNode is reference to a boolean which can be changes in the
    function, which indicates whether a proposal node is found.
    \return an iterator to a container of pointers to SPSnodes.
    \post haveNode will have been changed to true if a node has been found
    for the proposal using the probabilities under proposal.
    */
    static SPSnodeListItr proposeChangeMCMCState (const MCMCProposal& proposal,
                            SPSnodeList& nodes,
                            size_t numLeaves, size_t numCherries,
                            gsl_rng* rgsl, bool& haveNode);

    /*! \brief Determines whether to split a node to get a new MCMC state.

    This method probabilistically accepts or rejects a single-step change
    in the histogram state represented by this Adaptive Histogram resulting
    from a single split of a target splittable leaf node.

    Q(m' | m) is the transition probability from state m to state m'.

    This proposal is accepted if u drawn from Uniform(0,1) is such that
    u < (posterior probability of m' x Q(m | m'))/(posterior probability
    of m x Q(m' | m).  In this implementation, natural logs are used to
    simplify calculation, ie a proposal is accepted if
    log(u) < log [(posterior probability of m' x Q(m | m'))/(posterior
    probability of m x Q(m' | m)].

    \param target is a pointer to the target node proposed for splitting.
    \param proposal is a reference to a proposal distribution object.
    \param logPrior is a reference to a log prior object.
    \param rgsl is a pointer to a uniform random number generator.
    \param numLeaves is the number of splittable leaves in the SPSnode tree.
    \param numCherries is the number of cherries inthe SPSnode tree.
    \param minPoints is the minimum number of points allowed in a box.
    \param logging an enum controlling whether decision making output is
    sent to the log file.
    \param s is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain
    \return true if the proposal is accepted (the target will be split),
    false otherwise.
    */
    bool decisionMCMCSplit(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        size_t numLeaves, size_t numCherries, size_t minPoints,
                        LOGGING_LEVEL logging, const std::string& s, int i) const;
                        
    /*! \brief Determines whether to merge a node to get a new MCMC state.

    This method probabilistically accepts or rejects a single-step change n the
    histogram state represented by this Adaptive Histogram resulting from
    merging the leaf children of a target cherry node back into the target.

    Q(m' | m) is the transition probability from state m to state m'.

    This proposal is accepted if u drawn from Uniform(0,1) is such that
    u < (posterior probability of m' x Q(m | m'))/(posterior probability
    of m x Q(m' | m).  In this implementation, natural logs are used to
    simplify calculation, ie a proposal is accepted if
    log(u) < log [(posterior probability of m' x Q(m | m'))/(posterior
    probability of m x Q(m' | m)].

    \param target is a pointer to the target node proposed for merging.
    \param proposal is a reference to a proposal distribution object.
    \param logPrior is a reference to a log prior object.
    \param rgsl is a pointer to a uniform random number generator.
    \param numLeaves is the number of splittable leaf nodes in the SPSnode tree.
    \param numCherries is the number of cherries in the SPSnode tree.
    \param minPoints is the minimum number of points allowed in a box.
    \param logging an enum controlling whether decision making output is
    sent to the log file.
    \param s is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain
    \return true if the proposal is accepted (the target will be merged),
    false otherwise.
    */
    bool decisionMCMCMerge(SPSnode* target,
                        const MCMCProposal& proposal,
                        const LogMCMCPrior& logPrior, gsl_rng* rgsl,
                        size_t numLeaves, size_t numCherries, size_t minPoints,
                        LOGGING_LEVEL logging, const std::string& s, int i) const;


    /*! \brief Changes the state of this Adaptive Histogram by splitting a node.

    This method carries out a move to a new state in the histogram MCMC state
    chain by splitting the target splittable leaf node.

    minPoints restricts the pointers to new leaves resulting from the split
    which can be put into the container of pointers to splittable leaves and
    cherry nodes.  A pointer to a new leaf can only be put into the container
    if splitting that leaf would result in the children having less than
    minPoints data points associated with them.

    \pre an SPSnode tree representing some histogram state, a container of
    pointers to the splittable leaf and cherry nodes of the tree in its
    current state, ordered with the leaves first followed by the cherries,
    the number of splittable leaf nodes and the number of cherry nodes.
    \post the SPSnode tree changed as a result of the split, the container of
    splittable leaf and cherry nodes updated for any change in state, and the
    number of splittable leaves and cherries updated similarly.

    \param target is a pointer to the target node proposed for splitting.
    \param nodes is a reference to a container of pointers to the leaf and
    cherry nodes of the SPSnode tree managed by this AdaptiveHistogram, which
    will be updated by the method.
    \param numLeaves is a reference to a variable storing the number of
    splittable leaves in the SPSnode tree, which will be updated by the method.
    \param numCherries is a reference to a variable storing the number of
    cherries in the SPSnode tree, which will be updated by the method.
    \param minPoints is the minimum number of points allowed in a box.
    \return true if has been a successful change in state, false otherwise.
    */
    bool changeStateForSplit(SPSnode* target, SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        size_t minPoints);

    /*! \brief Changes the state of this Adaptive Histogram by merging cherry.

    This method carries out a move to a new state in the histogram MCMC state
    chain by merging the target cherry node.

    \pre an SPSnode tree representing some histogram state, a container of
    pointers to the splittable leaf and cherry nodes of the tree in its
    current state, ordered with the leaves first followed by the cherries,
    the number of splittable leaf nodes and the number of cherry nodes.
    \post the SPSnode tree changed as a result of the merge, the container of
    splittable leaf and cherry nodes updated for any change in state, and the
    number of splittable leaves and cherries updated similarly.

    \param target is a pointer to the target node proposed for merging.
    \param nodes is a reference to a container of pointers to the splittable
    leaf and cherry nodes of the SPSnode tree managed by this AdaptiveHistogram,
    which will be updated by the method.
    \param numLeaves is a reference to a variable storing the number of
    splittable  leaf nodes in the SPSnode tree, which will be updated
    by the method.
    \param numCherries is a reference to a variable storing the number of
    cherries in the SPSnode tree, which will be updated by the method.
    \param minPoints is the minimum number of points allowed in a box, which is
    needed to be able to tell which, if any, of the merged nodes children need
    to be taken out of the nodes container on the merge.
    \return true if has been a successful change in state, false otherwise.
    */
    bool changeStateForMerge(SPSnode* target, SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        size_t minPoints);
                        
    /*! \brief Method to do checking for whether a node is splittable.

    Decides whether a node is splittable based on checking volume 
	and number of points that would result in child nodes on split.

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
	
	If the node has already been split, the test will use the actual numbers
	of points in the children; if the node is a leaf (ie not split) then
	the test will consider the number of points that would go to the 
	each child if it were to be split.

    \param spn is a pointer to the target node.
    \param volChecking indicates whether volume is being checked
    \param minVol is the minimum volume allowed to be tested for.
    \param minChildPoints is the minimum number of points that there would be
    in the children if the node were to be split.
    \return true if has been a test conditions satisfied, false otherwise.
    */
    static bool checkNodeCountForSplit(const SPSnode * const spn,
                bool volChecking, double minVol, size_t minChildPoints);




    public:

    /*! \brief Default constructor

    By default, only counts are maintained in subpaving this manages, rather
    than all available stats.
    */
    explicit AdaptiveHistogram();

    /*! \brief Initialised constructor.

    Initialised with parameter controlling whether all available
    statistics be maintained in the SPSnode tree managed by this
    AdaptiveHistogram (true for all stats, false for counts only).
    */
    explicit AdaptiveHistogram  (bool as);

    /*! \brief Initialised constructor.

    Initialised  with domain box.
    By default, only counts are maintained as stats in the in the SPSnode tree
    managed by this AdaptiveHistogram.

    Ideal constructor when the support domain of data is known a priori
    or has been transformed to a known domain but splitting criteria have
    not been determined a priori.
    */
    explicit AdaptiveHistogram(ivector& v, bool as = false);

    /*! \brief  Copy constructor.
    */
    AdaptiveHistogram(const AdaptiveHistogram& other);

    /*! \brief Copy assignment operator.
    */
    AdaptiveHistogram& operator=(const AdaptiveHistogram& rhs);

    /*! \brief Overloaded addition operator.

    Makes a new histogram by adding this and rhs together.
    The subpaving will be the union of the subpavings of this and rhs.
    Data is reinserted into the new histogram so that the counts in each box
    of the subpaving are exactly right for that subpaving.  This is in contrast
    to the way that the AdaptiveHistogramCollator works (divides a count evenly
    in two when apportioning it between bisected boxes).
    holdAllStats will be set to the logical and of the values for this and rhs.
    */
    AdaptiveHistogram operator+(const AdaptiveHistogram& rhs);


    //! Destructor
    ~AdaptiveHistogram();



    /*! \brief Return a pointer to the SPSnode this manages.
    */
    SPSnode* getSubPaving() const;
    
    //src_trunk_0701
    /*! \brief Return the label for this.
	*/
    int getLabel() const;


    /*! \brief Gets the mean from the root box of the paving this manages.
    */
    rvector getRootPavingMean() const;


    /*! \brief Gets variance covariance vector from root box of rootpaving.
    */
    RealVec getRootPavingVarCovar() const;

    /*! \brief Gets count in the root paving.
    */
    size_t getRootCounter() const;

    /*! \brief Gets number of leaf nodes in the root paving.
    */
    size_t getRootLeaves() const;

    /*! \brief Gets the sum of leaf count over volume in root paving.
    */
    real getRootSumLeafCountOverVol() const;


    /*! \brief get the EMP part of the COPERR score.

    */
    real getEMPScoreCOPERR() const;

    /*! \brief get the EMP part of the AIC score.

    */
    real getEMPScoreAIC() const;


    /*! \brief get the COPERR score.

    \param pen is the penalty function object
    \param verbose option: true for extra console output (default false)
    */
    real getScoreCOPERR(const PenObj& pen, bool verbose = false) const;

    /*! \brief get AIC score.

    \param pen is the penalty function object
    \param verbose option: true for extra console output (default false)
    */
    real getScoreAIC(const PenObj& pen, bool verbose = false) const;

    /*! \brief get the PEN value.

    Get the value of the PEN given penalty function object pen.

    \param pen is the penalty function object
    \param deltaLeaf the number of additional leaves (can be negative) to
    calulate the PEN with:  allows changes in histogram shape to be evaluated.
    */
    real getPENValue(const PenObj& pen, int deltaLeaf = 0) const;

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

	//src_trunk_0701
    /*! \brief Get whether this has a subpaving to manage.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return true if this has a subpaving to manage.
	false otherwise.*/
    bool hasSubPaving() const;
    
    //src_trunk_0701
    	/*! \brief Get the box of the subpaving managed by this.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return copy of the box of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
	cxsc::ivector getRootBox() const;

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

    /*! \brief Insert a single data point into AdaptiveHistogram object.

    A method which  attempts to insert a datapoint into this AdaptiveHistogram
    object's dataCollection and associate the datapoint with one of the
    leaves of the tree pointed to by this's rootPaving.
    If the datapoint falls outside the boundaries of the rootBox it will
    not be inserted and a message will be printed to standard output.

    \param newdata the datapoint to be inserted, an rvector.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).

    \pre this AdaptiveHistogram object must have a rootPaving pointing to
    an initialised SPSnode; this SPSnode may already have data associated
    with it.
    \post If the datapoint is within the boundaries of the rootBox it
    will have been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection.
    */
    void insertOne(rvector newdata,  const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging=NOLOG);


    
    /** @name Insert 1-d data from txt file into AdaptiveHistogram object.

    A group of overloaded versions of a functions which read
    in lines of data representing doubles or integers from a txt file.

    Expects one integer or double per line.  No input validation or checking.
    Rejects lines with input that cannot be converted to integer or double.
    Ignore any extra 'words' on the line (anything after first white space)

    For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    12.04.  A string "30 abc 12.04.0006" will be read as  30, and a string
    "abc 30 12.04.0006" will be rejected.

    \param s the name of the txt file to read data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \pre A file with filename s in the same directory as main() or with the
    filename incorporating directory location.
	\param headerlines is number of headerlines to skip before reading 
	data.  Defaults to 0.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \return true if at least some data inserted, false if txt file cannot
    be read, is empty, or if something failed in after data is read and
    when it is being inserted.
    */
    //@{

    /** All data is associated with the root paving, no splitting. */
    bool insertOneDimDataFromTxt(const std::string& s,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging = NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertOneDimDataFromTxt(s, sn, headerlines, logging);
    }

    /** Adaptive splitting with each data point inserted. */
    bool insertOneDimDataFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging = NOLOG);
    //@}



    /** @name Insert doubles in a txt file into AdaptiveHistogram object.

    A group of overloaded versions of a functions which read
    in lines of data representing doubles from a txt file.

    Expects one integer or double per line.  No input validation or checking.
    Can deal with doubles with no decimal point.
    Rejects lines with input that cannot be converted to double.
    Ignore any extra 'words' on the line (anything after first white space)

    For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    12.04.  A string "30 abc 12.04.0006" will be read as  30, and a string
    "abc 30 12.04.0006" will be rejected.

    \deprecated Function kept for backwards compatibility, effectively replaced
    by insertOneDimDataFromTxt.

    \param s the name of the txt file to read data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \pre A file with filename s in the same directory as main() or with the
    filename incorporating directory location.
	\param headerlines is number of headerlines to skip before reading 
	data.  Defaults to 0.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \return true if at least some data inserted, false if txt file cannot
    be read, is empty, or if something failed in after data is read and
    when it is being inserted.
    */
    //@{

    /** All rvectors are associated with the root paving, no splitting. */
    bool insertDoublesFromTxt(const std::string& s,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertOneDimDataFromTxt(s, sn, headerlines, logging);
    }

    /** Adaptive splitting with each data point inserted. */
    bool insertDoublesFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG)
    {
        return insertOneDimDataFromTxt(s, boolTest, logging);
    }
    //@}


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
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertRvectorsFromTxt(s, sn, headerlines, logging);
    }


    /** Adaptive splitting with each data point inserted. */
    bool insertRvectorsFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);
    //@}


    /** @name Insert all rvectors from a container of rvectors.

    \param rvec is the container of rvectors to get data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \return true if at least some data inserted, false if there are no
    labeled points with the specified label in the RSSample object.
    */
    //@{

    /*! All rvectors are associated with the root paving, no splitting. */
    bool insertFromRVec(const RVecData& rvec, LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertFromRVec(rvec, sn, logging);
    }

    /*! Adaptive splitting with each data point inserted. */
    bool insertFromRVec(const RVecData& rvec, const SplitDecisionObj& boolTest,
                             LOGGING_LEVEL logging=NOLOG);

    //@}

    /** @name Insert a set number of rvectors from a container.

    The overloaded versions which take a random number generator as a parameter
    can be used to take successive samples from the same container.
    Otherwise the random number generator is created and destroyed during
    the scope of the function and repeating the identical function call will
    produce an identical sample from the container.

    \param samplesize the size of the sample to draw.
    \param gsl_rng * rgsl is a random number generator.
    \param seed is a seed for a random number generator.
    \param rvec the container to get data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \return true if at least some data inserted, false otherwise.
    */
    //@{

    /** All rvectors are associated with the root paving, no splitting,
    random number generator supplied. */
    bool insertSampleFromRVec(size_t samplesize, gsl_rng * rgsl,
            const RVecData& rvec, LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRVec(samplesize, rgsl, rvec,
                                    sn, logging);
    }

    /** All rvectors are associated with the root paving, no splitting,
    seed for creating a random number generator supplied. */
    bool insertSampleFromRVec(size_t samplesize, int seed,
            const RVecData& rvec, LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRVec(samplesize, seed, rvec,
                                    sn, logging);
    }
    /** All rvectors are associated with the root paving, no splitting,
    no random number generator supplied, default will be created. */
    bool insertSampleFromRVec(size_t samplesize,
            const RVecData& rvec, LOGGING_LEVEL logging=NOLOG)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRVec(samplesize, rvec, sn, logging);
    }

    /** Adaptive splitting with each data point inserted,
    random number generator supplied.*/
    bool insertSampleFromRVec(size_t samplesize, gsl_rng * rgsl,
            const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging=NOLOG);
    /** Adaptive splitting with each data point inserted,
    seed for creating random number generator supplied.*/
    bool insertSampleFromRVec(size_t samplesize, int seed,
            const RVecData& rvec, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging=NOLOG);
    /** Adaptive splitting with each data point inserted,
    no random number generator supplied, default will be created.*/
    bool insertSampleFromRVec(size_t samplesize, const RVecData& rvec,
            const SplitDecisionObj& boolTest, LOGGING_LEVEL logging=NOLOG);

    //@}



    /** @name Insert all rvectors from an RSSample object.

    Insert rvectors from the labeled point Samples of an RSSample
    object where the label for the point matches the label specified
    as a argument to the function.

    \param rss the RSSample object to get data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file.
    \param label the label for the labeled points in the rss.Samples
    that we want to use.
    \return true if at least some data inserted, false if there are no
    labeled points with the specified label in the RSSample object.
    */
    //@{

    /** All rvectors are associated with the root paving, no splitting. */
    bool insertFromRSSample(const RSSample& rss,
                LOGGING_LEVEL logging, int label = 0)
    {
        SplitNever sn; // a dummy split decision object
        return insertFromRSSample(rss, sn, logging, label);
    }

    /** Adaptive splitting with each data point inserted. */
    bool insertFromRSSample(const RSSample& rss,
                            const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging, int label = 0);
    //@}


    /** @name Insert a set number of rvectors from an RSSample object.

    Insert a set number of rvectors from the labeled point Samples of an
    RSSample object where the label for the point matches the label specified
    as a argument to the function.  Sampling is random with replacement.

    The overloaded versions which take a random number generator as a parameter
    can be used to take successive samples from the same RSSample object.
    Otherwise the random number generator is created and destroyed during
    the scope of the function and repeating the identical function call will
    produce an identical sample from the RSSample object.

    \param samplesize the size of the sample to draw.
    \param gsl_rng * rgsl is a random number generator.
    \param seed is a seed for a random number generator.
    \param rss the RSSample object to get data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param label is the label for the labeled points in the rss.Samples
    which we want to sample from.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \return true if at least some data inserted, false if there are no
    labeled points with the specified label in the RSSample object.
    */
    //@{

    /** All rvectors are associated with the root paving, no spliting,
    random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize, gsl_rng * rgsl,
            const RSSample& rss, LOGGING_LEVEL logging, int label = 0)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rgsl, rss, sn,
                                        logging, label);
    }
    /** All rvectors are associated with the root paving, no spliting,
    seed for creating a random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, LOGGING_LEVEL logging, int label = 0)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, seed, rss, sn,
                                        logging, label);
    }
    /** All rvectors are associated with the root paving, no spliting,
    no random number generator supplied, default will be created. */
    bool insertSampleFromRSSample(size_t samplesize,
            const RSSample& rss, LOGGING_LEVEL logging, int label = 0)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rss, sn,
                                        logging, label);
    }

    /** Adaptive splitting with each data point inserted,
    random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize, gsl_rng * rgsl,
            const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging, int label = 0);
    /** Adaptive splitting with each data point inserted,
    seed for creating random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging, int label = 0);
    /** Adaptive splitting with each data point inserted,
    no random number generator supplied, default will be created.*/
    bool insertSampleFromRSSample(size_t samplesize, const RSSample& rss,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging, int label = 0);

    //@}


    /** @name prioritySplit methods.

    These methods takes a histogram and progressively split using a priority
    queue to determine which node to split first. Splitting continues until
    some criteria applying either to individual nodes or to the histogram
    as a whole is satisfied, or there are no more splittable nodes.

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
    \return true if the priority split was successful, false otherwise.
    */
    //@{
    /** Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, minChildPoints, 0.0, 
										maxLeafNodes); }

    /** Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, 0, minVolB, maxLeafNodes); }

    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, 0, 0.0, maxLeafNodes); }

    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
							 double minVolB, size_t maxLeafNodes);

    /** With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl, 
							 size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, minChildPoints, 0.0, rgsl,
									  maxLeafNodes); }

    /** With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl, size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, 0, minVolB, rgsl, maxLeafNodes); }

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes)
    { return prioritySplit(compTest, he, logging, 0, 0.0, rgsl, maxLeafNodes); }

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes);

    //@}

    /*! \brief Priority merge to reduce number of leaves in histogram.

    This method takes a histogram where where all the data is associated with
    multiple nodes and progressively merges the children of sub-terminal leaves
    using a priority queue to determine which node to merge first.
    Merging continues until some criteria applying either to individual
    nodes or to the histogram as a whole is satisfied, or the histogram only
    has one bin.

    If more than one node is equally 'small', on the basis of the node
    comparison compTest used, then a random choice is made between all equally
    small nodes to find the node which will be merged.

    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to determine
    when to stop merging.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)
    \return true if the priority merge was successful, false otherwise.
    */
    bool priorityMerge(const NodeCompObj& compTest, const HistEvalObj& he,
                        LOGGING_LEVEL logging=NOLOG);

    /*! \brief Merge a multileaf histogram up to just root box.

    No prioritistion, just brute force
    */
    bool mergeUp();

    /*! \brief Split a histogram  to a specified shape.

    \param instruction specifies the required shape, eg "3, 3. 2, 1"
    */
    bool splitToShape(std::string instruction);


    /*! \brief Outputting MCMC samples from histogram state space.

    The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  In this implementation the
    the Metropolis-Hastings algorithm is used on the SPSnode tree managed by
    this Adaptive Histogram to generate samples from the histogram state
    probability density.

    MCMC requires a prior distribution over the state space and a likelihood of
    being in a particular state given the data, from which can be found (up to
    proportionality) a posterior distribution proportional to the
    likelihood x prior.

    A change in state can take place by either splitting a leaf node (bisecting
    the box represesented by that leaf and sending the data associated with the
    box down to the two new children) or absorbing two sibling leaf nodes back
    into their parent so that that parent becomes a leaf of the tree.  The
    parent node of two sibling leaf nodes is referred to as a 'cherry' and the
    reabsorbtion of the sibling leaf children of a cherry is referred to as
    'merging' a cherry.

    minPoints restricts the possible states by not allowing the chain to
    include a state where a leaf node has less than minPoints data points
    associated with it, unless that child leaf node has 0 points and its sibling
    has all the parent's points and number of parent's points >= minPoints.
    The final condition allows a node to be split when all the data goes to
    just one child provided that the number of points in the node >= minPoints
    so that the process can 'home in' on small peaks of data.

    When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can collate thinned out samples from the chain of states
    generated from the Markov Chain process and average them, outputting
    txt file represenations of the collation and the average.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.

    \param loops how many states in the chain to generate altogether,
    (ie including burnin).
    \param burnin how many states to consider as burnin, ie the number of
    changes in state required for the process to forget its starting state.
    \param thinout what step size to use in sampling for averaging.
    If thinout > 0, an average histogram is created by collating
    the first state in the chain after burnin states together with samples
    taken from the chain with thinout stepsize.  Txt and .dot file representions
    of the sampling process and the collation and the average are output.
    \param proposal is a reference to the proposal distribution object.
    \param logPrior is a reference to the prior distribution object.
    \param minPoints is the minimum number of points to allow in a box
    represented by a leaf in the SPnode tree (defaults to 0).
    \param logging an enum controlling whether histogram creation output is
    sent to log files (defaults to no logging).  TXT gives logging to a
    txt file, TXTANDGRAPH gives txt file logging and graphs,
    LOGSAMPLES does a txt log file for samples only, LOGANDGRAPHSAMPLES does
    a txt log file and dot graphs for samples only
    \return true if MCMC successfully carried out for required number of
    states in the chain (loops through the MCMC routine), false otherwise.
    */
    bool MCMC(unsigned int loops, unsigned int burnin, unsigned int thinout,
                MCMCProposal& proposal, LogMCMCPrior& logPrior,
                size_t minPoints = 0, LOGGING_LEVEL logging=NOLOG);

		/*! \brief Generating MCMC samples from histogram state space.

    The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  In this implementation the
    the Metropolis-Hastings algorithm is used on the SPSnode tree managed by
    this Adaptive Histogram to generate samples from the histogram state
    probability density.

    MCMC requires a prior distribution over the state space and a likelihood of
    being in a particular state given the data, from which can be found (up to
    proportionality) a posterior distribution proportional to the
    likelihood x prior.

    A change in state can take place by either splitting a leaf node (bisecting
    the box represesented by that leaf and sending the data associated with the
    box down to the two new children) or absorbing two sibling leaf nodes back
    into their parent so that that parent becomes a leaf of the tree.  The
    parent node of two sibling leaf nodes is referred to as a 'cherry' and the
    reabsorbtion of the sibling leaf children of a cherry is referred to as
    'merging' a cherry.

    minPoints restricts the possible states by not allowing the chain to
    include a state where a leaf node has less than minPoints data points
    associated with it, unless that child leaf node has 0 points and its sibling
    has all the parent's points and number of parent's points >= minPoints.
    The final condition allows a node to be split when all the data goes to
    just one child provided that the number of points in the node >= minPoints
    so that the process can 'home in' on small peaks of data.

    When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can collate thinned out samples from the chain of states
    generated from the Markov Chain process and average them, outputting
    txt file represenations of the collation and the average.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.

    \param samples is a reference to a container to add 
    AdaptiveHistogram samples to.
    \param loops how many states in the chain to generate altogether,
    (ie including burnin).
    \param burnin how many states to consider as burnin, ie the number of
    changes in state required for the process to forget its starting state.
    \param thinout what step size to use in sampling for averaging.
    If thinout > 0, samples are the first state in the chain after 
    burnin states together with samples taken from the chain with 
    thinout stepsize.
    \param proposal is a reference to the proposal distribution object.
    \param logPrior is a reference to the prior distribution object.
    \param minPoints is the minimum number of points to allow in a box
    represented by a leaf in the SPnode tree (defaults to 0).
    \param logging an enum controlling whether histogram creation output is
    sent to log files (defaults to no logging).  TXT gives logging to a
    txt file, TXTANDGRAPH gives txt file logging and graphs,
    LOGSAMPLES does a txt log file for samples only, LOGANDGRAPHSAMPLES does
    a txt log file and dot graphs for samples only
    \return the reference to the container of samples if MCMC 
    successfully carried out for required number of states in the chain
    empty container otherwise.
    */
	std::vector < AdaptiveHistogram >& MCMCsamples(
						std::vector < AdaptiveHistogram >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging);


    //docs and make private?
    void publicOutputMCMCStateSample(int ci, int i, bool confirm = false);

    //docs
    void publicLogMCMCSample(std::string s, int i);

    /*! \brief Changes the state of this Adaptive Histogram using MCMC process.

    This method proposes and probabilistically accepts a single-step change
    in the histogram state represented by this Adaptive Histogram, i.e a change
    in the SPSnode tree managed by the histogram and representing the partion of
    the data space into histogram 'bins'.  A single-step change in
    state is the change in state resulting from a single split of a leaf node
    or merge of a cherry node.

    Proposals are made by selecting a leaf or cherry at random.  If a cherry is
    chosen the proposed change in state is to merge the two leaf children of
    that cherry; if a leaf is chosen the proposal is to split (bisect) the leaf.

    minPoints > 0 restricts the leaves which can be selected for change: if
    splitting a leaf would result in the children having less than minPoints
    data points associated with them, the leaf cannot be selected
    for the proposed change unless that child leaf node has 0 points and its
    sibling has all the parent's points (>= minPoints). The final condition
    allows a node to be split when all the data goes to just one child provided
    that the number of points in the node >= minPoints so that the process can
    'home in' on small peaks of data..  Note however that minPoints has no
    effect on the cherries which can be proposed for merging.  A leaf whose
    split would result in both children having at least minPoints points in it,
    or where the split would give all the data to one child, is referred
    to as a 'splittable leaf'.

    Thus the implementation of the MCMC algorithm for generating a new state
    in the chain needs to maintain separately the overall state of the tree
    and the set of <b>splittable leaf nodes</b>.

    The log-likelihood of a state given the data is given by
    sum over leaves of (counts in leaf x -ln(count in leaf / (n x vol of leaf)))
    where n is the total number of data points in the histogram.

    The posterior distribution is proportional to the prior x likelihood.

    The Metropolis-Hastings algorithm also requires a proposal density which
    depends on the current state m to generate a proposed state m'.

    Q(m' | m) is the transition probability from state m to state m'.

    This proposal is accepted if u drawn from Uniform(0,1) is such that
    u < (posterior probability of m' x Q(m | m'))/(posterior probability
    of m x Q(m' | m).  In this implementation, natural logs are used to
    simplify calculation, ie a proposal is accepted if
    log(u) < log [(posterior probability of m' x Q(m | m'))/(posterior
    probability of m x Q(m' | m)].

    If the proposal is accepted the state is changed to m', otherwise it
    stays at m.

    \pre an SPSnode tree representing some histogram state, a container of
    pointers to the splittable leaf and cherry nodes of the tree in its
    current state, ordered with the leaves first followed by the cherries,
    the number of splittable leaf nodes and the number of cherry nodes.
    \post the SPSnode tree in some state which is either in the same state as
    the pre-state or in a new state reachable by a single-step change
    (i.e. through the split of one splittable leaf or the merge of one cherry
    of the pre-state), the container of splittable leaf and cherry nodes
    updated for any change in state, and the number of splittable leaves
    and cherries updated similarly.

    \note that the container of splittable leaf nodes and cherries is
    maintained separately from the overall state of the tree to save having to
    repeatedly assess whether a leaf can be split (given minPoints).  The
    number of splittable leaves and cherries is maintained for convenience to
    avoid repeatedly counting nodes of different types in the container.

    \param nodes is a reference to a container of pointers to the leaf and
    cherry nodes of the SPSnode tree managed by this AdaptiveHistogram, which
    will be updated if the method results in change in tree state.
    \param numLeaves is a reference to a variable storing the number of leaves
    in the SPSnode tree, which will be updated if the method results in change
    in tree state.
    \param numCherries is a reference to a variable storing the number of
    cherries in the SPSnode tree, which will be updated if the method results
    in a change in tree state.
    \param proposal is a reference to a proposal distribution object.
    \param logPrior is a reference to a log prior object.
    \param minPoints is the minimum number of points allowed in a box.
    \param rgsl is a pointer to a uniform random number generator.
    \param logging an enum controlling whether histogram creation output is
    sent to log files (defaults to no logging).  TXT gives logging to a
    txt file, TXTANDGRAPH gives txt file logging and graphs.
    \param s is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain.
    \return true if there has been a successful proposal and successful
    probabilistic acceptance/rejection of the change in state, false otherwise.
    */
    bool changeMCMCState (SPSnodeList& nodes, size_t& numLeaves,
                            size_t& numCherries,
                            MCMCProposal& proposal, LogMCMCPrior& logPrior,
                            size_t minPoints, gsl_rng* rgsl,
                            LOGGING_LEVEL logging, std::string s, int i);

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

    /*! \brief Output the subpaving managed by this to a txt file.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then node contribution to EMP
    under COPERR, the change that would result in the EMP under COPERR if the
    node were split, the node contribution to EMP under AIC, the change that
    would result in the EMP under AIC if the node were split, and then the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param s the name of the txt file to send output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputToTxtTabsWithEMPs(const std::string& s,
                                bool confirm = false) const;


    /*! \brief Output details of full sample (from root) to txt tile.

    Format is a mixture of alpha and  numeric data.

    \param s the name of the txt file to send output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputRootToTxt(const std::string& s, bool confirm = false) const;
	 
	 //--src_trunk_0701
	 	/*! \brief Change this so that the subpaving it manages is
	the union of this's subpaving and the subpaving of that of a 
	PiecewiseConstantFunction.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer or if the subpaving managed by \a other
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving of this has no box
	or if the subpaving of \a other has no box.
	
	Throws an IncompatibleDimensions_Error if the
	subpaving boxes of this and \a other are not identical.

	There will be no change in this if the subpaving of \other is
	everywhere less split than the subpaving of this.
	
	\param other is the %PiecewiseConstantFunction to make the union against.
	\pre Both this and \a other have subpavings with boxes to manage.
	\pre The boxes of the subpavings of this and \a other are the same. 
	\post the subpaving managed by this has the shape that is the
	union of its shape before the operation and the shape of 
	the subpaving managed by \a other.  \a other is unchanged.      */
 	void reshapeToUnion(const PiecewiseConstantFunction& other);

	/*! \brief Change this so that the subpaving it manages is
	as close as possible to the union of this's subpaving 
	and the subpaving of that of a 
	PiecewiseConstantFunction.
	
	If \a other has a subpaving that is more split than
	the subpaving managed by this at any node, this will
	 not exactly follow the shape of \a other if the resulting
	 nodes would not splittable according to 
	 SPSnode::isSplittableNode(size_t minChildPoints).  If any node 
	 cannot be split to follow the shape of \a other due to 
	 \a minChildPoints, a message will be printed to std::cerr.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer or if the subpaving managed by \a other
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving of this has no box
	or if the subpaving of \a other has no box.
	
	Throws an IncompatibleDimensions_Error if the
	subpaving boxes of this and \a other are not identical.

	There will be no change in this if the subpaving of \other is
	everywhere less split than the subpaving of this.
	
	\param other is the %PiecewiseConstantFunction to make the union against.
	\param minChildPoints is the minumum child points to use
	to check if this can be split in order to follow \a other.
	\pre Both this and \a other have subpavings with boxes to manage.
	\pre The boxes of the subpavings of this and \a other are the same. 
	\post the subpaving managed by this has the shape that is as
	close as possible to the
	union of its shape before the operation and the shape of 
	the subpaving managed by \a other, given \a minChildPoints.
	\a other is unchanged.      */
 	void reshapeToUnion(const PiecewiseConstantFunction& other,
						size_t minChildPoints);
	//--src_trunk_0701
	 
	 
	 /*! \brief Clear the histogram's data and counters
	 */
	 void makeEmpty();
	 
	 /** @name Distribution-free Likelihood Estimation
	 \param RVecData the observed data
	 or
	 \param RSSample the observed data
	 \param dx is a user-defined dx - see AHABC meeting notes
	 \param WeightHist weight for the histogram in model 1
	 \param WeightsPM weights for model 0
 	 \param wt is the mass added to the histogram to ensure positive
	  probability everywhere

	 \return (log)-likelihood estimate
	 */
	//@{
	/** Get the estimated log likelihood from RSSample. */
   real getEstLogLikelihoodFromRSSample(RSSample& labSampledData,
										double dx, double wt)
	{	double WeightHist = 1.0;
		std::map<rvector, double, std::less<rvector> > WeightsPM;
		return getEstLogLikelihoodFromRSSample(labSampledData,
										dx, wt, WeightHist, WeightsPM); }	
	real getEstLogLikelihoodFromRSSample(RSSample& labSampledData,
										double dx, double wt,  double WeightHist,
						std::map<rvector, double, std::less<rvector> >& WeightsPM);
		
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
	
	/** Get the IAE for a finite gaussian mixture distribution using interval 
	 * techniques. */
	 cxsc::interval getFinMixIntervalIAE(FinMix& mixt, double tol, int deg);
	
	/** Get the IAE of a bivariate gaussian/Levy 2D/Rosen 2D distribution. */
	cxsc::real get2DIAE(taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval));
	
	/** Get the IAE for a uniform mixture distribution. */ 	
	cxsc::real getUnifIAE(AdaptiveHistogram& myPart, double weight, std::vector<int> holesLoc);
	
	/** Get the IAE for a uniform distribution. */ 	
	cxsc::real getUnifIAE();
	
	/** Get the IAE for mapped function. */
	//cxsc::real getMappedFunctionIAE(RealMappedSPnode nodeEst);
	
	/** Get the IAE for a laplace distribution with mu=0 and b=1 using interval 
	 * techniques. */
	 cxsc::interval getLaplaceIntervalIAE(double tol, int deg);
	 
	 	/** Get the IAE for a laplace distribution with mu=0 and b=1 using interval 
	 * techniques. */
	 cxsc::interval getLognormalIntervalIAE(double tol, int deg);
	//@}

	/*! \brief Find the coverage of the boxes in an AdaptiveHistogram.
	*/
	void findDensityRegion(double cov, double weightPM, vector<SPSnode*> & covNodes,
									string covFileName);
									
	
	// Jenny addition for Gloria's convergence work
	/*! \brief Append current state of histogram to a txt log file.

	Format is a tab-delimited file of numeric data.
	Output is plain: just vols, counters, and boxes.

	\param s the name of the txt file to send output to.
	\param i the number of pass (ie, 0, 1, 2, 3 etc) in process
	*/
	void outputLogPlain(const std::string& s, const int i) const;
	
	// gloria addition for makingADHtoMapped
	/*! \brief Returns a RealMappedSPnode from the AdaptiveHistogram object.
	*/
	//cxsc::real getMappedIAE(RealMappedSPnode& nodeEst, ivector pavingBox) const;
	
	
		/** @name prioritySplit and getting the states using index from a distr
    */
    //@{
    /** Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, size_t maxLeafNodes,
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, minChildPoints, 0.0, 
										maxLeafNodes, States, Sampled); }

    /** Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, 0, minVolB, maxLeafNodes, 
										 States, Sampled); }

    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, 0, 0.0, maxLeafNodes,
										 States, Sampled); }

    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
							 double minVolB, size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled);

    /** With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl, 
							 size_t maxLeafNodes,
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, minChildPoints, 0.0, rgsl,
									  maxLeafNodes, States, Sampled); }

    /** With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl, size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, 0, minVolB, rgsl, maxLeafNodes, 
										 States, Sampled); }

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, std::vector<size_t>& Sampled)
    { return prioritySplitGet(compTest, he, logging, 0, 0.0, rgsl, maxLeafNodes, 
									    States, Sampled); }

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitGet(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes, 
                      std::vector<AdaptiveHistogram>& States, 
                      std::vector<size_t>& Sampled);
    //@}
	

	//gat41
	/** @name prioritySplit and getting likelihood*prior (which is proportional 
	    to posterior)
    */
    //@{
    /** Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, size_t maxLeafNodes,
                      std::vector<real>& Posterior, LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, minChildPoints, 0.0, 
										maxLeafNodes, Posterior, logPrior); }

    /** Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, size_t maxLeafNodes, std::vector<real>& Posterior,
                      LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, 0, minVolB, maxLeafNodes, 
    Posterior, logPrior); }

    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes, 
                      std::vector<real>& Posterior, LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, 0, 0.0, maxLeafNodes, Posterior,
    logPrior); }

    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
							 double minVolB, size_t maxLeafNodes, 
							 std::vector<real>& Posterior, LogMCMCPrior& logPrior);

    /** With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl, 
							 size_t maxLeafNodes, std::vector<real>& Posterior,
							 LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, minChildPoints, 0.0, rgsl,
									  maxLeafNodes, Posterior, logPrior); }

    /** With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl, size_t maxLeafNodes, 
                      std::vector<real>& Posterior, LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, 0, minVolB, rgsl, maxLeafNodes, 
    Posterior, logPrior); }

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes, 
                      std::vector<real>& Posterior, LogMCMCPrior& logPrior)
    { return prioritySplitMCMC(compTest, he, logging, 0, 0.0, rgsl, maxLeafNodes, 
    Posterior, logPrior); }

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitMCMC(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes, std::vector<real>& Posterior,
							 LogMCMCPrior& logPrior);
    //@}

	// gloria addition for total variation as a stopping criteria in prioritySplit
    /** @name prioritySplitWithTotalVar methods.
    */
    //@{
    /** Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, size_t maxLeafNodes, int StopVal,
                      std::vector<AdaptiveHistogram>& HistAtValley, int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, minChildPoints, 0.0, 
										maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, size_t maxLeafNodes,int StopVal,
                      std::vector<AdaptiveHistogram>& HistAtValley, int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, 0, minVolB, 
							maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes, int StopVal,
                      std::vector<AdaptiveHistogram>& HistAtValley, int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, 0, 0.0, 
								maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
							 double minVolB, size_t maxLeafNodes, int StopVal,
							 std::vector<AdaptiveHistogram>& HistAtValley, int simNum);

    /** With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl, 
							 size_t maxLeafNodes, int StopVal, 
							 std::vector<AdaptiveHistogram>& HistAtValley, int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, minChildPoints, 0.0, rgsl,
									  maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl, size_t maxLeafNodes, 
                      int StopVal, std::vector<AdaptiveHistogram>& HistAtValley,
                      int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, 0, minVolB, rgsl, 
								maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes, 
                      int StopVal, std::vector<AdaptiveHistogram>& HistAtValley,
                      int simNum)
    { return prioritySplitWithTotalVar(compTest, he, logging, 0, 0.0, rgsl, 
						maxLeafNodes, StopVal, HistAtValley, simNum); }

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitWithTotalVar(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes, int StopVal, 
							 std::vector<AdaptiveHistogram>& HistAtValley, int simNum);
        //@}
   
   // gloria addition for prioritySplit with switches
    /** @name prioritySplitWithSwitches methods.
    */
    //@{
    /** Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, size_t maxLeafNodes, int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, minChildPoints, 0.0, 
										maxLeafNodes, removeBox); }

    /** Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, size_t maxLeafNodes, int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, 0, minVolB, 
							maxLeafNodes, removeBox); }

    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes, int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, 0, 0.0, 
								maxLeafNodes, removeBox); }

    /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
							 double minVolB, size_t maxLeafNodes, int removeBox);

    /** With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl, 
							 size_t maxLeafNodes, int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, minChildPoints, 0.0, rgsl,
									  maxLeafNodes, removeBox); }

    /** With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl, size_t maxLeafNodes, 
                      int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, 0, minVolB, rgsl, 
								maxLeafNodes, removeBox); }

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes, 
                      int removeBox)
    { return prioritySplitWithSwitches(compTest, he, logging, 0, 0.0, rgsl, 
						maxLeafNodes, removeBox); }

    /** With random number generator. All other parameters supplied.*/
    bool prioritySplitWithSwitches(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes, int removeBox);
        //@}
   
   /*! Check if we should stop splitting using total variation distance
    */
    bool checkStopCrit(double StopCritCurrent, double StopCritPrevious, int& Prev); 
    
    //src_truk_0701
    	/*! \brief Swap the contents of this and another histogram.
	 
	\post After the swap \a adh will manage the subpaving that this used
	to manage, and this will manage the subpaving that \a adh used to 
	managed, and the values of the other data members of this
	and \a adh will also be swapped.  */
	void swap(AdaptiveHistogram& adh); // throw()

/** @name prioritySplit methods which will ouput the IAE for mapped functions at each split.
*/
  //@{
   /** minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplitMappedIAE(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t minChildPoints, 
											double minVolB, size_t maxLeafNodes,
  										PiecewiseConstantFunction& nodeEst,
  										std::vector<real>& vecIAE, std::vector<int> sequence);
    
    /** Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplitMappedIAE(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, size_t maxLeafNodes, 
                      PiecewiseConstantFunction& nodeEst,
                      std::vector<real>& vecIAE, std::vector<int> sequence)
    { return prioritySplitMappedIAE(compTest, he, logging, 0, 0.0, maxLeafNodes, nodeEst, vecIAE, sequence); }

   /** With random number generator. All other parameters supplied.*/
    bool prioritySplitMappedIAE(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl,
							 size_t maxLeafNodes, PiecewiseConstantFunction& nodeEst,
							 std::vector<real>& vecIAE, std::vector<int> sequence);

    /** With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplitIAEMappedIAE(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl, size_t maxLeafNodes, 
                      PiecewiseConstantFunction& nodeEst,
                      std::vector<real>& vecIAE, std::vector<int> sequence)
    { return prioritySplitMappedIAE(compTest, he, logging, 0, 0.0, rgsl, maxLeafNodes, nodeEst, vecIAE, sequence); }
     
  //@}
}; // end of AdaptiveHistogram class declarations


/*! a class for histogram-related exceptions.
*/
class HistException : public std::exception
{
   std::string s;
   public:
   HistException(std::string ss);
   ~HistException () throw ();
   virtual const char* what() const throw();
};



/*! \brief A class for the histogram description.

Describes a histogram as a ordered collection of leaf node depths.
Root has depth 0.

\todo complete this class and use instead of strings in describing
histograms.

*/

class HistDescription {
private:
    /*! The depth string, a string of comma separated levels.
    */
    std::string depthString;

    /*! A flag for whether the instruction string is good.
    */
    mutable bool goodString;


public:
    /*! \brief  Default constructor.
    */
    explicit HistDescription() : depthString(""), goodString(true) {}

    /*! \brief  Constructor from string.
    */
    explicit HistDescription(std::string str) : depthString(str),
                                                goodString(true) {}

    /*! \brief  Copy constructor.
    */
    explicit HistDescription(const HistDescription& other)
                : depthString(other.depthString),
                  goodString(other.goodString) {}

    /*! \brief Copy assignment operator.
    */
    HistDescription& operator=(const HistDescription& rhs);

    /*! \brief  Set the good string flag.
    */
    void setIsGood(bool b)
    { goodString = b; }

    /*! \brief  Get the first level in the description.
    */
    int peekFirst() const;

    /*! \brief  Strip off the first level in the description.
    */
    bool popFirst();

    /*! \brief Get the goodstring flag
    */
    bool getIsGood() const
    { return goodString; }

    /*! \brief Get the depthString
    */
    std::string getDepthString() const
    { return depthString; }

};



// ----------  declarations of non-member functions ----------------------

/*! \brief Output the HistDescription object.
*/
std::ostream & operator<<(std::ostream &os,
                                    const HistDescription& hd);

/*! \brief Comparison operator for the histogram description.
*/
bool operator<(const HistDescription& lhs,
                            const HistDescription& rhs);


/*! \brief Output the contents of an AdaptiveHistogram object.

Verbose output for an AdaptiveHistogram object, including all boxes
(not just leaves), data, and summary statistics.
*/
std::ostream & operator<<(std::ostream &os, const AdaptiveHistogram& adh);

//=====================end of non-member functions declarations=================

} // end of namespace subpavings

/*! A specialisation of std::swap for AdaptiveHistogram types.*/
namespace std
{
	template <>
	void swap (subpavings::AdaptiveHistogram & a1, 
			subpavings::AdaptiveHistogram & a2);
}

#endif

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

/*! \file      adaptivehistogramcollator.hpp
\brief AdaptiveHistogramCollator declarations.
*/

#ifndef ___ADAPTIVEHISTCOLL_HPP__
#define ___ADAPTIVEHISTCOLL_HPP__

#include "adaptivehistogram.hpp"
#include "sptypes.hpp"
#include "histevalobj.hpp"

namespace subpavings {


/*! \brief A wrapper or manager for a CollatorSPSnode.

AdaptiveHistogramCollator class objects manage \link subpavings::CollatorSPnode
CollatorSPnode \endlink objects for the purpose of collating information
from a number of AdaptiveHistogram objects.

The AdaptiveHistogramCollator's CollatorSPnode tree represents the subpaving
that is the union of all the subpavings associated with each AdaptiveHistogram
in the collation.  Each node in the tree has a data member which is a container
structure holding one value for each collated histogram.  For a collation of
histograms the container holds the <b>normalised height</b>, for each collated histogram,
of the histogram bin represented by the box of that node.

(The normalised height associated with a bin which is represented by the box
of a leaf node of a tree managed by an AdaptiveHistogram object is the number
of data points associated with that bin divided by (the total number of
data points in the histogram x the volume of the bin).  Thus the areas
(heights x volumes) of the bins sum to 1.

Since the tree represents the union of the subpavings associated with each
AdaptiveHistogram in the collation, the collation tree will have at least as
many and usually more bins than any of the collated histograms.

Each collated histogram will be represented in the summaries in the order in
which it was added to the collation.  Eg, the heights of the bins of the first
histogram to be collated will be first (index [1]) in the summary container.

If the collated AdaptiveHistograms have been properly formed and added to the
collation, the sum, over all the leaf nodes of the collation, of the volume of
the box of the leaf node multiplied by the values of the leaf node summary
corresponding a particular collated histogram will be 1 for each collated
histogram.

*/

class AdaptiveHistogramCollator {
private:

    /*! \brief Private initialised constructor.

    Initialised  with pointer to subpaving.
    */
    explicit AdaptiveHistogramCollator(CollatorSPnode * spn);

    /*! \brief Pointer to the root CollatorSPnode.

    A CollatorSPSnode is a binary tree representation of information from
    a number of subpavings.

    The summary information held by the CollatorSPnode for an SPSnode is the
    count/(total count in tree * volume), ie the normalised histogram height
    for a bin corresponding the box of that SPSnode.
    */
    CollatorSPnode* rootCollator;


    /*! \brief Private method for making a collation histogram from RVecData.

    \param samplesize is the number of datapoints in each sample to draw.
    \param numberSamples is the number samples to draw (number of histograms
    to make and include in this collator).
    \param rv is the container of rvectors to draw samples from.
    \param pavingBox is the box to be used as the root box for all histograms
    \param indImmedSplit is an indicator for immediate splitting,
            =1 if we split as data comes (no priority queue)
            =0 otherwise (priority queue splitting).
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to control
    when to stop splitting.
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \return true if the collation process successful, false otherwise.
    */
    bool collateFromRVec(size_t samplesize, size_t numberSamples,
                const RVecData rv, ivector pavingBox,
                int indImmedSplit, const SplitDecisionObj& boolTest,
                const NodeCompObj& compTest, const HistEvalObj& he,
                size_t minChildPoints, double minVolB);


    /*! \brief Private method for making a collation histogram from RSSample.

    \param samplesize is the number of datapoints in each sample to draw.
    \param numberSamples is the number samples to draw (number of histograms
    to make and include in this collator).
    \param rss is the RSSample object to draw samples from.
    \param pavingBox is the box to be used as the root box for all histograms
    \param indImmedSplit is an indicator for immediate splitting,
            =1 if we split as data comes (no priority queue)
            =0 otherwise (priority queue splitting).
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to control
    when to stop splitting.
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \param label is the label for the labeled points in the rss sample which
    we want to draw our samples from.
    \return true if the collation process successful, false otherwise.
    */
    bool collateFromRSSample(size_t samplesize, size_t numberSamples,
                const RSSample rss, ivector pavingBox,
                int indImmedSplit,
                const SplitDecisionObj& boolTest, const NodeCompObj& compTest,
                const HistEvalObj& he, size_t minChildPoints, double minVolB,
                int label);

	/*! \brief Add current state of collation to a log file.

    \param s is the name of the file to log to.
    \param i is a number representing the index of this state in a sequence.
    */
    void outputLog(const std::string& s, const int i) const;

public:
    /*! \brief Default constructor.
    */
    explicit AdaptiveHistogramCollator();

    /*! \brief initialised constructor.
    Initialised with an AdaptiveHistogram object.
    */
    explicit AdaptiveHistogramCollator(const AdaptiveHistogram& adh);

    /*! \brief  Copy constructor
    */
    AdaptiveHistogramCollator
        (const AdaptiveHistogramCollator& other);


    /*! \brief Copy assignment operator.
    */
    AdaptiveHistogramCollator&
        operator=(const AdaptiveHistogramCollator& rhs);

    
    
    /*! \brief Destructor.
    */
    ~AdaptiveHistogramCollator();

	
	/*! Addition operator.

    Addition gives a histogram collator managing a tree which represents a
    subpaving which is the union of the subpavings represented by the operand
    collators.  The summary for each node in the tree contains all the values
    from the summaries of the corresponding nodes in the trees managed by
    the operand AdaptiveHistogramCollators.

    */
    AdaptiveHistogramCollator operator+(const
                AdaptiveHistogramCollator& rhs) const;

    /*! Incremental or inplace addition operator.

    Addition gives this histogram collator an expaned tree which represents a
    subpaving which is the union of the subpavings represented by the this's
    original tree and the tree represented by the rhs collator.  The summary
    for each node in the expandd tree contains all the values
    from the summaries of the original nodes and the summaries of the
    corresponding nodes in the tree managed by the rhs collator.
    */
    AdaptiveHistogramCollator& operator+=(const
                AdaptiveHistogramCollator& rhs);

	/*! Subtraction operator.

    Subtraction gives a histogram collator managing a tree which represents a
    subpaving which is the union of the subpavings represented by the operand
    collators.  The summary for each node in the tree contains all the values
    from the summary of the corresponding node in the trees managed by the left
    hand side operand and the negation of the values from the summary of the
    corresponding node in the trees managed by the right hand side operand.
    */
    AdaptiveHistogramCollator operator-(const
                AdaptiveHistogramCollator& rhs) const;


     /*! \brief Return a pointer to the CollatorPSnode this manages.
    */
    CollatorSPnode* getSubPaving() const;


    
    /*! An AdaptiveHistogramCollator which is the average over this collation.

    Makes and returns an AdaptiveHistogramCollator which is the average over
    the collation of histograms represented by this.  The tree managed
    by the average has structure exactly the same as the tree managed by this
    and one value in the summary for each node where that value is
    the average of the summary of the corresponding node in this.

    This method can only be performed on summaries which do not include negated

    */
    AdaptiveHistogramCollator makeAverage() const;


    /*! \brief Add an AdaptiveHistogram object to the data collation.

    Attempts to add an AdaptiveHistogram object into the collation
    of AdaptiveHistogram information.

    \param adh the AdaptiveHistogram to be included in the collation.
    \post This will include summary data from the adh.  */
    void addToCollation(const AdaptiveHistogram& adh);

    /*! \brief Add the negation of an AdaptiveHistogram object to data collation.

    Attempts to add an AdaptiveHistogram object into the collation
    of AdaptiveHistogram information.

    \param adh the AdaptiveHistogram to be negated and included in the collation.
    \post This will include negated summary data from the adh.  */
    void addNegationToCollation(const AdaptiveHistogram& adh, double c);


    /*! \brief Collate samples from container of rvectors, immediate splitting.

    A method to find the collated histograms from a specified
    number of samples of specifed size from a given container of rvectors.

    Sampling from the container for each histogram is with
    replacement and different samples can are drawn independently: one
    data point in the container can be included in more than
    one histogram and more than once in any one histogram.

    Splitting takes place immediately as the data is inserted.

    \param samplesize is the size of each sample.
    \param numberSamples is the number of samples to draw.
    \param rv is the container of rvectors to draw samples from.
    \param pavingBox is the box to use for all subpaving created.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \return True if data was found and the collation created.
    */
    bool collateFromRVecSplitNow(size_t samplesize, size_t numberSamples,
                            const RVecData rv, ivector pavingBox,
                            const SplitDecisionObj& boolTest);


    /*! \brief Collate samples from container of rvectors, priority queue split.

    A method to find the collated histograms from a specified
    number of samples of specifed size from a given container of rvectors.

    Sampling from the container's samples for each histogram is with
    replacement and different samples can are drawn independently: one
    data point in the container's samples can be included in more than
    one histogram and more than once in any one histogram.

    The histograms are created with priority queue splitting once data is
    initially inserted into the root box of each histogram.

    \param samplesize is the size of each sample.
    \param numberSamples is the number of samples to draw.
    \param rv is the container of rvectors to draw samples from.
    \param pavingBox is the box to use for all subpaving created.
    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to control
    when to stop splitting.
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \return True if data was found and average created.  */
    bool collateFromRVecSplitPQ(size_t samplesize, size_t numberSamples,
                            const RVecData rv, ivector pavingBox,
                            const NodeCompObj& compTest,
                            const HistEvalObj& he,
                            size_t minChildPoints, double minVolB);

    /*! \brief Collate samples from an RSSample object, immediate splitting.

    A method to find the average histogram from a specified
    number of samples of specifed size from a given RSSample object.

    Sampling from the RSSample object's samples for each histogram is with
    replacement and different samples can are drawn independently: one
    data point in the RSSample object's samples can be included in more than
    one histogram and more than once in any one histogram.

    Splitting takes place immediately as the data is inserted.

    \param samplesize is the size of each sample.
    \param numberSamples is the number of samples to draw.
    \param rss is the RSSample object to draw samples from.
    \param pavingBox is the box to use for all subpaving created.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param label is the label of the labeled points in the RSSample object
    which we want to use to take samples from.
    \return True if data was found and collation created.
    */
    bool collateFromRSSampleSplitNow(size_t samplesize, size_t numberSamples,
                            const RSSample rss, ivector pavingBox,
                            const SplitDecisionObj& boolTest, int label = 0);

    /*! \brief Collate samples from RSSample object, priority queue splitting.

    A method to find the average histogram from a specified
    number of samples of specifed size from a given RSSample object.

    Sampling from the RSSample object's samples for each histogram is with
    replacement and different samples can are drawn independently: one
    data point in the RSSample object's samples can be included in more than
    one histogram and more than once in any one histogram.

    The histograms are created with priority queue splitting once data is
    initially inserted into the root box of each histogram.

    \param samplesize is the size of each sample.
    \param numberSamples is the number of samples to draw.
    \param rss is the RSSample object to draw samples from.
    \param pavingBox is the box to use for all subpaving created.
    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to control
    when to stop splitting.
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \param label is the label of the labeled points in the RSSample object
    which we want to use to take samples from.
    \return True if data was found and collation created.
    */
    bool collateFromRSSampleSplitPQ(size_t samplesize, size_t numberSamples,
                            const RSSample rss, ivector pavingBox,
                            const NodeCompObj& compTest,
                            const HistEvalObj& he,
                            size_t minChildPoints, double minVolB,
                            int label = 0);


    /*! \brief Get the number of Adaptive Histogram objects collated.
    */
    size_t getNumberCollated() const;
    
    /*! \brief Get number of leaf nodes.
     */
     size_t getNumLeaves();


    /*! \brief Make a .dot graph file from collated histogram structure.

    Makes a simple .dot graph from the histogram using node names and the
    .png image for this graph.

    \pre a constructed histogram
    \post a .dot file and a .png in the same directory as the program creating
    it was run in.
    */
    bool outputGraphDot() const;

    /*! \brief Output average normalised histogram over collation to a txt file.

    This method does not make the average histogram directly but, for each
    leaf node in the collated tree, calculates and outputs the average of the
    summary associated with that leaf node.

    Output tab delimited data on the average to a text file.
    Outputs the normalised average histogram bins and heights.

    \param s the name of the file to send the output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputAverageToTxtTabs(const std::string& s,
                                bool confirm = false) const;
    

    /*! \brief Output the collated information to a txt file.

    Output tab delimited data on the collation to a text file.

    \param s the name of the file to send the output to.
    \param confirm is a boolean controlling whether confirmation goes to
    console output (defaults to false).
    */
    void outputToTxtTabs(const std::string& s,
                            bool confirm = false) const;

	/*! \brief Add current state of collation to a log file.

    \param s is the name of the file to log to.
    \param i is a number representing the index of this state in a sequence.
    */
    void publicOutputLog(const std::string& s, const int i) const;
	 
	 //this was originally here - but removed - and now inserted again
	 /*! \brief Output the accumulated data over the collation to a txt file.

    Output tab delimited data on the average to a text file.
    Outputs the accumulation over the collation of data summarised in this.

    \param s the name of the file to send the output to.
    */
    void outputAccumulationToTxtTabs(const std::string& s) const;

    /*! \brief Get the estimated log likelihood from RSSample. 
	 */
    real getEstLogLikelihoodFromRSSample(RSSample& labSampledData,
										double dx, double wt)
	{	double WeightHist = 1.0;
		std::map<rvector, double, std::less<rvector> > WeightsPM;
		return getEstLogLikelihoodFromRSSample(labSampledData,
										dx, wt, WeightHist, WeightsPM); }
	
	real getEstLogLikelihoodFromRSSample(RSSample& labSampledData,
										double dx, double wt,  double WeightHist,
						std::map<rvector, double, std::less<rvector> >& WeightsPM);

	/*! \brief Make a marginalised version of this histogram collator.

	Marginalises to take out the unwanted dimensions and adjust summaries
	so that overall sum of (node vol x accumulated summaries) 
	is the same as for this.
	
	\reqDims is a vector of the dimensions to include in marginal.
	\return An AdaptiveHistogramCollator managing a subpaving which is
	the marginalised version of the subpaving managed by this.
	\pre \a reqDims must be compatible with current dimensions.
	\note allowed dimensions start at 1, ie dimensions to
	marginalise on can include 1, 2, ... #dimensions of this.
	\post returned histogram will have one summary value for each node and 
	have sum of 
	(node vol x accumulated summaries) = that for this.
	*/
    AdaptiveHistogramCollator marginalise(
								const std::vector<int>& reqDims) const;
								
	/*! \brief Find the coverage value for boxes.
	*/
	void findDensityRegion(double cov, double weightPM, vector<CollatorSPnode*> & covNodes,
									string covFileName);

	/*! \brief get leaf levels string
	*/
   std::string getLeafLevelsString() const;

	// not in trunk
    /*! \brief Get the sum of the variances for an area-related scalar summary.

    Treats the collation as representing a sample of histograms.

    The scalar summary is defined so that the variance for a histogram in the
    sample is the square of the sum of the areas of difference between
    a histogram sample and the average histogram over the sample.

    'Area', for a particular box, is the height of the bin multiplied by the
    box volume.

    \return The sum of the variances of the area-related scalar summary
    over each histogram collated.
    */
    real getSumVarianceAreaScalar() const;

	// not in trunk
    /*! \brief Get the sample variance for an area-related scalar summary.

    Treats the collation as representing a sample of histograms.

    The scalar summary is defined so that the variance for a histogram in the
    sample is the square of the sum of the areas of difference between
    a histogram sample and the average histogram over the sample.

    'Area', for a particular box, is the height of the bin multiplied by the
    box volume.

    The sample variance is:
    (sum of the variances over each histogram collated)/(number collated - 1).

    The sample variance is only defined in the number collated is > 1.

    \return The sample variance of the area-related scalar summary
    over the sample held in the collation.
    */
    real getSampleVarianceAreaScalar() const;

	// not in trunk
    /*! \brief Get the sum of the variances for scalar summary total height.

    Treats the collation as representing a sample of histograms.

    The scalar summary is the sum of the heights of the bins over the histogram.

    \return The sum of the variances of the total height scalar summary
    over each histogram collated.
    */
    real getSumVarianceTotalHeightScalar() const;

	// not in trunk
    /*! \brief Get the sample variance for scalar summary total height.

    Treats the collation as representing a sample of histograms.

    The scalar summary is the sum of the heights of the bins over the histogram.

    The sample variance is:
    (sum of the variances over each histogram collated)/(number collated - 1).

    The sample variance is only defined in the number collated is > 1.

    \return The sample variance of the total height scalar summary
    over the sample held in the collation.
    */
    real getSampleVarianceTotalHeightScalar() const;
    
		// Jenny addition for Gloria's convergence work
		// make a collator that is the differences of each element in a sample to the sample average
		AdaptiveHistogramCollator makeDifferencesToAverage() const;

		// Jenny addition for Gloria's convergence work
		// take a container and return the same container, which has been
		// cleared (if necessary) and re-filled with 
		// L1-distances-to-average, one for each histogram in collation
		RealVec& getL1DistancesToAverage(RealVec& container) const;
		
		/** @name Get the IAE of a distribution
		Get the integrated absolute error of the specified distribution.
		\return the integrated absolute error for this realization
		*/
		//@{
		/** Get the IAE for a uniform (mixture) distribution. */ 
		cxsc::real getUnifIAE(AdaptiveHistogram& myPart, std::vector<int> holesLoc, 
									double weight);
		
		/** Get the IAE for the standard uniform distribution. 
		 * Note: later could add arguments for specific boundaries. */
		 cxsc::real getUnifIAE();
		 
		 /** Get the IAE for a finite gaussian mixture distribution using interval 
		* techniques. */
	 cxsc::interval getFinMixIntervalIAE(FinMix& mixt, double tol, int deg);
	
		/** Get the IAE for a laplace distribution with mu=0 and b=1 using interval 
		* techniques. */
		cxsc::interval getLaplaceIntervalIAE(double tol, int deg);
	 
	 	/** Get the IAE for a laplace distribution with mu=0 and b=1 using interval 
		* techniques. */
		cxsc::interval getLognormalIntervalIAE(double tol, int deg);
 
		//@}

		//gloria addition
		/*! \brief Get the IAE between an AdaptiveHistogramCollator and a RealMappedSPnode.
		 */
		//cxsc::real getMappedIAE(RealMappedSPnode& nodeEst, ivector pavingBox) const;
		
			/** @name Get the Scheffe Tournament Winner
	    Exhaustive and bisection search.
		\param listScheffeSet the list of set of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecScheffeSet the current vector of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecWinnerVec the updated vector of vectors of tournament winners
      \param vecDeltaWinnerVec the updated vector of vectors of delta of winners		
		 
	*/
	//@{
	/** Get the Scheffe set from subpavings */
	void getHistScheffeSet(
	std::vector< std::vector< std::set<CollatorSPnode*, std::less<CollatorSPnode*> > > > & vecScheffeSetVec);
	
	void getHistYatSet(
	std::vector< std::set<CollatorSPnode*, std::less<CollatorSPnode*> > > & vecYatSet);
	
	void getHistScheffeWinner( 
	std::vector< std::vector< std::set<CollatorSPnode*, std::less<CollatorSPnode*> > > > & vecScheffeSetVec, 
	std::vector< std::vector<int> > & vecWinnerVec,
	std::vector< std::vector<double> > & vecDeltaWinnerVec);
	
	 double getNodesDelta(
				std::set<CollatorSPnode*, std::less<CollatorSPnode*> > & YatSet, 
				int thisTheta, size_t sizeColl);
	
	double getNodesMaxDelta(
			vector< set<CollatorSPnode*, less<CollatorSPnode*> > > & vecYatSet, 
			int thisTheta);
				
	void getMinDistEst(vector<double> & maxDelta, 	vector< set<CollatorSPnode*, less<CollatorSPnode*> > > & vecYatSet);

	//@}
		
};

// non-member functions
/*! \brief Output the contents of an AdaptiveHistogramCollator object.

Verbose output for an AdaptiveHistogram object.
*/

std::ostream & operator<<(std::ostream &os,
            const AdaptiveHistogramCollator& adhc);

} // end of namespace subpavings

#endif



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

/*! \file      adaptivehistogramvcollator.hpp
\brief AdaptiveHistogramVCollator declarations.
*/

#ifndef ___ADAPTIVEHISTVCOLL_HPP__
#define ___ADAPTIVEHISTVCOLL_HPP__

#include "adaptivehistogram.hpp"
#include "adaptivehistogramvalidation.hpp"
#include "sptypes.hpp"
#include "histevalobjval.hpp"
#include <set>

namespace subpavings {

/*! \brief A wrapper or manager for a CollatorSPVnode.

AdaptiveHistogramVCollator class objects manage 
\link subpavings::CollatorSPVnode CollatorSPVnode \endlink objects for the 
purpose of collating information from a number of AdaptiveHistogramVal objects.

The AdaptiveHistogramVCollator's CollatorSPVnode tree represents the subpaving
that is the union of all the subpavings associated with each 
AdaptiveHistogramVal in the collation.  Each node in the tree has two data 
members: a container structure holding one value for each collated histogram's 
training data and a double corresponding to the empirical measure of the 
histogram's validation data.  For a collation of histograms the container holds
the <b>normalised height</b> for each collated histogram of the histogram bin 
represented by the box of that node.

(The normalised height associated with a bin which is represented by the box
of a leaf node of a tree managed by an AdaptiveHistogramVal object is the number
of data points associated with that bin divided by (the total number of
data points in the histogram x the volume of the bin).  Thus the areas
(heights x volumes) of the bins sum to 1.

Since the tree represents the union of the subpavings associated with each
AdaptiveHistogramVal in the collation, the collation tree will have at least as
many and usually more bins than any of the collated histograms.

Each collated histogram will be represented in the summaries in the order in
which it was added to the collation.  Eg, the heights of the bins of the first
histogram to be collated will be first (index [1]) in the summary container.

If the collated AdaptiveHistogramVals have been properly formed and added to the
collation, the sum, over all the leaf nodes of the collation, of the volume of
the box of the leaf node multiplied by the values of the leaf node summary
corresponding a particular collated histogram will be 1 for each collated
histogram.

*/

//==============AdaptiveHistogramVCollator class declarations================//
class AdaptiveHistogramVCollator {

//----------------private members-----------------------------------------------
private:

    /*! \brief Private initialised constructor.

    Initialised  with pointer to subpaving.
    */
    explicit AdaptiveHistogramVCollator(CollatorSPVnode * spn);

    /*! \brief Pointer to the root CollatorSPVnode.

    A CollatorSPVnode is a binary tree representation of information from
    a number of subpavings.

    The summary information held by the CollatorSPVnode for an SPSVnode is the
    training data count/(total training data in tree * volume), ie the 
	 normalised histogram height for a bin corresponding the box of that SPSVnode.
	 
	 The Vemp information held by the CollatorSPVnode for an SPSVnode is the
	 validation data count/total validation data.
    */
    CollatorSPVnode * rootVCollator;

//----------------public members-----------------------------------------------
public:
    /*! \brief Default constructor.
    */
    explicit AdaptiveHistogramVCollator();

    /*! \brief initialised constructor.
    Initialised with an AdaptiveHistogramVal object.
    */
    explicit AdaptiveHistogramVCollator(const AdaptiveHistogramValidation& adh, 
                                        int whatSum);

    /*! \brief  Copy constructor
    */
    AdaptiveHistogramVCollator
        (const AdaptiveHistogramVCollator& other);
	
	/*! \brief Copy constructor for the subtracted ADHVC
	*/ 	  
	 AdaptiveHistogramVCollator(
			const AdaptiveHistogramVCollator& other, int toSubtract);	  

    /*! \brief Copy assignment operator.
    */
    AdaptiveHistogramVCollator&
        operator=(const AdaptiveHistogramVCollator& rhs);

    /*! Addition operator.

    Addition gives a histogram vcollator managing a tree which represents a
    subpaving which is the union of the subpavings represented by the operand
    collators.  The summary for each node in the tree contains all the values
    from the summaries of the corresponding nodes in the trees managed by
    the operand AdaptiveHistogramVCollators.

    */
    AdaptiveHistogramVCollator& operator+=(const
                AdaptiveHistogramVCollator& rhs);

    /*! Subtraction operator.

    Subtraction gives a histogram collator managing a tree which represents a
    subpaving which is the union of the subpavings represented by the operand
    collators.  The summary for each node in the tree contains all the values
    from the summary of the corresponding node in the trees managed by the left
    hand side operand and the negation of the values from the summary of the
    corresponding node in the trees managed by the right hand side operand.
    */
    AdaptiveHistogramVCollator operator-(const
                AdaptiveHistogramVCollator& rhs) const;


    /*! \brief Destructor.
    */
    ~AdaptiveHistogramVCollator();


    /*! An AdaptiveHistogramVCollator which is the average over this collation.

    Makes and returns an AdaptiveHistogramVCollator which is the average over
    the collation of histograms represented by this.  The tree managed
    by the average has structure exactly the same as the tree managed by this
    and one value in the summary for each node where that value is
    the average of the summary of the corresponding node in this.

    This method can only be performed on summaries which do not include negated

    */
    AdaptiveHistogramVCollator makeAverage() const;


    /*! \brief Return a pointer to the CollatorSPVnode this manages.
    */
    CollatorSPVnode* getSubPaving() const;

    /*! \brief Get the number of Adaptive Histogram objects collated.
    */
    size_t getNumberCollated() const;


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

    Outputs tab delimited data on the average to a text file.
    Outputs the normalised average histogram bins and heights.

    \param s the name of the file to send the output to.
    */
    void outputAverageToTxtTabs(const std::string& s) const;

    /*! \brief Output the accumulated data over the collation to a txt file.

    Output tab delimited data on the average to a text file.
    Outputs the accumulation over the collation of data summarised in this.

    \param s the name of the file to send the output to.
    */
    void outputAccumulationToTxtTabs(const std::string& s) const;
	 void outputDifferenceToTxtTabs(const std::string& s) const;


    /*! \brief Output the collated information to a txt file.

    Output tab delimited data on the collation to a text file.

    \param s the name of the file to send the output to.
    */
    void outputToTxtTabs(const std::string& s) const;
	 void outputToTxtTabs(const std::string& s, int whichColl) const;

    /*! \brief Add an AdaptiveHistogramVal object to the collation.		
    */
    void addToCollationWithVal(const AdaptiveHistogramValidation& adh, 
											int whatSum,
											size_t & agg);
	 
	 /*! \brief  Get a CollatorSPVnode pointer to the corresponding SPSVnode that was split 
	   
			\param chosenLargest the current splitted node
			\return splitCollNode the current splitted node as a CollatorSPVnode	  	
    */
	 void getSplitNodePtr(CollatorSPVnode* &splitCollNode, SPSVnode * const spn);
    
    /*! \brief Get the Yatracos set.
     
      The Yatracos set is obtained by pairwise comparisons of the heights of 
      the histogram at each leaf box and unions of leaf boxes for each theta 
		(number of splits). For all i, j pairs, the leaf box/unions of leaf boxes 
      at the i-th split goes into the Yatracos set if its height is higher than 
		the height at leaf box/unions of leaf boxes at the j-th split.
      
      The winning leaf boxes/union of leaf boxes are stored in a set of 
		\link subpavings::CollatorSPVnode CollatorSPVnode \endlink. The sets are 
		then stored in a sorted list. Any repetitive elements are removed from the 
      list by using the STL unique() algorithm.
      
		\param spiltCollNode the current splitted node
      \param current list of pointers to sets of 
				  \link subpavings::CollatorSPVnode CollatorSPVnode \endlink
      \param vecRowYatSet the current row of the growing Yatracos matrix
		\param vecColYatSet the current column of the growing Yatracos matrix
		
      \return an updated list of pointers to sets of 
					\link subpavings::CollatorSPVnode CollatorSPVnode \endlink
		\return an updated row and column vector of the current edges
   */

	void getYatracosClassAll(CollatorSPVnode * const splitCollNode,
	std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecRowYatSet,
	std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecColYatSet,
	std::list< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & listYatSet );
 
    /*! \brief Get the delta value for each Yatracos element at the specified 
	             split number.
   
      Get the delta value for each element of the list of pointers to the 
      Yatracos leaf boxes at the specified theta.
         
      The delta value is the absolute difference between the empirical measure 
      of the training data and the validation data.
      
      \param YatSet the set of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param theta the split number
      \param currentStage the most current split number
      \return the delta value		
    */
    double getNodesDelta(
				std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > & YatSet, 
				int thisTheta);
   double getNodesMaxDelta(
				std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecYatSet, 
				int thisTheta);


   /*! \brief Get the maximum delta values for all thetas for the Yatracos Class.
   
      Get the maximum delta value at each theta over the list of Yatracos leaf boxes for all thetas.
         
      The delta value is the absolute difference between the empirical measure 
      of the training data and the validation data.
      
      \param listYatSet the list of set of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecMaxDeltaVec the current vector of vectors to the maximum delta values for each theta 
      */
   void getYatracosDelta(std::list< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & listYatSet, 
   std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecRowYatSet, 
   std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecColYatSet, 
   std::vector< std::vector<double> > & vecMaxDeltaVec);
  
  
  /*! \brief Get the maximum delta values for all thetas for the Yatracos Class at the end.
   
      Get the maximum delta value at each theta over the list of Yatracos leaf boxes for all thetas.
         
      The delta value is the absolute difference between the empirical measure 
      of the training data and the validation data.
      
      \param listYatSet the list of set of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecMaxDeltaVec the current vector of vectors to the maximum delta values for each theta 
      */
   void getYatracosDeltaEnd(
   std::list< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & listYatSet, 
   std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecRowYatSet, 
   std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecColYatSet, 
   std::vector<double> & vecMaxDeltaVec);
  
  
	/*! \brief Get the minimum distance estimate from the MDE method. 
		\param method: 1. Scheffe tournament
							  2. Minimum distance estimate
		\param winner vector of winners
		\param loser  vector for losers
		\param deltaWinner vector for deltas of winners
		\param deltaLoser vector for deltas of losers
	*/
	void getBisectionSearchEstimate(int method, std::vector<size_t> & winner, 
	                            std::vector<size_t> & loser,
									    std::vector<double> & deltaWinner,
										 std::vector<double> & deltaLoser);	

   /*! \brief Get the thetas that gives the minimum distance.		
		Get the thetas that gives the minimum distance.
		Write definition for the min.  dist. theta  here.
		
      \param vecMinDistTheta the current vector of vectors of the thetas that gives the minimum distance for each theta 
      \param vecDeltaMaxVec vector of vectors of the maximum delta values
      \param n sample size
      \return vecMinDistTheta the updated vector of vectors of the thetas that gives the minimum distance for each theta
   */
   void getMinDistTheta(
   std::vector< std::vector<int> > & vecMinDistTheta, 
   std::vector< std::vector<double> > & vecDeltaMaxVec, int n);
	
   bool getMinDelta(int maxCheck, std::vector< std::vector<double> > & vecDeltaMaxVec);
	
  
   /*! \brief Get the infimum delta value for all thetas.
		
		Get the infimum delta value at each theta over the list of Yatracos leaf boxes for all thetas.
         
      This may also be used as a stopping critera.
      \param vecInfDelta the current vector of vectors to the infimum delta values for each theta 
      \param vecMaxDeltaVec a vector of vectors to the maximum delta values for each theta 
      \param n sample size
      \return vecInfDelta a vector of vectors to the infimum delta values for each theta
   */
  std::vector<double> getInfDelta(std::vector<double> & vecInfDelta, 
                   std::vector< std::vector<double> > & vecDeltaMaxVec, int n);

	/*! \brief Get the root box.	 
		Get the root box of the AdaptiveHistogramVCollator object.
	*/
	ivector getRootBox();
	
	
	/** @name Get the Scheffe Tournament Winner
	    Exhaustive and bisection search.
		\param listScheffeSet the list of set of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecScheffeSet the current vector of pointers to sets of \link subpavings::CollatorSPVnode CollatorSPVnode \endlink (pointers to the leaf boxes)
      \param vecWinnerVec the updated vector of vectors of tournament winners
      \param vecDeltaWinnerVec the updated vector of vectors of delta of winners		
		 
	*/
	//@{
	/** Get Scheffe Tournmanet Winner (exhaustive). */
	void getScheffeSetAll(CollatorSPVnode * const splitNode,
		vector< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & vecScheffeSet,
		list< set<CollatorSPVnode*, less<CollatorSPVnode*> > > & listScheffeSet);
		
   /** Get the winner and delta values for all pairwise comparisons for 
	 * the Scheffe set. */
   void getScheffeWinner( 
		std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecScheffeSet, 
		std::vector< std::vector<int> > & vecWinnerVec,
		std::vector< std::vector<double> > & vecDeltaWinnerVec);

	/** Get the Scheffe set from subpavings */
	void getHistScheffeSet(
	std::vector< std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > > & vecScheffeSetVec);
	
	void getHistYatSet(
	std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & vecYatSet);
	
	void getHistScheffeWinner( 
	std::vector< std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > > & vecScheffeSetVec, 
	std::vector< std::vector<int> > & vecWinnerVec,
	std::vector< std::vector<double> > & vecDeltaWinnerVec);
	
	void getMinDistEst(std::vector<double> & maxDelta, std::vector< std::set<CollatorSPVnode*, std::less<CollatorSPVnode*> > > & veYatSet);
	

	//@}
	
   // computes minimal subpavings from sibling subpavings
   // a subpaving is minimal if it has no sibling leaves
   // a minimal subpaving is created by discarding sibling leaves
   // and create summary data for new parent from children
   // warning: nodeReunite would not normally be used with
   // CollatorSPVnodes but is in the base class and is
   // reimplemented to try do it appropriately for this
   // derived class should it be needed.
   // This function is untested.
	void makeMinimal();	
	
										
	/*! \brief get leaf levels string
	*/
   std::string getLeafLevelsString() const;
   
   size_t getTotalNodes();

}; // end of declaring functions of AdaptiveHistogramVCollator class

//=========End of AdaptiveHistogramVCollator class declarations================//
} // end of namespace subpavings

//===============Non-member function declarations=============================// 
/*! \brief Output the contents of an AdaptiveHistogramVCollator object.

Verbose output for an AdaptiveHistogram object.
*/
std::ostream & operator<<(std::ostream &os,
            const subpavings::AdaptiveHistogramVCollator& adhc);

/*! \brief Find if double is negative.
*/
double isNegative(double d);

/*! \brief Output all boxes in collator to text file
*/
void outputAllNodesToTxtTabs(const std::string& s, const subpavings::AdaptiveHistogramVCollator& adhc); 
//===============Non-member function declarations=============================//



#endif

 

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
\brief AdaptiveHistogram declarations.
*/

#ifndef ___ADAPTIVEHIST_HPP__
#define ___ADAPTIVEHIST_HPP__

#include "piecewise_constant_function.hpp"
#include "realmappedspnode.hpp"
#include "spsnode.hpp"
#include "splitdecisionobj.hpp"
#include "nodecompobj.hpp"
#include "histmcmcobjs.hpp"
#include "histpenalty.hpp"
#include "histevalobj.hpp"
#include "spsnode_measure_obj.hpp"



// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

#include <gsl/gsl_rng.h>        // to know about the gsl random number generator

namespace subpavings {

// a class for comparison between spsnodes
class MyCompare
{
    const NodeCompObj& myNC;

    public:
    MyCompare(const NodeCompObj& nc);

    bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const;

};

/*! \brief A wrapper or manager for an SPSnode aka StatsSubPaving in
conjunction with massive amounts of sample data.

Here sample data is multi-dimensional point-valued data in a cxsc::rvector
container.  The %AdaptiveHistogram class
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

The %AdaptiveHistogram class uses the C-XSC library class rvector
for sample data points.  rvectors can have 1 or many dimensions.

Each %AdaptiveHistogram has an integer label.

@todo adjust references to holdAllStats etc to refer to getter methods.
   
@todo clean up nomenclature
*/

class AdaptiveHistogram {

    public:
	
	/*! \brief Types for evaluating priority split queues. */ 
	class PrioritySplitQueueEvaluator {
		
		public:
			PrioritySplitQueueEvaluator(SPSNodeMeasure& m, 
						real cs, size_t maxL);
						
			PrioritySplitQueueEvaluator(SPSNodeMeasure& m, 
						size_t cs, size_t maxL);
			
			PrioritySplitQueueEvaluator(SPSNodeMeasure& m, 
						int cs, size_t maxL);
			
			PrioritySplitQueueEvaluator(SPSNodeMeasure& m, size_t maxL);
			
			SPSNodeMeasure& getMeasurer() const;
			
			real getCritStop() const;
			
			size_t getMaxLeaves() const;
			
			bool getUsingCritStop() const;
			
			void setMaxLeaves(size_t ml);
			
			void setUsingCritStop(bool b);
			
			
			
		private:
			SPSNodeMeasure& measurer;
			real critStop;
			size_t maxLeaves;
			bool usingCritStop;
		
	};
	
	/*! \brief Type for collecting change of state information. */ 
	class ChangeOfStateInformation {
		
		friend class AdaptiveHistogram;
		
		public:
		
			virtual ~ChangeOfStateInformation() {}
		
		private:	
			/*! \brief Notify this of change in log-posterior from last change. */			
			virtual void notifyDeltaPi(real dp) = 0;
			/*! \brief Notify this of a split. */
			virtual void notifySplit(const SPnode * const spn) = 0;
			/*! \brief Notify this of a merge. */
			virtual void notifyMerge(const SPnode * const spn) = 0;
		
	};

    /*! \brief Default constructor

    By default, only counts are maintained in subpaving this manages, rather
    than all available stats.  The default label is 0;
    */
    AdaptiveHistogram();

    /*! \brief Initialised constructor.

    \param as indicator controlling whether all available
    statistics be maintained in the SPSnode tree managed by this
    AdaptiveHistogram (true for all stats, false for counts only).
	\param lab value for the label for this (defaults to 0).
    */
    explicit AdaptiveHistogram  (bool as, int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with domain box.
    By default, only counts are maintained as stats in the SPSnode tree
    managed by this.  Default label is 0.
    
    Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

    Ideal constructor when the support domain of data is known a priori
    or has been transformed to a known domain but splitting criteria have
    not been determined a priori.
    */
    explicit AdaptiveHistogram(const ivector& v, bool as = false);

	/*! \brief Initialised constructor.

    Initialised  with domain box, label and indicator for whether 
	all stats are maintained in the in the SPSnode tree
    managed by this or only counts.
    
    Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

    Ideal constructor when the support domain of data is known a priori
    or has been transformed to a known domain but splitting criteria have
    not been determined a priori, and a specific label is to be specified.
    */
    AdaptiveHistogram(const ivector& v, bool as, int lab);

    /*! \brief  Copy constructor.
    */
    AdaptiveHistogram(const AdaptiveHistogram& other);

    /*! \brief Copy assignment operator.
    */
    AdaptiveHistogram& operator=(AdaptiveHistogram rhs);

    /*! \brief Overloaded add-to-self operator.

    This becomes the result of adding this and \a rhs together using a 
    non-minimal union of subpavings.
	
    If before the operation, one of either this or \a rhs have no subpaving,
    after the operation this will have the non-null subpaving.  If  
	before the operation, both this or rhs have no subpaving,
    after the operation this will have no subpaving. 
     
	holdAllStats will be unchanged as a result of this operation.
	label will be unchanged as a result of this operation.
    
    \param rhs the histogram to add to this.
    \return the result of adding \a rhs to this.
	\pre If this and \a rhs both have a subpaving with a box, then those
	boxes must have equal dimensions and box side lengths.
	\post This will have a subpaving that is the non-minimal 
	union of the orignal subpavings of this and \a rhs,
    and the data collection will be the union of the original
    data collection and the data collection of \a rhs, and the 
    the counts in each part of the subpaving will be correct, and
    (if held) the optional statistics can also calculated correctly.
    holdAllStats will be unchanged.
    */
	AdaptiveHistogram& operator+=(const AdaptiveHistogram& rhs);

    /*! \brief Overloaded addition operator.

    Return the result of adding another histogram to this.
    	
    If before the operation, one of either this or \a rhs have no subpaving,
    the result will have the non-null subpaving.  If  
	before the operation, both this and \a rhs have no subpaving,
    the result will have no subpaving.
     
	\param rhs the histogram to add to this.
    
    \return the result of adding this  and \a rhs together.
	\pre If this and \a rhs both have a subpaving with a box, then those
	boxes must have equal dimensions and box side lengths.
	\post The histogram returned have a subpaving that is the 
    non-minimal union of the subpavings 
    of this and \a rhs and a data collection that is the 
    union of the data collections of this and \a rhs.
    The holdAllStats indicator of the returned histogram 
    will be the logical OR of the holdAllStats indicators
	of this and \a rhs.  The label of the returned histogram
	will be the label of this before the operation if that
	is the same as the label of \a rhs, and but will be 0 if
	this and \a rhs had different labels before the operation.*/
    const AdaptiveHistogram operator+(const AdaptiveHistogram& rhs) const;


    //! Destructor
    ~AdaptiveHistogram();



    /*! \brief Return a pointer to the SPSnode this manages.
    \todo This is bad bad bad and should not happen.  Change when 
	need for it for testing is over...
	*/
    SPSnode* getSubPaving() const;
	
	/*! \brief Return the label for this.
	*/
    int getLabel() const;
	
	/*! \brief Return a copy of the data collection this manages.
    */
    BigDataCollection getDataCollection() const;
	
	
    /*! \brief get the value of holdAllStats field.

    This determines whether the histrogram's rootPaving will maintain all
    available stats (true) or just the counts (false).
    */
    bool getHoldAllStats() const;

    /*! \brief Get whether this has a subpaving to manage.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return true if this has a subpaving to manage.
	false otherwise.*/
    bool hasSubPaving() const;
	
	/*! \brief Get the box of the subpaving managed by this.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return copy of the box of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
	cxsc::ivector getRootBox() const;
	
	/*! \brief get the volume of root box of the subpaving this manages.

    \return volume of the root paving box.
    \pre hasSubPaving() == true.*/
	real getRootVolume() const;
   
   
	/*! \brief get the dimensions of the subpaving this manages.

    \return 0 if this does not have a subpaving, else returns the
    dimensions of the subpaving.*/
    int getDimensions() const;
    
    /*! \brief Gets sample mean for data held.
	
	This calculates the mean for the
	data and returns it as a d*d-dimensional rvector, 
	where d is the dimension of the data.

	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return If means are held (see holdAllStats) and there is at 
	least one data point,
	then return the mean. Otherwise return an rvector
	of cxsc::SignalingNaN values.
	\pre hasSubPaving() == true.	     */
    rvector getRootPavingMean() const;


    /*! \brief Gets sample variance covariance vector for data held.
	
	This calculates the sample variance-covariance matrix for the
	data and returns it as a d*d-dimensional vector of reals, 
	where d is the dimension of the data.

	cov(i,j) is at index i*d+j in the returned RealVec (indices 
	from 0 to d*d-1, where d is the dimension of the data).
	
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return If variance-covariances are held (see holdAllStats) 
	and there are at least two data points (see getRootCounter())
	then return the sample variance-covariance matrix as a vector.
	Otherwise return a RealVec of cxsc::SignalingNaN values.
	\pre hasSubPaving() == true.	     */
    RealVec getRootPavingVarCovar() const;
	
	/*! \brief Gets sample mean for data held in the dataCollection.
	
	This calculates the maean for the
	data in the data collection and returns it as a d*d-dimensional rvector, 
	where d is the dimension of the data.
	* 
	If means are held (see holdAllStats) 
	then this should give the same result as getMean(),
	but if this does not hold all stats then getMean() will give a
	collection of NaN values whereas this method recalculates
	the variance-covariances directly from the data collection.


	\return If there is at 
	least one data point,
	then return the mean. Otherwise return an rvector
	of cxsc::SignalingNaN values.*/
    rvector getDataCollectionMean() const;


    /*! \brief Gets sample variance covariance vector for data held
	 * in the data collection.
	
	This calculates the sample variance-covariance matrix for the
	data in the data collection
	and returns it as a d*d-dimensional vector of reals, 
	where d is the dimension of the data.
	
	If variance-covariances are held (see holdAllStats) 
	then this should give the same result as getVarCov(),
	but if this does not hold all stats then getVarCov() will give a
	collection of NaN values whereas this method recalculates
	the variance-covariances directly from the data collection.

	cov(i,j) is at index i*d+j in the returned RealVec (indices 
	from 0 to d*d-1, where d is the dimension of the data).
	
	\return If there are at least two data points
	then return the sample variance-covariance matrix as a vector.
	Otherwise return a RealVec of cxsc::SignalingNaN values.*/
    RealVec getDataCollectionVarCovar() const;

    /*! \brief Gets count in the rootpaving in the root paving.
     
    Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer.
	
	\return the total count of data in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    size_t getRootCounter() const;

    /*! \brief Gets number of leaf nodes in the root paving.
    
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return the total number of leaves in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    size_t getRootLeaves() const;
	
	/*! \brief Gets number of cherry nodes in the root paving.
    
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return the total number of cherries in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    size_t getRootCherries() const;
	
	
	/*! \brief Gets the total depth over all the leaves in
	the root.
    
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return the total depth of of all the leaves in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    unsigned long int getRootTotalLeafDepth() const;

    /*! \brief Gets the sum of leaf count over volume in root paving.
    
    Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer.  
	
	\return the sum over all leaves of the subpaving managed of the 
	leaf counter over the leaf volume.
    \pre hasSubPaving() == true.	     */
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
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
	\param prec the precision for output formatting. ie, number
	of decimal places.    */
    void outputLog(const std::string& s, int i, int prec = 5) const;
	
	/*! \brief Append current state of histogram to a txt log file.

    Format is a tab-delimited file of numeric data.
    Output is plain: just vols, counters, and boxes.

    \param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
	\param prec the precision for output formatting. ie, number
	of decimal places.    */
	void outputLogPlain(const std::string& s, 
										int i, int prec = 5) const;

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
	\return true if the point was successfully associated with one of 
	the rootPaving's leaves, false otherwise (ie is outside the 
	bounds of the box the paving represents).
	
	Throws the following exceptions:
	<ul>
	<li>Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.</li>
	<li>Throws a NoBox_Error if the subpaving that this manages
	has no box.</li>
	<li>Throws an IncompatibleDimensions_Error if the \a newdata has 
	different dimensions to the subpaving that this manages.</li>
	</ul>

	\return true if the data point is within the box of the subpaving
	managed by this.
    \pre this must have a rootPaving pointing to
    an SPSnode with a box; this SPSnode may already have data associated
    with it.
    \post If return value is true the data point \a newdata 
    will have been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data point will have been incremented by one, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).
    */
    bool insertOne(rvector newdata,  const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging=NOLOG);


    /** @name Insert rvectors in a txt file into AdaptiveHistogram object.

    Reads in lines of data representing rvectors from a txt file.
    The dimensions of the data may given or may be deduced 
	from existing data in this or from the data in the file.
	
	Expects one line per rvector with the elements separated by 
	white space (space or tabs), or by commas.
	
	Input files can be checked to try to ensure data conforms with 
	some given dimensions and/or has no illegal characters.  
	* 
	Illegal characters are any characters other than those in the string
	"eE+-.,0123456789 \t".
	* 
	Note that this allows numbers in scientific format, and with space,
	tab, or comma delimiters.  However, allowing the characters 'e'
	and 'E' as part of scientific format, and the delimiters,
	also allows the possibility that these are erroneously present.  
	The level of checking for illegal characters is not rigorous and will
	not pick up the inadvertent or erroneous presence of these non-numeric
	characters.  This includes the presence of delimiters before the first 
	numeric character, which may result in 0.0 being inserted into the
	first dimension(s) of the rvector read from the relevant data line.
	* 
	<strong>Ordinary checking</strong> means that:<ul>
	<li>The first valid line is checked against the 
	given data dimensions if any.  Oherwise the first valid data line
	is used to find data dimensions.</li>
	<li> If data dimensions are given and the dimensions found in the 
	first valid data line are less than those given dimensions, an error
	message is given and reading of all data is aborted. Any data 
	in excess of the given dimensions will however just be ignored.    
	<li> The first valid data line, as described above, is the 
	first line after the given headers not containing illegal characters.
	All lines after that are also checked to ensure that they do not 
	contain illegal characters.    
	</ul>
	Note that under ordinary checking,   data lines after the first 
	valid line are not checked to ensure that they conform to 
	expected dimensions: data in excess of the given dimensions will 
	be ignored and 'missing' dimensions will be padded with 0.0's.    
	 
	<strong>Paranoid checking</strong> means that expected 
	data dimensions must be given and checking is carried out as for 
	Ordinary checking plus further checking:<ul>
	<li> All lines after the first valid data line are also checked to
	to check that the number of blocks of numbers is equal to the
	expected dimensions.  </li> 
	</ul>
	* 
	<strong>Fast checking</strong> effectively means no checking (fast
	reading of data).  Data is assumed to conform to the given dimensions:
	data in excess of the given dimensions will be ignored and
	all data lines which are less than the given dimensions will result in
	an rvector "padded" with 0.0 in the additional dimensions.  
	Illegal characters in the input, including initial delimiters before
	the first numeric values, may be interpreted as part of a
	real vector, or as delimiters, or ignored.  No guarantees are given 
	about the resulting data read in if the data does not conform to the
	values given for the parameters for headerlines and data dimensions
	or if the data is corrupted by illegal characters. 
	 
	While reading of the file continues (has not been aborted, for 
	example because the dimension of the data does not seem to be as 
	expected) , input lines which do not pass checks are rejected and logged
	but the rest of the file will continue to be processed.
	
    Blank lines are ignored.

    If no data dimension is specified by the user (\a dim = -1), 
	then the expected dimension \a dim will be taken as the dimension of the
	first valid data line read in from the file.

	If this has no subpaving to manage, then a subpaving to manage will
	be created by the method.  The subpaving created will be of the
	same dimensions as the expected dimensions of the data and will
	have a box big enough to contain all the data of the expected
	dimensions.
	
	Throws an invalid_argument exception if \a dim is supplied by the
	user and \a dim < 1.
	
	For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    3-dimensional rvector, a string "12.04 1.00005E-10 -30.0006" will be read as a
    3-dimensional rvector, a string "-30.0006" will be read as a
    1-dimensional rvector and a string "30 20" will be will be read as a
    2-dimensional rvector 

	\param s the name of the txt file to read data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
	\param headerlines is number of headerlines to skip before reading 
	data.
	\param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)
    \param dim expected data dimension.  
	
   \return true if at least some data was read in, else false (if
    the file could not be opened, if the file contained no data, 
    or if no valid data was found in the file).
    \pre A file with filename s (s can include path and name), if 
	\a dim is supplied then \a dim >= 1. 
    \pre If hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If this has no subpaving to manage 
	prior to the insertion, a  subpaving with an appropriate box will
	have been made to fit all the data of the expected dimensions.*/
   //@{

    /*! \brief Ordinary level checking, no splitting, 1-dimensional data. 

	\deprecated Use insertRvectorsFromTxtOrd with \a dim = 1 or one of the
	other insertRvectorsFromTxt variants. */
    bool insertOneDimDataFromTxt(const std::string& s,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging = NOLOG);

    /*! \brief Ordinary level checking, adaptive splitting 
	 * with each data point inserted, 1-dimensional data. 

	\deprecated Use insertRvectorsFromTxtOrd with \a dim = 1 or one of the
	other insertRvectorsFromTxt variants. */
    bool insertOneDimDataFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging = NOLOG);
 
	/*! \brief Ordinary level checking, no splitting, no dimensions specified. 
	
	 \deprecated Use insertRvectorsFromTxtOrd or one of the
	other insertRvectorsFromTxt variants. */
    bool insertRvectorsFromTxt(const std::string& s,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);

    /*! \brief Ordinary level checking, adaptive splitting with each data point 
	 * inserted, no dimensions specified. 

	 \deprecated Use insertRvectorsFromTxtOrd or one of the
	other insertRvectorsFromTxt variants. */
     bool insertRvectorsFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);
    
	/*! \brief Ordinary level checking, no splitting, dimensions specified. 
    
	\deprecated Use insertRvectorsFromTxtOrd or one of the
	other insertRvectorsFromTxt variants. */
    bool insertRvectorsFromTxt(const std::string& s,
								int dim,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);


    /*! \brief Ordinary level checking, adaptive splitting with each data point 
	 * inserted, dimensions specified. 
	
	 \deprecated Use insertRvectorsFromTxtOrd or one of the
	other insertRvectorsFromTxt variants. */
    bool insertRvectorsFromTxt(const std::string& s,
                    const SplitDecisionObj& boolTest,
					int dim,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);

	/*! \brief Ordinary level checking, no splitting, no dimensions specified. */
	bool insertRvectorsFromTxtOrd(const std::string& s,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);


    /*! \brief Ordinary level checking, adaptive splitting with each data point 
	 * inserted, no dimensions specified. */
    bool insertRvectorsFromTxtOrd(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);
    
	/*! \brief Ordinary level checking, no splitting, dimensions specified. */
    bool insertRvectorsFromTxtOrd(const std::string& s,
								int dim,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);


    /*! \brief Ordinary level checking, adaptive splitting with each data point 
	 * inserted, dimensions specified. */
    bool insertRvectorsFromTxtOrd(const std::string& s,
                    const SplitDecisionObj& boolTest,
					int dim,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);

	//new
	/*! \brief Ordinary level checking, no splitting, selected required dimensions. */
    bool insertRvectorsFromTxtOrd(const std::string& s,
								const std::vector < int >& reqDims,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);

	//new
    /*! \brief Ordinary level checking, adaptive splitting with each data point 
	 * inserted, , selected required dimensions. */
    bool insertRvectorsFromTxtOrd(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::vector < int >& reqDims,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);

	// paranoid level checking
	/*! \brief Paranoid level checking, no splitting, dimensions specified. */
    bool insertRvectorsFromTxtParanoid(const std::string& s,
								int dim,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);


    /*! \brief Paranoid level checking, adaptive splitting with each data point 
	 * inserted, dimensions specified. */
    bool insertRvectorsFromTxtParanoid(const std::string& s,
                    const SplitDecisionObj& boolTest,
					int dim,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);
	//new
 	/*! \brief Paranoid level checking, no splitting, selected required dimensions. */
    bool insertRvectorsFromTxtParanoid(const std::string& s,
								const std::vector < int>& reqDims,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);

	//new
    /*! \brief Paranoid level checking, adaptive splitting with each data point 
	 * inserted, selected required dimensions. */
    bool insertRvectorsFromTxtParanoid(const std::string& s,
                    const SplitDecisionObj& boolTest,
					const std::vector < int>& reqDims,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);
 
	// fast - minimal - level checking
	/*! \brief Fast level checking, no splitting, dimensions specified. */
    bool insertRvectorsFromTxtFast(const std::string& s,
								int dim,
								const std::size_t headerlines = 0,
                                LOGGING_LEVEL logging=NOLOG);


    /*! \brief Fast level checking, adaptive splitting with each data point 
	 * inserted, dimensions specified. */
    bool insertRvectorsFromTxtFast(const std::string& s,
                    const SplitDecisionObj& boolTest,
					int dim,
					const std::size_t headerlines = 0,
                    LOGGING_LEVEL logging=NOLOG);

    
    //@}
	
	
	/** @name Insert rvectors from a vector of vectors of doubles
	into AdaptiveHistogram object.

	Reads in vectors of doubles, representing rvectors, from a std::vector.
    The dimensions of the rvector can be given or (if \a dim is not given) 
	deduced from from the 
	first vector in the input vector \a inputData.
   
    The rest of the data then is expected to have the same dimension as
	that given or found from the frst vector in the input vector.  
	Any data with dimensions less than the expected dimensions will be 
	rejected.  If data has more than the expected dimensions, only the
	values on the first expected dimensions will be read in (the rest 
	will be ignored).
	   
    Expects one inner vector per rvector.

    Throws an invalid_argument exception if \a dim is supplied by the
	user and \a dim < 1.
	
	\param inputData is the collection of vectors of doubles to be input.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
	\param dim is the dimensions of the data to expect.  If \a dim = -1 then
	the method will try to establish the dimensions from the first 
	data point in the \a inputData.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)

    \return true if at least some data was input, else false (if
    no valid data was found in the \a inputData).
    \pre If \a dim is supplied then \a dim >= 1. 
    \pre If hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If this has no subpaving to manage 
	prior to the insertion, a  subpaving with an appropriate box will
	have been made to fit all the data of the expected dimensions.*/

   //@{

    /*! \brief All rvectors are associated with the root paving, no splitting. */
    bool insertRvectorsFromVectorOfVecDbls(
								const std::vector < VecDbl >& inputData,
								LOGGING_LEVEL logging=NOLOG);

    /*! \brief Adaptive splitting with each data point inserted. */
    bool insertRvectorsFromVectorOfVecDbls(
					const std::vector < VecDbl > & inputData,
                    const SplitDecisionObj& boolTest,
					LOGGING_LEVEL logging=NOLOG);
    
	/*! \brief All rvectors are associated with the root paving, no splitting, 
	\a dims given */
    bool insertRvectorsFromVectorOfVecDbls(
								const std::vector < VecDbl >& inputData,
								int dim,
								LOGGING_LEVEL logging=NOLOG);


    /*! \brief Adaptive splitting with each data point inserted, 
	\a dims given. */
    bool insertRvectorsFromVectorOfVecDbls(
					const std::vector < VecDbl > & inputData,
                    const SplitDecisionObj& boolTest,
					int dim,
					LOGGING_LEVEL logging=NOLOG);
   
	//@}
    
    
    /** @name Insert rvectors from a container of rvectors.
	
	If \a checkDims == false, appends the entire contents of
	\a rvec to the end of \a data unless the first element in \a rvec
	is an empty (0-d) rvector.
	
	If \a checkDims == true, checks the dimensions of each 
	element in \a inputData and
	only inserts those elements where with dimension equal to
	that of the first element in \a inputData .
	
	Throws an illegal_argument exception if the first element
	in \a inputData is an empty rvector.
	
	\param inputData is the container of rvectors we want to transfer
	to \a theData.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param checkDims is an indicator for whether the dimensions
	each of the elements of \a inputData are compared to the dimensions
	of the first (true) or not (false). 
	\param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    
	\return true if at least some data was input, else false (if
    no valid data was found in the \a inputData).
    \pre If hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If this has no subpaving to manage 
	prior to the insertion, a  subpaving with an appropriate box will
	have been made to fit all the data of the expected dimensions.*/

    //@{

    /*! All rvectors are associated with the root paving, no splitting. */
    bool insertFromRVec(const RVecData& inputData, 
							LOGGING_LEVEL logging=NOLOG);

    /*! Adaptive splitting with each data point inserted. */
    bool insertFromRVec(const RVecData& inputData, 
							const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging=NOLOG);

	/*! All rvectors are associated with the root paving, no splitting
	\a checkDims specified. */
    bool insertFromRVec(const RVecData& inputData, 
							bool checkDims, 
							LOGGING_LEVEL logging=NOLOG);

    /*! Adaptive splitting with each data point inserted,
	\a checkDims specified */
    bool insertFromRVec(const RVecData& inputData, 
							const SplitDecisionObj& boolTest,
							bool checkDims,
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
    \pre if hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If the subpaving managed had no 
	box prior to the insertion, a box will have been made to fit all
	the data.*/
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

    Insert rvectors from the labeled point samples of an RSSample
    object.
	* 
	If a label parameter is not specified, insert only points from
	sample where the point label matches the label for this.  Otherwise
	insert only points where the point label matches \a label.

    \param rss the RSSample object to get data from.
	\param label the label to use to select labelled points from \a rss 
	to insert.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file.
    \pre if hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If the subpaving managed had no 
	box prior to the insertion, a box will have been made to fit all
	the data.*/
   //@{

    /** All rvectors are associated with the root paving, no splitting,
	 inserts points with label \a lab. */
	bool insertFromRSSample(const RSSample& rss, int lab,
                LOGGING_LEVEL logging)
    {
        SplitNever sn; // a dummy split decision object
        return insertFromRSSample(rss, lab, sn, logging);
    }
	
	 /** All rvectors are associated with the root paving, no splitting,
	 inserts only points with same label as this. */
    bool insertFromRSSample(const RSSample& rss,
                LOGGING_LEVEL logging)
    {
        int lab = label; // label for this
		SplitNever sn; // a dummy split decision object
        return insertFromRSSample(rss, lab, sn, logging);
    }

    /** Adaptive splitting with each data point inserted,
	 inserts points with label \a lab. */
    bool insertFromRSSample(const RSSample& rss, int lab,
                            const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging);
	
	/** Adaptive splitting with each data point inserted,
	 inserts only points with same label as this. */
	bool insertFromRSSample(const RSSample& rss,
                            const SplitDecisionObj& boolTest,
                            LOGGING_LEVEL logging)
	{
        int lab = label; // label for this
		return insertFromRSSample(rss, lab, boolTest, logging);
    }
    //@}


    /** @name Insert a set number of rvectors from an RSSample object.

    Insert a set number of rvectors from the labeled point samples of an
    RSSample object.  Sampling is random with replacement.

	If a label parameter is not specified, insert only points from
	sample where the point label matches the label for this.  Otherwise
	insert only points where the point label matches \a label.

    The overloaded versions which take a random number generator as a parameter
    can be used to take successive samples from the same RSSample object.
    Otherwise the random number generator is created and destroyed during
    the scope of the function and repeating the identical function call will
    produce an identical sample from the RSSample object.

    \param samplesize the size of the sample to draw.
    \param label the label to use to select labelled points from \a rss 
	to insert.
    \param gsl_rng * rgsl is a random number generator.
    \param seed is a seed for a random number generator.
    \param rss the RSSample object to get data from.
    \param boolTest is a reference to an object providing a function
    operator determining whether to split a node when a data point arrives.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \pre if hasSubPaving() == true, then the dimensions of the data
	must match the dimensions of the root paving box.  
    \post If return value is true at least one data point will have 
    been put into the AdaptiveHistogram's dataCollection and
    also associated with one of the rootPaving's leaves via an iterator
    to dataCollection, the counts of all parts of the paving containing
	the data points added will have been incremented by the number of 
	new data points they contain, and values
	held to calculate other statistics will have been adjusted if 
	maintained (see holdAllStats).  If the subpaving managed had no 
	box prior to the insertion, a box will have been made to fit all
	the data.*/
   //@{

    /** All rvectors are associated with the root paving, no spliting,
    random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize,
			gsl_rng * rgsl,
            const RSSample& rss, int lab, LOGGING_LEVEL logging)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rgsl, rss, lab, sn,
                                        logging);
    }
	 /** All rvectors are associated with the root paving, no spliting,
    random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize, gsl_rng * rgsl,
            const RSSample& rss, LOGGING_LEVEL logging)
    {
		int lab = label; // label for this
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rgsl, rss, lab, sn,
                                        logging);
    }
    /** All rvectors are associated with the root paving, no spliting,
    seed for creating a random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, int lab, LOGGING_LEVEL logging)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, seed, rss, lab, sn,
                                        logging);
    }
	/** All rvectors are associated with the root paving, no spliting,
    seed for creating a random number generator supplied. */
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, LOGGING_LEVEL logging)
    {
		int lab = label; // label for this
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, seed, rss, lab, sn,
                                        logging);
    }
    /** All rvectors are associated with the root paving, no spliting,
    no random number generator supplied, default will be created. */
    bool insertSampleFromRSSample(size_t samplesize, 
            const RSSample& rss, int lab, LOGGING_LEVEL logging)
    {
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rss, lab, sn,
                                        logging);
    }
	/** All rvectors are associated with the root paving, no spliting,
    no random number generator supplied, default will be created. */
    bool insertSampleFromRSSample(size_t samplesize,
            const RSSample& rss, LOGGING_LEVEL logging)
    {
        int lab = label; // label for this
        SplitNever sn; // a dummy split decision object
        return insertSampleFromRSSample(samplesize, rss, lab, sn,
                                        logging);
    }

    /** Adaptive splitting with each data point inserted,
    random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize,  
			gsl_rng * rgsl,
            const RSSample& rss, int lab,
			const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging);
    /** Adaptive splitting with each data point inserted,
    random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize, gsl_rng * rgsl,
            const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
	{
		int lab = label; // label for this
        return insertSampleFromRSSample(samplesize, rgsl, rss, lab,
						boolTest, logging);
	}
    /** Adaptive splitting with each data point inserted,
    seed for creating random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, int lab, 
			const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging);
    /** Adaptive splitting with each data point inserted,
    seed for creating random number generator supplied.*/
    bool insertSampleFromRSSample(size_t samplesize, int seed,
            const RSSample& rss, const SplitDecisionObj& boolTest,
            LOGGING_LEVEL logging)
	{
		int lab = label; // label for this
        return insertSampleFromRSSample(samplesize, seed, rss, lab,
						boolTest, logging);
	}
    /** Adaptive splitting with each data point inserted,
    no random number generator supplied, default will be created.*/
    bool insertSampleFromRSSample(size_t samplesize, 
								const RSSample& rss,
								int lab,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging);
	/** Adaptive splitting with each data point inserted,
    no random number generator supplied, default will be created.*/
    bool insertSampleFromRSSample(size_t samplesize, const RSSample& rss,
                                const SplitDecisionObj& boolTest,
                                LOGGING_LEVEL logging)
	{
		int lab = label; // label for this
        return insertSampleFromRSSample(samplesize, rss, lab,
						boolTest, logging);
	}
    //@}

	
    /** @name prioritySplit methods.

    These methods takes a histogram and progressively splits using a priority
    queue to determine which node to split first. Splitting continues until
    some criteria applying either to individual nodes or to the histogram
    as a whole is satisfied, or there are no more splittable nodes.

    Nodes are considered to be splittable if they satisfy two criteria:
    First, their volume is greater than the minimum volume specified for the
    histogram as a whole through \a minVolB.
    Second if both prospective children would have at least \a
    minChildPoints data points associated with them.

    If more than one node is equally 'large', on the basis of the node
    comparison compTest used, then a random choice is made between all equally
    large nodes to find the node which will be split.

    The random number generator used for random selection between equally
    'large' nodes uses a default seed to ensure that results can be
    replicated.  If you are looking at distributions of results across
    mulitple histograms, supply the random number generator to the priority
    queue to ensure that each histogram will make different random choices.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this does not satisfy
	the criteria implied by \a minVolB and \a minChildPoints, ie if
	a cherry not satisfying the implied min volume and child points
	rules has been split.
	
	Throws an std::logic_error if the split becomes muddled because of 
	some failure within the logic of the algorithm itself.
	 
	Aborts if there are no splittable leaves left (or none at the start).

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
    \param rgsl is a pointer to a gsl random number generator.
    \return true if \a he is satisfied at the end of the operation,
	false otherwise (will return false if the split could not be started
	or had to abort before \a he was satisfied because there were no
	more splittable nodes)
	\pre hasSubPaving() == true.
	\pre The state of this is legal with respect to \a minChildPoints and minVolB
	ie this does not include cherries that should not have been split
	according to these criteria.  	
    \post if the method returned true, this will be in a state such
	that \a he is satisfied;
	if the method returned false splitting could not start or had to 
	be aborted before \a he was satisfied.    */
    //@{
    /** %HistEvalObj to stop splitting. 
	 * Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints)
    { return prioritySplit(compTest, he, logging, minChildPoints, 0.0); }

    /** %HistEvalObj to stop splitting. 
	Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB)
    { return prioritySplit(compTest, he, logging, 0, minVolB); }

    /** %HistEvalObj to stop splitting. 
	Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging)
    { return prioritySplit(compTest, he, logging, 0, 0.0); }

    /** %HistEvalObj to stop splitting. 
	minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB);

    /** %HistEvalObj to stop splitting. 
	With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl)
    { return prioritySplit(compTest, he, logging, minChildPoints, 0.0, rgsl); }

    /** %HistEvalObj to stop splitting. 
	With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl)
    { return prioritySplit(compTest, he, logging, 0, minVolB, rgsl); }

    /** %HistEvalObj to stop splitting. 
	With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging, gsl_rng * rgsl)
    { return prioritySplit(compTest, he, logging, 0, 0.0, rgsl); }

    /** %HistEvalObj to stop splitting. 
	With random number generator. All other parameters supplied.*/
    bool prioritySplit(const NodeCompObj& compTest, const HistEvalObj& he,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl);


    //@}

	/** @name prioritySplit methods.

    These methods takes a histogram and progressively splits using a priority
    queue to determine which node to split first. Splitting continues until
    there are \a maxLeaves leaves, or there are no more splittable nodes.

	\note
	These overloadings of the prioritySplit function are supplied because 
	it is much more efficient to check for number of leaves directly than
	by using a function object.

    Nodes are considered to be splittable if they satisfy two criteria:
    First, their volume is greater than the minimum volume specified for the
    histogram as a whole through \a minVolB.
    Second if both prospective children would have at least \a
    minChildPoints data points associated with them.

    If more than one node is equally 'large', on the basis of the node
    comparison compTest used, then a random choice is made between all equally
    large nodes to find the node which will be split.

    The random number generator used for random selection between equally
    'large' nodes uses a default seed to ensure that results can be
    replicated.  If you are looking at distributions of results across
    mulitple histograms, supply the random number generator to the priority
    queue to ensure that each histogram will make different random choices.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this does not satisfy
	the criteria implied by \a minVolB and \a minChildPoints, ie if
	a cherry not satisfying the implied min volume and child points
	rules has been split.
	
	Throws an std::logic_error if the split becomes muddled because of 
	some failure within the logic of the algorithm itself.
	 
	Aborts if there are no splittable leaves left (or none at the start).

    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param maxLeaves is the number of leaves to have for the
	final histogram.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging).
    \param minChildPoints is the minimum number of points any prospective child
    must have for a leaf node to be splittable.
    \param minVolB is a multiplier applied to (log n)^2/n to give the the
    minimum volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n is not splittable.  Important with AIC or COPERR.
    \param rgsl is a pointer to a gsl random number generator.
    \return true if \a he is satisfied at the end of the operation,
	false otherwise (will return false if the split could not be started
	or had to abort before \a he was satisfied because there were no
	more splittable nodes)
	\pre hasSubPaving() == true.
	\pre The state of this is legal with respect to \a minChildPoints and minVolB
	ie this does not include cherries that should not have been split
	according to these criteria. 
	\post if the method returned true, this will have \a maxLeaves leaves;
	if the method returned false splitting could not start or had to 
	be aborted before this got \a maxLeaves leaves.    */
    //@{
    /** \a maxLeaves to stop splitting. 
	 * Only minChildPoints supplied, no minvolB, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints)
    { return prioritySplit(compTest, maxLeaves, logging, minChildPoints, 0.0); }

    /** \a maxLeaves to stop splitting. 
	Only minVolB supplied, no minChildPoints, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      double minVolB)
    { return prioritySplit(compTest, maxLeaves, logging, 0, minVolB); }

    /** \a maxLeaves to stop splitting. 
	Neither minVolB nor minChildPoints supplied, no random number generator. */
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging)
    { return prioritySplit(compTest, maxLeaves, logging, 0, 0.0); }

    /** \a maxLeaves to stop splitting. 
	minVolB and minChildPoints supplied but no random number generator.*/
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB);

    /** \a maxLeaves to stop splitting. 
	With random number generator. Only minChildPoints supplied, no minvolB.*/
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, gsl_rng * rgsl)
    { return prioritySplit(compTest, maxLeaves, logging, minChildPoints, 0.0, rgsl); }

    /** \a maxLeaves to stop splitting. 
	With random number generator. Only minVolB supplied, no minChildPoints. */
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      double minVolB, gsl_rng * rgsl)
    { return prioritySplit(compTest, maxLeaves, logging, 0, minVolB, rgsl); }

    /** \a maxLeaves to stop splitting. 
	With random number generator. Neither minVolB nor minChildPoints supplied. */
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging, gsl_rng * rgsl)
    { return prioritySplit(compTest, maxLeaves, logging, 0, 0.0, rgsl); }

    /** \a maxLeaves to stop splitting. 
	With random number generator. All other parameters supplied.*/
    bool prioritySplit(const NodeCompObj& compTest, size_t maxLeaves,
                      LOGGING_LEVEL logging,
                      size_t minChildPoints, double minVolB, gsl_rng * rgsl);

	
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

	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.

    Throws an UnfulfillableRequest_Error if the subpaving that this manages
	is a leaf.

	Throws an std::logic_error if the merge becomes muddled.

    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param he is an instance of a class which provides a function to determine
    when to stop merging.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)
    \return true if the priority merge was successful, false otherwise.
	\pre hasSubPaving() == true.    */
    bool priorityMerge(const NodeCompObj& compTest, const HistEvalObj& he,
                        LOGGING_LEVEL logging=NOLOG);
	
	/*! \brief Priority merge to reduce number of leaves in histogram.

    This method takes a histogram where where all the data is associated with
    multiple nodes and progressively merges the children of sub-terminal leaves
    using a priority queue to determine which node to merge first.
    Merging continues until histogram has only \a minLeaves leaves or
	merging has had to be aborted.

    If more than one node is equally 'small', on the basis of the node
    comparison compTest used, then a random choice is made between all equally
    small nodes to find the node which will be merged.

	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.

    Throws an UnfulfillableRequest_Error if the subpaving that this manages
	is a leaf.

	Throws an std::logic_error if the merge becomes muddled.

    \param compTest is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param minLeaves is the number of leaves to aim for.
    \param logging an enum controlling whether histogram creation output is
    sent to a log file (defaults to no logging)
    \return true if the priority merge was successful, false otherwise.
	\pre hasSubPaving() == true.    */
    bool priorityMerge(const NodeCompObj& compTest, size_t minLeaves,
                        LOGGING_LEVEL logging=NOLOG);

    /*! \brief Merge a multileaf histogram up to just root.
    
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.

	Throws an std::logic_error if the merge becomes muddled.

    No prioritisation, just brute force
    \return true if the histogram was merged up to a single leaf.
	\pre hasSubPaving() == true.   */
    bool mergeUp();

	/*! \brief Get summary information on non-empty
	box numbers and volumes in this histogram.

	\return A pair, where the first value in the pair is the
	total number of non-empty leaf boxes in the subpaving managed 
	by this, and 
	the second value is the proportion of the total volume
	of the root box of this that is in the non-empty leaf boxes
	of the subpaving managed by this.
	\pre This must have a subpaving with a box to managed.*/		
	std::pair<size_t, cxsc::real> getNonEmptyBoxSummary() const;


	/*! \brief Find the coverage value for a data point.
	
	The coverage value is 
	1 - (sum of density of all the boxes with heights >  
	the height of the box where the data point is).
	
	If the point is not in the histogram at all, coverage = 0;
	If the point is in the lowest box in the histogram, 
	coverage = count lowest box / total count;
	If the point is in the highest box of the histogram, coverage = 1
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an IncompatibleDimensions_Error if the given point does not
	have the same dimensions as the subpaving that this manages.
	
    \param pt the point to find coverage for
	\return coverage for the point given.	
	\pre hasSubPaving() == true and the dimensions of the point
	and the root paving match.*/
	double findCoverage(const rvector& pt) const;

	/*! \brief Find the empirical density for a data point.
	
	The empirical density is the relative density of the histogram
	at the box the given data point is in.
	
	If the point is not in the histogram at all, empirical density = 0;
	If the point is in the some box in the histogram with height d
	(d = count in box / (total count x box volume), then d is the 
	empirical density.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an IncompatibleDimensions_Error if the given point does not
	have the same dimensions as the subpaving that this manages.

	\param pt the point to find empirical density for.
	\return the empirical density at the point.
	\pre hasSubPaving() == true and the dimensions of the point
	and the root paving match.		*/
	double findEmpiricalDensity(const rvector& pt) const;
	
	//new AHABC
	/*! \brief Get whether the root box for this contains a given point.
	 
	Does the support of the histogram include the point \a pt?
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an IncompatibleDimensions_Error if the given point does not
	have the same dimensions as the subpaving that this manages.

	\param pt the point to check.
	\return true if the root box of the subpaving this manages contains
	\a pt, false otherwise.
	\pre hasSubPaving() == true and the dimensions of the point
	and the root paving match.		*/
	bool histBoxContains(const rvector& pt) const;

    /*! \brief Split a histogram  to a specified shape.
    
   	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving box is empty.

	Prints a message to the standard error output if the instruction
	could not be carried out. 

    \param instruction specifies the required shape, eg "3, 3, 2, 1"
	\return true if the split was successful, false otherwise
    \pre hasSubPaving() == true.*/
    bool splitToShape(std::string instruction);

    
	/** @name Get a PiecewiseConstantFunction as the average
	of samples from a Markov-Chain Monte Carlo process starting from this.
	
	\note It is \b important to note that this is changed during the MCMC
	operation and that this operation
	may be implemented in such a way that, to improve efficiency,
	the state of this at the end of the operation is
	not necessarily equivalent to the final state in the chain.

    The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  The
    the Metropolis-Hastings algorithm is used on the SPSnode tree managed by
    this Adaptive Histogram to generate PiecewiseConstantFunction
	samples from the histogram state
    probability density.  This method returns the average over the samples

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

    \internal
	When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws a std::invalid_argument exception if \a burnin > \a loops.
	 
	Throws a std::logic_error exception if the histogram given as the
	starting point of the chain is "illegal" with respect to the given
	\a minPoints, ie has cherries which should not have been split
	according to this criterion.

	Throws a std::logic_error exception if the histogram becomes
	corrupted during the evolution of the chain.

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
    \return A %PiecewiseConstantFunction representing the average of
	the samples collected from the chain started from this.
    \pre hasSubPaving() == true, and \a loops >= \a burnin.
	\post This has the state at least as split as any state reached by
	the chain in the MCMC process.*/
	//@{
	/*! \brief MCMC average method using a random number generator
	set and seeded internally, using the gsl default generator and
	seed.*/
	PiecewiseConstantFunction MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging);
					
	/*! \brief MCMC average method using a random number generator
	set internally using the gsl default generator, and with seed
	set by the user.
	
	\param seed is the seed to use for the random number generator.*/
	PiecewiseConstantFunction MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					long unsigned int seed);
	
	/*! \brief MCMC average method using a random number generator
	supplied by the user.
	
	\param rgsl is the random number generator to use.*/
	PiecewiseConstantFunction MCMC(
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    MCMCProposal& proposal, LogMCMCPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					gsl_rng * rgsl);
	//@}
	
	
	/*! @name Generating MCMC samples from histogram state space.

    \note It is \b important to note that this is changed during the MCMC
	operation and that this operation
	may be implemented in such a way that, to improve efficiency,
	the state of this at the end of the operation is
	not necessarily equivalent to the final state in the chain.

	The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  In this implementation the
    the Metropolis-Hastings algorithm is used on the SPSnode tree managed by
    this Adaptive Histogram to generate PiecewiseConstantFunction
	samples from the histogram state
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

    \internal
	When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.

    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws a std::invalid_argument exception if \a burnin > \a loops.

	Throws a std::logic_error exception if the histogram given as the
	starting point of the chain is "illegal" with respect to the given
	\a minPoints, ie has cherries which should not have been split
	according to this criterion.
	 
	Throws a std::logic_error exception if the histogram becomes
	corrupted during the evolution of the chain.
	
    \param samples is a reference to a container to add 
    PiecewiseConstantFunction samples to.
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
    (the container is unchanged if MCMC process is not successful).
    \pre hasSubPaving() == true.
	\pre The state of this is legal with respect to \a minPoints, ie
	this does not include cherries that should not have been split
	according to this criterion.  
	\pre \a loops >= \a burnin, and 
	(\a loops-\a burnin+1)/ \a thinout < 1  when \a thinout > 0.
	\post if the MCMC process is successful \a samples contains the\
	 samples collected from the chain (\a samples is unchanged if 
	 the MCMC process is not successful).*/
	//@{
	/*! \brief MCMC samples method using a random number generator
	set and seeded internally, using the gsl default generator and
	seed.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging);
	
	/*! \brief MCMC samples method using a random number generator
	set internally using the gsl default generator, and with seed
	set by the user.
	
	\param seed is the seed to use for the random number generator.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						long unsigned int seed);

	/*! \brief MCMC samples method using a random number generator
	supplied by the user.
	
	\param rgsl is the random number generator to use.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamples(
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						MCMCProposal& proposal, LogMCMCPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						gsl_rng * rgsl);
	//@}
	
	//docs
	/*! @name Generating MCMC samples from histogram state space using 
	an Independent Metropolis Hastings chain.

	/note This is an \b experimental routine that has not been fully tested
	and where the implemententation possibly requires further work. It is 
	not recommended for use in anything other than an experimental setting.

    The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  In this implementation the
    the Independent Metropolis-Hastings algorithm is used on the SPSnode
	tree managed by
    this %AdaptiveHistogram to generate PiecewiseConstantFunction
	samples from the histogram state
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

    \internal
	When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.

    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws a std::invalid_argument exception if \a burnin > \a loops.

	Throws a std::logic_error exception if the histogram given as the
	starting point of the chain is "illegal" with respect to the given
	\a minPoints, ie has cherries which should not have been split
	according to this criterion.
	 
	Throws a std::logic_error exception if the histogram becomes
	corrupted during the evolution of the chain.
	
    \param samples is a reference to a container to add 
    PiecewiseConstantFunction samples to.
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
    (the container is unchanged if MCMC process is not successful).
    \pre hasSubPaving() == true.
	\pre The state of this is legal with respect to \a minPoints, ie
	this does not include cherries that should not have been split
	according to this criterion.  
	\pre \a loops >= \a burnin.
	\post This has the state at least as split as any state reached by
	the chain in the MCMC process.
	\post if the MCMC process is successful \a samples contains the\
	 samples collected from the chain (\a samples is unchanged if 
	 the MCMC process is not successful).*/
	//@{
	/*! \brief MCMC samples method using a random number generator
	set and seeded internally, using the gsl default generator and
	seed.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging);
	
	/*! \brief MCMC samples method using a random number generator
	set internally using the gsl default generator, and with seed
	set by the user.
	
	\param seed is the seed to use for the random number generator.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						long unsigned int seed);

	/*! \brief MCMC samples method using a random number generator
	supplied by the user.
	
	\param rgsl is the random number generator to use.*/
	std::vector < PiecewiseConstantFunction >& MCMCsamplesIMH(
						size_t maxLeaves,
						std::vector < PiecewiseConstantFunction >& samples, 
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, LOGGING_LEVEL logging,
						gsl_rng * rgsl);
	//@}
	

    /*! @name Change the state of this Adaptive Histogram using MCMC process
	 * and capture changes in an information object.

	\internal This method should really be private, not exposed.  I exposed
	it so that the auto-mcmc objects could use it.  Better design 
	might have either auto-mcmc objects as friends of %AdaptiveHistogram
	(but that's not very nice either) or create a adaptation of the 
	%AdaptiveHistogram that is a friend of %AdaptiveHistogram and which
	exposes an adaptation of this method (but not other private methods 
	and properties) and have the auto-mcmc methods use that adaptation.
	
	\internal To improve efficiency this method runs this 
	%AdaptiveHistogram and the RealMappedSPnode rooted at \a rmsp
	side-by-side.  The RealMappedSPnode rooted at \a rmsp represents 
	a virtual state of this that will eventually be reflected in samples,
	averages etc.  For efficiency this itself may be more split than
	the virtual state (ie merge changes proposed and accepted may be 
	implemented only in \a rmsp, not this).  The actual state of this and
	the RealMappedSPnode rooted at \a rmsp may differ both before and
	after the operation. Both before and after, the leaf nodes in 
	\a nodes should be equivalent to the splittable leaf nodes in
	\a leaves, and the cherry nodes in \a nodes should be equivalent
	to the cherry nodes in \a cherries, but the %SPS tree managed by 
	this may have more actual leaf nodes (both unsplittable and 
	splittable) than are in \a nodes and 
	more actual splittable leaf nodes than are in \a leaves, 
	and more actual cherry nodes than
	are in \a nodes (or \a cherries)
	
	\todo Improve design for this method, which is exposed so that 
	auto-mcmc types can use it.
	
	\note Changes in this %AdaptiveHistogram are 'mirrored' by changes in
	the RealMappedSPnode rooted at \a rmsp.  The RealMappedSPnode rooted
	at \a rmsp represents the virtual state of this.  This operation
	may be implemented in such a way that, to improve efficiency,
	the state of this both before and after the operation is actually
	different than that of the the RealMappedSPnode rooted at \a rmsp.

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
    effect on the cherries which can be proposed for merging. 
	
	minVol > 0.0 also restricts the leaves that can be selected for a change.
	If the volume of the leaf is < 2*minVol (so that the children would have
	vol < minVol) then the leaf is not splittable.
	
	A leaf with volume >= 2*minVol whose
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

    \pre This manages an SPSnode tree.  \a nodes is a container of
    pointers to some of the splittable leaf and cherry nodes of 
	this %SPSnode tree,
	ordered with the leaves first followed by the cherries.  The nodes
	in \nodes may represent a subset of the actual nodes of the %SPSnode
	tree managed by this.  
	\a rmsp is a %RealMappedSPnode which is as split or less split than
	the %SPSnode tree managed by this, all of whose leaf nodes are in 
	\a leaves and all of whose cherry nodes are in \a cherries.   
	
    \post the SPSnode tree in some state which is either in the same state as
    the pre-state or in a new state reachable by a single split
	from the state before the change.  The container \a nodes 
	is updated to reflect the effect on the splittable nodes in \a nodes
	that results from a split change proposed 
	and accepted in the course of the operation or the effect on the 
	cherries in \a nodes of 
	a merge change proposed and accepted in the course of the operation.
	\a numLeaves and \a numCherries are updated to reflect changes in
	the contents of \a nodes.  
	\a rmsp is updated to reflect the effect on \a rmsp
	of a split or merge change equivalent to an accepted split
	or merge change that takes place in the course of the operation, and
	\a leaves and \a cherries are updated to reflect the change in
	\a rmsp.  Results may be logged to \a sHist as directed by the 
	parameter \a logging.

    \internal The container of splittable leaf nodes and cherries is
    maintained separately from the overall state of the tree to save having to
    repeatedly assess whether a leaf can be split (given minPoints).  The
    number of splittable leaves and cherries is maintained for convenience to
    avoid repeatedly counting nodes of different types in the container.
	Similarly both \a rmsp and the containers of leaves and cherries in 
	\a rmsp are passed separately in order to avoid having to repeatedly
	re-collect the leaves and cherries.

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
	\param minVol is the minimum volume to check for to be able to split a node.
    \param rmsp is the root of a %RealMappedSPnode that represents
	the virtual state of this.
	\param leaves is a container of the leaf nodes of \a rmsp.
	\param cherries is a container of the cherry nodes of \a rmsp.
	\param info is a reference to an object to record some of the 
	effects on this of the change of state that takes place in the 
	course of the operation.
	\param rgsl is a pointer to a random number generator.
    \param logging an enum controlling whether histogram creation output is
    sent to log files (defaults to no logging).  TXT gives logging to a
    txt file, TXTANDGRAPH gives txt file logging and graphs.
    \param sHist is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain.
    \return Change in log-posterior as a result of the change in state.    */

	//@{
	/*! \brief Version of method using \a minVol.*/
	AdaptiveHistogram::ChangeOfStateInformation& changeMCMCState(
						SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						AdaptiveHistogram::ChangeOfStateInformation& info,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i);
						
	AdaptiveHistogram::ChangeOfStateInformation& changeMCMCState(
						SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						real minVol,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						AdaptiveHistogram::ChangeOfStateInformation& info,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i);
	
	bool changeMCMCStateIMH(
						size_t maxLeaves,
                        const MCMCPartitionGenerator& partitioner,
						const LogMCMCIMHPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode& rmsp,
						size_t& numLeaves,
						real& unscaledLogLik,
						real& logPosterior,
						gsl_rng* rgsl, 
						LOGGING_LEVEL logging,
                        const std::string& sDec,
						int i,
						const string& failureLogFilename);
	
	//docs
	/** @name Get a PiecewiseConstantFunction as the average
	of samples from a Markov-Chain Monte Carlo process starting from this.

    The leaves of the SPSnode tree represent the partition of the data space
    (the root box of the tree).  A histogram state is a particular partition
    of the root box which will be represented by a particular tree number
    and disposition of nodes of the tree.

    The Markov-Chain Monte Carlo process considers possible histogram states,
    given data, as a probability distribution.  The
    the Metropolis-Hastings algorithm is used on the SPSnode tree managed by
    this Adaptive Histogram to generate PiecewiseConstantFunction
	samples from the histogram state
    probability density.  This method returns the average over the samples

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

    \internal
	When minPoints > 0, proposals are effectively drawn from set of leaf and
    cherry nodes which does not include any leaf which, if split, would have
    a child whose number of points is < minPoints, unless that child leaf node
    has 0 points and its sibling has all the parent's points and number of
    parent's points >= minPoints.  Thus the implementation
    needs to distinguish between the overall state of the tree and the
    set of <b>splittable leaf nodes</b>.

    The method can log the process, including the components of the calculation
    for each change in state and .dot graphs for each state in the chain.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws a std::invalid_argument exception if \a burnin > \a loops.
	 
	Throws a std::logic_error exception if the histogram given as the
	starting point of the chain is "illegal" with respect to the given
	\a minPoints, ie has cherries which should not have been split
	according to this criterion.

	Throws a std::logic_error exception if the histogram becomes
	corrupted during the evolution of the chain.

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
    \return A %PiecewiseConstantFunction representing the average of
	the samples collected from the chain started from this.
    \pre hasSubPaving() == true, and \a loops >= \a burnin.
	\post This has the state at least as split as any state reached by
	the chain in the MCMC process.*/
	//@{
	/*! \brief MCMC average method using a random number generator
	set and seeded internally, using the gsl default generator and
	seed.*/
	PiecewiseConstantFunction MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging);
					
	/*! \brief MCMC average method using a random number generator
	set internally using the gsl default generator, and with seed
	set by the user.
	
	\param seed is the seed to use for the random number generator.*/
	PiecewiseConstantFunction MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					long unsigned int seed);
	
	/*! \brief MCMC average method using a random number generator
	supplied by the user.
	
	\param rgsl is the random number generator to use.*/
	PiecewiseConstantFunction MCMC_IMH(
					size_t maxLeaves,
					unsigned int loops, unsigned int burnin,
                    unsigned int thinout,
                    LogMCMCIMHPrior& logPrior,
                    size_t minPoints, LOGGING_LEVEL logging,
					gsl_rng * rgsl);
	//@}

	
	/*! @name Find maximum log-posteriors using two priorty split functions.
	 * 
	 Do a number of priority queue splits with two priority queue functions,
	 operating one after the other, and find the maximum log-posterior points 
	 reached by each queue overall.
	 
	 The PrioritySplitQueueEvaluator \a psqe controls the first part of
	 each queue and the %PrioritySplitQueueEvaluator \a psqePosterior
	 controls the second part of each queue.  Each 
	 %PrioritySplitQueueEvaluator will specify a maximum number of
	 leaves.  One queue using just the
	 \a psqePosterior controlled queue is always run.  
	 
	 The parameter 
	 \a checkMaxStep effectively controls the number of different 
	 versions of the \a psqe controlled queue that are run:  each one is
	 run to cover a number of interations (splits) that increments by 
	 \a checkMaxStep each time, until the overall maximum number of 
	 leaves specified in \a psqe is reached or exceeded.  Each \a psqe
	 controlled queue thus started is monitored and the 
	 maximum \a support-adjusted log-posterior point in the 
	 last checkMaxStep splis is 
	 recorded and used as the start point for a \a psqePosterior
	 controlled queue that runs until the maximum number of leaves
	 specified in \a pqsePosterior is reached (or until there
	 are no more splittable nodes left).  If \a checkMaxStep = 0, no
	 \a psqe controlled queues are run (just the first queue using
	 just \a psqePosterior).  
	 
	 A support-adjusted log-posterior is unnormalised and is
	 calculated by adjusting the unnormalised log-posterior for the 
	 total empty box volume over the whole of the root box (ie, tries
	 to use the actual support as represented by boxes with at 
	 least some points in them).
	   
	 The effect is thus to run the \a psqe controlled queue in blocks
	 and within each block, find the maximum support-adjusted log-posterior
	 point, and use that as the start point for a \pqsePosterior
	 controlled queue, and record the (un-normalised) log-posterior for
	 the combined two-phase queue for the block. This approach is used 
	 as a kind of trial-and-error method of finding an appropriate 
	 number of leaves to run the first phase (typically the carving
	 phase) for when we have very little information to go on.   
	  
	 For each overall queue, including the
	 one just controlled by \a pqsePosterior, the maximum log-posterior
	 point is found and recorded using \a carvedLaunchPoints (the 
	 number of leaves the first \a psqe controlled queue runs to - for
	 the queue that only uses the \a psqePosterior queue, this is 
	 taken as the number of leaves in this at the start of the 
	 operation), 
	 \a maxPosteriorPoints (the number of leaves that the \a psqePosterior
	 controlled queue runs to) and \a maxPosteriors (the (unnormalised)
	 value of the maximum log-posterior found)). 
	
	\a emptyBoxVolsVec and \a posteriorSupportVec collect information 
	from the total or overall \a pqse phase, ie the whole \a psqe 
	controlled queue run
	until the maximum number of leaves specified in \a psqe is reached.
	\a emptyBoxVolsVec is used to store the total volume of the empty 
	boxes in this for each loop (split) in the \a psqe phase, and 
	\a posteriorSupportVec is used to store the (un-normalised)
	support-adjusted log-posterior of this after each split in the
	\a psqe phase.
	
	\a stopOnMaxPosterior controls how long the \a psqePosterior queue
	is run for after
	an (apparent) maximum posterior point is found.  This
	is only relevant to the queues that include a \a psqe controlled
	phase (the first, \a pqsePosterior only controlled, queue
	should be run until the maximum leaves specified in \a pqsePosterior
	is reached or until there are no more splittable leaves).  If 
	\a stopOnMaxPosterior is false then each \pqsePosterior controlled 
	phase of a queue with an initial phase controlled by \a pqse is also
	run until the maximum leaves specified in \a pqsePosterior
	is reached (or until there are no more splittable leaves). If 
	\a stopOnMaxPosterior is true then each \pqsePosterior controlled 
	queue with an initial phase controlled by \a pqse is still
	run for some distance past a point that appears to
	represent a maximum (to check that this is not just a local
	maximum), but not necessarily until the maximum leaves specified 
	in \a pqsePosterior is reached.
	
	In each queue, if the log-posteriors are still going up at the 
	end of the queue, ie when the maximum leaves specified in 
	\a pqsePosterior is reached, then the final log-posterior is
	recorded.  
	 
	\internal \a shiftCatalan is left over from Gloria's implementation
	- I don't think that it is doing what she actually wanted.  I am 
	using it to shift the prior to the 1- leaf point (log prior = 0) 
	even if the actual number of leaves in this at the start is  1.
	
	\param psqe an object controlling the first phase of each queue 
	(typically, the 'carving' phase).
	\param psqePosterior an object controlling the second phase of 
	each queue (typically, the 'seb' phase).
	\param logging an enum to control logging of the process.
	\param minChildPoints the minimum points that each child of a node
	that has any points at all must have for that node to be considered
	to be splittable.   
	\param minVol the minimum volume that any leaf node should have. 
	\param logPrior the prior object to use for calculating 
	log-posterior values.
	\param emptyBoxVolsVec a container in which to store the total
	empty box volume resulting at each point (iteration, split)
	in the \a psqe queue, including that at the start of
	the operation. 
	\param posteriorSupportVec a container in which to store the 
	support-adjusted log-posterior at each point (iteration, split)
	in the \a psqe queue, including that at the start of
	the operation. The support-adjusted log-posterior is the 
	unnormalised log-posterior adjusted for the empty box volume.  
	\param checkMaxStep is the size of the blocks (number of 
	iterations, splits in each block) in the \a psqe controlled
	queue.
	\param stopOnMaxPosterior is an indicator to control how long the
	\a pqsePosterior queue is run for after an apparent maximum is 
	found for every queue started in a block of the \a psqe controlled
	queue.  If this is false, each \a pqsePosterior queue is run until
	the maximum number of leaves specified in \a pqsePosterior is reached
	(or until there are no more splittlable leaves).   
	\param carvedLaunchPoints a container in which to record the number
	of leaves in the support-adjusted maximum log-posterior point for
	each block in the \a psqe controlled queue (the first entry will
	always be number of leaves in the initial state of this at the 
	start of the operation).
	\param maxPosteriorPoints a container in which to record the number
	of leaves in the maximum log-posterior point for each combined queue
	run.  The the first entry will be the number of leaves in the 
	maximum posterior point of the queue using only \a pqsePosterior, 
	subsequent entries (present if checkMaxStep > 0) record the 
	number of leaves in the max log-posterior point reached by the
	\pqsePosterior controlled queue started from the maximum 
	support-adjusted log-posterior point in each block.
	\param maxPosteriors a container in which to record the value of 
	the unnormalised log-posterior at the maximum log-posterior point 
	for each combined queue
	run.  The the first entry will be the the unnormalised log-posterior
	in the queue using only \a pqsePosterior, 
	subsequent entries (present if checkMaxStep > 0) record the 
	maximum unnormalised log-posterior reached by the
	\pqsePosterior controlled queue started from the maximum 
	support-adjusted log-posterior point in each block.
	\param rgsl a random number generator to use in the operation
	\param shiftCatalan a boolean to control whether the prior
	is shifted back.  If true, then the prior is shifted so that 
	the prior for the initial state of this is the prior of a state
	with 1 leaf and the prior for a state with n-(m-1) leaves is the
	prior otherwise used for a state with n leaves, where m is the 
	number of leaves in this at the start of the operation.
	\return true if the operation completed without running out of
	splittable nodes in the \a psqe controlled queue, false otherwise.  
	*/
	//@{
	
	/*! \brief Version taking \a minVol parameter. */
	bool prioritySplitWithSupportPosteriorMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
								LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								int checkMaxStep, // if 0 no carving
								bool stopOnMaxPosterior,
								std::vector<size_t>& carvedLaunchPoints,
								std::vector<size_t>& maxPosteriorPoints,
								std::vector<real>& maxPosteriors,
                                gsl_rng * rgsl,
                                bool shiftCatalan);
	
	/*! \brief Version without \a minVol parameter */	
	bool prioritySplitWithSupportPosteriorMaxLik(
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
                                bool shiftCatalan);
	//@}
	
														
	/*! @name Do a priority queue split using two priority queue 
	 functions, operating one after the other. 
	 * 
	 Typically used to get this into
	 the state previously found using 
	 prioritySplitWithSupportPosteriorMaxLik or some operation that
	 manipulates the results of prioritySplitWithSupportPosteriorMaxLik
	 to find other states of interest. 
	 
	 The PrioritySplitQueueEvaluator \a psqe controls the first part of
	 each queue and the %PrioritySplitQueueEvaluator \a psqePosterior
	 controls the second part of each queue.  Each 
	 %PrioritySplitQueueEvaluator will specify a maximum number of
	 leaves.  
	 
	\param psqe an object controlling the first phase of each queue 
	(typically, the 'carving' phase).
	\param psqePosterior an object controlling the second phase of 
	each queue (typically, the 'seb' phase).
	\param logging an enum to control logging of the process.
	\param minChildPoints the minimum points that each child of a node
	that has any points at all must have for that node to be considered
	to be splittable.   
	\param minVol the minimum volume that any leaf node should have. 
	\param logPrior the prior object to use for calculating 
	log-posterior values.
	\param posteriorVec a container in which to record the 
	unnormalised log-posterior in each state of the queue.  The 
	first entry is the unnormalised log-posterior at the start 
	of the operation.
	\param loglikVec a container in which to record the 
	log-likelihood in each state of the queue.  The 
	first entry is the log-likelihood at the start 
	of the operation.
	\param rgsl a random number generator to use in the operation
	\param shiftCatalan a boolean to control whether the prior
	is shifted back.  If true, then the prior is shifted so that 
	the prior for the initial state of this is the prior of a state
	with 1 leaf and the prior for a state with n-(m-1) leaves is the
	prior otherwise used for a state with n leaves, where m is the 
	number of leaves in this at the start of the operation.
	\return true if the operation completed without running out of
	splittable nodes, false otherwise.  
	*/
	//@{
	
	/*! \brief Version taking \a minVol parameter. */
	bool prioritySplitMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
					            LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& posteriorVec,
								std::vector<real>& loglikVec,
                                gsl_rng * rgsl,
                                bool shiftCatalan);
	
	/*! \brief Version without \a minVol parameter. */
	bool prioritySplitMaxLik(
								const PrioritySplitQueueEvaluator& psqe,
								const PrioritySplitQueueEvaluator& psqePosterior,
					            LOGGING_LEVEL logging,
                                size_t minChildPoints,
								LogMCMCPrior& logPrior,
                                std::vector<real>& posteriorVec,
								std::vector<real>& loglikVec,
                                gsl_rng * rgsl,
                                bool shiftCatalan);
	//@}
								
	/*! @name Run a single priority queue recording the 
	 * empty volumes, support-adjusted log-posterior, and log-posterior .
	 * 
	 The PrioritySplitQueueEvaluator \a psqe will specify a maximum 
	 number of leaves to run the queue to. Typically \a pqse will
	 contain some sort of support-carving queue function.  
	 	 
	 A support-adjusted log-posterior is unnormalised and is
	 calculated by adjusting the unnormalised log-posterior for the 
	 total empty box volume over the whole of the root box (ie, tries
	 to use the actual support as represented by boxes with at 
	 least some points in them).
	   
	 \a emptyBoxVolsVec, \a posteriorSupportVec, and \a posteriorVec
	 collect information from the queue.
	
	\param psqe an object controlling the queue 
	(typically, using a 'carving' type priority function ).
	\param logging an enum to control logging of the process.
	\param minChildPoints the minimum points that each child of a node
	that has any points at all must have for that node to be considered
	to be splittable.   
	\param minVol the minimum volume that any leaf node should have. 
	\param logPrior the prior object to use for calculating 
	log-posterior values.
	\param emptyBoxVolsVec a container in which to store the total
	empty box volume resulting at each point (iteration, split)
	in the queue, including that at the start of
	the operation. 
	\param posteriorSupportVec a container in which to store the 
	support-adjusted log-posterior at each point (iteration, split)
	in the queue, including that at the start of
	the operation. The support-adjusted log-posterior is the 
	unnormalised log-posterior adjusted for the empty box volume.  
	\param posteriorVec a container in which to store the 
	log-posterior at each point (iteration, split)
	in the queue, including that at the start of
	the operation.
	\param rgsl a random number generator to use in the operation
	\param shiftCatalan a boolean to control whether the prior
	is shifted back.  If true, then the prior is shifted so that 
	the prior for the initial state of this is the prior of a state
	with 1 leaf and the prior for a state with n-(m-1) leaves is the
	prior otherwise used for a state with n leaves, where m is the 
	number of leaves in this at the start of the operation.
	\return true if the operation completed without running out of
	splittable nodes in the queue, false otherwise.  
	*/
	//@{
	
	/*! \brief Version taking \a minVol parameter. */
	bool prioritySplitWithSupportPosterior(
								const PrioritySplitQueueEvaluator& psqe,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints,
								double minVol, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								std::vector<real>& posteriorVec,
								gsl_rng * rgsl,
                                bool shiftCatalan);
	
	/*! \brief Version without \a minVol parameter. */
	bool prioritySplitWithSupportPosterior(
								const PrioritySplitQueueEvaluator& psqe,
                                LOGGING_LEVEL logging,
                                size_t minChildPoints, 
                                LogMCMCPrior& logPrior,
                                std::vector<real>& emptyBoxVolsVec,
                                std::vector<real>& posteriorSupportVec,
								std::vector<real>& posteriorVec,
								gsl_rng * rgsl,
                                bool shiftCatalan);
	//@}

	//new
	/*! @name Run a single priority queue recording the 
	 * log-posterior and log-likelihood.
	 * 
	 The PrioritySplitQueueEvaluator \a psqe will specify a maximum 
	 number of leaves to run the queue to. Typically \a pqse will
	 contain some sort of SEB queue function.  
	 	 
	    
	 \a posteriorVec and \a loglikVec
	 collect information from the queue.
	
	\param psqe an object controlling the queue.
	\param logging an enum to control logging of the process.
	\param minChildPoints the minimum points that each child of a node
	that has any points at all must have for that node to be considered
	to be splittable.   
	\param minVol the minimum volume that any leaf node should have. 
	\param logPrior the prior object to use for calculating 
	log-posterior values.
	\param posteriorVec a container in which to store the 
	log-posterior at each point (iteration, split)
	in the queue, including that at the start of
	the operation.
	\param loglikVec a container in which to store the 
	log-likelihood at each point (iteration, split)
	in the queue, including that at the start of
	the operation.
	\param rgsl a random number generator to use in the operation
	\param shiftCatalan a boolean to control whether the prior
	is shifted back.  If true, then the prior is shifted so that 
	the prior for the initial state of this is the prior of a state
	with 1 leaf and the prior for a state with n-(m-1) leaves is the
	prior otherwise used for a state with n leaves, where m is the 
	number of leaves in this at the start of the operation.
	\return true if the operation completed without running out of
	splittable nodes in the queue, false otherwise.  
	*/
	//@{
	
	/*! \brief Version taking \a minVol parameter. */
	bool prioritySplitWithPosterior(
						const PrioritySplitQueueEvaluator& psqe,
						LOGGING_LEVEL logging,
						size_t minChildPoints, 
						LogMCMCPrior& logPrior,
						std::vector<real>& posteriorVec,
						std::vector<real>& loglikVec, 
						gsl_rng * rgsl,
						bool shiftCatalan);
	
	/*! \brief Version without \a minVol parameter. */
	bool prioritySplitWithPosterior(
						const PrioritySplitQueueEvaluator& psqe,
						LOGGING_LEVEL logging,
						size_t minChildPoints,
						double minVol, 
						LogMCMCPrior& logPrior,
						std::vector<real>& posteriorVec, 
						std::vector<real>& loglikVec,
						gsl_rng * rgsl,
						bool shiftCatalan);
	//@}

	/*! \brief Get the log-likelihood for this histogram.
	 * 
	 * The log-likelihood is the 
	 * sum over leaves of ln(leaf count/(leaf vol * total pts in histogram))
	 * 
	 * What to do if the histogram has no points in?  Could return
	 * Nan, or 0.   I have decided to return 0.  
	 * 
	 * \return The log-likelihood for this histogram.  Returns 0 if 
	 * this contains no points.  
	 * \pre This must have a subpaving to manage.	 */
	cxsc::real getLogLikelihood() const;
	
	
    /*! \brief Make a .dot graph file from histogram structure.

    Makes a simple .dot graph from the histogram using node names and 
	also makes the .png image for this graph.

    \post a .dot file and a .png in the same directory as the program creating
    it was run in.    */
    bool outputGraphDot() const;

	/*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then the description of the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param os is a reference to the stream to output the histogramm to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \return a reference to the given stream.
    */
    std::ostream & outputToStreamTabs(std::ostream & os, int prec = 5) const;

    /*! \brief Output the subpaving managed by this to a txt file.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then the description of the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param s the name of the txt file to send output to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
	\param confirm is a boolean controlling whether confirmation goes to
    console output.     */
	//@{
	void outputToTxtTabs(const std::string& s, int prec = 5) const;

	void outputToTxtTabs(const std::string& s, int prec, bool confirm) const;

	//@}

    /*! \brief Output the subpaving managed by this to a txt file.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then node contribution to EMP
    under COPERR, the change that would result in the EMP under COPERR if the
    node were split, the node contribution to EMP under AIC, the change that
    would result in the EMP under AIC if the node were split, and then the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param s the name of the txt file to send output to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \param confirm is a boolean controlling whether confirmation goes to
    console output.
    */
	//@{
    void outputToTxtTabsWithEMPs(const std::string& s,
                                int prec = 5) const;
								
	void outputToTxtTabsWithEMPs(const std::string& s,
                                int prec, bool confirm) const;

	//@}

    /*! \brief Output details of full sample (from root) to txt tile.

    Format is a mixture of alpha and  numeric data.

    \param s the name of the txt file to send output to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \param confirm is a boolean controlling whether confirmation goes to
    console output.
    */
	//@{
    void outputRootToTxt(const std::string& s, int prec = 5) const;

    void outputRootToTxt(const std::string& s, 
									int prec, bool confirm) const;


	//@}
	
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

	
	/*! Gets the L1 distance between this and another subpaving.

	The L1 distance is defined as the sum of the absolute values
	of the differences in 'area' represented by the leaf nodes
	of this and the other paving.  The 'area' represented
	by a leaf node is the proportion of the total count of the
	tree the leaf is part of that is in that leaf node.  
	
	Throws the following exceptions:
	<ul>
	<li>Throws a NullSubpavingPointer_Error if the subpaving 
	that this manages is a NULL pointer.</li>
	<li>Throws a NullSubpavingPointer_Error if the subpaving 
	that \a other  manages is a NULL pointer.</li>
	<li>Throws a NoBox_Error if either this's subpaving or 
	\a other's subpaving have no box.</li>
	<li>Throws a IncompatibleDimensions_Error if the dimensions
	and sizes of the boxes of the subpavings managed by this
	and \a other are not the same. </li>
	</ul>

	\param other the histogram to calculate the distance against
	the L1 distance against.
	\pre hasSubPaving() == true and \a other->hasSubPaving() == true
	and both have boxes of the same size and dimensions.*/
	cxsc::real getL1Distance(const AdaptiveHistogram& other) const;

    /*! Set holdAllStats
	 * 
	 * Changes the holdAllStats indicator to the required value
	 * and also ensures that the paving is changed appropriately, ie
	 * that paving properties and statistics held are altered.
	 * 
	 * Throws a NullSubpavingPointer_Error if the subpaving 
	that this manages is a NULL pointer.
	 * 
	 * \param setTo The value to set the member holdAllStats to;
	 * \pre hasSubPaving() == true.
	 * \post if setTo is true, paving will hold all stats for the data
	 * if setTo is false, paving will not hold all stats for the data.
	 * */
	void setHoldAllStats(bool setTo);
	
	
	/*! Set the label for this.
	 *\param lab the label for this.*/
	void setLabel(int lab);
	
	/*! Clear all the data in this.
	 * 
	 * Clears all the data from this's dataCollection and also 
	 * from the subpaving managed.
	 * 
	 * Throws a NullSubpavingPointer_Error if the subpaving 
	that this manages is a NULL pointer.

	 * \pre hasSubPaving() == true.
	 * \post dataCollection is empty and the subpaving has no data 
	 * (getRootCounter() == 0) 
	 * */
	void clearAllHistData();

	/*! \brief Swap the contents of this and another histogram.
	 
	\post After the swap \a adh will manage the subpaving that this used
	to manage, and this will manage the subpaving that \a adh used to 
	managed, and the values of the other data members of this
	and \a adh will also be swapped.  */
	void swap(AdaptiveHistogram& adh); // throw()

	/*! \brief Get a string summary of this histogram's properties.
		
	A string description of this. Includes the address of the subpaving
	managed but not details of that subpaving.
	
	\return the string summary.		*/
	std::string stringSummary() const;
	
	/*! \brief Get a string summary of the properties of this 
	which can be used for checking and debugging.
	
	Includes details of the subpaving that this manages including 
	addresses of the values associated with leaf nodes of the paving.
	
	\return the string summary.		*/
	std::string doubleCheckStringSummary() const;


	private:
	
		class NodePtrMeasurePair {
				
			public:
				NodePtrMeasurePair(SPSnode *p, real m);
			
				/* Comparison using measure.*/
				bool operator<(const NodePtrMeasurePair& rhs) const;
				
				std::string toString() const;
				
				SPSnode * nodePtr;
				cxsc::real measure;
				
			private:
				NodePtrMeasurePair();
			
		};
		
		typedef
		std::multiset<AdaptiveHistogram::NodePtrMeasurePair>
		PriorityQueueT;
		typedef PriorityQueueT::iterator PriorityQueueItrT;
		typedef PriorityQueueT::const_iterator PriorityQueueConstItrT;
	
		
		/*! \brief Inner type for evaluating priority split queues. */ 
		class PrivatePrioritySplitQueueEvaluator {
			
			public:
				PrivatePrioritySplitQueueEvaluator(const 
					PrioritySplitQueueEvaluator& psqe);
							
				cxsc::real measure(const SPSnode * const spn) const;
				
				bool stopQueueQuery(const PriorityQueueT& pq,
										size_t numLeaves) const;
				
				bool stopQueueQuery(size_t numLeaves) const;
				
				SPSNodeMeasure& measurer;
				real critStop;
				size_t maxLeaves;
				bool usingCritStop;
			
		};
	

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
    \post !hasSubPaving() when the function was entered, then
    rootPaving is pointed to a new SPSnode object whose root node has a box
    tailored to contain all the data read in.
    \post The data in theData has been put into the AdaptiveHistogram's
    dataCollection and also associated with  rootPaving's leaves via
    iterators to dataCollection.
    \return true if at least some data was successfully 
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
    \pre If hasSubPaving() == true then the dimensions of the data and
	the root box of the subpaving match.
    \post The data in theData which is within the boundaries of the rootBox
    has been put into the AdaptiveHistogram's dataCollection and also
    associated with  rootPaving's leaves via iterators to dataCollection.
    \return number of datapoints for which insertion is successful, 
	ie the number associated with the leaves of the actual subpaving.
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
    probability of acceptance of a proposed state under MCMC.

    \param s the name of the txt file to send output to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param deltaL is the log change in likelihood.
    \param deltaP is the log change in prior probability.
    \param deltaQ is the log change in transition probability.
    \param deltaPi is the log change in posterior probability.
    \param randChange is the draw from the Uniform (0,1) distribution which is
    compared to the changes in combined posterior and transition probabilities
    to determine whether the proposed change takes place.    */
    static void logMCMCDeltas(const std::string& s, int i,
                        real deltaL, real deltaP, real deltaQ, real deltaPi,
                        double randChange);

	/*! \brief Send a collection of changes in MCMC probabilities to log file.

    The probabilities being logged are the changes in the components of the
    probability of acceptance of a proposed state under MCMC.
	 
	This version logs more information about the target node.

    \param s the name of the txt file to send output to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param nodeType an indictator for the type of node (eg leaf, cherry).
	\param accepted an indicator for whether the move was accepted.
	\param nCount number of points in the node concerned.    
	\param deltaL is the log change in likelihood.
    \param deltaP is the log change in prior probability.
    \param deltaQ is the log change in transition probability.
    \param deltaPi is the log change in posterior.
    \param randChange is the draw from the Uniform (0,1) distribution which is
    compared to the changes in combined posterior and transition probabilities
    to determine whether the proposed change takes place.    */
    static void logMCMCDeltasAugmented(std::string s, int i,
                            int nodeType, int accepted, size_t nCount, 
							real deltaL, real deltaP, real deltaQ, real deltaPi,
                            double randChange);

	/*! \brief Send a collection of values for IMH MCMC to log file.

    The values being logged are the components of the
    probability of acceptance of a proposed state under IMH MCMC.
	 
	\param s the name of the txt file to send output to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param currentLeaves the current number of leaves.
	\param proposedLeaves the proposed number of leaves.
	\param deltaL is the log change in likelihood.
    \param deltaP is the log change in prior probability.
    \param deltaQ is the log change in transition probability.
    \param logAcceptance is the log acceptance probability.
    \param randChange is the draw from the Uniform (0,1) distribution which is
    compared to the changes in combined posterior and transition probabilities
    to determine whether the proposed change takes place.    
	\param accepted an indicator for whether the move is accepted. */
	static void logMCMCIMH(std::string s, int i,
                            size_t currentLeaves,
							size_t proposedLeaves,
							real deltaL,
							real deltaP,
							real deltaQ,
							real logAcceptance,
                            double randChange,
							bool accepted);

    /*! \brief Put header in a log file for MCMC.

    \param s is the name of the file to log the final state to.
    \param i is a number indicating the state number in a chain or series
    of changes in state.
    \param proposal a reference to a proposal function object.
    \param logPrior a reference to a log prior function object.
    */
    void MCMCStartLogFile(const std::string& s, int i, 
                                    const MCMCProposal& proposal,
                                    const LogMCMCPrior& logPrior);

	/*! \brief Put header in a log file for MCMC.

    \param s is the name of the file to log the final state to.
    \param logPrior a reference to a log prior function object.
    */
    void MCMCStartLogFile(const std::string& s, 
                                    const LogMCMCIMHPrior& logPrior);

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
    void MCMCLogFinalState(const std::string& s, int i);

	/*! Internal method used for independent Metropolis Hastings MCMC sampling*/
	std::vector < PiecewiseConstantFunction >& _MCMCsamplesIMH(
						std::vector < PiecewiseConstantFunction >& samples, 
						bool average,
						size_t maxLeaves,
						unsigned int loops, 
						unsigned int burnin,
						unsigned int thinout,
						LogMCMCIMHPrior& logPrior,
						size_t minPoints, 
						LOGGING_LEVEL logging,
						gsl_rng * rgsl);


	/*! Internal method used for MCMC sampling*/
	std::vector < PiecewiseConstantFunction >& _MCMCsamples(
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
						gsl_rng * rgsl);

    /*! @name Change the state of this Adaptive Histogram using MCMC process.

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
	\param minVol is the minimum volume to check for to be able to split a node.
    \param rmsp is the root of a %RealMappedSPnode that represents
	the virtual state of this.
	\param leaves is a container of the leaf nodes of \a rmsp.
	\param cherries is a container of the cherry nodes of \a rmsp.
	\param rgsl is a pointer to a random number generator.
    \param logging an enum controlling whether histogram creation output is
    sent to log files (defaults to no logging).  TXT gives logging to a
    txt file, TXTANDGRAPH gives txt file logging and graphs.
    \param sHist is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain.
    \return Change in log-posterior as a result of the change in state.    */

	//@{
	/*! \brief Version of method using \a minVol.*/
	cxsc::real changeMCMCState(SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						real minVol,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i);
	
	/*! \brief Version of method \b not using \a minVol.*/
	cxsc::real changeMCMCState(SPSnodeList& nodes,
                        size_t& numLeaves, size_t& numCherries,
                        MCMCProposal& proposal, LogMCMCPrior& logPrior,
                        size_t minPoints,
						RealMappedSPnode::ListPtrs& leaves,
						RealMappedSPnode::ListPtrs& cherries,
						gsl_rng* rgsl, LOGGING_LEVEL logging,
                        const std::string& sHist, int i);
	//@}

	
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

    \param target is a pointer to the target node proposed for splitting.
    \param proposal is a reference to a proposal distribution object.
    \param logPrior is a reference to a log prior object.
    \param rgsl is a pointer to a uniform random number generator.
    \param numLeaves is the number of splittable leaves in the SPSnode tree.
    \param numCherries is the number of cherries inthe SPSnode tree.
	\param nodes is a collection of splittable leaf and cherry nodes.
	\param realNumLeaves is the number actual leaves in the SPSnode tree
		(splittable or not).
    \param minPoints is the minimum number of points allowed in a box.
    \param volChecking is an indicator for whether node volume is 
	being checked to determine whether it is splittable.
    \param minVol is the minimum volume to check for if volChecking is true.
    \param addLeftChild reference to an indicator for whether split will result
	in the addition of the newLeftChild to the proposable nodes. The value
	of this indicator may be modified during the operation.
	\param addRightChild reference to an indicator for whether split will result
	in the addition of new right child to the proposable nodes. The value
	of this indicator may be modified during the operation.
	\param deductParent reference to an indicator for whether split will result
	in the loss of the parent of the target from the cherries. The value
	of this indicator may be modified during the operation.
	\param leftChildCount reference to a value for the number of points
	that will go to the left child if the node is split. The value
	of this parameter will be updated during the operation.
	\param rightChildCount reference to a value for the number of points
	that will go to the right child if the node is split. The value
	of this parameter will be updated during the operation.\
	\param deltaPi reference to a value for the change in the 
	log-posterior (log likelihood + log prior)
	that would result from the split. The value
	of this parameter will be updated during the operation.
	\param logging an enum controlling whether decision making output is
    sent to the log file.
    \param s is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain.
	\return true if the proposal is accepted (the target will be split),
    false otherwise.
	\post \a addLeftChild, \a addRightChild may be changed
	indicating whether split operation would result in addition of these
	nodes to the list of proposable nodes.  \a addParent may
	be changed indicating whether merge will affect cherry nodes.
	\a leftChildCount and 
	\a rightChildCount and \a deltaPi may be updated during the operation.*/
    bool decisionMCMCSplitNEW(SPSnode* target,
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
						const std::string& s, int i) const;
	

	/*! \brief Determines whether to merge a node to get a new MCMC state.

    This method probabilistically accepts or rejects a single-step change n the
    histogram state represented by this Adaptive Histogram resulting from
    merging the leaf children of a target cherry node back into the target.

    \param target is a pointer to the target node proposed for merging.
    \param proposal is a reference to a proposal distribution object.
    \param logPrior is a reference to a log prior object.
    \param rgsl is a pointer to a uniform random number generator.
    \param numLeaves is the number of splittable leaf nodes in the SPSnode tree.
    \param numCherries is the number of cherries in the SPSnode tree.
    \param leaves is a collection of RealMappedSPnode leaves, which
	are 'shadowing' the state of this.
	\param minPoints is the minimum number of points allowed in a box.
	\param volChecking is an indicator for whether node volume is 
	being checked to determine whether it is splittable.
    \param minVol is the minimum volume to check for if volChecking is true.
    \param deductLeftChild reference to value indicating if merge will
	result in loss of left child from the list of proposable nodes.
	\param deductRightChild reference to value indicating if merge will
	result in loss of right child from the list of proposable nodes.
    \param addParent reference to value indicating if merge will
	result in addition of parent of target to the list of cherries.
	\param deltaPi reference to a value for the change in the
	log-posterior (log likelihood + log prior)
	that would result from the merge. The value
	of this parameter will be updated during the operation.
	\param logging an enum controlling whether decision making output is
    sent to the log file.
    \param s is the name of the filename to send logging output to.
    \param i is an integer for keeping track of the index for this link in
    a Markov Chain
    \return true if the proposal is accepted (the target will be merged),
    false otherwise.
    \post \a deductLeftChild, \a deductRightChild may be changed
	indicating whether merge operation would result in loss of these
	nodes from the list of proposable nodes.    \a addParent may
	be changed indicating whether merge will affect cherry nodes.
	The value of \a deltaPi will be updated.*/
    bool decisionMCMCMergeNEW(SPSnode* target,
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
						const std::string& s, int i) const;

	
	/*! \brief Changes the state of this Adaptive Histogram by splitting a node.

    This method carries out a move to a new state in the histogram MCMC state
    chain by splitting the target splittable leaf node.

   \pre an SPSnode tree representing some histogram state, a container of
    pointers to the splittable leaf and cherry nodes of the tree in its
    current state, ordered with the leaves first followed by the cherries,
    the number of splittable leaf nodes and the number of cherry nodes.
    \post the SPSnode tree changed as a result of the split, the container of
    splittable leaf and cherry nodes updated for any change in state, and the
    number of splittable leaves and cherries updated similarly.
	
	Throws an UnfulfillableRequest_Error if a split is attempted that
	would result in a box being split beyond the limits to which 
	a floating point number can be represented, so that as far as the
	machine is concerned, the child box == the parent box.  A message
	is also printed out to error output stream cerr. 

    \param target is a pointer to the target node proposed for splitting.
    \param nodes is a reference to a container of pointers to the leaf and
    cherry nodes of the SPSnode tree managed by this AdaptiveHistogram, which
    will be updated by the method.
    \param numLeaves is a reference to a variable storing the number of
    splittable leaves in the SPSnode tree, which will be updated by the method.
    \param numCherries is a reference to a variable storing the number of
    cherries in the SPSnode tree, which will be updated by the method.
    \param addLeftChild indicator for whether split will result
	in the addition of the newLeftChild to the proposable nodes.
	\param addRightChild indicator for whether split will result
	in the addition of new right child to the proposable nodes. 
	\param deductParent indicator for whether split will result
	in the loss of the parent of the target from the cherries. 
	\post the target will be split and \a nodes, \a numLeaves,
	\a numCherries updated as appropriate.*/
    void changeStateForSplitNEW(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries,
						int addLeftChild, int addRightChild,
						int deductParent);
	
	/*! \brief Changes the state of a shadow RealMappedSubPaving
	by splitting a node.
	
	The shadow is the `real' state of this.  Changes in real state
	as a result of
	decisions in the Monte Carlo Markov Chain have to be 
	reflected in the shadow.

    Splits the leaf node in \a leaves that corresponds to \a spsTarget
	and gives the new children `unnormalised height' values based on the
	number of points that will be in them in the histogram (given
	by \a leftChildCount and \a rightChildCount). The
	`unnormalised height' is the count/node volume (ie is not
	normalised for total count of points in the entire histogram).

    \param spsTarget is a pointer to the target node proposed for splitting.
    \param leaves is a reference to a container of pointers to the leaf
	nodes of the shadow, which will be updated by the operation.
	\param cherries is a reference to a container of pointers to the cherry
	nodes of the shadow, which may be updated by the operation.
    \param leftChildCount the number of points
	that will go to the left child as a result of the split. 
	\param rightChildCountthe number of points
	that will go to the right child as a result of the split.
	\param info is a reference to the change of state information 
	object to be updated for this change of state.
    \post the node corresponding to \a target in the shadow will be
	split and the contents of the containers \a leaves and \a cherries 
	will be updated for the split, and \a info will be updated.*/
    static void changeStateForSplitRMSP(
						const SPSnode * const spsTarget,
                        RealMappedSPnode::ListPtrs& leaves, 
						RealMappedSPnode::ListPtrs& cherries, 
						const size_t leftChildCount,
						const size_t rightChildCount,
						AdaptiveHistogram::ChangeOfStateInformation& info);

	
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
    \param deductLeftChild indicator for whether merge will
	result in loss of left child from the list of proposable nodes.
	\param deductRightChild indicator for whether merge will
	result in loss of right child from the list of proposable nodes.
    \param addParent indicator for whether merge will
	result in addition of parent of target to list of cherries.
    \post the target will be merged and \a nodes, \a numLeaves,
	\a numCherries updated as appropriate.*/
    void changeStateForMergeNEW(SPSnode* target,
                        SPSnodeList& nodes, size_t& numLeaves,
                        size_t& numCherries,
						int deductLeftChild, int deductRightChild,
						int addParent);
	
	/*! \brief Changes the state of a shadow RealMappedSubPaving
	by merging a node.
	
	The shadow is the `real' state of this.  Changes in real state
	as a result of
	decisions in the Monte Carlo Markov Chain have to be 
	reflected in the shadow.

    Merges the cherry node in \a cherries that corresponds to \a spsTarget.

    \param spsTarget is a pointer to the target node proposed for splitting.
    \param leaves is a reference to a container of pointers to the leaf
	nodes of the shadow, which will be updated by the operation.
	\param cherries is a reference to a container of pointers to the cherry
	nodes of the shadow, which may be updated by the operation.
	\param info is a reference to the change of state information 
	object to be updated for this change of state.
    \post the node corresponding to \a target in the shadow will be
	merges and the contents of the containers \a leaves and \a cherries 
	will be updated for the merge, and \a info will be updated.*/
    static void changeStateForMergeRMSP(
						const SPSnode * const spsTarget,
                        RealMappedSPnode::ListPtrs& leaves, 
						RealMappedSPnode::ListPtrs& cherries,
						AdaptiveHistogram::ChangeOfStateInformation& info);

	/*! \brief Calculate change in log likelihood for a prospective
	change in state.
	* 
	\internal This method is not delegated to the node concerned because it is 
	more efficient for the histogram to calculate given information
	already gathered.
	
	@todo I'd still like to make this a static method of SPSnode class
	rather than the histogram class: fits better there.  Do this when 
	Gloria's stuff is all sorted out ... ?

    \param leftCount is the number of points prospectively going to 
	a new left child of a node.
    \param rightCount is the number of points prospectively going to 
	a new right child of a node.
    \return the change in log likelihood that would result from the split.*/
	static real getSplitChangeLogLik(
				size_t leftCount, size_t rightCount);

	/*! \brief Internal method to set up priority queue.*/
	std::multiset< SPSnode*, subpavings::MyCompare>&
			_setupPrioritySplit(
				std::multiset< SPSnode*, MyCompare>& pq,
				size_t minChildPoints, double minVol);

	/*! \brief Internal method to do priority split loop.*/
	bool _prioritySplitLoop(std::multiset< SPSnode*, MyCompare>& pq,
							size_t n,
                            size_t minChildPoints, double minVol,
                            gsl_rng * rgsl);

	/*! \brief Internal method to set up priority merge.*/
	std::multiset< SPSnode*, MyCompare>&
			_setupPriorityMerge(
						std::multiset< SPSnode*, MyCompare>& pq);

	/*! \brief Internal method for priority merge loop.*/
	bool _priorityMergeLoop(
							std::multiset< SPSnode*, MyCompare>& pq,
							size_t n,
							gsl_rng * rgsl);

	/*! \brief Internal method to launch priority queue split
	 * to find max posterior point, with logging of posterior.*/
	bool _launchPrioritySplitWithPosterior(
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
                    gsl_rng * rgsl_step);

	/*! \brief Internal method to launch priority queue split, with
	 * logging of posterior.*/
	bool _prioritySplitWithPosterior(
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
					gsl_rng * rgsl);


	/*! \brief Internal method to set up priority queue for split, for
	 * routines that measure empty volume in root box.*/
	cxsc::real _setupPrioritySplitWithEmptyVolMeasure(
					PriorityQueueT& pq,
					SPSNodeMeasure& measurer,
					size_t minChildPoints,
					double minVol) const;
	
	/*! \brief Internal method to reset priority queue for split, for
	 * routines that measure empty volume in root box.*/
	static void _resetPrioritySplitWithEmptyVolMeasure(
		const PriorityQueueT& pqOld,
		PriorityQueueT& pqNew,
		SPSNodeMeasure& measurer);
					
	
	/*! \brief Internal method to do a loop in priority queue for split, for
	 * routines that measure empty volume in root box.
	 
	\internal \a loopEmptyBoxVolume and \a deltaL are updated by reference.*/
	bool _prioritySplitLoopWithEmptyVolMeasure(
					PriorityQueueT& pq,
					SPSNodeMeasure& measurer,
					size_t n, // only needed for emps
					size_t minChildPoints,
					double minVol, 
					real& loopEmptyBoxVolume, real& deltaL, 
					gsl_rng * rgsl);


	/*! \brief Find the coverage value for a data point.
	
	The coverage value is 
	1 - (sum of density of all the boxes with heights >  
	the height of the box where the data point is).
	
	If the point is not in the histogram at all, coverage = 0;
	If the point is in the lowest box in the histogram, 
	coverage = count lowest box / total count;
	If the point is in the highest box of the histogram, coverage = 1
	
	\param pt the point to find coverage for.
	\return the coverage at the point.
	\pre hasSubPaving() == true and the dimensions of the point and the 
	root paving match.  		*/
	double _coverage(const rvector& pt) const;

	/*! \brief Find the empirical density for a data point.
	
	The empirical density is the relative density of the histogram
	at the box the given data point is in.
	
	If the point is not in the histogram at all, empirical density = 0;
	If the point is in the some box in the histogram with height d
	(d = count in box / (total count x box volume), then d is the 
	empirical density.
	
	\param pt the point to find empirical density for.
	\return the empirical density at the point.
	\pre hasSubPaving() == true and the dimensions of the point and the 
	root paving match.
		*/
	double _empiricalDensity(const rvector& pt) const;
	
	/*! \brief Calculate values for variance-covariance matrix
	 * from the big data collection.*/
	RealVec calculateVarCovarFromBigData() const;
	
	/*! \brief Calculate values for mean
	 * from the big data collection.*/
	cxsc::rvector calculateMeanFromBigData() const;
	
	/*! \brief Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/*! \brief Handle exceptions thrown in splitting root to a specific shape. */
	void handleSplitToShapeError(SPSnode& spn);
	
	/*! \brief Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief a constant for padding a box if it is tailor-made for data.

    The padding if the size of the root box is obtained from the min
    and max of the data that is fully fed in.
    */
    static const real padding;
	
	/*! A label for this.*/
	int label;
	
    /*! \brief Pointer to the root node of the subpaving tree.

    An SPSnode is a binary tree representation of a subpaving, designed for
    processing statistical data.
    */
    SPSnode* rootPaving;

    /*! \brief A container for all sample data passed to this.

    The sample that has come in thus far.
    */
    BigDataCollection dataCollection;

    /*! \brief Controls whether all available statistics are maintained in
    the rootPaving.  If set to false (default) only counts are maintained.
    */
    bool holdAllStats;

    
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

	

}; // end of AdaptiveHistogram class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an AdaptiveHistogram object.

	Verbose output for an AdaptiveHistogram object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, const subpavings::AdaptiveHistogram& adh);

} // end namespace subpavings

/*! A specialisation of std::swap for AdaptiveHistogram types.*/
namespace std
{
	template <>
	void swap (subpavings::AdaptiveHistogram & a1, 
			subpavings::AdaptiveHistogram & a2);
}

#endif



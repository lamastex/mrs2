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
\brief AdaptiveHistogramCollator declarations.
*/

#ifndef ___NEW_ADAPTIVEHISTCOLL_HPP__
#define ___NEW_ADAPTIVEHISTCOLL_HPP__


#include "adaptivehistogram.hpp"
#include "collatorspnode.hpp"
#include "sptypes.hpp" // includes a collection of collator node pointers - do we still need?

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

	public:
		/*! \brief Default constructor.
		*/
		AdaptiveHistogramCollator();

		/*! \brief initialised constructor.
		Initialised with an AdaptiveHistogram object.
		
		Throws a NullSubpavingPointer_Error if the pointer to the
		statistical subpaving managed by \a adh is NULL.		*/
		explicit AdaptiveHistogramCollator(const AdaptiveHistogram& adh);

		/*! \brief  Copy constructor
		
		Throws a NullSubpavingPointer_Error if the pointer to the
		collator managed by \a other is NULL.		*/
		AdaptiveHistogramCollator
			(const AdaptiveHistogramCollator& other);


		/*! \brief Copy assignment operator.
		*/
		AdaptiveHistogramCollator&
			operator=(AdaptiveHistogramCollator rhs);



		/*! \brief Destructor.
		*/
		~AdaptiveHistogramCollator();

		//new
		static AdaptiveHistogramCollator importCollator(
					const std::string& s);

		/*! \brief Return a pointer to the CollatorPSnode this manages.
		*/
		CollatorSPnode* getSubPaving() const;
		
		/*! \brief Return true if there is nothing in the collator this manages.
		
		\return True if this' rootCollator has collated no histograms, false
		otherwise.
		*/
		bool isEmptyCollation() const;

		/*! \brief Get the box of the subpaving managed by this.
	
		\note with the present constructors, it is possible for
		this to have a subpaving but for the subpaving to have no box.

		\return copy of the box of the subpaving managed by this.
		\pre !getSubPaving()->isEmpty().*/
		cxsc::ivector getRootBox() const;
	   
		/*! \brief get the dimensions of the subpaving this manages.

		\return 0 if this does not have a subpaving with a box,
		else returns the dimensions of the subpaving.*/
		int getDimensions() const;
		

		/*! Incremental or inplace addition operator.

		Addition gives this histogram collator an expanded tree which represents a
		subpaving which is the union of the subpavings represented by the this's
		original tree and the tree represented by the rhs collator.  The summary
		for each node in the expandd tree contains all the values
		from the summaries of the original nodes and the summaries of the
		corresponding nodes in the tree managed by the rhs collator.
		*/
		AdaptiveHistogramCollator& operator+=(const
					AdaptiveHistogramCollator& rhs);

		/*! Addition operator.

		Addition gives a histogram collator managing a tree which represents a
		subpaving which is the union of the subpavings represented by the operand
		collators.  The summary for each node in the tree contains all the values
		from the summaries of the corresponding nodes in the trees managed by
		the operand AdaptiveHistogramCollators.

		*/
		const AdaptiveHistogramCollator operator+(const
					AdaptiveHistogramCollator& rhs) const;
					
		/*! Incremental or inplace addition operator.

		Addition gives this histogram collator an expanded tree which represents a
		subpaving which is the union of the subpaving represented by the this's
		original tree and the tree represented by rhs as a collator.  The summary
		for each node in the expanded tree contains all the values
		from the summaries of the original nodes and the summaries of the
		corresponding nodes in the tree managed by the collator representation of rhs.
		
		* 
		\warning If an exception is thrown during the addition process, 
		this may be left in an incoherent state.  Users can make 
		a 'backup copy' of a collator before addition 
		if they want to be able to return to the state before 
		the failed addition process.
		*/
		AdaptiveHistogramCollator& operator+=(const
					AdaptiveHistogram& rhs);

		/*! Addition operator.

		Addition gives a histogram collator managing a tree which represents a
		subpaving which is the union of the subpavings represented by this
		collator and the collator representation of rhs.
		The summary for each node in the tree contains all the values
		from the summaries of the corresponding nodes in the trees for this
		collator and the collator representation of rhs.
		*/
		const AdaptiveHistogramCollator operator+(const
					AdaptiveHistogram& rhs) const;

		
		/*! An AdaptiveHistogramCollator which is the average over this collation.

		Makes and returns an AdaptiveHistogramCollator which is the average over
		the collation of histograms represented by this.  The tree managed
		by the average has structure exactly the same as the tree managed by this
		and one value in the summary for each node where that value is
		the average of the summary of the corresponding node in this.

		Throws an UnfulfillableRequest_Error if this has nothing
		collated.
		
		\return An AdaptiveHistogramCollation which is the average of the 
		given collation.
		\pre !isEmptyCollation().
		*/
		const AdaptiveHistogramCollator makeAverage() const;

		
		/*! Make a normalised version of this collation.

		Normalises this collation so that the sum over all the leaf
		nodes of the volume of the box represented by the leaf node and 
		the total summary value of the leaf node (ie, height) is 1
		
		All summaries in the returned collation contain only one 
		summary value.
		
		Throws a UnfulfillableRequest_Error if this has nothing collated.
		
		Throws an std::logic_error if the subpaving managed by this
		has no 'area', ie getTotalAbsValueTimesVol() == 0.
			
		\return The normalised version of this.
		\pre !isEmptyCollation() and getTotalAbsValueTimesVol() != 0.
		*/
		const AdaptiveHistogramCollator makeNormalised() const;
		
		/*! \brief Make a marginalised version of this histogram collator.

		Marginalises to take out the given dimensions and adjust summaries
		so that for each element in the collation, 
		sum of (node vol x node value) 
		is the same as before marginalisation, and hence that the 
		overall sum of (node vol x accumulated values for node)
		is the same as before marginalisation.
				
		\note allowed dimensions start at 1, ie dimensions to
		marginalise on can include 1, 2, ... dimensions of this.
				
		Throws a UnfulfillableRequest_Error if this has nothing collated.
		
		Throws an std::invalid_argument if the required dimensions
		\a reqDim is empty or contains dimensions outside the 
		range of the dimensions of this.
		
		\param reqDims is a vector of the dimensions to include in marginal.
		\return An AdaptiveHistogramCollator managing a subpaving which is
		the marginalised version of the subpaving managed by this.
		\pre \a reqDims must be compatible with current dimensions
			and !isEmptyCollation().
		\post returned histogram will have the same size of collation 
		as before and have sum of 
		(node vol x accumulated summaries) equal to that for this.
		*/
		const AdaptiveHistogramCollator makeMarginal(
									const std::vector<int>& reqDims) const;
		
		
		/*! \brief Find the coverage value for a data point.
		
		The coverage value is 
		1 - (sum of density of all the boxes with heights >  
		the height of the box where the data point is).
		
		Height of a box is taken as the total summary value
		and density is taken as normalised product of 
		height and volume of box.
		
		If the point is not in the histogram at all, coverage = 0;
		If the point is in the lowest box in the histogram, 
		coverage = count lowest box / total count;
		If the point is in the highest box of the histogram, coverage = 1
		
		\warning Coverage only makes sense for collations of positive
		heights (i.e. proper histograms).  The results are unpredictable
		if collations have somehow got negative heights. 
		
		Throws a UnfulfillableRequest_Error if this has nothing collated.
		
		Throws a IncompatibleDimensions_Error if the dimensions
		of this and \a pt are not equal.
				
		\param pt the point to find coverage for
		\return coverage for the point given.	
		\pre !isEmptyCollation() and dimensions of \a pt 
		and this must match.*/
		double findCoverage(const rvector& pt) const;
		
		//new AHABC
		/* fill covNodes with all nodes which represent
		 * the region of this which covers cov of the total histogram density
		 */
		std::vector< const CollatorSPnode* > & findDensityRegion(
								std::vector< const subpavings::CollatorSPnode* > & covNodes,
								double cov) const;

		/*! \brief Find the empirical density for a data point.
		
		The empirical density is the relative density of the histogram
		at the box the given data point is in.
		
		If the point is not in the histogram at all, empirical density = 0;
		If the point is in the some leaf box in the histogram with 'normalised' height d
		(d = total summary for box / sum over leaves of (total summary x box volume)),
		then d is the empirical density.
		
		\warning empirical density only makes sense for collations of positive
		heights (i.e. proper histograms).  The results are unpredictable
		if collations have somehow got negative heights. 
		
		Throws a UnfulfillableRequest_Error if this has nothing collated.
		
		Throws a IncompatibleDimensions_Error if the dimensions
		of this and \a pt are not equal.
		
		\param pt the point to find empirical density for.
		\return the empirical density at the point.
		\pre !isEmptyCollation() and the dimensions of the point and the 
		root paving match.			*/
		double findEmpiricalDensity(const rvector& pt) const;
		
		//new AHABC
		/*! \brief Get whether the root box for this contains a given point.
	 
		Does the support of the collator include the point \a pt?
		
		Throws a UnfulfillableRequest_Error if this has nothing collated.
		
		Throws an IncompatibleDimensions_Error if the given point does not
		have the same dimensions as the subpaving that this manages.

		\param pt the point to check.
		\return true if the root box of the collator subpaving this 
		manages contains \a pt, false otherwise.
		\pre !isEmptyCollation() and the dimensions of the point
		and the root paving match.		*/
		bool histCollatorBoxContains(const rvector& pt) const;

		/*! \brief Split subpaving managed by this to a specified shape.
    
		Used for testing.
		
		Throws a NullSubpavings_Error if the subpaving that this manages
		is a NULL pointer.
		
		Throws a NoBox_Error if the subpaving box is empty.

		Prints a message to the standard error output if the instruction
		could not be carried out. 

		\param instruction specifies the required shape, eg "3, 3, 2, 1"
		\return true if the split was successful, false otherwise
		\pre !isEmptyCollation().*/
		bool splitToShape(std::string instruction);


		
		/*! @name Get a collection of the L1 distance values
		for each element of the collation against the average
		element in the collation.
		 
		The L1 distance for a element against the average is 
		defined over all the leaves of the collator subpaving
		that this manages.  The L1 distance 
		for an element in that collation is the sum,
		over the leaves, of the absolute differences 
		between the 'height' value for the element for 
		the leaf and the 'height' value for the average for the
		leaf, multiplied by the volume
		of the leaf.
		
		Throws an UnfulfillableRequest_Error if this has nothing collated.
		
		\param container a reference to a container to use to 
		store the L1 distance values.  Any contents of the 
		given container will be discarded before new values
		are added.  
		\return An ordered collection of the L1 distance values
		for each element of the collation against the average
		element in the collation, in the same order as the 
		elements are added to this collation .
		\pre !isEmptyCollation().*/
		//@{
		RealVec getL1DistancesToAverage() const;
		
		
		RealVec& getL1DistancesToAverage(RealVec& container) const;
		//@}
		
		/*! @name Get a collection of the L1 distance values
		for each element of this against the average
		element over another collation.
		 
		The L1 distance for a element of the this against the average 
		for another collation is 
		defined over all the leaves of a non-minimal union
		the collator subpavings managed by this and and the other.
		The L1 distance for an element in this is the sum,
		over these leaves, of the absolute differences 
		between the 'height' value for the element for
		the leaf and the 'height' value for the average of 
		the other for the leaf, multiplied by the volume
		of the leaf.
		* 
		Throws the following exceptions:
		<ul>
		<li>Throws a NullSubpavingPointer_Error if the pointer
		to the collator subpaving managed by this or \a other is NULL.</li>
		<li>Throws an UnfulfillableRequest_Error if this has nothing 
		collated.</li>
		<li>Throws an IncompatibleDimensions_Error if the boxes of 
		the subpaving managed by this and the subpaving managed by 
		\a other do not have
		the same dimensions and size.</li>\
		</ul>
		
		\param other a pointer to another collator node:  L1 distances
		are calculated for all the elements of this
		against the average of the colallation pointed to by \a other. 
		\param container a reference to container to use to 
		store the L1 distance values.  Any contents of the 
		given container will be discarded before new values
		are added.  
		\return An ordered collection of the L1 distance values
		for each element of the collation against the average
		element in \a other, in the same order as the 
		elements are added to this collation .
		\pre !other->isEmptyCollation()other and if
		!isEmptyCollation() (this has something collated)
		the boxes for this and \a other should be
		 of equal sizes and dimensions*/
		//@{
		RealVec getL1DistancesToAverage(
				const AdaptiveHistogramCollator& other) const;
		

		// take a container and return the same container, which has been
		// cleared (if necessary) and re-filled with 
		// L1-distances to average-of-other, one for each histogram in collation
		RealVec& getL1DistancesToAverage(
					RealVec& container,
					const AdaptiveHistogramCollator& other) const;
		//@}
			
		/*! @name Get a collection of the L1 distance values
		for each element of the collation against an AdaptiveHistogram.
		 
		The L1 distance for a element of this against a
		statistical subpaving is 
		defined over all the leaves of a non-minimal union
		of the subpaving managed by this and the given statistical
		subpaving.  The L1 distance for an element in this is the sum,
		over these leaves, of the absolute differences 
		between the 'height' value for that element for
		the leaf and the 'height' value of the statistical subpaving
		for the leaf (ie counter/volume normalised by total count 
		in the whole paving), multiplied by the volume
		of the leaf.
		
		Throws a NullSubpavingPointer_Error if the pointer to 
		the subpaving managed by this is NULL or if the pointer to 
		the subpaving managed by \a adh is NULL.
		
		\param adh A pointer to a statistical subpaving to 
		calculate distances against.
		\param container a reference to container to use to 
		store the L1 distance values.  Any contents of the 
		given container will be discarded before new values
		are added.  
		\return An ordered collection of the L1 distance values
		for each element of the collation against \a adh,
		in the same order as the 
		elements are added to this collation.
		\pre the pointers to the subpavings managed by 
		this and by \a adh should be non-NULL.*/
		//@{			
		RealVec getL1Distances(const AdaptiveHistogram& adh) const;
		
		// take a container and return the same container, which has been
		// cleared (if necessary) and re-filled with 
		// L1-distances to the adaptive histogram, one for each histogram in collation
		RealVec& getL1Distances(RealVec& container,
								const AdaptiveHistogram& adh) const;
		//@}

		/*! \brief Add an AdaptiveHistogram object to the collation.

		Attempts to add an AdaptiveHistogram object into the collation
		of AdaptiveHistogram information.
		
		\warning If an exception is thrown during the addition process, 
		this may be left in an incoherent state.  Users can make 
		a 'backup copy' of a collator before adding to the collation 
		if they want to be able to return to the state before 
		the failed addition process.

		\param adh the AdaptiveHistogram to be included in the collation.
		\pre if !isEmptyCollation(), \a adh must have
		the same dimensions as this. 
		\post This will include summary data from the adh.  */
		void addToCollation(const AdaptiveHistogram& adh);

		/*! \brief Add a collection of AdaptiveHistogram objects
		to the collation.

		\warning If an exception is thrown during the addition process, 
		this may be left in an incoherent state.  Users can make 
		a 'backup copy' of a collator before adding to the collation 
		if they want to be able to return to the state before 
		the failed addition process.

		\param samples the collection of AdaptiveHistogram objects
		to be included in the collation.
		\pre everything in \a samples must have the same dimensions as
		each other.  If !isEmptyCollation(), everything
		in \a samples must have the same dimensions as this. 
		\post This will include summary data from the elements of 
		\a samples.  */
		void addToCollation(
					const std::vector < AdaptiveHistogram >& samples);

		/*! \brief Get the number of Adaptive Histogram objects collated.
		*/
		size_t getNumberCollated() const;


		/*! \brief Make a .dot graph file from collated histogram structure.

		Makes a simple .dot graph from the histogram using node names and the
		.png image for this graph.

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
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output.
		*/
		//@{
		void outputAverageToTxtTabs(const std::string& s,
									int prec = 5) const;
									
		void outputAverageToTxtTabs(const std::string& s,
									int prec, bool confirm) const;
		//@}

		/*! \brief Output the collated information to a txt file.

		Output tab delimited data on the collation to a text file.

		\param s the name of the file to send the output to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output..
		*/
		//@{
		void outputToTxtTabs(const std::string& s,
								int prec = 5) const;
		
		void outputToTxtTabs(const std::string& s,
								int prec, bool confirm) const;
		//@}

		/*! \brief Output the subpaving managed by this to a given stream.

		Format is a tab-delimited file of numeric data starting with nodeName, then
		the node box volume, then the node summary, then the description of the
		node box as a tab-delimited list of interval upper and lower bounds.

		\param os is a reference to the stream to output the histogramm to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\return a reference to the given stream.
		*/
		std::ostream & outputToStreamTabs(std::ostream & os,
													int prec = 5) const;

		//new
		/*! \brief Export a description of the collation in a format
		that can be read in again to remake it.
		
		If isEmptyCollation() then the file \a s created will be empty;

		\param s is the name of the file to export to.
		\param prec the precision for output formatting. ie, number
		of decimal places.*/
		void exportCollator(
							const std::string& s, 
							int prec = 5) const;


		/*! \brief Add current state of collation to a log file.

		\param s is the name of the file to log to.
		\param i is a number representing the index of this state in a sequence.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		*/
		void outputLog(const std::string& s, 
									const int i, int prec = 5) const;
		
				
		void swap(AdaptiveHistogramCollator& adh); // throw()

		private:

		/*! \brief Average this. */
		void _average();
		
		/*! \brief Normalise this. */
		void _normalise();
		
		/*! \brief Marginalise this using given required Dimensions.*/
		void _marginalise(const std::vector<int>& reqDims);

		/*! \brief Find the coverage value for a data point.
		
		The coverage value is 
		1 - (sum of density of all the boxes with heights >  
		the height of the box where the data point is).
		
		Height of a box is taken as the total summary value
		and density is taken as normalised product of 
		height and volume of box.
		
		If the point is not in the histogram at all, coverage = 0;
		If the point is in the lowest box in the histogram, 
		coverage = count lowest box / total count;
		If the point is in the highest box of the histogram, coverage = 1
		
		\param pt the point to find coverage for.
		\return the coverage at the point.
		\pre !isEmptyCollation() and the dimensions of the point and the 
		root paving match
			*/
		double _coverage(const rvector& pt) const;

		/*! \brief Find the empirical density for a data point.
		
		The empirical density is the relative density of the histogram
		at the box the given data point is in.
		
		If the point is not in the histogram at all, empirical density = 0;
		If the point is in the some leaf box in the histogram with 'normalised' height d
		(d = total summary for box / sum over leaves of (total summary x box volume)),
		then d is the empirical density.
		
		\param pt the point to find empirical density for.
		\return the empirical density at the point.
		\pre !isEmptyCollation() and the dimensions of the point and the 
		root paving match.
			*/
		double _empiricalDensity(const rvector& pt) const;
		
		/*! \brief Handle exceptions thrown in splitting root to a specific shape. */
		void handleSplitToShapeError(CollatorSPnode& spn);

		/*! \brief Handle exceptions thrown in constructors. */
		void constructor_error_handler();
	
		// data members
		/*! \brief Pointer to the root CollatorSPnode.

		A CollatorSPSnode is a binary tree representation of information from
		a number of subpavings.

		The summary information held by the CollatorSPnode for an SPSnode is the
		count/(total count in tree * volume), ie the normalised histogram height
		for a bin corresponding the box of that SPSnode.
		*/
		CollatorSPnode* rootCollator;

		
	};

	// non-member functions
	/*! \brief Output the contents of an AdaptiveHistogramCollator object.

	Verbose output for an AdaptiveHistogram object.
	*/

	std::ostream & operator<<(std::ostream &os,
				const AdaptiveHistogramCollator& adhc);

}

/*! A specialisation of std::swap for AdaptiveHistogramCollator types.*/
namespace std
{
	template <>
	void swap (subpavings::AdaptiveHistogramCollator & a1, 
			subpavings::AdaptiveHistogramCollator & a2); // throw ()
	
}

#endif



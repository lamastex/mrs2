/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief FunctionEstimatorKDE declarations.
*/

#ifndef ___FUNCTIONESTIMATOR_KDE_HPP__
#define ___FUNCTIONESTIMATOR_KDE_HPP__


#include "piecewise_constant_function.hpp"
#include "realmappedspnode.hpp"
#include "spnode.hpp"
#include "real_pointwise_function_estimate.hpp"

#include <set>      // to use the stl::multiset container
#include <queue> 


namespace subpavings {


/*! \brief A wrapper or manager for an RealMappedSPnode
tree to be used for function estimation based on a 
kernel density estimate (KDE) function available pointwise.
* 
\internal Constructing using a RealPointwiseFunctionEstimator means
that the function is decoupled from the method
of using function to put a value (real, or interval) on a box in the 
estimate.  Different RealPointwiseFunctionEstimators can be used
(which might have the same function to estimate) but do it in a different
way, eg using density at midpoint, or some form of averaging over
a random sample of points in the box, etc. This also has the 
advantage that a  priority queue using \a fobj will be guaranteed
to be consistent with the values that are actually put on any new
leaf nodes formed in the priority queue splitting.*/

class FunctionEstimatorKDE {

    public:
	
	/*! \brief Initialised constructor.

    Initialised  with domain box.
    
	Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

    Ideal constructor when the support domain of the function is
	set a priori.
	\param v The box to use for the subpaving to be managed.
	\param f An estimator for the kde to be represented by this.
	\param lab The label for this (defaults to 0).    
	\post The estimator constructed has subpaving that consists 
	of single leaf node (the root) with a box like \a v) 
	and the range on that single node is the 
	the real image of \a v under the function represented by \a f.*/
    FunctionEstimatorKDE(const ivector& v,
							const RealPointwiseFunctionEstimator& f,
							int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with a subpaving.
    
	\param spn A subpaving to copy as the subpaving to be managed.
	\param f An estimator for the kde to be represented by this.
	\param lab The label for this (defaults to 0).    
	\pre \a spn has a box.
	\post The estimator constructed has label \a lab, 
	a subpaving that that is a copy 
	of spn and the range on each box in that subpaving is 
	the real mid-image of that box
	under the function represented by \a f.*/
    FunctionEstimatorKDE(const SPnode& spn,
							const RealPointwiseFunctionEstimator& f, 
							int lab = 0);

	/*! \brief  Copy constructor.
    */
    FunctionEstimatorKDE(const FunctionEstimatorKDE& other);

	
    //! Destructor
    ~FunctionEstimatorKDE();

	/*! \brief Get the reference to the function object used by this.
	
	\return the function object reference used by this.*/
    const RealPointwiseFunctionEstimator& getFobjReference() const;

	/*! \brief Get the label.
	
	\return the label for this.*/
    int getLabel() const;

	/*! \brief Set the label.*/
    void setLabel(int lab);

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
   
	/*! \brief get the dimensions of the subpaving this manages.

    \return 0 if this does not have a subpaving, else returns the
    dimensions of the subpaving.*/
    int getDimensions() const;
	
	/*! \brief get volume of the root box of the subpaving this manages.

    \return volume of the root box of the subpaving this manages, or
	0.0 if this has no subpaving.*/
    cxsc::real getDomainVolume() const;
    
    /*! \brief Gets number of leaf nodes in the root paving.
    
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return the total number of leaves in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    size_t getRootLeaves() const;

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
	
	/*! \brief Make a PiecewiseConstantFunction out of this estimator.
	  
	 \return A %PiecewiseConstantFunction with the same subpaving
	 as this with the same function values mapped onto the subpaving
	 and the same label as this. */
	subpavings::PiecewiseConstantFunction 
						makePiecewiseConstantFunction() const;



	
	/** @name prioritySplit methods.

    These methods takes an estimator and progressively split using a 
	priority queue of splittable nodes 
	to determine which node to split first.  The ordering for 
	the queue is referred to here as the measure.
	
	\note
	These versions of the prioritySplit function are supplied because 
	it is much more efficient to check for number of leaves directly than
	by using a function object.
	
	The measure of a node 
	used is the total variation between the function 
	estimate on the node and the function estimates on the node's 
	prospective children.  Pieces with the largest
	total variation will be split first.  
	
	Note that this measure <b>does not 
	provide a globally optimal</b> ordering for the queue in the sense
	that it will not necessarily provide the best function estimate
	for a given number of leaves.  
	
	Splitting continues until there are \a maxLeaves leaves,
	or there are no more splittable nodes.

    Each node in the subpaving managed by this decides for itself 
	whether it is splittable, using isSplittableNode().

    If more than one splittable node is equally 'large', on the basis of the
	measure  used, then a random choice is made between all equally
    large nodes to find the node which will be split.

     The seed for the random number generator used for random
	selection between equally
    'large' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple estimates, supply the random number generator seed,
	or your own random number generator,
	to the priority queue to ensure that each estimator
	will make different random choices
	each time; use of the internally set seed will give the same
	results each time.
        
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this is not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Throws an std::logic_error if the split becomes muddled because of 
	some failure within the logic of the algorithm itself.
	 
	Aborts if there are no splittable leaves left (or none at the start).

    \param measure is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioritise splitting.
    \param maxMeasure is the maximum measure value to control splitting.  
	Splitting is stopped when the smallest measure in the queue 
	is > \a maxMeasure.
    \param maxLeaves is number of leaves to aim for in the estimator.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the pseudo-random number generator.
	\param rgsl is a pseudo-random number generator.
    \return true if the priority split was successful, false otherwise
	(returns false if the split could not be started , ie none of the 
	leaves of the initial state is splittable, or if the split
	had to be aborted).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then this has \a maxLeaves leaves;
	if the method returned false then splitting had to be 
	aborted before \a maxLeaves leaves was reached. */
	//@{
	/*! \brief Version where the \a measure is supplied, with a 
	pseudo-random number generator seed.*/
	bool prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	/*! \brief Version where the \a measure is supplied, and
	the \a maxMeasure, with a 
	pseudo-random number generator seed.*/
	bool prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	/*! \brief Version where the \a measure is supplied, with a 
	pseudo-random number generator.*/
	bool prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl);

		/*! \brief Version where the \a measure is supplied, and
	the \a maxMeasure, with a 
	pseudo-random number generator.*/
	bool prioritySplit(
						const RealMappedSPnode::Measurer& measure,
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl);

	/*! \brief Version using the default \a measure, the total variation,
	and a pseudo-random number generator seed.*/
	bool prioritySplit(size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
						
	/*! \brief Version where the default \a measure, the total 
	variation, is used, with 
	the \a maxMeasure supplied.*/
	bool prioritySplit(
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	/*! \brief Version using the default \a measure, the total variation,
	and a pseudo-random number generator is supplied.*/
	bool prioritySplit(
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl);
	
	/*! \brief Version using the default \a measure, the total variation,
	and a \a maxMeasure and a pseudo-random number generator are supplied.*/
	bool prioritySplit( 
						real maxMeasure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						gsl_rng * rgsl);
	
	//@}

	
		
	/*! \brief Split an estimator to a specified shape.
    
   	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving box is empty.

	Prints a message to the standard error output if the instruction
	could not be carried out. 

    \param instruction specifies the required shape, eg "3, 3, 2, 1"
	\return true if the split was successful, false otherwise
    \pre hasSubPaving() == true.
	\post this has subpaving with shape specified by \a instruction
	and the value mapped to each node in the subpaving is the 
	real mid-image of that node's box under the function represented
	by this.*/
    bool splitToShape(std::string instruction);

			
		
	/*! \brief  Get the total integral of the function as estimated
	by this.
	
	The integral is calculated as the sum over
	all the leaves of the subpaving managed by this of the absolute
	value of the real range
	on the leaf multiplied by the volume of the 
	box represented by the leaf.
	
	\return total area between function and 0 for the function as
	estimated by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalIntegralOfRealEstimate() const;
	
	/*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited data giving details of leaf nodes.

    \param os is a reference to the stream to output the estimator to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \return a reference to the given stream.   */
    std::ostream & outputToStreamTabs(std::ostream & os, int prec = 5) const;

	/*! @name Output the subpaving managed by this to a txt file.

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

	/*! @name Output details of full sample (from root) to txt file.

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
	
	/*! @name Output all nodes of the subpaving managed by this to
	a given stream.

    Format is a tab-delimited data giving details of all nodes.

    \param os is a reference to the stream to output to.
	\param prec the precision for output formatting. ie, number
	of decimal places, defaulting to 5.
    \return a reference to the given stream.   */
    std::ostream & outputRootToStreamTabs(std::ostream & os,
													int prec = 5) const;
	
	/*! \brief Append current state of estimator to a txt log file.

    Format is a tab-delimited file of numeric data.
    Output includes node contributions to unscaled EMP under COPERR and AIC
    and the changes in EMP that would result from splitting the node.

    \param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    */
    void outputLog(const std::string& s, int i, int prec = 5) const;
	
	
    /*! \brief Get a string summary of this estimator's properties.
		
	A string description of this. Includes the address of the subpaving
	managed but not details of that subpaving.
	
	\return the string summary.		*/
	std::string stringSummary() const;
	
	
	private:
		
		
		class NodePtrMeasurePair {
			
			public:
				explicit NodePtrMeasurePair(RealMappedSPnode *p, real m);
			
				/* Comparison using measure.*/
				bool operator<(const NodePtrMeasurePair& rhs) const;
				
				std::string toString() const;
			
				RealMappedSPnode * nodePtr;
				cxsc::real measure;
				
			private:
				NodePtrMeasurePair();
			
		};
		
		/*Private inner class to get the total variation, 
		if a node were to be split.*/
		class TotalVariationSplitMeasurer 
					: public RealMappedSPnode::Measurer
		{
			public:
				TotalVariationSplitMeasurer(const RealPointwiseFunctionEstimator& mf);
			
				cxsc::real operator()(
						const RealMappedSPnode * const rmspn) const;
				
			private:
			TotalVariationSplitMeasurer();
			const RealPointwiseFunctionEstimator& measurefobj;
		};

		/*Private inner class to get the increase in 
		total variation if a split node were to be
		merged.*/
		class TotalVariationMergeMeasurer 
					: public RealMappedSPnode::Measurer
		{
			public:
				TotalVariationMergeMeasurer();
			
				cxsc::real operator()(
						const RealMappedSPnode * const rmspn) const;
				
			private:
			
		};
		
		// multiset for the queues
		typedef
		std::multiset<FunctionEstimatorKDE::NodePtrMeasurePair>
		PriorityQueueT;
		typedef PriorityQueueT::iterator PriorityQueueItrT;
	
		typedef
		std::queue< RealMappedSPnode* > NodeQueueT;
		

		
		
	
	/*No argument constructor is private and not implemented.*/
    explicit FunctionEstimatorKDE();

	 /*Copy assignment operator is private an not implemented.    */
    FunctionEstimatorKDE& operator=(const FunctionEstimatorKDE& rhs);


	/*Return a pointer to the RealMappedSPnode 
	this manages. */
    RealMappedSPnode* getSubPaving() const;
	
	
	PriorityQueueT& _setupPrioritySplitQueue(PriorityQueueT& pq,
						const RealMappedSPnode::Measurer& measure);
	
	bool _prioritySplitQueueLoop(
					const RealMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					gsl_rng * rgsl);
	
	/*Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.    */
    void outputLogStart(const std::string& s) const;

	
							
	/*Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/* Handle exceptions thrown in splitting root to a specific shape. */
	void handleSplitToShapeError(RealMappedSPnode& spn);
	
	/*Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief Pointer to the root node of the subpaving tree.

    */
    RealMappedSPnode* rootPaving;

    /*! \brief A representation of the kde function to estimate.   */
    const RealPointwiseFunctionEstimator& fobj;
	
	/*! The label.*/
	int label;

    

}; // end of FunctionEstimatorKDE class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an FunctionEstimatorKDE object.

	Verbose output for an FunctionEstimatorKDE object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, const subpavings::FunctionEstimatorKDE& fek);

} // end namespace subpavings



#endif


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
\brief FunctionEstimatorInterval declarations.
*/

#ifndef __FUNCTIONESTIMATOR_INTERVAL_HPP__
#define __FUNCTIONESTIMATOR_INTERVAL_HPP__

#include "sp_check_visitor.hpp"
#include "piecewise_constant_function.hpp"
#include "intervalmappedspnode.hpp"
#include "sp_value_visitor.hpp"
#include "fei_evalobj.hpp"

#include "mappedFobj.hpp"



#include <set>      // to use the stl::multiset container
#include <queue> 


#include <gsl/gsl_rng.h>        // to know about the gsl random number generator

namespace subpavings {



/*! \brief A wrapper or manager for an IntervalMappedSPnode
tree to be used for function estimation using interval enclosures
over subpaving boxes.
* 
The box of the root of the subpaving managed by 
a %FunctionEstimatorInterval can be thought of as the domain
and the boxes of the leaves of the subpaving can be thought of as
the pieces in the partition of that domain 
for the current function estimate.  The interval on each piece
is an interval enclosure of the function to be estimated
over that piece. 
* 
* @todo should this class? also have the ability to simulate in the mrs
* way, ie squeeze etc: only this one knows about the top of the interval
* or the fobj target.  But we could also use a piecewise constant function?
*/

class FunctionEstimatorInterval {
	
	
	public:
	
	/*! \brief Initialised constructor.

    Initialised  with domain box.
    
	Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

    Ideal constructor when the support domain of the function is
	set a priori.
	\param v The box to use for the subpaving to be managed.
	\param f The function to be estimated by this.
	\param lab The label for this (defaults to 0).    
	\post The estimator constructed has subpaving that consists 
	of single leaf node (the root) with a box like \a v) 
	and the range on that single node is the 
	the interval image of \a v under the function represented by \a f.*/
    FunctionEstimatorInterval(const ivector& v,
										const MappedFobj& f, 
										int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with a subpaving.
    
	\param spn A subpaving to copy as the subpaving to be managed.
	\param f The function to be estimated by this.
	\param lab The label for this (defaults to 0).    
	\pre \a spn has a box.
	\post The estimator constructed has label \a lab, 
	a subpaving that that is a copy 
	of spn and the range on each box in that subpaving is 
	the interval image of that box
	under the function represented by \a f.*/
    FunctionEstimatorInterval(const SPnode& spn, 
								const MappedFobj& f, 
								int lab = 0);


	/*! \brief  Copy constructor.
    */
    FunctionEstimatorInterval(const FunctionEstimatorInterval& other);

   
    //! Destructor
    ~FunctionEstimatorInterval();

	/*! \brief Get the reference to the function object used by this.
	
	\return the function object reference used by this.*/
    const MappedFobj& getFobjReference() const;

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
   
	/*! \brief Get diameter of the interval image of the root box
	of the subpaving managed by this.
	
	\return copy of the box of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
	cxsc::real getRootRangeDiameter() const;
   
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
	 as this with.  The function values mapped onto each node of the
	 for the returned object are the real mid-images under the function
	 represented by this of the box associated with the node. 
	 The %PiecewiseConstantFunction created will have the same label as 
	 this.*/
	subpavings::PiecewiseConstantFunction 
						makePiecewiseConstantFunction() const;

	/*! \brief Create an estimate of the function described by fobj using a 
	depth-first brute force approach.
	
	The over all aim to create an estimate of the described by fobj
	such that each of the leaf nodes of the subpaving managed by this
	has a box associated with it
	that meets the <b>interval image tolerance requirement</b> specified
	by \a tolerance.  

	The interval image tolerance requirement
	is that the diameter of the interval image of the box associated
	with any leaf node of of the subpaving managed by this
	should be less than or equal to the tolerance.
	 
	But in some cases this can result in boxes being split beyond the limits
	of the real number screen.  If a split of a leaf node that does not
	meet the interval image tolerance requirement would result in
	the child nodes' volumes being < the value of the cxsc::MinReal (the 
	smallest representable real number) then that leaf will not be split.

	This means that after the operation, the subpaving managed by this
	may have leaf nodes where the interval image tolerance requirement 
	is not met. 
	
	This operation works depth first, ie working from a leaf root
	it will successively split the root then the root's left child
	and that child's, left child, etc.  It will only start
	splitting the right child of any node already split when all the
	descendents of the left child meet the interval image tolerance
	requirement (or they cannot 
	be split further).  This may cause the machine to run out of memory
	before the process is completed.   
	 
	\note If this operation is applied to this when the subpaving managed
	by this already has leaf nodes which exceed the interval
	image tolerance requirement, those leaf nodes will be left unaltered
	by this operation.
	 
	\param tolerance describes the tolerance to be used in making the 
			estimate.
	\pre hasSubPaving() == true.
	\pre \a tol >= cxsc::MinReal.
	\post the leaf nodes of the subpaving managed by this have a
	range which is the interval image of the box of the node under
	the function represented by this and either the interval image tolerance
	requirement is met for each leaf node or that leaf node is not
	splittable.*/
	void bruteForceEstimate(cxsc::real tolerance);
	
	/*! \brief Create an estimate of the described by fobj using a 
	breadth-first brute force approach with a maximum number of leaves.
	
	The over all aim to create an estimate of the described by fobj
	such that each of the leaf nodes of the subpaving managed by this
	meets whatever requirments are specified
	in \a checker.  If this is not possible given \a maxLeaves, 
	an estimate with at most \a maxLeaves leaves is created which 
	is approximately as close as possible to the requirments of \a checker.
	It is only approximately as close as possible
	because although the method works across the current leaves of the
	same depth or level first before
	descending (ie level by level),
	if the number of leaves reaches \a maxLeaves mid-way
	through a 'level' the process will stop, leaving the remainder of
	the leaves of that level unvisited. 

	After the operation, the subpaving managed by this
	may have leaf nodes which do not meet the requirements 
	specified in \a checker because the maximum number of leaves was
	was reached before those nodes were split. 
	
	This operation works breadth first, ie working from a leaf root
	it will successively split the root then one of the new children
	and then the other, and only then the children of the children.
	The child split first is randomly selected to avoid a strict
	left to right (or right to left) traversal of each level.
		 
	\note If this operation is applied to this when the subpaving managed
	by this already has leaf nodes which exceed the requirements of 
	\a checker, those leaf nodes will be left unaltered
	by this operation.
	 
	\param checker an instance of type %SPCheckVisitor 
	which determines if a leaf node needs to be split further.
	\param maxLeaves the maximum number of leaves the function estimate
	should have at the end of the process.
	\param seed a seed for the random number generator used to determine
	which of a node's children is revisited first.  
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post the leaf nodes of the subpaving managed by this have a
	range which is the interval image of the box of the node under
	the function represented by this and either the 
	estimate has \a maxLeaves leaves and the requirements
	specified by \a checker are \b not met for all leaf nodes or the 
	estimate has less than or equal to \a maxLeaves leaves and
	the requirements
	specified by \a checker are met for all leaf nodes.*/
	
	void breadthFirstBruteForceEstimate(
					const SPCheckVisitor& checker,
					size_t maxLeaves,
					long unsigned int seed = 1234);
	
	void breadthFirstBruteForceEstimate(
					const SPCheckVisitor& nodeChecker,
					const FEIEvalObj& fe,
					long unsigned int seed = 1234);

	void breadthFirstBruteForceEstimate(
					size_t maxLeaves,
					long unsigned int seed = 1234);
	
	void breadthFirstBruteForceEstimate(
					const FEIEvalObj& fe,
					long unsigned int seed = 1234);

	/*! \brief Tighten up the range enclosures of the estimate.

	\pre This must have a subpaving to manage.
	\post The interval range on any non-leaf node of the subpaving managed
	by this will be the interval hull of the ranges of its children.*/
	void hullPropagation();
	
	/*! \brief prioritySplitOnGain.

    These methods takes an estimator and progressively split using a 
	priority queue of splittable nodes 
	to determine which node to split first.  The ordering for 
	the queue is referred to here as the measure.
	
	The measure used is the reduction in the 'area' of the the interval
	band on each piece of the current function estimate resulting 
	from a division of the piece into two halves.  The area of a piece
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  The reduction of this
	area resulting from a split of a piece is (area of band
	for current piece) - (sum of areas of bands for 
	the two halves of the piece that would result from dividing
	that piece in half). Pieces with the largest 'gain' measured in
	this way will be split first. Note that <b>this does not 
	provide a globally optimal</b> ordering for the queue: this ordering
	only looks one step ahead and does not take account of large 
	potential gains that might result from further bisections of the
	current pieces even if the first bisection does not provide a
	large immediate gain.  In particular, if the 
	function over some piece is completely symmetrical in two consecutive
	dimensions then the potential gain on bisection on the first of
	these will be calculated as 0 (because the interval enclosures 
	of both halves of the current piece will be the same as the interval
	enclosure of the current piece itself).  The bivariate Gaussian
	on a symmetric domain box provides a good example of this 
	completely non-globally optimal behaviour.  
	   
	Splitting continues until
    some criteria specified by the function evaluator \a fe
	(applying either to individual nodes or to the estimator
    as a whole) is satisfied, or there are no more splittable nodes.

    Each node in the subpaving managed by this decides for itself 
	whether it is splittable, using isSplittableNode().

    If more than one splittable node is equally 'large', on the basis of the
	measure  used, then a random choice is made between all equally
    large nodes to find the node which will be split.

    The seed for the random number generator used for random
	selection between equally
    'large' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple histograms, supply the random number generator seed
	to the priority queue to ensure that each estimator
	will make different random choices
	each time: use of the internally set seed will give the same
	results each time.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this does not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Throws an std::logic_error if the split becomes muddled because of 
	some failure within the logic of the algorithm itself.
	 
	Aborts if there are no splittable leaves left (or none at the start).

    \param fe is the function evaluator, an instance of a class 
	which provides a criterion to determine when to stop splitting.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.  
    \param seed is a seed for the random number generator.
    \return true if the priority split was successful, false otherwise
	(returns false if the split could not be started , ie none of the 
	leaves of the initial state is splittable, or if the split
	had to be aborted).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then the state of this satisfies
	\a fe; if the method returned false then splitting had to be 
	aborted before \a fe was satisfied. */
	bool prioritySplitOnGain(	const FEIEvalObj& fe,
							LOGGING_LEVEL logging,
							long unsigned int seed = 1234);

	/*! \brief prioritySplitOnGain.

    These methods takes an estimator and progressively split using a 
	priority queue of splittable nodes 
	to determine which node to split first.  The ordering for 
	the queue is referred to here as the measure.
	
	\note
	This overloading of the prioritySplitOnGain function is supplied because 
	it is much more efficient to check for number of leaves directly than
	by using a function object.
	 
	The measure used is the reduction in the 'area' of the the interval
	band on each piece of the current function estimate resulting 
	from a division of the piece into two halves.  The area of a piece
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  The reduction of this
	area resulting from a split of a piece is (area of band
	for current piece) - (sum of areas of bands for 
	the two halves of the piece that would result from dividing
	that piece in half). Pieces with the largest 'gain' measured in
	this way will be split first. Note that <b>this does not 
	provide a globally optimal</b> ordering for the queue: this ordering
	only looks one step ahead and does not take account of large 
	potential gains that might result from further bisections of the
	current pieces even if the first bisection does not provide a
	large immediate gain.  In particular, if the 
	function over some piece is completely symmetrical in two consecutive
	dimensions then the potential gain on bisection on the first of
	these will be calculated as 0 (because the interval enclosures 
	of both halves of the current piece will be the same as the interval
	enclosure of the current piece itself).  The bivariate Gaussian
	on a symmetric domain box provides a good example of this 
	completely non-globally optimal behaviour.  
	   
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
    multiple histograms, supply the random number generator seed
	to the priority queue to ensure that each estimator
	will make different random choices
	each time: use of the internally set seed will give the same
	results each time.
    
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this does not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Throws an std::logic_error if the split becomes muddled because of 
	some failure within the logic of the algorithm itself.
	 
	Aborts if there are no splittable leaves left (or none at the start).

    \param maxLeaves is number of leaves to aim for in the final 
	estimator.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
    \return true if the priority split was successful, false otherwise
	(returns false if the split could not be started , ie none of the 
	leaves of the initial state is splittable, or if the split
	had to be aborted).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then this has \a maxLeaves leaves;
	if the method returned false then splitting had to be 
	aborted before the number of leaves reached \a maxLeaves. */
	bool prioritySplitOnGain(	size_t maxLeaves,
							LOGGING_LEVEL logging,
							long unsigned int seed = 1234);
	
	
	/** @name prioritySplit methods.

    These methods takes an estimator and progressively split using a 
	priority queue of splittable nodes 
	to determine which node to split first.  The ordering for 
	the queue is referred to here as the measure.
	
	The default measure used is the 'area' of the the interval band
	on each piece of the current function estimate, where the area
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  This is effectively the 
	difference of the Reimann integral to the top of the interval 
	enclosure of the box against the Reimann integral to the 
	bottom of the interval enclosure of the box (the
	<b>'Reimann Difference'</b>.  Pieces with the largest
	areas will be split first.  
	
	Note that the Reimann Difference <b>does not 
	provide a globally optimal</b> ordering for the queue in the sense
	that it will not necessarily provide the tightest function enclosure
	for a given number of leaves.  This measure may tend to prioritise
	the bisection of pieces with large volumes but small interval 
	enclosures (ie large pieces where the function is almost flat) over
	pieces where the Reimann Difference is smaller but the function
	is sufficiently variable within the piece for there to be the 
	potential for considerable tightening of the enclosure by 
	bisection.
	
	Splitting continues until
    some criteria specified by the function evaluator \a fe
	(applying either to individual nodes or to the estimator
    as a whole) is satisfied, or there are no more splittable nodes.

    Each node in the subpaving managed by this decides for itself 
	whether it is splittable, using isSplittableNode().

    If more than one splittable node is equally 'large', on the basis of the
	measure  used, then a random choice is made between all equally
    large nodes to find the node which will be split.

     The seed for the random number generator used for random
	selection between equally
    'large' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple estimates, supply the random number generator seed
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
    comparing spsnodes, to order the nodes to prioitise splitting.
    \param fe is the function evaluator, an instance of a class 
	which provides a criterion to determine when to stop splitting.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
    \return true if the priority split was successful, false otherwise
	(returns false if the split could not be started , ie none of the 
	leaves of the initial state is splittable, or if the split
	had to be aborted).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then the state of this satisfies
	\a fe; if the method returned false then splitting had to be 
	aborted before \a fe was satisfied. */
	//@{
	
	/*! \brief Measure for ordering of priority queue supplied by the user,
	and seed for the random number generator supplied by the user.*/	
	bool prioritySplit( const IntervalMappedSPnode::Measurer& measure,
						const FEIEvalObj& fe,
                        LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	
	/*! \brief  Uses the default measure (the Reimann difference
	for each piece of the function estimate) 
	and a seed for the random number generator supplied by the user.  */
	bool prioritySplit( const FEIEvalObj& fe,
                        LOGGING_LEVEL logging,
						long unsigned int seed = 1234);

	//@}
	
	/** @name prioritySplit methods.

    These methods takes an estimator and progressively split using a 
	priority queue of splittable nodes 
	to determine which node to split first.  The ordering for 
	the queue is referred to here as the measure.
	
	\note
	These overloadings of the prioritySplit function are supplied because 
	it is much more efficient to check for number of leaves directly than
	by using a function object.
	
	The default measure used is the 'area' of the the interval band
	on each piece of the current function estimate, where the area
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  This is effectively the 
	difference of the Reimann integral to the top of the interval 
	enclosure of the box against the Reimann integral to the 
	bottom of the interval enclosure of the box (the
	<b>'Reimann Difference'</b>.  Pieces with the largest
	areas will be split first.  
	
	Note that the Reimann Difference <b>does not 
	provide a globally optimal</b> ordering for the queue in the sense
	that it will not necessarily provide the tightest function enclosure
	for a given number of leaves.  This measure may tend to prioritise
	the bisection of pieces with large volumes but small interval 
	enclosures (ie large pieces where the function is almost flat) over
	pieces where the Reimann Difference is smaller but the function
	is sufficiently variable within the piece for there to be the 
	potential for considerable tightening of the enclosure by 
	bisection.
	
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
    multiple estimates, supply the random number generator seed
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
    \param maxLeaves is number of leaves to aim for in the estimator.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
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
	/*! \brief Version where \a measure is supplied by the user.*/	
	bool prioritySplit(
						const IntervalMappedSPnode::Measurer& measure,
						size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);

	/*! \brief Version using the default \a measure, the 'Reimann
	Difference'.*/
	bool prioritySplit(size_t maxLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	//@}

	/** @name priorityMergeOnLoss methods.
	 	
    These methods takes an estimator and progressively merge using a 
	priority queue of mergeable nodes 
	to determine which node to merge first.  The ordering for 
	the queue is referred to here as the measure.
	
	The measure used is the increase in the 'area' of the the interval
	band on each cherry piece of the current function estimate resulting 
	from a merge of two halves of the piece.  This is seen as a loss
	in the sense of being a loss in the tightness of the function
	enclosure, or an increase in the Reimann Difference.  The area of a piece
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  The loss of tightness
	resulting from a merge of two halves of a piece is (area of band
	for merged piece) - (sum of areas of bands for 
	the two halves of the piece). Pieces with the smallest 'loss' measured in
	this way will be merged first. Note that <b>this does not 
	provide a globally optimal</b> ordering for the queue.
	   
	Merging continues until
    some criteria specified by the function evaluator \a fe
	(applying either to individual nodes or to the estimator
    as a whole) is satisfied, or there are no more mergeable nodes.

    If more than one mergeable node is equally 'small', on the basis of the
	measure  used, then a random choice is made between all equally
    small nodes to find the node which will be merged.

     The seed for the random number generator used for random
	selection between equally
    'small' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple estimates, supply the random number generator seed
	to the priority queue to ensure that each estimator
	will make different random choices
	each time; use of the internally set seed will give the same
	results each time.
        
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this is not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Aborts if there are no mergeable cherries left.

    \param fe is the function evaluator, an instance of a class 
	which provides a criterion to determine when to stop merging.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
    \return true if the priority merge was successful, false otherwise
	(returns  false if the merge could not be started , ie there were
	no cherries in the initial state, or if the merge
	had to be aborted before \a fe is satisfied).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then the state of this satisfies
	\a fe; if the method returned false then merging had to be 
	aborted before \a fe was satisfied. */
	//@{
	
	/*! \brief  Uses the loss in tightness of the function enclosure
	as measured by the Reimann Difference on a merge
	as the priority queue ordering measure 
	and a seed for the random number generator supplied by the user.  */
	bool priorityMergeOnLoss(
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	/*! \brief  Uses the loss in tightness of the function enclosure
	as measured by the Reimann Difference on a merge
	as the priority queue ordering measure 
	and a seed for the random number generator supplied by the user.  */
	bool priorityMergeOnLoss(
						size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	//@}

	/** @name priorityMerge methods.
	 	
    These methods takes an estimator and progressively merge using a 
	priority queue of mergeable nodes 
	to determine which node to merge first.  The ordering for 
	the queue is referred to here as the measure.
	
	The default measure used is the 'area' of the the interval band
	on each piece of the current function estimate, where the area
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  This is effectively the 
	difference of the Reimann integral to the top of the interval 
	enclosure of the box against the Reimann integral to the 
	bottom of the interval enclosure of the box (the
	<b>'Reimann Difference'</b>.  Children of cherry node pieces 
	with the smallest areas will be merged first.  
	
	Note that the Reimann Difference <b>does not 
	provide a globally optimal</b> ordering for the queue.
	
	Merging continues until
    some criteria specified by the function evaluator \a fe
	(applying either to individual nodes or to the estimator
    as a whole) is satisfied, or there are no more mergeable nodes.

    If more than one mergeable node is equally 'small', on the basis of the
	measure  used, then a random choice is made between all equally
    small nodes to find the node which will be merged.

     The seed for the random number generator used for random
	selection between equally
    'small' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple estimates, supply the random number generator seed
	to the priority queue to ensure that each estimator
	will make different random choices
	each time; use of the internally set seed will give the same
	results each time.
        
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this is not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Aborts if there are no mergeable cherries left.

    \param measure is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise merging.
    \param fe is the function evaluator, an instance of a class 
	which provides a criterion to determine when to stop merging.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
    \return true if the priority merge was successful, false otherwise
	(returns  false if the merge could not be started , ie there were
	no cherries in the initial state, or if the merge
	had to be aborted before \a fe is satisfied).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then the state of this satisfies
	\a fe; if the method returned false then merging had to be 
	aborted before \a fe was satisfied. */

	//@{
	
	/*! \brief Measure for ordering of priority queue supplied by the user
	and a seed for the random number generator supplied by the user.  */
	bool priorityMerge(const IntervalMappedSPnode::Measurer& measure,
						const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	
	/*! \brief  Uses the Reimann Difference of the cherries
	as the priority queue ordering measure 
	and a seed for the random number generator supplied by the user.  */
	bool priorityMerge(const FEIEvalObj& fe,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	//@}

	/*! \brief Measure for ordering of priority queue supplied by the user
	and a seed for the random number generator supplied by the user.  */
	bool priorityMerge(const IntervalMappedSPnode::Measurer& measure,
						size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);
	
	/** @name priorityMerge methods.
	 	
    These methods takes an estimator and progressively merge using a 
	priority queue of mergeable nodes 
	to determine which node to merge first.  The ordering for 
	the queue is referred to here as the measure.
	
	The default measure used is the 'area' of the the interval band
	on each piece of the current function estimate, where the area
	is the diameter of the interval image of the piece multiplied 
	by the area of the piece (leaf) box.  This is effectively the 
	difference of the Reimann integral to the top of the interval 
	enclosure of the box against the Reimann integral to the 
	bottom of the interval enclosure of the box (the
	<b>'Reimann Difference'</b>.  Children of cherry node pieces 
	with the smallest areas will be merged first.  
	
	Note that the Reimann Difference <b>does not 
	provide a globally optimal</b> ordering for the queue.
	
	Merging continues until
    there are at most \a minLeaves nodes, 
	or there are no more mergeable nodes.

    If more than one mergeable node is equally 'small', on the basis of the
	measure  used, then a random choice is made between all equally
    small nodes to find the node which will be merged.

     The seed for the random number generator used for random
	selection between equally
    'small' nodes can be specified by the user or set by this.
	If you are looking at distributions of results across
    multiple estimates, supply the random number generator seed
	to the priority queue to ensure that each estimator
	will make different random choices
	each time; use of the internally set seed will give the same
	results each time.
        
    Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
	
	Throws an std::logic_error if the state of this is not 'legal',
	ie if this contains cherries that do not pass isSplittableNode().
	
	Aborts if there are no mergeable cherries left.

    \param measure is an instance of a class providing a function for
    comparing spsnodes, to order the nodes to prioitise merging.
    \param minLeaves is the number of leaves to aim for.
    \param logging an enum controlling whether estimator creation output is
    sent to a log file.
    \param seed is a seed for the random number generator.
    \return true if the priority merge was successful, false otherwise
	(returns  false if the merge could not be started , ie there were
	no cherries in the initial state, or if the merge
	had to be aborted before there were \a minLeaves).
	\pre hasSubPaving() == true.
	\pre The state of this is legal ie this does not include cherries
	that should not have been split according to isSplittableNode().
	\post if the method returned true, then this has \a minLeaves 
	leaves; if the method returned false then merging had to be 
	aborted before the number of leaves was reduced to \a minLeaves. */

	//@{
	
	
	/*! \brief  Uses the Reimann Difference of the cherries
	as the priority queue ordering measure 
	and a seed for the random number generator supplied by the user.  */
	bool priorityMerge(size_t minLeaves,
						LOGGING_LEVEL logging,
						long unsigned int seed = 1234);

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
	interval image of that node's box under the function represented
	by this.*/
    bool splitToShape(std::string instruction);

	/*! \brief  Get the total "area" of the interval bands
	estimating the function represented by this.
	
	The "area" of the interval bands is calculated as the sum over
	all the leaves of the subpaving managed by this of the diameter
	of the interval on the leaf multiplied by the volume of the 
	box represented by the leaf.
	
	\return total area of the interval band estimating the function.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalAreaOfIntervalBand() const;

	/*! \brief  Get the maximum "area" of the interval band
	over all the leaves of the subpaving managed by this.
	
	The "area" of the interval band on a leaf is the diameter
	of the interval on the leaf multiplied by the volume of the 
	box represented by the leaf.
	
	\return maximum area of the interval band over all 
	the leaves of the subpaving managed by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getMaximumAreaOfIntervalBand() const;

	/*! \brief  Get the maximum interval tolerance over
	all the leaves of the subpaving managed by this.
	
	For any leaf, the interval tolerance is the maximum distance between 
	the ends of the interval on the leaf and the image of the 
	midpoint of the leaf's box under the function being estimated
	by this.
	
	\return maximum interval tolerance over
	all the leaves of the subpaving managed by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getMaximumIntervalTolerance() const;
	
	/*! \brief  Get the maximum interval diameter over
	all the leaves of the subpaving managed by this.
	
	For any leaf, the interval diameter is the distance between 
	the ends of the interval on the leaf.
	
	\return maximum interval diameter over
	all the leaves of the subpaving managed by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getMaximumIntervalDiameter() const;

    
    /*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited data giving details of leaf nodes.

    \param os is a reference to the stream to output the histogramm to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \return a reference to the given stream.   */
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

    

    /*! \brief Output details of full sample (from root) to txt file.

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
	
	/*! \brief Output all nodes of the subpaving managed by this to
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
				explicit NodePtrMeasurePair(IntervalMappedSPnode *p, real m);
			
				/* Comparison using measure.*/
				bool operator<(const NodePtrMeasurePair& rhs) const;
				
				std::string toString() const;
			
				IntervalMappedSPnode * nodePtr;
				cxsc::real measure;
				
			private:
				NodePtrMeasurePair();
			
		};
		
		/*Private inner class to get the gain, 
		in terms of interval band area reduction,
		if a node were to be split.*/
		class IntervalBandAreaGainOnSplitMeasurer 
					: public IntervalMappedSPnode::Measurer
		{
			public:
				IntervalBandAreaGainOnSplitMeasurer(const MappedFobj& mf);
			
				cxsc::real operator()(
						const IntervalMappedSPnode * const imspn) const;
				
			private:
			IntervalBandAreaGainOnSplitMeasurer();
			const MappedFobj& measurefobj;
		};

		/*Private inner class to get the increase in 
		interval band area if a split node were to be
		merged.*/
		class IntervalBandAreaLossOnMergeMeasurer 
					: public IntervalMappedSPnode::Measurer
		{
			public:
				IntervalBandAreaLossOnMergeMeasurer();
			
				cxsc::real operator()(
						const IntervalMappedSPnode * const imspn) const;
				
			private:
			
		};
		
		// multiset for the queues
		typedef
		std::multiset<FunctionEstimatorInterval::NodePtrMeasurePair>
		PriorityQueueT;
		typedef PriorityQueueT::iterator PriorityQueueItrT;
	
		typedef
		std::queue< IntervalMappedSPnode* > NodeQueueT;
		
	
	/*No argument constructor is private and not implemented.*/
    FunctionEstimatorInterval();

	 /*Copy assignment operator is private an not implemented.    */
    FunctionEstimatorInterval& operator=(const FunctionEstimatorInterval& rhs);

	/*Return a pointer to the IntervalMappedSPnode 
	this manages. */
    IntervalMappedSPnode* getSubPaving() const;
	
	/*Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.    */
    void outputLogStart(const std::string& s) const;

	PriorityQueueT& _setupPrioritySplitQueue(
						PriorityQueueT& pq,
						const IntervalMappedSPnode::Measurer& measure);
						
	bool _prioritySplitQueueLoop(
					const IntervalMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					const SPValueVisitor<cxsc::interval>& estimator,
					gsl_rng * rgsl);

	PriorityQueueT& _setupPriorityMergeQueue(
				PriorityQueueT& pq,
				const IntervalMappedSPnode::Measurer& measure);

	bool _priorityMergeQueueLoop(
					const IntervalMappedSPnode::Measurer& measure,
					PriorityQueueT& pq,
					gsl_rng * rgsl);
	
	/* internal method for breadth first queue setup */
    NodeQueueT& _setupBreadthFirstQueue(
					NodeQueueT& nq,
					const SPCheckVisitor& nodeChecker);

	/* internal method for breadth first queue loop */
    void _breadthFirstQueueLoop(
					const SPCheckVisitor& nodeChecker,
					NodeQueueT& nq,
					const SPValueVisitor<cxsc::interval>& estimator,
					gsl_rng * r);
	
	/*Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/* Handle exceptions thrown in splitting root to a specific shape. */
	void handleSplitToShapeError(IntervalMappedSPnode& spn);
	
	/*Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief Pointer to the root node of the subpaving tree.

    */
    IntervalMappedSPnode* rootPaving;

    /*! \brief A representation of the function to estimate.   */
    const MappedFobj& fobj;
	
	/*! The label.*/
	int label;

    

}; // end of FunctionEstimatorInterval class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an FunctionEstimatorInterval object.

	Verbose output for an FunctionEstimatorInterval object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, const subpavings::FunctionEstimatorInterval& fei);

} // end namespace subpavings



#endif


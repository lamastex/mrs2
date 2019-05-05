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
\brief PiecewiseConstantFunction declarations.
*/

#ifndef ___PIECEWISE_FUNCTION_HPP__
#define ___PIECEWISE_FUNCTION_HPP__

#include "realmappedspnode.hpp"

#include <vector>

#include <gsl/gsl_rng.h>        // to know about the gsl random number generator


namespace subpavings {

	class AdaptiveHistogram;
	class AdaptiveHistogramValidation;

/*! \brief A wrapper or manager for an RealMappedSPnode
tree representing a piecewise constant function.

*/

class PiecewiseConstantFunction {

    public:
	
	/*! \brief No argument constructor.

    \post This has no subpaving and label 0.*/
    explicit PiecewiseConstantFunction();

	
	/*! \brief Initialised constructor.

    Initialised  with domain box and label.
	 
	Constructor to be used if piecewise constant function is to be
	subsequently formed by, for example, splitToShape and allocateRanges.
    
	Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

	\param v The box to use for the subpaving to be managed.
	\param lab The label for this (defaults to 0).    
	\post This has subpaving that consists 
	of single leaf node (the root) with a box like \a v) 
	and the range on that single node is 0.0.*/
    explicit PiecewiseConstantFunction(const ivector& v,
										int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with a subpaving.
    
	\param spn A subpaving to copy as the subpaving to be managed.
	\param lab The label for this (defaults to 0).    
	\pre \a rmspn has a box.
	\post The piecewise constant function constructed has label \a lab, 
	and a subpaving that is a copy 
	of \a rmspn.*/
    PiecewiseConstantFunction(const RealMappedSPnode& rmspn,
							int lab = 0);

	/*! \brief Initialised constructor.

    Initialised  with an AdaptiveHistogram.
    
	\param adh A subpaving to copy as the subpaving to be managed.
	\pre \a adh has a subpaving.
	\post The piecewise constant function constructed has a subpaving 
	that is a copy of the subpaving managed by adh and label the 
	same as the label for \a adh.  The value mapped onto each
	node of this is equal to the count/(volume) of the equivalent node
	of the subpaving managed by \a adh \b after normalising for the 
	total count in the whole histogram, ie the getTotalIntegral() = 1.0.*/
    PiecewiseConstantFunction(const AdaptiveHistogram& adh);
    
    //gat41
    PiecewiseConstantFunction(const AdaptiveHistogramValidation& adh);
    
	/*! \brief  Copy constructor.
    */
    PiecewiseConstantFunction(const PiecewiseConstantFunction& other);

    //! Destructor
    ~PiecewiseConstantFunction();

	/*! \brief Copy assignment operator. */
    PiecewiseConstantFunction& operator=(PiecewiseConstantFunction rhs);


	/*! \brief Get the label.
	
	\return the label for this.*/
    int getLabel() const;

	/*! \brief Set the label.*/
    void setLabel(int lab);
	
	/*! \brief Get a copy of the subpaving managed by this.
	
	\return a copy of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
    const RealMappedSPnode getCopySubPaving() const;
	
	/*! \brief Get whether this has a subpaving to manage.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return true if this has a subpaving to manage.
	false otherwise.*/
    bool hasSubPaving() const;
	
	/*! \brief Get whether the subpaving managed by this 
	 has negative range values.
	
	\return true if the subpaving managed by this has negative ranges,
	false otherwise. Returns false if this has no subpaving to manage.*/
    bool hasNegativePiecewiseConstantValues() const;
	
	/*! \brief Get whether the subpaving managed by this 
	 has infinite range values.
	
	\return true if the subpaving managed by this has infinite ranges,
	false otherwise. Returns false if this has no subpaving to manage.*/
    bool hasInfinitePiecewiseConstantValues() const;
	
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
    
	\return the total number of leaves in the subpaving managed by this.
	Returns 0 if this has no subpaving to manage.*/
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

	
	/*! \brief Split this to a specified shape.
	
	\note that ranges for each node in the revised shape will have to
	be set explicitly, using allocateRanges().
    
   	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving box is empty.

	Prints a message to the standard error output if the instruction
	could not be carried out. 

    \param instruction specifies the required shape, eg "3, 3, 2, 1"
	\return true if the split was successful, false otherwise
    \pre hasSubPaving() == true.
	\post this has subpaving with shape specified by \a instruction 
	and the function value mapped to each node in the subpaving
	is the same as the value mapped to the root box before the operation.*/
    bool splitToShape(std::string instruction);
	
	/*! \brief Recursively allocate a collection of ranges to this and children.
    
    Allocation order is this, left child with remainder of allocation, 
    * right child with remainder.*/
    void allocateRanges(const std::vector< cxsc::real >& rangesToAllocate);
    
	/*! \brief Change this so that the subpaving it manages is
	the union of this's subpaving and the subpaving of another
	%PiecewiseConstantFunction.
	
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
	
	/*! \brief Return a %PiecewiseConstantFunction
	that has subpaving that is
	the union of this's subpaving and the subpaving of another
	%PiecewiseConstantFunction.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer or if the subpaving managed by \a other
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving of this has no box
	or if the the subpaving of \a other has no box.
	
	Throws an IncompatibleDimensions_Error if the
	subpaving boxes of this and \a other are not identical.

	\param other is the %PiecewiseConstantFunction to make the union against.
	\return A %PiecewiseConstantFunction that has a subpaving that is the
	union of the shape of the subpaving managed by this and the shape of 
	the subpaving managed by \a other.  \a this and other are unchanged.
 	\pre Both this and \a other have subpavings with boxes to manage.
	\pre The boxes of the subpavings of this and \a other are the same.*/ 
	PiecewiseConstantFunction makeShapeToUnion(
							const PiecewiseConstantFunction& other);
	
	//src_trubk_0701
	PiecewiseConstantFunction makeShapeToUnionJ(
							const PiecewiseConstantFunction& other) const;


	/*! \brief  Get the total integral of the piecewise constant 
	function represented by this.
	
	The integral is calculated as the sum over
	all the leaves of the subpaving managed by this of the absolute
	value of the real range
	on the leaf multiplied by the volume of the 
	box represented by the leaf.
	* 
	If any of the values on the leaf nodes of the subpaving managed by
	this is Infinity, the value returned will be Infinity. 
	
	\return total area between function and 0 for the function as
	represented by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalIntegral() const;
	
	/*! \brief  Get the integral of a specified box or union of boxes of the piecewise constant 
	function represented by this.
	
	The integral is calculated as the sum over
	all the leaves of the subpaving managed by this of the absolute
	value of the real range on the leaf multiplied by the volume of the 
	box represented by the leaf.
	* 
	If any of the values on the leaf nodes of the subpaving managed by
	this is Infinity, the value returned will be Infinity. 
	
	\return total area between function and 0 for the function as
	represented by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalIntegralForScheffeElement(ivector& box, cxsc::real vol) const;

	/*! \brief  Get the total integrated absolute error (IAE) between
	this and another %PiecewiseConstantFunction.
	
	The IAE is the total of the absolute value of the differences
	in the integrals of this and the 
	other %PiecewiseConstantFunction \a pcf.
	
	\note that this method does not check that either this or \a rmsp
	have the same total integrals nor does it check that either 
	is normalised.
	
	If any of the values on the leaf nodes of the subpaving managed by
	this or \a pcf is Infinity, the value returned will be Infinity. 
	
	This method is symmetric, ie getIAE(pcf) == pcf.getIAE(*this).
	
	\param pcf the %PiecewiseConstantFunction against which to calculate
	the IAE.
	\return total absolute area of the difference between function
	represented by this and function represented by \a pcf.
	\pre this and \a pcf have the same domain box.
	\pre Both this and \a pcf must have a subpaving to manage.*/
	cxsc::real getIAE(const PiecewiseConstantFunction& pcf) const;
	
	//src_trunk_0701
		/*! \brief Get a 'log likelihood' using 
	function values from this and counts from \a spn.
	
	If \b any of the pieces of this have negative values (including -ve infinity)
	then the value returned will be cxsc::SignalingNaN.
	
	Otherwise, if any of the pieces of this <b>where \a adh has points</b> are
	have zero values then the value returned will be -cxsc::Infinity.  
	
	Otherwise, if any of the pieces of this <b>where \a adh has points</b> are
	infinite then the value returned will be cxsc::Infinity.  
	* 
	\note If testing the results, it seems to be more reliable to use
	gsl_isnan and gsl_isinf (both in gsl_math.h) than to use cxsc::IsSignalingNaN()
	and cxsc::IsInfinity().  Also note that both positive and negative 
	cxsc::Infinity will give gsl_isinf true, but negative cxsc::Infinity 
	will also give true for test < 0.0.
	  
	\note This is not normalised before the calculation:  the user is 
	responsible for normalising if that is required. Similarly, any
	pieces of this that have non-positive values will simply be ignored
	in the calculation:  the user is responsible for ensuring that 
	there are no non-positive values if this check is required.
	 
	\param adh is a histogram used to supply the information 
	about counts to be associated with each piece of this
	in calculating the log-likelihood.
	\return the sum over the pieces of this with positive values
	in the intersection of the shape of the partition of this
	and the partion of \a adh of
	the product of the 
	log of the value on the piece of this and the
	count on the corresponding bin of \a adh.  Returns 0 if \a adh
	has no data in it. Return values can include cxsc::SignalingNaN,
	-cxsc::Infinity, and cxsc::Infinity.
	\pre this and \a adh have subpavings to manage and both supavings
	have the same box.*/
	cxsc::real getLogLikelihood(
					const AdaptiveHistogram& adh) const;
					
	//src_trunk_0701	
	/*! \brief  Get maximum value of the constant for any piece in this.
	
	\return maximum value of the constant on any piece of this.
	\pre This has a subpaving to manage.	*/
	cxsc::real getMaxPiecewiseConstant() const;
	
	
	/*! \brief Addition to self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to add to this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of adding the
	values on the equivalent nodes of this before the operation
	and \a rhs.  The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator+= (
							const PiecewiseConstantFunction& rhs);

	/*! \brief Addition operator.
	 
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to add to this.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of adding the
	values on the equivalent nodes of this and \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const PiecewiseConstantFunction operator+ (
							const PiecewiseConstantFunction& rhs) const;
	
	/*! \brief Self-scalar addition operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param add the value to add to this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of adding \a add to the
	values on the equivalent nodes of this.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator+= (const cxsc::real& add);
	
	/*! \brief Scalar addition operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param add the value to add to this.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of adding \a add
	to the values on the equivalent nodes of this. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const PiecewiseConstantFunction operator+ (const cxsc::real& add) const;
	

	/*! \brief Subtraction from self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to subtract from this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of subtracting the
	values on the equivalent nodes of \a rhs from the 
	values on the equivalent nodes of this before the operation.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator-= (
							const PiecewiseConstantFunction& rhs);
    
	/*! \brief Subtraction operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to subtract from this.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of subtracting the
	values on the equivalent nodes of \a rhs from the values on 
	the equivalent nodes of this.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const PiecewiseConstantFunction operator- (
							const PiecewiseConstantFunction& rhs) const;
	
	/*! \brief Self-scalar subtraction operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param sub the value to subtract from this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of subtracting \a sub
	from the values on the equivalent nodes of this.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator-= (const cxsc::real& sub);
	
	/*! \brief Scalar subtraction operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param sub the value to subtract from this.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that is the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of subtracting \a sub
	from the values on the equivalent nodes of this. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const PiecewiseConstantFunction operator- (const cxsc::real& sub) const;


	/*! \brief Multiplication of self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to multiply this by.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of multiplying the
	values on the equivalent nodes of this by the  
	values on the equivalent nodes of \a rhs before the operation.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator*= (
							const PiecewiseConstantFunction& rhs);
	
    /*! \brief Multiplication operator.
    
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to multiply this by.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of multiplying the
	values on the equivalent nodes of this by the values on 
	the equivalent nodes of \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const PiecewiseConstantFunction operator* (
							const PiecewiseConstantFunction& rhs) const;
	
	/*! \brief Self-scalar multiplication operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param mult the value to multiply this by.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of multiplying the
	values on the equivalent nodes of this by \a mult.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator*= (const cxsc::real& mult);
	
	/*! \brief Scalar multiplication operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param mult the value to multiply this by.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of multiplying the
	values on the equivalent nodes of this by \a mult. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const PiecewiseConstantFunction operator* (const cxsc::real& mult) const;
	
	
	/*! \brief Division of self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\note If any of the nodes of the subpaving managed by \a rhs have
	a value of 0.0 on them, then after the division this will
	have nodes with infinite values on them (the result of dividing
	by 0.0).
	
	\param rhs the object to divide this by.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of dividing the
	values on the equivalent nodes of this by the  
	values on the equivalent nodes of \a rhs before the operation.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator/= (
							const PiecewiseConstantFunction& rhs);
	
	/*! \brief Division operator.
    
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\note If any of the nodes of the subpaving managed by \a rhs have
	a value of 0.0 on them, then the %PiecewiseConstantFunction returned
	will manage a subpaving that has nodes with infinite values on them
	(the result of dividing by 0.0).
	
	\param rhs the object to divide this by.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of dividing the
	values on the equivalent nodes of this by the values on 
	the equivalent nodes of \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const PiecewiseConstantFunction operator/ (
							const PiecewiseConstantFunction& rhs) const;
	
	/*! \brief Self-scalar division operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param div the value to divide this by.
	\pre This has a subpaving to manage.
	\pre \a div != 0.0.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of dividing the
	values on the equivalent nodes of this by \a div.
	The label of this is unchanged.*/     
    PiecewiseConstantFunction& operator/= (const cxsc::real& div);
	
	/*! \brief Scalar division operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param div the value to multiply this by.
	\return A %PiecewiseConstantFunction that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of dividing the
	values on the equivalent nodes of this by \a div. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.
	\pre \a div != 0.0.	*/
	const PiecewiseConstantFunction operator/ (const cxsc::real& div) const;
	
	/*! Normalise this.

	Normalises this so that the sum over all the leaf
	nodes of the subpaving managed by this 
	of the product of the volume of the box represented by
	the leaf node and the constant function value on that box
	is 1.0.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws an std::runtime_error if the subpaving managed by this
	has no 'area', ie getTotalIntegral() <= 0.
		
	\pre This has a subpaving to manage and
	getTotalIntegral() > 0.0.
	\post This has the same subpaving as before the operation
	and getTotalIntegral() == 1.0	*/
	void normalise();

	/*! Make a normalised version of this.

	Normalises this so that the sum over all the leaf
	nodes of the subpaving managed by this 
	of the product of the volume of the box represented by
	the leaf node and the constant function value on that box
	is 1.0.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws an std::runtime_error if the subpaving managed by this
	has no 'area', ie getTotalIntegral() <= 0.
		
	\return The normalised version of this.
	\pre This has a subpaving to manage and
	getTotalIntegral() > 0.0.	*/
	const PiecewiseConstantFunction makeNormalised() const;
		
	/*! \brief Make a marginalised version of this.

	Marginalises to take out the given dimensions and adjusts
	the constant values for each node of the subpaving
	managed by this so that the node vol x node value
	is the same as before marginalisation, and hence that the 
	overall sum of (node vol x valu3) over all leaf nodes
	(getTotalIntegral()) is the same as before marginalisation.
			
	\note allowed dimensions start at 1, ie dimensions to
	marginalise on can include 1, 2, ... dimensions of this.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving box is empty.
		
	Throws an std::invalid_argument if the required dimensions
	\a reqDim is empty or contains dimensions outside the 
	range of the dimensions of this.
	
	\param reqDims is a vector of the dimensions to include in marginal.
	\return A PiecewiseConstantFunction managing a subpaving which is
	the marginalised version of the subpaving managed by this.
	\pre This has a subpaving to manage
	and \a reqDims must be compatible with current dimensions.
	\post returned estimator will have sum over leaf nodes of
	(node vol x value) equal to that for this.
	*/
	const PiecewiseConstantFunction makeMarginal(
								const std::vector<int>& reqDims) const;
	
	/* docs 
	 * Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
*/
	const PiecewiseConstantFunction makeSlice(
					const std::vector < int >& sliceDims,
					const std::vector < cxsc::real >& slicePts) const;
	
	/*! \brief Find the coverage value for a data point.
		
	The coverage value refers to the subpaving managed by this and is
	1 - (sum of "area" of all leaf nodes with value >  
	the value of the leaf node whose box contains \a pt)/
	(sum of "area" of all leaf nodes) .
	
	If the point is not in the root box of the subpaving managed by 
	this, coverage = 0.0;
	If the point is in the box of the node associated with
	the lowest value in this, 
	coverage = area lowest node / total area;
	If the point is in the box of the node associated with
	the highest value in this, coverage = 1.0.
	 
	\note The coverage value returned will be in the interval [0,1] 
	whether this is normalised or not (ie it is a adjusted if necessary
	to be a proportion of the total integral).
		
	\warning Coverage only makes sense for piecewise constant functions with 
	non-negative values associated with all boxes (eg density estimates).  
	
	Throws a IncompatibleDimensions_Error if the dimensions
	of this and \a pt are not equal.
			
	\param pt the point to find coverage for
	\return coverage for the point given.	
	\pre This must have a subpaving to manage. 
	\pre Dimensions of \a pt and this must match.
	\pre hasNegativePiecewiseConstantValues() = false.	
	\pre hasInfinitePiecewiseConstantValues() = false.	*/
	cxsc::real findCoverage(const rvector& pt) const;
	
	/*! \brief Find the pointwise extension of the piecewise constant
	 function for a given data point.
	
	The pointwise extension is the value associated with the leaf node
	of the subpaving managed by this
	that has the box containing \a pt.
	 
	\note the pointwise extension is not adjusted if this is not
	normalised (getTotalIntegral() != 1.0): it is just the value
	associated with the leaf node containing \a pt no matter what
	the total integral of the function over the domain is.
	
	If the point is not in the rootbox of the subpaving managed by this
	at all, the pointwise extension is 0;
	Throws a IncompatibleDimensions_Error if the dimensions
	of this and \a pt are not equal.
	
	\param pt the point to find empirical density for.
	\return the empirical density at the point.
	\pre This must have a subpaving to manage. 
	\pre Dimensions of \a pt and this must match.			*/
	cxsc::real pointwiseExtension(
				const rvector& pt) const;
		
	/*! Gets the L1 distance between this and another 
	%PiecewiseConstantFunction.

	The L1 distance is defined as the sum of the absolute values
	of the differences in 'area' represented by this and \a other
	over the union of their partitions (leaf nodes in subpaving
	managed).
	
	Throws the following exceptions:
	<ul>
	<li>Throws a NullSubpaving_Error if either this or 
	by \a other have no subpaving.</li>
	<li>Throws a IncompatibleDimensions_Error if the dimensions
	and sizes of the root boxes of this and \a other are not the same. 
	</ul>
	 
	\note this will not attempt to adjust for any difference in
	total integral between this and \a other: the L1 distance 
	is simply taken as the difference between the 'areas' of
	the leaf boxes.  
	
	\note If this or \a other manages a subpaving in which there
	are leaf nodes with infinite value then the L1 distance between
	 them will be infinite (cxsc::Infinity).  

	\param other the %PiecewiseConstantFunction to calculate
	the L1 distance against.
	\pre Both this and \a other must have subpavings to manage and 
	the root boxes of those subpavings must be the same.
	\post this will be unchanged.       */
	cxsc::real getL1Distance(
				const PiecewiseConstantFunction& other) const;
	
	//20160904 - for KL computations				
	PiecewiseConstantFunction makeSmearZeroValues(
										cxsc::real totalSmear) const;	
	void smearZeroValues(cxsc::real totalSmear);									
	
	/*! @name Output the nodes constituting coverage region \a cov
	to a stream.
	
	The coverage region is the smallest subset of the boxes of the leaf nodes
	of the subpaving managed by this such that the sum of the
	"areas" of those leaf boxes is >= cov * total sum of the "areas"
	of the leaf nodes of the subpaving managed by this.
	
	ie it is the elements of the partition whose integral is at least
	cov% of the integral of the whole when the elements are chosen 
	in order from "tallest" first downwards.  
	
	\pre This must have a subpaving to manage.
	\pre 0.0 =< cov <= 1.0.
	\pre The subpaving managed by this must have no negative ranges,
	 ie hasNegativePiecewiseConstantValues() = false.
	\pre The subpaving managed by this must have no infinite ranges,
	 ie hasInfinitePiecewiseConstantValues() = false.*/
	//@{
	void outputCoverageRegion(	std::ostream & os,
											cxsc::real cov,
											int prec) const;
	
	void outputCoverageRegion(	std::ostream & os,
											cxsc::real cov) const;
	//@}
	
	/*! @name Output the nodes constituting coverage region \a cov
	to a file called \a covFileName.
	
	The coverage region is the smallest subset of the boxes of the leaf nodes
	of the subpaving managed by this such that the sum of the
	"areas" of those leaf boxes is >= cov * total sum of the "areas"
	of the leaf nodes of the subpaving managed by this.
	
	ie it is the elements of the partition whose integral is at least
	cov% of the integral of the whole when the elements are chosen 
	in order from "tallest" first downwards.  
	
	\pre This must have a subpaving to manage.
	\pre 0.0 =< cov <= 1.0.
	\pre The subpaving managed by this must have no negative ranges,
	 ie hasNegativePiecewiseConstantValues() = false.
	\pre The subpaving managed by this must have no infinite ranges,
	 ie hasInfinitePiecewiseConstantValues() = false.*/
	//@{
	void outputCoverageRegion(const std::string& covFileName,
							cxsc::real cov,
							int prec,
							bool confirm = true) const;
	
	void outputCoverageRegion(const std::string& covFileName,
							cxsc::real cov,
							bool confirm = true) const;
	//@}
	
	/*! @name Simulate some data from this.

    Data is simulated from this according to a uniform mixture
	distribution where each leaf box of the subpaving
	managed by this a component in the mixture with weight 
	equal to contribtuion to the total integral over this.
	The result should be equivalent to simulating the 
	by, for each data point, choosing a leaf box
	randomly with probability according to the relative 
	"areas" (value x volume)
	of the leaves and then randomly selecting a data point
	within that box (ie uniform probability of selection over the entire
	box).

    \param container is a reference to a container for the data.
	\param numberToSimulate is the number of data points to simulate.
	\return the reference to \a container.
	\pre This has a subpaving to manage.
	\post \a container contains any data it had before the operation 
	and then \a numberToSimulate data points simulated from this.*/
	//@{
    
	/*! \brief Simulator taking a random number generator argument.
	
	\param r is a random number generator to use for the simulation.*/
    RVecData& simulateData(RVecData& container,
					size_t numberToSimulate, 
					gsl_rng * r) const;
	
	/*! \brief Simulator taking a random number generator seed argument.
	
	Simulation uses the gsl_rng_mt19937 generator and seed \a seed.
	
	\param seed is a random number generator seed to use for the simulation.*/
    RVecData& simulateData(RVecData& container,
					size_t numberToSimulate, 
					long unsigned int seed) const;
	
	/*! \brief Simulator using default random number generator.
	
	Simulation uses the gsl_rng_mt19937 generator and seed = 1234.*/
    RVecData& simulateData(RVecData& container,
					size_t numberToSimulate) const;
	//@}
		
    /*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited data giving details of leaf nodes.

    \param os is a reference to the stream to output the histogramm to.
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
	
	/*! \brief Make a .dot graph file.

    Makes a simple .dot graph from the %piecewise constant function
	using node names and also makes the .png image for this graph.

    \post a .dot file and a .png in the same directory as the program creating
    it was run in.    */
    void outputGraphDot() const;
	
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
	
	void swap(PiecewiseConstantFunction& pcf); // throw()
	
	private:
		
	
	/*Return a pointer to the RealMappedSPnode 
	this manages. */
    RealMappedSPnode* getSubPaving() const;
	
	/*Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.    */
    void outputLogStart(const std::string& s) const;

	/* internal version of normalise */
    void _normalise();
	
	/* internal version of marginalise */
	void _marginalise(
						const std::vector<int>& reqDims);

	/* internal version of coverage for a point */
	cxsc::real _coverage(const rvector& pt) const;

	/* Fill \a covNodes with all nodes which represent
	 * the region of this which covers \a cov % of the total
	 * area of this.*/
	subpavings::RealMappedSPnode::ConstPtrs& findCoverageRegion(
			subpavings::RealMappedSPnode::ConstPtrs& covNodes,
			cxsc::real cov) const;

							
	/*Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/* Handle exceptions thrown changes to root node. */
	void handleSPError(RealMappedSPnode& spn);
	
	/*Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief Pointer to the root node of the subpaving tree.    */
    RealMappedSPnode* rootPaving;

	
	/*! The label.*/
	int label;

    

}; // end of PiecewiseConstantFunction class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an PiecewiseConstantFunction object.

	Verbose output for an PiecewiseConstantFunction object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, 
					const subpavings::PiecewiseConstantFunction& pcf);
								

} // end namespace subpavings

/*! A specialisation of std::swap for PiecewiseConstantFunction types.*/
namespace std
{
	template <>
	void swap (subpavings::PiecewiseConstantFunction & f1, 
			subpavings::PiecewiseConstantFunction & f2); // throw ()
	
}

#endif


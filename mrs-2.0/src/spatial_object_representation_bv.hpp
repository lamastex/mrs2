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
\brief SpatialObjectRepresentationBV declarations.
*/

#ifndef ___SPATIAL_OBJECT_BV_HPP__
#define ___SPATIAL_OBJECT_BV_HPP__

#include "booleanvaluemappedspnode.hpp"
#include "booleanvalue.hpp"

#include <vector>


namespace subpavings {


/*! \brief A wrapper or manager for an BooleanValueMappedSPnode
tree representing a spatial object.

*/

class SpatialObjectRepresentationBV {

    public:
	
	/*! \brief No argument constructor.

    \post This has no subpaving and label 0.*/
    SpatialObjectRepresentationBV();

	
	/*! \brief Initialised constructor.

    Initialised  with domain box and label.
	 
	Constructor to be used if spatial object representation is to be
	subsequently formed by, for example, splitToShape and allocateRanges.
    
	Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

	\param v The box to use for the subpaving to be managed.
	\param lab The label for this (defaults to 0).    
	\post This has subpaving that consists 
	of single leaf node (the root) with a box like \a v) 
	and the range on that single node is false.*/
    explicit SpatialObjectRepresentationBV(const ivector& v,
										int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with a subpaving.
    
	\param spn A subpaving to copy as the subpaving to be managed.
	\param lab The label for this (defaults to 0).    
	\pre \a bmspn has a box.
	\post The spatial object representation constructed has label \a lab, 
	and a subpaving that is a copy 
	of \a bmspn.*/
    explicit SpatialObjectRepresentationBV(const BooleanValueMappedSPnode& bmspn,
							int lab = 0);

	
	/*! \brief  Copy constructor.
    */
    SpatialObjectRepresentationBV(const SpatialObjectRepresentationBV& other);

    //! Destructor
    ~SpatialObjectRepresentationBV();

	/*! \brief Copy assignment operator. */
    SpatialObjectRepresentationBV& operator=(SpatialObjectRepresentationBV rhs);


	/*! \brief Get the label.
	
	\return the label for this.*/
    int getLabel() const;

	/*! \brief Set the label.*/
    void setLabel(int lab);
	
	/*! \brief Get a copy of the subpaving managed by this.
	
	\return a copy of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
    const BooleanValueMappedSPnode getCopySubPaving() const;
	
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
    void allocateRanges(const std::vector< bool >& rangesToAllocate);
    
	/*! \brief Recursively allocate a collection of ranges to this and children.
    
    Allocation order is this, left child with remainder of allocation, 
    * right child with remainder.*/
    void allocateRanges(const std::vector< BooleanMappedValue >& rangesToAllocate);
    
	/*! \brief Change this so that the subpaving it manages is
	the union of this's subpaving and the subpaving of another
	%SpatialObjectRepresentationBV.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer or if the subpaving managed by \a other
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving of this has no box
	or if the subpaving of \a other has no box.
	
	Throws an IncompatibleDimensions_Error if the
	subpaving boxes of this and \a other are not identical.

	There will be no change in this if the subpaving of \other is
	everywhere less split than the subpaving of this.
	
	\param other is the %SpatialObjectRepresentationBV to make the union against.
	\pre Both this and \a other have subpavings with boxes to manage.
	\pre The boxes of the subpavings of this and \a other are the same. 
	\post the subpaving managed by this has the shape that is the
	union of its shape before the operation and the shape of 
	the subpaving managed by \a other.  \a other is unchanged.      */
 	void reshapeToUnion(const SpatialObjectRepresentationBV& other);
	
	/*! \brief Return a %SpatialObjectRepresentationBV
	that has subpaving that is
	the union of this's subpaving and the subpaving of another
	%SpatialObjectRepresentationBV.
	
	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer or if the subpaving managed by \a other
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving of this has no box
	or if the the subpaving of \a other has no box.
	
	Throws an IncompatibleDimensions_Error if the
	subpaving boxes of this and \a other are not identical.

	\param other is the %SpatialObjectRepresentationBV to make the union against.
	\return A %SpatialObjectRepresentationBV that has a subpaving that is the
	union of the shape of the subpaving managed by this and the shape of 
	the subpaving managed by \a other.  \a this and other are unchanged.
 	\pre Both this and \a other have subpavings with boxes to manage.
	\pre The boxes of the subpavings of this and \a other are the same.*/ 
	SpatialObjectRepresentationBV makeShapeToUnion(
							const SpatialObjectRepresentationBV& other) const;


	/*! \brief  Get the total volume of the of the spatial object
	represented by this.
	
	The volume is calculated as the sum over
	all the leaves of the subpaving managed by this with 
	value true of the leaf box volume.
	* 
	\return total volume of the spatial object represented by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalVolume() const;
	
	/*! \brief Union to self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to add to this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of ORing the
	values on the equivalent nodes of this before the operation
	and \a rhs.  The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator+= (
							const SpatialObjectRepresentationBV& rhs);

	/*! \brief Union operator.
	 
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to add to this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of ORing the
	values on the equivalent nodes of this and \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const SpatialObjectRepresentationBV operator+ (
							const SpatialObjectRepresentationBV& rhs) const;
	
	/*! \brief Self-scalar OR operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to OR with this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of ORing \a val with the
	values on the equivalent nodes of this.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator+= (bool val);
	
	/*! \brief Scalar addition operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to OR with this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of ORing \a val
	with the values on the equivalent nodes of this. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const SpatialObjectRepresentationBV operator+ (bool val) const;
	

	/*! \brief XOR (symmetric set difference) against self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to do XOR with this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of doing XOR
	of the values on the equivalent nodes of \a rhs against 
	values on the equivalent nodes of this before the operation.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator-= (
							const SpatialObjectRepresentationBV& rhs);
    
	/*! \brief XOR (symmetric set difference) operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to subtract from this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of doing XOR
	of the values on the equivalent nodes of \a rhs
	against the values on the equivalent nodes of this.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const SpatialObjectRepresentationBV operator- (
							const SpatialObjectRepresentationBV& rhs) const;
	
	/*! \brief Self-scalar XOR (symmetric set difference) operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to use to do set difference against this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of doing the
	XOR of \a val
	against the values on the equivalent nodes of this.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator-= (bool val);
	
	/*! \brief Scalar XOR (symmetric set difference) operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to subtract from this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that is the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of doing XOR
	of \a val
	against the values on the equivalent nodes of this. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const SpatialObjectRepresentationBV operator- (bool val) const;


	/*! \brief Intersection of self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to intersect with this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of ANDing the
	values on the equivalent nodes of this with the  
	values on the equivalent nodes of \a rhs before the operation.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator*= (
							const SpatialObjectRepresentationBV& rhs);
	
    /*! \brief Intersection operator.
    
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to intersect with this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of ANDing the
	values on the equivalent nodes of this with the values on 
	the equivalent nodes of \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const SpatialObjectRepresentationBV operator* (
							const SpatialObjectRepresentationBV& rhs) const;
	
	/*! \brief Self-scalar ANd operator.
   
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param mult the value to AND with this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of ANDing the
	values on the equivalent nodes of this with \a val.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator*= (bool val);
	
	/*! \brief Scalar AND operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param mult the value to AND with this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of ANDing the
	values on the equivalent nodes of this with \a val. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.*/
	const SpatialObjectRepresentationBV operator* (bool val) const;
	
	
	/*! \brief Set difference against self operator.
	
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	
	\param rhs the object to do set difference against this.
	\pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.
	\post this manages a subpaving that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of the subpaving is the result of 
	set differencing the
	values on the equivalent nodes of this with the  
	values on the equivalent nodes of \a rhs before the operation.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator/= (
							const SpatialObjectRepresentationBV& rhs);
	
	/*! \brief  Set difference operator.
    
	Throws a NullSubpavingPointer_Error if either this or \a rhs 
	do not have a subpaving to manage.
	
	Throws an IncompatibleDimensions_Error if the root boxes of 
	this and \a rhs are not identical. 
	 
	\note Note that if this and \a rhs have the same label the result
	of the operation will have that label, but if this and \a rhs
	do not have the same label, the result will get the default label (0).
	
	\param rhs the object to set difference against this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that is the non-minimal union
	of the subpavings managed by this and \a rhs and the value
	on each node of that subpaving is the result of 
	set differencing the
	values on the equivalent nodes of this with the values on 
	the equivalent nodes of \a rhs.  If
	the labels of this and \a rhs are the same, the returned object
	will have the same label as this; otherwise the returned object's
	label will be the default label (0).     
    \pre Both this and \a rhs have subpavings to manage.
	\pre This and \a rhs have identical root boxes.*/
	const SpatialObjectRepresentationBV operator/ (
							const SpatialObjectRepresentationBV& rhs) const;
	
	/*! \brief Self-scalar set difference operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to set difference against this.
	\pre This has a subpaving to manage.
	\post this manages a subpaving that is the same as 
	the subpaving managed by this before the operation and the value
	on each node of the subpaving is the result of
	set differencing the
	values on the equivalent nodes of this with \a val.
	The label of this is unchanged.*/     
    SpatialObjectRepresentationBV& operator/= (bool val);
	
	/*! \brief Scalar division operator.
    
	Throws a NullSubpavingPointer_Error if this  
	does not have a subpaving to manage.
	
	\param val the value to set difference against this.
	\return A %SpatialObjectRepresentationBV that manages a subpaving
	that the same as the subpaving managed by this and the value
	on each node of that subpaving is the result of
	set differencing the
	values on the equivalent nodes of this with \a val. 
	The label of this is unchanged.
    \pre This has a subpaving to manage.	*/
	const SpatialObjectRepresentationBV operator/ (bool val) const;
	
	
	/*! @name Get a slice of this.
	
	\note The \a sliceDims are given as 1, 2, 3, etc, ie are numbered
	from 1 upwards <strong> not from 0 upwards<\strong>.
	
	Make a slice this on at the point jointly specified by \a sliceDims
	and \a slicePts.  For example, if this has a 3-dimensional 
	root box [-1,1]x[-1,1]x[-1,1], ie dimensions {1,2,3}
	and \a sliceDims = {2,3} and 
	\a slicePts is (0.0,0.5) then we are slicing at point 0.0 on 
	dimension 2, point 0.5 on dimension 3. The slice returned 
	would then have only one-dimensional boxes on dimensions 
	{1,2,3}\{2,3} = {1},
	each box being one that did contain point 0.0 on 
	dimension 2, point 0.5 on dimension 3.  The range will be 
	unchanged).  Any boxes 
	associated with this that do not contain 
	the point 0.0 on dimension 2 and 0.5 on dimension 3 will not be
	represented in the returned object.  
	 
	\param sliceDims is a vector of dimensions to slice on, indexed 
	from 1 onwards.
	\param slicePts is a vector of points to slice on, assumed to 
	correspond to the dimensions in \a sliceDims, ie the ith value
	in \a sliceDims gives the dimension for the ith point in
	\a slicePts.
	\param sliceFilename is the name of file to use to capture
	the boxes of this that are used in the slice by outputting 
	them to the file named \a sliceFilename.  Defaults to 
	the empty string "".  If \a sliceFilename is the empty string
	"" no boxes will be captured to file.
	\return The slice of this on the specified dimensions and points.
	\pre This must have a subpaving to manage.
	\pre sliceDims must contain only valid dimensions, ie dimensions
	in this, with minimum dimension 1. */
	//@{
	
	const SpatialObjectRepresentationBV makeSlice(
					const std::vector < int >& sliceDims,
					const std::vector < cxsc::real >& slicePts,
					const std::string sliceFilename = "") const;
	
	const SpatialObjectRepresentationBV makeSlice(
					const std::vector < int >& sliceDims,
					const std::vector < double >& slicePts,
					const std::string sliceFilename = "") const;
	
	//@}
	
	
	
	/*! \brief Find the pointwise extension of the piecewise constant
	 function for a given data point.
	
	The pointwise extension is the value associated with the leaf node
	of the subpaving managed by this
	that has the box containing \a pt.
	 
	If the point is not in the rootbox of the subpaving managed by this
	at all, the pointwise extension is false;
	Throws a IncompatibleDimensions_Error if the dimensions
	of this and \a pt are not equal.
	
	\param pt the point to find pointwise extension for.
	\return the value at the point.
	\pre This must have a subpaving to manage. 
	\pre Dimensions of \a pt and this must match.			*/
	bool pointwiseExtension(
				const rvector& pt) const;
		
		
    /*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited data giving details of leaf nodes.

    \param os is a reference to the stream to output the histogramm to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \return a reference to the given stream.   */
    std::ostream & outputToStreamTabs(std::ostream & os, int prec = 5) const;

	/*! @name Output this representation of the object to a txt file.
	 * 
	Output only includes pieces comprising the object (ie leaf nodes 
	with value true in the 
	subpaving managed by this)

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node value (true), then the description of the
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

	/*! @name Output details of whole representation (from root) to txt file.

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
	
	void swap(SpatialObjectRepresentationBV& pcf); // throw()
	
	private:
		
	
	/*Return a pointer to the BooleanValueMappedSPnode 
	this manages. */
    BooleanValueMappedSPnode* getSubPaving() const;
	
	/*Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.    */
    void outputLogStart(const std::string& s) const;

	/*Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/* Handle exceptions thrown changes to root node. */
	void handleSPError(BooleanValueMappedSPnode& spn);
	
	/*Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief Pointer to the root node of the subpaving tree.    */
    BooleanValueMappedSPnode* rootPaving;

	
	/*! The label.*/
	int label;

    

}; // end of SpatialObjectRepresentationBV class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an SpatialObjectRepresentationBV object.

	Verbose output for an SpatialObjectRepresentationBV object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, 
					const subpavings::SpatialObjectRepresentationBV& pcf);

} // end namespace subpavings

/*! A specialisation of std::swap for SpatialObjectRepresentationBV types.*/
namespace std
{
	template <>
	void swap (subpavings::SpatialObjectRepresentationBV & s1, 
			subpavings::SpatialObjectRepresentationBV & s2); // throw ()
	
}

#endif


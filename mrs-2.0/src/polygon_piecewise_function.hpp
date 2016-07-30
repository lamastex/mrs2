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
\brief PolygonPiecewiseConstantFunction declarations.
*/

#ifndef ___POLYGON_PIECEWISE_FUNCTION_HPP__
#define ___POLYGON_PIECEWISE_FUNCTION_HPP__

#include "realmappedshapenode.hpp"
#include "piecewise_constant_function.hpp"

#include <vector>


namespace subpavings {

	
/*! \brief A wrapper or manager for a collection of 
\link RealMappedShapeNode RealMappedShapeNodes\endlink
representing a polygon 'paving' with real values
mapped to each polygon-shaped piece in the 'paving'.
 
A %PolygonPiecewiseConstantFunction is a representation 
of a piecewise-constant function on a domain that is 
\b not necessarily an axis-aligned hyper-rectangle 
and where the pieces are \b not necessarily axis-aligned
hyper-rectangles.

A %PolygonPiecewiseConstantFunction may be formed by 
transforming a PiecewiseConstantFunction (ie 
a representation of a piecewise-constant
function where the pieces are the pieces of an axis-aligned
regular paving).  The transformations can include rotation, 
reflection, scaling, and translation, ie transformations that
preserve angles. 
* 
A %PolygonPiecewiseConstantFunction also hold its area pieces,
describing the 'floor area' or domain of the function it 
represents (ie have no 'height' or range).  The area pieces 
will describe disjoint areas if this represents a function over
a domain of disjoint parts, but where possible the different
pieces in the domain will have been joined together so that 
the domain is described by as few area pieces as possible.

*/

class PolygonPiecewiseConstantFunction {

    public:
	
	    /*! \brief Constructor with a transformation matrix.
		
		This is constructed by applying a volume and angle-preserving
		back-transformation to the boxes in the PiecewiseConstantFunction
		\a pcf, then scaling as specifed in \a scale,
		then shifting as specified in \a shift.
		
		The back-transformation
		is specified by the matrix represented by \a backtransform.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The polygon-subpaving managed by this 
		will have as many pieces as \a pcf and each 
		piece will be the result of 
		applying the back-transformation matrix \a backtransform
		to the equivalent piece in \a pcf, scaling
		it and shifting it, and then
		giving it a constant value equal to the contstant value on
		the equivalent piece in \a pcf divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the domain that results from scaling the shapes in 
		this).
		
		\param pcf the PiecewiseConstantFunction to use to provide the 
		subpaving structure for this and the constant values to be adjusted 
		for volume scaling to give the values on the pieces of this.
		\param backtransform A representation of a volume preserving
		and angle preserving back-transformation matrix specified
		in row major order, ie 
		the inner vectors of \a backtransform are the rows of 
		the back-transformation matrix.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\pre \a backtransform should represent a square matrix 
		of dimension d where d is also the dimension of 
		the subpaving managed by \a pcf. 
		\pre \a backtransform should be a volume-preserving 
		transformation matrix. 
		\pre \a scale and \a shift should both also be of dimension d. 
		\pre \a scale should contain only strictly positive values. 
		\post This will have the same total integral as \a pcf.*/
        PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift);
		
	    /*! \brief Constructor with no transformation matrix.
		
		This is constructed by scaling \pcf as specifed in \a scale,
		then shifting as specified in \a shift.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The polygon-subpaving managed by this 
		will have as many pieces as \a pcf and each 
		piece will be the result of 
		scaling and shifting the equivalent piece in \a pcf and then
		giving it a constant value equal to the contstant value on
		the equivalent piece in \a pcf divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the domain that results from scaling the shapes in 
		this).
		
		\param pcf the PiecewiseConstantFunction to use to provide the 
		subpaving structure for this and the constant values to be adjusted 
		for volume scaling to give the values on the pieces of this.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\pre \a scale and \a shift should both also be of dimension d
		where d is the dimension of the subpaving managed by \a pcf. 
		\pre \a scale should contain only strictly positive values. 
		\post This will have the same total integral as \a pcf.*/
        PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < double >& scale,
				const std::vector < double >& shift);
		
		/*! \brief Constructor with a transformation matrix and
		a coverage value.
		
		This is constructed by taking the \a cov coverage pieces of 
		\a pcf and applying a volume and angle-preserving
		back-transformation to them, then scaling them as specifed in \a scale,
		then shifting them as specified in \a shift.
		
		The coverage pieces for a given value of \a cov is one of the
		the (possibly non-unique): smallest
		sets of pieces of \a pcf such that the sum of the
		"areas" of those pieces is >= cov * total sum of the "areas"
		of the pieces of \a pcf.  The "area" of a piece of a 
		%PiecewiseConstantFunction is the product of the volume of 
		the box associated with the piece and the value on the piece.
				
		The back-transformation
		is specified by the matrix represented by \a backtransform.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The polygon-subpaving managed by this 
		will have as many pieces as the smallest size of the set
		of pieces of \a pcf such that the set covers \a cov
		of the total integral of \a pcf,  and each 
		piece will be the result of 
		applying the back-transformation matrix \a backtransform
		to the equivalent piece in \a pcf, scaling
		it and shifting it, and then
		giving it a constant value equal to the contstant value on
		the equivalent piece in \a pcf divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the domain that results from scaling the shapes in 
		this).
		
		\param pcf the PiecewiseConstantFunction to use to provide the 
		subpaving structure for this and the constant values to be adjusted 
		for volume scaling to give the values on the pieces of this.
		\param backtransform A representation of a volume preserving
		and angle preserving back-transformation matrix specified
		in row major order, ie 
		the inner vectors of \a backtransform are the rows of 
		the back-transformation matrix.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\param cov The minimum proportion of the total integral of
		\a pcf to be represented in the pieces of \a pcf used to form
		this:  the pieces used to form this will be one of the
		smallest possible sets of pieces of \a pcf representing
		at least \a cov of the total integral of \a pcf.
		\pre \a backtransform should represent a square matrix 
		of dimension d where d is also the dimension of 
		the subpaving managed by \a pcf. 
		\pre \a backtransform should be a volume-preserving 
		transformation matrix. 
		\pre \a scale and \a shift should both also be of dimension d. 
		\pre \a scale should contain only strictly positive values. 
		\post This will have integral >= \a cov multiplied by 
		the total integral of \a pcf.*/
        PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real cov);
				
		/*! \brief Constructor with no transformation matrix but
		specifying a coverage value.
		
		This is constructed by taking the \a cov coverage pieces of 
		\a pcf and then scaling these pieces as specifed in \a scale,
		then shifting them as specified in \a shift.
		
		The coverage pieces for a given value of \a cov is one of the
		the (possibly non-unique): smallest
		sets of pieces of \a pcf such that the sum of the
		"areas" of those pieces is >= cov * total sum of the "areas"
		of the pieces of \a pcf.  The "area" of a piece of a 
		%PiecewiseConstantFunction is the product of the volume of 
		the box associated with the piece and the value on the piece.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The polygon-subpaving managed by this 
		will have as many pieces as the smallest size of the set
		of pieces of \a pcf such that the set covers \a cov
		of the total integral of \a pcf, and each 
		piece will be the result of 
		scaling and shifting the equivalent piece in \a pcf and then
		giving it a constant value equal to the contstant value on
		the equivalent piece in \a pcf divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the domain that results from scaling the shapes in 
		this).
		
		\param pcf the PiecewiseConstantFunction to use to provide the 
		subpaving structure for this and the constant values to be adjusted 
		for volume scaling to give the values on the pieces of this.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\param cov The minimum proportion of the total integral of
		\a pcf to be represented in the pieces of \a pcf used to form
		this:  the pieces used to form this will be one of the
		smallest possible sets of pieces of \a pcf representing
		at least \a cov of the total integral of \a pcf.
		\pre \a scale and \a shift should both also be of dimension d
		where d is the dimension of the subpaving managed by \a pcf. 
		\pre \a scale should contain only strictly positive values.
		\post This will have integral >= \a cov multiplied by 
		the total integral of \a pcf.*/
        PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real cov);

		PolygonPiecewiseConstantFunction(
						const PolygonPiecewiseConstantFunction& other);
		
		//! Destructor
		~PolygonPiecewiseConstantFunction();

		/*! \brief Copy assignment operator. */
		PolygonPiecewiseConstantFunction& operator=(PolygonPiecewiseConstantFunction rhs);

		/*! \brief get the dimensions of the subpaving this manages.

		\return 0 if this does not have a subpaving, else returns the
		dimensions of the subpaving.*/
		int getDimensions() const;
		

		/*! @name Output the pieces of this to a txt file.

		Format is tab-delimited file of numeric data starting with nodeName, then
		the node value, then the vertices of the box.

		\param s the name of the txt file to send output to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output.     */
		//@{
		void outputToTxtTabs(const std::string& s, int prec = 5) const;

		void outputToTxtTabs(const std::string& s, int prec, bool confirm) const;

		//@}

		/*! @name Output the pieces of this to a stream.

		Format is tab-delimited file of numeric data starting with nodeName, then
		the node value, then the vertices of the box.

		\param os the name of the txt file to send output to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output.     */
		std::ostream& outputToStreamTabs(std::ostream& os, int prec = 5) const;
		
		/*! @name Output the area pieces of this to a txt file.
		
		The area pieces are the coordinates representing the 'floor area'
		or domain of this, ie have no 'height' or range.

		Format is tab-delimited file of numeric data starting with 
		nodeName (X), then the node value (0.0), then the vertices of the box.

		\param s the name of the txt file to send output to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output.     */
		//@{
		void outputAreaToTxtTabs(const std::string& s, int prec = 5) const;

		void outputAreaToTxtTabs(const std::string& s, int prec, bool confirm) const;

		//@}

		/*! @name Output the area pieces of this to a stream.

		The area pieces are the coordinates representing the 'floor area'
		or domain of this, ie have no 'height' or range.

		Format is tab-delimited file of numeric data starting with 
		nodeName (X), then the node value (0.0), then the vertices of the box.

		\param os the name of the txt file to send output to.
		\param prec the precision for output formatting. ie, number
		of decimal places.
		\param confirm is a boolean controlling whether confirmation goes to
		console output.     */
		std::ostream& outputAreaToStreamTabs(std::ostream& os, int prec = 5) const;
					
		/*! \brief Swap this and another %PolygonPiecewiseConstantFunction.*/			
		void swap(PolygonPiecewiseConstantFunction& pcf); // throw()
	
	private:
	
		// no-argument constructor 
		PolygonPiecewiseConstantFunction();
		
		void fillNodes(
				const RealMappedSPnode::ConstPtrs& pieces,
				const std::vector < cxsc::ivector >& boxes, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift);
		
		void fillNodes(
				const RealMappedSPnode::ConstPtrs& pieces,
				const std::vector < cxsc::ivector >& boxes, 
				const std::vector < double >& scale,
				const std::vector < double >& shift);
	
		/*Handle exceptions thrown in constructors. */
		void constructor_error_handler();
		
		/*! \brief Pointers to the nodes representing the polygon-pieces
		 * of this.    */
		RealMappedShapeNode::Ptrs nodes;
		
		/*! \brief Pointers to the nodes representing the largest areas 
		 * of the polygon-pieces of this.    */
		RealMappedShapeNode::Ptrs areaNodes;
	   

}; // end of PolygonPiecewiseConstantFunction class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an PolygonPiecewiseConstantFunction object.

	Verbose output for an PolygonPiecewiseConstantFunction object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, 
					const subpavings::PolygonPiecewiseConstantFunction& ppcf);

} // end namespace subpavings

/*! A specialisation of std::swap for PolygonPiecewiseConstantFunction types.*/
namespace std
{
	template <>
	void swap (subpavings::PolygonPiecewiseConstantFunction & p1, 
			subpavings::PolygonPiecewiseConstantFunction & p2); // throw ()
	
}

#endif


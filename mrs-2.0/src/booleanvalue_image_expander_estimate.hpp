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

#ifndef ___BOOLEANVALUE_IMAGE_EXPANDERESTIMATE_HPP__
#define ___BOOLEANVALUE_IMAGE_EXPANDERESTIMATE_HPP__


/*! \file
\brief declarations for BooleanValueImageExpanderEstimator

A type that visits \linkSPnode SPnodes\endlink 
to make a BooleanMappedValue estimate for a function image.

*/


#include "sp_expand_visitor.hpp"
#include "booleanvalue.hpp"

#include "interval.hpp"

namespace subpavings {

    class MappedFobj;

    class BooleanValueImageExpanderEstimator 
				: public SPExpandVisitor<BooleanMappedValue> {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function to be estimated.
			\param tol describes the tolerance to be used in making the 
			estimate.
			\pre \a tol >= cxsc::MinReal.*/
            BooleanValueImageExpanderEstimator(const MappedFobj& f, 
											const cxsc::interval& crit,
											cxsc::real tol);
			
			~BooleanValueImageExpanderEstimator();

			/*! \brief The visit operation.
			
			Checks whether the interval image 
			under the function whose boolean image is to be estimated
			of the box of \a spn
			is contained in, overlaps with, or is completely outside
			the interval criterion of this.
			 
			If the image is outside, returns false;
			If the image is inside, returns true;
			If the image is not completely inside or outside but
			the maximum width of the box of \a spn is < the tolerance of
			this, returns true;
			Else (the image is not completely inside or outside and
			the maximum width of the box of \a spn is >= the tolerance of
			this), splits \a spn and returns true.
			 
			\param spn a pointer to an SPnode to be visited.
			\return BooleanMappedValue(true)
			if there is some overlap between the interval
			image under the function whose boolean image is to be estimated
			of the box of \a spn and the criterion interval of this.
			\post \a spn will have been expanded if the return value 
			is BooleanMappedValue(true)
			 and the maximum width of the box of \a spn
			is >= the tolerance of this.*/
            BooleanMappedValue visit(SPnode * spn) const;

		private:

			BooleanValueImageExpanderEstimator();
			
			const MappedFobj& fobj;
			
			const cxsc::real& supCriterion;
			const cxsc::real& infCriterion;
			
            const cxsc::real tolerance;


            

    };
    // end of BooleanValueImageExpanderEstimator class

} // end namespace subpavings

#endif

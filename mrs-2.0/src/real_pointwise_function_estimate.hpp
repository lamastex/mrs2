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

#ifndef __REAL_POINTWISE_FUNCTION_ESTIMATE_HPP__
#define __REAL_POINTWISE_FUNCTION_ESTIMATE_HPP__


/*! \file
\brief Declarations for RealPointwiseFunctionEstimator

A type that can make a function estimate based a function that may 
only be available pointwise.
* 
This is a type of SPValueVisitor.
* 
\todo need an equivalent for a type to estimate a function for which a
an interval function can be specied - then use this rather than 
the mapped fobject itself in the FunctionEstimatorReal and
FunctionEstimatorInterval types - again decouple function and method
of using function to put a value (real, or interval) on a box in the 
estimate.
*/

/*! \internal
 * This type decouples specific ways to make the function estimate from 
 * the ability to make a function estimate in general.  
 */

#include "sp_value_visitor.hpp"
#include "cxsc.hpp"

#include <vector>

//#include <gls/gsl_qrng.h>

namespace subpavings {

	//namespace kde;
    //class subpavings::kde::TypeKDE;

    class RealPointwiseFunctionEstimator : public SPValueVisitor <cxsc::real> {

        
        public:
	
			
			virtual ~RealPointwiseFunctionEstimator(){}

			virtual cxsc::real operator()(const std::vector < cxsc::real >& x) const = 0;
				
			virtual cxsc::real operator()(const cxsc::rvector& x) const = 0;
			
			virtual cxsc::real operator()(const cxsc::ivector& x) const = 0;

                       

    };
    // end of RealPointwiseFunctionEstimator class

} // end namespace subpavings

#endif

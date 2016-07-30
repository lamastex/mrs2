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
\brief A type for calculating a kernel density estimate given a
*  multivariate data point.

*/

#ifndef __TYPEKDE_HPP__
#define __TYPEKDE_HPP__

	#include "cxsc.hpp"

	#include <vector>

namespace subpavings {
	namespace kde {
		
		/*! \brief A type for calculating a kernel density estimate given a
		*  multivariate data point.		*/
		class TypeKDE {
			
			public :
			
				
				virtual ~TypeKDE(){}
				
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of cxsc::real values.*/
				virtual cxsc::real kde(const std::vector < cxsc::real >& x) const = 0;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of double values.*/
				virtual cxsc::real kde(const std::vector < double >& x) const = 0;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a cxsc::rvector value.*/
				virtual cxsc::real kde(const cxsc::rvector& x) const = 0;
				
				
		}; // end class
	}

} // end namespace


#endif

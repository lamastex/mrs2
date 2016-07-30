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
\brief Declarations for a class for making a kernel density 
* estimate using normal reference rule.
*/

#ifndef __NORMKDE_HPP__
#define __NORMKDE_HPP__

	#include "type_kde.hpp"
	#include "cxsc.hpp"

	#include <vector>
	#include <string>
	#include <cstddef> // size_t

namespace subpavings {
	namespace kde {
		
		/*! \brief A class for making a kernel density 
		* estimate using normal reference rule 
		* (Scott, 'Multivariate Density Estimation: Theory, Practice, and Visualization',
		* 1992, pge 152).

		*/
		class NormKDE  : public TypeKDE {
			
			public :
			
				/*! \brief Constructor using reals.
				 * 
				 * \param dx the collection of sample data from which
				 *  to create the kernel density estimator */
				NormKDE(const std::vector < std::vector < real > >& dx);
				
				/*! \brief Constructor using doubles.
				 * 
				 * \param dx the collection of sample data from which
				 *  to create the kernel density estimator */
				NormKDE(const std::vector < std::vector < double > >& dx);
				
				/*! \brief Destructor. */
				~NormKDE(){}
				
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of cxsc::real values. 
				 * 
				 * \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const std::vector < cxsc::real >& x) const;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of double values. 
				 * 
				 * \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const std::vector < double >& x) const;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a cxsc::rvector value. 
				 * 
				 * \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const cxsc::rvector& x) const;
				
				/*! \brief Estimate bandwidths.
				 * 
				 * Uses the normal reference rule to estimate bandwidths.
				 * 
				 * \param logFilename the name of a file to which to send
				 * logging output used to record bandwidths and diagnostic
				 * information.
				 * \return the collection of bandwidths used.  These 
				 * are the entries on the diagonal of the diagonal
				 * bandwidth matrix created by this routine.*/
				std::vector < cxsc::real > df(	
									const std::string& logFilename) const ;
				
			private :

				cxsc::real _kde(const std::vector < cxsc::real >& rx) const;
								
				
				void rescaleData() const;
				
				std::vector < std::vector < real > > data_x;
				size_t dim;
				size_t data_num;
				cxsc::real cont; 
				mutable std::vector < cxsc::real > bandwidths;
				mutable cxsc::real productBandwidths;
				mutable std::vector < std::vector < real > > rescaled_x;
				
		}; // end class
	}

} // end namespace


#endif

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

#ifndef ___REAL_KDE_MID_ESTIMATE_HPP__
#define ___REAL_KDE_MID_ESTIMATE_HPP__


/*! \file
\brief declarations for RealKDEMidEstimator

A type that can visit \linkSPnode SPnodes\endlink 
to make an estimate for a KDE function.

The estimate supplied is the kernel density estimate evaluated at the 
mid-point of the box associated with the node.
*/

#include "real.hpp"
#include "type_kde.hpp"
#include "real_pointwise_function_estimate.hpp"

//#include <gls/gsl_qrng.h>

namespace subpavings {

	
	/*! \brief An estimator for a Kernel Density Estimate that estimates
	the kde on a box as the kde evaluated 
	at the mid-point of the box.*/
    class RealKDEMidEstimator : 	public RealPointwiseFunctionEstimator {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a kde to be estimated.*/
			RealKDEMidEstimator(const subpavings::kde::TypeKDE& f);
			
			~RealKDEMidEstimator();

			/*! \brief The visit operation.
			
			\param spn a pointer to an SPnode to be visited.
			\return the real kernel density estimate 
			evaluated at the midpoint of the box associated with the node 
			pointed to by \a spn.  */
            cxsc::real visit(SPnode * spn) const;
			
			
			/*! \brief Evaluate the kernel density estimate at single point.
			
			\param x the point at which to evaluate the kernel density,
			given as a vector of reals rather than an rvector.
			\return the real kernel density estimate 
			evaluated at the rvector corresponding to \a x.  */
            cxsc::real operator()(const std::vector < cxsc::real >& x) const;
				
			/*! \brief Evaluate the kernel density estimate at single point.
			
			\param x the point at which to evaluate the kernel density.
			\return the real kernel density estimate 
			evaluated at \a x.  */
            cxsc::real operator()(const cxsc::rvector& x) const;
			
			/*! \brief Estimate of the function value 
			over a box \a x.
			
			\param x the box over which to estimate the kernel density.
			\return the real kernel density estimate 
			evaluated at the mid-point of \a x.  */
			cxsc::real operator()(const cxsc::ivector& x) const;


		private:

			RealKDEMidEstimator();
			
			const subpavings::kde::TypeKDE& fobj;
			
			#if(0)
			const size_t n;
			
			gsl_qrng * q;
			
			
					if (!num) throw std::invalid_argument(
				"RealKDEMidEstimator(...) : num = 0");
					q = gsl_qrng_alloc (gsl_qrng_sobol, di,);



				if ((NULL != q) gsl_qrng_free (q);
				q = NULL;

			/* the average of the real kernel density estimate 
			evaluated at \a n points in the box associated with the node 
			pointed to by \a spn.  The points chosen will always include
			the midpoint and will also include \a n-1 points from
			a quasi-random sequence in the box.*/
            

			#endif
			
                       

    };
    // end of RealKDEMidEstimator class

} // end namespace subpavings

#endif

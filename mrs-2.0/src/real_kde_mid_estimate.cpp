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
\brief definitions for RealKDEMidEstimator


*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "real_kde_mid_estimate.hpp"

#include "spnode.hpp"


#include <stdexcept>

using namespace cxsc;

//#define MYDEBUG
//#define MYDEBUG_MIN

#ifdef NDEBUG
	#undef MYDEBUG
	#undef MYDEBUG_MIN
#endif

#if defined (MYDEBUG) || defined (MYDEBUG_MIN)
#define MYDEBUG_MIN
#include <iostream>
#endif

namespace subpavings {


    RealKDEMidEstimator::RealKDEMidEstimator(const subpavings::kde::TypeKDE& f)
                : fobj(f)
	{}
		
	
	RealKDEMidEstimator::~RealKDEMidEstimator()
	{}


    cxsc::real RealKDEMidEstimator::visit(SPnode * spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in RealKDEMidEstimator::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
        ivector box = spn->getBox();
		
		#ifdef MYDEBUG
			std::cout << "this box is " << box << std::endl;
        #endif
		
		real thisRange = this->operator()(box);
		
		#ifdef MYDEBUG
			std::cout << "this range at midpoint is " << thisRange << std::endl;
		#endif
				
		return thisRange;
    }
	



	cxsc::real RealKDEMidEstimator::operator()(const std::vector < cxsc::real >& x) const
	{
		return fobj.kde(x);
	}
		
	cxsc::real RealKDEMidEstimator::operator()(const cxsc::rvector& x) const
	{
		return fobj.kde(x);
	}

	// return the mid-image
	cxsc::real RealKDEMidEstimator::operator()(const cxsc::ivector& x) const
	{
		int lb = Lb(x);
		int dim = lb + Ub(x) + 1;
		std::vector < cxsc::real > rv(dim);
		for (int i = 0; i < dim; ++i) rv[i] = mid(x[lb+i]); 
		return fobj.kde(rv);
			
	}
} // end namespace subpavings

#if(0)
    cxsc::real RealKDEMidEstimator::visit(SPnode * spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in RealKDEMidEstimator::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
        ivector box = spn->getBox();
		
		#ifdef MYDEBUG
			std::cout << "this box is " << box << std::endl;
        #endif
		
		real thisRange(0.0);
		
		{
			int lb = Lb (ivec);

			std::vector < cxsc::real > x(dim);
			for (size_t i = 0; i < dim; ++i) x[i] = cxsc::mid(ivec[lb+i]);
			
			thisRange += fobj.kde(x);
		}
		
		#ifdef MYDEBUG
			std::cout << "this range at midpoint is " << thisRange << std::endl;
		#endif
		
		/* if n > 1 do quasi random points as well */
		if (n > 1) {
			
			/* want the sequence to be quasi random for this box itself */
			gsl_qrng_init (q);
			
			for (size_t i = 0; i < n-1; ++i) {
				double v[dim];

				gsl_qrng_get (q, v);
				std::vector < cxsc::real > x(v, v+dim);
				thisRange += fobj.kde(x);
			}
		}
			
		return thisRange/n;
    }
#endif

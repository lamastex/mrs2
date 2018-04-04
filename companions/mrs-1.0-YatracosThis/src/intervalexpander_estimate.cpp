/*
* Copyright (C) 2010, 2011, 2012 Jennifer Harlow
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
\brief definitions for IntervalExpanderEstimator


*/




#include "intervalexpander_estimate.hpp"

#include "mappedFobj.hpp"

#include "spnode.hpp"

#include "toolz.hpp" // for realVolume

#include <stdexcept>


using namespace cxsc;

//#define MYDEBUG
//#define MYDEBUG_MIN

#if defined (MYDEBUG) || defined (MYDEBUG_MIN)
#define MYDEBUG_MIN
#include <iostream>
#endif

namespace subpavings {


    IntervalExpanderEstimator::IntervalExpanderEstimator(const MappedFobj& f,
														cxsc::real tol)
                : fobj(f), tolerance(tol) 
	{
		if (tolerance < cxsc::MinReal) 
			throw std::invalid_argument(
				"IntervalExpanderEstimator::IntervalExpanderEstimator(MappedFobj&, cxsc::real) : tol < cxsc::MinReal");
	}
		
	
	IntervalExpanderEstimator::~IntervalExpanderEstimator() {}


    cxsc::interval IntervalExpanderEstimator::visit(SPnode * spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in IntervalExpanderEstimator::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
        // check if we need to split
        ivector box = spn->getBox();
		
		#ifdef MYDEBUG
			std::cout << "this box is " << box << std::endl;
        #endif
		
		interval thisRange = fobj(box);
		real thisMidImage = fobj.imageMid(box);
		
		#ifdef MYDEBUG
			std::cout << "this midImage is " << thisMidImage << std::endl;
		#endif
		
		#ifdef MYDEBUG
			std::cout << "this range is " << thisRange << std::endl;
		#endif
		
		#ifdef MYDEBUG
			std::cout << "this tolerance is " << tolerance << std::endl;
		#endif
		
		// split if so
		if (max(Sup(thisRange) - thisMidImage, thisMidImage - Inf(thisRange)) > tolerance) {
			spn->nodeExpand();

		}
		#ifdef MYDEBUG_MIN
		if (!(max(Sup(thisRange) - thisMidImage, thisMidImage - Inf(thisRange)) > tolerance)) {
		
			std::cout << "no need to expand further: interval image <= tolerance " << std::endl;
		}
		#endif
		
		return thisRange;
    }

    // end of IntervalExpanderEstimator class

} // end namespace subpavings



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
\brief definitions for BooleanValueImageExpanderEstimator


*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "booleanvalue_image_expander_estimate.hpp"

#include "mappedFobj.hpp"

#include "spnode.hpp"

#include "toolz.hpp" // for maxdiam

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


    BooleanValueImageExpanderEstimator::BooleanValueImageExpanderEstimator(
												const MappedFobj& f,
												const cxsc::interval& crit,
												cxsc::real tol)
                : fobj(f), supCriterion(Sup(crit)), infCriterion(Inf(crit)),
					tolerance(tol) 
	{
		if (tolerance < cxsc::MinReal) 
			throw std::invalid_argument(
				"BooleanValueImageExpanderEstimator::BooleanValueImageExpanderEstimator(MappedFobj&, cxsc::real) : tol < cxsc::MinReal");
	}
		
	
	BooleanValueImageExpanderEstimator::~BooleanValueImageExpanderEstimator() {}


    BooleanMappedValue BooleanValueImageExpanderEstimator::visit(SPnode * spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in BooleanValueImageExpanderEstimator::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
        // check if we need to split
        ivector box = spn->getBox();
		
		#ifdef MYDEBUG
			std::cout << "this box is " << box << std::endl;
        #endif
		
		#ifdef MYDEBUG
			std::cout << "this tolerance is " << tolerance << std::endl;
		#endif
		
		interval thisRange = fobj(box);
		
		#ifdef MYDEBUG
			std::cout << "this range is " << thisRange << std::endl;
		#endif
		
		
		bool retValue = false;
		
		if ( !(Sup(thisRange) < infCriterion) 
						&& !(Inf(thisRange) > supCriterion) ) { // overlap
			retValue = true;
			#ifdef MYDEBUG_MIN
				std::cout << "overlap, returning true" << std::endl;
			
			#endif
			
			if ( (Sup(thisRange) > supCriterion) 
				|| (Inf(thisRange) < infCriterion) ) { //indet
				
				int maxdiamcomp;
				double maxDiam = MaxDiam (box, maxdiamcomp);
				// check if box max width is < tolerance, expand if not
				if (!(maxDiam < tolerance)) {
					spn->nodeExpand();
					#ifdef MYDEBUG_MIN
						std::cout << "and max width >= tolerance " << std::endl;
					
					#endif
				}
				
			}
			// else thisRange is inside criterion				
		}
		
		
		return BooleanMappedValue(retValue);
    }

    // end of BooleanValueImageExpanderEstimator class

} // end namespace subpavings



/*
* Copyright (C)2012 Jennifer Harlow
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
\brief definitions for BooleanValueImageEstimator


*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "booleanvalue_image_estimate.hpp"

#include "mappedFobj.hpp"

#include "spnode.hpp"



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


    BooleanValueImageEstimator::BooleanValueImageEstimator(const MappedFobj& f,
											const cxsc::interval& crit)
                : fobj(f), supCriterion(Sup(crit)), infCriterion(Inf(crit))
	{}
		
	
	BooleanValueImageEstimator::~BooleanValueImageEstimator() {}


    BooleanMappedValue BooleanValueImageEstimator::visit(SPnode * spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in BooleanValueImageEstimator::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
        ivector box = spn->getBox();
		
		#ifdef MYDEBUG
			std::cout << "this box is " << box << std::endl;
        #endif
		
		interval thisRange = fobj(box);
			
		#ifdef MYDEBUG
			std::cout << "this range is " << thisRange << std::endl;
		#endif
		
		return BooleanMappedValue( !(Sup(thisRange) < infCriterion) 
				&& !(Inf(thisRange) > supCriterion) );
    }
	
	// end of BooleanValueImageEstimator class

} // end namespace subpavings



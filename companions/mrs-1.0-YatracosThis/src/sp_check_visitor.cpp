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
\brief definitions for SPCheckVisitor type and concrete subtypes.


*/




#include "sp_check_visitor.hpp"

#include "spnode.hpp"

#include "toolz.hpp"

#include <stdexcept>


using namespace cxsc;

//#define MYDEBUG
//#define MYDEBUG_MIN

#if defined (MYDEBUG) || defined (MYDEBUG_MIN)
#define MYDEBUG_MIN
#include <iostream>
#endif

namespace subpavings {


    // SplittableCheck class
	SplittableCheck::SplittableCheck() {}
		
	
	SplittableCheck::~SplittableCheck() {}


    void SplittableCheck::visit(const SPnode * const spn) const
    {
		result  = spn->isSplittableNode();
    }
	
	bool SplittableCheck::getResult() const
    {
		return !result;
	} 

    // end of SplittableCheck class
	
	
	// IntervalImageToleranceCheck class
	IntervalImageToleranceCheck::IntervalImageToleranceCheck(const MappedFobj& f,
														cxsc::real tol)
                : fobj(f), tolerance(tol) 
	{
		if (tolerance < cxsc::MinReal) 
			throw std::invalid_argument(
				"IntervalImageToleranceCheck::IntervalImageToleranceCheck(MappedFobj&, cxsc::real) : tol < cxsc::MinReal");
	}
		
	
	IntervalImageToleranceCheck::~IntervalImageToleranceCheck() {}


    void IntervalImageToleranceCheck::visit(const SPnode * const spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in IntervalImageToleranceCheck::visit, for " << spn->getNodeName() << std::endl;
		#endif
		
		if (spn->isSplittableNode()) {
		
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
			
			result = true;
			if (max(Sup(thisRange) - thisMidImage, 
				thisMidImage - Inf(thisRange) )
							> tolerance) {
				result = false;

			}
			
		}
		#ifdef MYDEBUG_MIN
			std::cout << "conclusion " << result << std::endl;
		
		#endif
		
    }
	
	bool IntervalImageToleranceCheck::getResult() const
    {
		return result;
	} 

    // end of IntervalImageToleranceCheck class
	

	
	// ReimannDiffToleranceCheck class
    ReimannDiffToleranceCheck::ReimannDiffToleranceCheck(const MappedFobj& f,
														cxsc::real tol)
                : fobj(f), tolerance(tol) 
	{
		if (tolerance < cxsc::MinReal) 
			throw std::invalid_argument(
				"ReimannDiffToleranceCheck::ReimannDiffToleranceCheck(MappedFobj&, cxsc::real) : tol < cxsc::MinReal");
	}
		
	
	ReimannDiffToleranceCheck::~ReimannDiffToleranceCheck() {}


    void ReimannDiffToleranceCheck::visit(const SPnode * const spn) const
    {
		#ifdef MYDEBUG
			std::cout << "in ReimannDiffToleranceCheck::visit, for " << spn->getNodeName() << std::endl;
		#endif
			
		result = true;
		if (spn->isSplittableNode()) {
		
			ivector box = spn->getBox();
			
			#ifdef MYDEBUG
				std::cout << "this box is " << box << std::endl;
			#endif
			
			interval thisRange = fobj(box);
			
			#ifdef MYDEBUG
				std::cout << "this range is " << thisRange << std::endl;
			#endif
			
			#ifdef MYDEBUG
				std::cout << "this tolerance is " << tolerance << std::endl;
			#endif
			
			if ((realVolume(box) * cxsc::diam(thisRange)) > tolerance) {
				result = false;

			}
			
		}
		#ifdef MYDEBUG_MIN
			std::cout << "conclusion " << result << std::endl;
		
		#endif
		
    }
	
	bool ReimannDiffToleranceCheck::getResult() const
    {
		return result;
	} 

    // end of ReimannDiffToleranceCheck class

} // end namespace subpavings



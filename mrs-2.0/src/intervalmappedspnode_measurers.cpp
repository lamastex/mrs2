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
\brief Definitions for types that are concrete subclasses of
* IntervalMappedSPnode::Measurer.
*/

#include "intervalmappedspnode_measurers.hpp"




using namespace subpavings;
using namespace std;


ReimannDiffMeasurer::ReimannDiffMeasurer() {}
			
cxsc::real ReimannDiffMeasurer::
			operator()(const IntervalMappedSPnode * const imspn) const
{
	return ( imspn->getIntervalAreaOnBox() );
}


IntervalToleranceMeasurer::IntervalToleranceMeasurer(const MappedFobj& f)
 : fobj(f) {}
			
cxsc::real IntervalToleranceMeasurer::
			operator()(const IntervalMappedSPnode * const imspn) const
{
	
	interval thisRange = imspn->getRange();
	real thisMidImage = fobj.imageMid(imspn->getBox());
			
	 return max( cxsc::abs(Sup(thisRange) - thisMidImage), 
				cxsc::abs(thisMidImage - Inf(thisRange) ) );
}


			
cxsc::real IntervalWidthMeasurer::
			operator()(const IntervalMappedSPnode * const imspn) const
{
	
	return cxsc::diam(imspn->getRange());
	
}

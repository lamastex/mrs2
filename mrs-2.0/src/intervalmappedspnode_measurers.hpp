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
\brief Declarations for types that are concrete subclasses of
* IntervalMappedSPnode::Measurer.
*/

#ifndef ___INTERVALMAPPEDSPNODE_MEASURER_HPP__
#define ___INTERVALMAPPEDSPNODE_MEASURER_HPP__

#include "intervalmappedspnode.hpp"
#include "mappedFobj.hpp"



namespace subpavings {


	/*! \brief A type for measuring nodes that uses the total
	 * area of the interval band on a node, ie the volume of the box
	 * multiplied by the diameter of the interval.
	*/

	class ReimannDiffMeasurer 
				: public IntervalMappedSPnode::Measurer
	{
		public:
			ReimannDiffMeasurer();
		
			cxsc::real operator()(
					const IntervalMappedSPnode * const imspn) const;
			
		private:
	};

	/*! \brief A type for measuring nodes that uses the 
	 * maximum distance between the ends of the interval 
	 * and real image of the midpoint of the box of the node
	 * under some function.
	*/

	class IntervalToleranceMeasurer 
				: public IntervalMappedSPnode::Measurer
	{
		public:
			
			explicit IntervalToleranceMeasurer(const MappedFobj& f);
		
			cxsc::real operator()(
					const IntervalMappedSPnode * const imspn) const;
			
		private:
			IntervalToleranceMeasurer();
			
			const MappedFobj& fobj;
	};

	/*! \brief A type for measuring nodes that uses the 
	 * width of the interval as the measure.
	*/

	class IntervalWidthMeasurer 
				: public IntervalMappedSPnode::Measurer
	{
		public:
			
			cxsc::real operator()(
					const IntervalMappedSPnode * const imspn) const;
		
	};



} // end namespace subpavings



#endif


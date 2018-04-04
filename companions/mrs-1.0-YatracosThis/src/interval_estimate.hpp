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

#ifndef ___INTERVALESTIMATE_HPP__
#define ___INTERVALESTIMATE_HPP__


/*! \file
\brief declarations for IntervalEstimator

A type that visits \linkSPnode SPnodes\endlink 
to make an estimate for a function.

The estimate supplied is the interval image of the box associated
with the node.
*/

#include "interval.hpp"

#include "sp_value_visitor.hpp"


namespace subpavings {

    class MappedFobj;

    class IntervalEstimator : public SPValueVisitor <cxsc::interval> {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function to be estimated.*/
			IntervalEstimator(const MappedFobj& f);
			
			~IntervalEstimator();

			/*! \brief The visit operation.
			
			\param spn a pointer to an SPnode to be visited.
			\return the interval image under the function to be estimated by this
			of the box associated with the node pointed
			to by \a spn.*/
            cxsc::interval visit(SPnode * spn) const;

		private:

			IntervalEstimator();
			
			const MappedFobj& fobj;
                       

    };
    // end of IntervalEstimator class

} // end namespace subpavings

#endif

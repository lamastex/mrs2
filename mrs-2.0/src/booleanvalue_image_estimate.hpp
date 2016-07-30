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

#ifndef ___BOOLEANVALUEIMAGE_ESTIMATE_HPP__
#define ___BOOLEANVALUEIMAGE_ESTIMATE_HPP__


/*! \file
\brief declarations for BooleanValueImageEstimator

A type that visits \linkSPnode SPnodes\endlink 
to make an estimate for a function.

The estimate supplied is a
BooleanMappedValue indicating whether the interval image of the box associated
with the node is at least partially overlaps the interval range criterion
of this.
*/

#include "sp_value_visitor.hpp"
#include "booleanvalue.hpp"
#include "interval.hpp"

namespace subpavings {

    class MappedFobj;

    class BooleanValueImageEstimator
						: public SPValueVisitor<BooleanMappedValue> {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function to be estimated.*/
			BooleanValueImageEstimator(const MappedFobj& f,
									const cxsc::interval& crit);
			
			~BooleanValueImageEstimator();

			/*! \brief The visit operation.
			
			\param spn a pointer to an SPnode to be visited.
			\return BooleanMappedValue indicating whether
			the interval image under the function to be estimated by this
			of the box associated with the node pointed
			to by \a spn is inside or partially inside (ie overlaps with)
			the interval criterion of this.*/
            BooleanMappedValue visit(SPnode * spn) const;

		private:

			BooleanValueImageEstimator();
			
			const MappedFobj& fobj;
			
			const cxsc::real& supCriterion;
			const cxsc::real& infCriterion;
                       

    };
    // end of BooleanValueImageEstimator class

} // end namespace subpavings

#endif

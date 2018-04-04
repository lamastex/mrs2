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

#ifndef ___INTERVALEXPANDERESTIMATE_HPP__
#define ___INTERVALEXPANDERESTIMATE_HPP__


/*! \file
\brief declarations for IntervalExpanderEstimator

A type that visits \linkSPnode SPnodes\endlink 
to make an estimate for a function.

The overall aim is to create an estimate 
of a function described by \a f such that, for any node visited
by this, each of the leaf nodes descended from that node 
(or the node itself if it is a leaf) has a box associated with it
that meets the <b>interval image tolerance requirement</b> of the
estimator.  

The interval image tolerance requirement of the estimator 
is defined by its tolerance, given by the constructor argument \a tol.
For any node visited by the estimator, the interval image tolerance requirement
is that the maximum distance between the real image of the mid point 
of the  box associated
with the node and the ends of the interval image of the box associated
with the node should be less than or equal to the tolerance.

(But in some cases this can result in boxes being split beyond the limits
of the real number screen.  For this (and possibly other reasons) the
node may be unwilling to be split and refuse to accept a visit from this).

If the node is visited by this, the node will be split 
if the maximum distance between the real image of the 
mid-point of the box associated
with the node and the ends of the interval image of the box associated
with the node is greater than the \a tolerance property of the 
esimator.

If a node where the interval image tolerance requirement 
is not met has not accepted the visit that node will not have 
been visited and will not have been expanded further. 

The estimator's visit operation should return the interval image under the 
function to be estimated of the box associated with the node visted.
*/

#include "interval.hpp"

#include "sp_expand_visitor.hpp"

namespace subpavings {

    class MappedFobj;

    class IntervalExpanderEstimator : public SPExpandVisitor<cxsc::interval> {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function to be estimated.
			\param tol describes the tolerance to be used in making the 
			estimate.
			\pre \a tol >= cxsc::MinReal.*/
            IntervalExpanderEstimator(const MappedFobj& f, cxsc::real tol);
			
			~IntervalExpanderEstimator();

			/*! \brief The visit operation.
			 
			Expands the node pointed to by \a spn if the  
			interval image tolerance requirement is not met and returns
			the interval image under the function to be estimated by this
			of the box associated with the node pointed
			to by \a spn.
			
			\param spn a pointer to an SPnode to be visited.
			\return the interval image under the function to be estimated by this
			of the box associated with the node pointed
			to by \a spn.*/
            cxsc::interval visit(SPnode * spn) const;

		private:

			IntervalExpanderEstimator();
			
			const MappedFobj& fobj;
            const cxsc::real tolerance;


            

    };
    // end of IntervalExpanderEstimator class

} // end namespace subpavings

#endif

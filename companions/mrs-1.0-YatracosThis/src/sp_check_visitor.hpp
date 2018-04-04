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

#ifndef __SP_CHECK_HPP__
#define __SP_CHECK_HPP__

#include "mappedFobj.hpp"
#include "real.hpp"

/*! \file
\brief declarations for SPCheckVisitor type and concrete subtypes.

A interface for a type that visits \linkSPnode SPnodes\endlink 
and returns some boolean check information for them.*/



namespace subpavings {

    class SPnode;

	class SPCheckVisitor {

        
        public:
	
			virtual ~SPCheckVisitor(){};
	
			/*! \brief The visit operation.*/
            virtual void visit(const SPnode * const spn) const = 0;
			
			/*! \brief Get the result of the visit operation.
			
			This should be true if the last node checked meets all 
			of the requirements of this, false otherwise.*/
            virtual bool getResult() const = 0;


    };
    // end of SPCheckVisitor class
	
	
	/*! \brief A type that visits \linkSPnodes SPnodes\endlink 
	to check if they are splittable.

	The result of the visit should be a boolean indicating whether
	the node should be split, according its own isSplittabeNode()
	method.*/
	class SplittableCheck : public SPCheckVisitor {

        
		public:

			/*! \brief Constructor.*/
			SplittableCheck();
			
			~SplittableCheck();

			/*! \brief The visit operation.
			 
			Checks the node pointed to by \a spn and records the result
			true if the node is splittable (ie if spn->isSplittableNode()).
			
			\param spn a pointer to an SPnode to be visited.
			\return true if \a spn is not splittable, false otherwise.*/
			void visit(const SPnode * const spn) const;
			
			/*! \brief Get the result of the visit operation.*/
			bool getResult() const;

		private:

			mutable bool result;


			

	};
	// end of SplittableCheck class
	
	/*! \brief A type that visits \linkSPnodes SPnodes\endlink 
	to check if they meet an interval image tolerance
	requirement.

	The interval image tolerance requirement
	is that the maximum distance between the real image of the mid point 
	of the  box associated
	with the node and the ends of the interval image of the box associated
	with the node should be less than or equal to the tolerance given
	by the constructor argument \a tol.

	A node may not meet the requirement but not be splittable: in this
	case the result of the visit should be true (ie the node is 
	deemed to meet the requirement because it is as close as
	it can be to meeting it).

	The visit operation should return true if 
	either the the node visited meets the
	interval image tolerance requirement or the node does not meet
	the requirement but is not splittable.  Return false otherwise.
	*/
	class IntervalImageToleranceCheck : public SPCheckVisitor {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function against which to check
			the interval image tolerance requirement.
			\param tol describes the tolerance to be used in the check.
			\pre \a tol >= cxsc::MinReal.*/
            IntervalImageToleranceCheck(const MappedFobj& f, cxsc::real tol);
			
			~IntervalImageToleranceCheck();

			/*! \brief The visit operation.
			 
			Checks the node pointed to by \a spn and records the
			result true if the  
			interval image tolerance requirement is met.
			
			\param spn a pointer to an SPnode to be visited.
			\return true if \a spn meets the interval image
			tolerance requirement or the node does not meet
			the requirement but is not splittable.  
			Return false otherwise.*/
            void visit(const SPnode * const spn) const;
			
			/*! \brief Get the result of the visit operation.*/
            bool getResult() const;

		private:

			IntervalImageToleranceCheck();
			
			const MappedFobj& fobj;
            const cxsc::real tolerance;
			
			mutable bool result;

    };
    // end of IntervalImageToleranceCheck class
	
	/*! \brief A type that visits \linkSPnodes SPnodes\endlink 
	to check if they meet a 'Reimann Difference' tolerance 
	requirement.

	The 'Reimann Difference' tolerance requirement
	is that the difference in the Reimann Integral taken to the top
	of the interval image of the box associated with the node against the 
	Reimann Integral taken to the bottom
	of the interval image of the box associated should be less than
	or equal to the tolerance given
	by the constructor argument \a tol. 
	* 
	ie the Reimann Difference can also be seen as the 'area' represented
	by the interval image of the box associated with the node, or the
	volume of the box associated with the node multiplied by the diameter
	of the interval enclosure of the box associated with the node.

	A node may not meet the requirement but not be splittable: in this
	case the result of the visit should be true (ie the node is 
	deemed to meet the requirement because it is as close as
	it can be to meeting it).

	The visit operation should return true if 
	either the the node visited meets the
	Reimann Difference tolerance requirement or the node does not meet
	the requirement but is not splittable.  Return false otherwise.
	*/
	class ReimannDiffToleranceCheck : public SPCheckVisitor {

        
        public:
	
			/*! \brief Constructor.
			
			\param f describes a function against which to check
			the interval image tolerance requirement.
			\param tol describes the tolerance to be used in the check.
			\pre \a tol >= cxsc::MinReal.*/
            ReimannDiffToleranceCheck(const MappedFobj& f, cxsc::real tol);
			
			~ReimannDiffToleranceCheck();

			/*! \brief The visit operation.
			 
			Checks the node pointed to by \a spn and stores the 
			result true if the  
			Reimann Difference tolerance requirement is met.
			
			\param spn a pointer to an SPnode to be visited.
			\return true if \a spn meets the Reimann Difference
			tolerance requirement or the node does not meet
			the requirement but is not splittable.  
			Return false otherwise.*/
            void visit(const SPnode * const spn) const;
			
			/*! \brief Get the result of the visit operation.*/
            bool getResult() const;

		private:

			ReimannDiffToleranceCheck();
			
			const MappedFobj& fobj;
            const cxsc::real tolerance;
			
			mutable bool result;

    };
    // end of ReimannDiffToleranceCheck class

} // end namespace subpavings

#endif

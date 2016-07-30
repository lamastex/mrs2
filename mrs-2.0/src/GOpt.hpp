/*
 * Copyright (C) 2005, 2006, 2007, 2008 Raazesh Sainudiin
 * Copyright (C) 2009 Jennifer Harlow
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
\brief Global Optimisation declarations.
Global optimisation using the C++ Toolbox for Verified Computing
See Hammer, R., Hocks, M., Kulisch, U., and Ratz, D (1995),
C++ Toolbox for Verified Computing, Springer, pp. 312-342

This routine uses a global variable.  This is not ideal but seems, on
balance, to be the most sensible way around the problems of having a
polymorphic implemention of AllGop, ie one that can take some function object
and use the appropriate method of that function object as the HTscalar_FctPtr
needed by AllGop.  Ideally, we would take an fObj as the parameter and use
its HTscalar_FctPtr method directly but to do that we have to bind 'this',
the function object, into the call to the method.  That is possible, but in
this case the signature of the HTscalar_FctPtr type, as defined for AllGop,
is HessType (*) (const HTvector&) and using the basic
stl functional binders, we can't bind to when one of the arguments is a
reference.  We could get around this by using the std::tr1::function binders,
which have wrappers for references, but this would mean that the process
would only compile and run where some implementation of tr1 is available and
this would be a severe limitiation.  If and when tr1 is 'officially'
implemented, we could clean all this up and get rid of the global.  We have
to use a pointer to an Fobj as the global because Fobj is an abstract class
and so cannot be declared on just as Fobj.  But in this only when the other
arguments
*/

#ifndef __GOPT_HPP__
#define __GOPT_HPP__



#include "Fobj.hpp"
#include "toolz.hpp"

/*! \brief This runs the global optimisation procedure AllGOp for global
    minimums and prints results.

 \param f a pointer to the Fobj we want to do global optimisation on
 \param search an initial search box
 \param t a tolerance for optimisation -- stops search when interval searched
        is below t in width
\param label the model or topology label
*/
void GOptMin(Fobj* f, ivector search, real t, const int label = 0);

/*! \brief This runs the global optimisation procedure AllGOp for global
    maximums and prints results.

\param f a pointer to the Fobj we want to do global optimisation on
\param search an initial search box
\param t a tolerance for optimisation -- stops search when interval searched
       is below t in width
\param label the model or topology label
*/
void GOptMax(Fobj* f, ivector search, real t, const int label = 0);

//! Function conforming to typedef HTscalar_FctPtr for global minimums
HessType funcHessMin(const HTvector& x);

//! Function conforming to typedef HTscalar_FctPtr for global maximums
HessType funcHessMax(const HTvector& x);

//! Print the results for mimimums
void printOutcomeMin(ivector& SearchInterval, real& Tolerance, imatrix& Opti,
                     intvector& Unique, int NumberOfOptis, interval& Minimum,
                     int Error);

//! Print the results for maximums
void printOutcomeMax(ivector& SearchInterval, real& Tolerance, imatrix& Opti,
                     intvector& Unique, int NumberOfOptis, interval& Maximum,
                     int Error);

// end of GOpt declarations
#endif

/*
 * Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
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

/*! \file toolz.hpp
\brief Declaration of various tools; functions and structs, for MRS.
*/

#ifndef __TOOLZ_HPP__
#define __TOOLZ_HPP__

#include <iostream>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp" // Include interval arithmetic package
#include "imath.hpp"    // Include interval standard functions
#include "rmath.hpp"    // Include real standard functions
#include "ivector.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_statistics.h>
#include <vector>
#include <functional>

#include <sys/types.h>
#include <unistd.h>

using namespace std;
using namespace cxsc;

//! Compute the sample mean of a double array using a recurrence.
double mean (const size_t ssize, const double *x);

//!Compute the Variance.
double var (const size_t ssize, double *x);

// !Compute the sample mean and variance.
long mean_var (const vector < real > &x, real & mean, real & var);

//! Compute MSE using gsl_stats_variance_with_fixed_mean.
double MSE (double exact, int ssize, double *x);

//! Find Minimum of two real types.
real Min(real x, real y);

//! Draw a real uniformly at random from the interval X in R.
real DrawUnifInterval(gsl_rng* rgsl, const interval& X);

//! Draw a vector uniformly at pseudo-random from unit box [0,1)^n_dimensions.
rvector DrawUnifUnitBox (gsl_rng* rgsl, const int n_dimensions);

//! Draw a vector uniformly at pseudo-random from a Box
rvector DrawUnifBox(gsl_rng* rgsl, const ivector& X);

//! Draw a vector uniformly at quasi-random from a Box
rvector DrawQZUnifBox(gsl_qrng* qrgsl, const ivector& X);

//! Draw a vector uniformly at quasi-random from a Box: trans-dimensional case.
rvector DrawQZUnifBoxV(double* v, const ivector& X);

//! Return the first dimension with maximal diameter.
int MaxDiamComp(ivector& iv);

//! Blow up a box by eps pivoted at FromZero and return it.
ivector BlowUpFromZero(ivector iv, real FromZero, real eps);

//! Stable Summation Routine of Kahan's
template<typename T>
struct kahan_sum
{
  T s,c,y,t;
  kahan_sum() : s(0.),c(0.),y(0.),t(0.){}
  T & operator()(  T & v,const T & i )
  {
    /* c is zero or close to it*/
    y=i-c;

    /* If s is big and y small, then low-order digits of y are lost during
    the summation into t*/
    t=s+y;

    /* (t - s) recovers the high-order part of y; subtracting y
      recovers -(low part of y)*/
    c=(t-s)-y;

    /* Algebriacally, c should always be zero. Beware eagerly optimising
      compilers! */
    s=t;

    /* Next time around, the lost low part will be added to y in a fresh
      attempt. */
    return s;
  }
};

//! Return the maximal diameter of box x.
//double MaxDiam (ivector& x, int& c);
double MaxDiam (ivector x, int& c);

/*! \brief Compute the intersection interval r of interval a and interval b.

\return 0 if the intersection is empty, else return 1.
*/
int Intersection (interval & r, const interval & a, const interval & b);

/*! \brief Compute the intersection box r of box a and box b.

\return 0 if the intersection is empty, else return dimension of r.
*/
int Intersection (ivector & r, const ivector & a, const ivector & b);

//! Bisect box x normal to direction "split" and  return the lower half.
ivector Lower (const ivector & x, int split);

//! Bisect box x normal to direction "split" and return the upper half.
ivector Upper (const ivector & x, int split);

//! Bisect box x normal to direction "split" and set box y to the lower half.
void Lower (const ivector & x, ivector & y, int split);

//! Bisect box x normal to direction "split" and set box y to the upper half.
void Upper (const ivector & x, ivector & y, int split);

/*! \brief Return the volume of box x.

  unclear what behaviour of vol() is if the interval vector x is uninitialised.
*/
double Volume(const ivector &x);

//src_trunk_0701
/*! \brief Return the volume of box x as a real.

  unclear what behaviour of vol() is if the interval vector x is uninitialised.
*/
real realVolume(const ivector &x);

/*! \brief Function to return log of Catalan number of k

In a histogram, k is the number of leaves-1, ie the number of splits.
Catalan(k) = 1/(k+1)Binomial(2k choose k) = (2k)!/((k+1)!k!)
ln(Catalan(k)) is 0 for k-0, k=1,
sum(from i=1 to i=k-2) of (ln(2k-i) - ln(k-i)) for integer k > 1.

\param k the integer to calculate the Catalan number for.
\return the natural log of the Catalan number for k.
*/
double lCk(const int k);

//--src_trunk_0701
/*! \brief String representation of rvector, formatted with brackets etc.
 * 
 * \param rv the rvector to describe in the output.
 * \return string representation of \a rv.*/
std::string toString(const cxsc::rvector& rv);

/*! \brief String representation of ivector, formatted with brackets etc.
 * 
 * \param iv the ivector to describe in the output.
 * \return string representation of \a iv.*/
std::string toString(const cxsc::ivector& iv);

/*! \brief String representation of interval, formatted with brackets etc.
 * 
 * \param ival the interval to describe in the output.
 * \return string representation of \a ival.*/
std::string toString(const cxsc::interval& ival);

/*! \brief Output formatted with brackets etc, suitable for human reading of rvectors.
 * 
 * \param out the stream to send the output to.
 * \param rv the rvector to describe in the output.
 * \return the stream with the new output.*/
std::ostream& prettyPrint(std::ostream& out, const cxsc::rvector& rv);

/*! \brief Output formatted with brackets etc, suitable for human reading of ivectors.
 * 
 * \param out the stream to send the output to.
 * \param iv the ivector to describe in the output.
 * \return the stream with the new output.*/
std::ostream& prettyPrint(std::ostream& out, const cxsc::ivector& iv);

/*! \brief Output formatted with brackets etc, suitable for human reading of intervals.
 * 
 * \param out the stream to send the output to.
 * \param ival the interval to describe in the output.
 * \return the stream with the new output.*/
std::ostream& prettyPrint(std::ostream& out, const cxsc::interval& ival);
//--src_trunk_0701

#endif

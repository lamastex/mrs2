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

/*! \file
\brief Definition of various tools; functions and structs, for MRS.
*/

#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <algorithm>// to use stl::algorithms
#include <stdexcept> 

#include "toolz.hpp"

using namespace std;

//functions

/*! Compute the arithmetic mean of a double array using the recurrence relation
mean_(n) = mean(n-1) + (x[n] - mean(n-1))/(n+1).
*/
double mean (const size_t ssize, const double *x)
{
  return gsl_stats_mean (x, 1, ssize);
 /* long double mean = 0;
  size_t i;

  for (i = 0; i < ssize; i++)
  {
    mean += (x[i] - mean) / (i + 1);
  }

  return mean;*/
}

/*! the non-parametric plug-in estimate of variance (1/N) sum (x_i - mean)^2
*/
double var (const size_t ssize, double *x)
{
  return gsl_stats_variance (x, 1, ssize)*((double) (ssize-1) / (double) ssize);
}

/*! Compute the arithmetic mean of a double array -- curruntly ugly/unstable

\todo using the recurrence relation
  in kahan_sum and the sample variance using
  a similar recurrence OR cxsc DotPrecision Accumulators.
*/
long mean_var (const vector < real > &x, real & mean, real & var)
{
  vector < real >::const_iterator it = x.begin ();
  real sum (0.0), sumsq (0.0);
  long count (0);
  for (; it != x.end (); ++it)
  {
    sum += *it;
    sumsq += *it * *it;
    count++;
  }
  mean = sum / count;
  var = sumsq / count - mean * mean;
  return count;
  //sum = std::accumulate(x.begin(),x.end(),0.,kahan_sum<real>());
}

// mean squared error via GSL
double MSE (double exact, int ssize, double *x)
{
  return gsl_stats_variance_with_fixed_mean (x, 1, ssize, exact);
}

// minimum of two reals
real Min (real x, real y)
{
  return (x < y) ? x : y;
}

// Draw a real uniformly at random from the interval X in R
real DrawUnifInterval (gsl_rng * rgsl, const interval & X)
{
  double r = gsl_rng_uniform (rgsl);
  return (Inf (X) + (r * diam (X)));
}

/*! return rvector representing point drawn uniformly pseudo-randomly from
  unit box with n_dimensions dimensions.
*/
rvector DrawUnifUnitBox (gsl_rng* rgsl, const int n_dimensions)
{
  rvector rand_vector(1, n_dimensions);
  for(int i=1; i<=n_dimensions; i++)
  {
    rand_vector[i] = gsl_rng_uniform(rgsl);
  }
  return rand_vector;
}

//Draw a real uniformly at random from the box X in R^n
rvector DrawUnifBox (gsl_rng * rgsl, const ivector & X)
{
  int i, a = Lb (X), z = Ub (X);
  rvector UX (a, z);
  for (i = a; i <= z; ++i)
  {
    UX[i] = DrawUnifInterval (rgsl, X[i]);
  }
  return UX;
}

// Draw a real uniformly at quasi-random from the box X in R^n
rvector DrawQZUnifBox (gsl_qrng * qrgsl, const ivector & X)
{
  int i, a = Lb (X), z = Ub (X);
  vector<double> v(z - a + 1);

  gsl_qrng_get (qrgsl, (&v[0]));
  rvector UX (a, z);
  for (i = a; i <= z; ++i)
  {
    UX[i] = Inf (X[i]) + (diam (X[i]) * v[i - a]);
  }
  return UX;
}

/*! Draw a real uniformly at quasi-random from the box X in R^n given a
quasi-random vector in hypercube; the quasi-random vector v may have
dimensionality greater than that of X, in which case the extra elements at
the end of v are not used.
*/
rvector DrawQZUnifBoxV (double *v, const ivector & X)
{
  int i, a = Lb (X), z = Ub (X);
  rvector UX (a, z);
  for (i = a; i <= z; ++i)
  {
    UX[i] = Inf (X[i]) + (diam (X[i]) * v[i - a]);
  }
  return UX;
}

// return the first dimension with maximal diameter
int MaxDiamComp (ivector & iv)
{
  int mc, k;
  rvector d (Lb (iv), Ub (iv));

  d = diam (iv);
  mc = Lb (iv);
  for (k = Lb (iv) + 1; k <= Ub (iv); ++k)
  {
    if (d[k] > d[mc])
      mc = k;
  }
  return mc;
}

// blow the box by eps
ivector BlowUpFromZero (ivector iv, real FromZero, real eps)
{
  int i, a = Lb (iv), z = Ub (iv);
  interval x;
  interval EPS = _interval (-eps, eps);
  ivector ov (a, z);
  for (i = a; i <= z; ++i)
  {
    x = iv[i] + EPS;
    ov[i] =
      _interval (pred (max (FromZero, Inf (x))),
      succ (max (FromZero, Sup (x))));
  }
  return ov;
}

//double MaxDiam (ivector& x, int& c)
double MaxDiam (ivector x, int& c)
{
  int i, a=Lb(x), z=Ub(x);
  real Diam = diam(x[a]); c = a;
  for (i=a+1;i<=z;i++)
    if (diam(x[i])>Diam) { Diam = diam(x[i]); c = i; }
    return _double(Diam);
}

// the intersection interval r of interval a and interval b.
// return 0 if the intersection is empty, else return 1
int Intersection (interval & r, const interval & a, const interval & b)
{
  // disjoint?
  if( (Inf(a) > Sup(b)) || (Inf(b) > Sup(a)) )
    return 0;
  else {r = a & b; return 1;}
}

// the intersection box r of box a and box b.
// return 0 if the intersection is empty, else return dimension of r.
int Intersection (ivector & r, const ivector & a, const ivector & b)
{
  int i, il = Lb(a), iu = Ub(a);
  int intersect = iu-il+2;
  if( (iu != Ub(b)) || (il != Lb(b)) )
  {
    cerr << "\n unequal index sets in Intersection\n";
    exit(0);
  }
  Resize(r, il, iu);
  for (i = il; i <= iu; i++)
  {
    // disjoint?
    if( (Inf(a[i]) > Sup(b[i])) || (Inf(b[i]) > Sup(a[i])) )
      return 0;
    else {r[i] = a[i] & b[i]; intersect--;}
  }
  return intersect;
}

// Bisect box x normal to direction "split" and return the lower half.
ivector Lower (const ivector & x, int split)
{
  ivector t = x;
  SetSup( t[split], mid(x[split]) );
  return t;
}

// Bisect box x normal to direction "split" and return the upper half.
ivector Upper (const ivector & x, int split)
{
  ivector t = x;
  SetInf( t[split], mid(x[split]) );
  return t;
}

// Bisect box x normal to direction "split" and set y to the lower half.
void Lower (const ivector & x, ivector & y, int split)
{
  Resize(y, Lb(x), Ub(x));
  y = x;
  SetSup( y[split], mid(x[split]) );
}

// Bisect box x normal to direction "split" and set box y to the upper half.
void Upper (const ivector & x, ivector & y, int split)
{
  Resize(y, Lb(x), Ub(x));
  y = x;
  SetInf( y[split], mid(x[split]) );
}

// to find the volume of an interval vector x
double Volume(const ivector &x)
{
  return _double( realVolume(x) );
}

// to find the volume of an interval vector x
real realVolume(const ivector &x)
{
  int low_index = Lb(x);
  int upp_index = Ub(x);
  
  if (upp_index < low_index) throw std::logic_error("Ub < Lb");

  real accVol = 1.0;

  for (int i= low_index; i<=upp_index; ++i)
  {

    accVol *= diam(x[i]);

  }                 //

  return accVol;
}

//Function to return log of Catalan number of k

//In a histogram, k is the number of leaves-1, ie the number of splits.
//Catalan(k) = 1/(k+1)Binomial(2k choose k) = (2k)!/((k+1)!k!)
//ln(Catalan(k)) is     0 for k=0, k=1,
//                      sum(from i=1 to i=k-2)of(ln(2k-i) - ln(k-i)), k>1
double lCk(const int k)
{
    double retValue = 0.0;
    if (k > 1) {
        for (int i = 0; i < k-1; ++i) {
            retValue += (log(2*k-i) - log(k-i));
        }
    }
    return retValue; // return 0 for k=0, k=1
}

// Output formatted with brackets etc, suitable for human reading of rvectors
std::ostream& prettyPrint(std::ostream& out, const cxsc::rvector& rv)
{
	out << "(";
	for(int i=Lb(rv); i<=Ub(rv); ++i)
	{
	out << " " << rv[i];
	}
	out << " )"; 

	return out;
}

// Output formatted with brackets etc, suitable for human reading of ivectors
std::ostream& prettyPrint(std::ostream& out, const cxsc::ivector& iv)
{
	for(int i=Lb(iv); i<=Ub(iv); ++i)
	{
	out << " ";
	prettyPrint(out, iv[i]);
	}


	return out;
}

// Output formatted with brackets etc, suitable for human reading of intervals
std::ostream& prettyPrint(std::ostream& out, const cxsc::interval& ival)
{
	out << "[ " << Inf(ival) << ", " << Sup(ival) << " ]";

	return out;
}

// Plain output (alternative to silly cxsc format),tab delimited, suitable for output to txt file
std::ostream& plainPrint(std::ostream& out, const cxsc::rvector& rv)
{
	int i = 0;
	for(i=Lb(rv); i<Ub(rv); ++i)
	{
		out << rv[i] << "\t";
	}
	out << rv[i];  // final one with no tab at the end 

	return out;
}

std::string toString(const cxsc::rvector& rv)
{
	std::ostringstream oss;
	prettyPrint(oss, rv);
	return oss.str();
}

std::string toString(const cxsc::ivector& iv)
{
	std::ostringstream oss;
	prettyPrint(oss, iv);
	return oss.str();
}

std::string toString(const cxsc::interval& ival)
{
	std::ostringstream oss;
	prettyPrint(oss, ival);
	return oss.str();
}

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

/*! \file      FLevy2D.cpp
\brief Implementation for example function class FLevy2D 
(Levy function, 2 dimensions).
*/

#include "FLevy2D.hpp"

/*
// Older tempalized version as a function object -- more convenient for testing 
template <class Arg, class Result>
class FLevy: public std::unary_function<Arg, Result>{
// Std normal located at -5; 1 dim
public:
    Result operator() (const Arg &X) const {
      int a=Lb(X);
      int i, z=Ub(X);
      Result isum, jsum, hh;
      isum = 0.0; jsum = 0.0;
      for (i = 1; i <= 5; i++) {
        isum = isum + double(i)*cos(double(i-1)*X[a] + double(i));
        jsum = jsum + double(i)*cos(double(i+1)*X[z] + double(i));
      }
      hh = isum * jsum +
           sqr(X[a] + 1.42513) +    // Avoid real con-
           sqr(X[z] + 0.80032);      // version error
        //return exp(-hh/10.0);
        return exp(-hh/TEMPERATURE); // TEMPERATURE = 1, 4, 40, 400, 4000
    }
};
*/

/*! We pass the global maximum as a parameter to ensure the containment of the
density form of the original Levy target as an "energy" function inside the
number screen.  We get the density from energy by exponentiating its negative.  
This allows ease of likelihood inference and can be circumvented in several 
ways.  For eg we can find the global min first for the Levy target energy and 
then use it as the global max parameter GlbMx to the density form (actually the
shape without the normalizing constant).  Here we have already computed the 
Global max parameter GlbMx using C-XSC Toolbox.
*/
FLevy2D::FLevy2D (real T, real GlbMx, real C1, real C2, 
                  real DomainLimit, bool LogPi)
:
Temperature (T), GlobalMax (GlbMx), Center1 (C1), Center2 (C2)
{
  setUsingLogDensity (LogPi);

  PriorType = 0;// Uniform PriorType is an inherited member from Fobj
  // set up the domain list
  ivector domain (1, 2);
  LabBox  Ldomain;
  for (int i = 1; i <= 2; i++)
  {
    domain[i] = interval (-DomainLimit, DomainLimit);
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);

}

interval FLevy2D::operator () (const LabBox & X) const
{
  n_interval_calls++;

  ivector Box = X.Box;
  int a = Lb (Box), z = Ub (Box);
  interval isum, jsum, hh;
  isum = 0.0;
  jsum = 0.0;
  for (int i = 1; i <= 5; i++)
  {
    isum = isum + double(i) * cos (double(i - 1) * Box[a] + double(i));
    jsum = jsum + double(i) * cos (double(i + 1) * Box[z] + double(i));
  }

                    // Avoid real con-
  hh = isum * jsum + sqr (Box[a] + Center1) +
                    // version error
    sqr (Box[z] + Center2);
    hh += GlobalMax;  
  // TEMPERATURE = 1, 4, 40, 400, 4000
  interval result = exp (-hh / Temperature);

  return (UsingLogDensity) ? ln (result) : result;
}

real FLevy2D::operator () (const LabPnt & X) const
{
  n_real_calls++;
  rvector Pnt = X.Pnt;
  int a = Lb (Pnt), z = Ub (Pnt);
  real isum, jsum, hh;
  isum = 0.0;
  jsum = 0.0;

  for (int i = 1; i <= 5; i++)
  {
    isum = isum + double (i) * cos (double (i - 1) * Pnt[a] + double (i));
    jsum = jsum + double (i) * cos (double (i + 1) * Pnt[z] + double (i));
  }
  // Avoid real conversion error
  hh = isum * jsum + sqr (Pnt[a] + Center1) +
       sqr (Pnt[z] + Center2);
  hh += GlobalMax;

  // TEMPERATURE = 1, 4, 40, 400, 4000
  real result = exp (-hh / Temperature);
  //result /= sqrt(Temperature);
  //result /= Temperature;
  assert (result < 1.0);
  return (UsingLogDensity) ? ln (result) : result;
}

// HessType operator()
// label is not used but needed to match Fobj signature
HessType FLevy2D::operator()(const HTvector& x, const int label) const
{

  HessType isum(2), jsum(2), hh(2);
  isum = 0.0;
  jsum = 0.0;

  for (int i = 1; i <= 5; i++)
  {
    isum = isum + double (i) * cos (double (i - 1) * x[1] + double (i));
    jsum = jsum + double (i) * cos (double (i + 1) * x[2] + double (i));
  }

  hh = isum * jsum + sqr (x[1] + Center1 + sqr (x[2] + Center2));
  hh = hh + GlobalMax;

  // TEMPERATURE = 1, 4, 40, 400, 4000
  hh = exp (-hh / Temperature);

  return (UsingLogDensity) ? ln (hh) : hh;
}

FLevy2D_Lkl_Tfrom1data::FLevy2D_Lkl_Tfrom1data (ivector & data, 
                                                real GlbMx, real C1, real C2, 
                                                interval DomainLimit, 
                                                bool LogPi)
:
Data (data), GlobalMax (GlbMx), Center1 (C1), Center2 (C2)
{
  setUsingLogDensity (LogPi);
  //PriorType = 0;// Uniform PriorType is an inherited member from Fobj
  //domain for the temperature parameter in the likelihood
  ivector domain (1, 1);
  LabBox Ldomain;
  for (int i = 1; i <= 1; i++)
  {
    domain[i] = DomainLimit;
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);

}

interval FLevy2D_Lkl_Tfrom1data::operator () (const interval & X) const
{
  int a = Lb (Data), z = Ub (Data);
  interval isum, jsum, hh;
  isum = 0.0;
  jsum = 0.0;
  for (int i = 1; i <= 5; i++)
  {
    isum =
      isum + double (i) * cos (double (i - 1) * Data[a] + double (i));
    jsum =
      jsum + double (i) * cos (double (i + 1) * Data[z] + double (i));
  }
                    // Avoid real con-
  hh = isum * jsum + sqr (Data[a] + Center1) +
                    // version error
    sqr (Data[z] + Center2);
  hh += GlobalMax;

   // TEMPERATURE = 1, 4, 40, 400, 4000
  interval result = exp (-hh / X);

  result /= sqr (X);
  //result /= X;
  return (UsingLogDensity) ? ln (result) : result;
}

real FLevy2D_Lkl_Tfrom1data::operator () (const real & X) const
{
  int a = Lb (Data), z = Ub (Data);
  // here we simply take the mid-point of the data to get thin data
  rvector realData = mid (Data);
  real isum, jsum, hh;
  isum = 0.0;
  jsum = 0.0;
  for (int i = 1; i <= 5; i++)
  {
    isum =
      isum + double (i) * cos (double (i - 1) * realData[a] +
      double (i));
    jsum =
      jsum + double (i) * cos (double (i + 1) * realData[z] +
      double (i));
  }

                    // Avoid real con-
  hh = isum * jsum + sqr (realData[a] + Center1) +
                    // version error
    sqr (realData[z] + Center2);
  hh += GlobalMax;

  // TEMPERATURE = 1, 4, 40, 400, 4000
  real result = exp (-hh / X);
  result /= sqr (X);
  return (UsingLogDensity) ? ln (result) : result;
}

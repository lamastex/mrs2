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

/*! \file      FRosenbrock.cpp
    \brief Implementation for example function class FRosenbrock (Rosenbrock function).
*/

#include "FRosenbrock.hpp"

// Constructor
// if LogPi is true, the log of this will be returned
FRosenbrock::FRosenbrock (int n_dimensions, real T, real H, real DomainLimit, 
                          bool LogPi):num_dim (n_dimensions), Tinverse (T),
Height (H)
{
  n_interval_calls = 0;
  n_real_calls = 0;

  setUsingLogDensity (LogPi);
  PriorType = 0;// Uniform PriorType is an inherited member from Fobj

  // set up the domain
  ivector domain (1, n_dimensions);
  LabBox  Ldomain;
  for (int i = 1; i <= n_dimensions; i++)
  {
    domain[i] = interval (-DomainLimit, DomainLimit);
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  //LabDomainList is a data member of the base Fobj class
  LabDomainList.push_back (Ldomain);
}

// operator() returning interval
interval FRosenbrock::operator () (const LabBox & X) const
{
  n_interval_calls++;

  int i, a = Lb (X.Box), z = Ub (X.Box);

  interval result = _interval (0.0);

  if (num_dim == 1)
  {
    result = (Height * sqr (1.0 - sqr (X.Box[a])) + sqr (X.Box[a] - 1.0));
  }
  else
  {
    for (i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (X.Box[i] - sqr (X.Box[i - 1])) +
        sqr (X.Box[i - 1] - 1.0));
    }
  }

  result = exp (-(Tinverse * result));
  //if(Sup(result) < 0.0) Sup(result) = 0.0;
  //if(Inf(result) < 0.0) Inf(result) = 0.0;
  return (UsingLogDensity) ? ln (result) : result;
}

// operator () returning real
real FRosenbrock::operator () (const LabPnt & X) const
{
  n_real_calls++;

  int i, a = Lb (X.Pnt), z = Ub (X.Pnt);

  real result = 0.0;

  if (num_dim == 1)
  {
    result = (Height * sqr (1.0 - sqr (X.Pnt[a])) + sqr (X.Pnt[a] - 1.0));
  }
  else
  {
    for (i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (X.Pnt[i] - sqr (X.Pnt[i - 1])) +
        sqr (X.Pnt[i - 1] - 1.0));
    }
  }

  result = exp (-(Tinverse * result));
  //if(result < 0.0) result = 0.0;
  //if(result < 0.0) result = 0.0;
  return (UsingLogDensity) ? ln (result) : result;

}

// operator() returning HessType
// label is not used but needed to match Fobj signature
HessType FRosenbrock::operator()(const HTvector& x, const int label) const
{

  HessType hh(num_dim);
  hh = 0.0;

  if (num_dim == 1)
  {
    hh = (Height * sqr (1.0 - sqr (x[1])) + sqr (x[1] - 1.0));
  }
  else
  {
    for (int i = 2; i <= num_dim; i++)
    {
      hh = hh + (Height * sqr (x[i] - sqr (x[i-1])) +
        sqr (x[i - 1] - 1.0));
    }
  }

  hh = exp (-(Tinverse * hh));

  return (UsingLogDensity) ? ln (hh) : hh;
}

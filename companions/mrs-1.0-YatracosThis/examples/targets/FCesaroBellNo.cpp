/* 
 * Copyright (C) 2005, 2006, 2007, 2008 Raazesh Sainudiin and Thomas York
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

#include "FCesaroBellNo.hpp"

// #define DOMAINLIMIT 10.0 // for now domain size is [-DOMAINLIMIT, DOMAINLIMIT]^d
// implementation of FCesaroBellNo class

// constructor for absolute shape of the integrand in Cesaro's integral for n-th Bell number
FCesaroBellNo::FCesaroBellNo (real Nn, interval DomainLimit, bool LogPi):
n(Nn)
{
  n_interval_calls = 0;
  n_real_calls = 0;
  setUsingLogDensity (LogPi);

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

interval
FCesaroBellNo::operator () (const LabBox & X)
     const
     {
       n_interval_calls++;
       interval x = X.Box[1];
       interval cx=cos(x);
       interval sx=sin(x);
       interval ecx=exp(cos(x));
       interval result=abs( exp(ecx * cos(sx)) * sin(ecx * sin(sx)) * sin(n * x) );
       return (UsingLogDensity) ? ln (result) : result;
     }

real
FCesaroBellNo::operator () (const LabPnt & X)
     const
     {
       n_real_calls++;
       real x = X.Pnt[1];
       real cx=cos(x);
       real sx=sin(x);
       real ecx=exp(cos(x));
       real result=abs( exp(ecx * cos(sx) ) * sin(ecx * sin(sx)) * sin(n * x) );
       return (UsingLogDensity) ? ln (result) : result;
     }

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

// example function object class to use with MRSampler class

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"
#include "imath.hpp"
#include "rmath.hpp"
#include "ivector.hpp"
#include <functional>
#include <assert.h>

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FShiryaev1D.hpp"

// implementation of FShiryaev1D class -- the stretched oscilating exponential

FShiryaev1D::FShiryaev1D (real A, real E, 
                          real DomainLimit, bool LogPi, int Prior)
:
aa (A),
cc(E)
{
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  if(PriorType!=0)
  {
    cerr << "Only uniform prior defined here... Reset to Uniform Prior\n";
    PriorType=0;
  }
  AA = _interval (aa);
  CC = _interval(cc);
  ivector domain (1, 1);
  LabBox Ldomain;
  domain[1] = interval (0.0000000000001, DomainLimit);
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);
}

interval
FShiryaev1D::operator () (const LabBox & X)
const
{
  n_interval_calls++;
  ivector Box = X.Box;
                    //, z=Ub(Box);
  size_t a = Lb (Box);
  //interval Nc = pow(AA, 1.0/BB) / GAMMA1PLUS1OVERBB;
  interval aXtoL = AA * pow (Box[a], BB);
  //interval result = 
  // Nc * exp(- aXtoL) * (1.0 + (CC * sin( aXtoL * tan( BB * PI) )));
  interval result =
    exp (-aXtoL) * (1.0 + (CC * sin (aXtoL * tan (BB * PI))));
  return (UsingLogDensity) ? ln (result) : result;
}

real
FShiryaev1D::operator () (const LabPnt & X) const
{
  n_real_calls++;
  rvector Pnt = X.Pnt;
                    //, z=Ub(Pnt);
  size_t a = Lb (Pnt);
  //real Nc = pow(aa, 1.0/bb) / Gamma1plus1overbb;
  real aXtoL = aa * pow (Pnt[a], bb);
  //real result = Nc * exp(- aXtoL) * 
  //              (1.0 + (cc * sin( aXtoL * tan( bb * Pi) )));
  real result =
    exp (-aXtoL) * (1.0 + (cc * sin (aXtoL * tan (bb * Pi))));
  //assert( result < 1.0 );
  return (UsingLogDensity) ? ln (result) : result;
}

FShiryaev1D_Lkl_aa_fromData::
FShiryaev1D_Lkl_aa_fromData (RSSample & d, real E, interval DomainInterval, 
                             bool LogPi, int Prior)
:
Data (d), 
cc(E)
{
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  if(PriorType!=0)
  {
    cerr << "Only uniform prior defined here... Reset to Uniform Prior\n";
    PriorType=0;
  }
  CC = _interval(cc);
                    //domain for the alpha parameter in the likelihood
  ivector domain (1, 1);
  LabBox Ldomain;
  domain[1] = DomainInterval;
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);
}

interval
FShiryaev1D_Lkl_aa_fromData::operator () (const LabBox & X) const
{
  n_interval_calls++;
                    //, z=Ub(Box);
  size_t a = Lb (X.Box);
  // the logarithm of data-independent term in the product likelihood
  interval lnNc = _real(Data.Samples.size()) * 
                  ln( pow (X.Box[a], 1.0 / BB) / GAMMA1PLUS1OVERBB );
  
  interval result = lnNc;
  vector<LabPnt>::const_iterator it = Data.Samples.begin();
  for(size_t i0=Lb(it->Pnt); it!=Data.Samples.end(); it++)
  {
    interval aXtoL = X.Box[a] * pow (_interval(it->Pnt[i0]), BB);
    //! \todo use dotprec accumulator here instead
    result += ( (-aXtoL) + ln(1.0 + (CC * sin (aXtoL * tan (BB * PI)))) );
  }
  return (UsingLogDensity) ? (result) : exp(result);
}

real FShiryaev1D_Lkl_aa_fromData::operator () (const LabPnt & X) const
{
  n_real_calls++;
  size_t a = Lb (X.Pnt);
  // the logarithm of data-independent term in the product likelihood
  real lnNc = _real(Data.Samples.size()) * 
              ln( pow (X.Pnt[a], 1.0 / bb) / Gamma1plus1overbb );
  real result = lnNc;
  vector<LabPnt>::const_iterator it = Data.Samples.begin();
  for(size_t i0=Lb(it->Pnt); it!=Data.Samples.end(); it++)
  {
    real aXtoL = X.Pnt[a] * pow (it->Pnt[i0], bb);
    //! \todo use dotprec accumulator here instead
    result += ( (-aXtoL) + ln(1.0 + (cc * sin (aXtoL * tan (bb * Pi)))) );
  }
  return (UsingLogDensity) ? (result) : exp(result);
}

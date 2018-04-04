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

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FWitchesHat.hpp"

// for now domain size is [-DOMAINLIMIT, DOMAINLIMIT]^d
// implementation of FWitchesHat class

// constructor for "witches hat": mixture of uniform and conical peak
// if LogPi is true, the log of this will be returned
FWitchesHat::FWitchesHat (int n_dimensions, real PMean, real PRadius, 
                          real PWeight, real UWeight, real DomainLimit, 
                          bool LogPi, int Prior)
:
PeakMean (PMean), 
PeakRadius (PRadius), 
PeakWeight (PWeight),
UniformWeight (UWeight)
{
  //  witches hat peak is at (PeakMean, PeakMean, ... , PeakMean)
  n_interval_calls = 0;
  n_real_calls = 0;
  setUsingLogDensity (LogPi);
  h[0]=0.0; h[1]=1.0; h[2]=0.95493; h[3]=0.95493; h[4]=1.01321; h[5]=1.13986;
  h[6]=1.35456; h[7]=1.69321; h[8]=2.21745; h[9]=3.03167; h[10]=4.31345;
  //{0.0,1.0,0.95493,0.95493,1.01321,1.13986,1.35456,
  //1.69321,2.21745,3.03167,4.31345};

  ivector domain (1, n_dimensions);
  LabBox Ldomain;
  DomainVolume = 1.0;
  Rtod = PeakRadius;
  PriorType = Prior;
  if(PriorType!=0)
  {
    cerr << "Only uniform prior defined here... Reset to Uniform Prior\n";
    PriorType=0;
  }
  for (int i = 1; i <= n_dimensions; i++)
  {
    domain[i] = interval (-DomainLimit, DomainLimit);
    DomainVolume *= 2.0 * DomainLimit;
    Rtod *= PeakRadius;
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);

}

// vector<LabBox> FWitchesHat::get_domain(){ return LabDomainList; }

interval
FWitchesHat::operator () (const LabBox & X)
const
{
  n_interval_calls++;
  int i, a = Lb (X.Box), z = Ub (X.Box);
  int ndim = 1 + z - a;
  if (ndim > 10)
  {
    cout << "> 10 dimensions. exiting. " << endl;
    exit (EXIT_FAILURE);
  }
  //   real DomainVolume = (2.0*DOMAINLIMIT);
  interval r = sqr (X.Box[a] - PeakMean);
  for (i = a + 1; i <= z; i++)
  {
    r += sqr (X.Box[i] - PeakMean);
    //  DomainVolume *= (2.0*DOMAINLIMIT);
  }
  r = sqrt (r);
  interval result = h[ndim] * (PeakRadius - r) / Rtod;

  if (Sup (result) < 0.0)
    Sup (result) = 0.0;
  if (Inf (result) < 0.0)
    Inf (result) = 0.0;
                    // mixture of uniform and wh peak
  result = UniformWeight / DomainVolume + PeakWeight * result;
  return (UsingLogDensity) ? ln (result) : result;
}

real
FWitchesHat::operator () (const LabPnt & X)
const
{
  n_real_calls++;
  int i, a = Lb (X.Pnt), z = Ub (X.Pnt);
  int ndim = 1 + z - a;
  if (ndim > 10)
  {
    cout << "> 10 dimensions. exiting. " << endl;
    exit (EXIT_FAILURE);
  }
  //  real DomainVolume = (2.0*DOMAINLIMIT);
  real r = sqr (X.Pnt[a] - PeakMean);
  for (i = a + 1; i <= z; i++)
  {
    r += sqr (X.Pnt[i] - PeakMean);
    //  DomainVolume *= (2.0*DOMAINLIMIT);
  }
  r = sqrt (r);
  real result = h[ndim] * (PeakRadius - r) / Rtod;

  if (result < 0.0)
    result = 0.0;
  // if(result > 0.0) 
  //   cout << "Pnt: " << X.Pnt << "  r: "  << r << "  peakhgt: " << result 
  //        << "  Rtod: " << Rtod << endl;
  
  // mixture of uniform and wh peak
  result = UniformWeight / DomainVolume + PeakWeight * result;
  return (UsingLogDensity) ? ln (result) : result;
}

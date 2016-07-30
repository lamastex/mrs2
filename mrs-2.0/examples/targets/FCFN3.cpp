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
/*! \file FCFN3.cpp
  \brief Trans-dimensional three-taxa Cavender-Farris-Neyman phylogenetic
  model example function object class to use with MRSampler class.
*/

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
#include "FCFN3.hpp"

// in case range enclosure is thin interval, lower Inf
#define FATTEN_THIN_INTERVAL_RE true
/*! constructor for the likelihood of the common branchlength under
  CFN model for a rooted and clocked 3 taxa tree
*/
// if LogPi is true, the log of this will be returned
FCFN3Star::FCFN3Star(int Countid, int Countnid, interval Domain, 
                     bool LogPi, int Prior)
:
Cid (Countid), Cnid (Countnid)
{
  n_interval_calls = 0;
  n_real_calls = 0;
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  //TotSites = real(Cid+Cnid);
  //f0 = real(Cid); // /TotSites;
  n_dimensions=1;
  ivector domain (1, n_dimensions);
  LabBox Ldomain;
  for (int i = 1; i <= n_dimensions; i++)
  {
    domain[i] = Domain;
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);
}

interval
FCFN3Star::operator () (const LabBox & X)
const
{
  n_interval_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Box) != 1 || Lb(X.Box) != 1)
  {
    cerr << "dimensions !=1 OR start-index != 1... exiting. " 
         << endl; exit (EXIT_FAILURE);
  }
  #endif
  interval em4t = (exp (-4.0 * X.Box[1])) ;
  interval result = Cid*ln( ((1.0 + (3.0*em4t)) / 8.0)) + 
                    Cnid*ln( ((1.0 - em4t) / 8.0)) ;
  // result = (TotSites*result);
  /*
  if (Sup (result) < 0.0) Sup (result) = 0.0;
  if (Inf (result) < 0.0) Inf (result) = 0.0;
  */
  return (UsingLogDensity) ? (result) : exp(result);
}

real
FCFN3Star::operator () (const LabPnt & X)
const
{
  n_real_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Pnt) != 1 || Lb(X.Pnt) != 1)
  {
    cout << "dimensions !=1 OR start-index != 1... exiting. " 
         << endl; exit (EXIT_FAILURE);
  }
  #endif
  real em4t = (exp (-4.0 * X.Pnt[1])) ;
  real result = Cid*ln( ((1.0 + (3.0*em4t)) / 8.0)) + 
                (Cnid)*ln( ((1.0 - em4t) / 8.0)) ;
  // result = (TotSites*result);
  return (UsingLogDensity) ? (result) : exp(result);
}

/*! constructor for the likelihood of the branchlengths under CFN model for an 
  unrooted 3 taxa tree
*/
// if LogPi is true, the log of this will be returned
FCFN3UnRooted::FCFN3UnRooted(int Countxxx, int Countxxy, int Countyxx, 
                             int Countxyx, 
                             interval Domain, bool LogPi, int Prior)
:
Cxxx (Countxxx), Cxxy (Countxxy), Cyxx (Countyxx), Cxyx (Countxyx)
{
  n_interval_calls = 0;
  n_real_calls = 0;
  n_hesstype_calls = 0;
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  // TotSites = real(Cxxx + Cxxy + Cyxx + Cxyx);
  Inf(PositiveProbInterval)=1e-100;
  Sup(PositiveProbInterval)=1.0;
  fxxx = real(Cxxx);// /TotSites;
  fxxy = real(Cxxy);// /TotSites;
  fyxx = real(Cyxx);// /TotSites;
  fxyx = real(Cxyx);// /TotSites;
  n_dimensions=3;
  ivector domain (1, n_dimensions);
  LabBox Ldomain;   
  /*! Sampling domain is a cube with sides given by the interval Domain */
  for (int i = 1; i <= n_dimensions; i++)
  {
    domain[i] = Domain;
  }
  Ldomain.Box = domain;
  Ldomain.L = 0;
  LabDomainList.push_back (Ldomain);
}

interval
FCFN3UnRooted::operator () (const LabBox & X)
const
{
  n_interval_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Box) != 3 || Lb(X.Box) != 1)
  {
    cout << "dimensions !=3 OR start-index != 1. exiting. " 
         << endl; exit (EXIT_FAILURE);
  }
  #endif
  interval em2t1plust2 = exp(-2.0*(X.Box[1]+X.Box[2]));
  interval em2t2plust3 = exp(-2.0*(X.Box[2]+X.Box[3]));
  interval em2t1plust3 = exp(-2.0*(X.Box[1]+X.Box[3]));

  interval lxxx = ((1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0);
  interval lxxy = ((1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0);
  interval lyxx = ((1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0);
  interval lxyx = ((1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0);

  if(Inf(lxxx)<=0.0) lxxx = lxxx & PositiveProbInterval;
  if(Inf(lxxy)<=0.0) lxxy = lxxy & PositiveProbInterval;
  if(Inf(lyxx)<=0.0) lyxx = lyxx & PositiveProbInterval;
  if(Inf(lxyx)<=0.0) lxyx = lxyx & PositiveProbInterval;

  //cout << lxxx << lxxy << lyxx << lxyx << endl; getchar();
  interval result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) + 
                    fxyx*ln(lxyx) ;
  // result = (TotSites*result);
  //result = exp(TotSites*result);
  //cout << X.Box[a] << '\t' << result << '\n'; getchar();
  /* //forcing probability intersections 
  if (Sup (result) < 0.0)
  Sup (result) = 0.0;
  if (Inf (result) < 0.0)
  Inf (result) = 0.0;
  */
  //return (UsingLogDensity) ? ln(result) : (result);
  return (UsingLogDensity) ? (result) : exp(result);
}

real
FCFN3UnRooted::operator () (const LabPnt & X)
const
{
  n_real_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Pnt) != 3 || Lb(X.Pnt) != 1)
  {
    cout << "dimensions !=3 OR start-index != 1. exiting. " << endl;
    exit (EXIT_FAILURE);
  }
  #endif
  real em2t1plust2 = exp(-2.0*(X.Pnt[1]+X.Pnt[2]));
  real em2t2plust3 = exp(-2.0*(X.Pnt[2]+X.Pnt[3]));
  real em2t1plust3 = exp(-2.0*(X.Pnt[1]+X.Pnt[3]));

  real lxxx = (1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0;
  real lxxy = (1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0;
  real lyxx = (1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0;
  real lxyx = (1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0;

  real result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) + fxyx*ln(lxyx) ;
  // result = (TotSites*result);
  //result = exp(TotSites*result);
  //cout << X.Box[a] << '\t' << result << '\n'; getchar();

  /* //forcing probability intersections
  if (Sup (result) < 0.0)
  Sup (result) = 0.0;
  if (Inf (result) < 0.0)
  Inf (result) = 0.0;
  */
  //return (UsingLogDensity) ? ln(result) : (result);
  return (UsingLogDensity) ? (result) : exp(result);
}

HessType 
FCFN3UnRooted::operator()(const HTvector& x, const int label)
const
{
  n_hesstype_calls++;
  #ifdef TESTDIMS
  if (Ub(x) != 3 || Lb(x) != 1)
  {
    cout << "dimensions !=3 OR start-index != 1. exiting. " 
         << endl; exit (EXIT_FAILURE);
  }
  #endif
  HessType em2t1plust2 = exp(-2.0*(x[1]+x[2]));
  HessType em2t2plust3 = exp(-2.0*(x[2]+x[3]));
  HessType em2t1plust3 = exp(-2.0*(x[1]+x[3]));

  HessType lxxx = ((1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0);
  HessType lxxy = ((1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0);
  HessType lyxx = ((1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0);
  HessType lxyx = ((1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0);

/*
  if(Inf(lxxx)<=0.0) lxxx = lxxx & PositiveProbInterval;
  if(Inf(lxxy)<=0.0) lxxy = lxxy & PositiveProbInterval;
  if(Inf(lyxx)<=0.0) lyxx = lyxx & PositiveProbInterval;
  if(Inf(lxyx)<=0.0) lxyx = lxyx & PositiveProbInterval;
*/
  //cout << fValue(lxxx) << fValue(lxxy) << fValue(lyxx) << fValue(lxyx) << endl; getchar();
  //HessType result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) +  fxyx*ln(lxyx) ;
  HessType result = power(lxxx,Cxxx) * power(lxxy,Cxxy) * power(lyxx,Cyxx) *  power(lxyx,Cxyx) ;
  // result = (TotSites*result);
  //result = exp(TotSites*result);
  //cout << X.Box[a] << '\t' << result << '\n'; getchar();
  /* //forcing probability intersections 
  if (Sup (result) < 0.0)
  Sup (result) = 0.0;
  if (Inf (result) < 0.0)
  Inf (result) = 0.0;
  */
  //return (UsingLogDensity) ? ln(result) : (result);
  //return (UsingLogDensity) ? (result) : exp(result);
  return result;
}


/*! constructor for the likelihood of the branchlength under CFN model for a 
  rooted and clocked 3 taxa tree
*/
// if LogPi is true, the log of this will be returned
FCFN3Rooted::FCFN3Rooted(int Countxxx, int Countxxy, int Countyxx, 
                         int Countxyx, 
                         interval Domain, bool LogPi, int Prior)
:
Cxxx (Countxxx), Cxxy (Countxxy), Cyxx (Countyxx), Cxyx (Countxyx)
{
  n_interval_calls = 0;
  n_real_calls = 0;
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  // TotSites = real(Cxxx + Cxxy + Cyxx + Cxyx);
  Inf(PositiveProbInterval)=1.e-300;
  Sup(PositiveProbInterval)=1.0;
  fxxx = real(Cxxx);// /TotSites;
  fxxy = real(Cxxy);// /TotSites;
  fyxx = real(Cyxx);// /TotSites;
  fxyx = real(Cxyx);// /TotSites;
  n_dimensions=2;

  ivector domain (1, n_dimensions);
  LabBox Ldomain;   
  /*! For each topology label, Sampling domain is a rectangle with sides 
    given by the interval Domain
  */
  for (int TopLab = 1; TopLab <= 3; TopLab++)
  {
    for (int i = 1; i <= n_dimensions; i++)
    {
      domain[i] = Domain;
    }
    Ldomain.Box = domain;
    Ldomain.L = TopLab;
    LabDomainList.push_back (Ldomain);
  }
}

// Volume of rooted tree boxes is implemented here and NOT inherited from Fobj
inline real
FCFN3Rooted::LabBoxVolume(const LabBox& LB)
{
  real volume = diam(LB.Box[1]);
  //return (volume * (volume + diam (LB.Box[2])));
  //same as inherited LabBoxVolume--CORRECT
  return (volume * diam (LB.Box[2]));
}

interval
FCFN3Rooted::operator () (const LabBox & X)
const
{
  //  cout << "interval FCFN3Rooted" << endl;
  n_interval_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Box) != 2 || Lb(X.Box) != 1)
  {
    cout << "dimensions !=2 OR start-index != 1. exiting. " 
         << endl; exit (EXIT_FAILURE);
  }
  #endif
  interval em2t1plust2;
  interval em2t2plust3;
  interval em2t1plust3;
  if (X.L == 1)     /*! Tree ((12)3)*/
  {
    interval B12 = 2.0 * X.Box[1];
    em2t1plust2 = exp(-2.0 * B12);
    em2t2plust3 = exp(-2.0*(B12 + X.Box[2]));
                    //exp(-2.0*(B12 + X.Box[2]));
    em2t1plust3 = em2t2plust3;
  }
  else if (X.L == 2)/*! Tree ((23)1)*/
  {
    interval B23 = 2.0 * X.Box[1];
    em2t1plust2 = exp(-2.0*(B23 + X.Box[2]));
    em2t2plust3 = exp(-2.0* B23);
                    //exp(-2.0*(B23 + X.Box[2]));
    em2t1plust3 = em2t1plust2;
  }
  else              /*! Tree ((13)2)*/
  {
    interval B13 = 2.0 * X.Box[1];
    em2t1plust2 = exp(-2.0*(B13 + X.Box[2]));
    em2t2plust3 = em2t1plust2;
    em2t1plust3 = exp(-2.0* B13);
  }
  // cout << "x.B0x\n" << X.L << '\n' << X.Box << '\n';
  // cout << em2t1plust2 << em2t2plust3 << em2t1plust3 
  // << endl << endl; getchar();

                    // & PositiveProbInterval;
  interval lxxx = ((1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0);
                    // & PositiveProbInterval;
  interval lxxy = ((1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0);
                    // & PositiveProbInterval;
  interval lyxx = ((1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0);
                    // & PositiveProbInterval;
  interval lxyx = ((1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0);
  //cout << lxxx << lxxy << lyxx << lxyx << endl; getchar();
  if(Inf(lxxx)<=0.0) lxxx = lxxx & PositiveProbInterval;
  if(Inf(lxxy)<=0.0) lxxy = lxxy & PositiveProbInterval;
  if(Inf(lyxx)<=0.0) lyxx = lyxx & PositiveProbInterval;
  if(Inf(lxyx)<=0.0) lxyx = lxyx & PositiveProbInterval;
  interval result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) + 
                    fxyx*ln(lxyx) ;
  //result = exp(TotSites*result);
  //cout << X.Box[a] << '\t' << result << '\n'; getchar();
  return (UsingLogDensity) ? (result) : exp(result);
}

real
FCFN3Rooted::operator () (const LabPnt & X)
const
{
  // cout << "real FCFN3Rooted" << endl;
  n_real_calls++;
  #ifdef TESTDIMS
  if (Ub(X.Pnt) != 2 || Lb(X.Pnt) != 1)
  {
    cout << "dimensions !=2 OR start-index != 1. exiting. " 
    << endl; exit (EXIT_FAILURE);
  }
  #endif
  real em2t1plust2;
  real em2t2plust3;
  real em2t1plust3;
  if (X.L == 1)     /*! Tree ((12)3)*/
  {
    real B12 = 2.0 * X.Pnt[1];
    em2t1plust2 = exp(-2.0 * B12);
    em2t2plust3 = exp(-2.0*(B12 + X.Pnt[2]));
                    //exp(-2.0*(B12 + X.Box[2]));
    em2t1plust3 = em2t2plust3;
  }
  else if (X.L == 2)/*! Tree ((23)1)*/
  {
    real B23 = 2.0 * X.Pnt[1];
    em2t1plust2 = exp(-2.0*(B23 + X.Pnt[2]));
    em2t2plust3 = exp(-2.0* B23);
                    //exp(-2.0*(B23 + X.Box[2]));
    em2t1plust3 = em2t1plust2;
  }
  else              /*! Tree ((13)2)*/
  {
    real B13 = 2.0 * X.Pnt[1];
    em2t1plust2 = exp(-2.0*(B13 + X.Pnt[2]));
    em2t2plust3 = em2t1plust2;
    em2t1plust3 = exp(-2.0* B13);
  }

  //cout << X.L    << '\n' << X.Pnt          << '\n'; getchar();
  real lxxx = (1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0;
  real lxxy = (1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0;
  real lyxx = (1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0;
  real lxyx = (1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0;

  real result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + 
                fyxx*ln(lyxx) + fxyx*ln(lxyx) ;
  // result = (TotSites*result);
  //result = exp(TotSites*result);
  //cout << X.Pnt    << '\t' << result << '\n'; getchar();
  return (UsingLogDensity) ? (result) : exp(result);
}

FCFN3::FCFN3(int Countxxx, int Countxxy, int Countyxx, int Countxyx, 
             interval Domain, bool LogPi, int Prior)
:
                    // , PriorType(Prior)
Cxxx (Countxxx), Cxxy (Countxxy), Cyxx (Countyxx), Cxyx (Countxyx)
{
  n_interval_calls = 0;
  n_real_calls = 0;
  setUsingLogDensity (LogPi);
  PriorType = Prior;
  Cid=Cxxx;
  Cnid=Cxxy + Cyxx + Cxyx;
  // TotSites = real(Cxxx + Cxxy + Cyxx + Cxyx);
  Inf(PositiveProbInterval)=1.e-300;
  Sup(PositiveProbInterval)=1.0;
  // f0 = real(Cid)// /TotSites;
  fxxx = real(Cxxx);// /TotSites;
  fxxy = real(Cxxy);// /TotSites;
  fyxx = real(Cyxx);// /TotSites;
  fxyx = real(Cxyx);// /TotSites;
  const int n_dimensions_of_label[] = {1,2,2,2,3};
  //! replace vol_labelled_domain with integral of prior over labelled domain
  //! use LabDomainPriorIntegralList
  //  vol_labelled_domain.reserve(5);
  for (int TopLab = 0; TopLab < 5; TopLab++)
  {
    ivector domain (1, n_dimensions_of_label[TopLab]);
    LabBox Ldomain; 
    /*! For each topology label, Sampling domain is a rectangle with 
      sides given by the interval Domain
    */
    //   vol_labelled_domain[TopLab]=1.0;
    for (int i = 1; i <= n_dimensions_of_label[TopLab]; i++)
    {
      domain[i] = Domain;
      //   vol_labelled_domain[TopLab] *= diam(Domain);
    }
    Ldomain.Box = domain;
    Ldomain.L = TopLab;
    LabDomainList.push_back (Ldomain);
  }
  // PriorType = 1; // 0 -> uniform prior, otherwise exponential
  SetupLabDomainPriorIntegralList();
  for(int i=0; i<5; i++)
  {
    cout << "label, domain vol1,2: " << i << "  " 
         << LabDomainPriorIntegralList[i] << endl;
  }

}

//! Volume of rooted tree boxes is implemented here and NOT inherited from Fobj
inline real
FCFN3::LabBoxVolume(const LabBox& LB)
{
  if(LB.L==0) return diam(LB.Box[1]);
  else if(LB.L==4)
  {
    real volume = 1.0;
    for (int i = 1; i <= VecLen(LB.Box); i++) volume *= diam (LB.Box[i]);
    return volume;
  }
  else
  {
    real volume = diam(LB.Box[1]);
    //return (volume * (volume + diam (LB.Box[2])));
                    //same as inherited LabBoxVolume--CORRECT
    return (volume * diam (LB.Box[2]));
  }
}

interval
FCFN3::operator () (const LabBox & X)
const
{ /*! \todo divide by volume of labelled domain -- volume of rooted trees 
    need transformation (specalise Fobj member function get_volume_LabBox)
   */
  // cout << "interval FCFN3" << endl;
  n_interval_calls++;
  interval result;
  if (X.L==0)       /*! Star Tree (123)*/
  {
    #ifdef TESTDIMS
    if (Ub(X.Box) != 1 || Lb(X.Box) != 1)
    {
      cerr << "dimensions !=1 OR start-index != 1... exiting. " 
           << endl; exit (EXIT_FAILURE);
    }
    #endif
    interval em4t = (exp (-4.0 * X.Box[1])) ;
    result = Cid*ln( ((1.0 + (3.0*em4t)) / 8.0)) + 
             (Cnid)*ln( ((1.0 - em4t) / 8.0)) ;
  }
  else      // un/rooted trees with > 1 distinct branch length
  {
    interval lxxx;
    interval lxxy;
    interval lyxx;
    interval lxyx;
    interval em2t1plust2;
    interval em2t2plust3;
    interval em2t1plust3;
    if (X.L==4)     /*! UnRooted Tree (1:,2:,3:)*/
    {
      #ifdef TESTDIMS
      if (Ub(X.Box) != 3 || Lb(X.Box) != 1)
      {
        cout << "dimensions !=3 OR start-index != 1. exiting. " 
             << endl; exit (EXIT_FAILURE);
      }
      #endif
      em2t1plust2 = exp(-2.0*(X.Box[1]+X.Box[2]));
      em2t2plust3 = exp(-2.0*(X.Box[2]+X.Box[3]));
      em2t1plust3 = exp(-2.0*(X.Box[1]+X.Box[3]));
      //cout << lxxx << lxxy << lyxx << lxyx << endl; getchar();
    }
    else            // Rooted trees
    {
      #ifdef TESTDIMS
      if (Ub(X.Box) != 2 || Lb(X.Box) != 1)
      {
        cout << "dimensions !=2 OR start-index != 1. exiting. " 
             << endl; exit (EXIT_FAILURE);
      }
      #endif
      if (X.L == 1) /*! Tree ((12)3)*/
      {
        interval B12 = 2.0 * X.Box[1];
        em2t1plust2 = exp(-2.0 * B12);
        em2t2plust3 = exp(-2.0*(B12 + X.Box[2]));
                    //exp(-2.0*(B12 + X.Box[2]));
        em2t1plust3 = em2t2plust3;
      }
                    /*! Tree ((23)1)*/
      else if (X.L == 2)
      {
        interval B23 = 2.0 * X.Box[1];
        em2t1plust2 = exp(-2.0*(B23 + X.Box[2]));
        em2t2plust3 = exp(-2.0* B23);
                    //exp(-2.0*(B23 + X.Box[2]));
        em2t1plust3 = em2t1plust2;
      }
      else          /*! Tree ((13)2)*/
      {
        interval B13 = 2.0 * X.Box[1];
        em2t1plust2 = exp(-2.0*(B13 + X.Box[2]));
        em2t2plust3 = em2t1plust2;
        em2t1plust3 = exp(-2.0* B13);
      }
    }
    lxxx = ((1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0);
    lxxy = ((1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0);
    lyxx = ((1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0);
    lxyx = ((1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0);
    if(Inf(lxxx)<=0.0) lxxx = lxxx & PositiveProbInterval;
    if(Inf(lxxy)<=0.0) lxxy = lxxy & PositiveProbInterval;
    if(Inf(lyxx)<=0.0) lyxx = lyxx & PositiveProbInterval;
    if(Inf(lxyx)<=0.0) lxyx = lxyx & PositiveProbInterval;
    result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) + fxyx*ln(lxyx) ;
  }
  //result = exp(TotSites*result);
  //cout << X.Box[a] << '\t' << result << '\n'; getchar();
  //return (UsingLogDensity) ? ln(result) : (result);
  // result = (TotSites*result);
  //	 cout << "interval result: " << result << endl;

  if(FATTEN_THIN_INTERVAL_RE)
  {
    if(UsingLogDensity)
    {
      Inf(result) -= 1.0e-5;
    }
    else
    {
      Inf(result) *= (1.0 - 1.0e-5);
    }
  }
  //
  //	 cout << "X.L: " << X.L << "  result: " << result << "  ln(vollabdom): " 
  //   <<  ln(vol_labelled_domain[X.L]) << endl;
  //	 result = result - ln(vol_labelled_domain[X.L]);
  // The following line has to do with the prior - could this go in
  // Fobj for the basic (unif and exponential) priors?
  result = result - ln(LabDomainPriorIntegralList[X.L]);
  return (UsingLogDensity) ? (result) : exp(result);
}

real
FCFN3::operator () (const LabPnt & X)
const
{
  //  cout << "real FCFN3, top" << endl;
  n_real_calls++;
  real result;
  if (X.L==0)       /*! Star Tree (123)*/
  {
    #ifdef TESTDIMS
    if (Ub(X.Pnt) != 1 || Lb(X.Pnt) != 1)
    {
      cerr << "dimensions !=1 OR start-index != 1... exiting. " 
           << endl; exit (EXIT_FAILURE);
    }
    #endif
    real em4t = (exp (-4.0 * X.Pnt[1])) ;
    result = Cid*ln( ((1.0 + (3.0*em4t)) / 8.0)) + 
             (Cnid)*ln( ((1.0 - em4t) / 8.0)) ;
  }
  else              // un/rooted trees with > 1 distinct branch length
  {
    real lxxx;
    real lxxy;
    real lyxx;
    real lxyx;
    real em2t1plust2;
    real em2t2plust3;
    real em2t1plust3;
    if (X.L==4)     /*! UnRooted Tree (1:,2:,3:)*/
    {
      #ifdef TESTDIMS
      if (Ub(X.Pnt) != 3 || Lb(X.Pnt) != 1)
      {
        cerr << "dimensions !=3 OR start-index != 1. exiting. " 
             << endl; exit (EXIT_FAILURE);
      }
      #endif
      em2t1plust2 = exp(-2.0*(X.Pnt[1]+X.Pnt[2]));
      em2t2plust3 = exp(-2.0*(X.Pnt[2]+X.Pnt[3]));
      em2t1plust3 = exp(-2.0*(X.Pnt[1]+X.Pnt[3]));
    }
    else            // Rooted trees
    {
      #ifdef TESTDIMS
      if (Ub(X.Pnt) != 2 || Lb(X.Pnt) != 1)
      {
        cout << "dimensions !=2 OR start-index != 1. exiting. " 
             << endl; exit (EXIT_FAILURE);
      }
      #endif
      if (X.L == 1) /*! Tree ((12)3)*/
      {
        real B12 = 2.0 * X.Pnt[1];
        em2t1plust2 = exp(-2.0 * B12);
        em2t2plust3 = exp(-2.0*(B12 + X.Pnt[2]));
                    //exp(-2.0*(B12 + X.Pnt[2]));
        em2t1plust3 = em2t2plust3;
      }
                    /*! Tree ((23)1)*/
      else if (X.L == 2)
      {
        real B23 = 2.0 * X.Pnt[1];
        em2t1plust2 = exp(-2.0*(B23 + X.Pnt[2]));
        em2t2plust3 = exp(-2.0* B23);
                    //exp(-2.0*(B23 + X.Pnt[2]));
        em2t1plust3 = em2t1plust2;
      }
      else          /*! Tree ((13)2)*/
      {
        real B13 = 2.0 * X.Pnt[1];
        em2t1plust2 = exp(-2.0*(B13 + X.Pnt[2]));
        em2t2plust3 = em2t1plust2;
        em2t1plust3 = exp(-2.0* B13);
      }
    }
    lxxx = ((1.0 + em2t1plust2 + em2t2plust3 + em2t1plust3) / 8.0);
    lxxy = ((1.0 + em2t1plust2 - em2t2plust3 - em2t1plust3) / 8.0);
    lyxx = ((1.0 - em2t1plust2 + em2t2plust3 - em2t1plust3) / 8.0);
    lxyx = ((1.0 - em2t1plust2 - em2t2plust3 + em2t1plust3) / 8.0);
    //	cout << lxxx << "  " << lxxy << "  " << lyxx << "  " << lxyx << "  " 
    //       << endl; getchar();
    //	cout << fxxx << "  " << fxxy << "  " << fyxx << "  " << fxyx << endl;
    result = fxxx*ln(lxxx) + fxxy*ln(lxxy) + fyxx*ln(lyxx) + fxyx*ln(lxyx) ;
  }
  //result = exp(TotSites*result);
  //	cout << X.Pnt << '\t' << result << '\n'; getchar();
  //return (UsingLogDensity) ? ln(result) : (result);
  // result = (TotSites*result);
  //	cout << "real result: " << result << endl;
  //	 cout << "X.L:  " << X.L << " result: " << result << "  ln(vollabdom): " 
  //        <<  ln(vol_labelled_domain[X.L]) << endl;
  // result = result - ln(vol_labelled_domain[X.L]);
  result = result - ln(LabDomainPriorIntegralList[X.L]);
  return (UsingLogDensity) ? (result) : exp(result);
}

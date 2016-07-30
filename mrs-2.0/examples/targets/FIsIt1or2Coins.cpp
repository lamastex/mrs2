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
/*! \file FIsIt1or2Coins.cpp
\brief Trans-dimensional posterior distribution over prob. of heads for possibly two coins
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
#include "FIsIt1or2Coins.hpp"

/*! constructor for the the toy example of a target which is a mixture of two densities
   on spaces of different dimensionality.
   There are two set of coin flips N1 flips of coin 1 with n1 heads,
   N2 flips of coin 2 with n2 heads.
   The target is a mixture of a model with 2 parameters, 
   q1 = prob. of head for coin 1, q2 = prob. of heads for coin2,
   and a model with 1 parameter q = the prob. of heads for both coins.
*/
FIsIt1or2Coins::FIsIt1or2Coins(int No1, int no1, int No2, int no2, bool LogPi):N1(No1), n1(no1), N2(No2), n2(no2)
{
    n_interval_calls = 0;
    n_real_calls = 0;
    setUsingLogDensity (LogPi);
    PriorType = 0; // unif. prior is default
    SmallestPositive = 1e-300;
    IgnoreBinomCoeffs = true; // ignore the binomial coefficients if set to true

    for(int Label = 0; Label <= 1; Label++){
      int n_dimensions = Label+1; // dimensionality of this labeled space
      ivector domain (1, n_dimensions); // vector of intervals, with indices 1...2
      LabBox Ldomain;
      for (int i = 1; i <= n_dimensions; i++){
        //if working with log ikelihood then domain should be a 
        //proper subset of [0,1] to ensure local Lipschitz condition
        if(LogPi) 
          domain[i] = interval (SmallestPositive, 1.0 - SmallestPositive);
        else domain[i] = interval (0.0, 1.0);
      }
      Ldomain.Box = domain;
      Ldomain.L = Label;
      cout << "Ldomain.L: " << Ldomain.L << endl;
      LabDomainList.push_back (Ldomain);
    }
    w1 = 1.0;
    w2 = 1.0;
    if (IgnoreBinomCoeffs){
      bc1=interval(1.0,1.0);
      bc2=interval(1.0,1.0);
    }
    else {
      bc1 = binary_coefficient(N1, n1); 
      bc2 = binary_coefficient(N2, n2);
    }
    // bc = binary_coefficient(N1 + N2, n1 + n2);
    rbc1 = 0.5*(Inf(bc1) + Sup(bc1));
    rbc2 = 0.5*(Inf(bc2) + Sup(bc2));
    // rbc = 0.5*(Inf(bc) + Sup(bc));
    cout << "end of FIsIt1or2Coins constructor. \n";
}

// vector<LabBox> FIsIt1or2Coins::get_domain(){ return LabDomainList; }

interval
FIsIt1or2Coins::operator () (const LabBox & X)
const
{
    // L(N1, n1, q1) = power(q1, n1) * power(1.0 - q1, N1 - n1) * bc1;
    n_interval_calls++;

    interval result;
    if(X.L == 0){
        interval prior_density(1.0, 1.0); // 1.0/(1-d volume) 
        if(UsingLogDensity){
          interval q = X.Box[Lb(X.Box)];
          interval qcomp = 1.0 - q; if(Inf(qcomp) <= 0.0){ SetInf(qcomp, SmallestPositive); }
          interval LogLikelihood = (n1 + n2)*ln(q) + (N1+N2 - (n1+n2))*ln(qcomp) + ln(bc1)+ln(bc2);
          result = ln(w1) + LogLikelihood + ln(prior_density);
        }
        else{
          interval q = X.Box[Lb(X.Box)];
          interval qcomp = 1.0 - q; if(Inf(qcomp) <= 0.0){ SetInf(qcomp, 0.0); }
          interval likelihood = power(q, n1+n2) * power(qcomp, N1+N2 - (n1+n2)) * bc1*bc2;
          result = w1*likelihood*prior_density;
        }
    }
    else if(X.L == 1){
        interval prior_density(1.0, 1.0); // 1.0/(2-d volume)
        interval q1 = X.Box[Lb(X.Box)];
        interval q2 = X.Box[Lb(X.Box) + 1];
        if(UsingLogDensity){
          interval q1comp = 1.0 - q1; if(Inf(q1comp) <= 0.0){ SetInf(q1comp, SmallestPositive); }
          interval q2comp = 1.0 - q2; if(Inf(q2comp) <= 0.0){ SetInf(q2comp, SmallestPositive); }
          interval LogLikelihood = n1*ln(q1) + (N1 - n1)*ln(q1comp) + ln(bc1);
          LogLikelihood += n2*ln(q2) + (N2 - n2)*ln(q2comp) + ln(bc2);
          result = ln(w2) + LogLikelihood + ln(prior_density);
        }
        else{
          interval q1comp = 1.0 - q1; if(Inf(q1comp) <= 0.0){ SetInf(q1comp, 0.0); }
          interval q2comp = 1.0 - q2; if(Inf(q2comp) <= 0.0){ SetInf(q2comp, 0.0); }
          interval likelihood = power(q1, n1) * power(q1comp, N1 - n1) * bc1;
          likelihood *= power(q2, n2) * power(q2comp, N2 - n2) * bc2;
          result = w2*likelihood*prior_density;
        }
    }
    //cout << result << '\n';getchar();
    return result;
}

real
FIsIt1or2Coins::operator () (const LabPnt & X)
    const
{    
    n_real_calls++;
    real result;
    if(X.L == 0){
        real prior_density = 1.0; // 1.0/(1-d volume) 
        real q = X.Pnt[Lb(X.Pnt)];
       
        if(UsingLogDensity){
            real LogLikelihood = (n1 + n2)*ln(q) + (N1+N2 - (n1+n2))*ln(1.0 - q) + ln(rbc1)+ln(rbc2);
            result = ln(w1) + LogLikelihood + ln(prior_density);
        }
        else{
            real likelihood = power(q, n1+n2) * power(1.0 - q, N1+N2 - (n1+n2)) * rbc1*rbc2;
            result = w1*likelihood*prior_density;
        }
    }
    else if(X.L == 1){
        real prior_density = 1.0; // 1.0/(2-d volume)
        real q1 = X.Pnt[Lb(X.Pnt)];
        real q2 = X.Pnt[Lb(X.Pnt) + 1];
        if(UsingLogDensity){
            real LogLikelihood = n1*ln(q1) + (N1 - n1)*ln(1.0 - q1) + ln(rbc1);
            LogLikelihood += n2*ln(q2) + (N2 - n2)*ln(1.0 - q2) + ln(rbc2);
            result = ln(w2) + LogLikelihood + ln(prior_density);
        }
        else{
            real likelihood = power(q1, n1) * power(1.0 - q1, N1 - n1) * rbc1;
            likelihood *= power(q2, n2) * power(1.0 - q2, N2 - n2) * rbc2;
            result = w2*likelihood*prior_density;
        }
    }
    return result;
}

interval
FIsIt1or2Coins::binary_coefficient(int N, int n)
{
    int kmax = (n <= N - n)? n: N - n;
    int nf = N; 
    int df = 1;
    interval nmrtr(1,1);
    interval dnmntr(1,1);
    for(int k = 1; k <= kmax; k++){
        nmrtr *= nf; nf--;
        dnmntr *= df; df++;
    }
    return nmrtr/dnmntr;
} 

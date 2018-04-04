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

#ifndef __FSHIRYAEV1D__
#define __FSHIRYAEV1D__

#define PI (_interval(3.14159265358979323846264338328))
#define Pi (3.14159265358979323846264338328)
//lambda parameter is fixed at 9/20 and we compute the corresponding 
// Gamma(1+1/lambda) for likelihood comps
#define BB  (_interval(9.0/20.0))
#define bb (9.0/20.0)
#define GAMMA1PLUS1OVERBB (_interval(2.47859397716744434365241145010))
#define Gamma1plus1overbb (2.47859397716744434365241145010)

//! one-dimensional Shiryaev density as a function object class
class FShiryaev1D: public Fobj
{
  real aa;
  real cc;
  interval AA;
  interval CC;
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  public:
    FShiryaev1D(real aa, real cc, real DomainLimit, bool LogPi, int Prior);
    interval operator()(const LabBox& X) const;
    real operator()(const LabPnt& X) const ;
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};

//! one-dimensional Shiryaev likelhood as a function object class
class FShiryaev1D_Lkl_aa_fromData: public Fobj
{
  RSSample & Data;  //!< Store the data to estimate parameter aa
  real cc;          //!< cc is assumed to be known
  interval CC;
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  public:
    FShiryaev1D_Lkl_aa_fromData(RSSample & d, real cc, 
                                interval DomainInterval, bool LogPi, int Prior);
    interval operator()(const LabBox& A) const;
    real operator()(const LabPnt& A) const ;
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};
#endif

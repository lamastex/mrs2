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

/*! \file      FLevy2D.hpp
\brief Declarations for example function class FLevy2D 
(Levy function, 2 dimensions).
*/

#ifndef __FLEVY2D__
#define __FLEVY2D__

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"
#include "imath.hpp"
#include "rmath.hpp"
#include "ivector.hpp"
#include <functional>

#include <gop.hpp>  // cxsc global optimisation

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"

//! two-dimensional Levy density as a function object class
class FLevy2D: public Fobj
{
  real Temperature, GlobalMax, Center1, Center2;

  //! Track number of interval function calls
  mutable int n_interval_calls;

  //! Track number of real function calls
  mutable int n_real_calls;
  public:
    //! Constructor
    FLevy2D(real Temperature, real GlobalMax, real Center1, real Center2, 
            real DomainLimit, bool LogPi);

    //! interval operator()
    interval operator()(const LabBox& X) const;

    //! real operator()
    real operator()(const LabPnt& X) const;

    //! HessType operator()
    HessType operator()(const HTvector& x, const int label = 0) const;

    //! get volume of a labeled box
    virtual real LabBoxVolume(const LabBox& LB)
    {
      return Fobj::LabBoxVolume(LB);
    }

    //!Get number of interval function calls
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

//! two-dimensional Levy likelihood as a function object class
//! \todo get rid of Fobj1D and make these inherit from Fobj as FCFN3Star
class FLevy2D_Lkl_Tfrom1data: public Fobj1D
{
  ivector Data;
  real GlobalMax, Center1, Center2;
  public:
    FLevy2D_Lkl_Tfrom1data(ivector& Data, real GlobalMax, real Center1, real Center2, interval DomainLimit, bool LogPi);
    interval operator()(const interval& T) const;
    real operator()(const real& T) const ;
};
#endif

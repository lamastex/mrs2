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

/*! \file      FRosenbrock.hpp
\brief Declarations for example function class FRosenbrock (Rosenbrock function).
*/

#ifndef __FROSENBROCK_HPP__
#define __FROSENBROCK_HPP__

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
                    // Include interval arithmetic package
#include "interval.hpp"
#include "imath.hpp"// Include interval standard functions
#include "rmath.hpp"// Include real standard functions
#include "ivector.hpp"
#include <functional>

#include <gop.hpp>  // cxsc global optimisation

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"

//! n-dimensional Rosenbrock density as a function object class
class FRosenbrock: public Fobj
{

  int num_dim;
  real Tinverse;
  real Height;

                    //!< Track number of interval function calls
  mutable int n_interval_calls;
                    //!< Track number of real function calls
  mutable int n_real_calls;

  public:
    //! Constructor
    FRosenbrock(int n_dimensions, real T, real H,real DomainLimit, bool LogPi);

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

};                  // end of FRosenbrock declarations
#endif

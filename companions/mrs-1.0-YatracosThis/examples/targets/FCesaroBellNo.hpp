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

#ifndef __FCESAROBELLNO_HPP__
#define __FCESAROBELLNO_HPP__

/*! \file      FCesaroBellNo.hpp
\brief Declarations for example function class FCesaroBellNo.

 Absolute integrand shape in Cesaro's integral for the n-the Bell 
 number as a function object class. 
*/
// example function object class to use with MRSampler class

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"		// Include interval arithmetic package
#include "imath.hpp"		// Include interval standard functions
#include "rmath.hpp"		// Include real standard functions
#include "ivector.hpp"
#include <functional>

using namespace std;
using namespace cxsc;


#include "SmallClasses.hpp"
#include "Fobj.hpp"

class FCesaroBellNo: public Fobj{
    real n;
    mutable int n_interval_calls;
    mutable int n_real_calls;
public:
    FCesaroBellNo(real Nn, interval DomainLimit, bool LogPi);
    interval operator()(const LabBox& X) const;
    real operator()(const LabPnt& X) const;
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    int get_interval_calls(){ return n_interval_calls; }
    int get_real_calls(){ return n_real_calls; }
};

#endif

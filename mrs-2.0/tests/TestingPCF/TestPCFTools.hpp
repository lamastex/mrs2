/*
* Copyright (C) 2011, 2012 Jennifer Harlow
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


/*! \file
\brief Testing PiecewiseConstantFunction tools
 */

#ifndef __TESTPCFTOOLS_HPP__
#define __TESTPCFTOOLS_HPP__


#include "sptypes.hpp"
#include "cxsc.hpp"

#include <vector> 

std::vector< real > makeRanges1(cxsc::real rootvol);
std::vector< real > makeRanges2(cxsc::real rootvol);
std::vector< real > makeRanges3(cxsc::real rootvol);
std::vector< real > makeRanges4(cxsc::real rootvol);
std::vector< real > makeRanges5(cxsc::real rootvol);
std::vector< real > makeRanges6(cxsc::real rootvol);
std::vector< real > makeRangesNegative(cxsc::real rootvol);
std::vector< real > makeRangesInfinite();
std::vector< real > makeRangesPositiveAndNegative();
std::vector< real > makeRangesPositiveAndNegativeInfinite();
std::vector< real > makeRangesInfiniteAndNegativeInfinite();
std::vector< real > makeRangesPositiveAndZero();
std::vector< real > makeRangesInfiniteAndNegative();
std::vector< real > makeRangesInfiniteAndZero();
std::vector< real > makeRangesNegativeAndZero();

std::vector< real > makeRangesArithmeticSpecial1(cxsc::real rootvol);
std::vector< real > makeRangesArithmeticSpecial2(cxsc::real rootvol);

subpavings::RVecData& getData1(subpavings::RVecData& data);
subpavings::RVecData& getData2(subpavings::RVecData& data);
subpavings::RVecData& getData2LeftOnly(subpavings::RVecData& data);
subpavings::RVecData& getData2RightOnly(subpavings::RVecData& data);
subpavings::RVecData& getData3(subpavings::RVecData& data);
subpavings::RVecData& getDataExtra1(subpavings::RVecData& data);
subpavings::RVecData& getDataExtra2(subpavings::RVecData& data);

bool checkSame(cxsc::real r1, const cxsc::real r2, int n); 

		
#endif

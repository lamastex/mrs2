/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief Definitions for multivariate Gaussian example function object class.

*/

#include "GaussianFobj.hpp"

#include <iostream>

using namespace cxsc;
using namespace std;

subpavings::GaussianFobj::GaussianFobj(){}

interval subpavings::GaussianFobj::operator()(
								const cxsc::ivector& ivec) const
{
	int lb = Lb (ivec);
	int ub = Ub (ivec);
	int dim = ub - lb + 1;
	
	cxsc::interval twoPi = interval(2.0)*cxsc::Pi();
	
	interval a = cxsc::pow(twoPi, cxsc::interval(dim/2.0));
	interval b = cxsc::power(ivec[lb],2);
	for (int i = lb+1; i <= ub; ++i) {
		b += cxsc::power(ivec[i],2);
	}
	b *= interval(-0.5);
	
	interval result = cxsc::interval(1.0)/a * cxsc::exp(b);
	
	return result;
}

real subpavings::GaussianFobj::operator()(const cxsc::rvector& r) const
{
	int lb = Lb (r);
	int ub = Ub (r);
	int dim = ub - lb + 1;
	
	real a = cxsc::pow(cxsc::Pi2_real, cxsc::real(dim/2.0));
	real b = cxsc::power(r[lb],2);
	for (int i = lb+1; i <= ub; ++i) {
		b += cxsc::power(r[i],2);
	}
	
	b*= real(-0.5);
	
	real result = cxsc::real(1.0)/a * cxsc::exp(b);

	return result;
	
}

std::string subpavings::GaussianFobj::getName() const
{
	return std::string("Gaussian");
}

subpavings::GaussianFobj::~GaussianFobj(){}

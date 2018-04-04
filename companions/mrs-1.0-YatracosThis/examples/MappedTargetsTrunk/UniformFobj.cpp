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
\brief Definitions for multivariate Uniform example function object class.

*/

#include "UniformFobj.hpp"
#include "toolz.hpp"

#include <iostream>

using namespace cxsc;
using namespace std;

subpavings::UniformFobj::UniformFobj(const cxsc::ivector& ivec)
	: domain( ivec ) {}

interval subpavings::UniformFobj::operator()(
								const cxsc::ivector& ivec) const
{
	int lb = Lb (ivec);
	int ub = Ub (ivec);
	int dim = ub - lb + 1;
	
	cxsc::interval result = 1.0/interval( realVolume(domain) );
	
	for (int i = lb; i <= ub; ++i) {
		if (Inf(ivec[i]) < Inf(domain[i]) || Sup(ivec[i]) > Sup(domain[i]) )
			SetInf(result, 0.0);
		if (Inf(ivec[i]) > Sup(domain[i]) || Sup(ivec[i]) < Inf(domain[i]) )
			SetSup(result, 0.0);
		
	}
	
	return result;
}

real subpavings::UniformFobj::operator()(const cxsc::rvector& r) const
{
	int lb = Lb (r);
	int ub = Ub (r);
	int dim = ub - lb + 1;
	
	cxsc::real result = 1.0/realVolume(domain);
	
	for (int i = lb; i <= ub && (result != 0.0); ++i) {
		if (r[i] < Inf(domain[i]) || r[i] > Sup(domain[i]) )
			result = 0.0;
		if (r[i] > Sup(domain[i]) || r[i] < Inf(domain[i]) )
			result = 0.0;
		
	}
	return result;
	
}

std::string subpavings::UniformFobj::getName() const
{
	return std::string("Uniform");
}

subpavings::UniformFobj::~UniformFobj(){}

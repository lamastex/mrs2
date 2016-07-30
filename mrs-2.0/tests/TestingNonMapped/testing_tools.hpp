/*
* Copyright (C) 2011 Jennifer Harlow
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
\brief Testing tools declarations
 */

#ifndef __TESTINGTOOLS_HPP__
#define __TESTINGTOOLS_HPP__

#include "sptypes.hpp"
#include "adaptivehistogram.hpp"
#include "spsnode.hpp"

#include <gsl/gsl_matrix.h>

#include <iostream>
#include <string>


bool checkFileLines(const std::string& s, std::size_t expectedLines);

subpavings::RVecData combineData(const subpavings::RVecData& v1, const subpavings::RVecData& v2);

subpavings::RVecData combineData(const subpavings::RVecData& v1, const subpavings::RVecData& v2,  const subpavings::RVecData& v3);

void outputADH(const std::string& s, const subpavings::AdaptiveHistogram& adh);

void outputSPS(const std::string& s, const subpavings::SPSnode& spn, int prec = 5);

void swapCheckOutput(const std::string& s, const subpavings::AdaptiveHistogram& adh);

void swapCheckOutput(const std::string& s, const subpavings::SPSnode& spn, bool append = false);

void checkOutput(const std::string& s, const subpavings::AdaptiveHistogram& adh); 

void doubleCheckOutput(const std::string& s, const subpavings::SPSnode& spn, bool append);

void doubleCheckOutput(const std::string& s, const subpavings::AdaptiveHistogram& adh);

void cov_calculate_gsl(gsl_matrix *r, gsl_matrix *m);		

subpavings::RealVec checkVarCov(const subpavings::RVecData& data); 

cxsc::rvector checkMean(const subpavings::RVecData& data); 

bool checkAllNaN(const subpavings::RealVec& vcov); 

bool checkAllNaN(const cxsc::rvector& rv);

bool checkSame(const cxsc::rvector& rv1, const cxsc::rvector& rv2, int n);

bool checkSame(const subpavings::RealVec& vec1, const subpavings::RealVec& vec2, int n);

bool checkSame(cxsc::real r1, const cxsc::real r2, int n) ;

#endif

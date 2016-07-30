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
\brief Definitions for classes for evaluating when to stop
 changing function estimators (using intervals).
*/

#include "fei_evalobj.hpp"
#include "functionestimator_interval.hpp"

using namespace subpavings;


/*Class for testing the number of leaves in estimator.
*/
FEICritLeaves_GTE::FEICritLeaves_GTE(size_t t) : test(t) {}

bool FEICritLeaves_GTE::operator()
			(const FunctionEstimatorInterval& fei) const
{
	return (fei.getRootLeaves() >= test);
}

/*Class for testing the number of leaves in estimator.
*/
FEICritLeaves_LTE::FEICritLeaves_LTE(size_t t) : test(t) {}

bool FEICritLeaves_LTE::operator()(
			const FunctionEstimatorInterval& fei) const
{
	return (fei.getRootLeaves() <= test);
}


/*Class for testing interval band over estimator.
*/
FEICritIntervalBand_LTE::FEICritIntervalBand_LTE(cxsc::real t) : test(t) {}

bool FEICritIntervalBand_LTE::operator()(
			const FunctionEstimatorInterval& fei) const
{
	return (fei.getTotalAreaOfIntervalBand() <= test);
}


/*Class to bale out of priority queue splitting.
*/
bool FEICritStopAll::operator()(
			const FunctionEstimatorInterval&) const
	{ return true; }
	



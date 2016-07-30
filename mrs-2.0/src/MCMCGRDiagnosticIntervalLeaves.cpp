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
\brief Definitions for an object to do Gelman-Rubin diagnostics for MCMC
using leaves as the scalar value.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCGRDiagnosticIntervalLeaves.hpp"

using namespace cxsc;
using namespace subpavings;
using namespace std;

			
// concrete ones only
const std::string MCMCGRDiagnosticIntervalLeaves::GRScalarFilename 
					= "GelmanRubinIntervalLeaves";
const std::string MCMCGRDiagnosticIntervalLeaves::GRWorkingCalcsFilename 
					= "GelmanRubinIntervalLeavesWorkingCalcs";
const std::string MCMCGRDiagnosticIntervalLeaves::scalarsName = "IvalLvs";
const std::string MCMCGRDiagnosticIntervalLeaves::baseScalarsColName
			= "Leaves";


MCMCGRDiagnosticIntervalLeaves::MCMCGRDiagnosticIntervalLeaves(
											cxsc::real t,
											size_t si,
											double percent,
											bool req)
	: MCMCGRDiagnosticInterval(t, si, percent, req)
{}

MCMCGRDiagnosticIntervalLeaves::MCMCGRDiagnosticIntervalLeaves(
											cxsc::real t,
											size_t si,
											size_t sm,
											double percent,
											bool req)
	: MCMCGRDiagnosticInterval(t, si, sm, percent, req)
{}

MCMCGRDiagnosticIntervalLeaves::~MCMCGRDiagnosticIntervalLeaves()
{}

std::string MCMCGRDiagnosticIntervalLeaves::getGRWorkingCalcsFilename() const
{
	return GRWorkingCalcsFilename;
}

std::string MCMCGRDiagnosticIntervalLeaves::getGRDiagnosticsFilename() const
{
	return GRScalarFilename;
}

std::string MCMCGRDiagnosticIntervalLeaves::getScalarsName() const
{
	return scalarsName;
}



// private


real MCMCGRDiagnosticIntervalLeaves::getScalarValue(
	const ChangeOfStateInformationAutoMCMC& info) const
{
	return 1.0*info.getCurrentLeaves();
}
					
std::string MCMCGRDiagnosticIntervalLeaves::getBaseScalarsColName() const
{
	return baseScalarsColName;
}

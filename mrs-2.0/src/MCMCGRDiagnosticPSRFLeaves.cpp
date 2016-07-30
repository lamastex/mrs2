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

#include "MCMCGRDiagnosticPSRFLeaves.hpp"

using namespace cxsc;
using namespace subpavings;
using namespace std;

			
// concrete ones only
const std::string MCMCGRDiagnosticPSRFLeaves::GRScalarFilename 
					= "GelmanRubinPSRFLeaves";
const std::string MCMCGRDiagnosticPSRFLeaves::GRWorkingCalcsFilename 
					= "GelmanRubinPSRFLeavesWorkingCalcs";
const std::string MCMCGRDiagnosticPSRFLeaves::scalarsName = "PSRFLvs";
const std::string MCMCGRDiagnosticPSRFLeaves::baseScalarsColName
			= "Leaves";

MCMCGRDiagnosticPSRFLeaves::MCMCGRDiagnosticPSRFLeaves(cxsc::real t,
														bool req)
	: MCMCGRDiagnosticPSRF(t, req)
{}

MCMCGRDiagnosticPSRFLeaves::~MCMCGRDiagnosticPSRFLeaves()
{}

std::string MCMCGRDiagnosticPSRFLeaves::getGRWorkingCalcsFilename() const
{
	return GRWorkingCalcsFilename;
}

std::string MCMCGRDiagnosticPSRFLeaves::getGRDiagnosticsFilename() const
{
	return GRScalarFilename;
}

std::string MCMCGRDiagnosticPSRFLeaves::getScalarsName() const
{
	return scalarsName;
}


// private


real MCMCGRDiagnosticPSRFLeaves::getScalarValue(
	const ChangeOfStateInformationAutoMCMC& info) const
{
	return 1.0*info.getCurrentLeaves();
}
					
					
std::string MCMCGRDiagnosticPSRFLeaves::getBaseScalarsColName() const
{
	return baseScalarsColName;
}					
					

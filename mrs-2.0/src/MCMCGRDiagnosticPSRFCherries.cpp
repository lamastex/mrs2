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
using cherries as the scalar value.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCGRDiagnosticPSRFCherries.hpp"

using namespace cxsc;
using namespace subpavings;
using namespace std;

			
// concrete ones only
const std::string MCMCGRDiagnosticPSRFCherries::GRScalarFilename 
					= "GelmanRubinPSRFCherries";
const std::string MCMCGRDiagnosticPSRFCherries::GRWorkingCalcsFilename 
					= "GelmanRubinPSRFCherriesWorkingCalcs";
const std::string MCMCGRDiagnosticPSRFCherries::scalarsName = "PSRFChrs";
const std::string MCMCGRDiagnosticPSRFCherries::baseScalarsColName
			= "Cherries";

MCMCGRDiagnosticPSRFCherries::MCMCGRDiagnosticPSRFCherries(cxsc::real t,
														bool req)
	: MCMCGRDiagnosticPSRF(t, req)
{}

MCMCGRDiagnosticPSRFCherries::~MCMCGRDiagnosticPSRFCherries()
{}

std::string MCMCGRDiagnosticPSRFCherries::getGRWorkingCalcsFilename() const
{
	return GRWorkingCalcsFilename;
}

std::string MCMCGRDiagnosticPSRFCherries::getGRDiagnosticsFilename() const
{
	return GRScalarFilename;
}

std::string MCMCGRDiagnosticPSRFCherries::getScalarsName() const
{
	return scalarsName;
}


// private

real MCMCGRDiagnosticPSRFCherries::getScalarValue(
	const ChangeOfStateInformationAutoMCMC& info) const
{
	return 1.0*info.getCurrentCherries();
}
					
std::string MCMCGRDiagnosticPSRFCherries::getBaseScalarsColName() const
{
	return baseScalarsColName;
}

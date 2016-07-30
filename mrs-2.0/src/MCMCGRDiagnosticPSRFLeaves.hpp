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
\brief Declarations for an object to do Gelman-Rubin diagnostics for MCMC
using leaves as the scalar value.
 */

#ifndef __MCMCGRAUTO_DIAG_PSRF_LEAVES_HPP__
#define __MCMCGRAUTO_DIAG_PSRF_LEAVES_HPP__

#include "MCMCGRDiagnosticPSRF.hpp"

namespace subpavings {
	
	class MCMCGRDiagnosticPSRFLeaves : public MCMCGRDiagnosticPSRF {
		
		public:

			/*! \brief Constructor.
			 
			\param t The tolerance to use for the 
			Gelman-Rubin convergence criteria.  
			\param req flag for whether this
			diagnostic must be satisfied for burn-in. Defaults
			to false (this does not have to be satisfied).  
						 
			\pre t > 0.0.*/
			explicit MCMCGRDiagnosticPSRFLeaves(cxsc::real t,
												bool req = false);
			
			~MCMCGRDiagnosticPSRFLeaves();
			
			std::string getGRWorkingCalcsFilename() const;

			std::string getGRDiagnosticsFilename() const;
			
			std::string getScalarsName() const;
			
		private:
		
			MCMCGRDiagnosticPSRFLeaves();
			
			MCMCGRDiagnosticPSRFLeaves(const MCMCGRDiagnosticPSRFLeaves& other);
			
			MCMCGRDiagnosticPSRFLeaves operator=(const MCMCGRDiagnosticPSRFLeaves& rhs);
			MCMCGRDiagnosticPSRFLeaves operator=(MCMCGRDiagnosticPSRFLeaves rhs);
			
			real getScalarValue(
				const ChangeOfStateInformationAutoMCMC& info) const;
				
			std::string getBaseScalarsColName() const;
				
			static const std::string GRScalarFilename;
			static const std::string GRWorkingCalcsFilename;
			
			static const std::string scalarsName;
			static const std::string baseScalarsColName;

			
	};
} // end subpavings
#endif 

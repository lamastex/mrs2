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
using cherries as the scalar value and interval method for diagnostic.
 */

#ifndef __MCMCGRAUTO_DIAG_INTERVAL_CHERRIES_HPP__
#define __MCMCGRAUTO_DIAG_INTERVAL_CHERRIES_HPP__

#include "MCMCGRDiagnosticInterval.hpp"

namespace subpavings {
	
	class MCMCGRDiagnosticIntervalCherries : public MCMCGRDiagnosticInterval {
		
		public:

			/*! \brief Constructor.
			
			The largest number of states used for the calculation of
			each within chain interval will be 1,000,000 or the 
			next integer multiple of \si above that figure.   
			
			\param t The tolerance to use for the 
			Gelman-Rubin convergence criteria. 
			\param si the size of the sampling interval to use for 
			recalculating the diagnostic, eg \a si = 1 will recalculate
			for every change of state, \a si = 1000 will recalculate
			every 1000 changes of state, etc.
			\param percent the percentage interval to calculate.
			\param req flag for whether this
			diagnostic must be satisfied for burn-in. Defaults
			to false (this does not have to be satisfied).  
						 
			\pre t > 0.0, 0.0 < percent < 1.0.*/
			MCMCGRDiagnosticIntervalCherries(cxsc::real t,
											size_t si,
											double percent,
											bool req = false);
			
			/*! \brief Constructor.
			 
			\param t The tolerance to use for the 
			Gelman-Rubin convergence criteria. 
			\param si the size of the sampling interval to use for 
			recalculating the diagnostic, eg \a si = 1 will recalculate
			for every change of state, \a si = 1000 will recalculate
			every 1000 changes of state, etc.
			\param sm the maximum number of states to be used for
			any one within-chain interval.  If this is not an integer
			multiple of \a si then the maximum number of states used in
			the calculations will be \a sm rounded up to the next
			integer multiple of \a si greater than \a sm.
			\param percent the percentage interval to calculate.
			\param req flag for whether this
			diagnostic must be satisfied for burn-in.  
						 
			\pre t > 0.0, 0.0 < percent < 1.0.*/
			MCMCGRDiagnosticIntervalCherries(cxsc::real t,
								size_t si,
								size_t sm,
								double percent,
								bool req = false);
			
			~MCMCGRDiagnosticIntervalCherries();
			
			std::string getGRWorkingCalcsFilename() const;

			std::string getGRDiagnosticsFilename() const;
			
			std::string getScalarsName() const;
			
		private:
		
			MCMCGRDiagnosticIntervalCherries();
			
			MCMCGRDiagnosticIntervalCherries(const MCMCGRDiagnosticIntervalCherries& other);
			
			MCMCGRDiagnosticIntervalCherries operator=(const MCMCGRDiagnosticIntervalCherries& rhs);
			MCMCGRDiagnosticIntervalCherries operator=(MCMCGRDiagnosticIntervalCherries rhs);
			
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

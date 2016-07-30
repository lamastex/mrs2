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
\brief Declarations for a type to calculate Gelman-Rubin 
potential scale reduction factor (PSRF) heuristic 
diagnostics.
 
See Brooks, Stephen P, and Gelman, Andrew (1998), 
'General methods for monitoring convergence of iterative simulations',
Journal of Computational and Graphical Statistics, vol. 7, no. 4, pp. 434--455
 */

#ifndef __MCMCGR_DIAG_PSRF_HPP__
#define __MCMCGR_DIAG_PSRF_HPP__


#include "MCMCGRAuto.hpp"  


#include <vector>

namespace subpavings {
	
	/*! \brief Constructor.
			 
	Type to do Gelman-Rubin potential scale reduction factor (PSRF)
	diagnostic calculations.
	 
	Calculations use the (approximate) last half of the scalar values
	from each chain only.  If the number of states in the chain is n,
	the number used for the calculations is n-(n/2) where the division 
	is integer division, ie if n is odd we effectively use the last 
	(n+1)/2 values and if n is even we use the last n/2 values.  */
	class MCMCGRDiagnosticPSRF : public MCMCGRAuto::Diagnostic {
		
		public:

			/*! \brief Constructor.
			 
			\param t The tolerance to use for the 
			Gelman-Rubin convergence criteria. 
			\param req flag for whether this
			diagnostic must be satisfied for burn-in.  
						 
			\pre t > 0.0.*/
			explicit MCMCGRDiagnosticPSRF(cxsc::real t,
											bool req);
			
			/*! \brief Destructor.*/
			virtual ~MCMCGRDiagnosticPSRF();
			
			/*! \brief Get whether this diagnostic
			is required for burn-in. */
			bool isRequired() const;
			
			/*! \brief Clean old calculations out.*/
			virtual void clean();
			
			/*! \brief Clean old calculations out and 
			 * initialise for new calculations.
			 * 
			 * \param ch The number of chains to initialise for.
			 * \param ch Suggested number of loops to reserve
			 * space for: this can choose to reserve more 
			 * but should not reserve less capacity than this.
			 * \param logs Indicator for whether this is to keep
			 * additional logs for calculations.*/
			virtual void initialise(size_t ch, 
									size_t res,
									bool logs = false);
			
			/*! \brief Initialise values for each chain.
			 * 
			 * \param infoObjs The collection of 
			 * ChangeOfStateInformationAutoMCMC objects that
			 * this can ask to initialise its scalar values 
			 * for each chain. 
			 * \pre infoObjs.size() == chains .*/
			virtual void initialiseChainValues(
				const std::vector < 
				ChangeOfStateInformationAutoMCMC >& infoObjs);
			
			/*! \brief Update calculations for each chain.
			 * 
			 * \param infoObjs The collection of 
			 * ChangeOfStateInformationAutoMCMC objects that
			 * this can ask for updates to its scalar values 
			 * for each chain. 
			 * \pre infoObjs.size() == chains .*/
			virtual void updateChainValuesInLoop(
				const std::vector < 
				ChangeOfStateInformationAutoMCMC>& infoObjs);
			
			/*! \brief Do convergence diagnostic calculations
			 * for a loop or iteration of the MCMC process, 
			 * ie using current values assuming all chains
			 * have been iterated through for this loop.
			 * 
			 \return 1 if convergence criteria satisfied, 0 otherwise*/
			virtual int calcDiagnosticsForLoop();
			
			/*! \brief Output results including the contents
			 * of \a sampledIndPtr.
			 * 
			 * \note Results go to filename
			 * returned by getGRDiagnosticsFilename().
			 *
			 * Outputs the results of the Gelman-Rubin
			 * convergence diagnostic calculations (W, B, estimated
			 * variance, estimated R, whether convergence
			 * criteria satisfied, and \a sampledIndPtr values).
			 * 
			 * \note This version used for debugging.
			 * 
			 * \param sampledIndPtr A pointer to a collection 
			 * of values, each one corresponding to a state
			 * for which this holds calcuations, indicating
			 *  whether sampling took place in that state. 
			 * \param precData The precision for the output.*/
			void outputResults(
					const std::string& filenameGRScalars,
					const RealVec& sampledInd,
					int precData = 5) const;
			
			/*! \brief Output results.
			 * 
			 * \note Results go to filename
			 * returned by getGRDiagnosticsFilename().
			 * 
			 * Outputs the results of the Gelman-Rubin
			 * convergence diagnostic calculations (W, B, estimated
			 * variance, estimated R, whether convergence
			 * criteria satisfied).
			 
			 * \param precData The precision for the output.*/
			void outputResults(
					const std::string& filenameGRScalars,
					int precData = 5) const;
			
			/*! \brief Output values used in the calculations.
			 * 
			 * \note Used for debugging.  Results go to filename
			 * returned by getGRWorkingCalcsFilename().
			 * 
			 * Outputs the intermediary values held and 
			 * calculations for the the Gelman-Rubin
			 * convergence diagnostic calculations.
			 * \param precData The precision for the output.*/
			void outputCalculations(
					const std::string& filenameGRWorkingCalcs,
					int precData = 5) const;
			
			/*! \brief Get a const reference to the scalars held by this.*/
			const std::vector < std::vector < cxsc::real > >& 
							getScalarsRef() const;
			
			/*! \brief Get the current estimated R held.
			 * 
			 * \note Used for debugging.
			 * 					 * 
			  * \return The current estimated R value.
			 * \pre This holds at least one estimated R.*/
			cxsc::real getCurrentRhatValue() const;

			/*! \brief Get the name used by this in output files 
			 * for working calculations.*/
			virtual std::string getGRWorkingCalcsFilename() const = 0;

			/*! \brief Get the name used by this in results files.*/
			virtual std::string getGRDiagnosticsFilename() const = 0;
			
			/*! \brief Get the name used by this for the scalar 
			 * type it is using.*/
			virtual std::string getScalarsName() const = 0;
			
			/*! \brief Get the tolerance used by this.*/
			virtual cxsc::real getTolerance() const;
			
			
		private:
							
			MCMCGRDiagnosticPSRF();
			
			/*! \internal Doing it this way means that the concrete MCMCGRDiagnosticPSRF
			 * types have got to 'know about' the 
			 * ChangeOfStateInformationAutoMCMC type,
			 * ie know what method of that object to 
			 * call to get their own scalar value.  This is not ideal
			 * but the alternative is a proper implementation of 
			 * the visitor pattern and that the 
			 * ChangeOfStateInformationAutoMCMC has to
			 * then know about each concrete MCMCGRDiagnosticPSRF type: seemed
			 * easier this way around but no doubt I will regret it.... */
			virtual cxsc::real getScalarValue(
				const ChangeOfStateInformationAutoMCMC& info) const = 0;
			
			std::vector <std::string >& getScalarColNames(
				std::vector <std::string >& colNames) const;
			
			virtual std::string getBaseScalarsColName() const = 0;
			
			const cxsc::real tol;
			
			const bool required;
			
			bool keepLogs;
			
			size_t chains;
			
			int rhatDiagnosticFlag;
			
			size_t states;
			
			size_t statesNotInCalcs;
			
			
			std::vector<std::string> runningSumColNames;
			std::vector<std::string> sampleVarianceColNames;
			
			
			static const std::string baseRunningSumColName;
			static const std::string baseSampleVarianceColName;
			static const std::string overallRunningSumColName;
			
			/* vector containing vector of scalars for each chain
			we can work out the average for each chain so far from this
			start with a running sum of 0.0 for each chain */
			std::vector < RealVec > scalars;
			
			/* vector containing one running sum of diagnostic for each chain
			we can work out the average for each chain so far from this
			start with a running sum of 0.0 for each chain */
			RealVec runningSum;

			/* vector containing one running sum of 
			squared diagnostic for each chain
			we can work out the average of the squared v's ie v^2
			for each chain so far from this
			start with a running sum of 0.0 for each chain.
			(Use a dotprecision for each running sum to keep accuracy 
			when accumulating products of reals) */
			VecDotPrec runningSumSquared;

			/* value of running sum of diagnostic over all chains
			we can work out the average v over all chains so far from this */
			real runningSumAllChains;
			
			// only use ifdef MYDEBUG_LOGS
				// keep a vector of all the overall running sums as well
				RealVec runningSumOverall;
				// keep a vector of the runningsums for each chain as well
				std::vector < RealVec > runningSumChains;
				// keep a vector of the sample variances for each chain as well
				std::vector < RealVec > sampleVariances;
				/* keep a vector of the flag for convergence
				 * (it's not a real, but easier to output it if we treat it like one) */
				RealVec rhatFlag;
			
			RealVec Ws; // to hold the Ws
			RealVec Bs; // to hold the Bs
			RealVec estVarV; // to hold the estimated var(v)
			RealVec rhat; // to hold the rhats

			/* need to accumulate sum over all chains of the square of 
			* the running sum of vs 
			* for each chain for this starting state */
			cxsc::real initialSumOfSquaresOfRunningSums;
			
			// variables used in loops
			
			/* we want to accumulate the sample variance of the scalar summary leaves
			 * for each chain up to the point reached in this loop */
			cxsc::real sumOfSampleVariancesOverChains;
			
			/* also accumulate sum over all chains of the square of 
			 * the running sum of vs 
			 * for each chain up to the point reached in this loop */
			cxsc::real sumOfSquaresOfRunningSums;
	
	}; // end of class



} // end subpavings
#endif 

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
interval heuristic diagnostics.
 
See Brooks, Stephen P, and Gelman, Andrew (1998), 
'General methods for monitoring convergence of iterative simulations',
Journal of Computational and Graphical Statistics, vol. 7, no. 4, pp. 434--455.
 */

#ifndef __MCMCGR_DIAG_INTERVAL_HPP__
#define __MCMCGR_DIAG_INTERVAL_HPP__


#include "MCMCGRAuto.hpp"  


#include <vector>
#include <set>

namespace subpavings {
	
	/*! \brief Constructor.
			 
	Type to do Gelman-Rubin potential scale reduction factor (PSRF)
	diagnostic calculations.
	 
	Calculations use the (approximate) last half of the scalar values
	from each chain only.  If the number of states in the chain is n,
	the number used for the calculations is n-(n/2) where the division 
	is integer division, ie if n is odd we effectively use the last 
	(n+1)/2 values and if n is even we use the last n/2 values.  */
	class MCMCGRDiagnosticInterval : public MCMCGRAuto::Diagnostic {
		
		public:

			/*! \brief Constructor.
			 * 
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
			diagnostic must be satisfied for burn-in.  
						 
			\pre t > 0.0, 0.0 < percent < 1.0.*/
			MCMCGRDiagnosticInterval(cxsc::real t,
								size_t si,
								double percent,
								bool req);
			
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
			MCMCGRDiagnosticInterval(cxsc::real t,
								size_t si,
								size_t sm,
								double percent,
								bool req);
			
			/*! \brief Destructor.*/
			virtual ~MCMCGRDiagnosticInterval();
			
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
							
			MCMCGRDiagnosticInterval();
			
			/*! \internal
			 * Doing it this way means that the concrete MCMCGRDiagnosticInterval
			 * types have got to 'know about' the 
			 * ChangeOfStateInformationAutoMCMC type,
			 * ie know what method of that object to 
			 * call to get their own scalar value.  This is not ideal
			 * but the alternative is a proper implementation of 
			 * the visitor pattern and that the 
			 * ChangeOfStateInformationAutoMCMC has to
			 * then know about each concrete MCMCGRDiagnosticInterval type: seemed
			 * easier this way around but no doubt I will regret it.... */
			virtual cxsc::real getScalarValue(
				const ChangeOfStateInformationAutoMCMC& info) const = 0;
			
			
			std::vector <std::string >& getScalarColNames(
				std::vector <std::string >& colNames) const;
				
			virtual std::string getBaseScalarsColName() const = 0;
			
			real getIntervalLength(
						const std::multiset < cxsc::real >& set) const;
			
			const cxsc::real tol;
			
			const size_t samplingInterval;
			
			const double p;
			
			const bool required;
			
			bool keepLogs;
			
			size_t chains;
			
			int rhatDiagnosticFlag;
			
			size_t states;
			
			size_t statesNotInCalcs;
			
			/* maximum states to use for any one within-chain interval - 
			 * need this for larger problems or states used gets too 
			 * big and calculations take too long, and states 
			 * for total chain which is union of states for each
			 * within chain interval gets way way way too large. */
			size_t maxStatesForCalcs;
			
			std::vector<std::string> intervalLengthsColNames;
			static const std::string baseIntervalLengthsColName;
			static const std::string overallIntervalLengthsColName;
			static const std::string interchainIntervalStatisticColName;
			
			/* vector containing vector of scalars for each chain
			we can work out the average for each chain so far from this
			start with a running sum of 0.0 for each chain */
			std::vector < RealVec > scalars;
			
			/* multiset of scalars for all chains */
			std::multiset < cxsc::real > overallSets;
			
			/* current lengths of intervals one for each chain*/
			std::vector < cxsc::real >  currentIntervalLengths;
			
			/* Vector of over-individual-chain summary stat (eg mean, whatever*/
			std::vector < cxsc::real > interchainIntervalStatistic;
		
			/* Vector of total sequence interval lengths*/
			std::vector < cxsc::real > overallIntervalLengths;
			
			/* vector containing states actually used for calculations */
			std::vector < size_t > calculationStates;
			
			RealVec rhat; // to hold the rhats
			
			// only use ifdef MYDEBUG_LOGS
				/* vector containing vector of all interval lengths for each chain */
				std::vector < std::vector < cxsc::real > > intervalLengths;
			
				RealVec rhatFlag;
			
	}; // end of class



} // end subpavings
#endif 

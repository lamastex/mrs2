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
\brief Declarations for a type to do MCMC with some Gelman-Rubin 
* heuristic automatic sampling rule.
 */

#ifndef __MCMCGRAUTO_HPP__
#define __MCMCGRAUTO_HPP__


#include "adaptivehistogram.hpp"  

#include "piecewise_constant_function.hpp"

#include "adaptivehistogram_changeofstateinfo_automcmc.hpp"

#include <vector>
namespace subpavings {
	
	/*! \brief Type for objects capable of doing MCMC with some
	 * automatic convergence/burn-in testing rule using Gelman-Rubin
	 * heuristics.*/ 
	class MCMCGRAuto {
		
		public:

			/*! \brief An enum to control 
			the behaviour of the sampler if the criteria for the 
			Gelman-Rubin convergence diagnostic
			ceases to be met when it has previously been met.   
			
			<ul>
			<li>NORMAL: just burn in once and then sample from there, even if we
			 move out of the tolerance band afterwards.</li>
			 <li>DISCRIMINATING: sample only when we are in the tolerance band, but if we 
			 move out, keep all samples and just wait to come back in again, then
			 sample again.</li>
			 <li>CAUTIOUS: throw away all existing samples if we move out of the tolerance
			 band, wait to come back in, and sample from scratch again.</li>
			 </ul>*/
			typedef enum{NORMAL = -1, DISCRIMINATING = 0, CAUTIOUS = 1}
							SAMPLING_TYPE;
			
			/*! @name Get a PiecewiseConstantFunction from MCMC with some
			 * automatic convergence/burnin testing rule using Gelman-Rubin
			 * heuristics.
			 
			 The %PiecewiseConstantFunction returned is the average over
			 \a samplesNeeded samples from the chains started at the
			 histograms in \a histPtrsVec.
			
			Samples are collected from all the chains: at each sampling
			point (as determined by the automatically determined
			burnin point and the \a thinout) each chain is sampled
			in order, until \a samplesNeeded samples have been collected.
			* 
			Throws a std::invalid_argument exception if the number of 
			histograms in \a histPtrsVec is less than 2.
			 
			Throws a std::logic_error if any of the histograms in 
			\a histPtrsVec has an illegal state with respect to minPoints, ie
			if any of the cherries should not have been split.
			
			Throws a std::runtime_error if the process fails.  This can
			happen if there there is no legal change of 
			state available to a histogram in \a histPtrsVec (eg starting
			state is a single leaf root node which cannot be split)
			or if the required number of samples is not collected
			before the process is aborted because the maximum number
			of loops \a maxLoops has been reached.
			
			Each use of this to doMCMCAuto should initiate an MCMC
			process guided by a unique sequence of random numbers
			but also so that the first use of this will always be associated
			with the same sequence, the second use with its own sequence,
			etc etc, irrespective of how many prngs are generated from
			each separate sequence (ie how many loops are called for
			in each use of doMCMCAuto).    
			
			\internal
			This can have its own prng, which is used to generate seeds
			for a prng to be used in each separate call on this to
			doMCMCAuto.  Thus each call will have a different MCMC
			process but the overall process over several calls to
			doMCMCAuto is replicable and the starting points of the 
			sequence of prngs used for each call are independent of
			the number of loops in each call.
						 
			\param histPtrsVec is a collection of pointers to histograms 
			from which to start the chains (there should be at least 2).
			\param maxLoops is the maximum number of loops of the MCMC
			process to do before failing (if the required samples 
			\a samplesNeeded has not been collected when every 
			histogram in \a histPtrsVec has gone through \a maxLoops changes 
			of state then a runtime_error will be thrown).
			\param samplesNeeded is the number of samples needed in total. 
			\param thinout determines how frequently samples are taken
			once the convergence has been achieved. \a thinout = 1 means
			every state is sampled, \a thinout = 2 means every other state is
			sampled, etc.  
			\param samplesMonitored the number of post-burnin samples 
			monitored using 
			the convergence diagnostics if \a samplingType = NORMAL
			(this parameter has no effect if \a samplingType != NORMAL). 
			\param minPoints is the minimum number of points to allow 
			in a leaf node of a histogram state in the chain if the 
			node has any points in at all. This parameter affects what
			is considered  to be a "splittable node" in a histogram.
			\param rhatFlagCounterThreshold is the number of diagnostic
			criteria that must be satisfied for convergence to be 
			deemed to have taken place.
			\param runID is a integer identifier for the run, or a 
			run number.
			\param scalarsFilename is a string to use to find the filenames
			to output traces leaves, logposteriors, etc for logging
			the MCMC process. If scalarsFilename = "" (empty string) 
			no traces or logs will be output.
			\param samplingType is a value controlling 
			the behaviour if the criteria for the 
			Gelman-Rubin convergence diagnostic
			ceases to be met when it has previously been met.
			\param logging controls logging of samples from the chain.
			Suggested value LOGSAMPLES.
			\param loggingInStateChanges controls logging of the 
			evolution of the chains, state by state. Suggested value
			NOLOG.  
			\return A PiecewiseConstantFunction that is the average
			over all the samples collected of the states sampled
			from the histograms (this %PiecewiseConstantFunction
			will be normalised, ie have integral 1.0).
			\pre Each histogram in histPtrsVec has a legal state with 
			respect to \a minPoints and \a minPoints also permits
			at least one change in state in every histogram in \a histPtrsVec, 
			and there are at least two histograms in \a histPtrsVec. 
			\post Each histogram in \a histPtrsVec is in the state reached 
			when the MCMC process was terminated.  */ 
			//@{
				
				
									
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID, LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType = NORMAL) const; 
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int samplesMonitored,
							int thinout, 
							size_t minPoints,
							int runID, LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType = NORMAL) const; 
		
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID, 
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout,
							int samplesMonitored, 
							size_t minPoints,
							int runID, 
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const; //gat41
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int samplesMonitored, 
							int thinout, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName, 
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType) const; 
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& h,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID, LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType = NORMAL) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& h,
							int maxLoops, 
							int samplesNeeded,
							int thinout,
							int samplesMonitored,  
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID, LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType = NORMAL) const; 
		
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID, 
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout,
							int samplesMonitored,   
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID, 
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							real minVol,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,   
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName, 
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType = NORMAL) const; 
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			virtual PiecewiseConstantFunction doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,  
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const; //gat41
			
			//@}
			
			/*! \brief Get the maximum flag threshold that makes sense
			for this.*/				
			virtual size_t getMaxFlagThreshold() const = 0;
			
			
			/*! \brief Get the time taken, in seconds, for the 
			last successful run.  
			* 
			If the last run was not successful, the value returned by
			this method will not necessarily record the time taken 
			for that last run.*/
			virtual double getLastRunTime() const = 0;
			
			/*! \brief Get the state at which burnin is first reached.  
			 
			If the sampling type is CAUTIOUS, this will be the state
			at which the last successful attempt at burnin was reached.*/
			virtual size_t getFirstBurninTime() const = 0;
			
			
			virtual ~MCMCGRAuto(){}
			
			// inner class
			/*! \brief An abstract type to calculate Gelman-Rubin 
			 * convergence diagnostics for an MCMC process.
			 * 
			 * Each concrete type must be able to calculate its
			 * \f$ \widehat{R} \f$ (from within-chain and between
			 * chain variances) and report whether this is within
			 * a specified tolerance of 1.0.  Each concrete type 
			 * may also be configured as required for burn-in or not.
			 * 
			 * Required for burn-in is used by objects assessing
			 * burn-in over a variety of diagnostics.  Requried
			 * for burn-in indicates that even if other
			 * diagnostics have satisfactory \f$ \widehat{R} \f$
			 * values, if this does not then an MCMCGRAuto should
			 * not deem that burn-in has been achieved.

		    */
			class Diagnostic {
				
				public:

					/*! \brief Destructor.*/
					virtual ~Diagnostic() {};
					
					/*! \brief Get whether this diagnostic
					is required for burn-in. */
					virtual bool isRequired() const = 0;
					
					/*! \brief Clean old calculations out.*/
					virtual void clean() = 0;
					
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
											bool logs = false) = 0;
					
					/*! \brief Initialise values for each chain.
					 * 
					 * \param infoObjs The collection of 
					 * ChangeOfStateInformationAutoMCMC objects that
					 * this can ask to initialise its scalar values 
					 * for each chain. 
					 * \pre infoObjs.size() == chains .*/
					virtual void initialiseChainValues(
						const std::vector < 
						ChangeOfStateInformationAutoMCMC >& infoObjs) = 0;
					
					/*! \brief Update calculations for each chain.
					 * 
					 * \param infoObjs The collection of 
					 * ChangeOfStateInformationAutoMCMC objects that
					 * this can ask for updates to its scalar values 
					 * for each chain. 
					 * \pre infoObjs.size() == chains .*/
					virtual void updateChainValuesInLoop(
						const std::vector < 
						ChangeOfStateInformationAutoMCMC>& infoObjs) = 0;
					
					/*! \brief Do convergence diagnostic calculations
					 * for a loop or iteration of the MCMC process, 
					 * ie using current values assuming all chains
					 * have been iterated through for this loop.
					 * 
					 \return 1 if convergence criteria satisfied, 0 otherwise*/
					virtual int calcDiagnosticsForLoop() = 0;
					
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
					virtual void outputResults(
							const std::string& filenameGRScalars,
							const RealVec& sampledInd,
							int precData = 5) const = 0;
					
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
					virtual void outputResults(
							const std::string& filenameGRScalars,
							int precData = 5) const = 0;
					
					/*! \brief Output values used in the calculations.
					 * 
					 * \note Used for debugging.  Results go to filename
					 * returned by getGRWorkingCalcsFilename().
					 * 
					 * Outputs the intermediary values held and 
					 * calculations for the the Gelman-Rubin
					 * convergence diagnostic calculations.
					 * \param precData The precision for the output.*/
					virtual void outputCalculations(
							const std::string& filenameGRWorkingCalcs,
							int precData = 5) const = 0;
					
					/*! \brief Get a const reference to the scalars held by this.*/
					virtual const std::vector < std::vector < cxsc::real > >& 
									getScalarsRef() const = 0;
					
					/*! \brief Get the current estimated R held.
					 * 
					 * \note Used for debugging.
					 * 					 * 
					  * \return The current estimated R value.
					 * \pre This holds at least one estimated R.*/
					virtual cxsc::real getCurrentRhatValue() const = 0;

					/*! \brief Get the name used by this in output files 
					 * for working calculations.*/
					virtual std::string getGRWorkingCalcsFilename() const = 0;

					/*! \brief Get the name used by this in results files.*/
					virtual std::string getGRDiagnosticsFilename() const = 0;
					
					/*! \brief Get the name used by this for the scalar 
					 * type it is using.*/
					virtual std::string getScalarsName() const = 0;
					
					/*! \brief Get the tolerance used by this.*/
					virtual cxsc::real getTolerance() const = 0;
					
									
		}; // end of inner class


		private :
	
			virtual PiecewiseConstantFunction _doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,  
							size_t minPoints,
							real minVol,
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const = 0; //gat41
			
			virtual PiecewiseConstantFunction _doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,  
							size_t minPoints,
							real minVol,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const = 0; //gat41
			
			virtual int getDefaultSamplesMonitored() const = 0;
			
		
		
	};
} // end subpavings
#endif // __MCMCGRAUTO_HPP__

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
\brief Declarations for a new version of an object to do MCMC with Gelman-Rubin 
* heuristic automatic sampling rule using flexible criteria for convergence
* diagnostics.
 */

#ifndef __MCMCGRAUTO_NEW_HPP__
#define __MCMCGRAUTO_NEW_HPP__


#include "MCMCGRAuto.hpp"  

#include <gsl/gsl_rng.h>        // to know about the gsl random number generator

namespace subpavings {
	
	/*! \brief This class is instantiated with one or more
	MCMCGRAuto::Diagnostic diagnostics objects.  Each one 
	calculates its own Gelman-Rubin values and determines
	whether its own tolerance criteria is met. The MCMCGRAutoNew
	object organises the MCMC process and determines whether some
	overall burn-in criteria across all the diagnostics it is using
	is satisfied (eg, whether some threshold number of diagnostic 
	criteria \a rhatFlagCounterThreshold has been met).  
	If this overall burn-in criteria is 
	satisfied, this object samples from the chains.  When sampling
	is completed this object returns the sample average.   
	* 
	If rhatFlagCounterThreshold is not specified as an argument to
	to  doMCMCGRAuto(), the default value is the number of diagnostic 
	criteria, ie all have to be 'true' for convergence 
	to be deemed to have taken place.
	* 
	Objects of this type use a LogCatalan prior and a 
	UniformSSMProposal (stay-split-merge base chain) with 
	probability of staying in the same state
	\f$ \sigma = 1.0E-6 \f$ for the MCMC routine.*/ 
	class MCMCGRAutoNew : public MCMCGRAuto {
		
		public:

			/*! \brief Constructor.
			 
			The constructor uses \a seed to initialise a prng to generate
			random seeds to be used for the doMCMCGRAuto process.
			
			The number of post-burnin samples monitored using the convergence
			diagnostics if NORMAL sampling is used (i.e., burnin once)
			is set to some small default value.  This just means that 
			output from the diagnostics that can be used to make trace
			plots extends a little way past the point where
			burnin is diagnosed.  
						 
			\param diagObj The Gelman-Rubin convergence diagnostic object.  
			\param seed seed to use for a pseudo random number generator
			(defaults to 1234).
						 
			 \pre t > 0.0.*/
			explicit MCMCGRAutoNew(const std::vector < MCMCGRAuto::Diagnostic* >& diagPtrs,
								unsigned long int seed = 1234);
			
			~MCMCGRAutoNew();
			
			size_t getMaxFlagThreshold() const;
			
			double getLastRunTime() const;
			
			size_t getFirstBurninTime() const;

			
		private:
		
			
			MCMCGRAutoNew();
			
			MCMCGRAutoNew(const MCMCGRAutoNew& other);
			
			MCMCGRAutoNew operator=(const MCMCGRAutoNew& rhs);
			MCMCGRAutoNew operator=(MCMCGRAutoNew rhs);
			
			PiecewiseConstantFunction _doMCMCGRAuto(
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
							LogMCMCPrior& logPrior) const; //gat41
			
			PiecewiseConstantFunction _doMCMCGRAuto(
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
							LogMCMCPrior& logPrior) const; //gat41
			
			int getDefaultSamplesMonitored() const;
			
			void clean() const;
			
			static std::string getScalarFilename(const std::string& scalarType,
								const std::string& runIDstr,
								const std::string& scalarsFileName);
			
			static std::string getPath(const std::string& scalarsFileName);
			
			static void outputLeafLogPost(std::vector < Size_tVec >* const leavesPtr,
						std::vector < RealVec >* const logpostPtr,
						const std::vector < std::string>& leavesColNames,
						const std::vector < std::string>& logpostColNames,
						const std::string& GR_Logpost_Filename,
						int precData);
						
			static void outputLeafLogPost(std::vector < Size_tVec >* const leavesPtr,
						std::vector < RealVec >* const logpostPtr,
						const std::vector < std::string>& leavesColNames,
						const std::vector < std::string>& logpostColNames,
						const std::string& GR_Logpost_Filename,
						size_t startPos, 
						int precData);
			
			void outputLeafLogPost(
						const std::string& filenameLeavesAndLogpost,
						int precData) const;
						
			void outputLeafLogPost(
						const std::string& filenameLeavesAndLogpost,
						size_t startPos, 
						int precData) const;
			
			std::vector < MCMCGRAuto::Diagnostic* > diagObjPtrs;
			
			gsl_rng * rgsl;
			
			mutable double time;
			// keep track of when we (first) reached burnin
			mutable size_t burntinReachedState; 
	
			
			/* one vector of log posteriors for each chain
			 for which we need to keep loglikelihoods and log priors as well*/
			mutable std::vector < RealVec > logpost;  

			mutable RealVec sampledInd;

			mutable std::vector< std::string > leavesColNames;
			mutable std::vector< std::string > logpostColNames;
						
			static const int defaultSamplesMonitored;
			static const std::string  baseLeavesColName;
			static const std::string  baseLogpostColName;
			static const std::string baseSequenceStateFilename;
			
			static const std::string samplesLogFilename;
			
			static const size_t logInterval;
	
	
			
	
	};
} // end subpavings
#endif // __MCMCGRAUTO_HPP__

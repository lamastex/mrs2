/*
* Copyright (C) 2012 Jennifer Harlow
*
*/


/*! \file
\brief Common routines for KDE tests

 */

#ifndef __TEST_DEN_COMMON_HPP__
#define __TEST_DEN_COMMON_HPP__

//#include "mcmc_kde.hpp"
//#include "norm_kde.hpp"
//#include "type_kde.hpp"
#include "mixture_mvn.hpp"
#include "piecewise_constant_function.hpp"
//#include "histmcmcobjs.hpp"

#include "sptypes.hpp"

#include "cxsc.hpp"


#include <vector>
#include <string>

#include <cstddef>

/*
subpavings::PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					subpavings::LogMCMCPrior& logPrior);				
			
subpavings::PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const std::string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					subpavings::LogMCMCPrior& logPrior);

subpavings::PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded,
					subpavings::LogMCMCPrior& logPrior);

subpavings::PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const std::string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					subpavings::LogMCMCPrior& logPrior);
*/

std::vector < real >& getPCFDensitiesCensor(const subpavings::PiecewiseConstantFunction& pcf,
						subpavings::RVecData& intPts,
						std::vector < real >& pcfDensities, int dim);

std::vector < real >& getPCFDensities(const subpavings::PiecewiseConstantFunction& pcf,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities);
						
std::vector < real >& getPCFDensities(const subpavings::PiecewiseConstantFunction& pcf,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities,
						std::vector < double >& timing);

std::vector < real >& getPCFDensitiesCensor(const subpavings::PiecewiseConstantFunction& pcf,
						std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities);
						
std::vector < real >& getPCFDensitiesCensor(const subpavings::PiecewiseConstantFunction& pcf,
						std::vector < std::vector < real > >& intPts,
						std::vector < real >& pcfDensities,
						std::vector < double >& timing);
						

	
#endif

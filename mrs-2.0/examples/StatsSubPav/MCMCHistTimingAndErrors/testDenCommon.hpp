/*
* Copyright (C) 2012 Jennifer Harlow
*
*/


/*! \file
\brief Common routines for KDE tests

 */

#ifndef __TEST_DEN_COMMON_HPP__
#define __TEST_DEN_COMMON_HPP__

#include "mcmc_kde.hpp"
#include "norm_kde.hpp"
#include "type_kde.hpp"
#include "mixture_mvn.hpp"
#include "piecewise_constant_function.hpp"

#include "cxsc.hpp"


#include <vector>
#include <string>

#include <cstddef>


subpavings::PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded);				
			
subpavings::PiecewiseConstantFunction doHistMCMC(int rep, size_t n, size_t minPoints,	
					const std::string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed);

subpavings::PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed,
					unsigned int thinout,
					unsigned int samplesNeeded);

subpavings::PiecewiseConstantFunction doHistMCMC(const cxsc::ivector& box,
					int rep, size_t n, size_t minPoints,	
					const std::string& histFilenameBase,
					const std::vector < std::vector < double > >& simdata,
					std::vector < double >& timingMake,
					long unsigned int seed);


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
						
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		subpavings::kde::MixtureMVN* mixMVNptr,
		unsigned int thinout,
		unsigned int samplesNeeded);
		
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		subpavings::kde::MixtureMVN* mixMVNptr,
		const std::vector < std::vector < int > >& vecSliceDims,
		const std::vector < std::vector < double > >& vecSlicePts);
		
void doDenEst(
		size_t dim,
		long unsigned int seed,
		int reps,
		size_t minPoints,
		size_t n,
		size_t intN,
		const string& histFilenameBase,
		const string& logFilenameBase,	
		subpavings::kde::MixtureMVN* mixMVNptr,
		const std::vector < std::vector < int > >& vecSliceDims,
		const std::vector < std::vector < double > >& vecSlicePts,
		unsigned int thinout,
		unsigned int samplesNeeded);
	
#endif

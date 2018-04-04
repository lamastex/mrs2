/*
* Copyright (C) 2012 Jennifer Harlow
*
*/


#ifndef __TESTDENTOOLS_HPP__
#define __TESTDENTOOLS_HPP__

/* tools for comparing mcmc and kde

 */



#include "mixture_mvn.hpp"
#include "cxsc.hpp"

#include <vector>
#include <string>
#include <cstddef>

					
real avLogDen(const std::vector < real >& densities);

real avAbsDiffDen(const std::vector < real >& densities1,
				const std::vector < real >& densities2);

						
std::vector < std::vector < real > >& getQuasiRandomPoints(
			const cxsc::ivector& box,
			std::vector < std::vector < real > >& intPts,
			size_t N);


std::string makeDir(const std::string& baseOutputDir, 
				const std::string& thisDir);

void outputResults(const std::string& logFilename,
			const std::vector < std::vector < double > >& timingMake,
			const std::vector < std::vector < double > >& timingIntDensities,
			const std::vector < std::vector < real > >& avLogDens,
			const std::vector < std::vector < real > >& avLogDenRatios,
			const std::vector < std::vector < real > >& estL1ErrorsQR,
			const std::vector < std::vector < size_t > >& leaves,
			const std::vector < std::vector < double > >& topt);
			
void outputResults(const std::string& logFilename,
			const std::vector < std::vector < double > >& timings,
			const std::vector < std::vector < real > >& avLogDens,
			const std::vector < std::vector < real > >& avLogDenRatios,
			const std::vector < std::vector < real > >& estL1ErrorsQR,
			const std::vector < std::vector < size_t > >& leaves);

#endif

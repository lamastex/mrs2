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
\brief Definitions for a type to do Gelman-Rubin diagnostics for MCMC.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCGRDiagnosticPSRF.hpp"

#include "sptools.hpp"

#include <fstream>  // input and output streams
//#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions

#define MYDEBUG // extra console output for what is happening in process
//#define MYDEBUG_CALCS // extra console output for calculations
//#define MYDEBUG_CALCS_EXTRA // even more, with getchar() for step by step debugging
//#define NDEBUG // uncomment this to turn off assertion checking and all for this module only
#define OLDCALCMETHOD
#ifdef NDEBUG // ie only allow defines for the others if we have not defined NDEBUG for no debugging
	#undef MYDEBUG
	#undef MYDEBUG_CALCS
	#undef MYDEBUG_CALCS_EXTRA

#endif

using namespace cxsc;
using namespace subpavings;
using namespace std;


// statics
const std::string MCMCGRDiagnosticPSRF::baseRunningSumColName = "vsSum_";
const std::string MCMCGRDiagnosticPSRF::baseSampleVarianceColName = "vsVar_";
const std::string MCMCGRDiagnosticPSRF::overallRunningSumColName = "OverallVsSum";





MCMCGRDiagnosticPSRF::MCMCGRDiagnosticPSRF(cxsc::real t, bool req)
	: tol(t), required(req),
		keepLogs(false),
		chains(0),
		rhatDiagnosticFlag(0),
		states(0),
		statesNotInCalcs(0),
		runningSumAllChains(0.0),
		initialSumOfSquaresOfRunningSums(0.0)
{
	if (!(tol > 0.0)) throw std::invalid_argument(
			"MCMCGRDiagnosticPSRF::MCMCGRDiagnostic(...) : t");
}

MCMCGRDiagnosticPSRF::~MCMCGRDiagnosticPSRF()
{}

bool MCMCGRDiagnosticPSRF::isRequired() const
{
	return required;
}

void MCMCGRDiagnosticPSRF::clean()
{
	vector<string>().swap(runningSumColNames);
	vector<string>().swap(sampleVarianceColNames);
	
	keepLogs = false;
	
	chains = 0;
	
	rhatDiagnosticFlag = 0;
	
	states = 0;
	statesNotInCalcs = 0;
	
	std::vector < RealVec >().swap( scalars );
	RealVec().swap( runningSum );
	RealVec().swap( Ws );
	RealVec().swap( Bs); 
	RealVec().swap( estVarV ); 
	RealVec().swap( rhat ); 
	
	RealVec().swap( runningSumOverall );
	std::vector < RealVec >().swap( runningSumChains );
	std::vector < RealVec >().swap( sampleVariances );
	RealVec().swap( rhatFlag );
	
	VecDotPrec().swap(runningSumSquared);

	runningSumAllChains = 0.0;
	
	initialSumOfSquaresOfRunningSums = 0.0;
}

void MCMCGRDiagnosticPSRF::initialise(size_t ch, size_t res, bool logs)
{
	clean();
	
	keepLogs = logs;
	
	chains = ch;
	
	rhatDiagnosticFlag = 0;
	
	states = 1;
	
	statesNotInCalcs = 0;
	
	runningSumColNames = vector<string>(chains);
	sampleVarianceColNames = vector<string>(chains);
	
	scalars = std::vector <RealVec > (chains);
	for (size_t ci = 0; ci < chains; ++ci) {
		scalars[ci].reserve(res);
	}
	
	/* vector containing one running sum of v's for each chain
	we can work out the average v for each chain so far from this
	start with a running sum of 0.0 for each chain */
	runningSum = RealVec (chains, cxsc::real(0.0));

	if (keepLogs) {
		// keep a vector of all the overall running sums as well
		runningSumOverall = RealVec(1, cxsc::real(0.0));
		runningSumOverall.reserve(res);
		// keep a vector of the runningsums for each chain as well
		runningSumChains = std::vector < RealVec >(chains);
		// keep a vector of the sample variances for each chain as well
		sampleVariances = std::vector < RealVec >(
					chains, RealVec(1, cxsc::real(0.0)) );
		
		for (size_t ci = 0; ci < chains; ++ci) {
			runningSumChains[ci].reserve(res);
			sampleVariances[ci].reserve(res);
		}
		
		/* keep a vector of the flag for convergence
		 * (it's not a real, but easier to output it if we treat it like one) */
		rhatFlag = RealVec(1, cxsc::real(0.0));
		rhatFlag.reserve(res);
		
	}

	Ws = RealVec(1, cxsc::real (0.0) ); // to hold the Ws
	Ws.reserve(res);
	Bs = RealVec(1, cxsc::real (0.0) ); // to hold the Bs
	Bs.reserve(res);
	estVarV = RealVec(1, cxsc::real (0.0) ); // to hold the estimated var(v)
	estVarV.reserve(res);
	rhat = RealVec(1, cxsc::real (0.0) ); // to hold the rhats
	rhat.reserve(res);

	/* vector containing one running sum of 
	squared v's for each chain
	we can work out the average of the squared v's ie v^2
	for each chain so far from this
	start with a running sum of 0.0 for each chain.
	(Use a dotprecision for each running sum to keep accuracy 
	when accumulating products of reals) */
	runningSumSquared = VecDotPrec(chains, cxsc::dotprecision(0.0));

	/* 
	runningSumAllChains and
	initialSumOfSquaresOfRunningSums  are done in clean() */
	if (keepLogs) {
		for (size_t ci = 0; ci < chains; ci++) {
			{
				std::ostringstream stm;
				stm << baseRunningSumColName << ci;
				runningSumColNames[ci] = stm.str();
			}
			
			{
				std::ostringstream stm;
				stm << baseSampleVarianceColName << ci;
				sampleVarianceColNames[ci] = stm.str();
			}
			
		} // end loop initialising stuff for each chain
	}
	
}

void MCMCGRDiagnosticPSRF::initialiseChainValues(
				const std::vector < ChangeOfStateInformationAutoMCMC>& infoObjs)
{
	for (size_t ci = 0; ci < chains; ++ci) {
		//get the last value
		real lastValue = getScalarValue(infoObjs[ci]);
		scalars.at(ci).push_back(lastValue);
		
		// update the running sum of v's for the chain, held in runningSum
		cxsc::real newRunningSum = runningSum.at(ci) 
									+ lastValue;
		runningSum.at(ci) = newRunningSum;
				
		// accumulate the square of the running sum of v's 
		initialSumOfSquaresOfRunningSums += newRunningSum*newRunningSum;
				
		/* update the running sum of squared v's over this chain
		 *  held in runningSumSquared as a dot precision */
		cxsc::accumulate( runningSumSquared[ci], lastValue, lastValue );
		
		// update  the overall running sum runningSumAllChains 
		runningSumAllChains += lastValue;
		
		if (keepLogs) {
			//sampleVariances.at(ci) was initialised to 0.0
			runningSumChains.at(ci).push_back (newRunningSum);
			// store the current runningSumAllChains as well
			runningSumOverall.back() += newRunningSum;
		}
	}
}


/* the overall running sum runningSumAllChains 
 * was initialised to 0.0 
 * and if keeping logs, runningSumOverall was initialised to contain one 0.0 
 * and similarly rhatFlagPtr was initialised to contain one 0.0*/

/* and we started the convergence statistics for chains with just one state in
 * with one 0.0 in each (Ws, Bs, estVarsVs, rhats)
 * when we initialised */


void MCMCGRDiagnosticPSRF::updateChainValuesInLoop(
		const std::vector < ChangeOfStateInformationAutoMCMC>& infoObjs)
{
	// do initial values for everything so far
	/* we want to accumulate the sample variance of the scalar summary
	 * for each chain up to the point reached in this loop */
	sumOfSampleVariancesOverChains = 0.0;
	
	/* also accumulate sum over all chains of the square of 
	 * the running sum of v's 
	 * for each chain up to the point reached in this loop */
	sumOfSquaresOfRunningSums = 0.0;
	
	++states;
	
	size_t newStatesNotInCalcs = states/2;
	#ifdef OLDCALCMETHOD
		newStatesNotInCalcs = 0;
	#endif
	
	
	int dropped = newStatesNotInCalcs - statesNotInCalcs;
	assert((dropped == 0) || (dropped == 1));
	statesNotInCalcs = newStatesNotInCalcs;
	
	/*On the first run, states will be 2 but statesForCalcs will be 1)*/
	size_t statesForCalcs = states - statesNotInCalcs;
	#ifdef MYDEBUG_CALCS
		cout << "\nstates for calcs = " << statesForCalcs << endl;
		cout << "dropped = " << dropped << endl;
		cout << "statesNotInCalcs = " << statesNotInCalcs << endl;
	#endif	
	
	for (size_t ci = 0; ci < chains; ++ci) {
		//get the last value
		real lastValue = getScalarValue(infoObjs[ci]);
		real droppedValue = 0.0;
		if (dropped) droppedValue = scalars.at(ci)[statesNotInCalcs-1];
		#ifdef MYDEBUG_CALCS
			cout << "\nchain = " << ci << endl;
			cout << "droppedValue = " << droppedValue << endl;
		#endif
		scalars.at(ci).push_back(lastValue);
		
		// update the running sum of v's for the chain, held in runningSum
		cxsc::real newRunningSum = runningSum.at(ci) 
						+ lastValue - droppedValue;
		runningSum.at(ci) = newRunningSum;
		
		// accumulate the square of the running sum of v's 
		sumOfSquaresOfRunningSums += newRunningSum*newRunningSum;
		
		/* update the running sum of squared v's over this chain
		 *  held in runningSumSquared as a dot precision */
		cxsc::accumulate( runningSumSquared[ci], 
					lastValue, lastValue );
		cxsc::accumulate( runningSumSquared[ci], 
					-droppedValue, droppedValue );
					
		// update  the overall running sum runningSumAllChains 
		runningSumAllChains += (lastValue - droppedValue);
		
		/* accumulate the sample variance for v's for this chain: 
		 * sample variance for the scalar summary v
		 * calculated as (sum of squares - n * square of averages)/(n-1)
		 * which equals (sum of squares - square of sums/n)/(n-1) */
		
		cxsc::real thisSampleVariance = 0.0;
		if (statesForCalcs > 1) {
			thisSampleVariance =  ( 1.0/(statesForCalcs - 1) )
				*( cxsc::rnd(runningSumSquared[ci])
				-  (newRunningSum*newRunningSum/(statesForCalcs * 1.0)) ) ;
		}
		
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "statesForCalcs = " << statesForCalcs << endl;
			cout << "runningSumSquared[ci] = " << runningSumSquared[ci] << endl;
			cout << "newRunningSum = " << newRunningSum << endl;
			cout << "thisSampleVariance = " << thisSampleVariance << endl;
		#endif
		
		sumOfSampleVariancesOverChains += thisSampleVariance;
		
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "now, sumOfSampleVariancesOverChains = " << sumOfSampleVariancesOverChains << endl;
		#endif
		
		if (keepLogs) {
			sampleVariances.at(ci).push_back( thisSampleVariance );
			runningSumChains.at(ci).push_back (newRunningSum);
		}
		
		#ifdef MYDEBUG_CALCS
		
			//check thisSampleVariance is correct, doing it the long way
			// scalars[ci] has the v_ij for each chain i
			
			real acc(0.0);
			{
				RealVecItr it = scalars.at(ci).begin();
				advance(it, statesNotInCalcs);
				assert( distance(it, scalars.at(ci).end()) == statesForCalcs);
				for ( ; it < scalars.at(ci).end(); ++it) {
					acc+= (*it);
				}
			}
			
			cxsc::real av = acc/(1.0*statesForCalcs);
			cxsc::dotprecision accDiffs(0.0);
			
			{
				RealVecItr it = scalars.at(ci).begin();
				advance(it, statesNotInCalcs);
				
				for ( ; it < scalars.at(ci).end(); ++it) {
					cxsc::real thisDiff = (*it) - av;
					// sum up the squares of the differences compared to overall average
					cxsc::accumulate(accDiffs, thisDiff, thisDiff);
				}
			}
			cxsc::real altVar = 0.0;
			if (statesForCalcs > 1) altVar = rnd(accDiffs)/( statesForCalcs- 1.0 );
			
			cout << "\nthisSampleVariance is\t" 
					<< sampleVariances.at(ci).back() << endl;
			cout << "and value calculated from basics is \t" << altVar << endl;
	
		#endif
	}
}


int MCMCGRDiagnosticPSRF::calcDiagnosticsForLoop()
{
	if (keepLogs) {
		// store the current runningSumAllChains as well
		runningSumOverall.push_back(runningSumAllChains);
	}

	// convergence diagnostics calculations for v's

	// the Ws: average, over chains, of sample variance of scalar value
	cxsc::real thisW = sumOfSampleVariancesOverChains/(chains * 1.0); 
	
	#ifdef MYDEBUG_CALCS_EXTRA
			cout << "and thisW = " << thisW << endl;
	#endif
	
				
	Ws.push_back(thisW); 
	// the Bs
	size_t statesForCalcs = states - statesNotInCalcs;
	
	cxsc::real thisB = (1.0/( (chains - 1) * statesForCalcs ) 
						* ( sumOfSquaresOfRunningSums 
						- (runningSumAllChains 
						* runningSumAllChains/(chains * 1.0)) ) );
	Bs.push_back(thisB); 
	
	#ifdef MYDEBUG_CALCS
		//check thisB is correct, doing it the long way
		// runningSumPtr has one running sum for each chain
		RealVec chainAverages;
		cxsc::real accRunningSums(0.0);
		for (RealVecItr it = runningSum.begin(); it < runningSum.end(); ++it) {
			cxsc::real thisChainRunningSum = (*it);
			cxsc::real thisChainAv = thisChainRunningSum/(statesForCalcs * 1.0);
			chainAverages.push_back(thisChainAv);
			accRunningSums+=thisChainRunningSum;
		}
		cxsc::real overallAv = accRunningSums/(statesForCalcs * chains * 1.0);
		cxsc::dotprecision accDiffs(0.0);
		for (RealVecItr it = chainAverages.begin(); it < chainAverages.end(); ++it) {
			cxsc::real thisDiff = (*it) - overallAv;
			// sum up the squares of the differences compared to overall average
			cxsc::accumulate(accDiffs, thisDiff, thisDiff);
		}
		cxsc::real altB = rnd(accDiffs)*( statesForCalcs/(chains - 1.0) );
		
		cout << "\nthisB for v's is\t" << thisB << endl;
		cout << "altB for v's is\t" << altB << endl;
		//assert(thisB == altB);
	
	#endif
	
	// the estimated var(v)
	cxsc::real thisVarV(0.0);
	if (statesForCalcs > 1) {
		thisVarV = statesForCalcs/(statesForCalcs-1.0) 
					* thisW + (1.0/statesForCalcs)*thisB;
	}
	#ifndef OLDCALCMETHOD
		if (statesForCalcs > 1) thisVarV +=thisB/(1.0*chains*statesForCalcs);
	#endif
	
	estVarV.push_back(thisVarV); 
	// the rhats
	cxsc::real thisRhat(0.0);
	// allow division by 0 if w = 0 when var does not
	if (thisW > 0.0 || thisVarV > 0.0) {
		thisRhat = thisVarV/thisW;
	}
	rhat.push_back(thisRhat); 
	
	#ifdef MYDEBUG_CALCS_EXTRA
			cout << "thisRhat = " << thisRhat << " - press any key " << endl;
			getchar();
	#endif
	
	if ( (thisRhat <= 1.0 + tol) 
				&& (thisRhat >= 1.0 - tol) ) {
		// if we have not been converged before on this scalar value
		if (!rhatDiagnosticFlag)  {
			#ifdef MYDEBUG
				cout << "\n" << getScalarsName() 
					<< " convergence test satisfied in state " 
					  << states << " (states in calcs = " 
					  << statesForCalcs << ")"<< endl;
			#endif
			// set the flag for this scalar value
			rhatDiagnosticFlag = 1;
		}
	} 
	else { // not converged on this scalar value
		// if we were okay on this scalar value before
		if (rhatDiagnosticFlag) {
			#ifdef MYDEBUG
				cout << "\n--------- Note: " << getScalarsName() 
				<< " convergence test now NOT satisfied in state " 
				  << states << endl;
		
			#endif
			
			rhatDiagnosticFlag = 0; // update the flag
			
		} 
	}
	
	if (keepLogs) {
		// store the flag as well, as a real, which is a fudge...
		rhatFlag.push_back(rhatDiagnosticFlag);
	}

	// end of checking diagnostic for v's

	return rhatDiagnosticFlag;
} // end calculations


void MCMCGRDiagnosticPSRF::outputResults(
				const std::string& filenameGRScalars,
				const RealVec& sampledInd,
				int precData) const
{
	/* sampledInd may contain more values than we monitored */
	RealVec::const_iterator it;
	advance(it, scalars.size());
	RealVec sampledIndTmp(sampledInd.begin(), it);
	
	std::vector < std::string > colNames;
	getScalarColNames(colNames);
	colNames.push_back("W");
	colNames.push_back("B");
	colNames.push_back("estVarV");
	colNames.push_back("rhat");
	
		colNames.push_back("rhatFlag");
		colNames.push_back("sampled?");
	
	std::vector < const RealVec* > data;
	addDataPtrs(data, scalars);
	data.push_back(&Ws);
	data.push_back(&Bs);
	data.push_back(&estVarV);
	data.push_back(&rhat);
	
		data.push_back(&rhatFlag);
		data.push_back(&sampledIndTmp);
	
	outputToFileVertical(data, colNames, 
	filenameGRScalars, precData);
} 

void MCMCGRDiagnosticPSRF::outputResults(
			const std::string& filenameGRScalars,
			int precData) const
{
	std::vector <std::string > colNames;
	getScalarColNames(colNames);
	colNames.push_back("W");
	colNames.push_back("B");
	colNames.push_back("estVarV");
	colNames.push_back("rhat");
	
	std::vector < const RealVec* > data;
	addDataPtrs(data, scalars);
	data.push_back(&Ws);
	data.push_back(&Bs);
	data.push_back(&estVarV);
	data.push_back(&rhat);
	
	outputToFileVertical(data, colNames, 
	filenameGRScalars, precData);
} 




void MCMCGRDiagnosticPSRF::outputCalculations(
				const std::string& filenameGRWorkingCalcs,
				int precData) const
{
	
	std::vector <std::string > colNames;
	
	getScalarColNames(colNames);
	
	/* output working calcs: all v's for each chain, 
	 * running sums for each chain, sample variances,
	 * overall running sums */
	colNames.insert(colNames.end(), runningSumColNames.begin(), runningSumColNames.end());
	colNames.insert(colNames.end(), sampleVarianceColNames.begin(), sampleVarianceColNames.end());
	colNames.push_back(overallRunningSumColName);
	
	std::vector < const RealVec* > data;
	data = addDataPtrs(data, scalars);
	data = addDataPtrs(data, runningSumChains);
	data = addDataPtrs(data, sampleVariances);
	data.push_back(&runningSumOverall);
	
	outputToFileVertical(data, colNames, filenameGRWorkingCalcs, precData);
	
}



const std::vector < std::vector < cxsc::real > >& 
			MCMCGRDiagnosticPSRF::getScalarsRef() const
{
	return scalars;
}

cxsc::real MCMCGRDiagnosticPSRF::getCurrentRhatValue() const
{
	if (rhat.empty())
		throw std::runtime_error(
			"MCMCGRDiagnosticPSRF::getCurrentRhatValue(size_t): no value to get");
	return rhat.back();
}

cxsc::real MCMCGRDiagnosticPSRF::getTolerance() const
{
	return tol;
}

// private 

std::vector <std::string >& MCMCGRDiagnosticPSRF::getScalarColNames(
	std::vector <std::string >& colNames) const
{	
	/* and start column names for each scalar-chain */
	for (size_t i = 0; i < chains; ++i) {
		
		std::ostringstream stm;
		stm << getBaseScalarsColName() << i;
		colNames.push_back(stm.str());
	}
	return colNames;
}

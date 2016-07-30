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

#include "MCMCGRDiagnosticInterval.hpp"

#include "sptools.hpp"

#include <numeric> // accumulate
#include <iterator> 
#include <fstream>  // input and output streams
//#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions

#define MYDEBUG // extra console output for what is happening in process
//#define MYDEBUG_CALCS // extra console output for calculations
//#define MYDEBUG_CALCS_EXTRA // even more, with getchar() for step by step debugging
//#define NDEBUG // uncomment this to turn off assertion checking and all for this module only
//#define OLDCALCMETHOD // use all states for calculation of Rhat values, not just second half of sequence
#ifdef NDEBUG // ie only allow defines for the others if we have not defined NDEBUG for no debugging
	#undef MYDEBUG
	#undef MYDEBUG_CALCS
	#undef MYDEBUG_CALCS_EXTRA

#endif

using namespace cxsc;
using namespace subpavings;
using namespace std;


// statics
const std::string MCMCGRDiagnosticInterval::baseIntervalLengthsColName 
			= "Len_";
const std::string MCMCGRDiagnosticInterval::overallIntervalLengthsColName 
			= "OallLen";
const std::string MCMCGRDiagnosticInterval::interchainIntervalStatisticColName 
			= "LenStat";


MCMCGRDiagnosticInterval::MCMCGRDiagnosticInterval(cxsc::real t, 
													size_t si,
													double percent,
													bool req)
	: 	tol(t), samplingInterval(si),
		p((1-percent)/2), required(req),
		keepLogs(false),
		chains(0),
		rhatDiagnosticFlag(0),
		states(0),
		statesNotInCalcs(0),
		maxStatesForCalcs(1000000)
{
	if (!(tol > 0.0)) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : t");
	if (!(percent > 0.0) || !(percent < 1.0)) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : percent");
	if (si == 0) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : si");
	if (maxStatesForCalcs%si) {
		maxStatesForCalcs = (maxStatesForCalcs/si + 1) * si;
	}
}


MCMCGRDiagnosticInterval::MCMCGRDiagnosticInterval(cxsc::real t, 
													size_t si,
													size_t ms,
													double percent,
													bool req)
	: 	tol(t), samplingInterval(si), 
		p((1-percent)/2), required(req),
		keepLogs(false),
		chains(0),
		rhatDiagnosticFlag(0),
		states(0),
		statesNotInCalcs(0),
		maxStatesForCalcs(ms)
{
	if (!(tol > 0.0)) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : t");
	if (!(percent > 0.0) || !(percent < 1.0)) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : percent");
	if (si == 0) throw std::invalid_argument(
			"MCMCGRDiagnosticInterval::MCMCGRDiagnostic(...) : si");
	if (maxStatesForCalcs%si) {
		maxStatesForCalcs = (maxStatesForCalcs/si + 1) * si;
	}
}

MCMCGRDiagnosticInterval::~MCMCGRDiagnosticInterval()
{}

bool MCMCGRDiagnosticInterval::isRequired() const
{
	return required;
}

void MCMCGRDiagnosticInterval::clean()
{
	vector<string>().swap(intervalLengthsColNames);
	
	keepLogs = false;
	
	chains = 0;
	
	rhatDiagnosticFlag = 0;
	
	states = 0;
	statesNotInCalcs = 0;
	
	std::vector < RealVec >().swap( scalars );
	
	RealVec().swap( rhat ); 
	
	std::multiset < cxsc::real >().swap( overallSets );
	
	std::vector < cxsc::real >().swap(currentIntervalLengths);
	
	std::vector < std::vector < cxsc::real > >().swap( intervalLengths );
	
	std::vector < size_t >().swap( calculationStates );
	
	std::vector < cxsc::real >().swap( overallIntervalLengths );
	
	std::vector < cxsc::real >().swap( interchainIntervalStatistic );
	
	RealVec().swap( rhatFlag );
	
}

void MCMCGRDiagnosticInterval::initialise(size_t ch, size_t res, bool logs)
{
	clean();
	
	keepLogs = logs;
	
	chains = ch;
	
	rhatDiagnosticFlag = 0;
	
	states = 1;
	
	statesNotInCalcs = 0;
	
	scalars = std::vector <RealVec > (chains);
	for (size_t ci = 0; ci < chains; ++ci) {
		scalars[ci].reserve(res);
	}
	
	currentIntervalLengths = vector < real >(chains, 0.0);
	
	if (keepLogs) {
		
		intervalLengthsColNames = vector<string>(chains);
		intervalLengths = std::vector <RealVec > (chains);
		for (size_t ci = 0; ci < chains; ++ci) {
			intervalLengths[ci].reserve(res);
			intervalLengths[ci].push_back(0.0);
			{
				std::ostringstream stm;
				stm << baseIntervalLengthsColName << ci;
				intervalLengthsColNames[ci] = stm.str();
			}
		}
		
		/* keep a vector of the flag for convergence
		 * (it's not a real, but easier to output it if we treat it like one) */
		rhatFlag = RealVec(1, cxsc::real(0.0));
		rhatFlag.reserve(res);
		
	}
	
	
	rhat = RealVec(1, cxsc::real (0.0) ); // to hold the rhats
	rhat.reserve(res);
	
	overallIntervalLengths.reserve(res);
	overallIntervalLengths.push_back(0);
	interchainIntervalStatistic.reserve(res);
	interchainIntervalStatistic.push_back(0);
	calculationStates.reserve(res);
	calculationStates.push_back(0);
	
	

}

void MCMCGRDiagnosticInterval::initialiseChainValues(
				const std::vector < ChangeOfStateInformationAutoMCMC>& infoObjs)
{
	for (size_t ci = 0; ci < chains; ++ci) {
		//get the last value
		real lastValue = getScalarValue(infoObjs[ci]);
		scalars.at(ci).push_back(lastValue);

	}
}


/* the overall running sum runningSumAllChains 
 * was initialised to 0.0 
 * and if keeping logs, runningSumOverall was initialised to contain one 0.0 
 * and similarly rhatFlagPtr was initialised to contain one 0.0*/

/* and we started the convergence statistics for chains with just one state in
 * with one 0.0 in each (Ws, Bs, estVarsVs, rhats)
 * when we initialised */


void MCMCGRDiagnosticInterval::updateChainValuesInLoop(
		const std::vector < ChangeOfStateInformationAutoMCMC>& infoObjs)
{
	++states; // states is effectively number of transitions + 1
	statesNotInCalcs = states/2;
	/* but limit states used for calcs */
	if ( (states - statesNotInCalcs) > maxStatesForCalcs ) {
		statesNotInCalcs = states - maxStatesForCalcs;
	}
	#ifdef OLDCALCMETHOD
		statesNotInCalcs = 0;
	#endif
		
	/*On the first run, states will be 2 but statesForCalcs will be 1)*/
	size_t statesForCalcs = states - statesNotInCalcs;
	
	//If we are updating
	if ( (states % samplingInterval == 0) 
			&& (statesForCalcs % samplingInterval == 0) ) {
		
		#ifdef MYDEBUG_CALCS
			cout << "\nsampling" << endl;
			cout << "\nstates for calcs = " << statesForCalcs << endl;
			cout << "statesNotInCalcs = " << statesNotInCalcs << endl;
		#endif	
		
		vector < real >(chains, 0.0).swap(currentIntervalLengths);
		
		//the overallSets container should be clear
		assert(overallSets.empty()); 
		
		for (size_t ci = 0; ci < chains; ++ci) {
			
			//get the last value
			real lastValue = getScalarValue(infoObjs[ci]);
			scalars.at(ci).push_back(lastValue);
			
			/* for each chain put the values to be used into the set*/
			std::vector < real >::iterator it = scalars[ci].begin();
			advance(it, statesNotInCalcs);
			// a new set from the scalars
			std::multiset < cxsc::real > thisSet( it, scalars[ci].end() );
			// and also add to the overall set
			overallSets.insert( it, scalars[ci].end() );
			
			#ifdef MYDEBUG_CALCS_EXTRA
				cout << "\n chain = " << ci << endl;
			#endif
			//get the interval length
			real len(0.0);
			if (statesForCalcs > 1) len = getIntervalLength(thisSet);
			currentIntervalLengths[ci] = len;
			
			if (keepLogs) intervalLengths[ci].push_back(len);
			
			#ifdef MYDEBUG_CALCS_EXTRA
				cout << "len = " << len << endl;
			#endif
			
		}
	}
	else { // not updating, get the scalar values only
		for (size_t ci = 0; ci < chains; ++ci) {
			
			//get the last value
			real lastValue = getScalarValue(infoObjs[ci]);
			scalars.at(ci).push_back(lastValue);
			
			// and put old values into log files
			if (keepLogs) intervalLengths[ci].push_back(
					currentIntervalLengths[ci]);
			
		}
	}
}


int MCMCGRDiagnosticInterval::calcDiagnosticsForLoop()
{
	//if we are updating
	if ( (states % samplingInterval == 0) 
			&& ((states - statesNotInCalcs) % samplingInterval == 0) ) {
		calculationStates.push_back(states);
		//get the interval length for the overall set
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "\noverall set = "<< endl;
		#endif
		real len(0.0);
		if (states - statesNotInCalcs > 1) len = getIntervalLength(overallSets);
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "overall len = " << len << endl;
		#endif
		
		//clear out the overallSets container
		std::multiset < cxsc::real >().swap( overallSets ); 
		
		
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "chain interval lens are:" << endl;
			ostream_iterator<real> out_it (cout,"\t");
			copy ( currentIntervalLengths.begin(),currentIntervalLengths.end(), out_it );
			cout << endl;
		#endif
		
		real meanChainIntervalLength(0.0);
		meanChainIntervalLength = std::accumulate(currentIntervalLengths.begin(),
		currentIntervalLengths.end(), meanChainIntervalLength)/(1.0*chains);
		
		/* keep the diagnostics */
		overallIntervalLengths.push_back(len);
		interchainIntervalStatistic.push_back(meanChainIntervalLength);
		
		real thisRhat(0.0);
		if (meanChainIntervalLength > 0) thisRhat = len/meanChainIntervalLength;
		
		rhat.push_back(thisRhat); 
		
		#ifdef MYDEBUG_CALCS_EXTRA
			cout << "meanChainIntervalLength = " << meanChainIntervalLength << endl;
			cout << "thisRhat = " << thisRhat << " - press any key " << endl;
			getchar();
		#endif
	
		if ( !((thisRhat > 1.0 + tol) 
					|| (thisRhat < 1.0 - tol)) ) {
			// if we have not been converged before on this scalar value
			if (!rhatDiagnosticFlag)  {
				#ifdef MYDEBUG
					cout << "\n" << getScalarsName() 
						<< " convergence test satisfied in state " 
						  << states << " (states in calcs = " 
						  << (states - statesNotInCalcs) << ")"<< endl;
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
	
		#ifdef MYDEBUG_CALCS
			cout << "\nsampling" << endl;
			cout << "\nstates = " << states << endl;
			cout << "statesNotInCalcs = " << statesNotInCalcs << endl;
		#endif	

	}
	else { // not updating, just replicate values
		calculationStates.push_back(calculationStates.back());
		overallIntervalLengths.push_back(overallIntervalLengths.back());
		interchainIntervalStatistic.push_back(interchainIntervalStatistic.back());
		
		rhat.push_back(rhat.back()); 
		
	}
	if (keepLogs) {
		// store the flag as well, as a real, which is a fudge...
		rhatFlag.push_back(rhatDiagnosticFlag);
	}
	// old value if not updated
	return rhatDiagnosticFlag;
} // end calculations


void MCMCGRDiagnosticInterval::outputResults(
				const std::string& filenameGRScalars,
				const RealVec& sampledInd,
				int precData) const
{
	/* sampledInd may contain more values than we monitored */
	RealVec::const_iterator it = sampledInd.begin();
	advance(it, scalars.size());
	RealVec sampledIndTmp(sampledInd.begin(), it);
	
	std::vector < std::string > colNames;
	getScalarColNames(colNames);
	colNames.push_back(interchainIntervalStatisticColName);
	colNames.push_back(overallIntervalLengthsColName);
	colNames.push_back("rhat");
	
		colNames.push_back("rhatFlag");
		colNames.push_back("sampled?");
	
	std::vector < const RealVec* > data;
	addDataPtrs(data, scalars);
	data.push_back(&interchainIntervalStatistic);
	data.push_back(&overallIntervalLengths);
	data.push_back(&rhat);
	
		data.push_back(&rhatFlag);
		data.push_back(&sampledIndTmp);
	
	outputToFileVertical(calculationStates, data, colNames, 
	filenameGRScalars, precData);
} 

void MCMCGRDiagnosticInterval::outputResults(
			const std::string& filenameGRScalars,
			int precData) const
{
	std::vector <std::string > colNames;
	getScalarColNames(colNames);
	colNames.push_back(interchainIntervalStatisticColName);
	colNames.push_back(overallIntervalLengthsColName);
	colNames.push_back("rhat");
	
	std::vector < const RealVec* > data;
	addDataPtrs(data, scalars);
	data.push_back(&interchainIntervalStatistic);
	data.push_back(&overallIntervalLengths);
	data.push_back(&rhat);
	
	outputToFileVertical(calculationStates, data, colNames, 
	filenameGRScalars, precData);
} 




void MCMCGRDiagnosticInterval::outputCalculations(
				const std::string& filenameGRWorkingCalcs,
				int precData) const
{
	
	std::vector <std::string > colNames;
	
	getScalarColNames(colNames);
	
	/* output working calcs: all v's for each chain, 
	 * running sums for each chain, sample variances,
	 * overall running sums */
	colNames.insert(colNames.end(), intervalLengthsColNames.begin(), intervalLengthsColNames.end());
	colNames.push_back(interchainIntervalStatisticColName);
	colNames.push_back(overallIntervalLengthsColName);
	
	std::vector < const RealVec* > data;
	data = addDataPtrs(data, scalars);
	data = addDataPtrs(data, intervalLengths);
	data.push_back(&interchainIntervalStatistic);
	data.push_back(&overallIntervalLengths);
	
	outputToFileVertical(calculationStates, data, colNames, 
						filenameGRWorkingCalcs, precData);
	
}



const std::vector < std::vector < cxsc::real > >& 
			MCMCGRDiagnosticInterval::getScalarsRef() const
{
	return scalars;
}

cxsc::real MCMCGRDiagnosticInterval::getCurrentRhatValue() const
{
	if (rhat.empty())
		throw std::runtime_error(
			"MCMCGRDiagnosticInterval::getCurrentRhatValue(size_t): no value to get");
	return rhat.back();
}

cxsc::real MCMCGRDiagnosticInterval::getTolerance() const
{
	return tol;
}

// private 

std::vector <std::string >& MCMCGRDiagnosticInterval::getScalarColNames(
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


/* Interval length is calculated as 
 * upper 1-(1-percent)/2 point - lower (1-percent)/2 point
 * 
 For percentile = (1-(1-percent)/2) and then ( (1-percent)/2 ),
 the method used to calculate the percentile point
 uses interpolation and is supposed to be equivalent
 to the method used for gsl_stats_quantile_from_sorted_data
 (see the GNU Scientific Library).
 
 if percentile = p, 0 < p < 1, then
 
 percentile points = (1-delta)*set[i] + delta*set[i+1]
 
 where 
 set[i] is the value in position i of the (ordered) set
 n = size of data set
 i = floor((n-1)p)
 delta = (n-1)p - i
 */
real MCMCGRDiagnosticInterval::getIntervalLength(
						const std::multiset < cxsc::real >& set) const
{
	size_t n = set.size();
	
	#ifdef MYDEBUG_CALCS_EXTRA
		cout << "set size (n) = " << n << endl;
	#endif
	real lb(0.0);
	real ub(0.0);
		
	double interpolPt = (n-1)*p;
	size_t i = static_cast<size_t>(std::floor(interpolPt));

	double delta = interpolPt - i;
	
	#ifdef MYDEBUG_CALCS_EXTRA
		cout << "interpolPt = " << interpolPt << ", i = " << i << endl;
		cout << "delta = " << delta << endl;
	#endif
		
	{
		std::multiset < cxsc::real >::const_iterator it = set.begin();
		advance(it, i);
		lb = (1-delta) * (*it);
		advance(it, 1);
		lb += delta * (*it);
	}
	{
		//use reverse iterator
		std::multiset < cxsc::real >::const_reverse_iterator it = set.rbegin();
		advance(it, i);
		ub = delta * (*it); // this is the 'top' end
		advance(it, 1);
		ub += (1-delta) * (*it); // this is the 'bottom' end
	}
	#ifdef MYDEBUG_CALCS_EXTRA
		cout << "lb = " << lb << "\tub = " << ub << endl;
	#endif
	assert(!(ub < lb));
	return (ub-lb); 
}

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
\brief Definitions for type collecting info for automcmc from adaptive histograms.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "adaptivehistogram_changeofstateinfo_automcmc.hpp"

using namespace cxsc;
using namespace std;

#include <climits>

namespace subpavings {


ChangeOfStateInformationAutoMCMC::ChangeOfStateInformationAutoMCMC(
							real lp,
							size_t cl, size_t cc, 
							unsigned long int tld)
	: logPosterior(lp), deltaPi(0.0), 
		currentLeaves(cl), currentCherries(cc),
		totalLeafDepth (tld)
	{}

ChangeOfStateInformationAutoMCMC::ChangeOfStateInformationAutoMCMC(
					real lp,
					const AdaptiveHistogram& adh)
	: logPosterior(lp), deltaPi(0.0),
		currentLeaves(0), currentCherries(0),
			totalLeafDepth (0)
{
	currentLeaves = adh.getRootLeaves();
	currentCherries = adh.getRootCherries();
	totalLeafDepth = adh.getRootTotalLeafDepth();
}

ChangeOfStateInformationAutoMCMC::ChangeOfStateInformationAutoMCMC(
					real lp,
					size_t cl, size_t cc, 
					const AdaptiveHistogram& adh)
	: logPosterior(lp), deltaPi(0.0),
	currentLeaves(cl), currentCherries(cc)
{
	totalLeafDepth = adh.getRootTotalLeafDepth();
}

/*Change in log posterior from last change. */			
real ChangeOfStateInformationAutoMCMC::getDeltaPi() const
{
	return deltaPi;
}

/* Get current log posterior from last change. */			
real ChangeOfStateInformationAutoMCMC::getCurrentLogPosterior() const
{
	return logPosterior;
}
			

/*Get current number of leaves. */
size_t ChangeOfStateInformationAutoMCMC::getCurrentLeaves() const
{
	return currentLeaves;
}

/*Get current number of cherries. */
size_t ChangeOfStateInformationAutoMCMC::getCurrentCherries() const
{
	return currentCherries;
}


/*Get current average leaf depth. */
real ChangeOfStateInformationAutoMCMC::getAverageLeafDepth() const
{
	return (totalLeafDepth*1.0)/currentLeaves;
}



/* Set change in log posterior from last change. */			
void ChangeOfStateInformationAutoMCMC::notifyDeltaPi(
								real dp)
{
	deltaPi = dp;
	logPosterior+=deltaPi;
}

/* Notify of split. */
void ChangeOfStateInformationAutoMCMC::notifySplit(
								const SPnode * const spn)
{
	++currentLeaves;
	if (!spn->hasLeafSibling()) ++currentCherries;
	
	size_t thisNodeDepth = 
			static_cast<unsigned long int>(spn->getNodeDepth());
	if (ULONG_MAX - (thisNodeDepth + 2) < totalLeafDepth) {
		throw std::runtime_error(
			"ChangeOfStateInformationAutoMCMC::notifySplit(...)");
	}
	totalLeafDepth += thisNodeDepth + 2;
}

/* Notify of merge. */
void ChangeOfStateInformationAutoMCMC::notifyMerge(
								const SPnode * const spn)
{
	--currentLeaves;
	if (!spn->hasLeafSibling()) --currentCherries;
	
	size_t thisNodeDepth = 
			static_cast<unsigned long int>(spn->getNodeDepth());
	
	totalLeafDepth -= (thisNodeDepth + 2);

}




}// end namespace subpavings

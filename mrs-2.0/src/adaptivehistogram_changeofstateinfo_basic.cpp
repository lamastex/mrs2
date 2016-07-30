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
\brief Definitions for type collecting basic info from adaptive histograms.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "adaptivehistogram_changeofstateinfo_basic.hpp"

using namespace cxsc;
using namespace std;



namespace subpavings {


ChangeOfStateInformationBasic::ChangeOfStateInformationBasic(
	size_t cl, size_t cc)
	: deltaPi(0.0), currentLeaves(cl), currentCherries(cc)
	{}

ChangeOfStateInformationBasic::ChangeOfStateInformationBasic(
					const AdaptiveHistogram& adh)
	: deltaPi(0.0), currentLeaves(0), currentCherries(0)
{
	currentLeaves = adh.getRootLeaves();
	currentCherries = adh.getRootCherries();
}

/*Change in log posterior from last change. */			
real ChangeOfStateInformationBasic::getDeltaPi() const
{
	return deltaPi;
}

/*Get current number of leaves. */
size_t ChangeOfStateInformationBasic::getCurrentLeaves() const
{
	return currentLeaves;
}

/*Get current number of cherries. */
size_t ChangeOfStateInformationBasic::getCurrentCherries() const
{
	return currentCherries;
}



/* Set change in log posterior from last change. */			
void ChangeOfStateInformationBasic::notifyDeltaPi(
								real dp)
{
	deltaPi = dp;
}

/* Notify of split. */
void ChangeOfStateInformationBasic::notifySplit(
								const SPnode * const spn)
{
	++currentLeaves;
	if (!spn->hasLeafSibling()) ++currentCherries;
	
}

/* Notify of merge. */
void ChangeOfStateInformationBasic::notifyMerge(
								const SPnode * const spn)
{
	--currentLeaves;
	if (!spn->hasLeafSibling()) --currentCherries;
}




}// end namespace subpavings

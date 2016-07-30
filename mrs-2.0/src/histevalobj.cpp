/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief Definitions for classes for evaluating when to stop changing histograms.
*/

#include "histevalobj.hpp"
#include "adaptivehistogram.hpp"

namespace subpavings {
	
	//base class destructor
	HistEvalObj::~HistEvalObj(){}

	/*Class for testing the number of bins of a histogram.
	*/
	CritLeaves_GTE::CritLeaves_GTE(size_t t) : test(t) {}

	/* True if the number of leaves is >= test. */
	bool CritLeaves_GTE::operator()
			(const AdaptiveHistogram * const adh) const
	{
		return (adh->getRootLeaves() >= test);
	}



	/* Class for testing the number of bins of a histogram.
	*/
	CritLeaves_LTE::CritLeaves_LTE(size_t t) : test(t) {}

	/* True if the number of leaves is <= test. */
	bool CritLeaves_LTE::operator()(const AdaptiveHistogram * const adh) const
	{
		return (adh->getRootLeaves() <= test);
	}


	/* Class for testing the count of the node with the smallest
	count in histogram's subpaving.
	*/
	CritSmallestCount_LTE ::CritSmallestCount_LTE(size_t t) : test(t) {}

	/* True if the smallest leaf count is <= test. */
	bool CritSmallestCount_LTE::operator()(const AdaptiveHistogram * const adh) const
	{
		return (((adh->getSubPaving())->getSmallestLeafCount()) <= test);
	}
	
	/* \brief Class for testing the count of the node with the largest
	count in histogram's subpaving.
	*/
	CritLargestCount_LTE::CritLargestCount_LTE(size_t t) : test(t) {}

	/* True if the largest leaf count is <= test. */
	bool CritLargestCount_LTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		return (((adh->getSubPaving())->getLargestLeafCount()) <= test);
	}

	
	
	/* Class for testing the count of the node with the largest
	count or the maximum number of leaves in histogram's subpaving.
	*/
	CritLargestCount_LTE_CritLeaves_GTE::CritLargestCount_LTE_CritLeaves_GTE(size_t tCount, size_t tLeaves) 
				: testCount(tCount), testLeaves(tLeaves) {}

	/* True if the largest leaf count is <= testCount OR 
	 * if the number of leaves is >= testLeaves. */
	bool CritLargestCount_LTE_CritLeaves_GTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		return 
		((((adh->getSubPaving())->getLargestLeafCount()) <= testCount)
		||
		(adh->getRootLeaves() >= testLeaves));
	}


	/* Class for testing the volume of the box with the smallest
	volume in the histogram's subpaving.
	*/
	CritSmallestVol_LTE::CritSmallestVol_LTE(double t) : test(t) {}

	/* True if the volume of the smallest leaf (by vol) is <= test. */
	bool CritSmallestVol_LTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		return ((adh->getSubPaving())->getSmallestLeafVol() <= test);
	}


	/* Class for testing the volume of the box with the largest
	volume in the histogram's subpaving.
	*/
	CritLargestVol_LTE::CritLargestVol_LTE(double t) : test(t) {}

	/* True if the volume of the largest leaf (by vol) is <= test. */
	bool CritLargestVol_LTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		return ((adh->getSubPaving())->getLargestLeafVol() <= test);
	}


	/* Class for testing the count of the node with the largest
	count or the maximum number of leaves in histogram's subpaving.
	*/
	CritLargestVol_LTE_CritLeaves_GTE::CritLargestVol_LTE_CritLeaves_GTE(
					double vol, size_t tLeaves) 
				: testVol(vol), testLeaves(tLeaves) {}

	/* True if the largest leaf volume is <= testVol OR 
	 * if the number of leaves is >= testLeaves. */
	bool CritLargestVol_LTE_CritLeaves_GTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		return 
		((((adh->getSubPaving())->getLargestLeafVol()) <= testVol)
		||
		(adh->getRootLeaves() >= testLeaves));
	}


	
	
	/* \brief Class for testing the change in COPERR score from splitting.

	\warning <b>Never</b> just use CritCOPERRChange on its own:  it could keep
	splitting for ever.  Use CritCOPERRChangeOrLeaves_GTE instead.
	*/
	CritCOPERRChange_GTE::CritCOPERRChange_GTE(
					const PenObj& p, double t): pen(p), test(t) {}

	/* True if the change in COPERR score from splitting best node to
	split >= test. */
	bool CritCOPERRChange_GTE::operator()(
					const AdaptiveHistogram * const adh) const
	{
		dotprecision best = adh->getSubPaving()->getBestSplitChangeEMPCOPERR(
								adh->getSubPaving()->getCounter());
		real change = rnd(best) + pen(adh, 1) - pen(adh, 0);

		return (change >= test);
	}



	/* Class for testing the change in AIC score from splitting.

	\warning <b>Never</b> just use critAICChange on its own:  it could keep
	splitting for ever.  Use critAICChangeOrLeaves_GTE instead.
	*/
	CritAICChange_GTE::CritAICChange_GTE(
				const PenObj& p, double t) : pen(p), test(t) {}

	/* True if the change in AIC score from splitting best node to
	split >= test. */
	bool CritAICChange_GTE::operator()(const AdaptiveHistogram * const adh) const
	{
		dotprecision best = adh->getSubPaving()->getBestSplitChangeEMPAIC();
		real change = rnd(best) + pen(adh, 1) - pen(adh, 0);

		return (change >= test);
	}

	

	/* Class for testing change in COPERR or number leaves from splitting.
	*/
	CritCOPERRChangeOrLeaves_GTE::CritCOPERRChangeOrLeaves_GTE(
						const PenObj& p, size_t tl, double ts)
						: pen(p), testLeaves(tl), testScore(ts)  {}

	/* True if the change in COPERR score from splitting best node to
	split >= testScore OR if number of leaves is >= testLeaves. */
	bool CritCOPERRChangeOrLeaves_GTE::operator()(
				const AdaptiveHistogram * const adh) const
	{
		CritLeaves_GTE critLeaves(testLeaves);
		CritCOPERRChange_GTE critScore(pen, testScore);

		return (critLeaves(adh) || critScore(adh));
	}



	/* Class for testing change in AIC or number leaves from splitting.
	*/
	CritAICChangeOrLeaves_GTE::CritAICChangeOrLeaves_GTE(const PenObj& p, size_t tl, double ts)
					: pen(p), testLeaves(tl), testScore(ts)  {}

	/* True if the change in AIC score from splitting best node to
	split >= testScore OR if number of leaves is >= testLeaves. */
	bool CritAICChangeOrLeaves_GTE::operator()(const AdaptiveHistogram * const adh) const
	{
		CritLeaves_GTE critLeaves(testLeaves);
		CritAICChange_GTE critScore(pen, testScore);

		return (critLeaves(adh) || critScore(adh));
	}


	/* Class for testing change in COPERR or largest count from splitting.

	\warning this stopping rule can get stuck: largest count may not stop splitting.
	*/
	CritCOPERRChangeOrLargestCount_LTE::CritCOPERRChangeOrLargestCount_LTE(
									const PenObj& p, size_t tc, double ts)
									: pen(p), testCount(tc), testScore(ts)  {}

	/* True if the change in COPERR score from splitting best node to
	split >= testScore OR if the largest count in any leaf if <= testCount. */
	bool CritCOPERRChangeOrLargestCount_LTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		CritLargestCount_LTE critCount(testCount);
		CritCOPERRChange_GTE critScore(pen, testScore);

		return (critCount(adh) || critScore(adh));
	}


	/* Class for testing change in AIC or largest count from splitting.

	\warning this stopping rule can get stuck: largest count may not stop splitting.
	*/
	CritAICChangeOrLargestCount_LTE::CritAICChangeOrLargestCount_LTE(const PenObj& p,
											size_t tc, double ts)
										: pen(p), testCount(tc), testScore(ts)  {}

	/* True if the change in AIC score from splitting best node to
	split >= testScore OR if the largest count in any node is <= testCount. */
	bool CritAICChangeOrLargestCount_LTE::operator()(const AdaptiveHistogram * const adh) const
	{
		CritLargestCount_LTE critCount(testCount);
		CritAICChange_GTE critScore(pen, testScore);

		return (critCount(adh) || critScore(adh));
	}


	/* Class for testing the change in COPERR score from merging.

	*/
	CritCOPERRMergeChange_GTE::CritCOPERRMergeChange_GTE(
										const PenObj& p, double t)
											: pen(p), test(t)  {}

	/* True if the change in COPERR score from merging best node to
	merge >= test. */
	bool CritCOPERRMergeChange_GTE::operator()(
					const AdaptiveHistogram * const adh) const
	{

		dotprecision best = adh->getSubPaving()->getBestSplitChangeEMPCOPERR(
								adh->getSubPaving()->getCounter());
		real change = rnd(best) + pen(adh, -1) - pen(adh, 0); // merge = -1 leaf

		return (change >= test);
	}


	/* Class for testing the change in AIC score from merging.

	True if the change in AIC score from merging best node to merge >= test.
	*/
	CritAICMergeChange_GTE::CritAICMergeChange_GTE(const PenObj& p, double t) : pen(p), test(t) {}

	/* True if the change in AIC score from merging best node to
	merge >= test. */
	bool CritAICMergeChange_GTE::operator()(const AdaptiveHistogram * const adh) const
	{
		dotprecision best = adh->getSubPaving()->getBestMergeChangeEMPAIC();

		real change = rnd(best) + pen(adh, -1) - pen(adh, 0); // merge = -1 leaf

		return (change >= test);
	}



	/* Class for testing change in COPERR or number leaves from merging.
	*/
	CritCOPERRMergeChangeOrLeaves_LTE::CritCOPERRMergeChangeOrLeaves_LTE(
					const PenObj& p, size_t tl, double ts)
							: pen(p), testLeaves(tl), testScore(ts)  {}

	/* True if the change in COPERR score from merging best node to
	merge >= testScore OR if number of leaves is <= testLeaves. */
	bool CritCOPERRMergeChangeOrLeaves_LTE::operator()(
					const AdaptiveHistogram * const adh) const
	{
		CritLeaves_LTE critLeaves(testLeaves);
		CritCOPERRMergeChange_GTE critScore(pen, testScore);

		return (critLeaves(adh) || critScore(adh));
	}


	/* Class for testing change in AIC or  number leaves from merging.

	True if the change in AIC score from merging best node to merge >= 0
	OR if number of leaves is <= test.
	*/
	CritAICMergeChangeOrLeaves_LTE::CritAICMergeChangeOrLeaves_LTE(
						const PenObj& p, size_t tl, double ts)
						: pen(p), testLeaves(tl), testScore(ts)  {}

	/* True if the change in AIC score from merging best node to
	merge >= testScore OR if number of leaves is <= testLeaves. */
	bool CritAICMergeChangeOrLeaves_LTE::operator()(
						const AdaptiveHistogram * const adh) const
	{
		CritLeaves_LTE critLeaves(testLeaves);
		CritAICMergeChange_GTE critScore(pen, testScore);

		return (critLeaves(adh) || critScore(adh));
	}



	/* Class to bale out of priority queue splitting.
	*/
	bool CritStopAll::operator()(const AdaptiveHistogram * const) const
		{ return true; }



} // end namespace subpavings





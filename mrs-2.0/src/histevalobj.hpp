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
\brief Declarations for classes for evaluating when to stop changing histograms.
*/

#ifndef ___HISTEVAL_HPP__
#define ___HISTEVAL_HPP__

// to use histogram penalty function objects
#include "histpenalty.hpp"

namespace subpavings {
	//! Forward class declarations
	class AdaptiveHistogram;


	/*! \brief A Virtual class providing a way to stop histogram changes.
	*/
	class HistEvalObj {

		public:

		/*! return true when splitting should stop. */
		virtual bool operator()(const AdaptiveHistogram * const adh) const = 0;
		
		virtual ~HistEvalObj();
	};


	/** @name Concrete classes derived from HistEvalObj

	These classes provide a stopping test for priority queue splitting or merging.
	The function operator () returns true when splitting or merging should stop.

	*/
	//@{


	/*! \brief Class for testing the number of bins of a histogram.
	*/
	class CritLeaves_GTE : public HistEvalObj
	{
		size_t test;

		CritLeaves_GTE(); // private default constructor

		public:

		explicit CritLeaves_GTE(size_t t);

		/*! True if the number of leaves is >= test. */
		bool operator()
			(const AdaptiveHistogram * const adh) const;
	};


	/*! \brief Class for testing the number of bins of a histogram.
	*/
	class CritLeaves_LTE : public HistEvalObj
	{
		size_t test;

		CritLeaves_LTE(); // private default constructor

		public:

		explicit CritLeaves_LTE(size_t t);

		/*! True if the number of leaves is <= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing the count of the node with the smallest
	count in histogram's subpaving.
	*/
	class CritSmallestCount_LTE : public HistEvalObj
	{
		size_t test;

		CritSmallestCount_LTE(); // private default constructor

		public:

		explicit CritSmallestCount_LTE(size_t t);

		/*! True if the smallest leaf count is <= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};
	
	
	/*! \brief Class for testing the count of the node with the largest
	count in histogram's subpaving.
	*/
	class CritLargestCount_LTE : public HistEvalObj
	{
		size_t test;

		CritLargestCount_LTE(); // private default constructor

		public:

		explicit CritLargestCount_LTE(size_t t);

		/*! True if the largest leaf count is <= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
		
	};
	
	
	/*! \brief Class for testing the count of the node with the largest
	count or the maximum number of leaves in histogram's subpaving.
	*/
	class CritLargestCount_LTE_CritLeaves_GTE : public HistEvalObj
	{
		size_t testCount;
		size_t testLeaves;

		CritLargestCount_LTE_CritLeaves_GTE(); // private default constructor

		public:

		CritLargestCount_LTE_CritLeaves_GTE(size_t tCount, size_t tLeaves);

		/*! True if the largest leaf count is <= testCount OR 
		 * if the number of leaves is >= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};

	/*! \brief Class for testing the volume of the box with the smallest
	volume in the histogram's subpaving.
	*/
	class CritSmallestVol_LTE : public HistEvalObj
	{
		double test;

		CritSmallestVol_LTE(); // private default constructor

		public:

		explicit CritSmallestVol_LTE(double t);

		/*! True if the volume of the smallest leaf (by vol) is <= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing the volume of the box with the largest
	volume in the histogram's subpaving.
	*/
	class CritLargestVol_LTE : public HistEvalObj
	{
		double test;

		CritLargestVol_LTE(); // private default constructor

		public:

		explicit CritLargestVol_LTE(double t);

		/*! True if the volume of the largest leaf (by vol) is <= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing the count of the node with the largest
	count or the maximum number of leaves in histogram's subpaving.
	*/
	class CritLargestVol_LTE_CritLeaves_GTE : public HistEvalObj
	{
		double testVol;
		size_t testLeaves;

		CritLargestVol_LTE_CritLeaves_GTE(); // private default constructor

		public:

		CritLargestVol_LTE_CritLeaves_GTE(double vol, size_t tLeaves);

		/*! True if the largest leaf volume is <= testVol OR 
		 * if the number of leaves is >= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};
	
	
	/*! \brief Class for testing the change in COPERR score from splitting.

	\warning <b>Never</b> just use CritCOPERRChange on its own:  it could keep
	splitting for ever.  Use CritCOPERRChangeOrLeaves_GTE instead.
	*/
	class CritCOPERRChange_GTE : public HistEvalObj
	{
		const PenObj& pen;
		double test;

		CritCOPERRChange_GTE(); // private default constructor

		public:

		CritCOPERRChange_GTE(const PenObj& p, double t);

		/*! True if the change in COPERR score from splitting best node to
		split >= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};


	/*! \brief Class for testing the change in AIC score from splitting.

	\warning <b>Never</b> just use critAICChange on its own:  it could keep
	splitting for ever.  Use critAICChangeOrLeaves_GTE instead.
	*/
	class CritAICChange_GTE : public HistEvalObj
	{
		const PenObj& pen;
		double test;

		CritAICChange_GTE(); // private default constructor

		public:

		CritAICChange_GTE(const PenObj& p, double t);

		/*! True if the change in AIC score from splitting best node to
		split >= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};

	/*! \brief Class for testing change in COPERR or number leaves from splitting.
	*/
	class CritCOPERRChangeOrLeaves_GTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testLeaves;
		double testScore;

		CritCOPERRChangeOrLeaves_GTE(); // private default constructor

		public:

		CritCOPERRChangeOrLeaves_GTE(const PenObj& p, size_t tl, double ts);

		/*! True if the change in COPERR score from splitting best node to
		split >= testScore OR if number of leaves is >= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;

	};

	/*! \brief Class for testing change in AIC or number leaves from splitting.
	*/
	class CritAICChangeOrLeaves_GTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testLeaves;
		double testScore;

		CritAICChangeOrLeaves_GTE(); // private default constructor

		public:

		CritAICChangeOrLeaves_GTE(const PenObj& p, size_t tl, double ts);

		/*! True if the change in AIC score from splitting best node to
		split >= testScore OR if number of leaves is >= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing change in COPERR or largest count from splitting.

	\warning this stopping rule can get stuck: largest count may not stop splitting.
	*/
	class CritCOPERRChangeOrLargestCount_LTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testCount;
		double testScore;

		CritCOPERRChangeOrLargestCount_LTE(); // private default constructor

		public:

		CritCOPERRChangeOrLargestCount_LTE(const PenObj& p,
											size_t tc, double ts);

		/*!True if the change in COPERR score from splitting best node to
		split >= testScore OR if the largest count in any leaf if <= testCount. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing change in AIC or largest count from splitting.

	\warning this stopping rule can get stuck: largest count may not stop splitting.
	*/
	class CritAICChangeOrLargestCount_LTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testCount;
		double testScore;

		CritAICChangeOrLargestCount_LTE(); // private default constructor

		public:

		CritAICChangeOrLargestCount_LTE(const PenObj& p,
											size_t tc, double ts);

		/*! True if the change in AIC score from splitting best node to
		split >= testScore OR if the largest count in any node is <= testCount. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};



	/*! \brief Class for testing the change in COPERR score from merging.

	*/
	class CritCOPERRMergeChange_GTE : public HistEvalObj
	{
		const PenObj& pen;
		double test;

		CritCOPERRMergeChange_GTE(); // private default constructor

		public:

		CritCOPERRMergeChange_GTE(const PenObj& p, double t);

		/*! True if the change in COPERR score from merging best node to
		merge >= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};

	/*! \brief Class for testing the change in AIC score from merging.

	True if the change in AIC score from merging best node to merge >= test.
	*/
	class CritAICMergeChange_GTE : public HistEvalObj
	{
		const PenObj& pen;
		double test;

		CritAICMergeChange_GTE(); // private default constructor

		public:

		CritAICMergeChange_GTE(const PenObj& p, double t);

		/*! True if the change in AIC score from merging best node to
		merge >= test. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};


	/*! \brief Class for testing change in COPERR or number leaves from merging.
	*/
	class CritCOPERRMergeChangeOrLeaves_LTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testLeaves;
		double testScore;

		CritCOPERRMergeChangeOrLeaves_LTE(); // private default constructor

		public:

		CritCOPERRMergeChangeOrLeaves_LTE(const PenObj& p,
										size_t tl, double ts);

		/*! True if the change in COPERR score from merging best node to
		merge >= testScore OR if number of leaves is <= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};


	/*! \brief Class for testing change in AIC or  number leaves from merging.

	True if the change in AIC score from merging best node to merge >= 0
	OR if number of leaves is <= test.
	*/
	class CritAICMergeChangeOrLeaves_LTE : public HistEvalObj
	{
		const PenObj& pen;
		size_t testLeaves;
		double testScore;

		CritAICMergeChangeOrLeaves_LTE(); // private default constructor

		public:

		CritAICMergeChangeOrLeaves_LTE(const PenObj& p,
										size_t tl, double ts);

		/*! True if the change in AIC score from merging best node to
		merge >= testScore OR if number of leaves is <= testLeaves. */
		bool operator()(const AdaptiveHistogram * const adh) const;
	};


	/*! \brief Class to bale out of priority queue splitting.
	*/
	class CritStopAll: public HistEvalObj
	{
		public:

		// use default constructor

		/*! True always */
		bool operator()(const AdaptiveHistogram * const) const;
	};

	//@}
} // end namespace subpavings

#endif



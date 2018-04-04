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
\brief Declarations for classes for evaluating when to stop
 changing function estimators (using intervals).
*/

#ifndef ___FEI_EVAL_HPP__
#define ___FEI_EVAL_HPP__

#include "real.hpp"

#include <cstdlib>



namespace subpavings {
	//! Forward class declarations
	class FunctionEstimatorInterval;


	/*! \brief An interface for a type providing a way to stop estimator changes.
	*/
	class FEIEvalObj {

		public:

		/*! return true when splitting should stop. */
		virtual bool operator()(
				const FunctionEstimatorInterval& fei) const = 0;
		
		virtual ~FEIEvalObj() {};
	};


	/** @name Concrete classes derived from FEIEvalObj

	These classes provide a stopping test for priority queue splitting or merging.
	The function operator () returns true when splitting or merging should stop.

	*/
	//@{


	/*! \brief Class for testing the number of leaves in estimator.
	*/
	class FEICritLeaves_GTE : public FEIEvalObj
	{
		public:

			explicit FEICritLeaves_GTE(size_t t);

			/*! True if the number of leaves is >= test. */
			bool operator()(
						const FunctionEstimatorInterval& fei) const;
				
		private: 
			size_t test;

			FEICritLeaves_GTE(); // private default constructor

		
	};


	/*! \brief Class for testing the number of leaves in estimator.
	*/
	class FEICritLeaves_LTE : public FEIEvalObj
	{
		public:

			explicit FEICritLeaves_LTE(size_t t);

			/*! True if the number of leaves is <= test. */
			bool operator()(
						const FunctionEstimatorInterval& fei) const;
			
		private: 
			size_t test;

			FEICritLeaves_LTE(); // private default constructor

		
	};

	/*! \brief Class for testing the interval band on the estimator.
	*/
	class FEICritIntervalBand_LTE : public FEIEvalObj
	{
		public:

			explicit FEICritIntervalBand_LTE(cxsc::real t);

			/*! True if the band is <= test. */
			bool operator()(
						const FunctionEstimatorInterval& fei) const;
			
		private: 
			cxsc::real test;

			FEICritIntervalBand_LTE(); // private default constructor

		
	};


	/*! \brief Class to bale out of priority queue splitting.
	*/
	class FEICritStopAll: public FEIEvalObj
	{
		public:

		// use default constructor

			/*! True always */
			bool operator()(
						const FunctionEstimatorInterval& fei) const;
	};

	//@}
} // end namespace subpavings

#endif



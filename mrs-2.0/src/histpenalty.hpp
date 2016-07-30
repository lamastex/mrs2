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
\brief Declaration of classes for histogram penalty objects declarations.
*/

#ifndef ___HISTPEN_HPP__
#define ___HISTPEN_HPP__

#include "real.hpp"

namespace subpavings {
	
	
	//! Forward class declarations
	class AdaptiveHistogram;


	/*! \brief Abstract class for objects with penalty function for histogram fit.

	The penalty function is used in a histogram fit scoring method such as
	COPERR or AIC.
	*/

	class PenObj
	{

		public:

		/*! \brief the penalty function.
		\param adh is the histogram object to calculate the penalty on.
		\param deltaLeaf change in number leaves to take into account.
		\return the penalty as a cxsc::real.
		*/
		virtual cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const = 0;
		
		virtual ~PenObj();
	};


	/** @name Classes derived from PenObj

	These classes provide a penalty  function is used in a histogram fit scoring
	method such as COPERR or AIC.

	*/
	//@{

	/*! \brief Penalty function as number of leaves in histogram
	*/
	class PenLeaves : public PenObj
	{
		double beta;

		public:

		PenLeaves();
		PenLeaves(double b);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
		int deltaLeaf) const;
	};

	/*! \brief Class for penalty function 1 for AIC.
	*/
	class PenAIC1 : public PenObj
	{
		double c;
		double alpha;
		double r;

		public:

		// default constructor has c = 1, alpha = 0.5, r = 2
		PenAIC1();

		PenAIC1(double cc, double aa, double rr);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const;
	};

	/*! \brief Class for penalty function 2 for AIC.
	*/
	class PenAIC2 : public PenObj
	{
		double c;
		double alpha;
		double r;

		// default constructor private
		PenAIC2();

		public:

		PenAIC2(double cc, double aa, double rr);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const;
	};

	/*! \brief Class for penalty function 3 for AIC.

	Same as PenAIC2 at present.
	*/
	class PenAIC3 : public PenObj
	{
		double c;
		double alpha;
		double r;

		// default constructor private
		PenAIC3();

		public:

		PenAIC3(double cc, double aa, double rr);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const;
	};

	/*! \brief Class for penalty function 4 for AIC.

	Penalty is just some multiple of logCatK.
	Should have c = 0.5 to stop splitting on one point with AIC EMP,
	or with c=2, penalty before taking logs is (CatK)^2).
	*/
	class PenAIC4 : public PenObj
	{
		double c;

		// default constructor private
		PenAIC4();

		public:

		PenAIC4(double cc);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const;
	};

	/*! \brief Class for penalty function 5 for AIC.

	Penalty is just some multiple of logCatK.
	Should have c = 1 to stop splitting on one point with AIC EMP.
	*/
	class PenAIC5 : public PenObj
	{
		double c;

		// default constructor private
		PenAIC5();

		public:

		PenAIC5(double cc);

		cxsc::real operator() (const AdaptiveHistogram * const adh,
									int deltaLeaf) const;
	};

	//@}

	} // end namespace subpavings

#endif



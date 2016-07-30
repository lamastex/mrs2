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
\brief Definitions of classes for histogram penalty objects declarations.
*/
#include "histpenalty.hpp"
#include "adaptivehistogram.hpp"
#include "toolz.hpp"

namespace subpavings {
	
	 PenObj::~PenObj(){}


	/* Penalty function as number of leaves in histogram
	*/
	PenLeaves::PenLeaves() : beta(1.0) {}
	PenLeaves::PenLeaves(double b) : beta(b) {}

	real PenLeaves::operator() (const AdaptiveHistogram * const adh,
	int deltaLeaf) const
	{
		real penalty;
		penalty = beta*(adh->getSubPaving()->getNumberLeaves()) + deltaLeaf;
		return penalty;
	}


	/* Class for penalty function 1 for AIC.
	*/
	// default constructor has c = 1, alpha = 0.5, r = 2
	PenAIC1::PenAIC1() : c(1.0), alpha (0.50), r(2.0) {}

	PenAIC1::PenAIC1(double cc, double aa, double rr)
			: c(cc), alpha (aa), r(rr) {}

	real PenAIC1::operator() (const AdaptiveHistogram * const adh,
								int deltaLeaf) const
	{
		dotprecision penalty(0.0);

		int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
		double logCatk= lCk(k);

		accumulate(penalty, c, logCatk); // pen = c*logCatk
		accumulate(penalty, alpha, k); // pen = c*logCatk + alpha*k
		accumulate(penalty, c*r, log(k+1.0));
		// now pen = c*logCatk + alpha*k + c*r*log(k+1)

		dotprecision temp(0.0);
		accumulate(temp, c*alpha*k, logCatk);
		accumulate(temp, r, log(k+1.0));
		// temp = c*alpha*k*logCatk + r*log(k+1)

		return (rnd(penalty) + 2*sqrt(rnd(temp)));
		// return c*logCatk + alpha*k + c*r*log(k+1)
		//                          + 2*sqrt(c*alpha*k*logCatk + r*log(k+1))
	}


	/* Class for penalty function 2 for AIC.
	*/
	PenAIC2::PenAIC2(double cc, double aa, double rr)
			: c(cc), alpha (aa), r(rr) {}

	real PenAIC2::operator() (const AdaptiveHistogram * const adh,
								int deltaLeaf) const
	{
		dotprecision penalty(0.0);
		int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
		double logCatk= lCk(k);
		size_t counter = adh->getRootCounter(); // total number points
		accumulate(penalty, c, logCatk); // pen = c*logCatk

		real s = adh->getRootSumLeafCountOverVol();

		accumulate(penalty, s/(1.0*counter), alpha);
		//now pen = c*logCatk + (alpha/counter)*sum(leaf counts/vols)
		accumulate(penalty, 2*r, log(k+1));
		//now pen = c*logCatk + (alpha/counter)*sum(leaf counts/vols) + 2*r*log(k+1)

		return rnd(penalty);
	}



	/* Class for penalty function 3 for AIC.

	Same as PenAIC2 at present.
	*/
	
	PenAIC3::PenAIC3(double cc, double aa, double rr)
			: c(cc), alpha (aa), r(rr) {}

	real PenAIC3::operator() (const AdaptiveHistogram * const adh,
								int deltaLeaf) const
	{
		dotprecision penalty(0.0);
		int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
		double logCatk= lCk(k);
		size_t counter = adh->getRootCounter(); // total number points
		accumulate(penalty, c, logCatk); // pen = c*logCatk

		real s = adh->getRootSumLeafCountOverVol();

		accumulate(penalty, s/(1.0*counter), alpha);
		//now pen = c*logCatk + (alpha/counter)*sum(leaf counts/vols)
		accumulate(penalty, 2*r, log(k+1));
		//now pen = c*logCatk + (alpha/counter)*sum(leaf counts/vols) + 2*r*log(k+1)

		return rnd(penalty);
	}


	/* Class for penalty function 4 for AIC.

	Penalty is just some multiple of logCatK.
	Should have c = 0.5 to stop splitting on one point with AIC EMP,
	or with c=2, penalty before taking logs is (CatK)^2).
	*/
		PenAIC4::PenAIC4(double cc) : c(cc) {}

	real PenAIC4::operator() (const AdaptiveHistogram * const adh,
								int deltaLeaf) const
	{
		dotprecision penalty(0.0);

		int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
		double logCatk= lCk(k);

		accumulate(penalty, c, logCatk); // pen = c*logCatk

		return (rnd(penalty));
		// return c*logCatk
	}


	/* Class for penalty function 5 for AIC.

	Penalty is just some multiple of logCatK.
	Should have c = 1 to stop splitting on one point with AIC EMP.
	*/
	PenAIC5::PenAIC5(double cc) : c(cc) {}

	real PenAIC5::operator() (const AdaptiveHistogram * const adh,
								int deltaLeaf) const
	{
		dotprecision penalty(0.0);

		int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
		
		accumulate(penalty, c*k, log(2.0)); // pen = c*k*log(2)

		return (rnd(penalty));
	}
} // end namespace subpavings





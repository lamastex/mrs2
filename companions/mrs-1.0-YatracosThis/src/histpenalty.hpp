/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
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

/*! \file      histpenalty.hpp
\brief Declaration of classes for histogram penalty objects declarations.
*/

#ifndef ___HISTPEN_HPP__
#define ___HISTPEN_HPP__

double lCk(const int k);

//using namespace std;
//using namespace cxsc;

namespace subpavings {
	
	//! Forward class declarations
	class AdaptiveHistogram;
	class AdaptiveHistogramValidation;
	
	
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
	    \return the penalty as a real.
	    */
	    virtual real operator() (const AdaptiveHistogram * const adh,
	                                int deltaLeaf) const = 0;
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
	
	    PenLeaves() : beta(1.0) {}
	    PenLeaves(double b) : beta(b) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
	    int deltaLeaf) const
	    {
	        real penalty;
	        penalty = beta*spLeaves(adh->getSubPaving()) + deltaLeaf;
	        return penalty;
	    }
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
	    PenAIC1() : c(1.0), alpha (0.50), r(2.0) {}
	
	    PenAIC1(double cc, double aa, double rr)
	            : c(cc), alpha (aa), r(rr) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
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
	
	    PenAIC2(double cc, double aa, double rr)
	            : c(cc), alpha (aa), r(rr) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
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
	
	    PenAIC3(double cc, double aa, double rr)
	            : c(cc), alpha (aa), r(rr) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
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
	
	    PenAIC4(double cc) : c(cc) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
	                                int deltaLeaf) const
	    {
	        dotprecision penalty(0.0);
	
	        int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
	        double logCatk= lCk(k);
	
	        accumulate(penalty, c, logCatk); // pen = c*logCatk
	
	        return (rnd(penalty));
	        // return c*logCatk
	    }
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
	
	    PenAIC5(double cc) : c(cc) {}
	
	    real operator() (const AdaptiveHistogram * const adh,
	                                int deltaLeaf) const
	    {
	        dotprecision penalty(0.0);
	
	        int k = adh->getRootLeaves() + deltaLeaf - 1; // leaves-1
	//        double logCatk= lCk(k);
	
	        accumulate(penalty, c*k, log(2.0)); // pen = c*k*log(2)
	
	        return (rnd(penalty));
	        // return c*k*log(2)
	    }
	};
	
	//@}
	
} // end of namespace subpavings
	
	#endif
	
	
	

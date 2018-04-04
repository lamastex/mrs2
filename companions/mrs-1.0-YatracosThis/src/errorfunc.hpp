/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
* Copyright (C) 2011 Gloria Teng
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

/*! \file      errorfunc.hpp
\brief L1-error function declarations.
*/

#ifndef ___ERRORFUNC_HPP__
#define ___ERRORFUNC_HPP__

#include "sptypes.hpp"
#include "adaptivehistogram.hpp"
//#include "adaptivehistogramvalidation.hpp"
#include <gsl/gsl_rng.h>        // to know about the gsl random number generator
#include <gsl/gsl_randist.h>    // we need these libraries to get the IAE for
										  // finite mixtures
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

#include "../examples/StatsSubPav/ExactInt/Int.h"
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"

#include "intervalw.h" //for interval routines
#include "ia_ad.h" //for interval routines

#include "itaylor.h" //for taylor integration routine

#include <vector>

using namespace subpavings;
using namespace cxsc;

//! Declarations of typedef taylor. This is used for 2D integrations.
typedef taylor::dim2taylor d2t;
typedef taylor::dim2taylor_vector d2tv;

namespace subpavings {

//! Forward class declarations
class AdaptiveHistogram;
class AdaptiveHistogramValidation;

} // end of namespace subpavings

/** @name Family of functions for the bivariate gaussian
*/ 
//@{
/** The integral of the absolute error for a bivariate gaussian. */
d2t BiGOP (d2tv X, interval fhat); 
//@}

/** @name Family of functions for the Levy 2D
*/ 
//@{
/** The integral of the absolute error for a Levy 2D. */
d2t LevyOP (d2tv X, interval fhat);
//@}

/** @name Family of functions for the Rosenbrock 2D
*/ 
//@{
/** The integral of the absolute error for a Rosenbrock 2D. */
d2t RosenOP (d2tv X, interval fhat); 
//@}



/*! \brief Get the mean of the data (this is used when building the regular histogram)
*/
//check if Jenny already has code for this 
double myMean(const RVecData& rvec);
/*! \brief get the standard deviation of the data (this is used when building the regular histogram)
*/
double myStd(const RVecData& rvec);
/*! \brief Gaussian probability density function. 
		Generates normal probability density values corresponding to X which is a 
		vector of doubles.
	 \param XX is a row vector of support points for which normal density 
	         values are required;
    \param SS is a row vector of standard deviations.
    \return NPD: row vector of normal probability density values.
*/
std::vector<double> gaussian(std::vector<double> & vecNPD, 
                             std::vector<double> & XX, double SS);
                             
//----the true density for gaussian mixtures----------------//
//for root finding routine
ia_ad f(const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S);

//for integration routine
typedef itaylor (*pfcn)(const itaylor &, vector<double>&, vector<double>&, vector<double>&);

itaylor integrand(const itaylor &x, vector<double>& W, vector<double>& M, vector<double>& S);                             
                             
//----------functions for finding roots----------------------------------//
//return the interval function?
intervalw F (const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S); 

//return the derivative of the function?
intervalw DF(const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S); 

//newton's routine
// need to modify this to suit fhat?
intervalw N (const intervalw &x, vector<double>& W, 
					vector<double>& M, vector<double>& S, double fhat);

//find the root using interval newton method
void findRoot(const intervalw &domain, vector<double>& W, vector<double>& M, 
					vector<double>& S, double fhat, vector<intervalw>& rootList);

//bisect the domain and decide which root-finding routine to use
void bisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootList, vector<double>& W, 
				vector<double>& M, vector<double>& S);

//----------functions for integration----------------------------------
interval riemannTerm(pfcn f, interval X, int Deg, vector<double>& W, 
				vector<double>& M, vector<double>& S);
				
interval integrate(pfcn f, interval X, int Deg, double Tol, vector<double>& W, 
				vector<double>& M, vector<double>& S);

interval getL1error(double fhat, interval& thisInt, int Deg, double TOL, 
							vector<double>& W, vector<double>& M, vector<double>& S);

//------For Laplace--------------------------------//
//for root finding routine
ia_ad LaplacePDF(const intervalw &x);

//for integration routine
typedef itaylor (*pfcnLaplace)(const itaylor &);

itaylor LaplaceIntegrand(const itaylor &x);                             
                             
//----------functions for finding roots----------------------------------//
//return the interval function?
intervalw LaplaceF (const intervalw &x); 

//return the derivative of the function?
intervalw LaplaceDF(const intervalw &x); 

//newton's routine
// need to modify this to suit fhat?
intervalw LaplaceNewton (const intervalw &x, double fhat);

//find the root using interval newton method
void LaplacefindRoot(const intervalw &domain, double fhat, vector<intervalw>& rootList);

//bisect the domain and decide which root-finding routine to use
void LaplaceBisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootList);

//----------functions for integration----------------------------------
interval LaplaceRiemannTerm(pfcn f, interval X, int Deg);
				
interval LaplaceIntegrate(pfcn f, interval X, int Deg, double Tol);

interval LaplaceGetL1error(double fhat, interval& thisInt, int Deg, double TOL);

//-------------end of functions for Laplace---------------------------//

//--------------For Standard LogNormal--------------------------------//
//for root finding routine
ia_ad LognormalPDF(const intervalw &x);

//for integration routine
typedef itaylor (*pfcnLognormal)(const itaylor &);

itaylor LognormalIntegrand(const itaylor &x);                             
                             
//----------functions for finding roots----------------------------------//
//return the interval function?
intervalw LognormalF (const intervalw &x); 

//return the derivative of the function?
intervalw LognormalDF(const intervalw &x); 

//newton's routine
// need to modify this to suit fhat?
intervalw LognormalNewton (const intervalw &x, double fhat);

//find the root using interval newton method
void LognormalfindRoot(const intervalw &domain, double fhat, vector<intervalw>& rootList);

//bisect the domain and decide which root-finding routine to use
void LognormalBisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootList);

//----------functions for integration----------------------------------
interval LognormalRiemannTerm(pfcn f, interval X, int Deg);
				
interval LognormalIntegrate(pfcn f, interval X, int Deg, double Tol);

interval LognormalGetL1error(double fhat, interval& thisInt, int Deg, double TOL);

/** @name Family of functions for finite mixtures.
*/
//@{
/** Structure for the mixture model.*/
struct FinMix {
        std::vector<double> W; //weights of components
        std::vector<double> M;  //means of components
        std::vector<double> S;  //standard deviations of components
        double fhat;
};
/** A finite mixture distribution is the sum of c gaussian components. 
    Its PDF (CDF) is obtained by summing the PDFs (CDF) of each component. The 
	 integrated absolute error is calculated using a quadrature routine from the GSL library.  
	\param x is the data point at which we take the PDF (CDF).
	\param x1 (only for the IAE at the boundary) is the infimum of the root box. 
	\param x2 (only for the IAE at the boundary) is the supremum of the root box.
	\param W is a vector of the weights of each component.
	\param M is a vector of the means of each component.
	\param S is a vector of the standard deviations of each component.
	\return the integrated absolute error.
*/
//@{
/** get the PDF of a finite mixture */
double FinMixPDF(double x, std::vector<double> &W, 
                          std::vector<double> &M, std::vector<double> &S);
/** get the CDF of a finite mixture */
double FinMixCDF(double, double, std::vector<double> &W,
                          std::vector<double> &M, 
                          std::vector<double> &S);
/** get the integrated absolute of a finite mixture of boxes in the root box */
double FinMixAbs(double x, void* params);
/** get the integrated absolute of a finite mixture at the boundaries */
dotprecision dpFinMixIAEBoun(double x1, double x2, FinMix& mixt) ;

/** Find the number of generated U(0,1) data that is less the weight of the 
    components.
   \param  u is an array that holds U(0,1) rv.
   \param  intp is always 0
   \param	n is the total number of data points
   \param	weight is the weights of the mixture
   \param	w is the weight component that we want to find the number of 
	         members of
*/
void findComp(vector<double> & u, int& intp, const int n, 
                        double* weight, int w);

/** Function to get cumulative sum
  	 \param weight is the weights of members in mixtures
    \param	w is the component of weight
*/
void cumsum(vector<double> & weight, double* w);
//@}

/** @name Family of functions for the regular histogram including error 
          calculations.
			Bandwidth method:
			0 = Scott's normal reference rule,
			1 = Wand's one-stage rule,
			2 = Wand's two-stage rule,
*/
//@{
/** Structure for a regular histogram. 
 * probably should make this into a class
*/
struct RegHist {
		std::vector<real> LowerBoxes;
		std::vector<real> UpperBoxes;
		std::vector<double> heights;
		cxsc::real binwidth;
};

/** Make a regular histogram, i.e. histogram of equal width. Only for 
    one-dimensional data. */
void makeRegularHist(RegHist& myRegHist, const RVecData& rvec, 
							ivector theBox, int bwmethod);
void makeRegularHist(RegHist& myRegHist, const RVecData& rvec, 
							ivector theBox, double bw);

/** Get psi (this is used when building the regular histogram). 
    Reference: Wand M.P. (1997), "Data-based choice of histogram bin width",
	 American Statistician 51, 59-64: Equations (2.2), (4.1)-(4.4).
*/
double psi(const RVecData& rvec, double g, double r);

/** Get the IAE for a uniform mixture distribution for regular histograms. */
real getRegHistUnifIAE(RegHist &myRegHist, 
                       AdaptiveHistogram &myPart, size_t n, double weight,
                       vector<int> holesLoc);

/** Get the IAE for a finite mixture distribution for regular histograms */
real getRegHistFinMixIAE(size_t n, RegHist &myRegHist, FinMix & mixt);

/** Get the interval IAE for a finite mixture distribution for regular histograms */
interval getRegHistFinMixIntervalIAE(size_t n, RegHist &myRegHist, FinMix & mixt, double tol, int deg);

/** Get the interval IAE for a finite mixture distribution for regular histograms */
interval getRegHistLaplaceIntervalIAE(size_t n, RegHist &myRegHist, double tol, int deg);

/** Get the interval IAE for a finite mixture distribution for regular histograms */
interval getRegHistLognormalIntervalIAE(size_t n, RegHist &myRegHist, double tol, int deg);

/** Output the regular histogram to .txt file*/
void outputRegHistToTxt(RegHist &myRegHist, std::string& s);
//@}



#endif

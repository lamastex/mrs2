/*
* Copyright (C) 2011, 2012 Jennifer Harlow
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
\brief Testing PiecewiseConstantFunction log likelihood
 */
#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
#include "adaptivehistogram.hpp"
#include "subpaving_exception.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>
#include <cfloat> // for DBL_EPSILON
#include <cmath>
#include <limits>

#include <gsl/gsl_math.h> //gsl_isnan(), isinf

using namespace cxsc;
using namespace std;
using namespace subpavings;



void testLogLik()
{
	cout << cxsc::ln(cxsc::MinReal) << endl; //cxsc::ln(0.0) error
	cout << std::log(0.0) << endl;
	
	double dinf = numeric_limits<double>::infinity();
	real rinf = cxsc::Infinity;
	real rneginf = -cxsc::Infinity;
	real rsnan = cxsc::SignalingNaN;
	real rqnan = cxsc::QuietNaN;
	
	cout << "riinf = cxsc::IsInfinity " << rinf << endl;
	cout << "rneginf = -cxsc::IsInfinity " << rneginf << endl;
	cout << "rneginf < 0.0 ? " << (rneginf < 0.0) << endl;
	
	//cant do cxsc::ln(rinf)
	
	cout << "cxsc::real(1.0/0.0)? " << cxsc::real(1.0/0.0) << endl;
	cout << "cxsc::IsInfinity(cxsc::real(1.0/0.0))? " << cxsc::IsInfinity(cxsc::real(1.0/0.0)) << endl;
	cout << "cxsc::IsInfinity(rinf)? " << cxsc::IsInfinity(rinf) << endl;
	cout << "cxsc::IsInfinity(rneginf)? " << cxsc::IsInfinity(rneginf) << endl;
	cout << "cxsc::IsInfinity(1.0/0.0)? " << cxsc::IsInfinity(1.0/0.0) << endl;
	cout << "cxsc::IsInfinity(dinf)? " << cxsc::IsInfinity(dinf) << endl;
	cout << "gsl_isinf(rinf)? " << gsl_isinf(_double(rinf)) << endl;
	cout << "gsl_isinf(rneginf)? " << gsl_isinf(_double(rneginf)) << endl;
	cout << "gsl_isinf(dinf)? " << gsl_isinf(_double(dinf)) << endl;
	cout << "gsl_isinf(-dinf)? " << gsl_isinf(_double(-dinf)) << endl;
	cout << "gsl_isinf(1/0.0)? " << gsl_isinf(1/0.0) << endl;
	cout << "gsl_isinf(-1/0.0)? " << gsl_isinf(-1/0.0) << endl;
	
	cout << "cxsc::exp(rneginf) " << cxsc::exp(rneginf) << endl;
	cout << "std::exp(-dinf) " << std::exp(-dinf) << endl;
	
	cxsc::dotprecision dotneginf(-cxsc::Infinity);
	
	cout << "dotneginf " << dotneginf << endl;
	cout << "cxsc::rnd(dotneginf) " << cxsc::rnd(dotneginf) << endl;
	cout << "cxsc::IsInfinity(cxsc::rnd(dotneginf) " 
			<< cxsc::IsInfinity(cxsc::rnd(dotneginf)) << endl;
	cout << "cxsc::exp(cxsc::rnd(dotneginf)) "  << cxsc::exp(cxsc::rnd(dotneginf)) << endl;
	cout << "gsl_isinf(_double(cxsc::rnd(dotneginf))) " 
			<< gsl_isinf(_double(cxsc::rnd(dotneginf))) << endl;
	cout << "cxsc::rnd(dotneginf)) < 0.0 ? " 
			<< (cxsc::rnd(dotneginf) < 0.0) << endl;
	cout << "cxsc::rnd(dotneginf)) > 0.0 ? " 
			<< (cxsc::rnd(dotneginf) > 0.0) << endl;
	
	cxsc::dotprecision dotinf(cxsc::Infinity);
	
	cout << "dotinf " << dotinf << endl;
	cout << "cxsc::rnd(dotinf) " << cxsc::rnd(dotinf) << endl;
	cout << "cxsc::IsInfinity(cxsc::rnd(dotinf) " 
			<< cxsc::IsInfinity(cxsc::rnd(dotinf)) << endl;
	cout << "gsl_isinf(_double(cxsc::rnd(dotinf))) " 
			<< gsl_isinf(_double(cxsc::rnd(dotinf))) << endl;
	cout << "cxsc::rnd(dotinf)) < 0.0 ? " 
			<< (cxsc::rnd(dotinf) < 0.0) << endl;
	cout << "cxsc::rnd(dotinf)) > 0.0 ? " 
			<< (cxsc::rnd(dotinf) > 0.0) << endl;
	// cant do cxsc::exp(cxsc::rnd(dotinf))
	
	
	cout << rsnan << endl;
	cout << rqnan << endl;
	
	cout << "cxsc::IsSignalingNaN(rsnan)? " << cxsc::IsSignalingNaN(rsnan) << endl;
	cout << "cxsc::IsQuietNaN(rqnan)? " << cxsc::IsQuietNaN(rqnan) << endl;
	
	cout << "cxsc::IsSignalingNaN(rqnan)? " << cxsc::IsSignalingNaN(rqnan) << endl;
	cout << "cxsc::IsQuietNaN(rsnan)? " << cxsc::IsQuietNaN(rsnan) << endl;
	
	// cannot do cxsc::exp(rqnan)
	// cannot do cxsc::exp(rsnan)
	
	
	
	
	int prec = 5; // default precision for output files
	// make a simple pcf and try log likelihoods
	int d = 2; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
	{
		cout << "\nTest log lik with no data in hist" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		cxsc::real loglik = pcf.getLogLikelihood(hist);
		
		cxsc::real shouldBe = 0.0;
		assert (loglik == shouldBe);
		cout << "loglik = " << loglik << endl;
	}
	
	{
		cout << "\nTest log lik with negative values in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesNegative(pcf.getDomainVolume());
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		RVecData data;
		data = getData2(data);
		bool successfulInsertion = hist.insertFromRVec(data);
		if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
		
		assert(successfulInsertion);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
		}
	}
	{
		cout << "\nTest log lik with some negative and some positive values in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesPositiveAndNegative();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		{
			RVecData data;
			data = getData2RightOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
		
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		{
			RVecData data;
			data = getData2LeftOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
		}
	}
	
	{
		cout << "\nTest log lik with infinite ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesInfinite();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		RVecData data;
		data = getData2(data);
		bool successfulInsertion = hist.insertFromRVec(data);
		if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
		
		assert(successfulInsertion);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isinf(_double(loglik)) && loglik > 0);
			cout << "loglik = " << loglik << endl;
		}
	}
	
	{
		cout << "\nTest log lik with infinite and negative ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesInfiniteAndNegative();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		{
			RVecData data;
			data = getData2LeftOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
			
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isinf(_double(loglik)) && loglik > 0.0);
			
		}
		{
			RVecData data;
			data = getData2RightOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
		
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
		}
	}
	{
		cout << "\nTest log lik with infinite and zero ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesInfiniteAndZero();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		RVecData data;
		data = getData2(data);
		bool successfulInsertion = hist.insertFromRVec(data);
		if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
		
		assert(successfulInsertion);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isinf(_double(loglik)) && loglik < 0.0);
			cout << "loglik = " << loglik << endl;
			cxsc::real lik = cxsc::exp(loglik);
			cout << "cxsc::exp(loglik) = " << lik << endl;
		}
	}
	
	{
		cout << "\nTest log lik with infinite and negative infinite ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesInfiniteAndNegativeInfinite();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		RVecData data;
		data = getData2(data);
		bool successfulInsertion = hist.insertFromRVec(data);
		if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
		
		assert(successfulInsertion);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
			
		}
	}
	
	{
		cout << "\nTest log lik with positive and negative infinite ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesPositiveAndNegativeInfinite();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		RVecData data;
		data = getData2(data);
		bool successfulInsertion = hist.insertFromRVec(data);
		if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
		
		assert(successfulInsertion);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
			
		}
	}
	
	{
		cout << "\nTest log lik with negative and zero ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesNegativeAndZero();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		{
			RVecData data;
			data = getData2LeftOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isinf(_double(loglik)) && loglik < 0.0);
			
		}
		
		{
			RVecData data;
			data = getData2RightOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isnan(_double(loglik)));
			cout << "loglik = " << loglik << endl;
			
		}
	}
	{
		cout << "\nTest log lik with positive and zero ranges in pcf" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRangesPositiveAndZero();
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		{
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		
		{
			RVecData data;
			data = getData2LeftOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
		
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			cxsc::real shouldBe = 0.0;
			assert (loglik == shouldBe);
		}
		{
			RVecData data;
			data = getData2RightOnly(data);
			bool successfulInsertion = hist.insertFromRVec(data);
			if (!successfulInsertion) cout << "unsuccessful insertion" << endl;
			
			assert(successfulInsertion);
		
			cxsc::real loglik = pcf.getLogLikelihood(hist);
			assert(gsl_isinf(_double(loglik)) && loglik < 0.0);
			cout << "loglik = " << loglik << endl;
			cxsc::real lik = cxsc::exp(loglik);
			cout << "cxsc::exp(loglik) = " << lik << endl;
		}
	}
	
	
	{
		cout << "\nTest log lik with simple hist" << endl;
		PiecewiseConstantFunction pcf(pavingBox);
		std::string split = "1,1";
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
		pcf.allocateRanges(ranges);
		
		AdaptiveHistogram hist(pavingBox);
		hist.splitToShape(split);
		
		cxsc::real loglik = pcf.getLogLikelihood(hist);
		
		//cxsc::real shouldBe = 0.0;
		//assert (loglik == shouldBe);
		cout << "Piecewise constant function is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
		cout << "loglik = " << loglik << endl;
	}
	
	#if(0)
	int prec = 5; // default precision for output files
		
	int d = 2; // dimension of the box to sample data from
	ivector pavingBox(d);
	interval pavingInterval(-2,2);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	try {
		
		cout << "\nTest making pcfs from histograms" << endl;
		
		
			
		try {
			cout << "\nTest constructing pcf from histogram with no subpaving (this should fail)" << endl;
			
			AdaptiveHistogram adh;
			PiecewiseConstantFunction pcf(adh);
			
			throw std::logic_error("Should not be able to get here");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe ) {
			std::string msg(nspe.what());
			cout << "\nFailed to construct pcf from histogram with no subpaving:\n" << msg << endl;
		}		
		{	
			cout << "\nTest constructing pcf from histogram (this should be okay)" << endl;
	
			std::string split1 = "3,4,4,2,2,4,4,3";
			bool holdAllStats = false;
			int lab = 2;
			AdaptiveHistogram hist1(pavingBox, holdAllStats, lab);
			
			assert (hist1.getLabel() == lab);
			hist1.splitToShape(split1);

			RVecData data1;
			data1 = getData1(data1);
			bool successfulInsertionHistFirst = hist1.insertFromRVec(data1);
			if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
			
			assert(successfulInsertionHistFirst);
			
			cout << "Histogram is " << endl;
			hist1.outputToStreamTabs(cout, prec);
			
			PiecewiseConstantFunction pcf(hist1);
			
			assert (pcf.getLabel() == lab);
			
			cxsc::real integral = pcf.getTotalIntegral();
			
			string s("PCFfromHist1.txt");
			pcf.outputToTxtTabs(s, prec, true);
			cout << "pcf.getTotalIntegral() = " << integral << endl;
			assert (integral == cxsc::real(1.0));
			cout << "Piecewise constant function is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
		}
		
		cout << "\nEnd of make pcf from histogram tests:\n" << endl;		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do make pcf from histogram tests:\n" << msg << endl;
		throw;
	}
	#endif
}		
	

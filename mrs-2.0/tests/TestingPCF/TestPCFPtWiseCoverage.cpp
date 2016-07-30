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
\brief Testing PiecewiseConstantFunction pointwise extension and coverage
 */

#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
#include "functionestimator_real.hpp"
#include "simpleFobj2.hpp"
#include "subpaving_exception.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>
#include <cfloat> // for DBL_EPSILON


using namespace cxsc;
using namespace std;
using namespace subpavings;


void TestPointWiseAndCoverage()
{
	int d = 2; // dimension of the box to sample data from
	ivector pavingBox(d);
	interval pavingInterval(0,1);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	int prec = 5;
	
	SimpleFobj2 fobj;

	std::string split = "2,3,4,4,2,3,3";
	
	
	
	try {
		
		cout << "\n\n Pointwise extensions:\n" << endl;
		
		try {
			
			cout << "\nTry with pcf with no subpaving (this should fail)" << endl;
			PiecewiseConstantFunction pcfNull;
			rvector pt(d);
			pt[1] = 1.0;
			pt[2] = 2.0; 
			cxsc::real r = pcfNull.pointwiseExtension(pt);
			throw std::logic_error("Should not be able to get here");
		}
		catch (NullSubpavingPointer_Error& nspe) {
			cout << "\nNullSubpavingPointer_Error: " << nspe.what() << endl;
		}
		
		{
		
			cout << "\nTest pointwise extension with negative values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			
			{
				rvector pt(d);
				pt[1] = 0.0;
				pt[2] = 1.0; // XL
				
				cxsc::real r = pcf.pointwiseExtension(pt);
				cxsc::real shouldBe(ranges[1]);
				assert ( r == shouldBe);
				cout << "Passed assert that pointwiseExtension(";
				prettyPrint(cout, pt);
				cout << ") = " << shouldBe << endl;
			}
			{
				rvector pt(d);
				pt[1] = 0.5;
				pt[2] = 1.0; // XR
				
				cxsc::real r = pcf.pointwiseExtension(pt);
				cxsc::real shouldBe(ranges[2]);
				assert ( r == shouldBe);
				cout << "Passed assert that pointwiseExtension(";
				prettyPrint(cout, pt);
				cout << ") = " << shouldBe << endl;
			}
							
		}
		{
		
			cout << "\nTest pointwise extension with infinite values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasInfinitePiecewiseConstantValues());
			
			{
				rvector pt(d);
				pt[1] = 0.0;
				pt[2] = 1.0; // XL
				
				cxsc::real r = pcf.pointwiseExtension(pt);
				cxsc::real shouldBe(ranges[1]);
				assert ( r == shouldBe);
				cout << "Passed assert that pointwiseExtension(";
				prettyPrint(cout, pt);
				cout << ") = " << shouldBe << endl;
			}
			{
				rvector pt(d);
				pt[1] = 0.5;
				pt[2] = 1.0; // XR
				
				cxsc::real r = pcf.pointwiseExtension(pt);
				cxsc::real shouldBe(ranges[2]);
				assert ( r == shouldBe);
				cout << "Passed assert that pointwiseExtension(";
				prettyPrint(cout, pt);
				cout << ") = " << shouldBe << endl;
			}
							
		}
		
		
		FunctionEstimatorReal fer(pavingBox, fobj);
		PiecewiseConstantFunction pcf;
		
		try {
			fer.splitToShape(split);
			
			pcf = fer.makePiecewiseConstantFunction();
			
			string s1 = "PCFfromSimpleFobj2Fei.txt";
			pcf.outputToTxtTabs(s1, prec, true);
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to make pcf from FunctionEstimatorReal:\n" << msg << endl;
			throw;
		}
		try {
			
			cout << "\nTry point with wrong dimensions (this should fail)" << endl;
			rvector pt(d+1);
			for (int k = 1; k <= d+1; ++k) pt[k] = 1.0;
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			throw std::logic_error("Should not be able to get here");
		}
		catch (IncompatibleDimensions_Error& ice) {
			cout << "\nIncompatibleDimensions_Error: " << ice.what() << endl;
		}
		
		
		{
			rvector pt(d);
			pt[1] = 1.0;
			pt[2] = 2.0; // not in box
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(0.0);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 0.0; // XLL
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(0.25);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 1.0;
			pt[2] = 1.0; // XRRR
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(2.625);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.5;
			pt[2] = 0.5; // XRRR
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(1.875);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.25;
			pt[2] = 0.75; // XLRRR
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(1.3125);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.25;
			pt[2] = 0.5; // XLRRL
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe(0.9375);
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do real pointwise extension tests:\n" << msg << endl;
		throw;
	}
	try {
		
		cout << "\n\n Coverage region:\n" << endl;
		
		try {
			
			cout << "\nTry with pcf with no subpaving (this should fail)" << endl;
			PiecewiseConstantFunction pcfNull;
			real cov = 0.0;
			pcfNull.outputCoverageRegion(cout, cov);
			throw std::logic_error("Should not be able to get here");
		}
		catch (NullSubpavingPointer_Error& nspe) {
			cout << "\nNullSubpavingPointer_Error: " << nspe.what() << endl;
		}
		
		try {
		
			cout << "\nTest coverage region with negative values in pcf (this should fail)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			
			real cov = 0.0;
			pcf.outputCoverageRegion(cout, cov);
			throw std::logic_error("Should not be able to get here");
							
		}
		catch (runtime_error& re) {
			cout << "\nruntime_error:\n" << re.what() << endl;
		}
		
		try {
		
			cout << "\nTest coverage region with infinite values in pcf (this should fail)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasInfinitePiecewiseConstantValues());
			
			real cov = 0.0;
			pcf.outputCoverageRegion(cout, cov);
			throw std::logic_error("Should not be able to get here");
		
		}
		catch (runtime_error& re) {
			cout << "\nruntime_error:\n" << re.what() << endl;
		}		
		
		FunctionEstimatorReal fer(pavingBox, fobj);
		
		try {
			fer.splitToShape(split);
			
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do splitToShape:\n" << msg << endl;
			throw;
		}
		PiecewiseConstantFunction pcf = fer.makePiecewiseConstantFunction();
		
		
		try {
		
			cout << "\nTest coverage region cov < 0.0 (this should fail)" << endl;
			real cov = -0.1;
			pcf.outputCoverageRegion(cout, cov);
			throw std::logic_error("Should not be able to get here");
			
							
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nstd::invalid_argument:\n" << msg << endl;
		}
		try {
		
			cout << "\nTest coverage region cov > 1.0 (this should fail)" << endl;
			real cov = 1.1;
			pcf.outputCoverageRegion(cout, cov);
			throw std::logic_error("Should not be able to get here");
			
							
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nstd::invalid_argument:\n" << msg << endl;
		}
		{
			real cov = 0.0;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.1;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.328124;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.328125;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.644531124;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.64453125;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 0.999999;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 1.0;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do coverage region tests:\n" << msg << endl;
		throw;
	}
	
	try {
		
		cout << "\n\n Coverage for a point:\n" << endl;
		
		FunctionEstimatorReal fer(pavingBox, fobj);
		
		try {
			
			cout << "\nTry with pcf with no subpaving (this should fail)" << endl;
			PiecewiseConstantFunction pcfNull;
			rvector pt(d);
			pt[1] = 1.0;
			pt[2] = 2.0; 
			cxsc::real r = pcfNull.findCoverage(pt);
			throw std::logic_error("Should not be able to get here");
		}
		catch (NullSubpavingPointer_Error& nspe) {
			cout << "\nNullSubpavingPointer_Error: " << nspe.what() << endl;
		}
		
		try {
		
			cout << "\nTest coverage with negative values in pcf (this should fail)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 1.0; // outside
			
			real cov = pcf.findCoverage(pt);
			
			throw std::logic_error("Should not be able to get here");
			
							
		}
		catch (std::runtime_error& re) {
			std::string msg(re.what());
			cout << "\nruntime_error:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nTest coverage with infinite values in pcf (this should fail)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasInfinitePiecewiseConstantValues());
			
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 1.0; // outside
			
			real cov = pcf.findCoverage(pt);
			
			throw std::logic_error("Should not be able to get here");
			
							
		}
		catch (std::runtime_error& re) {
			std::string msg(re.what());
			cout << "\nruntime_error:\n" << msg << endl;
		}
		
		
		try {
			fer.splitToShape(split);
			
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do splitToShape:\n" << msg << endl;
			throw;
		}
		PiecewiseConstantFunction pcf = fer.makePiecewiseConstantFunction();
		
		try {
			
			cout << "\nTry point with wrong dimensions (this should fail)" << endl;
			rvector pt(d+1);
			for (int k = 1; k <= d+1; ++k) pt[k] = 1.0;
			
			cxsc::real r = pcf.findCoverage(pt);
			throw std::logic_error("Should not be able to get here");
		}
		catch (IncompatibleDimensions_Error& ice) {
			cout << "\nIncompatibleDimensions_Error: " << ice.what() << endl;
		}
		{
			rvector pt(d);
			pt[1] = -1.0;
			pt[2] = -1.0; // outside
			
			real cov = pcf.findCoverage(pt);
			real shouldBe = 0.0;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.75;
			pt[2] = 0.5; // XRRR
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.5;
			pt[2] = 0.5; // XRRL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.328125;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.25;
			pt[2] = 0.75; // XLRRR
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.5625;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.25;
			pt[2] = 0.5; // XLRRL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.64453125;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.5;
			pt[2] = 0.0; // XRL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.703125;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 0.5; // XLRL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.890625;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 0.0; // XLL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0 - 0.9375;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do findCoverage tests:\n" << msg << endl;
		throw;
	}
	
	try {
		
		cout << "\n\nTest pointwiseExtension, coverageRegion, findCoverage when pcf is not normalised:\n" << endl;
		
		PiecewiseConstantFunction pcf(pavingBox);
		
		std::string split = "1,1";
		
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRanges1(0.5*pcf.getDomainVolume());
		
		pcf.allocateRanges(ranges);
					
		ostringstream oss;
		pcf.outputRootToStreamTabs(oss, prec);
		string pcfRoot = oss.str();
		cout << "Unnormalised pcf is " << endl;
		cout << pcfRoot << endl;
		
		{
			rvector pt(d);
			pt[1] = 0.5;
			pt[2] = 0.0; // XLRRL
			
			cxsc::real r = pcf.pointwiseExtension(pt);
			cxsc::real shouldBe = ranges[2];
			assert ( r == shouldBe);
			cout << "Passed assert that pointwiseExtension(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << "\n" << endl;
		} 
		
		{
			real cov = 0.75;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		{
			real cov = 1.0;
			cout << "Coverage region for cov = " << cov << " is" << endl;
			pcf.outputCoverageRegion(cout, cov);
			
		}
		
		{
			rvector pt(d);
			pt[1] = 0.0;
			pt[2] = 0.0; // XL
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 0.25;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		{
			rvector pt(d);
			pt[1] = 0.5;
			pt[2] = 0.0; // XR
			
			real cov = pcf.findCoverage(pt);
			cout << "findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << cov << endl;
			real shouldBe = 1.0;
			assert(cov == shouldBe);
			cout << "Passed assert that findCoverage(";
			prettyPrint(cout, pt);
			cout << ") = " << shouldBe << endl;
		}
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do unnormalised pointwise and coverage tests:\n" << msg << endl;
		throw;
	}	
	cout << "\nEnd of pointwise and coverage tests:\n" << endl;	
		
		
}


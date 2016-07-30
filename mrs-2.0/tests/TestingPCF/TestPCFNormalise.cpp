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
\brief Testing PiecewiseConstantFunction normalising
 */

#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
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


void testNormalise()
{
	int prec = 5;
	try {
		try {
			cout << "\ncall makeNormalised on pcf with no subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf;
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to call makeNormalised on pcf with no subpaving:\n" << msg << endl;
		}
		
		
		int d = 2; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(0,1);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		try {
			cout << "\nmakeNormalised with pcf with no values (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			assert(pcf.getRootLeaves() == 1);
			assert(pcf.getTotalIntegral() == 0.0);
			
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::runtime_error& re) {
			std::string msg(re.what());
			cout << "\nFailed to call makeNormalised on pcf with no values:\n" << msg << endl;
		}
		
		try {
			
			cout << "\nTest normalise with infinite values in pcf (this should fail)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasInfinitePiecewiseConstantValues());
			assert (pcf.getTotalIntegral() == cxsc::Infinity);
			
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::runtime_error& re) {
			std::string msg(re.what());
			cout << "\nFailed to call makeNormalised on pcf with infinite values:\n" << msg << endl;
		}
		{
			
			cout << "\nTest normalise with negative values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			cxsc::real intUnnorm = pcf.getTotalIntegral();
			assert (intUnnorm != 1.0);
			
			cout << "unnormalised pcf with negative values has integral " << intUnnorm << " and is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			cxsc::real intNorm = norm.getTotalIntegral();
			assert (intNorm == 1.0);
			cout << "normalised pcf with negative values has integral " << intNorm << " and is " << endl;
			norm.outputRootToStreamTabs(cout, prec);
							
		}
		
		
		
		{
			
			cout << "\nTest normalise with positive values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			cxsc::real intUnnorm = pcf.getTotalIntegral();
			assert (intUnnorm != 1.0);
			
			cout << "unnormalised pcf has integral " << intUnnorm << " and is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			cxsc::real intNorm = norm.getTotalIntegral();
			assert (intNorm == 1.0);
			cout << "normalised pcf has integral " << intNorm << " and is " << endl;
			norm.outputRootToStreamTabs(cout, prec);
							
		}
		
		{
			
			cout << "\nTest normalise with more complicated pcf (this should be okay)" << endl;
			int d1 = 5; // dimension of the box to sample data from
			ivector pavingBox1(d1);
			interval pavingInterval1(-2.5, 2.5);
			for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;
		
			PiecewiseConstantFunction pcf(pavingBox1);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			cxsc::real intUnnorm = pcf.getTotalIntegral();
			assert (intUnnorm != 1.0);
			PiecewiseConstantFunction norm = pcf.makeNormalised();
			cxsc::real intNorm = norm.getTotalIntegral();
			assert (intNorm == 1.0);
			cout << "passed normalised asserts" << endl;
							
		}
		cout << "\nEnd of normalise tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test normalise:\n" << msg << endl;
		throw;
	}
		
}


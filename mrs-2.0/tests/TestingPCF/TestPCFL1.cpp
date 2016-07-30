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
\brief Testing PiecewiseConstantFunction L1 distances
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



void testL1()
{
	int prec = 5;
	cout << "\n\nTest L1 distance" << endl;

	try {
		try {
			cout << "\ncall getL1distance for pcf with no subpaving (this should fail)" << endl;
			
			int d = 1; // dimension of the box to sample data from
			interval pavingInterval(-2,2);
			ivector pavingBox(d);
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			PiecewiseConstantFunction pcfNull;
			PiecewiseConstantFunction pcf(pavingBox);
			
			cxsc::real dis = pcfNull.getL1Distance(pcf);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do getL1distance for pcf with no subpaving:\n" << msg << endl;
		}
		try {
			cout << "\ncall getL1distance against pcf with no subpaving (this should fail)" << endl;
			
			int d = 1; // dimension of the box to sample data from
			interval pavingInterval(-2,2);
			ivector pavingBox(d);
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			PiecewiseConstantFunction pcfNull;
			PiecewiseConstantFunction pcf(pavingBox);
			
			cxsc::real dis = pcf.getL1Distance(pcfNull);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do getL1distance against pcf with no subpaving:\n" << msg << endl;
		}
		try {
			cout << "\ncall getL1distance for two pcf with incompatible dimensions (this should fail)" << endl;
			
			int d = 3; // dimension of the box to sample data from
			interval pavingInterval(-2,2);
			ivector pavingBox1(d);
			for(int k=1; k <= d; k++) pavingBox1[k] = pavingInterval;
			ivector pavingBox2(d+1);
			for(int k=1; k <= d+1; k++) pavingBox2[k] = pavingInterval;
			
			PiecewiseConstantFunction pcf1(pavingBox1);
			PiecewiseConstantFunction pcf2(pavingBox2);
			
			cxsc::real dis = pcf1.getL1Distance(pcf2);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do getL1distance for two pcf with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\ncall getL1distance for two pcf with incompatible dimensions (this should fail)" << endl;
			
			int d = 2; // dimension of the box to sample data from
			interval pavingInterval1(-2,2);
			interval pavingInterval2(-1,2);
			ivector pavingBox1(d+1);
			for(int k=1; k <= d+1; k++) pavingBox1[k] = pavingInterval1;
			ivector pavingBox2(d+1);
			for(int k=1; k <= d; k++) pavingBox2[k] = pavingInterval1;
			pavingBox2[d+1] = pavingInterval2;
			
			PiecewiseConstantFunction pcf1(pavingBox1);
			PiecewiseConstantFunction pcf2(pavingBox2);
			
			cxsc::real dis = pcf1.getL1Distance(pcf2);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do getL1distance for two pcf with incompatible dimensions:\n" << msg << endl;
		}
		{
			cout << "\ncall getL1distance against self (this should be okay)" << endl;
			
			int d = 1; // dimension of the box to sample data from
			interval pavingInterval(-2,2);
			ivector pavingBox(d);
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
						
			cxsc::real dis = pcf.getL1Distance(pcf);
			
			cxsc::real shouldBe(0.0);
			assert(dis == shouldBe);
			
			cout << "Passed assert that distance is " << shouldBe << endl;
			
		}
		
		
		int d = 3; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		{
			cout << "\ncall getL1distance against copy of self (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf1(pavingBox);
			std::string split = "1,1";
			pcf1.splitToShape(split);
			std::vector< real > ranges = makeRanges1(pcf1.getDomainVolume());
			pcf1.allocateRanges(ranges);
			
			PiecewiseConstantFunction pcf2(pcf1);
						
			cxsc::real dis = pcf1.getL1Distance(pcf2);
			
			cxsc::real shouldBe(0.0);
			assert(dis == shouldBe);
			
			cout << "Passed assert that distance is " << shouldBe << endl;
			
		}
		{
			cout << "\ncall getL1distance against a pcf with no value (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf1(pavingBox);
			std::string split = "1,1";
			pcf1.splitToShape(split);
			std::vector< real > ranges = makeRanges1(pcf1.getDomainVolume());
			pcf1.allocateRanges(ranges);
			
			PiecewiseConstantFunction pcf2(pavingBox);
						
			cxsc::real dis1 = pcf1.getL1Distance(pcf2);
			
			cxsc::real shouldBe = pcf1.getTotalIntegral();
			assert(dis1 == shouldBe);
			
			cout << "Passed assert that distance pcf1.getL1Distance(pcf2) is " << shouldBe << endl;
			
			cxsc::real dis2 = pcf2.getL1Distance(pcf1);
			
			assert(dis2 == shouldBe);
			
			cout << "Passed assert that distance pcf2.getL1Distance(pcf1) is " << shouldBe << endl;
			
		}
		{
			cout << "\ncall getL1distance against a pcf with negative values (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf1(pavingBox);
			std::string split = "1,1";
			pcf1.splitToShape(split);
			std::vector< real > ranges1 = makeRanges1(pcf1.getDomainVolume());
			pcf1.allocateRanges(ranges1);
			
			PiecewiseConstantFunction pcf2(pavingBox);
			pcf2.splitToShape(split);
			std::vector< real > ranges2 = makeRangesNegative(pcf2.getDomainVolume());
			pcf2.allocateRanges(ranges2);
						
			cxsc::real dis1 = pcf1.getL1Distance(pcf2);
			
			cxsc::real shouldBe = (cxsc::abs(cxsc::abs(ranges1[1]) - cxsc::abs(ranges2[1]))
						+cxsc::abs(cxsc::abs(ranges1[2]) - cxsc::abs(ranges2[2])))/(cxsc::abs(ranges1[1]) + cxsc::abs(ranges1[2]));
			
			cout << "Distance should be " << shouldBe << endl;
			
			assert(dis1 == shouldBe);
			
			cout << "Passed assert that distance pcf1.getL1Distance(pcf2) is " << shouldBe << endl;
			
			cxsc::real dis2 = pcf2.getL1Distance(pcf1);
			
			assert(dis2 == shouldBe);
			
			cout << "Passed assert that distance pcf2.getL1Distance(pcf1) is " << shouldBe << endl;
			
		}
		{
			cout << "\ncall getL1distance against a pcf with infinite values (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf1(pavingBox);
			std::string split = "1,1";
			pcf1.splitToShape(split);
			std::vector< real > ranges1 = makeRanges1(pcf1.getDomainVolume());
			pcf1.allocateRanges(ranges1);
			
			PiecewiseConstantFunction pcf2(pavingBox);
			pcf2.splitToShape(split);
			std::vector< real > ranges2 = makeRangesInfinite();
			pcf2.allocateRanges(ranges2);
						
			cxsc::real dis1 = pcf1.getL1Distance(pcf2);
			
			cxsc::real shouldBe = cxsc::Infinity;
			
			cout << "Distance should be " << shouldBe << endl;
			
			assert(dis1 ==shouldBe);
			
			cout << "Passed assert that distance pcf1.getL1Distance(pcf2) is " << shouldBe << endl;
			
			cxsc::real dis2 = pcf2.getL1Distance(pcf1);
			
			assert(dis2 == shouldBe);
			
			cout << "Passed assert that distance pcf2.getL1Distance(pcf1) is " << shouldBe << endl;
			
		}
		
		{
			cout << "\ncall getL1distance more complex (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf1(pavingBox);
			std::string split1 = "3,4,4,2,2,4,4,3";
			pcf1.splitToShape(split1);
			std::vector< real > ranges1 = makeRanges3(0.5*pcf1.getDomainVolume());
			pcf1.allocateRanges(ranges1);
			
			PiecewiseConstantFunction pcf2(pavingBox);
			std::string split2 = "2,3,4,4,2,3,3";
			pcf2.splitToShape(split2);
			std::vector< real > ranges2 = makeRanges4(pcf2.getDomainVolume());
			pcf2.allocateRanges(ranges2);
			
			cxsc::real shouldBe(10.0/4.0);
						
			cxsc::real dis1 = pcf1.getL1Distance(pcf2);
			assert(dis1 == shouldBe);
			cout << "Distance pcf1.getL1Distance(pcf2) is " << shouldBe << endl;
			
			cxsc::real dis2 = pcf2.getL1Distance(pcf1);
			
			assert(dis2 == dis1);
			
			cout << "Passed assert that distance pcf1.getL1Distance(pcf2) == pcf2.getL1Distance(pcf1) " << endl;
			
		}
		cout << "\nEnd of L1 tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test L1:\n" << msg << endl;
		throw;
	}
		
}


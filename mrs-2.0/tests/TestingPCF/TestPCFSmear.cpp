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
\brief Testing PiecewiseConstantFunction smearing zero values
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


void testSmear()
{
	int prec = 5;
	try {
		try {
			cout << "\ncall smear on pcf with no subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf;
			
			real totalSmear(1.0);
			
			pcf.smearZeroValues(totalSmear);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to call smear on pcf with no subpaving:\n" << msg << endl;
		}
		
		
		int d = 3; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		try {
			cout << "\nmakeMarginal with pcf with no values (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			real totalSmear(1.0);
			
			pcf.smearZeroValues(totalSmear);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::runtime_error& re) {
			std::string msg(re.what());
			cout << "\nFailed to call smear on pcf with no values:\n" << msg << endl;
		}
		try {
			cout << "\ncall smear with total smear = 0 (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			
			real totalSmear(0.0);
			
			pcf.smearZeroValues(totalSmear);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call smear with total smear = 0:\n" << msg << endl;
		}
		try {
			cout << "\ncall smear with total smear >= integral = 1.0(this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "3,4,4,2,2,4,4,3";;
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			
			real totalSmear = pcf.getTotalIntegral();
			
			pcf.smearZeroValues(totalSmear);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call smear with total smear >= integral (this should fail):\n" << msg << endl;
		}
		try {
			cout << "\ncall smear with total smear >= integral > 1.0(this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			
			real totalSmear = pcf.getTotalIntegral();
			
			pcf.smearZeroValues(totalSmear);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call smear with total smear >= integral (this should fail):\n" << msg << endl;
		}
		
		{
			cout << "\ncall smear with negative ranges (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "2,2,2,2";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges6(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			cout << "Before smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			real totalSmear = 0.1*pcf.getTotalIntegral();
			
			pcf.smearZeroValues(totalSmear);
			cout << "After smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
		}
		
		{
			cout << "\ncall smear with positive ranges (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			cout << "Before smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			real totalSmear = 0.1*pcf.getTotalIntegral();
			
			pcf.smearZeroValues(totalSmear);
			cout << "After smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
		}
		
		
		{
			
			cout << "\nTest smear with more complicated pcf (this should be okay)" << endl;
			int d1 = 5; // dimension of the box to sample data from
			ivector pavingBox1(d1);
			interval pavingInterval1(-2.5, 2.5);
			for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;
		
			PiecewiseConstantFunction pcf(pavingBox1);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			cout << "Before smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			real totalSmear = 0.1*pcf.getTotalIntegral();
			
			pcf.smearZeroValues(totalSmear);
			cout << "After smear, integral is " << pcf.getTotalIntegral() 
							<< " and pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
							
		}
		cout << "\nEnd of smear tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test marginal:\n" << msg << endl;
		throw;
	}
		
}



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
\brief Testing PiecewiseConstantFunction IAE calculations
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


void testIAE()
{
	int prec = 5;
	cout << "\n\nTest IAE" << endl;

	try {
		
		int d = 2; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		int label1 = 1;
		int label2 = 2;
		
		PiecewiseConstantFunction pcf1(pavingBox, label1);
		std::string split1 = "3,4,4,2,2,4,4,3";
		pcf1.splitToShape(split1);
		std::vector< real > ranges1 = makeRanges3(0.5*pcf1.getDomainVolume());
		pcf1.allocateRanges(ranges1);
		
		string s1("pcfArithmetic1.txt");
		pcf1.outputToTxtTabs(s1, prec, true);
				
		ivector pavingBoxW(d+1);
		for(int k=1; k <= d+1; k++) pavingBoxW[k] = pavingInterval;
		PiecewiseConstantFunction pcfWrongD(pavingBoxW, label1);
		
		ivector pavingBoxI(d);
		interval pavingIntervalW(-2,1);
		pavingBoxI[1] = pavingIntervalW;
		for(int k=2; k <= d; k++) pavingBoxI[k] = pavingInterval;
		PiecewiseConstantFunction pcfWrongI(pavingBoxI, label1);
		
		PiecewiseConstantFunction pcf2(pavingBox, label2);
		std::string split2 = "2,3,4,4,2,3,3";
		pcf2.splitToShape(split2);
		std::vector< real > ranges2 = makeRanges5(pcf2.getDomainVolume());
		pcf2.allocateRanges(ranges2);
		
		string s2("pcfArithmetic2.txt");
		pcf2.outputToTxtTabs(s2, prec, true);
		
		string additionLeafString("3,4,4,3,4,4,2,4,4,3");
		
		
		// IAE fail tests
		
		
		try {
			cout << "\nIAE with this as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			cxsc::real IAE = pcfNull.getIAE(pcf1);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do IAE with this as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nIAE with operand as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			cxsc::real IAE = pcf1.getIAE(pcfNull);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do IAE with operand as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nIAE with operand with incompatible dimensions (this should fail)" << endl;
			
			cxsc::real IAE = pcf1.getIAE(pcfWrongD);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do IAE with operand with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			cxsc::real IAE = pcf1.getIAE(pcfWrongI);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do IAE with operand with incompatible dimensions:\n" << msg << endl;
		}
		
		
		// IAE tests (should succeed)
		{
			
			
			string before;
			string after;
				
			{
				cout << "\nIAE with pcf with no values " << endl;
				
				PiecewiseConstantFunction rhs(pavingBox);
				
				real intBefore = pcf1.getTotalIntegral();
				assert(rhs.getTotalIntegral() == real(0.0));
				
				cxsc::real IAE = pcf1.getIAE(rhs);
				
				assert(IAE == intBefore);
				cout << "Passed asserts for IAE with pcf with no values" << endl;
				
			}
			{
				cout << "\nIAE for pcf with no values " << endl;
				
				PiecewiseConstantFunction lhs(pavingBox);
				
				real intBefore = pcf1.getTotalIntegral();
				assert(lhs.getTotalIntegral() == real(0.0));
				
				cxsc::real IAE = lhs.getIAE(pcf1);
				
				assert(IAE == intBefore);
				cout << "Passed asserts for IAE for pcf with no values" << endl;
				
			}
			{
				cout << "\nIAE for pcf against self " << endl;
				
				cxsc::real IAE = pcf1.getIAE(pcf1);
				
				assert(IAE == 0.0);
				cout << "Passed asserts for IAE for pcf against self" << endl;
				
			}
			{
				cout << "\nIAE for pcf against copy of self " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				cxsc::real IAE = pcf1.getIAE(pcf);
				
				assert(IAE == 0.0);
				cout << "Passed asserts for IAE for pcf against copy of self" << endl;
				
			}
			{
				cout << "\nIAE with pcf with values " << endl;
				
				cxsc::real IAE1_12 = pcf1.getIAE(pcf2);
				
				cout << "pcf1.getIAE(pcf2) = " << IAE1_12 << endl;
				
				PiecewiseConstantFunction pcf11(pcf1);
				pcf11 -= pcf2;
				assert(pcf11.getTotalIntegral() == IAE1_12);
				
				cout << "\nDifference is " << endl;
				pcf11.outputRootToStreamTabs(cout, prec);
				
				cxsc::real IAE1_21 = pcf2.getIAE(pcf1);
				assert(IAE1_12 == IAE1_21);
				
				cout << "\npcf2.getIAE(pcf1) = " << IAE1_21 << endl;
				
				PiecewiseConstantFunction pcf22(pcf2);
				pcf22 -= pcf1;
				assert(pcf22.getTotalIntegral() == IAE1_12);
				
				cout << "Passed asserts for subtraction to self with pcf with values" << endl;
									
			}
			
		}
		
		

		cout << "\nEnd of IAE tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test IAE:\n" << msg << endl;
		throw;
	}
		
}



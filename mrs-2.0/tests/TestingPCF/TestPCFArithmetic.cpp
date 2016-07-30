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
\brief Testing PiecewiseConstantFunction arithmetic
 */

#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
#include "subpaving_exception.hpp"
#include "intervalmappedspnode.hpp"
#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>
#include <cfloat> // for DBL_EPSILON


using namespace cxsc;
using namespace std;
using namespace subpavings;


void testArithmetic()
{
	int prec = 5;
	cout << "\n\nTest arithmetic" << endl;

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
		
		// addition fail tests
		try {
			cout << "\nAddition to self as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcfNull += pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do addition to self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nAddition to self with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcf1 += pcfNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do addition to self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nAddition to self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 += pcfWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do addition to self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nAddition to self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 += pcfWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do addition to self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nAddition with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			PiecewiseConstantFunction temp = pcfNull + pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do addition with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nAddition of self as null subpaving with scalar add(this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			cxsc::real add(0.0);
			
			pcfNull += add;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do addition as self as null subpaving with scalar add:\n" << msg << endl;
		}
		
		// subtraction fail tests
		try {
			cout << "\nSubtraction from self as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcfNull -= pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do subtraction from self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcf1 -= pcfNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do subtraction from self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 -= pcfWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do subtraction from self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 -= pcfWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do subtraction from self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			PiecewiseConstantFunction temp = pcfNull - pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do subtraction with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self as null subpaving with scalar sub(this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			cxsc::real sub(0.0);
			
			pcfNull -= sub;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do subtraction from self as null subpaving with scalar sub:\n" << msg << endl;
		}
		
		// multiplication fail tests
		try {
			cout << "\nMultiplication of self as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcfNull *= pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do multiplication of self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nMultiplication of self with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcf1 *= pcfNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to multiplication of self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nMultiplication of self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 *= pcfWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do multiplication of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 *= pcfWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do multiplication of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nMultiplication with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			PiecewiseConstantFunction temp = pcfNull * pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do multiplication with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nMultiplication of self as null subpaving with scalar mult(this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			cxsc::real mult(0.0);
			
			pcfNull *= mult;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do multiplication as self as null subpaving with scalar mult:\n" << msg << endl;
		}
		
		// division fail tests
		try {
			cout << "\nDivision of self as null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcfNull /= pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division of self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision of self with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			pcf1 /= pcfNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to division of self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision of self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 /= pcfWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do division of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			pcf1 /= pcfWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do division of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nDivision with null subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			
			PiecewiseConstantFunction temp = pcfNull / pcf1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision of self as null subpaving with scalar div (this should fail)" << endl;
			
			PiecewiseConstantFunction pcfNull;
			cxsc::real div(1.0);
			
			pcfNull /= div;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division as self as null subpaving with scalar div:\n" << msg << endl;
		}
		try {
			cout << "\nDivision of self with scalar div = 0.0 (this should fail)" << endl;
			
			cxsc::real div(0.0);
			
			pcf1 /= div;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to do division as self with scalar div = 0.0:\n" << msg << endl;
		}
		
		// addition tests (should succeed)
		{
			{
				cxsc::real scalar(0.0);
				
				cout << "\nAddition to self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf += scalar;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			{
				cxsc::real scalar(2.0);
				
				cout << "\nAddition to self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf += scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before + pcf.getDomainVolume()*scalar == after);
				
			}
			{
				int scalar = 2;
				
				cout << "\nAddition to self with integer scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf += scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before + pcf.getDomainVolume()*scalar == after);
				
			}
			string before;
			string after;
				
			{
				cout << "\nAddition to self with pcf with no values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				PiecewiseConstantFunction rhs(pavingBox);
				
				real intBefore = pcf.getTotalIntegral();
				assert(rhs.getTotalIntegral() == real(0.0));
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf += rhs;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getTotalIntegral() == intBefore);
				assert(before == after);
				assert(pcf.getLabel() == label1);
				cout << "Passed asserts for addition to self with pcf with no values" << endl;
				
			}
			{
				cout << "\nAddition to self with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				pcf += pcf2;
			
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getTotalIntegral() == intBefore + pcf2.getTotalIntegral());
				assert(pcf.getLabel() == label1);
				
				string s("pcfAddition.txt");
				pcf.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for addition to self with pcf with values" << endl;
									
			}
			{
				cout << "\nAddition with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf = pcf1 + pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert(pcf.getTotalIntegral() == pcf1.getTotalIntegral() + pcf2.getTotalIntegral());
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == 0);
				cout << "Passed asserts for addition with pcf with values  " << endl;
				
			}
			{
				cout << "\nAddition with pcf with values and same label " << endl;
				
				PiecewiseConstantFunction pcfT(pcf1);
				pcfT.setLabel(label2);
				
				PiecewiseConstantFunction pcf = pcfT + pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert(pcf.getTotalIntegral() == pcfT.getTotalIntegral() + pcf2.getTotalIntegral());
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == label2);
				cout << "Passed asserts for addition with pcf with values and same label " << endl;
				
			}
		}	
		
		// subtraction tests (should succeed)
		{
			{
				cxsc::real scalar(0.0);
				
				cout << "\nSubtraction from self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf -= scalar;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			{
				cxsc::real scalar(2.0);
				
				cout << "\nSubtraction from self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf -= scalar;
				
				real after = pcf.getTotalIntegral();
				
				//This only works because all ranges were positive and all go negative
				assert(cxsc::abs(before - pcf.getDomainVolume()*scalar) == after);
				
			}
			{
				int scalar = 2;
				
				cout << "\nSubtraction from self with integer scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf -= scalar;
				
				real after = pcf.getTotalIntegral();
				
				//This only works because all ranges were positive and all go negative
				assert(cxsc::abs(before - pcf.getDomainVolume()*scalar) == after);
				
			}
			string before;
			string after;
				
			{
				cout << "\nSubtraction from self with pcf with no values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				PiecewiseConstantFunction rhs(pavingBox);
				
				real intBefore = pcf.getTotalIntegral();
				assert(rhs.getTotalIntegral() == real(0.0));
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf -= rhs;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getTotalIntegral() == intBefore);
				assert(before == after);
				assert(pcf.getLabel() == label1);
				cout << "Passed asserts for subtraction to self with pcf with no values" << endl;
				
			}
			{
				cout << "\nSubtraction from self with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				pcf -= pcf2;
			
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getLabel() == label1);
				
				string s("pcfSubtraction.txt");
				pcf.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for subtraction to self with pcf with values" << endl;
									
			}
			{
				cout << "\nSubtraction with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf = pcf1 - pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == 0);
				cout << "Passed asserts for subtraction with pcf with values  " << endl;
				
			}
			{
				cout << "\nSubtraction with pcf with values and same label " << endl;
				
				PiecewiseConstantFunction pcfT(pcf1);
				pcfT.setLabel(label2);
				
				PiecewiseConstantFunction pcf = pcfT - pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == label2);
				cout << "Passed asserts for subtraction with pcf with values and same label " << endl;
				
			}
		}
		
		// multiplication tests (should succeed)
		{
			{
				cxsc::real scalar(0.0);
				
				cout << "\nMultiplication of self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				assert(pcf.getTotalIntegral() > cxsc::real(0.0));
				
				pcf *= scalar;
				
				assert(pcf.getTotalIntegral() == cxsc::real(0.0));
				
			}
			{
				cxsc::real scalar(1.0);
				
				cout << "\nMultiplication of self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf *= scalar;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			
			{
				cxsc::real scalar(2.0);
				
				cout << "\nMultiplication of self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf *= scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before*scalar == after);
				
			}
			{
				int scalar = 2;
				
				cout << "\nMultiplication of self with integer scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf *= scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before*scalar == after);
				
			}
			string before;
			string after;
				
			{
				cout << "\nMultiplication of self with pcf with no values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				PiecewiseConstantFunction rhs(pavingBox);
				
				real intBefore = pcf.getTotalIntegral();
				assert(rhs.getTotalIntegral() == real(0.0));
				
				pcf *= rhs;
				
				assert(pcf.getTotalIntegral() == real(0.0));
				assert(pcf.getLabel() == label1);
				cout << "Passed asserts for multiplication to self with pcf with no values" << endl;
				
			}
			{
				cout << "\nMultiplication of self with pcf with value 1.0 " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				RealMappedSPnode rmspn(pavingBox, real(1.0));
				PiecewiseConstantFunction rhs(rmspn);
				
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf *= rhs;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getTotalIntegral() == intBefore);
				assert(before == after);
				assert(pcf.getLabel() == label1);
				cout << "Passed asserts for multiplication to self with pcf with value 1.0" << endl;
				
			}
			{
				cout << "\nMultiplication of self with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				pcf *= pcf2;
			
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getLabel() == label1);
				
				string s("pcfMultiplication.txt");
				pcf.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for multiplication to self with pcf with values" << endl;
									
			}
			{
				cout << "\nMultiplication with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf = pcf1 * pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == 0);
				cout << "Passed asserts for multiplication with pcf with values  " << endl;
				
			}
			{
				cout << "\nMultiplication with pcf with values and same label " << endl;
				
				PiecewiseConstantFunction pcfT(pcf1);
				pcfT.setLabel(label2);
				
				PiecewiseConstantFunction pcf = pcfT * pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == label2);
				cout << "Passed asserts for multiplication with pcf with values and same label " << endl;
				
			}
		}
		
		// Division tests (should succeed)
		{
			
			{
				cxsc::real scalar(1.0);
				
				cout << "\nDivision of self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf /= scalar;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			
			{
				cxsc::real scalar(2.0);
				
				cout << "\nDivision of self with scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf /= scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before/scalar == after);
				
			}
			{
				int scalar = 2;
				
				cout << "\nDivision of self with integer scalar " << scalar << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				
				real before = pcf.getTotalIntegral();
				
				pcf /= scalar;
				
				real after = pcf.getTotalIntegral();
				
				assert(before/scalar == after);
				
			}
			string before;
			string after;
				
			
			{
				cout << "\nDivision of self with pcf with value 1.0 " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				RealMappedSPnode rmspn(pavingBox, real(1.0));
				PiecewiseConstantFunction rhs(rmspn);
				
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				pcf /= rhs;
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getTotalIntegral() == intBefore);
				assert(before == after);
				assert(pcf.getLabel() == label1);
				cout << "Passed asserts for division to self with pcf with value 1.0" << endl;
				
			}
			{
				cout << "\nDivision of self with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf(pcf1);
				real intBefore = pcf.getTotalIntegral();
				
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				pcf /= pcf2;
			
				{
					ostringstream oss;
					pcf.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(pcf.getLabel() == label1);
				
				string s("pcfDivision1.txt");
				pcf.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for division to self with pcf with values" << endl;
									
			}
			{
				cout << "\nDivision of self with pcf with values that include 0's" << endl;
				
				PiecewiseConstantFunction pcf(pcf2);
				real intBefore = pcf.getTotalIntegral();
				
				pcf /= pcf1;
			
				cout << "After division pcf is " << endl;
				pcf.outputRootToStreamTabs(cout, prec);
				
				std::vector < int > reqDims;
				reqDims.push_back(1);
				PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			
				cout << "marg is" << endl;
				marg.outputRootToStreamTabs(cout, prec);
			
				assert(pcf.hasInfinitePiecewiseConstantValues());
				assert(pcf.getTotalIntegral() == cxsc::Infinity);
				
				assert(pcf.getLabel() == label2);
				
				string s("pcfDivision2.txt");
				pcf.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for division to self with pcf with values that include 0's" << endl;
									
			}
			{
				cout << "\nDivision with pcf with values " << endl;
				
				PiecewiseConstantFunction pcf = pcf1 / pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == 0);
				cout << "Passed asserts for division with pcf with values  " << endl;
				
			}
			{
				cout << "\nDivision with pcf with values and same label " << endl;
				
				PiecewiseConstantFunction pcfT(pcf1);
				pcfT.setLabel(label2);
				
				PiecewiseConstantFunction pcf = pcfT / pcf2;
				
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (pcf.getLeafLevelsString()) == additionLeafString);
				assert(afterT == after);
				assert(pcf.getLabel() == label2);
				cout << "Passed asserts for division with pcf with values and same label " << endl;
				
			}
		}
		cout << "\nEnd of arithmetic tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test arithmetic:\n" << msg << endl;
		throw;
	}
		
}



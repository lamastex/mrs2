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
\brief Testing SpatialObjectRepresentationBV arithmetic
 */

#include "TestSORTools.hpp"

#include "spatial_object_representation_bv.hpp"
#include "subpaving_exception.hpp"
#include "booleanvaluemappedspnode.hpp"
#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>
#include <cfloat> // for DBL_EPSILON


using namespace cxsc;
using namespace std;
using namespace subpavings;


void testBVArithmetic()
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
		
		SpatialObjectRepresentationBV sor1(pavingBox, label1);
		std::string split1 = "3,4,4,2,2,4,4,3";
		sor1.splitToShape(split1);
		std::vector< bool > ranges1 = makeRanges3();
		sor1.allocateRanges(ranges1);
		
		string s1("sorArithmetic1.txt");
		sor1.outputToTxtTabs(s1, prec, true);
				
		ivector pavingBoxW(d+1);
		for(int k=1; k <= d+1; k++) pavingBoxW[k] = pavingInterval;
		SpatialObjectRepresentationBV sorWrongD(pavingBoxW, label1);
		
		ivector pavingBoxI(d);
		interval pavingIntervalW(-2,1);
		pavingBoxI[1] = pavingIntervalW;
		for(int k=2; k <= d; k++) pavingBoxI[k] = pavingInterval;
		SpatialObjectRepresentationBV sorWrongI(pavingBoxI, label1);
		
		SpatialObjectRepresentationBV sor2(pavingBox, label2);
		std::string split2 = "2,3,4,4,2,3,3";
		sor2.splitToShape(split2);
		std::vector< bool > ranges2 = makeRanges5();
		sor2.allocateRanges(ranges2);
		
		string s2("sorArithmetic2.txt");
		sor2.outputToTxtTabs(s2, prec, true);
		
		string unionLeafString("3,4,4,3,4,4,2,4,4,3");
		
		// union fail tests
		try {
			cout << "\nUnion with self as null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sorNull += sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do union with self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nUnion with self with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sor1 += sorNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do union with self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nUnion with self with incompatible dimensions (this should fail)" << endl;
			
			sor1 += sorWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do union with self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nUnion with self with incompatible dimensions (this should fail)" << endl;
			
			sor1 += sorWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do union with self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nUnion with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			SpatialObjectRepresentationBV temp = sorNull + sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do union with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nUnion of self as null subpaving with scalar add(this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			bool add =false;
			
			sorNull += add;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do union as self as null subpaving with scalar add:\n" << msg << endl;
		}
		
		// subtraction (XOR symmetric set difference) tests
		try {
			cout << "\nXOR of self as null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sorNull -= sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do XOR of self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nXOR of self with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sor1 -= sorNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to XOR of self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nXOR of self with incompatible dimensions (this should fail)" << endl;
			
			sor1 -= sorWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do XOR of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			sor1 -= sorWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do XOR of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nXOR with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			SpatialObjectRepresentationBV temp = sorNull / sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do XOR with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nXOR of self as null subpaving with scalar div (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			bool div = false;
			
			sorNull /= div;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do XOR as self as null subpaving with scalar div:\n" << msg << endl;
		}
		
		
		// intersection fail tests
		try {
			cout << "\nIntersection of self as null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sorNull *= sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do intersection of self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nIntersection of self with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sor1 *= sorNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to intersection of self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nIntersection of self with incompatible dimensions (this should fail)" << endl;
			
			sor1 *= sorWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do intersection of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nSubtraction from self with incompatible dimensions (this should fail)" << endl;
			
			sor1 *= sorWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do intersection of self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nIntersection with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			SpatialObjectRepresentationBV temp = sorNull * sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do intersection with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nIntersection of self as null subpaving with scalar mult(this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			bool mult = true;
			
			sorNull *= mult;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do intersection as self as null subpaving with scalar mult:\n" << msg << endl;
		}
		
		// division (set difference) fail tests
		try {
			cout << "\nDivision from self as null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sorNull /= sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division from self as null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision from self with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			sor1 /= sorNull;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division from self with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision from self with incompatible dimensions (this should fail)" << endl;
			
			sor1 /= sorWrongD;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do division from self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nDivision from self with incompatible dimensions (this should fail)" << endl;
			
			sor1 /= sorWrongI;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::IncompatibleDimensions_Error& ice) {
			std::string msg(ice.what());
			cout << "\nFailed to do division from self with incompatible dimensions:\n" << msg << endl;
		}
		try {
			cout << "\nDivision with null subpaving (this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			
			SpatialObjectRepresentationBV temp = sorNull / sor1;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division with null subpaving:\n" << msg << endl;
		}
		try {
			cout << "\nDivision from self as null subpaving with scalar sub(this should fail)" << endl;
			
			SpatialObjectRepresentationBV sorNull;
			bool sub = false;
			
			sorNull /= sub;
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to do division from self as null subpaving with scalar sub:\n" << msg << endl;
		}
		
		
		// union tests (should succeed)
		{
			{
				bool scalar = false;
				
				cout << "\nUnion with self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor += scalar;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			{
				bool scalar = true;
				
				cout << "\nUnion to self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				real before = sor.getTotalVolume();
				
				sor += scalar;
				
				real after = sor.getTotalVolume();
				
				assert(sor.getDomainVolume() == after);
				
			}
			
			
			string before;
			string after;
				
			{
				cout << "\nUnion with self with sor with no values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				SpatialObjectRepresentationBV rhs(pavingBox);
				
				real volBefore = sor.getTotalVolume();
				assert(rhs.getTotalVolume() == real(0.0));
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor += rhs;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(sor.getTotalVolume() == volBefore);
				assert(before == after);
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for union with self with sor with no values" << endl;
				
			}
			{
				cout << "\nUnion with self with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				cout << "\nVol before " << volBefore<< endl;
				
				sor += sor2;
			
				real volAfter = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				cout << "\nVol after " << volAfter<< endl;
				
				assert(volAfter >= volBefore);
				assert(sor.getLabel() == label1);
				
				string s("sorUnion.txt");
				sor.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for union with self with sor with values" << endl;
									
			}
			{
				cout << "\nUnion with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor = sor1 + sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor = sor1 + sor2 = " << sor.getTotalVolume() << endl;
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == 0);
				cout << "Passed asserts for union with sor with values  " << endl;
				
			}
			{
				cout << "\nUnion with sor with values and same label " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				sorT.setLabel(label2);
				
				SpatialObjectRepresentationBV sor = sorT + sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == label2);
				cout << "Passed asserts for union with sor with values and same label " << endl;
				
			}
		}	
		
		
		// XOR (-, symmetric set difference) tests (should succeed)
		{
			
			{
				bool scalar = false;
				
				cout << "\nXOR of self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor -= scalar;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			
			{
				bool scalar = true;
				
				cout << "\nXOR of self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				real before = sor.getTotalVolume();
				
				sor -= scalar;
				
				real after = sor.getTotalVolume();
				
				assert(after == sor.getDomainVolume() - before);
				
			}
			
			string before;
			string after;
				
			
			{
				cout << "\nXOR of self with sor with value false " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				BooleanValueMappedSPnode nmspn(pavingBox, false);
				SpatialObjectRepresentationBV rhs(nmspn);
				
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor -= rhs;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(sor.getTotalVolume() == volBefore);
				assert(before == after);
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for XOR to self with sor with value 1.0" << endl;
				
			}
			{
				cout << "\nXOR of self with sor with value true " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				BooleanValueMappedSPnode nmspn(pavingBox, true);
				SpatialObjectRepresentationBV rhs(nmspn);
				
				real volBefore = sor.getTotalVolume();
				
				sor -= rhs;
				
				assert(sor.getTotalVolume() == sor.getDomainVolume() - volBefore);
				assert(sor.getRootLeaves() == sor1.getRootLeaves());
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for XOR to self with sor with value 1.0" << endl;
				
			}
			{
				cout << "\nXOR of self with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				sor -= sor2;
			
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor /= sor2 = " << sor.getTotalVolume() << endl;
				
				assert(sor.getLabel() == label1);
				
				string s("sorXOR.txt");
				sor.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for XOR to self with sor with values" << endl;
									
			}
			
			{
				cout << "\nXOR with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor = sor1 - sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == 0);
				cout << "Passed asserts for XOR with sor with values  " << endl;
				
			}
			{
				cout << "\nXOR with sor with values and same label " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				sorT.setLabel(label2);
				
				SpatialObjectRepresentationBV sor = sorT - sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == label2);
				cout << "Passed asserts for XOR with sor with values and same label " << endl;
				
			}
		}
		
		
		// intersection tests (should succeed)
		{
			{
				bool scalar = false;
				
				cout << "\nIntersection of self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				assert(sor.getTotalVolume() > cxsc::real(0.0));
				
				sor *= scalar;
				
				assert(sor.getTotalVolume() == cxsc::real(0.0));
				
			}
			{
				bool scalar = true;
				
				cout << "\nIntersection of self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor *= scalar;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			
			
			string before;
			string after;
				
			{
				cout << "\nIntersection of self with sor with no values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				SpatialObjectRepresentationBV rhs(pavingBox);
				
				real volBefore = sor.getTotalVolume();
				assert(rhs.getTotalVolume() == real(0.0));
				
				sor *= rhs;
				
				assert(sor.getTotalVolume() == real(0.0));
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for intersection to self with sor with no values" << endl;
				
			}
			{
				cout << "\nIntersection of self with sor with value false " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				BooleanValueMappedSPnode bmspn(pavingBox, false);
				SpatialObjectRepresentationBV rhs(bmspn);
				
				sor *= rhs;
				
				assert(sor.getTotalVolume() == 0.0);
				assert(sor.getRootLeaves() == sor1.getRootLeaves());
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for intersection to self with sor with value 1.0" << endl;
				
			}
			{
				cout << "\nIntersection of self with sor with value true " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				BooleanValueMappedSPnode rmspn(pavingBox, true);
				SpatialObjectRepresentationBV rhs(rmspn);
				
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor *= rhs;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(sor.getTotalVolume() == volBefore);
				assert(before == after);
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for intersection to self with sor with value 1.0" << endl;
				
			}
			{
				cout << "\nIntersection of self with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				sor *= sor2;
			
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor1 *= sor2 = " << sor.getTotalVolume() << endl;
				
				assert(sor.getLabel() == label1);
				assert(sor.getTotalVolume() <= sor1.getTotalVolume());
				
				string s("sorIntersection.txt");
				sor.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for intersection to self with sor with values" << endl;
									
			}
			{
				cout << "\nIntersection with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor = sor1 * sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor *= sor2 = " << sor.getTotalVolume() << endl;
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == 0);
				cout << "Passed asserts for intersection with sor with values  " << endl;
				
			}
			{
				cout << "\nIntersection with sor with values and same label " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				sorT.setLabel(label2);
				
				SpatialObjectRepresentationBV sor = sorT * sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == label2);
				cout << "Passed asserts for intersection with sor with values and same label " << endl;
				
			}
		}
		
		// set difference tests / (should succeed)
		{
			{
				bool scalar = false;
				
				cout << "\nSet difference against self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				string before;
				string after;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor /= scalar;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(before == after);
				
			}
			{
				bool scalar = true;
				
				cout << "\nSet difference against self with scalar " << scalar << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				
				real before = sor.getTotalVolume();
				
				sor /= scalar;
				
				real after = sor.getTotalVolume();
				
				assert(after == 0.0);
				
			}
			
			string before;
			string after;
				
			{
				cout << "\nSet difference against self with sor with no values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				SpatialObjectRepresentationBV rhs(pavingBox);
				
				real volBefore = sor.getTotalVolume();
				assert(rhs.getTotalVolume() == real(0.0));
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
				
				sor /= rhs;
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				assert(sor.getTotalVolume() == volBefore);
				assert(before == after);
				assert(sor.getLabel() == label1);
				cout << "Passed asserts for set difference to self with sor with no values" << endl;
				
			}
			{
				cout << "\nSet difference against self with sor with values " << endl;
				
				SpatialObjectRepresentationBV sor(sor1);
				real volBefore = sor.getTotalVolume();
				
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					before = oss.str();
				}
			
				sor /= sor2;
			
				{
					ostringstream oss;
					sor.outputRootToStreamTabs(oss, prec);
					after = oss.str();
				}
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor1 -= sor2 = " << sor.getTotalVolume() << endl;
				
				assert(sor.getLabel() == label1);
								
				string s("sorSetDifference.txt");
				sor.outputToTxtTabs(s, prec, true);
				cout << "Passed asserts for set difference against self with sor with values" << endl;
									
			}
			{
				cout << "\nSet difference against sor with values " << endl;
				
				SpatialObjectRepresentationBV sor = sor1 / sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
				
				cout << "\nVol sor1 = " << sor1.getTotalVolume() << endl;
				cout << "\nVol sor2 = " << sor2.getTotalVolume() << endl;
				cout << "\nVol sor = sor1 - sor2 = " << sor.getTotalVolume() << endl;
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == 0);
				cout << "Passed asserts for set difference with sor with values  " << endl;
				
			}
			{
				cout << "\nSet difference against sor with values and same label " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				sorT.setLabel(label2);
				
				SpatialObjectRepresentationBV sor = sorT / sor2;
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string afterT = oss.str();
							
				assert( (sor.getLeafLevelsString()) == unionLeafString);
				assert(afterT == after);
				assert(sor.getLabel() == label2);
				cout << "Passed asserts for set difference with sor with values and same label " << endl;
				
			}
		}
		
		{
			{
				cout << "\nXOR compared to union and intersection " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				
				
				
				SpatialObjectRepresentationBV sorXOR = sorT - sor2;
				
				SpatialObjectRepresentationBV sorUnion = sor1 + sor2;
				SpatialObjectRepresentationBV sorIntersection = sor1 * sor2;
				
				assert( sorXOR.getTotalVolume() == sorUnion.getTotalVolume() - sorIntersection.getTotalVolume());
				cout << "Passed asserts for XOR compared to union and intersection " << endl;
				
			}
			{
				cout << "\nset diff compared to original and intersection " << endl;
				
				SpatialObjectRepresentationBV sorT(sor1);
				
				
				
				SpatialObjectRepresentationBV sorSD = sorT / sor2;
				
				SpatialObjectRepresentationBV sorIntersection = sor1 * sor2;
				
				assert( sorSD.getTotalVolume() == sor1.getTotalVolume() - sorIntersection.getTotalVolume());
				cout << "Passed asserts for set diff compared to original and intersection " << endl;
				
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



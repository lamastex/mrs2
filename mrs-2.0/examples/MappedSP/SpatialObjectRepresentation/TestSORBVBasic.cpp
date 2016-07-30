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
\brief Testing SpatialObjectRepresentationBV basic
 */

#include "TestSORTools.hpp"

#include "spatial_object_representation_bv.hpp"
#include "subpaving_exception.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>


using namespace cxsc;
using namespace std;
using namespace subpavings;


// test basics for sor
void testBVBasic()
{
	int prec = 5; // default precision for output files
	
					
	try {
		cout << "\nSOR constructed with default box" << endl;
		
		cxsc::ivector uselessBox;

		SpatialObjectRepresentationBV temp(uselessBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct SOR with default box:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nSOR constructed with thin interval box" << endl;
		
		cxsc::ivector thinBox(1.0,2.0);

		SpatialObjectRepresentationBV temp(thinBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct SOR with thin interval box:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nConstruct SOR from spn subpaving with no box (this should fail)" << endl;
		
		BooleanValueMappedSPnode spn;

		SpatialObjectRepresentationBV sor(spn);
		
		throw std::logic_error("Should not be able to get here");
		
	}
	catch (subpavings::NoBox_Error& nbe) {
		cout << "\nException:\n" << nbe.what() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct from spn test:\n" << msg << endl;
		throw;
	}
	
	
	
	try {
		cout << "\nDefault SOR" << endl;
		
		SpatialObjectRepresentationBV sor;
		
		cout << sor.stringSummary() << endl;
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed default constructor test:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTrue volume on default SOR" << endl;
		
		SpatialObjectRepresentationBV sor;
		
		assert(sor.getRootLeaves() == 0);
		
		real volume = sor.getTotalVolume();
		throw std::logic_error("Should not be able to get here");
		
	}
	catch (NullSubpavingPointer_Error& nsp) {
		cout << "\nException:\n" << nsp.what() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed volume on default SOR:\n" << msg << endl;
		throw;
	}
	
	
	
	int d = 2; // dimension
	ivector pavingBox(d);
	interval pavingInterval(0,1);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
	try {
		cout << "\nConstruct just with box, default label" << endl;
		
		SpatialObjectRepresentationBV sor(pavingBox);
		
		assert(sor.getRootLeaves() == 1);
		assert(sor.getLabel() == 0);
		
		ostringstream oss;
		sor.outputRootToStreamTabs(oss, prec);
		string sorRoot = oss.str();
		
		cout << "sor root is:\n" << sorRoot << endl;
		
		real volume = sor.getTotalVolume();
		assert(volume == 0.0);
		cout << "Volume is " << volume << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct just with box test:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nConstruct just with box and label" << endl;
		
		int lab = 1;
		SpatialObjectRepresentationBV sor(pavingBox, lab);
		
		assert(sor.getRootLeaves() == 1);
		assert(sor.getLabel() == lab);
		
		ostringstream oss;
		sor.outputRootToStreamTabs(oss, prec);
		string sorRoot = oss.str();
		
		cout << "sor root is:\n" << sorRoot << endl;
		
		real volume = sor.getTotalVolume();
		assert(volume == 0.0);
		cout << "Volume is " << volume << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct just with box and label:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct just with box and split to shape" << endl;
		
		SpatialObjectRepresentationBV sor(pavingBox);
		
		std::string split = "1,2,2";
		
		sor.splitToShape(split);
		
		ostringstream oss;
		sor.outputRootToStreamTabs(oss, prec);
		string sorRoot = oss.str();
		assert(sor.getRootLeaves() == 3);
		cout << "sor root is:\n" << sorRoot << endl;
		
		real volume = sor.getTotalVolume();
		assert(volume == 0.0);
		cout << "Volume is " << volume << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailedto construct just with box and split to shape:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct just with box and split to shape and allocate ranges" << endl;
		
		SpatialObjectRepresentationBV sor(pavingBox);
		
		std::string split = "1,1";
		
		sor.splitToShape(split);
		std::vector< bool > ranges = makeRanges1();
		
		sor.allocateRanges(ranges);
					
		ostringstream oss;
		sor.outputRootToStreamTabs(oss, prec);
		string sorRoot = oss.str();
		assert(sor.getRootLeaves() == 2);
		cout << "sor root is:\n" << sorRoot << endl;
		
		real volume = sor.getTotalVolume();
		assert(volume == 0.5);
		cout << "Volume is " << volume << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct just with box and split to shape and allocate ranges:\n" << msg << endl;
		throw;
	}
	
	
	try {
		cout << "\nCopy constructor" << endl;
		
		int lab = 4;
		SpatialObjectRepresentationBV sor1(pavingBox, lab);
		assert(sor1.getLabel() == lab);
		
		std::string split = "1,1";
		
		sor1.splitToShape(split);
		std::vector< bool > ranges = makeRanges1();
		
		sor1.allocateRanges(ranges);
					
		string sorRoot1;
		{
			ostringstream oss;
			sor1.outputRootToStreamTabs(oss, prec);
			sorRoot1 = oss.str();
		}
		
		SpatialObjectRepresentationBV sor2(sor1);
		assert(sor2.getLabel() == lab);
		
		string sorRoot2;
		{
			ostringstream oss;
			sor2.outputRootToStreamTabs(oss, prec);
			sorRoot2 = oss.str();
		}
		
		assert(sorRoot1 == sorRoot2);
		cout << "Summmary of sor1 is " << endl;
		cout << sor1.stringSummary() << endl;
		cout << "Summmary of sor2 made as a copy of sor1 is " << endl;
		cout << sor2.stringSummary() << endl;
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to copy constructor test:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest swap" << endl;
		
		int label1 = 1;
		SpatialObjectRepresentationBV sor1(pavingBox, label1);
		std::string split1 = "1,1";
		sor1.splitToShape(split1);
		std::vector< bool > ranges1 = makeRanges1();
		sor1.allocateRanges(ranges1);
		
		int d2 = 1; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval(0,1);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;
	
		int label2 = 2;
		SpatialObjectRepresentationBV sor2(pavingBox2, label2);
		std::string split2 = "2,2,2,2";
		sor2.splitToShape(split2);
		std::vector< bool > ranges2 = makeRanges2();
		sor2.allocateRanges(ranges2);
		
		string sorRoot11;
		string sorRoot21;
		
		{
			ostringstream oss;
			sor1.outputRootToStreamTabs(oss, prec);
			sorRoot11 = oss.str();
		}
		{
			ostringstream oss;
			sor2.outputRootToStreamTabs(oss, prec);
			sorRoot21 = oss.str();
		}
		
		cout << "Summmary of sor1 before swap is " << endl;
		cout << sor1.stringSummary() << endl;
		cout << "Summmary of sor2 before swap is " << endl;
		cout << sor2.stringSummary() << endl;
		
		std::swap(sor1, sor2);
		
		string sorRoot12;
		string sorRoot22;
		
		{
			ostringstream oss;
			sor1.outputRootToStreamTabs(oss, prec);
			sorRoot12 = oss.str();
		}
		{
			ostringstream oss;
			sor2.outputRootToStreamTabs(oss, prec);
			sorRoot22 = oss.str();
		}
		
		assert(sorRoot11 == sorRoot22);
		assert(sorRoot21 == sorRoot12);
		
		cout << "Passed swap assert" << endl;
		
		cout << "Summmary of sor1 after swap is " << endl;
		cout << sor1.stringSummary() << endl;
		cout << "Summmary of sor2 after swap is " << endl;
		cout << sor2.stringSummary() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do swap test\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest assignment operator" << endl;
		
		int label1 = 1;
		SpatialObjectRepresentationBV sor1(pavingBox, label1);
		std::string split1 = "1,1";
		sor1.splitToShape(split1);
		std::vector< bool > ranges1 = makeRanges1();
		sor1.allocateRanges(ranges1);
		
		int d2 = 1; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval(0,1);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;
	
		int label2 = 2;
		SpatialObjectRepresentationBV sor2(pavingBox2, label2);
		std::string split2 = "2,2,2,2";
		sor2.splitToShape(split2);
		std::vector< bool > ranges2 = makeRanges2();
		sor2.allocateRanges(ranges2);
		
		string sorRoot11;
		string sorRoot21;
		
		{
			ostringstream oss;
			sor1.outputRootToStreamTabs(oss, prec);
			sorRoot11 = oss.str();
		}
		{
			ostringstream oss;
			sor2.outputRootToStreamTabs(oss, prec);
			sorRoot21 = oss.str();
		}
		
		string sor1Summary1 = sor1.stringSummary();
		cout << "Summmary of sor1 before assignment is " << endl;
		cout << sor1Summary1 << endl;
		cout << "Summmary of sor2 before assignment is " << endl;
		cout << sor2.stringSummary() << endl;
		
		sor2 = sor1;
		assert(sor2.getLabel() == label1);
		
		string sorRoot22;
		
		{
			ostringstream oss;
			sor2.outputRootToStreamTabs(oss, prec);
			sorRoot22 = oss.str();
		}
		
		assert(sorRoot11 == sorRoot22);
		
		string sor1Summary2 = sor1.stringSummary();
		assert(sor1Summary1 == sor1Summary2);
		cout << "Passed assignment asserts" << endl;
		cout << "Summmary of sor1 after assignment is " << endl;
		cout << sor1Summary2 << endl;
		cout << "Summmary of sor2 after assignment is " << endl;
		cout << sor2.stringSummary() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do assignment test\n" << msg << endl;
		throw;
	}
	cout << "\nEnd of basic tests:\n" << endl;	
}


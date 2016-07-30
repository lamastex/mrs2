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
\brief Testing PiecewiseConstantFunction basic
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


// test basics for pcf
void testBasic()
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
					
	try {
		cout << "\nPCF constructed with default box" << endl;
		
		cxsc::ivector uselessBox;

		PiecewiseConstantFunction temp(uselessBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct PCF with default box:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nPCF constructed with thin interval box" << endl;
		
		cxsc::ivector thinBox(1.0,2.0);

		PiecewiseConstantFunction temp(thinBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct PCF with thin interval box:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nConstruct PCF from spn subpaving with no box (this should fail)" << endl;
		
		RealMappedSPnode spn;

		PiecewiseConstantFunction pcf(spn);
		
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
		cout << "\nDefault PCF" << endl;
		
		PiecewiseConstantFunction pcf;
		
		cout << pcf.stringSummary() << endl;
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed default constructor test:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nIntegral on default PCF" << endl;
		
		PiecewiseConstantFunction pcf;
		
		assert(pcf.getRootLeaves() == 0);
		
		real integral = pcf.getTotalIntegral();
		throw std::logic_error("Should not be able to get here");
		
	}
	catch (NullSubpavingPointer_Error& nsp) {
		cout << "\nException:\n" << nsp.what() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed integral on default PCF:\n" << msg << endl;
		throw;
	}
	
	
	
	int d = 2; // dimension of the box to sample data from
	ivector pavingBox(d);
	interval pavingInterval(0,1);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
	try {
		cout << "\nConstruct just with box, default label" << endl;
		
		PiecewiseConstantFunction pcf(pavingBox);
		
		assert(pcf.getRootLeaves() == 1);
		assert(pcf.getLabel() == 0);
		
		ostringstream oss;
		pcf.outputRootToStreamTabs(oss, prec);
		string pcfRoot = oss.str();
		
		cout << "pcf root is:\n" << pcfRoot << endl;
		
		real integral = pcf.getTotalIntegral();
		assert(integral == 0.0);
		cout << "Integral is " << integral << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct just with box test:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nConstruct just with box and label" << endl;
		
		int lab = 1;
		PiecewiseConstantFunction pcf(pavingBox, lab);
		
		assert(pcf.getRootLeaves() == 1);
		assert(pcf.getLabel() == lab);
		
		ostringstream oss;
		pcf.outputRootToStreamTabs(oss, prec);
		string pcfRoot = oss.str();
		
		cout << "pcf root is:\n" << pcfRoot << endl;
		
		real integral = pcf.getTotalIntegral();
		assert(integral == 0.0);
		cout << "Integral is " << integral << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct just with box and label:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct just with box and split to shape" << endl;
		
		PiecewiseConstantFunction pcf(pavingBox);
		
		std::string split = "1,2,2";
		
		pcf.splitToShape(split);
		
		ostringstream oss;
		pcf.outputRootToStreamTabs(oss, prec);
		string pcfRoot = oss.str();
		assert(pcf.getRootLeaves() == 3);
		cout << "pcf root is:\n" << pcfRoot << endl;
		
		real integral = pcf.getTotalIntegral();
		assert(integral == 0.0);
		cout << "Integral is " << integral << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailedto construct just with box and split to shape:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct just with box and split to shape and allocate ranges" << endl;
		
		PiecewiseConstantFunction pcf(pavingBox);
		
		std::string split = "1,1";
		
		pcf.splitToShape(split);
		std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
		
		pcf.allocateRanges(ranges);
					
		ostringstream oss;
		pcf.outputRootToStreamTabs(oss, prec);
		string pcfRoot = oss.str();
		assert(pcf.getRootLeaves() == 2);
		cout << "pcf root is:\n" << pcfRoot << endl;
		
		real integral = pcf.getTotalIntegral();
		assert(integral == 1.0);
		cout << "Integral is " << integral << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct just with box and split to shape and allocate ranges:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nCopy constructor" << endl;
		
		int lab = 4;
		PiecewiseConstantFunction pcf1(pavingBox, lab);
		assert(pcf1.getLabel() == lab);
		
		std::string split = "1,1";
		
		pcf1.splitToShape(split);
		std::vector< real > ranges = makeRanges1(pcf1.getDomainVolume());
		
		pcf1.allocateRanges(ranges);
					
		string pcfRoot1;
		{
			ostringstream oss;
			pcf1.outputRootToStreamTabs(oss, prec);
			pcfRoot1 = oss.str();
		}
		
		PiecewiseConstantFunction pcf2(pcf1);
		assert(pcf2.getLabel() == lab);
		
		string pcfRoot2;
		{
			ostringstream oss;
			pcf2.outputRootToStreamTabs(oss, prec);
			pcfRoot2 = oss.str();
		}
		
		assert(pcfRoot1 == pcfRoot2);
		cout << "Summmary of pcf1 is " << endl;
		cout << pcf1.stringSummary() << endl;
		cout << "Summmary of pcf2 made as a copy of pcf1 is " << endl;
		cout << pcf2.stringSummary() << endl;
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to copy constructor test:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest swap" << endl;
		
		int label1 = 1;
		PiecewiseConstantFunction pcf1(pavingBox, label1);
		std::string split1 = "1,1";
		pcf1.splitToShape(split1);
		std::vector< real > ranges1 = makeRanges1(pcf1.getDomainVolume());
		pcf1.allocateRanges(ranges1);
		
		int d2 = 1; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval(0,1);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;
	
		int label2 = 2;
		PiecewiseConstantFunction pcf2(pavingBox2, label2);
		std::string split2 = "2,2,2,2";
		pcf2.splitToShape(split2);
		std::vector< real > ranges2 = makeRanges2(pcf2.getDomainVolume());
		pcf2.allocateRanges(ranges2);
		
		string pcfRoot11;
		string pcfRoot21;
		
		{
			ostringstream oss;
			pcf1.outputRootToStreamTabs(oss, prec);
			pcfRoot11 = oss.str();
		}
		{
			ostringstream oss;
			pcf2.outputRootToStreamTabs(oss, prec);
			pcfRoot21 = oss.str();
		}
		
		cout << "Summmary of pcf1 before swap is " << endl;
		cout << pcf1.stringSummary() << endl;
		cout << "Summmary of pcf2 before swap is " << endl;
		cout << pcf2.stringSummary() << endl;
		
		std::swap(pcf1, pcf2);
		
		string pcfRoot12;
		string pcfRoot22;
		
		{
			ostringstream oss;
			pcf1.outputRootToStreamTabs(oss, prec);
			pcfRoot12 = oss.str();
		}
		{
			ostringstream oss;
			pcf2.outputRootToStreamTabs(oss, prec);
			pcfRoot22 = oss.str();
		}
		
		assert(pcfRoot11 == pcfRoot22);
		assert(pcfRoot21 == pcfRoot12);
		
		cout << "Passed swap assert" << endl;
		
		cout << "Summmary of pcf1 after swap is " << endl;
		cout << pcf1.stringSummary() << endl;
		cout << "Summmary of pcf2 after swap is " << endl;
		cout << pcf2.stringSummary() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do swap test\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest assignment operator" << endl;
		
		int label1 = 1;
		PiecewiseConstantFunction pcf1(pavingBox, label1);
		std::string split1 = "1,1";
		pcf1.splitToShape(split1);
		std::vector< real > ranges1 = makeRanges1(pcf1.getDomainVolume());
		pcf1.allocateRanges(ranges1);
		
		int d2 = 1; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval(0,1);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;
	
		int label2 = 2;
		PiecewiseConstantFunction pcf2(pavingBox2, label2);
		std::string split2 = "2,2,2,2";
		pcf2.splitToShape(split2);
		std::vector< real > ranges2 = makeRanges2(pcf2.getDomainVolume());
		pcf2.allocateRanges(ranges2);
		
		string pcfRoot11;
		string pcfRoot21;
		
		{
			ostringstream oss;
			pcf1.outputRootToStreamTabs(oss, prec);
			pcfRoot11 = oss.str();
		}
		{
			ostringstream oss;
			pcf2.outputRootToStreamTabs(oss, prec);
			pcfRoot21 = oss.str();
		}
		
		string pcf1Summary1 = pcf1.stringSummary();
		cout << "Summmary of pcf1 before assignment is " << endl;
		cout << pcf1Summary1 << endl;
		cout << "Summmary of pcf2 before assignment is " << endl;
		cout << pcf2.stringSummary() << endl;
		
		pcf2 = pcf1;
		assert(pcf2.getLabel() == label1);
		
		string pcfRoot22;
		
		{
			ostringstream oss;
			pcf2.outputRootToStreamTabs(oss, prec);
			pcfRoot22 = oss.str();
		}
		
		assert(pcfRoot11 == pcfRoot22);
		
		string pcf1Summary2 = pcf1.stringSummary();
		assert(pcf1Summary1 == pcf1Summary2);
		cout << "Passed assignment asserts" << endl;
		cout << "Summmary of pcf1 after assignment is " << endl;
		cout << pcf1Summary2 << endl;
		cout << "Summmary of pcf2 after assignment is " << endl;
		cout << pcf2.stringSummary() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do assignment test\n" << msg << endl;
		throw;
	}
	cout << "\nEnd of basic tests:\n" << endl;	
}


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
\brief Testing PiecewiseConstantFunction - making a random shape
 */

#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
#include "MCMCPartitionGenerator.hpp"
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


// test random for pcf
void testRandom()
{
	int prec = 5; // default precision for output files
	
	int d = 2; // dimension of the box
	ivector pavingBox1(d);
	ivector pavingBox2(d);
	interval pavingInterval1(0,1);
	interval pavingInterval2(-1,1);
	for(int k=1; k <= d; k++) {
		pavingBox1[k] = pavingInterval1;
		pavingBox2[k] = pavingInterval2;
	}
	
	long unsigned int seed = 1234;
	
	MCMCPartitionGenerator partitioner(seed);
					
	{
		cout << "\nrandom PCF with 1 piece" << endl;
		
		long unsigned int numLeaves = 1;
		
		PiecewiseConstantFunction pcf(numLeaves, partitioner, pavingBox1);
		
		assert(pcf.getRootLeaves() == numLeaves);
		
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "pcf root is:\n" << pcfRoot << endl;
		}
		pcf.normalise();
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "normalised pcf root is:\n" << pcfRoot << endl;
		}
		
		
	}
	
	{
		cout << "\nrandom PCF with 2 pieces" << endl;
		
		long unsigned int numLeaves = 2;
		
		PiecewiseConstantFunction pcf(numLeaves, partitioner, pavingBox1);
		
		assert(pcf.getRootLeaves() == numLeaves);
		
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "pcf root is:\n" << pcfRoot << endl;
		}
		pcf.normalise();
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "normalised pcf root is:\n" << pcfRoot << endl;
		}
	}
	
	{
		cout << "\nrandom PCF with 10 pieces" << endl;
		
		long unsigned int numLeaves = 10;
		
		PiecewiseConstantFunction pcf(numLeaves, partitioner, pavingBox1);
		
		assert(pcf.getRootLeaves() == numLeaves);
		
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "pcf root is:\n" << pcfRoot << endl;
		}
		pcf.normalise();
		{
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			string pcfRoot = oss.str();
			
			cout << "normalised pcf root is:\n" << pcfRoot << endl;
		}
	}
	
	try {
		cout << "\nrandom PCF with 1048576 pieces" << endl;
		
		long unsigned int numLeaves = 1048576;
		
		int d10 = 10; // dimension of the box
		ivector pavingBox10(d10);
		interval pavingIntervalSmall(0,1);
		for(int k=1; k <= d10; k++) {
			pavingBox10[k] = pavingIntervalSmall;
		}
		
		PiecewiseConstantFunction pcf(numLeaves, partitioner, pavingBox10);
		
		assert(pcf.getRootLeaves() == numLeaves);
		
		cout << "total integral is " << pcf.getTotalIntegral() << endl;
		pcf.normalise();
		cout << "after normalise, total integral is " << pcf.getTotalIntegral() << endl;
		
	}
	catch (subpavings::UnfulfillableRequest_Error& ure) {
		cout << "could not do that one: error is " << ure.what() << endl;
	}
	
	try {
		cout << "\nrandom PCF with 1048576 pieces, interval [-5, 5]" << endl;
		
		long unsigned int numLeaves = 1048576;
		
		int d10 = 10; // dimension of the box
		ivector pavingBox10(d10);
		interval pavingIntervalSym(-5.0,5.0);
		for(int k=1; k <= d10; k++) {
			pavingBox10[k] = pavingIntervalSym;
		}
		
		PiecewiseConstantFunction pcf(numLeaves, partitioner, pavingBox10);
		
		assert(pcf.getRootLeaves() == numLeaves);
		
		cout << "total integral is " << pcf.getTotalIntegral() << endl;
		pcf.normalise();
		cout << "after normalise, total integral is " << pcf.getTotalIntegral() << endl;
		
	}
	catch (subpavings::UnfulfillableRequest_Error& ure) {
		cout << "could not do that one: error is " << ure.what() << endl;
	}
	
	
	cout << "\nEnd of  random pcf tests:\n" << endl;	
}


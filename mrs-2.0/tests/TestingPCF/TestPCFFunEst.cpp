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
\brief Testing PiecewiseConstantFunction from function estimates
 */

#include "TestPCFTools.hpp"

#include "functionestimator_interval.hpp"
#include "functionestimator_real.hpp"
#include "piecewise_constant_function.hpp"
#include "subpaving_exception.hpp"
#include "simpleFobj1.hpp"
#include "simpleFobj2.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>
#include <cfloat> // for DBL_EPSILON


using namespace cxsc;
using namespace std;
using namespace subpavings;


void testMakeFromFunctionEstimates()
{
	int d = 2; // dimension of the box to sample data from
	ivector pavingBox(d);
	interval pavingInterval(0,1);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
	string feiRoot;
	string ferRoot;
	string pcf1Root;
	string pcf2Root;

	SimpleFobj2 fobj;

	std::string split = "2,3,4,4,2,3,3";
	
	int prec = 5;

	try {
		
		cout << "\n\nMake a piecewise constant function from a function estimator interval:\n" << endl;
		int lab = 1;
		FunctionEstimatorInterval fei(pavingBox, fobj, lab);
		assert(fei.getLabel() == lab);
		
		try {
			fei.splitToShape(split);
			cout << "String summmary of FunctionEstimatorInterval is: " << endl;
			cout << fei.stringSummary() << endl;
			
			ostringstream oss;
			fei.outputRootToStreamTabs(oss, prec);
			feiRoot = oss.str();
			cout << "\nfei root is:\t" << feiRoot << endl;
			
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do splitToShape:\n" << msg << endl;
		}
		
		try {
			PiecewiseConstantFunction pcf = fei.makePiecewiseConstantFunction();
			
			cout << "String summmary of PiecewiseConstantFunction is: " << endl;
			cout << pcf.stringSummary() << endl;
			
			assert(pcf.getLabel() == lab);
			assert(pcf.getRootLeaves() == fei.getRootLeaves());
			
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			pcf1Root = oss.str();
			
			cout << "pcf root is:\t" << pcf1Root << endl;
			
			string s1 = "PCFfromSimpleFobj2Fei.txt";
			pcf.outputToTxtTabs(s1, prec, true);
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to .makePiecewiseConstantFunction:\n" << msg << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make a piecwise constant function from a function estimator real:\n" << msg << endl;
		throw;
	}
	try {
		
		cout << "\n\nMake a piecewise constant function from a function estimator real:\n" << endl;
		
		int lab = 2;
		FunctionEstimatorReal fer(pavingBox, fobj, lab);
		assert(fer.getLabel() == lab);
				
		try {
			fer.splitToShape(split);
			cout << "String summmary of FunctionEstimatorReal is: " << endl;
			cout << fer.stringSummary() << endl;
			
			ostringstream oss;
			fer.outputRootToStreamTabs(oss, prec);
			ferRoot = oss.str();
			cout << "\nfer root is:\t" << ferRoot << endl;
			
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do splitToShape:\n" << msg << endl;
		}
		
		try {
			PiecewiseConstantFunction pcf = fer.makePiecewiseConstantFunction();
			
			cout << "String summmary of PiecewiseConstantFunction is: " << endl;
			cout << pcf.stringSummary() << endl;
			
			assert(pcf.getLabel() == lab);
			assert(pcf.getRootLeaves() == fer.getRootLeaves());
			
			ostringstream oss;
			pcf.outputRootToStreamTabs(oss, prec);
			pcf2Root = oss.str();
			
			cout << "pcf root is:\t" << pcf2Root << endl;
			assert(pcf2Root == pcf1Root);
			assert(pcf2Root == ferRoot);
			
			cout << "Passed asserts on state of pcf" << endl;
			
			string s1 = "PCFfromSimpleFobj2Fer.txt";
			pcf.outputToTxtTabs(s1, prec, true);
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to .makePiecewiseConstantFunction:\n" << msg << endl;
		}
		cout << "\nEnd of function estimator tests:\n" << endl;	

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make a piecewise constant function from a function estimator real:\n" << msg << endl;
		throw;
	}
	
}


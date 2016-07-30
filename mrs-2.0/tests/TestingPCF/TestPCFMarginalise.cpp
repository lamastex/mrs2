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
\brief Testing PiecewiseConstantFunction marginalising
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


void testMarginalise()
{
	int prec = 5;
	try {
		try {
			cout << "\ncall makeMarginal on pcf with no subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf;
			
			std::vector < int > reqDims;
			reqDims.push_back(1);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to call makeMarginal on pcf with no subpaving:\n" << msg << endl;
		}
		
		
		int d = 3; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		{
			cout << "\nmakeMarginal with pcf with no values (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			assert(pcf.getRootLeaves() == 1);
			std::vector < int > reqDims;
			reqDims.push_back(1);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			assert(pcf.getTotalIntegral() == marg.getTotalIntegral());
			assert(marg.getDimensions() == 1);
			
		}
		
		try {
			cout << "\ncall makeMarginal with empty req dims (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > reqDims;
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeMarginal on pcf with empty req dims:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeMarginal with invalid req dims (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > reqDims;
			reqDims.push_back(1);
			reqDims.push_back(4);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeMarginal on pcf with invalid req dims:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeMarginal with invalid req dims (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > reqDims;
			reqDims.push_back(0);
			reqDims.push_back(3);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeMarginal on pcf with invalid req dims:\n" << msg << endl;
		}
		
		
		{
			std::vector < int > reqDims;
			reqDims.push_back(1);
			reqDims.push_back(3);
			
			cout << "\nTest marginal with negative values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			assert(marg.getDimensions() == reqDims.size());
			assert(pcf.getTotalIntegral() == marg.getTotalIntegral());
							
		}
		{
			std::vector < int > reqDims;
			reqDims.push_back(1);
			reqDims.push_back(3);
			
			cout << "\nTest marginal with infinite values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			assert(pcf.hasInfinitePiecewiseConstantValues());
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
						
			assert(marg.getDimensions() == reqDims.size());
			assert(pcf.getTotalIntegral() == marg.getTotalIntegral());
		}
		
		{
			std::vector < int > reqDims;
			reqDims.push_back(3);
			reqDims.push_back(2);
			
			cout << "\nTest marginal with positive values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			assert(marg.getDimensions() == reqDims.size());
			assert(pcf.getTotalIntegral() == marg.getTotalIntegral());
							
		}
		{
			std::vector < int > reqDims;
			reqDims.push_back(3);
			reqDims.push_back(2);
			reqDims.push_back(1);
			
			cout << "\nTest marginal with all dims required (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			assert(marg.getDimensions() == reqDims.size());
			assert(pcf.getTotalIntegral() == marg.getTotalIntegral());
			
			string pcfRoot;
			{
				ostringstream oss;
				pcf.outputRootToStreamTabs(oss, prec);
				pcfRoot = oss.str();
			}
			string margRoot;
			{
				ostringstream oss;
				marg.outputRootToStreamTabs(oss, prec);
				margRoot = oss.str();
			}
			assert (pcfRoot == margRoot);
		}
		
		{
			
			cout << "\nTest marginal with more complicated pcf (this should be okay)" << endl;
			int d1 = 5; // dimension of the box to sample data from
			ivector pavingBox1(d1);
			interval pavingInterval1(-2.5, 2.5);
			for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;
		
			PiecewiseConstantFunction pcf(pavingBox1);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			std::vector < int > reqDims;
			reqDims.push_back(3);
			reqDims.push_back(2);
			
			PiecewiseConstantFunction marg = pcf.makeMarginal(reqDims);
			assert(marg.getDimensions() == reqDims.size());
			cout << "pcf.getTotalIntegral() = "  << pcf.getTotalIntegral() << endl; 
			cout << "marg.getTotalIntegral() = "  << marg.getTotalIntegral() << endl; 
			
			bool isNotEqual = cxsc::abs(pcf.getTotalIntegral() - marg.getTotalIntegral()) 
			> DBL_EPSILON * cxsc::max(pcf.getTotalIntegral(), marg.getTotalIntegral());
			assert(!isNotEqual);
		
			cout << "passed marginalisation asserts" << endl;
			
			string filename("MargTest");
			{
				ostringstream oss;
				oss << filename << "Before.txt";
				pcf.outputToTxtTabs(oss.str(), prec, true);
			}
			string margRoot;
			{
				ostringstream oss;
				oss << filename << "AfterMargOn" << reqDims[0] << "," << reqDims[1] << ".txt";
				marg.outputToTxtTabs(oss.str(), prec);
			}
							
		}
		cout << "\nEnd of marginalise tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test marginal:\n" << msg << endl;
		throw;
	}
		
}



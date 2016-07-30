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
\brief Testing PiecewiseConstantFunction slices
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

void testSlice()
{
	cout << "\n\nTest slice" << endl;

	int prec = 5;
	try {
		try {
			cout << "\ncall makeSlice on pcf with no subpaving (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf;
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe) {
			std::string msg(nspe.what());
			cout << "\nFailed to call makeSlice on pcf with no subpaving:\n" << msg << endl;
		}
		
		
		int d = 3; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		{
			cout << "\nmakeSlice with pcf with no values (this should be okay)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			assert(pcf.getRootLeaves() == 1);
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			cout << "pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			cout << "slice is " << endl;
			slc.outputRootToStreamTabs(cout, prec);
			
			assert(pcf.getTotalIntegral() == slc.getTotalIntegral());
			assert(slc.getDimensions() == d - sliceDims.size());
			
		}
		
		try {
			cout << "\ncall makeSlice with empty slice dims (this should fail)" << endl;
			
			std::vector < int > sliceDims;
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction pcf(pavingBox);
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with empty slice dims:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with empty slice pts (this should fail)" << endl;
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			std::vector < cxsc::real > slicePts;
			
			PiecewiseConstantFunction pcf(pavingBox);
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with empty slicePts:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with slice dims = all dims (this should fail)" << endl;
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			sliceDims.push_back(2);
			sliceDims.push_back(3);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction pcf(pavingBox);
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with slice dims = all dims:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with incompatible slice dims and slice pts (this should fail)" << endl;
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			sliceDims.push_back(2);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction pcf(pavingBox);
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice incompatible slice dims and slice pts:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with invalid slice dims (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(4);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with invalid slice dims:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with invalid slice dims (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(0);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with invalid slice dims:\n" << msg << endl;
		}
		
		try {
			cout << "\ncall makeSlice with invalid slice pts (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(3.0));
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with invalid slice pts:\n" << msg << endl;
		}
		try {
			cout << "\ncall makeSlice with invalid slice pts (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			sliceDims.push_back(2);
			sliceDims.push_back(3);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(-3.0));
			
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with invalid slice pts:\n" << msg << endl;
		}
		
		try {
			cout << "\ncall makeSlice with duplicate slice dimensions (this should fail)" << endl;
			
			PiecewiseConstantFunction pcf(pavingBox);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(2);
			sliceDims.push_back(2);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
						
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to call makeSlice on pcf with invalid slice pts:\n" << msg << endl;
		}
		
		
		{
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			sliceDims.push_back(3);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			
			cout << "\nTest slice with negative values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesNegative(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasNegativePiecewiseConstantValues());
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			assert(slc.getDimensions() == d - sliceDims.size());
							
		}
		{
			
			cout << "\nTest slice with infinite values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRangesInfinite();
			pcf.allocateRanges(ranges);
			
			assert(pcf.hasInfinitePiecewiseConstantValues());
			
			{
				std::vector < int > sliceDims;
				sliceDims.push_back(1);
				std::vector < cxsc::real > slicePts;
				slicePts.push_back(cxsc::real(0.0));
				
				PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
				assert(slc.getDimensions() == d - sliceDims.size());
			}
			{
				std::vector < int > sliceDims;
				sliceDims.push_back(2);
				std::vector < cxsc::real > slicePts;
				slicePts.push_back(cxsc::real(0.0));
				
				PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
				assert(slc.getDimensions() == d - sliceDims.size());
			}					
		}
		{
			std::vector < int > sliceDims;
			sliceDims.push_back(3);
			sliceDims.push_back(2);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			
			cout << "\nTest slice with positive values in pcf (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			assert(slc.getDimensions() == d - sliceDims.size());
							
		}
		
		{
			
			cout << "\nTest slice with simple example (this should be okay)" << endl;
			PiecewiseConstantFunction pcf(pavingBox);
			std::string split = "1,1";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges1(pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			assert(pcf.getTotalIntegral() == 1.0);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0)); // XR
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			
			cout << "pcf has integral " << pcf.getTotalIntegral() << " and is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			
			cout << "sliced pcf has integral " << slc.getTotalIntegral() << " and is " << endl;
			slc.outputRootToStreamTabs(cout, prec);
			
			//cout << "(ranges[2]/(ranges[1]+ranges[2])) = " << (ranges[2]/(ranges[1]+ranges[2])) << endl;
			
			assert(slc.getDimensions() == d - sliceDims.size());
			assert(slc.getTotalIntegral() == 0.5*(ranges[2]/(ranges[1]+ranges[2])));				
		
			
							
		}
		
		
		{
			
			cout << "\nTest slice with more complicated pcf (this should be okay)" << endl;
			int d1 = 3; // dimension of the box to sample data from
			ivector pavingBox1(d1);
			interval pavingInterval1(-2.5, 2.5);
			for(int k=1; k <= d1; k++) pavingBox1[k] = pavingInterval1;
		
			PiecewiseConstantFunction pcf(pavingBox1);
			std::string split = "3,4,4,2,2,4,4,3";
			pcf.splitToShape(split);
			std::vector< real > ranges = makeRanges3(0.5*pcf.getDomainVolume());
			pcf.allocateRanges(ranges);
			
			std::vector < int > sliceDims;
			sliceDims.push_back(1);
			sliceDims.push_back(2);
			std::vector < cxsc::real > slicePts;
			slicePts.push_back(cxsc::real(0.0));
			slicePts.push_back(cxsc::real(0.0));
			
			
			PiecewiseConstantFunction slc = pcf.makeSlice(sliceDims, slicePts);
			assert(slc.getDimensions() == d1 - sliceDims.size());
			cout << "pcf.getTotalIntegral() = "  << pcf.getTotalIntegral() << endl; 
			cout << "slc.getTotalIntegral() = "  << slc.getTotalIntegral() << endl; 
			
			
			string filename("SliceTest");
			{
				ostringstream oss;
				oss << filename << "Before.txt";
				pcf.outputToTxtTabs(oss.str(), prec, true);
			}
			{
				ostringstream oss;
				oss << filename << "AfterSliceOn" << sliceDims[0]
						<< "(" << _double(slicePts[0]) << ")," << sliceDims[1]
						<< "(" << _double(slicePts[1]) << ").txt";
				slc.outputToTxtTabs(oss.str(), prec, true);
			}
							
		}
		cout << "\nEnd of slice tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test slice:\n" << msg << endl;
		throw;
	}
		
}


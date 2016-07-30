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
\brief Testing PiecewiseConstantFunction making from histograms
 */
#include "TestPCFTools.hpp"

#include "piecewise_constant_function.hpp"
#include "adaptivehistogram.hpp"
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



void testMakeFromHist()
{
	int prec = 5; // default precision for output files
		
	int d = 2; // dimension of the box to sample data from
	ivector pavingBox(d);
	interval pavingInterval(-2,2);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	try {
		
		cout << "\nTest making pcfs from histograms" << endl;
		
		
			
		try {
			cout << "\nTest constructing pcf from histogram with no subpaving (this should fail)" << endl;
			
			AdaptiveHistogram adh;
			PiecewiseConstantFunction pcf(adh);
			
			throw std::logic_error("Should not be able to get here");
		}
		catch (subpavings::NullSubpavingPointer_Error& nspe ) {
			std::string msg(nspe.what());
			cout << "\nFailed to construct pcf from histogram with no subpaving:\n" << msg << endl;
		}		
		{	
			cout << "\nTest constructing pcf from histogram (this should be okay)" << endl;
	
			std::string split1 = "3,4,4,2,2,4,4,3";
			bool holdAllStats = false;
			int lab = 2;
			AdaptiveHistogram hist1(pavingBox, holdAllStats, lab);
			
			assert (hist1.getLabel() == lab);
			hist1.splitToShape(split1);

			RVecData data1;
			data1 = getData1(data1);
			bool successfulInsertionHistFirst = hist1.insertFromRVec(data1);
			if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
			
			assert(successfulInsertionHistFirst);
			
			cout << "Histogram is " << endl;
			hist1.outputToStreamTabs(cout, prec);
			
			PiecewiseConstantFunction pcf(hist1);
			
			assert (pcf.getLabel() == lab);
			
			cxsc::real integral = pcf.getTotalIntegral();
			
			string s("PCFfromHist1.txt");
			pcf.outputToTxtTabs(s, prec, true);
			cout << "pcf.getTotalIntegral() = " << integral << endl;
			assert (integral == cxsc::real(1.0));
			cout << "Piecewise constant function is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
		}
		
		{	
			cout << "\nTest constructing pcf from histogram with smear(this should be okay)" << endl;
	
			std::string split1 = "3,4,4,2,2,4,4,3";
			bool holdAllStats = false;
			int lab = 2;
			AdaptiveHistogram hist1(pavingBox, holdAllStats, lab);
			
			assert (hist1.getLabel() == lab);
			hist1.splitToShape(split1);

			RVecData data1;
			data1 = getData1(data1);
			bool successfulInsertionHistFirst = hist1.insertFromRVec(data1);
			if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
			
			assert(successfulInsertionHistFirst);
			
			
			PiecewiseConstantFunction pcf(hist1);
			
			assert (pcf.getLabel() == lab);
			
			cxsc::real integralBefore = pcf.getTotalIntegral();
			
			cout << "pcf is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			cout << "pcf.getTotalIntegral() = " << integralBefore << endl;
			assert (integralBefore == cxsc::real(1.0));
			
			cxsc::real totalSmear((1.0/hist1.getRootCounter()));
			
			pcf.smearZeroValues(totalSmear);
			
			cxsc::real integralAfter = pcf.getTotalIntegral();
			
			cout << "pcf after smear is " << endl;
			pcf.outputRootToStreamTabs(cout, prec);
			cout << "pcf.getTotalIntegral() = " << integralBefore << endl;
			assert (integralBefore == cxsc::real(1.0));
		}
		
		cout << "\nEnd of make pcf from histogram tests:\n" << endl;		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do make pcf from histogram tests:\n" << msg << endl;
		throw;
	}
}		
	

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
special output for logging intermediate states
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


void testArithmeticSpecial()
{
	int prec = 5;
	cout << "\n\nTest arithmetic special" << endl;

	try {
		
		int d = 1; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		PiecewiseConstantFunction pcf1(pavingBox);
		std::string split1 = "3,3,2,3,3,3,3";
		pcf1.splitToShape(split1);
		std::vector< real > ranges1 = 
					makeRangesArithmeticSpecial1(pcf1.getDomainVolume());
		pcf1.allocateRanges(ranges1);
		
		cout << "Total integral for pcf 1 is " << pcf1.getTotalIntegral() << endl;
		
		string s1("pcfArithmeticSpecial1.txt");
		pcf1.outputToTxtTabs(s1, prec, true);
				
		
		PiecewiseConstantFunction pcf2(pavingBox);
		std::string split2 = "3,4,4,3,3,4,4,4,4,2";
		pcf2.splitToShape(split2);
		std::vector< real > ranges2 = 
					makeRangesArithmeticSpecial2(pcf2.getDomainVolume());
		pcf2.allocateRanges(ranges2);
		
		cout << "Total integral for pcf 2 is " << pcf2.getTotalIntegral() << endl;
				
		string s2("pcfArithmeticSpecial2.txt");
		pcf2.outputToTxtTabs(s2, prec, true);
		

		cout << "\nAddition " << endl;
				
		PiecewiseConstantFunction pcf = pcf1 + pcf2;
		
		string s("pcfArithmeticSpecial1.txt");
		pcf1.outputToTxtTabs(s, prec, true);
		
		cout << "\nEnd of arithmetic tests:\n" << endl;	
	}
		
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to test arithmetic special:\n" << msg << endl;
		throw;
	}
		
}



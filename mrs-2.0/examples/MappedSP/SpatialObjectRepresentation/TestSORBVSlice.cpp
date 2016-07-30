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
\brief Testing SpatialObjectRepresentationBV slice
 */

#include "TestSORTools.hpp"

//#include "functionimage_estimator_booleanvalue.hpp"


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


// test slice for sor
void testBVSlice()
{
	interval pavingInterval(-1,1);
	int d = 3;
	ivector pavingBox(d);
	for (int k = 1; k <= d; ++k) pavingBox[k] = pavingInterval;
	
	int prec = 5; // default precision for output files
	
	cout << "\nConstruct just with box and split to shape and allocate ranges" << endl;
	
	SpatialObjectRepresentationBV sor(pavingBox);
	
	std::string split = "3,3,3,3,3,3,3,3";
	
	sor.splitToShape(split);
	std::vector< bool > ranges = makeRanges4();
	
	sor.allocateRanges(ranges);
	
	{			
		ostringstream oss;
		sor.outputRootToStreamTabs(oss, prec);
		string sorRoot = oss.str();
		cout << "sor root is:\n" << sorRoot << endl;
	
		real volume = sor.getTotalVolume();
		cout << "Volume is " << volume << endl;
	}
	
	vector < int > sliceDims(1,3);  // dim 3
	
	vector < real > slicePts(1,-0.5);  
	
	SpatialObjectRepresentationBV slice = sor.makeSlice(sliceDims, slicePts);
	
	{			
		ostringstream oss;
		slice.outputRootToStreamTabs(oss, prec);
		string sliceRoot = oss.str();
		cout << "slice root is:\n" << sliceRoot << endl;
	
		real volume = slice.getTotalVolume();
		cout << "Slice volume is " << volume << endl;
	}
	
	cout << "\nEnd of slice tests:\n" << endl;	
}


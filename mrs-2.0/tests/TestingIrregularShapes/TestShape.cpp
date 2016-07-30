/*
* Copyright (C) 2011 Jennifer Harlow
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
\brief Testing RealMappedShapeNode and PolygonPiecewiseConstantFunction
 */

#include "polygon_piecewise_function.hpp"
#include "piecewise_constant_function.hpp"
#include "realmappedshapenode.hpp"
#include "realmappedspnode.hpp"
#include "subpaving_exception.hpp"

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iterator>  
#include <cfloat> // DBL_EPSILON
#include <cassert> // assert

using namespace cxsc;
using namespace std;
using namespace subpavings;

std::vector < real > makeRanges1(cxsc::real rootvol);		
std::vector < real > makeRanges2(cxsc::real rootvol);		
std::vector < real > makeRanges3(cxsc::real rootvol);	
std::vector < real > makeRanges4(cxsc::real rootvol);
std::vector < real > makeRanges5(cxsc::real rootvol);
void outputNode(const std::string& s, 
			const RealMappedShapeNode& rmsn, const int prec);	




int main()
{
	
	ivector pavingBox1;
	ivector pavingBox2;
	ivector pavingBox3;
	int prec = 5;
	
	{
		int d = 1; 
		pavingBox1 = ivector(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox1[k] = pavingInterval;
	}
	{
		int d = 2; 
		pavingBox2 = ivector(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox2[k] = pavingInterval;
	}

	{
		int d = 3; 
		pavingBox3 = ivector(d);
		interval pavingInterval(-2,2);
		for(int k=1; k <= d; k++) pavingBox3[k] = pavingInterval;
	}

	std::vector < std::vector < double > > 
					backtransform1(1, std::vector < double >(1, 1.0));
	std::vector < std::vector < double > > 
					backtransform2(2, std::vector < double >(2, 0.0));
	std::vector < std::vector < double > > 
					backtransform3(3, std::vector < double >(3, 0.0));
	
	{ // first row of backtransform2
		double mydoubles[] = {.707, .707}; 
		backtransform2[0].assign (mydoubles,mydoubles+2);
	}
	{ // second row of backtransform2
		double mydoubles[] = {-0.707, .707}; 
		backtransform2[1].assign (mydoubles,mydoubles+2);
	}	
	{ // first row of backtransform3
		double mydoubles[] = {.36, .48, -0.8}; 
		backtransform3[0].assign (mydoubles,mydoubles+3);
	}
	{ // second row of backtransform3
		double mydoubles[] = {-0.8, .6, 0.0}; 
		backtransform3[1].assign (mydoubles,mydoubles+3);
	}
	{ // third row of backtransform3
		double mydoubles[] = {.48, .64, .6}; 
		backtransform3[2].assign (mydoubles,mydoubles+3);
	}
	
	/* backtransform2 is 45deg rotation clockwise, ie the reverse of a
	 * 45deg rotation counter clockwise */
	/* backtransform3 is a ~74deg rotation around (-1/3, 2/3, 2/3) */

	std::vector < double > scale1(1, 3.0);
	std::vector < double > scale2(2, 0.0);
	std::vector < double > scale3(1, 0.0);
	
	{ // scale 2
		double mydoubles[] = {1.5, 2.0}; 
		scale2.assign (mydoubles,mydoubles+2);
	}
	{ // scale 3
		double mydoubles[] = {1.5, 2.0, 2.5}; 
		scale3.assign (mydoubles,mydoubles+3);
	}
	
	std::vector < double > shift1(1, -1.0);
	std::vector < double > shift2(2, 0.0);
	std::vector < double > shift3(1, 0.0);
	
	{ // shift 2
		double mydoubles[] = {1.0, 1.5}; 
		shift2.assign (mydoubles,mydoubles+2);
	}
	{ // shift 3
		double mydoubles[] = {1.0, 1.5, 2.0}; 
		shift3.assign (mydoubles,mydoubles+3);
	}

	std::string split1 = "3,4,4,2,2,4,4,3";
	std::string split2 = "2,3,4,4,2,3,3";
	std::string split3 = "1,2,3,3";
	std::string split4 = "2,2,2,2";
	std::string split5 = "3,3,3,3,3,3,3,3";
	
	#if(0)
	{
		cout << "\n\nConstruct 1-d rmsp and turn into irregular shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, scale1, shift1);
		
		outputNode("ShapeScaleShift1.txt", shape, prec);
		
		assert(shape.getDimension() == rmsp.getDimension());

		
	}
	
	{
		cout << "\n\nConstruct 2-d rmsp and turn into irregular shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, scale2, shift2);
		
		outputNode("ShapeScaleShift2.txt", shape, prec);

		assert(shape.getDimension() == rmsp.getDimension());

	}
	
	{
		cout << "\n\nConstruct 3-d rmsp and turn into irregular shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, scale3, shift3);
		
		outputNode("ShapeScaleShift3.txt", shape, prec);

		assert(shape.getDimension() == rmsp.getDimension());
		
	}
	
	
	{
		cout << "\n\nConstruct 1-d rmsp and turn into irregular shape"
				<< " - backtransform, scale and shift" << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, backtransform1, scale1, shift1);
		
		outputNode("ShapeTransformScaleShift1.txt", shape, prec);

		assert(shape.getDimension() == rmsp.getDimension());

	}
	
	{
		cout << "\n\nConstruct 2-d rmsp and turn into irregular shape"
				<< " - backtransform, scale and shift" << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, backtransform2, scale2, shift2);
		
		outputNode("ShapeTransformScaleShift2.txt", shape, prec);

		assert(shape.getDimension() == rmsp.getDimension());
		
	}
	
	{
		cout << "\n\nConstruct 3-d rmsp and turn into irregular shape"
				<< " - backtransform, scale and shift" << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		cout << "\nrmsp is"  << endl;
		rmsp.leavesOutputTabs(cout);
		cout << "\nmake shape"  << endl;
		
		RealMappedShapeNode shape(rmsp, backtransform3, scale3, shift3);
		
		outputNode("ShapeTransformScaleShift3.txt", shape, prec);

		assert(shape.getDimension() == rmsp.getDimension());
		
	}
	
	{
		cout << "\n\nConstruct 1-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 only " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		std::vector < double > scale(1, 1.0);
		std::vector < double > shift(1, 0.0);
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	{
		cout << "\n\nConstruct 1-d pcf and turn into irregular polygon shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale1, shift1);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	
	{
		cout << "\n\nConstruct 2-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 only " << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		std::vector < double > scale(2, 1.0);
		std::vector < double > shift(2, 0.0);
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	{
		cout << "\n\nConstruct 2-d pcf and turn into irregular polygon shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale2, shift2);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	
	{
		cout << "\n\nConstruct 3-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 only " << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		std::vector < double > scale(3, 1.0);
		std::vector < double > shift(3, 0.0);
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	{
		cout << "\n\nConstruct 3-d pcf and turn into irregular polygon shape"
				<< " - scale and shift only " << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split3);
		std::vector< real > ranges = makeRanges3(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		cout << "\nmake shape"  << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale3, shift3);
		
		shape.outputToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	#endif
	
	#if(1)
	{
		cout << "\n\nConstruct 1-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split1);
		std::vector< real > ranges = makeRanges1(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		cout << "\npcf integral is " << (pcf.getTotalIntegral())  << endl;
		
		real cov(1.0);
		std::vector < double > scale(1, 1.0);
		std::vector < double > shift(1, 0.0);
		
		
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	
	{
		cout << "\n\nConstruct 1-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split1);
		std::vector< real > ranges = makeRanges1(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.80);
		std::vector < double > scale(1, 1.0);
		std::vector < double > shift(1, 0.0);
		
		
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	#endif
	#if(1)
	{
		cout << "\n\nConstruct 1-d pcf and turn into irregular polygon shape"
				<< " - scale and shift and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox1);
    
		rmsp.splitRootToShape(split1);
		std::vector< real > ranges = makeRanges1(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.80);
		
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale1, shift1, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	#endif
	
	{
		cout << "\n\nConstruct 2-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split4);
		std::vector< real > ranges = makeRanges4(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.80);
		std::vector < double > scale(2, 1.0);
		std::vector < double > shift(2, 0.0);
		
		
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	{
		cout << "\n\nConstruct 3-d pcf and turn into irregular polygon shape"
				<< " - scale=1 and shift=0 and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split5);
		std::vector< real > ranges = makeRanges5(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.50);
		std::vector < double > scale(3, 1.0);
		std::vector < double > shift(3, 0.0);
		
		
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale, shift, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	
	{
		cout << "\n\nConstruct 2-d pcf and turn into irregular polygon shape"
				<< " - scale and shift and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox2);
    
		rmsp.splitRootToShape(split4);
		std::vector< real > ranges = makeRanges4(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.80);
				
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale2, shift2, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	{
		cout << "\n\nConstruct 3-d pcf and turn into irregular polygon shape"
				<< " - scale and shift and coverage " << endl;
		RealMappedSPnode rmsp(pavingBox3);
    
		rmsp.splitRootToShape(split5);
		std::vector< real > ranges = makeRanges5(rmsp.nodeRealVolume());
		rmsp.allocateRanges(ranges);
		
		PiecewiseConstantFunction pcf(rmsp);
		
		cout << "\npcf is"  << endl;
		pcf.outputRootToStreamTabs(cout);
		
		real cov(0.50);
				
		cout << "\nmake shape with coverage " << cov << endl;
		
		PolygonPiecewiseConstantFunction shape(pcf, scale3, shift3, cov);
		
		shape.outputToStreamTabs(cout);
		
		cout << "\nArea vertices are" << endl;
		shape.outputAreaToStreamTabs(cout);
		
		assert(shape.getDimensions() == rmsp.getDimension());

		
	}
	
	
	cout << "\n\nend of test program" << endl;
    return 0;

} // end of test program

std::vector< real > makeRanges1(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((13.0/13.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((7.0/13.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((6.0/13.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((1.0/13.0)*(8.0/rootvol)));//XLLL
	result.push_back(cxsc::real((5.0/13.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real((3.0/13.0)*(16.0/rootvol)));//XLLRL
	result.push_back(cxsc::real((2.0/13.0)*(16.0/rootvol)));//XLLRR
	result.push_back(cxsc::real((1.0/13.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real((6.0/13.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/13.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((5.0/13.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((4.0/13.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((3.0/13.0)*(16.0/rootvol)));//XRRLL
	result.push_back(cxsc::real((1.0/13.0)*(16.0/rootvol)));//XRRLR
	result.push_back(cxsc::real((1.0/13.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges2(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((2.0/2.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real(0.0));//XLL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real(0.0));//XLRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((1.0/2.0)*(16.0/rootvol)));//XLRRL
	result.push_back(cxsc::real(0.0));//XLRRR
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real(0.0));//XRL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real(0.0));//XRRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges3(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((3.0/3.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((1.0/3.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/3.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((1.0/3.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real(0.0));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges4(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((10.0/10.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((6.0/10.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((5.0/10.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((1.0/10.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real((4.0/10.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((3.0/10.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((1.0/10.0)*(4.0/rootvol)));//XRR
	
	
	return result;
	
}

std::vector< real > makeRanges5(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((10.0/10.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((4.0/10.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((1.0/10.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((0.5/10.0)*(8.0/rootvol)));//XLLL
	result.push_back(cxsc::real((0.5/10.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real((3.0/10.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XLRL
	result.push_back(cxsc::real((2.0/10.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((6.0/10.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((2.0/10.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XRLL
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XRLR
	result.push_back(cxsc::real((4.0/10.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((3.0/10.0)*(8.0/rootvol)));//XRRR
	
	
	return result;
	
}

void outputNode(const std::string& s, 
			const RealMappedShapeNode& rmsn, const int prec)
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		rmsn.leavesOutputTabs(os, prec);
		os.close();
		cout << s << " output to file" << endl;
	}
	else cout << "Error opening file " << s << endl;
}


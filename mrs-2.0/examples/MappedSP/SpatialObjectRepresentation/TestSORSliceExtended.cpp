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
\brief Testing SpatialObjectRepresentationBV from function estimates
with slices for 3-D sphere images

 */

#include "TestSORTools.hpp"

#include "functionimage_estimator_booleanvalue.hpp"
#include "spatial_object_representation_bv.hpp"
#include "subpaving_exception.hpp"
#include "SphereFobj.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>


using namespace cxsc;
using namespace std;
using namespace subpavings;


SpatialObjectRepresentationBV makeImage(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				const std::string& fnamestart,
				int prec=5 );
	
void imageSlice(const SpatialObjectRepresentationBV& sor,
				const vector < vector < int > >& sliceDimsVec,
				const vector < vector < real > >& slicePtsVec,
				int d,
				real tolerance,
				const std::string& fnamestart,
				int prec);

void testSORSphereSlice()
{
	cout << "\n\nMake spatial object representations from a function image estimator:\n" << endl;
			
	cxsc::interval crit(1.0,cxsc::sqrt(2.0));
	interval pavingInterval(-3,3);
	string fnamestart("SORfromSphere");
	int prec = 5;
	
	int dimsArr[] = {3};
	vector<int> dims (dimsArr, dimsArr + sizeof(dimsArr) / sizeof(int) );
	
	real tolArr[] = {0.05};
	vector<real> tols (tolArr, tolArr + sizeof(tolArr) / sizeof(real) );
	
	vector < vector < int > > sliceDimsVec;
	{
		vector < int > sliceDims(1,3);  // dim 3
	
		sliceDimsVec.push_back(sliceDims);
		sliceDimsVec.push_back(sliceDims);
		sliceDimsVec.push_back(sliceDims);
		sliceDimsVec.push_back(sliceDims);
	}
	vector < vector < real > > slicePtsVec;
	{
		vector < real > slicePts(1,0.0);  
		slicePtsVec.push_back(slicePts);
	}
	{
		vector < real > slicePts(1,0.5);  
		slicePtsVec.push_back(slicePts);
	}
	{
		vector < real > slicePts(1,1.0);  
		slicePtsVec.push_back(slicePts);
	}
	{
		vector < real > slicePts(1,1.25);  
		slicePtsVec.push_back(slicePts);
	}
	
	
	for (int i = 0; i < dims.size(); ++i) {
		
		int d = dims[i];
		ivector pavingBox(d);
		
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		cout << "dim = " << d << endl;
		
		for (int j = 0; j < tols.size(); ++ j) {
			
			real tolerance = tols[j];
			
			cout << "tolerance = " << tolerance << endl;
			
			SphereFobj fobj;
					
			SpatialObjectRepresentationBV sor = makeImage(pavingBox,
						fobj,
						crit,
						tolerance,
						fnamestart,
						prec);

			
			// now do slice
			imageSlice(sor,
				sliceDimsVec,
				slicePtsVec,
				d,
				tolerance,
				fnamestart,
				prec) ;
			
			
			
			
		}
	}
	
}


SpatialObjectRepresentationBV makeImage(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				const std::string& fnamestart,
				int prec) 
{
			
			FunctionImageEstimatorBooleanValue fei(pavingBox, fobj, crit);
						
			fei.bruteForceEstimate(tolerance);
			
			SpatialObjectRepresentationBV sor = fei.makeSpatialObjectRepresentationBV();
			
			cout << "sor.getRootLeaves() = " << sor.getRootLeaves() << endl;
			
			ostringstream oss;
			oss << fnamestart << "_ForSlice_" << _double(tolerance) 
						<< "_" << VecLen(pavingBox) << "D.txt";
						
			string s1 = oss.str();
			sor.outputToTxtTabs(s1, prec, true);
			
			return sor;
			
}

void imageSlice(const SpatialObjectRepresentationBV& sor,
				const vector < vector < int > >& sliceDimsVec,
				const vector < vector < real > >& slicePtsVec,
				int d,
				real tolerance,
				const std::string& fnamestart,
				int prec) 
{

	for (int i = 0; i < sliceDimsVec.size(); ++i) {
		
		vector <int> sliceDims = sliceDimsVec[i];
		vector <real> slicePts = slicePtsVec[i];
		
		int dims = sor.getDimensions();
		
		string slicePtsStr;
		{
			ostringstream oss;
			for (int j = 0; j < dims; ++j) {
				for (int k = 0; k < sliceDims.size(); ++k) {
					if (j == sliceDims[k]-1)   oss << "_" << _double(slicePts[k]); 
					else   oss << "_x"; 
				}
			}
			slicePtsStr = oss.str();
		}
		
		SpatialObjectRepresentationBV slice = sor.makeSlice(sliceDims, slicePts);
		
		cout << "slice on " << slicePtsStr << " has getRootLeaves() = " 
			<< sor.getRootLeaves() << endl;
		{	
			ostringstream oss;
			oss << fnamestart << "_" << _double(tolerance) 
						<< "_" << d << "D_Slice"
						<< slicePtsStr << ".txt";
						
			slice.outputToTxtTabs(oss.str(), prec, true);
		}
		
	}

}

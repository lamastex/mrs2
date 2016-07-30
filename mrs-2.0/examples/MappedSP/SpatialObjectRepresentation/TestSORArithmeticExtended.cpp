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
with arithmetic for sphere images

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


SpatialObjectRepresentationBV* makeImagePtr(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				int cshift,
				const std::string& fnamestart,
				int prec=5 );
	
void imageArithmetic(const std::vector < SpatialObjectRepresentationBV* >& images,
				int d,
				real tolerance,
				const std::string& fnamestart,
				int prec);

void testSORSphereArithmetic()
{
	cout << "\n\nMake spatial object representations from a function image estimator:\n" << endl;
			
	cxsc::interval crit(1.0,cxsc::sqrt(2.0));
	interval pavingInterval(-3,3);
	string fnamestart("SORfromSphere");
	int prec = 5;
	
	// values to use for shifting centre
	real cValsArray[] = {0.7, 0.0, 0.7};
	vector<real> cVals (cValsArray, cValsArray + sizeof(cValsArray) / sizeof(real) );
	
	
	int dimsArr[] = {2, 3};
	vector<int> dims (dimsArr, dimsArr + sizeof(dimsArr) / sizeof(int) );
	
	real tolArr[] = {0.05, 0.03};
	vector<real> tols (tolArr, tolArr + sizeof(tolArr) / sizeof(real) );
	
	int cArr[] = {0, 1};
	vector<int> centreShifts (cArr, cArr + sizeof(cArr) / sizeof(int) );
	
	for (int i = 0; i < dims.size(); ++i) {
		
		int d = dims[i];
		ivector pavingBox(d);
		
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		cout << "dim = " << d << endl;
		
		for (int j = 0; j < tols.size(); ++ j) {
			
			real tolerance = tols[j];
			
			cout << "tolerance = " << tolerance << endl;
			
			std::vector < SpatialObjectRepresentationBV* > images;
					
			for (int c = 0; c < centreShifts.size(); ++ c) {
			
				int cshift = centreShifts[c];
				
				rvector centre(d);
				for(int k=1; k <= d; k++) { 
					
					if (cshift) {
						
						centre[k] = cVals[(k-1)%cVals.size()];
						
					}
					else centre[k] = 0.0;
				}
				
				cout << "centre = ";
				prettyPrint(cout, centre);
				cout << endl;
				
				
				SphereFobj fobj(centre);
					
				images.push_back( makeImagePtr(pavingBox,
							fobj,
							crit,
							tolerance,
							cshift,
							fnamestart,
							prec) );

			}
			
			// now do arithmetic
			imageArithmetic(images,
				d,
				tolerance,
				fnamestart,
				prec);
			
			//delete the images on the heap
			for (vector < SpatialObjectRepresentationBV* >::iterator it = images.begin();
					it < images.end();
					++it) {
					
				delete (*it);
				*it = NULL;
			}
			
		}
	}
	
}


SpatialObjectRepresentationBV* makeImagePtr(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				int cshift,
				const std::string& fnamestart,
				int prec) 
{
			
			FunctionImageEstimatorBooleanValue fei(pavingBox, fobj, crit);
						
			fei.bruteForceEstimate(tolerance);
			
			SpatialObjectRepresentationBV sor = fei.makeSpatialObjectRepresentationBV();
			
			cout << "sor.getRootLeaves() = " << sor.getRootLeaves() << endl;
			
			ostringstream oss;
			oss << fnamestart << "_BF" << _double(tolerance) 
						<< "_" << VecLen(pavingBox) << "D";
			if (cshift ) oss << "_shifted";
			oss << ".txt";
						
			string s1 = oss.str();
			sor.outputToTxtTabs(s1, prec, true);
			
			return new SpatialObjectRepresentationBV(sor);
			
}

void imageArithmetic(const std::vector < SpatialObjectRepresentationBV* >& images,
				int d,
				real tolerance,
				const std::string& fnamestart,
				int prec) 
{

	for (int i = 0; i < images.size()-1; ++i) {
		
		for (int j = i+1; j < images.size(); ++j) {
		
			{ // union
				SpatialObjectRepresentationBV sor = (*images[i]) + (*images[j]);
		
				cout << "union " << (i+1) << " and " << (j+1) 
					<< " has getRootLeaves() = " << sor.getRootLeaves() << endl;
				
				ostringstream oss;
				oss << fnamestart << "_Union" << _double(tolerance) 
							<< "_" << d << "D_"
							<< (i+1) << "_" << (j+1) << ".txt";
							
				sor.outputToTxtTabs(oss.str(), prec, true);
			}
			{ // XOR, or symmetric set difference
				SpatialObjectRepresentationBV sor = (*images[i]) - (*images[j]);
		
				ostringstream oss;
				oss << fnamestart << "_XOR" << _double(tolerance) 
							<< "_" << d << "D_"
							<< (i+1) << "_" << (j+1) << ".txt";
							
				sor.outputToTxtTabs(oss.str(), prec, true);
			}
			{ // union
				SpatialObjectRepresentationBV sor = (*images[i]) * (*images[j]);
		
				ostringstream oss;
				oss << fnamestart << "_Intersect" << _double(tolerance) 
							<< "_" << d << "D_"
							<< (i+1) << "_" << (j+1) << ".txt";
							
				sor.outputToTxtTabs(oss.str(), prec, true);
			}
			{ // set difference
				SpatialObjectRepresentationBV sor = (*images[i]) / (*images[j]);
		
				ostringstream oss;
				oss << fnamestart << "_SetDiff" << _double(tolerance) 
							<< "_" << d << "D_"
							<< (i+1) << "_" << (j+1) << ".txt";
							
				sor.outputToTxtTabs(oss.str(), prec, true);
			}
			
		}
		
	}
}

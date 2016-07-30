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


void makeImageBV(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				real cshift,
				const std::string& fnamestart,
				int prec=5 );
	

void testBVMakeFromFunctionImageEstimates()
{
	{
		int d = 2; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(-3,3);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		SphereFobj fobj;

		std::string split = "2,3,4,4,2,3,3";
		
		int prec = 5;

		try {
			
			cout << "\n\nMake a spatial object representation from a function image estimator:\n" << endl;
			int lab = 1;
			cxsc::interval crit(1.0,cxsc::sqrt(2.0));
			FunctionImageEstimatorBooleanValue fei(pavingBox, fobj, crit, lab);
			assert(fei.getLabel() == lab);
			
			try {
				fei.splitToShape(split);
				cout << "String summmary of FunctionImageEstimatorBooleanValue is: " << endl;
				cout << fei.stringSummary() << endl;
				
				ostringstream oss;
				fei.outputRootToStreamTabs(oss, prec);
				string feiRoot = oss.str();
				cout << "\nfei root is:\t" << feiRoot << endl;
				
								
			}
			catch (std::exception& ee) {
				std::string msg(ee.what());
				cout << "\nFailed to do splitToShape:\n" << msg << endl;
			}
			
			try {
				SpatialObjectRepresentationBV sor = fei.makeSpatialObjectRepresentationBV();
				
				cout << "String summmary of SpatialObjectRepresentationBV is: " << endl;
				cout << sor.stringSummary() << endl;
				
				assert(sor.getLabel() == lab);
				assert(sor.getRootLeaves() == fei.getRootLeaves());
				
				ostringstream oss;
				sor.outputRootToStreamTabs(oss, prec);
				string sorRoot = oss.str();
				
				cout << "sor root is:\t" << sorRoot << endl;
				
				string s1 = "SORfromSphere1_2D.txt";
				sor.outputToTxtTabs(s1, prec, true);
								
			}
			catch (std::exception& ee) {
				std::string msg(ee.what());
				cout << "\nFailed to makeSpatialObjectRepresentationBV:\n" << msg << endl;
			}
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to make a spatial object representation from a function estimator real:\n" << msg << endl;
			throw;
		}
	}	
	
	cout << "\n\nMake spatial object representations from a function image estimator:\n" << endl;
			
	cxsc::interval crit(1.0,cxsc::sqrt(2.0));
	interval pavingInterval(-3,3);
	string fnamestart("SORfromSphere");
	int prec = 5;
	
	int dimsArr[] = {2, 3};
	vector<int> dims (dimsArr, dimsArr + sizeof(dimsArr) / sizeof(int) );
	
	real tolArr[] = {0.05, 0.01};
	vector<real> tols (tolArr, tolArr + sizeof(tolArr) / sizeof(real) );
	
	real cArr[] = {0, 0.5};
	vector<real> centreShifts (cArr, cArr + sizeof(cArr) / sizeof(real) );
	
	for (int i = 0; i < dims.size(); ++i) {
		
		int d = dims[i];
		ivector pavingBox(d);
		
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
		
		cout << "dim = " << d << endl;
		
		for (int j = 0; j < tols.size(); ++ j) {
			
			real tolerance = tols[j];
			
			cout << "tolerance = " << tolerance << endl;
					
			for (int c = 0; c < centreShifts.size(); ++ c) {
			
				real cshift = centreShifts[c];
				
				cout << "centre shift = " << cshift << endl;
				
				if (cshift > 0.0) {
					rvector centre(d);
					for(int k=1; k <= d; k++) centre[k] = cshift;
					
					SphereFobj fobj(centre);
					
					for(int k=1; k <= d; k++) {
						Inf(pavingBox[k])+= centre[k];
						Sup(pavingBox[k])+= centre[k];
					}
					
					makeImageBV(pavingBox,
								fobj,
								crit,
								tolerance,
								cshift,
								fnamestart,
								prec);

				}
				else {
					SphereFobj fobj;
					makeImageBV(pavingBox,
								fobj,
								crit,
								tolerance,
								cshift,
								fnamestart,
								prec);
				}

			}
		}
	}
	
}


void makeImageBV(const cxsc::ivector& pavingBox,
				const MappedFobj& fobj,
				const cxsc::interval& crit,
				real tolerance,
				real cshift,
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
			if (cshift > 0.0) oss << "_shifted" << _double(cshift);
			oss << ".txt";
						
			string s1 = oss.str();
			sor.outputToTxtTabs(s1, prec, true);
			
}

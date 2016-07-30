/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief Testing Sphere function object
 */


#include "SphereFobj.hpp"
#include "toolz.hpp"

#include <stdexcept> // throwing exceptions
#include <iostream> 


using namespace cxsc;
using namespace subpavings;
using namespace std;

void doTest(int dd, const MappedFobj& fobj);						
						
int main(int argc, char* argv[])
{
	try {
		
		vector <int> dims;
		dims.push_back(2);
		//dims.push_back(3);
		
		for (size_t i = 0; i < dims.size(); ++i) {
			int dd = dims[i];
		
			rvector centre(dd);
			real centrePt(1.0);
			for(int i=1; i <= dd; i++) centre[i] = centrePt+1.0;
			
			cout << "\nTest with default centre = origin\n" << endl;
			{
				SphereFobj fobj;
				doTest(dd, fobj);
			}
			
			cout << "\nTest with default centre = ";
			prettyPrint(cout, centre);
			cout << "\n" << endl;
			{
				SphereFobj fobj(centre);
				//doTest(dd, fobj);
			}
					
		}
		
		return 0;
	}
	catch (std::exception& e) {
		cout << "Exception:\n" << e.what() << endl;
		throw;
	}
	catch (...) {
		cout << "Unknown exception" << endl;
		throw;
	}
	
}

void doTest(int dd, const MappedFobj& fobj)
{
	{
		ivector ivec(dd);
		rvector rvec(dd);
		real rr1(-0.8);
		real rr2(0.4);
		interval ii1(-1.5,1.1);
		interval ii2(-1.3,1.2);
		for(int k=1; k <= dd; k+=2) {
			ivec[k] = ii1;
			rvec[k] = rr1;
		}
		for(int k=2; k <= dd; k+=2) {
			ivec[k] = ii2;
			rvec[k] = rr2;
		}
		
		cout << "\n\nd = " << dd << endl;
		
		cout << "ivec = " << (toString(ivec)) << endl;
		interval i_image = fobj(ivec);
		
		cout << "interval image of ivec = " << (toString(i_image)) << endl;
		
		cout << "rvec = " << toString(rvec) << endl;
		real r_image = fobj(rvec);
		cout << "real image of rvec = " << r_image << endl;
	}
	{
		ivector ivec(dd);
		rvector rvec(dd);
		real rr1(6.0);
		real rr2(0.4);
		interval ii1(1.5,6.0);
		interval ii2(-1.3,1.2);
		for(int k=1; k <= dd; k+=2) {
			ivec[k] = ii1;
			rvec[k] = rr1;
		}
		for(int k=2; k <= dd; k+=2) {
			ivec[k] = ii2;
			rvec[k] = rr2;
		}
		
		cout << "\n\nd = " << dd << endl;
		
		cout << "ivec = " << (toString(ivec)) << endl;
		interval i_image = fobj(ivec);
		cout << "interval image of ivec = " << (toString(i_image)) << endl;
		
		cout << "rvec = " << toString(rvec) << endl;
		real r_image = fobj(rvec);
		cout << "real image of rvec = " << r_image << endl;
	}
	{
		ivector ivec(dd);
		rvector rvec(dd);
		real rr1(2.5);
		real rr2(5.0);
		interval ii1(-6,-5.5);
		interval ii2(-1.3,1.2);
		for(int k=1; k <= dd; k+=2) {
			ivec[k] = ii1;
			rvec[k] = rr1;
		}
		for(int k=2; k <= dd; k+=2) {
			ivec[k] = ii2;
			rvec[k] = rr2;
		}
		
		cout << "\n\nd = " << dd << endl;
		
		cout << "ivec = " << (toString(ivec)) << endl;
		interval i_image = fobj(ivec);
		cout << "interval image of ivec = " << (toString(i_image)) << endl;
		
		cout << "rvec = " << toString(rvec) << endl;
		real r_image = fobj(rvec);
		cout << "real image of rvec = " << r_image << endl;
	}
	
			

}

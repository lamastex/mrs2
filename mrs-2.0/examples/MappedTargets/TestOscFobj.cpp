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
\brief Testing Oscillating function object
 */


#include "oscFobj1.hpp"
#include "toolz.hpp"

#include <stdexcept> // throwing exceptions
#include <iostream> 


using namespace cxsc;
using namespace subpavings;
using namespace std;

						
						
int main(int argc, char* argv[])
{
	try {
		
		vector <int> dims;
		dims.push_back(1);
		
		vector <interval > domains;
		
		domains.push_back(interval(0.0, 1.0));
		domains.push_back(interval(0.0, 0.5));
		domains.push_back(interval(0.5, 1.0));
		domains.push_back(interval(0.0, 0.25));
		domains.push_back(interval(0.25, 0.5));
		domains.push_back(interval(0.5, 0.75));
		domains.push_back(interval(0.75, 1.0));
		domains.push_back(interval(0.75, 0.875));
		domains.push_back(interval(0.875, 1.0));
		
		OscFobj oscF;
		
		for (size_t i = 0; i < dims.size(); ++i) {
			int dd = dims[i];
		
			for (size_t j= 0; j < domains.size(); ++j) {
				
				ivector box(dd);
			
				for(int k=1; k <= dd; k++) { box[k] = domains[j]; }
			
				interval range = oscF(box);
				
				cout << "Box ";
				prettyPrint(cout,box);
				cout << "\nrange is " << range;
				cout << "\twidth is " << (Sup(range) - Inf(range) ) << endl;
				cout << "\tarea is " << ((Sup(range) - Inf(range))*(Sup(box) - Inf(box)) ) << endl;
				
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


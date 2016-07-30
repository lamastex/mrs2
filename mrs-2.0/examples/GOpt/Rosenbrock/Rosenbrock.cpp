/* 
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin
 * Copyright (C) 2009 Jennifer Harlow
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

/*! \file examples/GOpt/Rosenbrock/Rosenbrock.cpp
\brief Global optimisation example using Rosenbrock function
Using an example function object class with GOpt
*/

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"		
#include "imath.hpp"		
#include "rmath.hpp"		
#include "intvector.hpp"
#include "ivector.hpp"
#include "rvector.hpp"
#include "imatrix.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <functional>
#include<algorithm>
#include<numeric>

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;
using namespace cxsc;

#include "toolz.hpp"
#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FRosenbrock.hpp"
#include "GOpt.hpp"


int
main (int argc, char **argv)
{
	ios::sync_with_stdio ();	// call this function so iostream works with stdio
	cout << SetPrecision (20, 15);	// Number of mantissa digits in I/O
	
	// set default values
	int n_dimensions = 2; 
	
	//Parameters specific to the Rosenbrock target
	real Tinverse = 1.0;
	real Height = 100.0;
	real RosenDomainLimit = 10.0;
	
	bool UseLogPi = false;
	bool use_f_scale = true;
	
	cout << "Tinverse = " << Tinverse << "\n Height =  " << Height 
             << endl; 
       //getchar();
	
	UseLogPi = false; // log scale won't work naively
	
	// make the function object
	FRosenbrock FRosen (n_dimensions, Tinverse, Height, 
                      RosenDomainLimit, UseLogPi);
	
	// ***************** global optimisation ****************
	// set up a search box
	ivector search (1, n_dimensions);
	for (int i = 1; i <= n_dimensions; i++)
	{
		search[i] = interval (-RosenDomainLimit, RosenDomainLimit);
	}
	
	real tolerance;	// set a tolerance
	tolerance = 1e-16;
	
	//minimums
	GOptMin(&FRosen, search, tolerance); // call GoptMin with pointer to FRosen
	
	
	//maximums
	GOptMax(&FRosen, search, tolerance); // call GoptMin with pointer to FRosen
	
	return 0;			// end main statement
}

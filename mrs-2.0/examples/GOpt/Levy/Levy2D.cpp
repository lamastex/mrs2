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

/*! \file examples/GOpt/Levy/Levy2D.cpp
\brief GlobalOptimisation example using 2 dimensional Levy function
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
#include <interval.hpp>		// Include interval arithmetic package
#include <imath.hpp>		// Include interval standard functions
#include <rmath.hpp>		// Include real standard functions
#include <intvector.hpp>
#include <ivector.hpp>
#include <rvector.hpp>
#include <imatrix.hpp>
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
#include "FLevy2D.hpp"
#include "GOpt.hpp"


int main (int argc, char **argv)
{
	ios::sync_with_stdio ();	// call this function so iostream works with stdio
	cout << SetPrecision (20, 15);	// Number of mantissa digits in I/O
	
	
	int n_dimensions = 2; 
	int n_boxes = 1000;
	int n_samples = 100;
	double Alb = 1.0;// partition until lower bound on Acceptance Prob. is > Alb
	unsigned theSeed = 0;
	
	
	bool UseLogPi = false;
	bool use_f_scale = true;
	
	if (argc >= 2)
	{
		sscanf (argv[1], "%i", &n_boxes);
		if (argc >= 3)
		{
			sscanf (argv[2], "%i", &n_samples);
			if (argc >= 4){
				sscanf (argv[3], "%ui", &theSeed);
			}
			if (argc >= 5)
				cout << "# Usage: MRS <n_boxes> <n_samples> <seed>; "
             << "extra arguments ignored.\n";
		}
		
		else cout << "# Usage: MRS <n_boxes> <n_samples> <seed>; "
              << "extra arguments ignored.\n";
	}
	
	cout << "# n_dimensions: " << n_dimensions << "  n_boxes: " 
       << n_boxes << "  n_samples: " 
		   << n_samples << "  rng_seed = " << theSeed  << endl; //getchar();
	
	//Parameters specific to the Levy target
	real Temperature = 10.0;
	real Center1 = 1.42513;
	real Center2 = 0.80032;
	real GlobalMax = 176.14;
	real DomainLimit = 10.0;	//0.999999999999999;
	UseLogPi = false; // log scale won't work naively
	
	// make the function object
	FLevy2D F_Levy_Temp_2D(Temperature, GlobalMax, 
                         Center1, Center2, DomainLimit, UseLogPi);
	
	
	// ***************** global optimisation ****************
	
	// set up a search box
	ivector search (1, n_dimensions);
	for (int i = 1; i <= n_dimensions; i++)
	{
		search[i] = interval (-DomainLimit, DomainLimit);
	}
	
	real tolerance;	// set a tolerance
	tolerance = 1e-8;
	
	// Minimums
  // call GoptMin with pointer to F_Levy_Temp_2D
	GOptMin(&F_Levy_Temp_2D, search, tolerance); 
	
	// Maximums
  // call GoptMax with pointer to F_Levy_Temp_2D
	GOptMax(&F_Levy_Temp_2D, search, tolerance); 
	
	return 0;			// end main statement
}

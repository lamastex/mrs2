/* 
 * Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
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

/*! \file examples/MooreRejSam/Levy/Levy2D.cpp
\brief MRS example using 2 dimensional Levy function
Using an example function object class with MRSampler class
*/

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"		// Include interval arithmetic package
#include "imath.hpp"		// Include interval standard functions
#include "rmath.hpp"		// Include real standard functions
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
#include "FLevy2D.hpp"
#include "MRSampler.hpp"

void
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, 
                 double Alb, unsigned seed, bool use_f_scale) 
{
	ofstream out ("MRS_Levy2D.samples");//file to store the i.i.d samples
	clock_t T1 = clock (), T2, T3;
	// Construct theSampler with the chosen target shape object FTG
	MRSampler theSampler (f, n_boxes, Alb, seed, (use_f_scale == 1));
	//MRSampler theSampler(f, n_boxes, 0, Alb, seed);  
	/* // To print out the partition of the domain
	cout << "Domain Partition: \n" ;
	theSampler.Print_Domain_Partition();
	theSampler.PrintBoxes(0);
	getchar(); // press enter to continue...
	*/
	
	// To realize a file output of the RangeDomainSet	
	ofstream os("Levy2DOutput.txt");         // Filename
	os << theSampler << endl;                   
	cout << "The output has been written to Levy2DOutput" << endl << endl;
	
	T2 = clock ();
	double Ptime = (double) (T2 - T1) / CLOCKS_PER_SEC;
	
	RSSample rs_sample;
	cout << "before Rej..SampleMany \n";
	cout << "n_samples: " << n_samples << endl;
	theSampler.RejectionSampleMany (n_samples, rs_sample);
	cout << "after Rej..SampleMany \n";
	double IntegralEstimate = _double (rs_sample.IntegralEstimate ());
	cout << "rs_sample IU, N, Nrs: " << rs_sample.EnvelopeIntegral << " " 
       << rs_sample.Samples.size() << " " << rs_sample.Samples.size() << endl;
	cout << "RSSampleMany, integral est: " << IntegralEstimate << endl;
	cout << "RSSampleMany mean: \n"; rs_sample.Mean ();
	rs_sample.Print(out);
	
	cout << "n interval function calls: " << f.get_interval_calls () << endl;
	cout << "n real function calls: " << f.get_real_calls () << endl;
	
	//----------------------------------------------------------------------------
	T3 = clock ();
	double Stime = (double) (T3 - T2) / CLOCKS_PER_SEC;
	cout << "# CPU Time (seconds). Partitioning: " << Ptime << "  Sampling: " 
       << Stime << "  Total: " << (Ptime + Stime) << endl;
	cout << "# CPU time (secods) per estimate: " 
       << (Ptime + Stime) / (double) (n_samples) << endl;
	
	
}

int
main (int argc, char **argv)
{
	ios::sync_with_stdio ();	// call this function so iostream works with stdio
	cout << SetPrecision (20, 15);	// Number of mantissa digits in I/O
	
	
	int n_dimensions = 2; 
	int n_boxes = 1000;
	int n_samples = 10000;
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
	
	cout << "# n_dimensions: " << n_dimensions << "  n_boxes: " << n_boxes 
       << "  n_samples: " << n_samples << "  rng_seed = " << theSeed  
       << endl; //getchar();
	
	//Parameters specific to the Levy target
	real Temperature = 40.0;
	real Center1 = 1.42513; 
	real Center2 = 0.80032; 
	real GlobalMax = 176.14;
	real DomainLimit = 10.0;	//0.999999999999999;
	UseLogPi = false; // log scale won't work naively
	
	// make the function object
	FLevy2D F_Levy_Temp_2D(Temperature, GlobalMax, 
                         Center1, Center2, DomainLimit, UseLogPi);
	
	// produce the samples
	ProduceMRSamples(F_Levy_Temp_2D, n_boxes, n_samples, 
                   Alb, theSeed, use_f_scale);
	

	return 0;			// end main statement
}

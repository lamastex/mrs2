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

/*! \file examples/MooreRejSam/Rosenbrock/Rosenbrock.cpp
\brief MRS example using Rosenbrock function
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
#include "FRosenbrock.hpp"
#include "MRSampler.hpp"

void
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, 
                 double Alb, unsigned seed, bool use_f_scale) 
{
	ofstream out ("MRS_Rosenbrock.samples");//file to store the i.i.d samples
	clock_t T1 = clock (), T2, T3;
	// Construct theSampler with the chosen target shape object FTG
	MRSampler theSampler (f, n_boxes, Alb, seed, (use_f_scale == 1));
	
        ///////////////////////////////////////////////////////////////////////
        //--comment this printing block if not needed
        // To print out the partition of the domain
        //cout << "Domain Partition: \n" ;
        //ofstream Partout ("MRS_RosenbrockDomain.txt"); //Filename
        //theSampler.Print_Domain_Partition(Partout);
        //cout << "The output has been written to MRS_IsIt1or2CoinsDomain.txt" 
        //     << endl << endl;
 
        // To realize a file output of the RangeDomainSet
        ofstream os("MRS_RosenbrockRangeDomainSet.txt");         // Filename
        os << theSampler << endl;                   
        cout << "output has been written to MRS_IsIt1or2CoinsRangeDomainSet.txt" 
             << endl << endl;
        //--end of printing block
        ///////////////////////////////////////////////////////////////////////
	// To print out the partition of the domain
	
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
	
	// set default values
	int n_dimensions = 2; 
	int n_boxes = 100;
	int n_samples = 50;
	double Alb = 1.0;// partition until lower bound on Acceptance Prob. is > Alb
	unsigned theSeed = 0;
	
	//Parameters specific to the Rosenbrock target
	real Tinverse = 1.0;
	real Height = 100.0;
	real RosenDomainLimit = 10.0;
	
	if (argc >= 2)
	{
		sscanf (argv[1], "%i", &n_dimensions);
		if (argc >= 3)
		{
			sscanf (argv[2], "%i", &n_boxes);
			
			if (argc >= 4)
				sscanf (argv[3], "%i", &n_samples);
			if (argc >= 5) 
				sscanf (argv[4], "%ui", &theSeed);
			if (argc >= 6)
				cout << "# Usage: MRS <n_dimensions> <n_boxes> <n_samples> <seed>; "
             << "extra arguments ignored.\n";
		}
		
		else cout << "# Usage: MRS <n_dimensions> <n_boxes> <n_samples> <seed>; "
              << "extra arguments ignored.\n";
		
	}
	
	cout << "# n_dimensions: " << n_dimensions << "  n_boxes: " << n_boxes 
       << "  n_samples: " << n_samples << "  rng_seed = " << theSeed  
       << endl; //getchar();
	
	bool UseLogPi = false; // log scale won't work naively
	bool use_f_scale = false;
	
	// make the function object
	FRosenbrock FRosen (n_dimensions, 
                      Tinverse, Height, RosenDomainLimit, UseLogPi);
	
	// produce the samples
	ProduceMRSamples(FRosen, n_boxes, n_samples, 
                   Alb, theSeed, use_f_scale);
	
	return 0;			// end main statement
}

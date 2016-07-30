/* 
 * Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
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
/*! \file      examples/MooreRejSam/PhyloPOT/PhyloPOT.cpp
\brief Example to use FPhyloPOT and MRSampler to produce rejection samples 
for trans-dimensional phylogenetic problems.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <string>

#include <math.h>
#include <time.h>

//CSXC includes
#include "interval.hpp"		// Include interval arithmetic package
#include "imath.hpp"		// Include interval standard functions
#include "rmath.hpp"		// Include real standard functions
#include "intvector.hpp"
#include "ivector.hpp"
#include "rvector.hpp"
#include "imatrix.hpp"

#include <stdio.h>
#include <stdlib.h>

#include "toolz.hpp"
#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FPhyloPOT.hpp"
#include "MRSampler.hpp"


using namespace std;
using namespace cxsc;


void
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, 
                 double Alb, unsigned seed, bool use_f_scale) 
{
	clock_t T1 = clock (), T2, T3;
  	// Construct theSampler with the chosen target shape object f
	MRSampler theSampler (f, n_boxes, Alb, seed, (use_f_scale == 1));
	T2 = clock ();
	double Ptime = (double) (T2 - T1) / CLOCKS_PER_SEC;
	
	//ofstream dout("PhyloPOTdomain.txt");
	//theSampler.Output_Domain_Partition(dout);
	
	// To realize a file output of the RangeDomainSet
	//ofstream os("PhyloPOTOutput.txt");         // Filename
	//os << theSampler << endl;                   // 
	//cout << "The output has been written to PhyloPOTOutput.txt" << endl << endl;
	
	ofstream out ("ChorB100_MRS_PhyloPOT.samples");//file to store the i.i.d samples
	
	
	RSSample rs_sample;
	
	T3 = clock ();
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

////	cout << "n GradType function calls: " << f.get_gradtype_calls () << endl;
	
	cout << "n real function calls: " << f.get_real_calls () << endl;
	
		
	//----------------------------------------------------------------------------
	double Stime = (double) (T3 - T2) / CLOCKS_PER_SEC;
	cout << "# CPU Time (seconds). Partitioning: " << Ptime << "  Sampling: " 
       << Stime << "  Total: " << (Ptime + Stime) << endl;
	cout << "# CPU time (seconds) per estimate: " 
       << (Ptime + Stime) / (double) (n_samples) << endl;
	
}

int main (int argc, char **argv)
{
	ios::sync_with_stdio ();	// call this function so iostream works with stdio
	cout << SetPrecision (20, 15);	// Number of mantissa digits in I/O
	
	// default value
	double Alb = 1.0;// partition until lower bound on Acceptance Prob. is > Alb
	
	// default value
	int n_boxes = 10000;
	
	// default value
        int CenFrm=1;

	// default value
	int n_samples = 10000;
	
	// default value
	unsigned theSeed = 1234;
	
	// default value
	bool UseLogPi = true;
	
	// default value
	bool use_f_scale = true;
	
	// default value
	int prior_type = 0;
	
	// parameters for FPhyloPOT
	// default value
	interval DomainCFN3(1e-10,1.0);
	
	// default value
	int tree_space = 3;	// number of taxa
        int character_space = 4; //number of nucleotide characters
	
	// check the input 
	if (argc >= 2)
	{
		sscanf (argv[1], "%i", &tree_space);
cout << "aaa" << endl;
		if (argc >= 3) 
			sscanf (argv[2], "%i", &character_space);
cout << "aaa" << endl;
		if (argc >= 4)
cout << "aaa" << endl;
			sscanf (argv[3], "%i", &CenFrm);
cout << "bbb" << endl;
		if (argc >= 5) 
			sscanf (argv[4], "%i", &n_boxes);
cout << "aaa" << endl;
		if (argc >= 6)
			sscanf (argv[5], "%i", &n_samples);
cout << "aaa" << endl;
		if (argc >= 7) 
			sscanf (argv[6], "%ui", &theSeed);
		if (argc >=8)
			cout << "# Usage: MRS <tree_space> <character_space> <centered_form> "
				<< " <n_boxes> <n_samples> <seed> "
           			<< "extra arguments ignored.\n";
	}
	
	bool CenteredForm;
        CenteredForm = CenFrm;

	cout << "# tree_space: " << tree_space << " # characters: " << character_space 
		<< " Centered Form " << CenteredForm << "  n_boxes: " << n_boxes 
       		<< "  n_samples: " << n_samples << "  rng_seed = " << theSeed  << endl; 
	
	
	FPhyloPOT FPhylo(tree_space, character_space, DomainCFN3, CenteredForm, UseLogPi, prior_type);
	
	ProduceMRSamples(FPhylo, n_boxes, n_samples, Alb, theSeed, use_f_scale);
	
	
	return 0;			
}

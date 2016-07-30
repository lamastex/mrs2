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
#include "FBinomialPartition.hpp"
#include "MRSampler.hpp"

void
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, 
                 double Alb, unsigned seed, bool use_f_scale) 
{
  ofstream out ("MRS_BinomialPartition.samples");//file to store the i.i.d samples
  clock_t T1 = clock (), T2, T3;
  // Construct theSampler with the chosen target shape object FTG
  MRSampler theSampler (f, n_boxes, Alb, seed, (use_f_scale == 1));
  //////////////////////////////////////////////////////////////////////////////////
  //--comment this printing block if not needed as this will SLOW DOWN due to I/O
/*     // To print out the partition of the domain
     //cout << "Domain Partition: \n" ;
     //ofstream Partout ("MRS_BinomialPartitionDomain.txt"); //Filename
     //theSampler.Print_Domain_Partition(Partout);
     //cout << "The output has been written to MRS_BinomialPartitionDomain.txt" 
     //     << endl << endl;
  
     // To realize a file output of the RangeDomainSet
     ofstream os("MRS_BinomialPartitionRangeDomainSet.txt");         // Filename
     os << theSampler << endl;                   
     cout << "The output has been written to MRS_BinomialPartitionRangeDomainSet.txt" 
          << endl << endl;
*/  //--end of printing block
  //////////////////////////////////////////////////////////////////////////////////
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
  ///---------------------------------------------------------------------------
  ///---To see the accepted samples uncomment this - I/O intensive outputcall 
  ///---best to comment it if you are only interested in posterior model probs 
  ///---and MAP estimates
  //rs_sample.Print(out);
  ///---------------------------------------------------------------------------
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

  int n_boxes = 1000000;
  int n_samples = 1;
  double Alb = 0.1;// partition until lower bound on Acceptance Prob. is > Alb
  unsigned theSeed = 0;
  bool UseLogPi = true;//now using the real power function by log in likelihood
  //bool UseLogPi = false;//now using the integer power function in likelihood  
  bool use_f_scale = true;
  int tmpULP= 1;
  int tmpint= 0;
  vector<int> Ns(0,2);
  vector<int> Ys(0,2);

  if (argc >= 2)
    {
      sscanf (argv[1], "%i", &n_boxes);
      if (argc >= 3)
	{
	  sscanf (argv[2], "%i", &n_samples);
	  if (argc >= 4){
            //cout << "theSeed: " << argv[3] << endl;
	    sscanf (argv[3], "%ui", &theSeed);	    
	  }
	  if (argc >= 5){
	    //cout << "N[0]: " << argv[4] << endl;
	    sscanf (argv[4], "%i", &tmpint);
            Ns[0]=tmpint; 
	  }
	  if (argc >= 6) {
	    //cout << "Y[0]: " << argv[5] << endl;
	    sscanf (argv[5], "%i", &tmpint);
            Ys[0]=tmpint; 
	  }
	  if (argc >= 7){
	    //cout << "N[1]: " << argv[6] << endl;
	    sscanf (argv[6], "%i", &tmpint);
            Ns[1]=tmpint; 
	  }
	  if (argc >= 8){
	    //cout << "Y[1]: " << argv[7] << endl;
	    sscanf (argv[7], "%i", &tmpint);
            Ys[1]=tmpint; 
	  }
	  if(argc >= 9){
	    sscanf (argv[8], "%i", &tmpULP);
	    UseLogPi = (bool)tmpULP;
	    //cout << "UseLogPi: " << argv[8] << endl;
	  }  
	  if (argc >= 10){
	   cout << "# Usage: " << argv[0] << "<n_boxes> <n_samples> "
           << "<seed> <N[0]=number of first coin tosses> <Y[0]=number of heads in first coin tosses>"
           << "<N[1]=number of second coin tosses> <Y[1]=number of heads in second coin tosses> <UseLogPi>" 
           << "; extra arguments ignored.\n";
	  }

	}
    }

  else { 
	 cout << "# Usage: " << argv[0] << "<n_boxes> <n_samples> "
         << "<seed> <N[0]=number of first coin tosses> <Y[0]=number of heads in first coin tosses>"
         << "<N[1]=number of second coin tosses> <Y[1]=number of heads in second coin tosses> <UseLogPi>" 
         << "; extra arguments ignored.\n";

  }
   
  cout << "  n_boxes: " << n_boxes << "  n_samples: " << n_samples << "  rng_seed = " 
       << theSeed  << endl; 
/*
       << "  N[0]=number of first coin tosses          : " << Ns[0] << endl
       << "  Y[0]=number of heads in first coin tosses : " << Ys[0] << endl
       << "  N[1]=number of second coin tosses         : " << Ns[1] << endl
       << "  Y[1]=number of heads in second coin tosses: " << Ys[1] << endl; //getchar();
*/
  // FBinomialPartition(int N[0], int Y[0], int N[1], int Y[1], bool LogPi);
  //FBinomialPartition FTG (Ns, Ys, UseLogPi);
  FBinomialPartition FTG (UseLogPi);
  ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);

  return 0;			// for the main statement
}

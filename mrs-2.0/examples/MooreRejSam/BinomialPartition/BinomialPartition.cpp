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
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, int M,
                 double Alb, unsigned seed, bool use_f_scale) 
{
  //ofstream out ("MRS_BinomialPartition.samples");//file to store the i.i.d samples
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

  /*std::vector<rvector> means = rs_sample.Mean();
  std::vector<rvector> means;
  std::vector<int> labels;
  std::vector<real> props;*/

  MeansLabelsProportions mlp = rs_sample.MeanLabelProportion();
  // more convenient output
      cout << "label proportions posterior_means" << endl; 
      std::vector<rvector>::const_iterator it;
      for(size_t i=0; i<mlp.means.size(); i++)
      {
        cout << mlp.labels[i] << " " << mlp.proportions[i] << " " << endl << mlp.means[i] << endl;
      }
    cout << "______________________ number of treatments, models = "<< M << " , " << mlp.means.size() << endl;
    cout << "model_labels number_of_parameters posterior_probability posterior_means" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
    if(M==1){//just 1 binomial trial y_1
      // model 0 = {{1}}
      cout << mlp.labels[0] << " "  << VecLen(mlp.means[0]) << " "  << mlp.proportions[0] << "\t"; 
      if (mlp.proportions[0]==0.0) cout << "NaN" << endl;
      else cout << mlp.means[0][1] << endl;
    }
    if(M==2){//just 2 binomial trials y_1,y_2
      //model 0 = {{1,2}}
      cout << mlp.labels[0] << "\t"  << VecLen(mlp.means[0]) << "\t"  << mlp.proportions[0] << "\t"; 
      if (mlp.proportions[0]==0.0) cout << "NaN\tNaN" << endl;
      else cout << mlp.means[0][1] << "\t" << mlp.means[0][1] << endl;
      //model 1 = {{1},{2}}
      cout << mlp.labels[1] << "\t"  << VecLen(mlp.means[1]) << "\t"  << mlp.proportions[1] << "\t"; 
      if (mlp.proportions[0]==0.0) cout << "NaN\tNaN" << endl;
      else cout  << mlp.means[1][1] << "\t" << mlp.means[1][2] << endl;
    }
    if(M==3){//just 3 binomial trials y_1,y_2,y_3
      //model 0 = {{1,2,3}}
      cout << mlp.labels[0] << "\t" << VecLen(mlp.means[0]) << "\t" << mlp.proportions[0] << "\t"; 
      if (mlp.proportions[0]==0.0) cout << "NaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[0][1] << "\t"
           << mlp.means[0][1] << "\t"
           << mlp.means[0][1] << endl;
      //model 1 = {{1},{2,3}}
      cout << mlp.labels[1] << "\t" << VecLen(mlp.means[1]) << "\t" << mlp.proportions[1] << "\t"; 
      if (mlp.proportions[1]==0.0) cout << "NaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[1][1] << "\t"
           << mlp.means[1][2] << "\t"
           << mlp.means[1][2] << endl;
      //model 2 = {{2},{1,3}}
      cout << mlp.labels[2] << "\t"  << VecLen(mlp.means[2]) << "\t"  << mlp.proportions[2] << "\t"; 
      if (mlp.proportions[2]==0.0) cout << "NaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[2][2] << "\t"
           << mlp.means[2][1] << "\t"
           << mlp.means[2][2] << endl;
      //model 3 = {{3},{1,2}}
      cout << mlp.labels[3] << "\t"  << VecLen(mlp.means[3]) << "\t"  << mlp.proportions[3] << "\t"; 
      if (mlp.proportions[3]==0.0) cout << "NaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[3][2] << "\t"
           << mlp.means[3][2] << "\t"
           << mlp.means[3][1] << endl;
      //model 4 = {{1},{2},{3}}
      cout << mlp.labels[4] << "\t"  << VecLen(mlp.means[4]) << "\t"  << mlp.proportions[4] << "\t"; 
      if (mlp.proportions[4]==0.0) cout << "NaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[4][1] << "\t"
           << mlp.means[4][2] << "\t"
           << mlp.means[4][3] << endl;
    }
    if(M==4){//just 4 binomial trials y_1,y_2,y_3,y_4
      //model label 0 = [{1, 2, 3, 4}] 
      cout << mlp.labels[0] << "\t"  << VecLen(mlp.means[0]) << "\t"  << mlp.proportions[0] << "\t"; 
      if (mlp.proportions[0]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[0][1] << "\t" 
           << mlp.means[0][1] << "\t" 
           << mlp.means[0][1] << "\t" 
           << mlp.means[0][1] << endl;
      // model label 1 = [{1}, {2, 3, 4}], 
      cout << mlp.labels[1] << "\t"  << VecLen(mlp.means[1]) << "\t"  << mlp.proportions[1] << "\t";
      if (mlp.proportions[1]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[1][1] << "\t" 
           << mlp.means[1][2] << "\t" 
           << mlp.means[1][2] << "\t" 
           << mlp.means[1][2] << endl;
      // model label 2 = [{2}, {1, 3, 4}], 
      cout << mlp.labels[2] << "\t"  << VecLen(mlp.means[2]) << "\t"  << mlp.proportions[2] << "\t";
      if (mlp.proportions[2]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[2][2] << "\t" 
           << mlp.means[2][1] << "\t" 
           << mlp.means[2][2] << "\t" 
           << mlp.means[2][2] << endl;
      // model label 3 = [{3}, {1, 2, 4}], 
      cout << mlp.labels[3] << "\t"  << VecLen(mlp.means[3]) << "\t"  << mlp.proportions[3] << "\t";
      if (mlp.proportions[3]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[3][2] << "\t" 
           << mlp.means[3][2] << "\t" 
           << mlp.means[3][1] << "\t" 
           << mlp.means[3][2] << endl;
      // model label 4 = [{4}, {1, 2, 3}], 
      cout << mlp.labels[4] << "\t"  << VecLen(mlp.means[4]) << "\t"  << mlp.proportions[4] << "\t";
      if (mlp.proportions[4]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[4][2] << "\t" 
           << mlp.means[4][2] << "\t" 
           << mlp.means[4][2] << "\t" 
           << mlp.means[4][1] << endl;
      // model label 5 = [{1, 2}, {3, 4}],
      cout << mlp.labels[5] << "\t"  << VecLen(mlp.means[5]) << "\t"  << mlp.proportions[5] << "\t";
      if (mlp.proportions[5]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[5][1] << "\t" 
           << mlp.means[5][1] << "\t" 
           << mlp.means[5][2] << "\t" 
           << mlp.means[5][2] << endl;
      // model label 6 = [{1, 3}, {2, 4}], 
      cout << mlp.labels[6] << "\t"  << VecLen(mlp.means[6]) << "\t"  << mlp.proportions[6] << "\t";
      if (mlp.proportions[6]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[6][1] << "\t" 
           << mlp.means[6][2] << "\t" 
           << mlp.means[6][1] << "\t" 
           << mlp.means[6][2] << endl;
      // model label 7 = [{1, 4}, {2, 3}], 
      cout << mlp.labels[7] << "\t"  << VecLen(mlp.means[7]) << "\t"  << mlp.proportions[7] << "\t";
      if (mlp.proportions[7]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[7][1] << "\t" 
           << mlp.means[7][2] << "\t" 
           << mlp.means[7][2] << "\t" 
           << mlp.means[7][1] << endl;
      // model label 8 = [{1}, {2}, {3, 4}], 
      cout << mlp.labels[8] << "\t"  << VecLen(mlp.means[8]) << "\t"  << mlp.proportions[8] << "\t";
      if (mlp.proportions[8]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[8][1] << "\t" 
           << mlp.means[8][2] << "\t" 
           << mlp.means[8][3] << "\t" 
           << mlp.means[8][3] << endl;
      // model label 9 = [{1}, {3}, {2, 4}],
      cout << mlp.labels[9] << "\t"  << VecLen(mlp.means[9]) << "\t"  << mlp.proportions[9] << "\t";
      if (mlp.proportions[9]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[9][1] << "\t" 
           << mlp.means[9][3] << "\t" 
           << mlp.means[9][2] << "\t" 
           << mlp.means[9][3] << endl;
      // model label 10 = [{1}, {4}, {2, 3}], 
      cout << mlp.labels[10] << "\t"  << VecLen(mlp.means[10]) << "\t"  << mlp.proportions[10] << "\t";
      if (mlp.proportions[10]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[10][1] << "\t" 
           << mlp.means[10][3] << "\t" 
           << mlp.means[10][3] << "\t" 
           << mlp.means[10][2] << endl;
      // model label 11 = [{2}, {3}, {1, 4}], 
      cout << mlp.labels[11] << "\t"  << VecLen(mlp.means[11]) << "\t"  << mlp.proportions[11] << "\t";
      if (mlp.proportions[11]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[11][3] << "\t" 
           << mlp.means[11][1] << "\t" 
           << mlp.means[11][2] << "\t" 
           << mlp.means[11][3] << endl;
      // model label 12 = [{2}, {4}, {1, 3}], 
      cout << mlp.labels[12] << "\t"  << VecLen(mlp.means[12]) << "\t"  << mlp.proportions[12] << "\t";
      if (mlp.proportions[12]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[12][3] << "\t" 
           << mlp.means[12][1] << "\t" 
           << mlp.means[12][3] << "\t" 
           << mlp.means[12][2] << endl;
      // model label 13 = [{3}, {4}, {1, 2}], 
      cout << mlp.labels[13] << "\t"  << VecLen(mlp.means[13]) << "\t"  << mlp.proportions[13] << "\t";
      if (mlp.proportions[13]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[13][3] << "\t" 
           << mlp.means[13][3] << "\t" 
           << mlp.means[13][1] << "\t" 
           << mlp.means[13][2] << endl;
      // model label 14 = [{1}, {2}, {3}, {4}], 
      cout << mlp.labels[14] << "\t"  << VecLen(mlp.means[14]) << "\t"  << mlp.proportions[14] << "\t";
      if (mlp.proportions[14]==0.0) cout << "NaN\tNaN\tNaN\tNaN" << endl;
      else cout 
           << mlp.means[14][1] << "\t" 
           << mlp.means[14][2] << "\t" 
           << mlp.means[14][3] << "\t" 
           << mlp.means[14][4] << endl;
     }
    cout << "-----------------------------------------------------------------------" << endl;
    cout << "______________________ end of model posteriors" << endl;

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
  int M = 4;
  vector<int> Ns; //(0,2);
  vector<int> Ys; //(0,2);

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
	    sscanf (argv[4], "%i", &tmpULP);
	    UseLogPi = (bool)tmpULP;
	    //cout << "UseLogPi: " << argv[8] << endl;
	  }
	  if (argc >= 6) {
	    sscanf (argv[5], "%i", &M);
            assert(M >=1 && M <= 4);
            Ns.push_back(M); // M=number of treatments 
            Ys.push_back(M); // M=number of treatments 
            assert(argc == 6+2*M);// make sure we have M ys and M ns
	  }
	  if (argc >= 6+2*M) {
	    for (int i=6; i < argc; ++i){
	      sscanf (argv[i], "%i", &tmpint); Ns.push_back(tmpint); 
	      sscanf (argv[1+i], "%i", &tmpint); Ys.push_back(tmpint); 
              ++i;
	    }
          }
	  }  
    }

  else { 
	 cout << "# Usage: " << argv[0] << "<n_boxes> <n_samples> "
         << "<seed> <UseLogPi> <NumberOfTreatments=M>"
         << " <N[1]=number of trials in first treatment> <Y[1]=number of successes in first treatment>" 
         << " <N[2]=number of trials in second treatment> <Y[2]=number of successes in second treatment>" 
         << " ... "
         << " <N[M]=number of trials in last treatment> <Y[M]=number of successes in last treatment>" 
         << "; extra arguments ignored.\n";
         cout << "Example usage for 4 treatments of pine seedling data:" << endl
              << "$./BinomialPartition 250000 100000 0 1 4 100 59 100 89 100 88 100 95 " << endl;

        exit(1);
  }
   
  cout << "  n_boxes: " << n_boxes << "  n_samples: " << n_samples << "  rng_seed: " << theSeed 
       << " UseLogPi: " << UseLogPi << " Number_of_Treatments=M: " << M << endl; 
  cout << "M N[1] N[2] ... N[M]" << endl;
  for (int i = 0; i < Ns.size(); i++) cout << Ns[i] << " "; cout << endl;
  cout << "M N[1] N[2] ... N[M]" << endl;
  for (int i = 0; i < Ys.size(); i++) cout << Ys[i] << " "; cout << endl;

/*
// pine seedling data
61 38 93 65 73 51 49 42
*/
FBinomialPartition FTG (Ns, Ys, UseLogPi);
//FBinomialPartition FTG (UseLogPi);
ProduceMRSamples(FTG, n_boxes, n_samples, M, Alb, theSeed, use_f_scale);

  return 0;			// for the main statement
}

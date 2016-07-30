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
#include "FCFN3.hpp"
#include "MRSampler.hpp"

void
ProduceMRSamples(Fobj & f, int n_boxes, int n_samples, 
                 double Alb, unsigned seed, bool use_f_scale) 
{
  ofstream out ("MRS_CFN3.samples");//file to store the i.i.d samples
  clock_t T1 = clock (), T2, T3;
  // Construct theSampler with the chosen target shape object FTG
  MRSampler theSampler (f, n_boxes, Alb, seed, (use_f_scale == 1));
  /* // To print out the partition of the domain
     cout << "Domain Partition: \n" ;
     theSampler.Print_Domain_Partition();
     theSampler.PrintBoxes(0);
     getchar(); // press enter to continue...
  */
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
  //rs_sample.Print(out);

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
  int n_samples = 10000000;
  double Alb = 0.001;// partition until lower bound on Acceptance Prob. is > Alb

  unsigned theSeed = 0;
  
  // theModelSpace options
  // 0: clocked star tree, 1 param. 
  // 1: clocked rooted, 2 param. 
  // 2: unclocked 3 param. 
  // 3: transdimensional.
  unsigned theModelSpace = 2;   


  // Mitochodrial data of Human, Chimpanzee and Gorilla Brown et al, 1982
  //according to Yang, Proc. Roy. Soc. Lon. B (2000) Table 1
  int Cxxx=762;
  int Cxxy=54;
  int Cyxx=38;
  int Cxyx=41;
  //the following two patterns have been incorrectly parsed in the Analysis of
  // Yang 2000 as well as Sainudiin and York 2009
  //int Cyxx=41;
  //int Cxyx=38;

  // other data sets
  //Cxxx=0; Cxxy=0; Cyxx=0; Cxyx=0;//no data -- so target is just the prior
  //Cxxx=300; Cxxy=80; Cyxx=65; Cxyx=55;//Fig4 Yang&Rannala 2005
  Cxxx=884; Cxxy=6; Cyxx=2; Cxyx=3;//Sainudiin and York - H,C,G - Brown et al data
  Cxxx=857; Cxxy=29; Cyxx=4; Cxyx=5;//Sainudiin and York - C,G,O - Brown et al data
  bool UseLogPi = true; int tmpULP; 
  bool use_f_scale = true;
  int PriorType = 0; // unif. prior is default
  double DomainUB = 10.0; // Domain is [1e-10, DomainUB] for each dimension

  if (argc >= 2)
    {
      sscanf (argv[1], "%ui", &theModelSpace);
      if (argc >= 3)
	{
	  sscanf (argv[2], "%i", &n_boxes);
	  if (argc >= 4){
	    sscanf (argv[3], "%i", &n_samples);
	  }
	  if (argc >= 5){
	    cout << "theSeed: " << argv[4] << endl;
	    sscanf (argv[4], "%ui", &theSeed);	    
	  }
	  if (argc >= 6){
	    cout << "Cxxx: " << argv[5] << endl;
	    sscanf (argv[5], "%i", &Cxxx); 
	  }
	  if (argc >= 7) {
	    cout << "Cxxy: " << argv[6] << endl;
	    sscanf (argv[6], "%i", &Cxxy);
	  }
	  if (argc >= 8){
	    cout << "Cyxx: " << argv[7] << endl;
	    sscanf (argv[7], "%i", &Cyxx);
	  }
	  if (argc >= 9){
	    cout << "Cxyx: " << argv[8] << endl;
	    sscanf (argv[8], "%i", &Cxyx);
	  }
	  if(argc >= 10){
	    cout << "UseLogPi: " << argv[9] << endl;
	    sscanf (argv[9], "%i", &tmpULP);
	    UseLogPi = (bool)tmpULP;
	    cout << "UseLogPi: " << argv[9] << endl;
	  }  
	  if(argc >= 11){
	    cout << "PriorType: " << argv[10] << endl;
	    sscanf (argv[10], "%i", &PriorType);
	    //  UseLogPi = (bool)tmpULP;
	    cout << "PriorType: " << argv[10] << endl;
	  }
	  if(argc >= 12){
	    cout << "DomainUB: " << argv[11] << endl;
	    sscanf (argv[11], "%lf", &DomainUB);
	    //  UseLogPi = (bool)tmpULP;
	    cout << "DomainUB: " << argv[11] << endl;
	  }
 
	  if (argc >= 13){
	    cout << "# Usage: " << argv[0] << "<theModelSpace> <n_boxes> <n_samples> "
           << "<seed> <Cxxx> <Cxxy> <Cyxx> <Cxyx> <UseLogPi> <PriorType> " 
           << "<DomainUB>; extra arguments ignored.\n";
	  }

	}
    }

  else { 
    cout << "# Usage: " << argv[0] << "<theModelSpace> <n_boxes> <n_samples> "
         << "<seed> <Cxxx> <Cxxy> <Cyxx> <Cxyx> <UseLogPi> <PriorType> "
         << "<DomainUB>; extra arguments ignored.\n";
  

    cout << "# theModelSpace: " << theModelSpace << "  n_boxes: " 
         << n_boxes << "  n_samples: " << n_samples << "  rng_seed = " 
         << theSeed  << endl; //getchar();
    
  }
   
  int Cin=Cxxx;
  int Cnin=Cxxy+Cyxx+Cxyx;
  cout << theSeed << "  " << UseLogPi << "   " << endl;
  //n_xxx, n_xxy, n_yxx, n_xyx
  cout << Cxxx << " " << Cxxy << " " << Cyxx << " " << Cxyx << endl;
  // CFN phylogenetic Models
  // domain for CFN model under rooted,clocked 3-taxa star tree
  interval DomainCFN3(1e-10, DomainUB); 
  switch (theModelSpace) {
  case 0: {//------CFN model on 3 taxa Star Trees
    //UseLogPi = true; Cin=762; Cnin=133; // Example of Brown et al, 1982.
    FCFN3Star FTG (Cin, Cnin, DomainCFN3, UseLogPi, PriorType);
    cout << "Cin, Cnin: " << Cin << "  " << Cnin << endl;
    
    // Ziheng Yang, Complexity of the Simplest Phylogenetic Model", 2000.
    cout << "MLE = " 
         << -0.25*log((4.0*(double(Cin)/(double(Cin+Cnin)))-1.0)/3.0) << endl;
    
    ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);
    break;}
  case 1: {//------CFN model on 3 taxa Rooted Trees with three topologies
    // UseLogPi = true; 
    FCFN3Rooted FTG (Cxxx, Cxxy, Cyxx, Cxyx, DomainCFN3, UseLogPi, PriorType);
    ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);
    break;}
  case 2: {//------CFN model on 3 taxa Unrooted Trees with one topology
    //UseLogPi = true; 
    FCFN3UnRooted FTG (Cxxx, Cxxy, Cyxx, Cxyx, DomainCFN3, UseLogPi, PriorType);
    //UseLogPi = false; FCFN3UnRooted FTG (76, 5, 4, 3, DomainCFN3, UseLogPi);
    ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);
    cout << "Cin, Cnin: " << Cin << "  " << Cnin << endl;
    break;}
  case 3: {//------Trans-dimensional CFN model on 5 labelled 3 taxa Trees, 
           //------ie 5 topologies [star+3rooted trees+1unrooted]
    //UseLogPi = true; 
    cout << "UseLogPi: " << UseLogPi << endl;
    FCFN3 FTG (Cxxx, Cxxy, Cyxx, Cxyx, DomainCFN3, UseLogPi, PriorType);
    ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);
    break;}
  default:{
    {cerr<< "contains unspecified theModelSpace = " << theModelSpace 
         << " [allowed values are 0,1,2,3] ! using theModelSpace = 3"; 
      theModelSpace=1;}
    //UseLogPi = true; 
    FCFN3Rooted FTG (Cxxx, Cxxy, Cyxx, Cxyx, DomainCFN3, UseLogPi, PriorType);
    ProduceMRSamples(FTG, n_boxes, n_samples, Alb, theSeed, use_f_scale);
    break;}
  }

  return 0;			// for the main statement
}

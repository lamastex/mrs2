/*
* Copyright (C) 2007-2015 Raazesh Sainudiin
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


/*! \file SimpleIO.cpp
\brief Estimate optimally smoothed MAP adaptive histogram estimate via CV
 */

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams

#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data
/* what's in dataprep.hpp
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
*/
//for gsl permutations
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include <limits>

/* heades in Gaussian Carver Example
#include "adaptivehistogram.hpp" 
#include "histevalobj.hpp"
#include "piecewise_constant_function.hpp"
*/
#include "carver_seb.hpp"

#include <time.h>   // clock and time classes
#include <iostream>  // input and output streams
#include <fstream>  // file streams

// to be able to manipulate strings as streams
#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions
#include <iterator>

#include <vector>
#include <map>

#include "testDenCommon.hpp" // to use density testing tools
#include "testDenTools.hpp"
#include "mixture_mvn.hpp" // to use MixtureMVN (Jenny's thesis)


using namespace cxsc;
using namespace std;
using namespace subpavings;
using namespace subpavings::kde;


int main(int argc, char ** argv) 
{
    // ------- prepare to generate some data for the tests -----------
    //const int n=10000;  // number to generate

    clock_t start, end;     // for timing
    double timeTaken;
    int reps=1;
    int dim = 2;
    long unsigned int seed =1234;
    string leafDepthFileName = "ldfn.txt"; // for leaf-depth encoded string specifying the shape of the SRP histogram
    string rangesFileName = "ranges.txt"; // for the range or height value of each leaf node
    string rootBoxFileName = "rootBox.txt";// root box of the histogram
    if (argc > 1) {
      dim = atoi(argv[1]);
    }
    if (argc > 2) {
      rootBoxFileName = argv[2];
    }
    if (argc > 3) {
      leafDepthFileName = argv[3];
    }
    if (argc > 4) {
      rangesFileName = argv[4];
    }
    else cout << "USAGE: ./SimpleIO 2 rootBox.txt ldfn.txt ranges.txt" << endl;
    // Mixture of Normals is made in either case
    long unsigned int MVNseed = seed+9876;// seed for data simulator

    std::ifstream ifs(leafDepthFileName.c_str());
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    // this replaces any end of line character by white space for compatibility with .splitToShape method
    for (int i = 0; i < str.length();i++) if (str[i] == '\n') str[i] = ' '; 
    ifs.close();
    //getchar();

    std::ifstream rfs(rangesFileName.c_str());
    std::vector< real > ranges;
    for (real a; rfs >> a;) ranges.push_back(a);
    rfs.close();

    std::ifstream bfs(rootBoxFileName.c_str());
    ivector  pavingBox(dim);
    for(int i=1; i <= dim; i++) {
       interval a; bfs >> a; pavingBox[i] = a;
    }
    bfs.close();

    // domain is a hypercube read from file
    //ivector pavingBox(dim);
    //interval pavingInterval(0,1);
    //for(int i=1; i <= dim; i++) pavingBox[i] = pavingInterval;
    // a container for the leaf boxes
    //    vector<ivector> Pboxes;
    //    size_t PartSize;
    
    PiecewiseConstantFunction pcf(pavingBox);

    bool successfulInstruction = pcf.splitToShape(str);

    //bool successfulInstruction = pcf.splitToShape("1,2,2");

    //std::vector < real > ranges(5); ranges[0]=5.0; ranges[1]=0.1; ranges[2]=0.2;ranges[3]=0.3; ranges[4]=0.4;
    //pcf.allocateRanges(ranges);

    //std::vector < real > ranges(3); ranges[0]=5.0; ranges[1]=0.1; ranges[2]=0.2;//ranges[3]=0.3; ranges[4]=0.4;
    pcf.allocateRangesToLeaves(ranges);

    cout << "Level string for new partition is "
             << pcf.getLeafLevelsString() << endl;
   //cout << endl << pcf << endl;
   pcf.outputToStreamTabs(cout);


    // containers for stuff we will be storing 
    size_t intN = 10;//00000;
    //MixtureMVN* mixMVNptr = makeMixture(dim , MVNseed); 
    MixtureMVN* mixMVNptr = makeStandard(dim , MVNseed); 
    

    clock_t starttime = clock();
        

    // make the below modular also!!!
    // a container for our histograms
    //std::vector< subpavings::AdaptiveHistogram* > hists;
    // a container for our PCFs of histograms
    //std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
    //getchar();

    clock_t endtime = clock();	
    double timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
    cout << "time to get prior-selected adaptive hist = " << timingStarts << endl;

    // L1 error calculations
    //getting number of leaves in optimal MAP estimate
    size_t NoOfLeaves = pcf.getRootLeaves();
    cout << "\nOptMAP estimate with " << NoOfLeaves << " leaves"<< endl; //getchar();
    //get quasi random points in the box 
    ivector box = pcf.getRootBox();
    std::vector < std::vector < real > > qrPts;
    getQuasiRandomPoints( box, qrPts, intN);
    // get points from true density 
    //std::vector < std::vector < real > > intPts;
    //mixMVNptr->prn(intPts, intN);
    //get MCMC histogram densities at the remaining integration points and log results
    std::vector < real > estDensitiesOptMAP_QR;
    getPCFDensities(pcf, qrPts, estDensitiesOptMAP_QR);
    
    //get true densities at the qr points 
    std::vector < real > trueIntPtDensities_QR;
    getTrueDensities(*mixMVNptr, qrPts, trueIntPtDensities_QR);
    for(int i=0; i < qrPts.size(); i++) { cout << qrPts[i]  << '\t' << trueIntPtDensities_QR[i] << '\t' << estDensitiesOptMAP_QR[i] << '\n';} cout << endl;
    real boxVol = realVolume(box);
    real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensitiesOptMAP_QR);

    cout << endl << "estimated L1 error = " << estL1_QR << endl;
    //to free all the contents of pcfs at the end
   //if (NULL != &pcf) delete *pcf;

    //delete the data generator 
    delete mixMVNptr;	

    return 0;

} // end of program


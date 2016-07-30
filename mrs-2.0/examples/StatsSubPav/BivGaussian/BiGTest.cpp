/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
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


/*! \file BiGTest.cpp
\brief Testing StatsSubPavings (aka SPSnodes) with Bivariate Gaussian data
 */

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams

#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data

using namespace cxsc;
using namespace std;
using namespace subpavings;


int main()
{
    // ------- prepare to generate some data for the tests -----------

    // set up a random number generator for bivariate gaussian rvs
    const gsl_rng_type * T;
    gsl_rng * r;

    int i;
    const int n=10000;  // number to generate
    double sigma_x=1;   // distribution parameter
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    //create a generator chosen by the environment variable GSL_RNG_TYPE

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    string samplesFileName; // for samples
    string outputFileName;// for output file
    ofstream oss;         // ofstream object
    oss << scientific;  // set formatting for input to oss
    oss.precision(5);

    double *x;
    double *y;

    x= new double[n];   // make x and y in dynamic memory
    y= new double[n];   // (so they must be freed later)

    double* itx;
    double* ity;


    // get n random variates chosen from the bivariate Gaussian
    // distribution with mean zero and given sigma_x, sigma_y.
    for (i = 0; i < n; i++)
    {
        gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                   rho, &x[i], &y[i]);

    }

    // free the random number generator
    gsl_rng_free (r);

    itx = &x[0];
    ity = &y[0];

    // create a name for the file to use
    samplesFileName = "bgSamples.txt";
    // output the sample data
    oss.open(samplesFileName.c_str());         // opens the file

    for(i=0; i<n; i++) {
        oss << (*itx) << "  " << (*ity);
        if (i<n-1) oss << endl; // new line if not final line
        itx++;
        ity++;
    }
    oss << flush;
    oss.close();

    cout << "Samples output to " << samplesFileName << endl << endl;

    clock_t start, end;     // for timing
    double timeTaken;

    bool successfulInsertion = false;
    bool successfulPQSplit = false;

    // ------ example to create one histogram with splitting value ----
    // --------------------entered by user ----------------------------

    // get a count of the lines in the txt file
    int dataCount = countLinesInTxt(samplesFileName);
    int myK = 0;

    // tell user how many lines there are in the file
    cout << "The file " << samplesFileName << " has " << dataCount
            << " lines in it" << endl << endl;
    // get a parameter for k
    cout << "Enter a parameter for your splitting criteria here please:  ";
    cin >> myK;
    cout << endl << endl; // myK has been input


    // make an Adaptive Histogram object with no specified box and, by default,
    // holdAllStats = false so that the underlying rootPaving managed by the
    // histogram will not maintain all available stats, only counts
    AdaptiveHistogram myHistFirst;

    start=clock();
    // clock running

    // make the function object to decide whether to split.
    // aim to get max myK data members in each box, default minimum
    // number of points allowed in each box is 0.
    SplitOnK splitK(myK);

    int dim = 2;
	size_t headerlines = 0;

    // insert the data on by one, checking whether to split on each insertion
    successfulInsertion = myHistFirst.insertRvectorsFromTxtOrd(samplesFileName,
            splitK, dim, headerlines, NOLOG); // no logging

    end=clock();

    timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
    cout << "Computing time : " <<timeTaken<< " s." << endl;

    // only do more if some data was fed in
    if(successfulInsertion) {

        // create a name for the file to output
        outputFileName = "BivGaussianFirst.txt";
        // To realize a file output
        myHistFirst.outputToTxtTabs(outputFileName);
    }

    // end of example for histogram with splitting value input by user

    successfulInsertion = false;
    successfulPQSplit = false;

    // example to create one histogram with pulse data and a priority
    // ---------- split to give a minimum number of bins -----------

    // make an Adaptive Histogram object with no specified box
    // or splitting value (and holdAllStats again defaults to false),
    // with the same data
    AdaptiveHistogram myHistSecond;

    start=clock();
    // clock running

    // put in the data in a 'pulse' with no splitting, ie one big box
    successfulInsertion = myHistSecond.insertRvectorsFromTxt(samplesFileName);

    if (successfulInsertion) {

        // set up function objects for a priority split

        // function object to compare nodes on count
        // ie split node with largest count first
        CompCount compCount;

        // function object to split until number of leaves is >= minLeaves
        size_t minLeaves = 50;
        CritLeaves_GTE critLeavesGTE(minLeaves);

        /* minimum points to use when splitting.
        A node will not be splittable if either child would then have
        < minPoints of data associated with it. */
        size_t minPoints = 1;

        // do the priority split
        successfulPQSplit = myHistSecond.prioritySplit(compCount,
                critLeavesGTE, NOLOG, minPoints); // no logging
    }

    end=clock();

    timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
    cout << "Computing time : " <<timeTaken<< " s." << endl;


    if(successfulPQSplit) { // only do more if split was successful

        // create a name for the file to output
        outputFileName = "BivGaussianSecond.txt";
        // To realize a file output
        myHistSecond.outputToTxtTabs(outputFileName);

    }

    delete x;   // free dynamic memory used for x and y
    delete y;

    return 0;

} // end of bivariate gaussian test program

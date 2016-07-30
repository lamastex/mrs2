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


/*! \file ScoreTest.cpp
\brief Testing StatsSubPavings (aka SPSnodes) with COPERR and AIC scoring data
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
    const int n=100;    // number to generate
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
    bool successfulPQMerge = false;

    // example to create one histogram with pulse data and a priority
    // ---------- split to give a minimum number of bins -----------

    cout << endl << endl;
    cout << "Start of first example:" << endl;
    cout << "Priority split on COPERR to give minimum number of bins" << endl;

    // Use default constructor for histogram.  By default,
    // holdAllStats = false so that the underlying rootPaving managed by the
    // myHistFirst will not maintain all available stats, only counts
    AdaptiveHistogram myHistFirst;

    start=clock();
    // clock running

    // put in the data in a 'pulse' with no splitting, ie one big box
    successfulInsertion = myHistFirst.insertRvectorsFromTxtOrd(samplesFileName);

    if (successfulInsertion) {

        // with min volume check
        // the minimum volume of a splittable node is minVolB(log n)^2/n
        double minVolB = 1.0;

        // prepare function objects to do a priority split

        // function object to compare nodes on cont. to EMP sum under COPERR
        // ie split node with EMP contribution (COPERR) first
        CompEMPSumChangeCOPERR nodeCompEMPCOPERR;

        // function object to split until number of leaves is >= minLeaves
        int minLeaves = 20;
        CritLeaves_GTE critLeavesGTE(minLeaves);

        // do the priority split
        /* minimum points minPoints to use when splitting defaults to zero.
        A node will not be splittable if either child would then have
        < minPoints of data associated with it. */
        successfulPQSplit = myHistFirst.prioritySplit(nodeCompEMPCOPERR,
              critLeavesGTE, NOLOG, minVolB); // no logging
    }

        end=clock();

        timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
        cout << "Computing time : " <<timeTaken<< " s." << endl;


    if(successfulPQSplit) { // only do more if split was successful

        // create a name for the file to output
        outputFileName = "ScoreTestFirst.txt";
        // To realize a file output
        myHistFirst.outputToTxtTabsWithEMPs(outputFileName);

        myHistFirst.outputGraphDot();
        // output will go to a dot file
        // this will then be made into a png image
        // console output will tell you where png image is

        // optional - print out the scores

        // create a penalty function object to use to give us the total
        // score (EMP + PEN) under COPERR
        PenLeaves penC;

        real scoreCOPERR = myHistFirst.getScoreCOPERR(penC, true); //verbose
        cout << "End of first example:" << endl;
        cout << "The total score (EMP+PEN) under COPERR is " << scoreCOPERR
                            << endl;
        real scoreAIC = myHistFirst.getEMPScoreAIC();
        cout << "The EMP score under AIC is " << scoreAIC << endl;
    }

    // example to create one histogram with pulse data and a priority
    // ---------- split to minimise COPERR score -----------


    successfulInsertion = false;
    successfulPQSplit = false;

    // example to create one histogram with pulse data and a priority
    // ---------- split to give a minimum number of bins -----------

    cout << endl << endl;
    cout << "Start of second example:" << endl;
    cout << "Priority split on AIC give minimum number of bins" << endl;

    // make an Adaptive Histogram object with no specified box
    // or splitting value (and holdAllStats again defaults to false),
    // with the same data
    AdaptiveHistogram myHistSecond;

    start=clock();
    // clock running

    // put in the data in a 'pulse' with no splitting, ie one big box
    successfulInsertion = myHistSecond.insertRvectorsFromTxtOrd(samplesFileName);

    if (successfulInsertion) {

        // prepare function objects to do a priority split

        // function object to compare nodes on cont. to EMP sum under COPERR
        // ie split node with EMP contribution (COPERR) first
        CompEMPSumChangeAIC nodeCompEMPAIC;

        // function object to split until number of leaves is >= minLeaves
        int minLeaves = 500;
        CritLeaves_GTE critLeavesGTE(minLeaves);

        /* minimum points minPoints to use when splitting defaults to zero.
        A node will not be splittable if either child would then have
        < minPoints of data associated with it.
        minVolB, the multiplier for (log n)^2/n to determine the minimum
        volume of a splittable node, also defaults to 0.0. */
        successfulPQSplit = myHistSecond.prioritySplit(nodeCompEMPAIC,
              critLeavesGTE, NOLOG); //no logging
    }

    end=clock();

    timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
    cout << "Computing time : " <<timeTaken<< " s." << endl;


    if(successfulPQSplit) { // only do more if split was successful

        // create a name for the file to output
        outputFileName = "ScoreTestSecond.txt";
        // To realize a file output
        myHistSecond.outputToTxtTabsWithEMPs(outputFileName);

        // optional - print out the scores
        // make a penalty function to give total (EMP + PEN) score under AIC
        double ca = 1.0;
        double alpha = 0.5;
        double ra = 2.0;
        PenAIC1 penA(ca, alpha, ra); // penalty for AIC

        real scoreCOPERR = myHistSecond.getEMPScoreCOPERR();
        cout << "End of second example:" << endl;
        cout << "The EMP score COPERR is " << scoreCOPERR << endl;
        real scoreAIC = myHistSecond.getScoreAIC(penA, true);
        cout << "The total score (EMP+PEN) under AIC is " << scoreAIC << endl;
    }

    // end of example to create one histogram with pulse data and a priority
    // ---------- split to minimise AIC score -----------

    // example to merge up histogram,prioritising merge on
    // ---------------- minimising increase in AIC score -------------

    if(successfulPQSplit) { // only do more if split was successful


        cout << endl << endl;
        cout << "Example to merge the second histogram upwards" << endl;

        // merge node with largest contribution to EMP sum first
        CompEMPSumChangeMergeAIC nodeCompEMPMergeAIC;

        int minLeaves = 0; // this will merge right back to the root box
        CritLeaves_LTE critLeavesLTE(minLeaves);

        successfulPQMerge = myHistSecond.priorityMerge(nodeCompEMPMergeAIC,
            critLeavesLTE, TXT); // with logging

        //try critAICMergeChangeOrLeaves_LTE instead of critLeaves_LTE

        if (successfulPQMerge) {
            // optional - print out the scores
            real scoreAIC = myHistSecond.getEMPScoreAIC();
            cout << "After merging:" << endl;
            cout << "the EMP score under AIC is " << scoreAIC << endl;
        }

    }

    cout << endl;
    cout << "End of examples" << endl << endl;

    // end of examples

    delete x;   // free dynamic memory used for x and y
    delete y;

    return 0;

} // end of COPERR and AIC scoring test program

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


/*! \file Averaging.cpp
\brief Testing CollatorSPSnodes with Bivariate Gaussian data
*/

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams

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

    // ----------------   example to create and ------------------
    //---------------- collate multiple histograms -------------------

    // make a box: the same box will be used by all histograms
    // so should be big enough for all of them
    int d = 2; // dimensions
    ivector pavingBox(d);
    interval dim1(-5,5);
    interval dim2(-5,5);
    pavingBox[1] = dim1;
    pavingBox[2] = dim2;

    // make a collation object, empty at present
    AdaptiveHistogramCollator coll;

    // the number of histograms to generate
    int numHist = 10;

    // for loop to generate histograms and add to collation
    for (int j=1; j<=numHist; j++) {

        //get n random variates chosen from the bivariate Gaussian
        // distribution with mean zero and given sigma_x, sigma_y.

        RVecData theData;   // a container for all the points generated

        // make a sample
        for (int i = 0; i < n; i++) {

            rvector thisrv(d);
            double x = 0;
            double y = 0;

            gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y,
                                    rho, &x, &y);
            thisrv[1] = x;
            thisrv[2] = y;

            // put points generated into container
            theData.push_back(thisrv);

        }  // data should be in theData


        // make an Adaptive Histogram object with a specified box.  By default,
        // holdAllStats = false so that the underlying rootPaving managed by the
        // myHistFirst will not maintain all available stats, only counts
        AdaptiveHistogram myHist(pavingBox);

        // find k, the maximum number of data members
        // to be allowed in each box of the histogram
        // as a function of j and n
        // applying SEB heuristics for k to satisfy k/n -> 0 as n -> +oo
        int k_int = (int(log2(double(n)))*2*j);

        cout << "Splitting with k = " << k_int << endl;

        bool successfulInsertion = false;
        bool successfulPQSplit = false;

        // make the function object to get max myK data members in each box
        SplitOnK splitK(k_int);

        // insert data into the histogram, splitting as we go, no logging
        successfulInsertion = myHist.insertFromRVec(theData, splitK, NOLOG);

        // only do more if some data was fed in
        if(successfulInsertion) {

            // create a name for the file to output
            string fileName = "BivGaussian";
            //convert j to a string
            std::ostringstream stm2;
            stm2 << j;
            // add the stringed j to the filename
            fileName += stm2.str();
            fileName += ".txt"; // and finish the filename

            // To realize a file output
            myHist.outputToTxtTabs(fileName);

            // add the histogram to the collection
            coll.addToCollation(myHist);

            // optional- create graph output
            // myHist.outputGraphDot();
        }

    } // end of for loop creating histograms

    // free the random number generator
    gsl_rng_free (r);

    string collfileName = "CollatorHistogram.txt";
    coll.outputToTxtTabs(collfileName); // output the collation to file

    // optional - create graph output - don't do for lots of leaves!
    //coll.outputGraphDot();

	//  Make an average
	string avgfileName = "AverageBG.txt";     // provide a filename

	AdaptiveHistogramCollator avColl = coll.makeAverage();
	avColl.outputToTxtTabs(avgfileName);  // output the average to file

    // ---- end of example to create and collate multiple histograms -----

    return 0;

} // end of averaging test program

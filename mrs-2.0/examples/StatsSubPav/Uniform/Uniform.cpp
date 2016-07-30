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


/*! \file Uniform.cpp
\brief Testing histograms with uniform distribution examples
*/

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams

#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data

using namespace cxsc;
using namespace std;
using namespace subpavings;

bool collateFromRVecSplitPQ(AdaptiveHistogramCollator& coll,
					size_t samplesize, size_t numberSamples,
					RVecData& allData,
                    ivector pavingBox, const NodeCompObj& compTest,
                    const HistEvalObj& he,
                    size_t minChildPoints, double minVolB);

int main()
{
    // ------- prepare to generate some data for the tests -----------

    // set up a random number generator for uniform rvs
    const gsl_rng_type * T;
    gsl_rng * r;

    //create a generator chosen by the environment variable GSL_RNG_TYPE

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    string fileName; // create a name for the file to use

/*
    double *x;
    double *y;

    x= new double[n];
    y= new double[n];

    double* itx;
    double* ity;

*/

    // ----------------   example to create and ------------------
    //---------------- collate multiple histograms -------------------

    // make a box: the same box will be used by all histograms
    // so should be big enough for all of them

    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
/*
    ivector pavingBox(2);
    interval pavingInterval1(-5,5);
    interval pavingInterval2(-5,5);
    pavingBox[1] = pavingInterval1;
    pavingBox[2] = pavingInterval2;
*/

    int totalDataPoints = 100000; // total points from random number generator

    RVecData allData;   // a container for all the points generated

    // make a simulated data set allData to sample from
    for (int i = 0; i < totalDataPoints; i++) {

        real thisreal = gsl_rng_uniform(r);
        rvector thisrv(thisreal);
        // put points generated into container

        allData.push_back(thisrv);

    }  // data to draw from should be in allData

    cout << "Simulated data set to sample from has been created "
         << endl << endl;



    // the number of histograms to generate
    // each histogram is from a parametrically bootstrapped sample of data
    // points from allData
    int numSamples = 10;
    const int samplesize=10000;  // number of samples to take for each

    // make a collation object, empty at present
    AdaptiveHistogramCollator coll;

    // set up function objects for priority queue splitting

    // node comparison on count k
    CompCount compCount;

    // and stopping when the largest count for any splittable node is
    // <= myK
    size_t myK = 100; // a value to cease to split the histograms on
    CritLargestCount_LTE critLargestCountLTE(myK);

    /* minPoints is the minumum number of points in any node.
    A node will not be splittable if splitting that node would give at least
    one child with < minPoints of data associated with it. */
    size_t minPoints = 0;

    /* minVolB is the multiplier for (log n)^2/n to determine the minumum
    volume for a splittable node.  A node with
    volume < minVolB(log n)^2/n will not be split. */
    double minVolB = 0.0;

    /* do the creation, collation and averaging of numSamples histograms
    using a priority queue to make the individual histograms */
    bool success = collateFromRVecSplitPQ(coll,
					samplesize, numSamples,
                    allData, pavingBox, compCount, critLargestCountLTE,
                    minPoints, minVolB);
    if (success) {
            string fileName = "CollatorHistogram.txt";
            coll.outputToTxtTabs(fileName); // output the collation to file

            //  Average the histograms
            fileName = "AverageHistogram.txt";     // provide a filename

            coll.outputAverageToTxtTabs(fileName);  // output the average to file
    }

    else    cout << "Failed to create collation over histograms" << endl;

    //Raaz raw it up -- bring the loop body of averageFromRVecSplitPQ [[chased to collation]] here
    // for loop to generate histograms and add to collation
/*
    for (int j=1; j<=numberSamples; j++) {

        // make an Adaptive Histogram object with a specified box
        AdaptiveHistogram myHist(pavingBox);
        //AdaptiveHistogram myHistCurr (pavingBox);
        //AdaptiveHistogram myHistProp (pavingBox);

        if (indImmedSplit == 1) { // doing immediate splitting

            successfulInsertion = myHist.insertSampleFromRVec(samplesize,
                rgsl, rv, boolTest, test);
        }

        if (indImmedSplit == 0) { // doing priority queue splitting

            successfulInsertion = myHist.insertSampleFromRVec(samplesize,
                rgsl, rv);

            bool successfulPQSplit;

            if (successfulInsertion) {

                successfulPQSplit = myHist.prioritySplit(compTest, stopTest,
                                        test);
            }

            successfulInsertion = successfulInsertion && successfulPQSplit;

        }

        // only do more if some data was fed in
        if(successfulInsertion) {

            // create a name for the file to output
            fileName = "Hist";
            //convert j to a string
            std::ostringstream stm2;
            stm2 << j;
            // add the stringed j to the filename
            fileName += stm2.str();
            fileName += ".txt"; // and finish the filename


            // To realize a file output
            myHist.outputToTxtTabs(fileName);

            // add the histogram to the collection represented by this
            addToCollation(myHist);

            countIn++; // increment the counter

        }


    }
*/
    gsl_rng_free (r);
    return 0;

} // end of bivariate gaussian test program

// Make a collated histogram from an RVec
// priority queue splitting
bool collateFromRVecSplitPQ(AdaptiveHistogramCollator& coll,
					size_t samplesize, size_t numberSamples,
					RVecData& allData,
                    ivector pavingBox, const NodeCompObj& compTest,
                    const HistEvalObj& he,
                    size_t minChildPoints, double minVolB)
{
	gsl_rng * rgsl = NULL;

    try {
		
		bool cancontinue = (!allData.empty());
		// cancontinue will be false if there was a problem getting data points
		// if cancontinue is true data should contain at least some data points

		bool retValue = false;

		if (cancontinue) {

			int countIn = 0; // track the number of histograms made and added

			// set up a random number sampler
			const gsl_rng_type * tgsl;

			// set the library variables *gsl_rng_default and
			// gsl_rng_default_seed to default environmental vars
			gsl_rng_env_setup();

			tgsl = gsl_rng_default; // make tgsl the default type
			rgsl = gsl_rng_alloc (tgsl); // set up with default seed

			std::string fileName; // a name for the files to use


			// for loop to generate histograms and add to collation
			for (size_t j=1; j<=numberSamples; j++) {

				bool successfulInsertion; // recognise successes

				// make an Adaptive Histogram object with a specified box
				AdaptiveHistogram myHist(pavingBox);
				
				successfulInsertion = myHist.insertSampleFromRVec(samplesize,
					rgsl, allData);

				bool successfulPQSplit;

				if (successfulInsertion) {

					successfulPQSplit = myHist.prioritySplit(compTest, he,
													NOLOG, minChildPoints, minVolB);
				}

				successfulInsertion = successfulInsertion && successfulPQSplit;

				// only do more if some data was fed in
				if(successfulInsertion) {

					// create a name for the file to output
					fileName = "Hist";
					//convert j to a string
					std::ostringstream stm2;
					stm2 << j;
					// add the stringed j to the filename
					fileName += stm2.str();
					fileName += ".txt"; // and finish the filename

					// To realize a file output
					myHist.outputToTxtTabs(fileName);

					// add the histogram to the collection represented by this
					coll.addToCollation(myHist);

					countIn++; // increment the counter
				}
			} // end of for loop creating histograms

			if (countIn == numberSamples) {

				retValue = true;
			}
			else { // did not add required number of histograms
				std::cerr << "Problem in collateFromRVec(): check "
						<< "console for error reports " << std::endl;
			}
			// free the random number generator
			gsl_rng_free (rgsl);
			rgsl = NULL;
			
		}

		return retValue;
	}
    catch (exception const&) {
        try {
			if (NULL != rgsl) {
				gsl_rng_free(rgsl);
				rgsl = NULL;
			} 
			// free the random number generator
		}
		catch (exception const& ee) {} // catch and swallow
		throw; // rethrow original exception
    }
}



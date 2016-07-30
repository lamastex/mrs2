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


/*! \file LevyTest.cpp
\brief Testing StatsSubPavings Levy2D data
*/

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams

#include "histall.hpp"  // headers for the histograms

#include <gsl/gsl_qrng.h>       // types needed by MRSampler.hpp
#include <gsl/gsl_randist.h>
#include "Fobj.hpp"     // to be able to use the Levy function objects
#include "FLevy2D.hpp"
#include "MRSampler.hpp"    // to be able to do MRS rejection sampling


using namespace cxsc;
using namespace std;
using namespace subpavings;

bool collateFromRSSampleSplitPQ(AdaptiveHistogramCollator& coll,
					size_t samplesize, size_t numberSamples,
					const RSSample rss,
                    ivector pavingBox, const NodeCompObj& compTest,
                    const HistEvalObj& he,
                    size_t minChildPoints, double minVolB, int label);

int main (int argc, char **argv)
{

    // example to average 10 samples from a 2-d Levy shape

    ios::sync_with_stdio (); // so iostream works with stdio
    cout << SetPrecision (20, 15);  // Number of mantissa digits in I/O


    int n_dimensions = 2;
    int n_boxes = 1000;
    int n_samples = 100000;
    double Alb = 1.0;// partition until lower bound on Acceptance Prob.>Alb
    unsigned theSeed = 12345;


    bool use_f_scale = true;

    cout << "# n_dimensions: " << n_dimensions << "  n_boxes: " << n_boxes
        << "  n_samples: " << n_samples << "  rng_seed = " << theSeed
        << endl; //getchar();

    //Parameters specific to the Levy target
    real Temperature = 40.0;
    real Center1 = 1.42513;
    real Center2 = 0.80032;
    real GlobalMax = 176.14;
    real DomainLimit = 10.0;    //0.999999999999999;
    bool UseLogPi = false; // log scale won't work naively

    // make the function object
    FLevy2D F_Levy_Temp_2D(Temperature, GlobalMax,
                        Center1, Center2, DomainLimit, UseLogPi);

    // create the sampler
    MRSampler theSampler (F_Levy_Temp_2D, n_boxes, Alb, theSeed,
                        (use_f_scale == 1));

    // produce the samples (n_sample samples should be produced)
    RSSample rs_sample; // object for the sample

    theSampler.RejectionSampleMany (n_samples, rs_sample);


    // make a box: the same box will be used by all histograms
    // so should be big enough for all of them, so use the function domain
    // set up the domain list
    ivector pavingBox(2);
    interval dim1(-DomainLimit, DomainLimit);
    interval dim2(-DomainLimit, DomainLimit);
    pavingBox[1] = dim1;
    pavingBox[2] = dim2;

    size_t samplesize = 10000; // number of samples to take from the RSSample

    // the number of histograms to generate
    size_t numSamples = 10;

    // make a collation object, empty at present
    AdaptiveHistogramCollator coll;

    // set up objects for priority queue splitting

    // node comparison using count of data points associated with nodes
    CompCount compCount;

    // stopping on smallest volume criteria for splittable nodes
    double vol = 0.05;
    CritSmallestVol_LTE critSmallestVol(vol);

    /* A node is not splittable if splitting that node would give at least
    one child with < minPoints of data associated with it.*/
    size_t minPoints = 0;

    /* minVolB is the multiplier for (log n)^2/n to determine the minumum
    volume for a splittable node where n is total points in subpaving.
    A node with volume < minVolB(log n)^2/n will not be splittable. */
    double minVolB = 0.0;
	
	int label = 0;

    bool success = collateFromRSSampleSplitPQ(coll, 
					samplesize, numSamples,
                    rs_sample, pavingBox, compCount,
                    critSmallestVol, minPoints, minVolB,
					label);

	if (success) {
            string fileName = "CollatorHistogram.txt";
            coll.outputToTxtTabs(fileName); // output the collation to file

            //  Average the histograms
            fileName = "AverageHistogram.txt";     // provide a filename

            coll.outputAverageToTxtTabs(fileName);  // output the average to file
    }

    else    cout << "Failed to create collation over histograms" << endl;

    return 0;

} // end of Levy test program

// Make a collated histogram from an RSSample
// priority queue splitting
bool collateFromRSSampleSplitPQ(AdaptiveHistogramCollator& coll,
					size_t samplesize, size_t numberSamples,
					const RSSample rss,
                    ivector pavingBox, const NodeCompObj& compTest,
                    const HistEvalObj& he,
                    size_t minChildPoints, double minVolB, int label)
{
	gsl_rng * rgsl = NULL;

    try {
		// container to put the rvectors into
		RVecData allData;

		//get the container of rvectors
		//use getRvectorsFromRSSample to put rvectors from labeled points in
		// rss.Samples into allData where the labeled point label matches label
		size_t numberFound = getRvectorsFromRSSample(allData, rss, label);

		bool cancontinue = (numberFound > 0);
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



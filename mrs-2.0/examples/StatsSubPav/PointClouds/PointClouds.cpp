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


/*! \file PointClouds.cpp
\brief Testing StatsSubPavings (aka SPSnodes) with Bivariate Gaussian / read-in data
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

using namespace cxsc;
using namespace std;
using namespace subpavings;


int main(int argc, char ** argv) 
{
    // ------- prepare to generate some data for the tests -----------
    int i;
    const int n=10000;  // number to generate
    double sigma_x=1;   // distribution parameter
    double sigma_y=1;   // distribution parameter
    double rho=0;       // x and y uncorrelated

    string samplesFileName; // for samples
    string pqFileNameBase("pcAdh_");
    string burstsFileBaseName;

    // create a name for the file to use
    //samplesFileName = "bgSamples.txt";
    //burstsFileBaseName = "dataFiles/A";
    //burstsFileBaseName = "dataFiles/B5K";
    //burstsFileBaseName = "dataFiles/B_n5000_d3_k3_m500_";
    //burstsFileBaseName = "dataFiles/B_n3000_d2_k2_m300_";
    //burstsFileBaseName = "dataFiles/B5_n3000_d2_k2_m300_";
    burstsFileBaseName = "dataFiles/tst/M3_B_n5000_d3_k3_m500_";
    if (argc == 2) {
      burstsFileBaseName = argv[1];
    }
    size_t B = 9;//20;//number of bursts
    //burstsFileBaseName = "dataFiles/B_n1000_d2_k2_m100_";
    samplesFileName = burstsFileBaseName+"0.txt";//"dataFiles/B2K0.txt";

    clock_t start, end;     // for timing
    double timeTaken;

    bool successfulInsertion = false;
    bool successfulPQSplit = false;
    // set to doe Carver PQ + SEBPB for posterior maximization
    bool CarvingMaxPosterior=true;// false means using SEB PQ with minChildPoints
    bool TempIterate=true;//true;
    // ------ example to create one histogram with splitting value ----
    // --------------------entered by user ----------------------------

    // get a count of the lines in the txt file
    int dataCount = countLinesInTxt(samplesFileName);
    size_t headerlines = 0;


    // tell user how many lines there are in the file
    cout << "The file " << samplesFileName << " has " << dataCount
            << " lines in it" << endl << endl;


    // set up prior distribution object
    //LogCatalanPrior logPrior;
    LogTemperaturePrior logPrior(10.0);// n=3000,d=2: 1/100 = undersmoother, 10 is oversmoothed/raggedy, (1/10,1.0) is ok

    // a container for our histograms
    std::vector< subpavings::AdaptiveHistogram* > hists;
    // a container for our PCFs of histograms
    std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
    // a container for likelihoods of small bursts given density under first big burst
    std::vector< cxsc::real > LklBursts;
    cxsc::rvector consecL1D = rvector(1,B);
    cxsc::rmatrix allPrsL1D = rmatrix(0,B,0,B);
    
    {
      //-----------------------------------------------------------------------
      // begin get data from first big burst
      //-----------------------------------------------------------------------
      AdaptiveHistogram adhA0Raw; // let it make its own root box
      //bool successfulInsertion = adhA0.insertFromRVec(simdata);
      // put in the data in a 'pulse' with no splitting, ie one big box
      bool insertFirstBigBurst=true;
      bool insertSmallBurstsAlso=true;
      size_t headerlines = 0;
      RVecData RawData;//container to keep the data from first burst or later bursts
      if(insertFirstBigBurst)
      {
        readRvectorsFromTxtOrd (RawData, samplesFileName, headerlines);
        //cout << "size of RawData = " << RawData.size() << endl; getchar(); 
        assert(RawData.size()>0); assert(RawData[0]==d);
      }
      if(insertSmallBurstsAlso)
      {
        //read in the smaller bursts
        std::string burstFileName;
        for (int burstI=1; burstI<=B; burstI++)
        {
          std::string burstHistFileName;
          {
            ostringstream oss;
            oss << burstsFileBaseName << burstI << ".txt";
            burstFileName = oss.str();
            ostringstream oss1;
            oss1 << burstsFileBaseName << "hist_" << burstI << ".txt";
            burstHistFileName = oss1.str();
          }
          // get a count of the lines in the txt file
          int dataCount = countLinesInTxt(burstFileName);
          // tell user how many lines there are in the file
          cout << "The file " << burstFileName << " has " << dataCount
              << " lines in it" << endl;
          //RVecData burstData;
          //readRvectorsFromTxtOrd (burstData, burstFileName, headerlines);
          //cout << "size of burstData = " << burstData.size() << endl;// getchar(); 
          //assert(burstData.size()>0); assert(burstData[0]==d);
          //successfulInsertion = false;
          //successfulInsertion = adhA0Raw.insertFromRVec(burstData);//insert transformed data
          //if (!successfulInsertion) throw std::runtime_error("Failed to insert burst data");
        readRvectorsFromTxtOrd (RawData, burstFileName, headerlines);
        //cout << "size of RawData = " << RawData.size() << endl; getchar(); 
        }
      }
      //successfulInsertion = adhA0Raw.insertRvectorsFromTxt(samplesFileName);
      successfulInsertion = false;
      successfulInsertion = adhA0Raw.insertFromRVec(RawData);//insert all raw data
      RawData.clear();
      adhA0Raw.setHoldAllStats(true);
      if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
      size_t n = adhA0Raw.getRootCounter();
      size_t d = adhA0Raw.getDimensions ();
      BigDataCollection rawDataCollection = adhA0Raw.getDataCollection();
      rvector mean = adhA0Raw.getRootPavingMean();//mean of data
      //cout << mean << endl; getchar();
      //for(size_t i=1; i<=d; i++) cout << mean[i] << endl; getchar();
      //cov(i,j) is at index i*d+j in the returned RealVec
      RealVec varcov = adhA0Raw.getRootPavingVarCovar();//var-covar of data
      //cout << varcov << endl; getchar();
      cxsc::rmatrix V = rmatrix(1,d,1,d);
      for(size_t i=0; i<d; i++) for(size_t j=0; j<d; j++) V[i+1][j+1]=varcov[i*d+j];
      //cout << V << endl; 
//cout <<"n = " << n << endl; getchar();
        
      RVecData transformedData;//container to keep the transformed data from first burst
      
      //RealVec r(d, cxsc::real(0.0));
      //int counter=0;
      for (BigDataConstItr it = rawDataCollection.begin();
                        it != rawDataCollection.end();
                        ++it) 
      {
        rvector rawdata = *it;        
        //counter++; getchar(); cout << counter << endl << rawdata << endl;
        // transformation by centered non-rotating rescaling        
        //cout << rawdata << endl; getchar();
        for (int i = 1; i <= d; ++i) {
            //cout << rawdata[i] << '\t' << mean[i] << '\t' << rawdata[i]-mean[i] << endl; getchar();
            rawdata[i] -= mean[i];
            rawdata[i] /= sqrt(V[i][i]);
        }// end i-loop
        //cout << rawdata << endl; getchar();
        transformedData.push_back(rawdata);
      } // end iteration through data

      adhA0Raw.clearAllHistData();//clear the raw data
      //-----------------------------------------------------------------------
      // end of get data from first big burst
      //-----------------------------------------------------------------------

      //-----------------------------------------------------------------------
      // begin getting best histogram for first big burst
      //-----------------------------------------------------------------------
      int rep=0;

      std::string postFileName;
      {
        ostringstream oss;
        oss << "Posteriors" << pqFileNameBase << "_d" << d << "_n" << n << "_r" << (rep+1) << ".txt";
        postFileName = oss.str();
      }
      std::string checkPostFileNameBase;
      {
        ostringstream oss;
        oss << "CheckPosterior" << pqFileNameBase << "_d" << d  << "_n" << n  << "_r" << (rep+1) << ".txt";
        checkPostFileNameBase = oss.str();
      }

      int precPQ = 5;

	bool stopOnMaxPosterior = true;//bool for carver PQ
	long unsigned int seed = 1234;
	int keep = 1;
	int chooseStarts = 10;
	size_t minPoints = 1;
	clock_t starttime = clock();
       
        
        //-----------------------------------------------------------------------
        // optimizing over smoothing parameter Temp - 
        //                             getting best histogram for first big burst
        //-----------------------------------------------------------------------
        // n=3000,d=2: 1/100 = undersmoother, 10 is oversmoothed/raggedy, (1/10,1.0) is ok
        double t_lo=0.001;//lowest temperature in search
        double t_hi=1.0;//10.0;//highest temperature in search
        int TempIterations=0;
        int MaxTempIterations=20;
        double t_Delta=(t_hi-t_lo)/double(MaxTempIterations-1);//change in temperature
        std::vector<double> Temperatures;
        for(int i=0; i<MaxTempIterations; i++) Temperatures.push_back(t_lo+(double(i))*t_Delta);
        //cout << Temperatures << endl; getchar();
        RealVec AvgHeldOutLkls;
        RealVec AvgHeldOutEmpiricalDeviations;
        cxsc::real CVgain=0.0;
        // first get a root box containing all points from first burst
        AdaptiveHistogram adhA0;//main hist object for first burst
        successfulInsertion = false;
        successfulInsertion = adhA0.insertFromRVec(transformedData);//insert transformed data
        if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
        //transformedData.clear();//keep the transformed data!!
      n = adhA0.getRootCounter();
      d = adhA0.getDimensions ();
//cout << "transformed data:  n = " << n << endl; getchar();

        size_t critSEB;
        size_t maxLeavesSEB;
        size_t maxLeavesCarving;
	unsigned long int seedStarts = seed+TempIterations;
if(TempIterate)
{
        const size_t N = adhA0.getRootCounter();//size of successfully inserted transformed data
        const size_t K = 10;//10; //K-fold CV
        const size_t KofN = N/K; //cout << KofN << endl; getchar();
        size_t nTrain;
        adhA0.clearAllHistData();//clear the transformed data to make space during CV
        RVecData transformedDataT;//container to keep the transformed Training data from first burst
        RVecData transformedDataV;//container to keep the transformed Validation data from first burst
        // set up for permutations
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_permutation * p = gsl_permutation_alloc (N);
        gsl_permutation * q = gsl_permutation_alloc (N);
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        //printf ("initial permutation:");  
        gsl_permutation_init (p);
        //gsl_permutation_fprintf (stdout, p, " %u"); printf ("\n"); getchar();
        do
        {
          LogTemperaturePrior logPrior(Temperatures[TempIterations]);//t_lo);
	  seedStarts += TempIterations;
          real HeldOutLkl = 0.0;
          real HeldOutEmpiricalDeviation = 0.0;
          // a container for our histograms at various temperatures
          AdaptiveHistogram adhA0cv(adhA0.getRootBox()); // make adh for CV with root box from adhA0
          for (int cvI=1; cvI<K; cvI++)//K-fold CV loop
          { 
            gsl_ran_shuffle (r, p->data, N, sizeof(size_t));
            successfulInsertion = false;
            std::vector< subpavings::AdaptiveHistogram* > histsT;
            adhA0cv.clearAllHistData();//clear the cv histogram before insertion
            adhA0cv.mergeUp();//Merge the possibly multileaf cv histogram up to just root. 
            transformedDataV.clear(); transformedDataT.clear();//clear Cv containers
            for(size_t i=0; i<KofN; i++) transformedDataV.push_back(transformedData[gsl_permutation_get(p,i)]);
            for(size_t i=KofN; i<N; i++) transformedDataT.push_back(transformedData[gsl_permutation_get(p,i)]);
            successfulInsertion = adhA0cv.insertFromRVec(transformedDataT);//insert transformed data
            if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
	    /* some guesses for max points in a node to stop posterior queue */
            nTrain = adhA0cv.getRootCounter(); //cout << "nTrain = " << nTrain << endl; getchar();
	    critSEB = static_cast<size_t>(std::log(static_cast<double>(nTrain)));//can be as low as 1
	    /* some guesses for maximum leaves we'll let SEB queue go to */
	    maxLeavesSEB = nTrain;// / critSEB; // integer division
	    maxLeavesCarving = maxLeavesSEB / 2; // integer division
	    SPSNodeMeasureVolMassMinus compCarving(nTrain);
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving(compCarving, maxLeavesCarving);
	    SPSNodeMeasureCount compSEB;
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB(compSEB, critSEB, maxLeavesSEB);
	    CarverSEB::findStartingPointsBest(adhA0cv, histsT, evaluatorCarving, evaluatorSEB, logPrior, 
						minPoints, chooseStarts, keep, stopOnMaxPosterior, 
						postFileName, checkPostFileNameBase, precPQ, seedStarts);
            PiecewiseConstantFunction pcfT(*histsT[0]);
            pcfT.smearZeroValues(0.0000001);
	    //assert (pcfT.getTotalIntegral() == cxsc::real(1.0));
            histsT[0]->clearAllHistData();//clear the data from training big burst
            histsT[0]->insertFromRVec(transformedDataV);//insert validation data
            HeldOutLkl += pcfT.getLogLikelihood(*histsT[0]);
            PiecewiseConstantFunction pcfV(*histsT[0]);
            //histsT[0]->clearAllHistData();//clear the data from validation big burst
            //histsT[0]->mergeUp();//merge up to root
            HeldOutEmpiricalDeviation += pcfT.getL1Distance(pcfV);
            //cout << HeldOutLkl << '\t' << HeldOutEmpiricalDeviation << endl; getchar();
            //to free all the contents of histsT at the end
            for (size_t i = 0; i < histsT.size(); ++i)
            {
              if (NULL != histsT[i]) delete histsT[i];
              histsT[i] = NULL;
            }
          }
          //cout << "Avg HeldOutLkl = " << HeldOutLkl/double(K) << '\n' 
          //     << "Avg HeldOutEmpiricalDeviation = " << HeldOutEmpiricalDeviation/double(K) << endl; getchar();
          TempIterations++;
          AvgHeldOutLkls.push_back(HeldOutLkl/double(K));
          AvgHeldOutEmpiricalDeviations.push_back(HeldOutEmpiricalDeviation/double(K));
        }
        while (TempIterations<MaxTempIterations);// && CVgain>0.1);
          cout << Temperatures << endl;
          cout << AvgHeldOutLkls << endl << AvgHeldOutEmpiricalDeviations << endl;
//          cout << "Temp Iteration Number " << TempIterations << endl; getchar();

        gsl_permutation_free (p);
        gsl_rng_free (r);
}

        double t_opt = 0.3;//0.3;//1.0;
if(TempIterate)
{
        const double lowest_double = -std::numeric_limits<double>::max();
        double AvgHeldOutLkls_opt = lowest_double;//-numeric_limits<double>::max();//-100000000000000000;
        cout << "lowest double = " << lowest_double << endl;
        for(int i=0; i < Temperatures.size(); i++)
        {
          if(AvgHeldOutLkls[i]>AvgHeldOutLkls_opt)
          {
            AvgHeldOutLkls_opt = _double(AvgHeldOutLkls[i]);
            t_opt=_double(Temperatures[i]);
          }
        }
//        cout << "optimal temperature and AvgHeldOutLkls are : " << t_opt << '\t' << AvgHeldOutLkls_opt << endl; getchar();
}
	    /* some guesses for max points in a node to stop posterior queue */
            if(adhA0.getRootCounter()==0)
            {  
              successfulInsertion = false;
              successfulInsertion = adhA0.insertFromRVec(transformedData);//insert all data
              if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
            }
	    critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));//can be as low as 1
	    /* some guesses for maximum leaves we'll let SEB queue go to */
	    maxLeavesSEB = n;//*adhA0.getDimensions();// / critSEB; // integer division
	    maxLeavesCarving = maxLeavesSEB;//*adhA0.getDimensions();///2; // integer division
	    SPSNodeMeasureVolMassMinus compCarving(n);
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving( compCarving, maxLeavesCarving);
	    SPSNodeMeasureCount compSEB;
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
            LogTemperaturePrior logPrior(t_opt);
            if(CarvingMaxPosterior)
            { 
				//find out more about this function
	    CarverSEB::findStartingPointsBest(adhA0, hists, evaluatorCarving, evaluatorSEB, logPrior, 
						minPoints, chooseStarts, keep, stopOnMaxPosterior, 
						postFileName, checkPostFileNameBase, precPQ, seedStarts);
            }
            else
            {
              // function object to compare nodes on count
              // ie split node with largest count first
              CompCount compCount;
              CompVol compVol;
              /* minimum points to use when splitting.
              A node will not be splittable if either child would then have
               < minPoints of data associated with it. */
              size_t minChildPoints = 10;//static_cast<size_t>(std::log(static_cast<double>(n)))
              //double minVolume = 0.01;
              //double minVolume = 0.001;
              double minVolume = 0.000001;
              cout << "minChildPoints, minVolume = " << minChildPoints << "\t" << minVolume << endl; getchar();
              adhA0.prioritySplit(compCount, maxLeavesSEB, NOLOG, minChildPoints, minVolume);
              //adhA0.prioritySplit(compVol, maxLeavesSEB, NOLOG, minVolume);
              hists.push_back(& adhA0);
              //AdaptiveHistogram* adhPtr = & adhA0;
              //hists.push_back(adhPtr);
            }
	clock_t endtime = clock();	
	double timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
	cout << "time to get starts = " << timingStarts << endl;

	assert(hists.size()==1);
        hists[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt", 6,true);
        // create a name for the file to output
        // To realize a file output
        //make a piecewise constant function of the best adaptive histogram density estimate
        //PiecewiseConstantFunction pcf0(*hists[0]);
        PiecewiseConstantFunction* pcfPtr = new PiecewiseConstantFunction(*hists[0]);
        //AdaptiveHistogram adh0(*hists[0]);//make adaptive histogram
        pcfs.push_back(pcfPtr);
        pcfs[0]->smearZeroValues(0.0000001);
	cxsc::real integral0 = pcfs[0]->getTotalIntegral();
        pcfs[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt"+"pcf", 6,true);
	cout << "pcfs[0]->getTotalIntegral() = " << integral0 << endl;
	assert (integral0 == cxsc::real(1.0));
        //make a copy of hists[0] if you want to keep it 
        //AdaptiveHistogram adhAburst(hists[0]);
        //hists[0]->clearAllHistData();//clear the data from first big burst
        //cout << hists.size() << endl; getchar();
      //-----------------------------------------------------------------------
      // end of getting best histogram for first big burst
      //-----------------------------------------------------------------------

        //read in the smaller bursts
        std::string burstFileName;
        for (int burstI=1; burstI<=B; burstI++)
        {
          size_t headerlines = 0;
          std::string burstHistFileName;
          {
            ostringstream oss;
            oss << burstsFileBaseName << burstI << ".txt";
            burstFileName = oss.str();
            ostringstream oss1;
            oss1 << burstsFileBaseName << "hist_" << burstI << ".txt";
            burstHistFileName = oss1.str();
          }
          // get a count of the lines in the txt file
          int dataCount = countLinesInTxt(burstFileName);
          // tell user how many lines there are in the file
          cout << "The file " << burstFileName << " has " << dataCount
              << " lines in it" << endl;
          RVecData burstData;
          RVecData TransformedBurstData;
          readRvectorsFromTxtOrd (burstData, burstFileName, headerlines);
          cout << "size of burstData = " << burstData.size() << endl;// getchar(); 
          assert(burstData.size()>0); assert(burstData[0]==d);
          //if(burstI==5) cout << endl << endl ; getchar();
          //transforming burstData in place 
          int counter=0;
          for (RVecDataItr it = burstData.begin(); it != burstData.end(); ++it) 
          {
            rvector rawdata = *it;        
            //counter++; getchar(); cout << counter << endl << rawdata << endl;
            // transformation by centered non-rotating rescaling        
            for (int i = 1; i <= d; ++i) {
                //cout << rawdata[i] << '\t' << mean[i] << '\t' << rawdata[i]-mean[i] << endl; getchar();
                rawdata[i] -= mean[i];
                rawdata[i] /= sqrt(V[i][i]);
            }// end i-loop
            //cout << rawdata << endl; getchar();
            TransformedBurstData.push_back(rawdata);
          } // end iteration through data
          //if(burstI==5) cout << endl << endl; getchar();

          hists[0]->clearAllHistData();//clear the data from previous burst
          hists[0]->outputToTxtTabs(burstHistFileName);
          ivector Rbox = hists[0]->getRootBox();
          cout << "root box is" << endl << Rbox << endl; //getchar(); 
          successfulInsertion = false;
          //successfulInsertion = hists[0]->insertFromRVec(burstData);//insert transformed data
          successfulInsertion = hists[0]->insertFromRVec(TransformedBurstData);//insert transformed data
          if (!successfulInsertion) throw std::runtime_error("Failed to insert burst data");



          hists[0]->outputToTxtTabs(burstHistFileName+"post");
          size_t m = hists[0]->getRootCounter();
          //Likelihood (probability under IID) of current small burst given density from first big burst 
          real lkl = pcfs[0]->getLogLikelihood(*hists[0]);
          cout << "size of small burst inserted and its likelihood = "<< m << '\t' << lkl << endl << endl;
          //if(burstI==5) cout << endl << endl; getchar();
          LklBursts.push_back(lkl);
          subpavings::PiecewiseConstantFunction* pcfPtr = new PiecewiseConstantFunction(*hists[0]);
          //subpavings::AdaptiveHistogram* adhPtr = new AdaptiveHistogram(*hists[0]);
          pcfs.push_back(pcfPtr);
          cxsc::real integralI = pcfPtr->getTotalIntegral();
        //pcf0.outputToTxtTabs(s, prec, true);
        cout << "pcfPtr->getTotalIntegral() = " << integralI << endl; //getchar();
        cout << "pcfPtr->getIAE(pcf0)          = " << pcfPtr->getIAE(*pcfs[0]) << endl; //getchar();
        //cout << "adh0.getL1Distance(*hists[0]) = " << adh0.getL1Distance(*hists[0]) << endl; //getchar();
        //cout << "adh0.getL1Distance(*adhPtr)   = " << adh0.getL1Distance(*adhPtr) << endl; //getchar();

	pcfPtr->outputToTxtTabs(burstHistFileName+"pcf", 6,true);
        }
/*
        for (size_t i = 0; i < pcfs.size(); ++i) 
        {
          cout << "pcfs[" << i << "]->getIAE(*pcfs[0])        " << pcfs[i]->getIAE(*pcfs[0]) << endl;
          cout << "pcfs[" << i << "]->getL1Distance(*pcfs[0]) " << pcfs[i]->getL1Distance(*pcfs[0]) << endl;
          if(i<pcfs.size()-1) cout << "pcfs[i]->getL1Distance(pcfs[i+1]) = " << pcfs[i]->getL1Distance(*pcfs[i+1]) << endl;
        }
*/      
        assert(B==pcfs.size());// check number of bursts
        for (size_t i = 0; i <= B; ++i) 
        {
          if(i!=B) consecL1D[i+1] = pcfs[i]->getL1Distance(*pcfs[i+1]);
          for (size_t j = 0; j <= i; ++j)
          {   
            allPrsL1D[i][j] = pcfs[i]->getL1Distance(*pcfs[j]);
            allPrsL1D[j][i] = allPrsL1D[i][j];
          } 
            //if(i<pcfs.size()-1) cout << "pcfs[i]->getL1Distance(pcfs[i+1]) = " << pcfs[i]->getL1Distance(*pcfs[i+1]) << endl;
        }
     }
     int prec = 6;
     cout << "\n likelihood of small bursts  (using cxsc formatting)" << endl;
     ostream_iterator<cxsc::real> out_it (cout,"\n");
     cout << cxsc::SetPrecision(prec+2, prec);
     copy ( LklBursts.begin(), LklBursts.end(), out_it );
     string PairwiseL1FileName;
     PairwiseL1FileName = burstsFileBaseName+"L1.txt";//"dataFiles/B2K0.txt";
     ofstream L1file;
     //L1file.open(PairwiseL1FileName.c_str(),ios::out | ios::app);
     L1file.open(PairwiseL1FileName.c_str(),ios::out);
     cout << "consecutive L1 distances " << endl;
     cout << consecL1D << endl;
     cout << "pairwise L1 distances " << endl;
     cout << allPrsL1D << endl << endl; 
     cout << allPrsL1D(1,B,1,B) << endl; 
     L1file << allPrsL1D(1,B,1,B) << endl; 
     L1file.close();


    //to free all the contents of pcfs at the end
    for (size_t i = 0; i < pcfs.size(); ++i) 
    {
        if (NULL != pcfs[i]) delete pcfs[i];
        pcfs[i] = NULL;
    }
if(CarvingMaxPosterior)
{ 
    //to free all the contents of hists at the end
    for (size_t i = 0; i < hists.size(); ++i) 
    {
        if (NULL != hists[i]) delete hists[i];
        hists[i] = NULL;
    }
}    
    return 0;

} // end of bivariate gaussian test program

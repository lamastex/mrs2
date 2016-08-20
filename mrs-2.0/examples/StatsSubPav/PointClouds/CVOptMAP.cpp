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


/*! \file CVOptMAP.cpp
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


// return a vector of the top k indices of a
void topk(RealVec a, vector<size_t> & indtop){
  multimap<real, size_t> m; // mapping from value to its index
  RealVecItr it;
  for (it = a.begin(); it != a.end(); ++it)
    m.insert(make_pair(*it, it - a.begin()));

  multimap<real, size_t>::reverse_iterator itm; // mapping from value to its index
  size_t indx=0;
  for (itm = m.rbegin(); itm != m.rend(); ++itm){
    //cout << itm->first <<" , "<< itm->second << endl;
    indtop[indx]=itm->second;
    indx++;
  } 
}

bool translateByMeanScaleByVar(RVecData & RawData);

bool selectPriorByLlkCV (RVecData & transformedData, size_t K, double t_lo, double t_hi,
				int LocalMaxTempIterations, int MaxTempIterations,
				RealVec & AvgHeldOutLkls, vector<double> & Temperatures, 
				double & t_opt, double & AvgHeldOutLkls_opt,
				size_t minPoints, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, int precPQ, 
				unsigned long int seedStarts);

bool selectPriorByLv1OutCV (RVecData & Data, double t_lo, double t_hi,
				int LocalMaxTempIterations, int MaxTempIterations,
				RealVec & Lv1OutCVScores, vector<double> & Temperatures, 
				double & t_opt, double & Lv1OutCVScores_opt,
				size_t minPoints, double minVolume, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, string burstsFileBaseName, 
				bool printHist, int precPQ, 
				bool CarvingMaxPosterior, 
				unsigned long int seedStarts);

bool optPQMCAdapHist (RVecData & transformedData, 
				double & t_opt, 
				vector< subpavings::AdaptiveHistogram* > & hists, 
				vector< subpavings::PiecewiseConstantFunction* > & pcfs,
				size_t minPoints, double minVolume, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, string burstsFileBaseName, 
				bool printHist, int precPQ, bool CarvingMaxPosterior, 
				unsigned long int seedStarts);

int main(int argc, char ** argv) 
{
    // ------- prepare to generate some data for the tests -----------
    //const int n=10000;  // number to generate

    clock_t start, end;     // for timing
    double timeTaken;

    RVecData Data;//container to keep the data from first burst or later bursts
    //default values
    bool simulateOrUseDataFile=false;
    bool transformData=false;
    bool priorSelectByCV=false;
    double t_opt = 0.5;//
    size_t dim = 1; //dimensions
    size_t n = 1000; //sample size
    int reps = 1; //replications - just used to change seed in each replicate simulation	
    size_t K = 5; //K-fold CV
    double t_lo=0.01;//lowest temperature in search 
	// make sure this is high enough to avoid
	//terminate called after throwing an instance of 'std::runtime_error'
  	//what():  findStartingPointsMaxLeafSpread(...): Not enough achievable starts

    double t_hi=5.0;//10.0;//highest temperature in search
    int LocalMaxTempIterations=10;
    int MaxTempIterations=5;//AAAAA
        int precPQ = 5;
	bool stopOnMaxPosterior = true;//bool for carver PQ
        // set to do Carver PQ + SEBPB for posterior maximization
        bool CarvingMaxPosterior=true;// false means using SEB PQ with minChildPoints
	long unsigned int seed = 7871234;
	int keep = 1;
	int chooseStarts = 10;//10;
	size_t minPoints = 1;
        double minVolume = 0.001;
/*
	int MaxTempIterations= atoi(argv[6]); //how many temperatures
	double t_lo=atof(argv[7]);//lowest temperature in search
	double t_hi=atof(argv[8]); //highest temperature in search
*/
    // create a name for the file to use
    string samplesFileName; // for samples
    string pqFileNameBase("dataCVOptMAP/PQ/pcAdh_");
    string burstsFileBaseName;
    burstsFileBaseName = "dataCVOptMAP/datasets/rosenbrock_d_1_n_E4";
    if (argc > 1) {
      burstsFileBaseName = argv[1];
    }
    simulateOrUseDataFile = (burstsFileBaseName == "s");
    if (argc > 2) {
       transformData = atoi(argv[2]) != 0;
    }
    if (argc > 3) {
       t_opt = atof(argv[3]);
       //cout << "t_opt = " << t_opt << endl; getchar();
       priorSelectByCV = (t_opt == 0.0);
    }
    if (argc > 4) {
       t_hi = atof(argv[4]);
    }
    if (argc > 5) {
       MaxTempIterations = atoi(argv[5]);
    }
    if (argc > 6) {
       n = atoi(argv[6]);
    }
    if (argc > 7) {
       dim = atoi(argv[7]);
    }
    if (argc > 8) {
       reps = atoi(argv[8]);// just using a seed for 1 replicate
    }
    if (argc > 9) {
       minVolume = atof(argv[9]);
    }
    if (argc > 10) {
       minPoints = atoi(argv[10]);
    }
    if (argc > 11) {
       chooseStarts = atoi(argv[11]);
    }
    // Mixture of Normals is made in either case
    long unsigned int MVNseed = 9876;// seed for data simulator
    MixtureMVN* mixMVNptr = makeMixture(dim , MVNseed); 
    if (simulateOrUseDataFile){
        stringstream ss; ss << "dataCVOptMAP/datasets/sim_" << dim << "_" << n;
        burstsFileBaseName = ss.str();
	//cout << "simulating data" << endl; getchar();
	long unsigned int dataseed = MVNseed+reps;// modified by reps seed
	mixMVNptr->resetPRNG(dataseed);
	cout << "\nGenerate " << n << " random values:" << endl;
	mixMVNptr->prn(Data, n);
        //cout << "size of Data = " << Data.size() << endl; getchar(); 
	//cout << "simulating data setup done" << endl; getchar();
    }
    else {
    	samplesFileName = burstsFileBaseName;//"dataFiles/B2K0.txt";
        // get a count of the lines in the txt file
        int dataCount = countLinesInTxt(samplesFileName);
        size_t headerlines = 0;
        cout << "The file " << samplesFileName << " has " << dataCount-1
             << " lines in it" << endl << endl;
        // put in the data in a 'pulse' with no splitting, ie one big box
        readRvectorsFromTxtOrd (Data, samplesFileName, headerlines);
        //bool successfulInsertion = adhA0.insertFromRVec(simdata);
        //bool successfulInsertion = adhA0Raw.insertRvectorsFromTxt(samplesFileName);
        cout << "size of Data = " << Data.size() << endl; getchar(); 
     }// end of getting data from file

     if(transformData){
       //bool successfullyTransformed=false;
       cout << "transforming data OK - press ENTER" << endl; getchar();
       bool successfullyTransformed=translateByMeanScaleByVar(Data);
       if (!successfullyTransformed) throw std::runtime_error("Failed to transform data");
     }
      string postFileName;
      postFileName = pqFileNameBase+"_Posterior";
      string checkPostFileNameBase;
      checkPostFileNameBase = pqFileNameBase+"_CheckPosterior";

      clock_t starttime = clock();
      // n=3000,d=2: 1/100 = undersmoother, 10 is oversmoothed/raggedy, (1/10,1.0) is ok
      std::vector<double> Temperatures;
      //cout << Temperatures << endl; getchar();
      RealVec CVScores;
      unsigned long int seedStarts = seed+12345;
        
      bool CVSuccessful=false;
      double CVScores_opt = 0.0;
        
      if(priorSelectByCV){
          cout << " doing automatic prior selection by CV " << endl;
	// this is a slower and generic K-fold CV method 
	// (held out Lkl is used as an example and should be replaced 
	//with the appropriate scoring rule for maximization over parameter in [t_lo, t_hi])
/*
	CVSuccessful = selectPriorByLlkCV (Data, K, t_lo, t_hi,
				LocalMaxTempIterations, MaxTempIterations,
				CVScores, Temperatures, 
				t_opt, CVScores_opt,
				minPoints, chooseStarts, keep, 
				stopOnMaxPosterior, postFileName, 
				checkPostFileNameBase, precPQ, seedStarts);
*/
	bool printCVHists=false;
	CVSuccessful = selectPriorByLv1OutCV (Data, t_lo, t_hi,
				LocalMaxTempIterations, MaxTempIterations,
				CVScores, Temperatures, 
				t_opt, CVScores_opt,
				minPoints, minVolume, chooseStarts, keep, 
				stopOnMaxPosterior, postFileName, 
				checkPostFileNameBase, burstsFileBaseName, 
				printCVHists, precPQ, 
				CarvingMaxPosterior, 
				seedStarts);
      }
      //cout << "optimal temperature and AvgHeldOutLkls are : " << t_opt << '\t' << AvgHeldOutLkls_opt << endl; getchar();

// make the below modular also!!!
    // a container for our histograms
    std::vector< subpavings::AdaptiveHistogram* > hists;
    // a container for our PCFs of histograms
    std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
    //getchar();
bool succPQMCopt = false;
bool printHist=true;
succPQMCopt = optPQMCAdapHist (Data, 
				t_opt, hists, pcfs,
				minPoints, minVolume, chooseStarts, keep, 
				stopOnMaxPosterior, postFileName, 
				checkPostFileNameBase, burstsFileBaseName, printHist, precPQ, 
				CarvingMaxPosterior, seedStarts);

	clock_t endtime = clock();	
	double timingStarts = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
	cout << "time to get prior-selected adaptive hist = " << timingStarts << endl;

    // L1 error calculations
    if (simulateOrUseDataFile){
    real lv1outCVScore = pcfs[0]->getLeave1OutCVScore(*hists[0]);
    //lv1outCVScore = tmpPcf.getTotalIntegral() + lv1outCVScore;
    cout << "leave-1-out CV summand = " << lv1outCVScore << endl;      
    /* Get L1-distance and KL-distance */
	cout << "Getting the L1-distance and KL-distance: " << endl;
	
	// get quasi random points in the box 
	ivector box = pcfs[0]->getRootBox();
    	std::vector < std::vector < real > > qrPts;
    	int intN = 10000000;
	getQuasiRandomPoints( box, qrPts, intN);
				
	std::vector < real > estDensities_QR;
	PiecewiseConstantFunction pcfSmeared 
					= pcfs[0]->makeSmearZeroValues(1/(1000000.0));
	getPCFDensities(pcfSmeared, qrPts, estDensities_QR);
	//cout << "qrPts.size() = " << qrPts.size() << endl 
	//	<< "estDensities_QR.size() = " << estDensities_QR.size() << endl; getchar();
	/*get true densities at the qr points */
	std::vector < real > trueIntPtDensities_QR;
	getTrueDensities(*mixMVNptr, qrPts, trueIntPtDensities_QR);
	/*approx L1 errors*/
	real boxVol = realVolume(box);
	real estL1_QR = boxVol * avAbsDiffDen(trueIntPtDensities_QR, estDensities_QR);
        //cout << "apprx L1 error = " << estL1_QR << endl; //getchar();
        cout << "optimal temperature, CVScores_opt, lv1outCV, apprx L1 error are : " << t_opt << '\t' << CVScores_opt << '\t' << lv1outCVScore << '\t' << estL1_QR << endl; //getchar();
        cout << "minPoints, minVolume, chooseStarts, keep are : " << minPoints << '\t' << minVolume << '\t' << chooseStarts << '\t' << keep << endl;
    }
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

        //delete the data generator 
        delete mixMVNptr;	

    return 0;

} // end of program

//-----------------------------------------------------------------------
// optimizing over smoothing parameter Temp - 
// getting best histogram over parameter in [t_lo, t_hi] by max/minmising 
// a CV based score
//-----------------------------------------------------------------------
bool selectPriorByLv1OutCV (RVecData & Data, double t_lo, double t_hi,
				int LocalMaxTempIterations, int MaxTempIterations,
				RealVec & Lv1OutCVScores, vector<double> & Temperatures, 
				double & t_opt, double & Lv1OutCVScores_opt,
				size_t minPoints, double minVolume, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, string burstsFileBaseName, 
				bool printHist, int precPQ, 
				bool CarvingMaxPosterior, 
				unsigned long int seedStarts)
{
  bool success=false;
  int TempIterations=0;
  std::vector<double> LocalTemperatures;
     do{// outer temp iterations
          if (TempIterations!=0){
            std::vector<size_t> indtop(Lv1OutCVScores.size());
            topk(Lv1OutCVScores, indtop);
            double tBest =Temperatures[indtop[0]];
            double t2ndBest =Temperatures[indtop[1]];
            double t3rdBest =Temperatures[indtop[2]];
            double tWorst = max(abs(tBest-t2ndBest),abs(tBest-t3rdBest));
            t_lo = max(t_lo,tBest-tWorst); t_hi = tBest+tWorst;
//            cout << t_lo << " , " << t_hi << " : "<< tBest << " , " << t2ndBest << " , " << t3rdBest << endl; getchar();
            if (t_hi-t_lo<0.001){ 
		//|| abs(Lv1OutCVScores[indtop[0]]-Lv1OutCVScores[indtop[1]])<0.00001) {
              	cout << "reaching temp values < 0.001 or Lv1OutCVScores diff < 0.00001"<<endl; 
		getchar(); 
		break;
            }
          }
          double t_Delta=(t_hi-t_lo)/double(LocalMaxTempIterations-1);
          LocalTemperatures.clear();
          for(int i=0; i<LocalMaxTempIterations; i++){ 
                        Temperatures.push_back(t_lo+(double(i))*t_Delta);
                        LocalTemperatures.push_back(t_lo+(double(i))*t_Delta);
          }
        int LocalTempIterations=0;
        do
        {
          //LogTemperaturePrior logPrior(Temperatures[TempIterations]);//t_lo);
          
          //LogTemperaturePrior logPrior(LocalTemperatures[LocalTempIterations]);//t_lo);
          double temperatureNow = LocalTemperatures[LocalTempIterations];
	  seedStarts += TempIterations;
	  // a container for our histograms
  	  std::vector< subpavings::AdaptiveHistogram* > hists;
  	  // a container for our PCFs of histograms
  	  std::vector< subpavings::PiecewiseConstantFunction* > pcfs;
  	  bool succPQMCopt = false;

	  //ostringstream strs; strs << temperatureNow;
          //string burstsFileBaseNameNow = burstsFileBaseName+"_"+strs.str();
	  succPQMCopt = optPQMCAdapHist (Data, temperatureNow, hists, pcfs,
				minPoints, minVolume, chooseStarts, keep, 
				stopOnMaxPosterior, postFileName, 
				checkPostFileNameBase, burstsFileBaseName, printHist, precPQ, 
				CarvingMaxPosterior, seedStarts);

    	  real lv1outCVScore = pcfs[0]->getLeave1OutCVScore(*hists[0]);
          Lv1OutCVScores.push_back(-1.0*lv1outCVScore);// -1.0* so we want to max to be min
          LocalTempIterations++;
    	  //to free all the contents of pcfs at the end
    	  for (size_t i = 0; i < pcfs.size(); ++i) 
    	    { if (NULL != pcfs[i]) delete pcfs[i]; pcfs[i] = NULL;}
      	  //to free all the contents of hists at the end
    	  if(CarvingMaxPosterior)
    	  { for (size_t i = 0; i < hists.size(); ++i) 
              { if (NULL != hists[i]) delete hists[i]; hists[i] = NULL;}
          }    
        }
        while (LocalTempIterations<LocalMaxTempIterations);// && CVgain>0.1);
      TempIterations++;
      }
      while (TempIterations<MaxTempIterations);// && CVgain>0.1);
////////////////////////////////////////////////////////////////////////////////////////////////
          cout << "Temperatures:" << endl;
          cout << Temperatures << endl;
          cout << "- Lv1OutCVScores: (looking for maximum of -1.0*Lv1OutCVScores)" << endl;
          cout << Lv1OutCVScores << endl << endl;
          cout << "Temp Iteration Number " << TempIterations << endl; //getchar();
          cout << "MaxTempIterations = " << MaxTempIterations << endl; //getchar();

        vector<size_t> indtop(Lv1OutCVScores.size());
        topk(Lv1OutCVScores, indtop);
        t_opt = _double(Temperatures[indtop[0]]);
        Lv1OutCVScores_opt = _double(Lv1OutCVScores[indtop[0]]);

        success=true;
        return success;
}

bool optPQMCAdapHist (RVecData & transformedData, 
				double & t_opt, 
				vector< subpavings::AdaptiveHistogram* > & hists, 
				vector< subpavings::PiecewiseConstantFunction* > & pcfs,
				size_t minPoints, double minVolume, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, string burstsFileBaseName, 
				bool printHist, int precPQ, bool CarvingMaxPosterior, 
				unsigned long int seedStarts)
{
        bool succPQMCopt = false;
        AdaptiveHistogram adhA0;//main hist object for first burst
	    /* some guesses for max points in a node to stop posterior queue */
        //    if(adhA0.getRootCounter()==0)
          //  {  
              bool successfulInsertion = false;
              successfulInsertion = adhA0.insertFromRVec(transformedData);//insert all data
              if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
              size_t n = adhA0.getRootCounter();
              size_t d = adhA0.getDimensions ();
            //}
	    size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));//can be as low as 1
	    /* some guesses for maximum leaves we'll let SEB queue go to */
	    size_t maxLeavesSEB = n;//*adhA0.getDimensions();// / critSEB; // integer division
	    size_t maxLeavesCarving = maxLeavesSEB/2;//*adhA0.getDimensions();///2; // integer division
	    SPSNodeMeasureVolMassMinus compCarving(n);
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving( compCarving, maxLeavesCarving);
	    SPSNodeMeasureCount compSEB;
	    AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB( compSEB, critSEB, maxLeavesSEB);

    // set up prior distribution object
    //LogCatalanPrior logPrior;
    // n=3000,d=2: 1/100 = undersmoother, 10 is oversmoothed/raggedy, (1/10,1.0) is ok
    LogTemperaturePrior logPrior(t_opt);

    if(CarvingMaxPosterior)
    { 
      //find out more about this function XXX minVolume = 0.01 - make smaller!!! PASS from outdide!!! TODO
	CarverSEB::findStartingPointsBest(adhA0, hists, evaluatorCarving, evaluatorSEB, 
						logPrior, 
				//minPoints, 0.1758, chooseStarts, keep, stopOnMaxPosterior, 
				minPoints, minVolume, chooseStarts, keep, stopOnMaxPosterior, 
				//minPoints, 0.001, chooseStarts, keep, stopOnMaxPosterior, 
				//minPoints, chooseStarts, keep, stopOnMaxPosterior, 
				"", "", precPQ, seedStarts);
				//postFileName, checkPostFileNameBase, precPQ, seedStarts);
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
//THINK!!! before setting these!!!
              size_t minChildPoints = 80;//static_cast<size_t>(std::log(static_cast<double>(n)))
              double minVol = 0.01;
              adhA0.prioritySplit(compCount, maxLeavesSEB, NOLOG, minChildPoints, minVol);
              //adhA0.prioritySplit(compVol, maxLeavesSEB, NOLOG, minVolume);
              hists.push_back(& adhA0);
              //AdaptiveHistogram* adhPtr = & adhA0;
              //hists.push_back(adhPtr);
    }
	//assert(hists.size()==1);

        //hists[0]->outputToTxtTabs(burstsFileBaseName+"hist_"+"0"+".txt", 6,true);
        // create a name for the file to output
        // To realize a file output
        //make a piecewise constant function of the best adaptive histogram density estimate
        //PiecewiseConstantFunction pcf0(*hists[0]);
        PiecewiseConstantFunction* pcfPtr = new PiecewiseConstantFunction(*hists[0]);
        //AdaptiveHistogram adh0(*hists[0]);//make adaptive histogram
        pcfs.push_back(pcfPtr);
        pcfs[0]->smearZeroValues(0.0000001);
	cxsc::real integral0 = pcfs[0]->getTotalIntegral();
        if(printHist) pcfs[0]->outputToTxtTabs(burstsFileBaseName+"_pcf.txt", 6,true);
	//cout << "pcfs[0]->getTotalIntegral() = " << integral0 << endl;
	//assert (integral0 == cxsc::real(1.0));
        //make a copy of hists[0] if you want to keep it 
        //AdaptiveHistogram adhAburst(hists[0]);
        //hists[0]->clearAllHistData();//clear the data from first big burst
        //cout << hists.size() << endl; getchar();
      //-----------------------------------------------------------------------
      // end of getting best histogram for first big burst
      //-----------------------------------------------------------------------

    succPQMCopt = true;
    return succPQMCopt;
}
//-----------------------------------------------------------------------
// optimizing over smoothing parameter Temp - 
// getting best histogram over parameter in [t_lo, t_hi] by max/minmising 
// a CV based score
//-----------------------------------------------------------------------
// this is a slower and generic K-fold CV method -- USE Leave1Out for Histograms!!!
// CAUTON: held out Lkl is used as an example and should be replaced 
//with the appropriate scoring rule for maximization over parameter in [t_lo, t_hi]
bool selectPriorByLlkCV (RVecData & transformedData, size_t K, double t_lo, double t_hi,
				int LocalMaxTempIterations, int MaxTempIterations,
				RealVec & AvgHeldOutLkls, vector<double> & Temperatures, 
				double & t_opt, double & AvgHeldOutLkls_opt,
				size_t minPoints, int chooseStarts, int keep, 
				bool stopOnMaxPosterior, string postFileName, 
				string checkPostFileNameBase, int precPQ, 
				unsigned long int seedStarts)
{
  RealVec AvgHeldOutEmpiricalDeviations;
  // first get a root box containing all points 
  AdaptiveHistogram adhA0;//main hist object 
  bool successfulInsertion = false;
  successfulInsertion = adhA0.insertFromRVec(transformedData);//insert transformed data
  if (!successfulInsertion) throw std::runtime_error("Failed to insert transformed data");
  //transformedData.clear();//keep the transformed data!!
  size_t n = adhA0.getRootCounter();
  size_t d = adhA0.getDimensions ();
  //cout << "transformed data:  n = " << n << endl; getchar();

  bool success=false;
  const size_t N = adhA0.getRootCounter();//size of successfully inserted transformed data
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
////////////////////////////////////////////////////////////////////////////////////////////////
        int TempIterations=0;
        std::vector<double> LocalTemperatures;
     do{// outer temp iterations
          if (TempIterations!=0){
            std::vector<size_t> indtop(AvgHeldOutLkls.size());
            topk(AvgHeldOutLkls, indtop);
            double tBest =Temperatures[indtop[0]];
            double t2ndBest =Temperatures[indtop[1]];
            double t3rdBest =Temperatures[indtop[2]];
            double tWorst = max(abs(tBest-t2ndBest),abs(tBest-t3rdBest));
            t_lo = max(0.00001,tBest-tWorst); t_hi = tBest+tWorst;
//            cout << t_lo << " , " << t_hi << " : "<< tBest << " , " << t2ndBest << " , " << t3rdBest << endl; getchar();
            if (t_hi-t_lo<0.001 || 
			abs(AvgHeldOutLkls[indtop[0]]-AvgHeldOutLkls[indtop[1]])<0.01) {
              cout << "reaching temp values < 0.001 or likl diff < 1.0"<<endl; getchar(); break;
            }
          }
          double t_Delta=(t_hi-t_lo)/double(LocalMaxTempIterations-1);
          LocalTemperatures.clear();
          for(int i=0; i<LocalMaxTempIterations; i++){ 
                        Temperatures.push_back(t_lo+(double(i))*t_Delta);
                        LocalTemperatures.push_back(t_lo+(double(i))*t_Delta);
          }
        int LocalTempIterations=0;
        do
        {
          //LogTemperaturePrior logPrior(Temperatures[TempIterations]);//t_lo);
          LogTemperaturePrior logPrior(LocalTemperatures[LocalTempIterations]);//t_lo);
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
	    size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(nTrain)));//can be as low as 1
	    /* some guesses for maximum leaves we'll let SEB queue go to */
	    size_t maxLeavesSEB = nTrain/2;// / critSEB; // integer division
	    size_t maxLeavesCarving = maxLeavesSEB / 3; // integer division
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
	    //if(pcfT.getTotalIntegral() != cxsc::real(1.0)) {cout << "433!!!"; getchar();}
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
          AvgHeldOutLkls.push_back(HeldOutLkl/double(K));
          AvgHeldOutEmpiricalDeviations.push_back(HeldOutEmpiricalDeviation/double(K));
          LocalTempIterations++;
        }
        while (LocalTempIterations<LocalMaxTempIterations);// && CVgain>0.1);
      TempIterations++;
      }
      while (TempIterations<MaxTempIterations);// && CVgain>0.1);
////////////////////////////////////////////////////////////////////////////////////////////////
          cout << Temperatures << endl;
          cout << AvgHeldOutLkls << endl << AvgHeldOutEmpiricalDeviations << endl;
          cout << "Temp Iteration Number " << TempIterations << endl; getchar();
          cout << "MaxTempIterations = " << MaxTempIterations << endl; getchar();

        vector<size_t> indtop(AvgHeldOutLkls.size());
        topk(AvgHeldOutLkls, indtop);
        t_opt = _double(Temperatures[indtop[0]]);
        AvgHeldOutLkls_opt = _double(AvgHeldOutLkls[indtop[0]]);

        gsl_permutation_free (p);
        gsl_rng_free (r);

        success=true;
        return success;
}


bool translateByMeanScaleByVar(RVecData & RawData) {
      bool successfullyTransformed=false;
      AdaptiveHistogram adhA0Raw; // let it make its own root box
      bool successfulInsertion = false;
      successfulInsertion = adhA0Raw.insertFromRVec(RawData);//insert all raw data
      RawData.clear();
      adhA0Raw.setHoldAllStats(true);
      if (!successfulInsertion) throw std::runtime_error("Failed to insert data");
      size_t n = adhA0Raw.getRootCounter();
      size_t d = adhA0Raw.getDimensions ();
      //cout << "size of RawData = " << n << endl; getchar();
      //cout << "dimn of RawData = " << d << endl; getchar();

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
      // doing in-place transformation
      //RVecData transformedData;//container to keep the transformed data from first burst
      
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
        //??transformedData.push_back(rawdata);
        RawData.push_back(rawdata);
      } // end iteration through data

      adhA0Raw.clearAllHistData();//clear the raw data
      
      successfullyTransformed=true;
      return successfullyTransformed;
}

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
/*! \file FBinomialPartition.cpp
\brief Trans-dimensional posterior distribution over m binomial partition 
  model example function object class to use with MRSampler class.
*/

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"
#include "imath.hpp"
#include "rmath.hpp"
#include "ivector.hpp"
#include <functional>

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FBinomialPartition.hpp"

/*! constructor for the the toy example of a target which is a binomial partition problem
   on spaces of different dimensionality.
   There are m sets of X_1 := (X_{1,1},X_{1,2},X_{1,n_1}),...,
   X_2 := (X_{2,1},X_{2,2},X_{2,n_2}),..., X_m := (X_{m,1},X_{m,2},X_{m,n_m}) 
   Bernoulli trials that are identical within each set and independent across the sets.     But only the sums of the Bernoulli trials Y_i := \sum_{j=1}^{n_i} X_{i,j}, i=1,...,m 
   are recorded. So effectively, we have m independent binomial trials each of known 
   sizes n_1,n_2,...,n_m with possibly different success probabilities 
   t_i := P(X_{i,j}=1).  We define a binomial partition model as follows.
   Let L={1,2,...,m} be the set of the m Binomial trial labels. 
   For each set partition of L with b blocks, where b=1,2,...,L, define a binomial 
   partition Model such that the success probability is identical between trials whose
   labels fall within a block but may be different between trials with labels in 
   distinct blocks. Thus there are m-th Bell number many partition models with different
   dimensions of parameter spaces.  Our goal is to sample exactly from the 
   transdimensional posterior distribution over the binomial partition models.
*/
FBinomialPartition::FBinomialPartition(vector<int> ns, vector<int> ys, bool LogPi) : n(ns), y(ys)
{
    n_interval_calls = 0;
    n_real_calls = 0;
    setUsingLogDensity (LogPi);
    PriorType = 0; // unif. (normalised Lebesgue measure) prior is default
    SmallestPositive = 1e-100;
    IgnoreBinomCoeffs = true; // ignore the binomial coefficients if set to true
    int M = n[0];//Trials m in y_1,...y_m
    //cout << "M = " << M << endl; getchar();
    //just output the data for now
    cout << "The data for number of successes/heads y_i's are:\n"; 
    for(int i=0;i<=M;i++) cout << y[i] << '\t';cout << "\n"; 
    cout << "The data for number of Bern trials are n_i's are:\n"; 
    for(int i=0;i<=M;i++) cout << n[i] << '\t';cout << "\n";
    
    vector<int> SixZeros(6,0);
    vector<interval> SixIZeros;
    vector<real> SixRZeros;
    for(int i=0;i<6;i++) SixIZeros.push_back(interval(0.0,0.0));
    for(int i=0;i<6;i++) SixRZeros.push_back(real(0.0));
    for(int i=0;i<6;i++) MbyBMatrix.push_back(SixZeros);
    MbyBMatrix[0][0] = 0;//m=0
    MbyBMatrix[1][0] = 1;//m=1, m-Bell=1
    MbyBMatrix[1][1] = 1;//m=1,b=1, m,b-Stirling2nd=1
    MbyBMatrix[2][0] = 2;//m=2, m-Bell=2
    MbyBMatrix[2][1] = 1;//m=2,b=1, m,b-Stirling2nd=1
    MbyBMatrix[2][2] = 1;//m=2,b=2, m,b-Stirling2nd=1
    MbyBMatrix[3][0] = 5;//m=3, m-Bell=5
    MbyBMatrix[3][1] = 1;//m=3,b=1, m,b-Stirling2nd=1
    MbyBMatrix[3][2] = 3;//m=3,b=2, m,b-Stirling2nd=3
    MbyBMatrix[3][3] = 1;//m=3,b=3, m,b-Stirling2nd=1
    MbyBMatrix[4][0] =15;//m=4, m-Bell=15
    MbyBMatrix[4][1] = 1;//m=4,b=1, m,b-Stirling2nd=1
    MbyBMatrix[4][2] = 7;//m=4,b=2, m,b-Stirling2nd=7
    MbyBMatrix[4][3] = 6;//m=4,b=3, m,b-Stirling2nd=6
    MbyBMatrix[4][4] = 1;//m=4,b=4, m,b-Stirling2nd=1
    MbyBMatrix[5][0] =52;//m=5, m-Bell=52
    MbyBMatrix[5][1] = 1;//m=5,b=1, m,b-Stirling2nd=1
    MbyBMatrix[5][2] =15;//m=5,b=2, m,b-Stirling2nd=15
    MbyBMatrix[5][3] =25;//m=5,b=3, m,b-Stirling2nd=25
    MbyBMatrix[5][4] =10;//m=5,b=4, m,b-Stirling2nd=10
    MbyBMatrix[5][5] = 1;//m=5,b=5, m,b-Stirling2nd=1
    for(int i=0;i<31;i++) Blocks.push_back(SixZeros);// Blocks = [[0,0,0,0,0,0],...] thera are 31 for m=5
    // Partns = [[0,0,0,0,0,0],...] there are 76 models m=1,2,3,4,5
    for(int i=0;i<76;i++) Partns.push_back(SixZeros);
    for(int i=0;i<76;i++) Succss.push_back(SixZeros);
    for(int i=0;i<76;i++) Trials.push_back(SixZeros);
    for(int i=0;i<76;i++) Failrs.push_back(SixZeros);
    for(int i=0;i<76;i++) Ibc.push_back(SixIZeros);
    for(int i=0;i<76;i++) Rbc.push_back(SixRZeros);
    //the indices of partition matrix is the model label for the model partition
    Blocks[0][0]=1;Blocks[0][1]=1;//{1} is encoded as [1,1,0,0,0,0], block of size 1 with trial label 1
    Blocks[1][0]=1;Blocks[1][1]=2;//{2} is encoded as [1,2,0,0,0,0], block of size 1 with trial label 2
    Blocks[2][0]=2;Blocks[2][1]=1;Blocks[2][2]=2;//{1,2}
    Blocks[3][0]=1;Blocks[3][1]=3;//{3}
    Blocks[4][0]=2;Blocks[4][1]=1;Blocks[4][2]=3;//{1,3}
    Blocks[5][0]=2;Blocks[5][1]=2;Blocks[5][2]=3;//{2,3}
    Blocks[6][0]=3;Blocks[6][1]=1;Blocks[6][2]=2;Blocks[6][3]=3;//{1,2,3}
    Blocks[7][0]=1;Blocks[7][1]=4;//{4}
    Blocks[8][0]=2;Blocks[8][1]=1;Blocks[8][2]=4;//{1,4}
    Blocks[9][0]=2;Blocks[9][1]=2;Blocks[9][2]=4;//{2,4}
    Blocks[10][0]=2;Blocks[10][1]=3;Blocks[10][2]=4;//{3,4}
    Blocks[11][0]=3;Blocks[11][1]=1;Blocks[11][2]=2;Blocks[11][3]=4;//{1,2,4}
    Blocks[12][0]=3;Blocks[12][1]=1;Blocks[12][2]=3;Blocks[12][3]=4;//{1,3,4}
    Blocks[13][0]=3;Blocks[13][1]=2;Blocks[13][2]=3;Blocks[13][3]=4;//{2,3,4}
    Blocks[14][0]=4;Blocks[14][1]=1;Blocks[14][2]=2;Blocks[14][3]=3;Blocks[14][4]=4;//{1,2,3,4}
     //{1,2} is encoded as [2,1,2,0,0,0], block of size 2 with trial labels 1 and 2
    if(M==1){//just 1 binomial trial y_1
      Partns[0][0]=1;Partns[0][1]=0;// model 0 = {{1}}
    }
    if(M==2){//just 2 binomial trials y_1,y_2
      Partns[0][0]=1;Partns[0][1]=2;//model 0 = {{1,2}}
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=1;//model 1 = {{1},{2}}
    }
    if(M==3){//just 3 binomial trials y_1,y_2,y_3
      Partns[0][0]=1;Partns[0][1]=6;//model 0 = {{1,2,3}}
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=5;//model 1 = {{1},{2,3}}
      Partns[2][0]=2;Partns[2][1]=1;Partns[2][2]=4;//model 2 = {{2},{1,3}}
      Partns[3][0]=2;Partns[3][1]=3;Partns[3][2]=2;//model 3 = {{3},{1,2}}
      Partns[4][0]=3;Partns[4][1]=0;Partns[4][2]=1;Partns[4][3]=3;//model 4 = {{1},{2},{3}}
    }
    if(M==4){//just 4 binomial trials y_1,y_2,y_3,y_4
      //model label 0 = [{1, 2, 3, 4}] 
      // model label 1 = [{1}, {2, 3, 4}], 
      // model label 2 = [{2}, {1, 3, 4}], 
      // model label 3 = [{3}, {1, 2, 4}], 
      // model label 4 = [{4}, {1, 2, 3}], 
      // model label 5 = [{1, 2}, {3, 4}],
      // model label 6 = [{1, 3}, {2, 4}], 
      // model label 7 = [{1, 4}, {2, 3}], 
      // model label 8 = [{1}, {2}, {3, 4}], 
      // model label 9 = [{1}, {3}, {2, 4}],
      // model label 10 = [{1}, {4}, {2, 3}], 
      // model label 11 = [{2}, {3}, {1, 4}], 
      // model label 12 = [{2}, {4}, {1, 3}], 
      // model label 13 = [{3}, {4}, {1, 2}], 
      // model label 14 = [{1}, {2}, {3}, {4}], 
      //model 0 = [{1, 2, 3, 4}] 
      Partns[0][0]=1;Partns[0][1]=14;
      // model 1 = [{1}, {2, 3, 4}], 
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=13;
      // model 2 = [{2}, {1, 3, 4}], 
      Partns[2][0]=2;Partns[2][1]=1;Partns[2][2]=12;
      // model 3 = [{3}, {1, 2, 4}], 
      Partns[3][0]=2;Partns[3][1]=3;Partns[3][2]=11;
      // model 4 = [{4}, {1, 2, 3}], 
      Partns[4][0]=2;Partns[4][1]=7;Partns[4][2]=6;
      // model 5 = [{1, 2}, {3, 4}],
      Partns[5][0]=2;Partns[5][1]=2;Partns[5][2]=10;
      // model 6 = [{1, 3}, {2, 4}], 
      Partns[6][0]=2;Partns[6][1]=4;Partns[6][2]=9;
      // model 7 = [{1, 4}, {2, 3}], 
      Partns[7][0]=2;Partns[7][1]=8;Partns[7][2]=5;
      // model 8 = [{1}, {2}, {3, 4}], 
      Partns[8][0]=3;Partns[8][1]=0;Partns[8][2]=1;Partns[8][3]=10;
      // model 9 = [{1}, {3}, {2, 4}],
      Partns[9][0]=3;Partns[9][1]=0;Partns[9][2]=3;Partns[9][3]=9;
      // model 10 = [{1}, {4}, {2, 3}], 
      Partns[10][0]=3;Partns[10][1]=0;Partns[10][2]=7;Partns[10][3]=5;
      // model 11 = [{2}, {3}, {1, 4}], 
      Partns[11][0]=3;Partns[11][1]=1;Partns[11][2]=3;Partns[11][3]=8;
      // model 12 = [{2}, {4}, {1, 3}], 
      Partns[12][0]=3;Partns[12][1]=1;Partns[12][2]=7;Partns[12][3]=4;
      // model 13 = [{3}, {4}, {1, 2}], 
      Partns[13][0]=3;Partns[13][1]=3;Partns[13][2]=7;Partns[13][3]=2;
      // model 14 = [{1}, {2}, {3}, {4}], 
      Partns[14][0]=4;Partns[14][1]=0;Partns[14][2]=1;Partns[14][3]=3;Partns[14][4]=7;
    }
    int NumModels=MbyBMatrix[M][0];// numer of models based on the number of trials M
    //cout << NumModels << endl; getchar();
    for(int m=0; m < NumModels; m++){
      int ModelDim = Partns[m][0];//number of dimensions (or blocks in partition) of model m
      ivector domain (1, ModelDim); // vector of intervals, with indices 1...2
      LabBox Ldomain;
      for (int i = 1; i <= ModelDim; i++){
        //if working with log ikelihood then domain should be a 
        //proper subset of [0,1] to ensure local Lipschitz condition
        //let it approach [0,1] as needed by your problem
        if(LogPi) {
          //getchar(); RHS is just 1 when SmallestPositive is too small
          domain[i] = interval(SmallestPositive, (1.0-SmallestPositive));
        }
        else domain[i] = interval (0.0, 1.0);
      }
      Ldomain.Box = domain;
      Ldomain.L = m;
      //cout << "Ldomain.L: " << Ldomain.L << endl;
      LabDomainList.push_back (Ldomain);
      //model prior weights are the same
      modelprior.push_back(1.0);
      Succss[m][0]= ModelDim;
      for(int Blki=1; Blki<= Succss[m][0]; Blki++){
        int BlkIndx = Partns[m][Blki];
        int BlkiSz = Blocks[BlkIndx][0];
        Succss[m][Blki]=0;Trials[m][Blki]=0;Failrs[m][Blki]=0;
        interval ibc(1.0,1.0);
        for(int j=1; j<= BlkiSz; j++){
          Trials[m][Blki] += n[ Blocks[BlkIndx][j]];
          Succss[m][Blki] += y[ Blocks[BlkIndx][j]];
          if(!(IgnoreBinomCoeffs)) 
            ibc *= binary_coefficient(n[Blocks[BlkIndx][j]],y[Blocks[BlkIndx][j]]);
        }
        Failrs[m][Blki] = Trials[m][Blki] - Succss[m][Blki];
        Ibc[m][Blki] = ibc;
        Rbc[m][Blki] = (0.5*(Inf(ibc) + Sup(ibc)));
/*
        cout << "m = " << m << " Blki = " << Blki << '\t'
             <<  Succss[m][Blki] << '\t' << Failrs[m][Blki] << '\t' << Trials[m][Blki]  
             << '\t' << Ibc[m][Blki] << '\t' << Rbc[m][Blki] << endl;
        getchar();
*/
      }
    }
    cout << "end of FBinomialPartition constructor. \n";

}

FBinomialPartition::FBinomialPartition(bool LogPi)
{// hard-coded default constructor..
    n_interval_calls = 0;
    n_real_calls = 0;
    setUsingLogDensity (LogPi);
    PriorType = 0; // unif. (normalised Lebesgue measure) prior is default
    SmallestPositive = 1e-100;
    IgnoreBinomCoeffs = true; // ignore the binomial coefficients if set to true
    //TO DO: send M and y_1,...y_m and n_1,...,n_m as input into constructor from command line by modifying
    //targets/FBinomialPartition.*pp and MooreRejSam/BinomialPartition.cpp 
    //TO DO: see previos TO DO, for now we hard-code the m in data y_1,...y_m and assign it to M
    //int M = 4;//Trials m in y_1,...y_m
    //TO DO: see previos TO DO, for now we hard-code the data y_1,...y_m
    //y.push_back(M);//y_0=m
    //TO DO: see previos TO DO, for now we hard-code the data n_1,...n_m
    //n.push_back(M);//n_0=m
    //int M = 2;//Trials m in y_1,...y_m
    //y.push_back(M);//y_0=m
    //n.push_back(M);//n_0=m
    //y.push_back(50);y.push_back(50);n.push_back(100);n.push_back(100);
    //y.push_back(1264);y.push_back(1542);n.push_back(1600);n.push_back(1600);
    //y.push_back(1264);y.push_back(1542);n.push_back(2000);n.push_back(2000);
    //y.push_back(1264);y.push_back(1542);n.push_back(5000);n.push_back(5000);
    //y.push_back(1264);y.push_back(1542);n.push_back(10000);n.push_back(10000);
    //y.push_back(1264);y.push_back(1542);n.push_back(250000);n.push_back(250000);
    
//---Cossoni and Veronese 1950 pine seedling mortality data, see Peter Green 1995
    int M = 4;//Trials m in y_1,...y_m
    y.push_back(M);//y_0=m
    n.push_back(M);//n_0=m
    //y.push_back(59);y.push_back(89);y.push_back(88);y.push_back(95);//LH,LD,SH,SD
    //n.push_back(100);n.push_back(100);n.push_back(100);n.push_back(100);

/* black robin data Rangatira unmanaged 1990-91 -- 1993-94 Euan Kennedy data 2007 PhD thesis*/
    //y.push_back(37);
    //y.push_back(35);y.push_back(42);y.push_back(38);//LH,LD,SH,SD
    //y.push_back(51);
    //n.push_back(51);
    //n.push_back(53);n.push_back(49);n.push_back(61);
    //n.push_back(65);
    //y.push_back(38);y.push_back(35);y.push_back(211);y.push_back(93);
    //n.push_back(61);n.push_back(53);n.push_back(302);n.push_back(114);
    y.push_back(38);y.push_back(65);y.push_back(51);y.push_back(42);
    n.push_back(61);n.push_back(93);n.push_back(73);n.push_back(49);

/*---2*2*2*2=16 data sets for the Bernoulli Partition Model <=> n_1=n_2=...=n_m
    int M = 4;//Trials m in y_1,...y_m
    y.push_back(M);//y_0=m
    n.push_back(M);//n_0=m
    y.push_back(1);y.push_back(1);y.push_back(1);y.push_back(1);
    y.push_back(0);y.push_back(0);y.push_back(0);y.push_back(0);
    n.push_back(1);n.push_back(1);n.push_back(1);n.push_back(1);
*/
/*
  //---Central Indian > 100mm Monsoon rain data between four groups of years:
  //---1951-1962, 1963-1975, 1976-1988, 1989-2000
    int M = 4;//Trials m in y_1,...y_m
    y.push_back(M);//y_0=m
    n.push_back(M);//n_0=m
    y.push_back(598);y.push_back(666);
    y.push_back(774);y.push_back(768);
    n.push_back(120000);n.push_back(130000);
    n.push_back(130000);n.push_back(120000);
*/
    //just output the hard-coded data for now
    cout << "The data for number of successes/heads y_i's are:\n"; 
    for(int i=0;i<=M;i++) cout << y[i] << '\t';cout << "\n"; 
    cout << "The data for number of Bern trials are n_i's are:\n"; 
    for(int i=0;i<=M;i++) cout << n[i] << '\t';cout << "\n";
    
    vector<int> SixZeros(6,0);
    vector<interval> SixIZeros;
    vector<real> SixRZeros;
    for(int i=0;i<6;i++) SixIZeros.push_back(interval(0.0,0.0));
    for(int i=0;i<6;i++) SixRZeros.push_back(real(0.0));
    for(int i=0;i<6;i++) MbyBMatrix.push_back(SixZeros);
    MbyBMatrix[0][0] = 0;//m=0
    MbyBMatrix[1][0] = 1;//m=1, m-Bell=1
    MbyBMatrix[1][1] = 1;//m=1,b=1, m,b-Stirling2nd=1
    MbyBMatrix[2][0] = 2;//m=2, m-Bell=2
    MbyBMatrix[2][1] = 1;//m=2,b=1, m,b-Stirling2nd=1
    MbyBMatrix[2][2] = 1;//m=2,b=2, m,b-Stirling2nd=1
    MbyBMatrix[3][0] = 5;//m=3, m-Bell=5
    MbyBMatrix[3][1] = 1;//m=3,b=1, m,b-Stirling2nd=1
    MbyBMatrix[3][2] = 3;//m=3,b=2, m,b-Stirling2nd=3
    MbyBMatrix[3][3] = 1;//m=3,b=3, m,b-Stirling2nd=1
    MbyBMatrix[4][0] =15;//m=4, m-Bell=15
    MbyBMatrix[4][1] = 1;//m=4,b=1, m,b-Stirling2nd=1
    MbyBMatrix[4][2] = 7;//m=4,b=2, m,b-Stirling2nd=7
    MbyBMatrix[4][3] = 6;//m=4,b=3, m,b-Stirling2nd=6
    MbyBMatrix[4][4] = 1;//m=4,b=4, m,b-Stirling2nd=1
    MbyBMatrix[5][0] =52;//m=5, m-Bell=52
    MbyBMatrix[5][1] = 1;//m=5,b=1, m,b-Stirling2nd=1
    MbyBMatrix[5][2] =15;//m=5,b=2, m,b-Stirling2nd=15
    MbyBMatrix[5][3] =25;//m=5,b=3, m,b-Stirling2nd=25
    MbyBMatrix[5][4] =10;//m=5,b=4, m,b-Stirling2nd=10
    MbyBMatrix[5][5] = 1;//m=5,b=5, m,b-Stirling2nd=1
    for(int i=0;i<31;i++) Blocks.push_back(SixZeros);// Blocks = [[0,0,0,0,0,0],...] thera are 31 for m=5
    // Partns = [[0,0,0,0,0,0],...] there are 76 models m=1,2,3,4,5
    for(int i=0;i<76;i++) Partns.push_back(SixZeros);
    for(int i=0;i<76;i++) Succss.push_back(SixZeros);
    for(int i=0;i<76;i++) Trials.push_back(SixZeros);
    for(int i=0;i<76;i++) Failrs.push_back(SixZeros);
    for(int i=0;i<76;i++) Ibc.push_back(SixIZeros);
    for(int i=0;i<76;i++) Rbc.push_back(SixRZeros);
    //the indices of partition matrix is the model label for the model partition
    Blocks[0][0]=1;Blocks[0][1]=1;//{1} is encoded as [1,1,0,0,0,0], block of size 1 with trial label 1
    Blocks[1][0]=1;Blocks[1][1]=2;//{2} is encoded as [1,2,0,0,0,0], block of size 1 with trial label 2
    Blocks[2][0]=2;Blocks[2][1]=1;Blocks[2][2]=2;//{1,2}
    Blocks[3][0]=1;Blocks[3][1]=3;//{3}
    Blocks[4][0]=2;Blocks[4][1]=1;Blocks[4][2]=3;//{1,3}
    Blocks[5][0]=2;Blocks[5][1]=2;Blocks[5][2]=3;//{2,3}
    Blocks[6][0]=3;Blocks[6][1]=1;Blocks[6][2]=2;Blocks[6][3]=3;//{1,2,3}
    Blocks[7][0]=1;Blocks[7][1]=4;//{4}
    Blocks[8][0]=2;Blocks[8][1]=1;Blocks[8][2]=4;//{1,4}
    Blocks[9][0]=2;Blocks[9][1]=2;Blocks[9][2]=4;//{2,4}
    Blocks[10][0]=2;Blocks[10][1]=3;Blocks[10][2]=4;//{3,4}
    Blocks[11][0]=3;Blocks[11][1]=1;Blocks[11][2]=2;Blocks[11][3]=4;//{1,2,4}
    Blocks[12][0]=3;Blocks[12][1]=1;Blocks[12][2]=3;Blocks[12][3]=4;//{1,3,4}
    Blocks[13][0]=3;Blocks[13][1]=2;Blocks[13][2]=3;Blocks[13][3]=4;//{2,3,4}
    Blocks[14][0]=4;Blocks[14][1]=1;Blocks[14][2]=2;Blocks[14][3]=3;Blocks[14][4]=4;//{1,2,3,4}
     //{1,2} is encoded as [2,1,2,0,0,0], block of size 2 with trial labels 1 and 2
    if(M==1){//just 1 binomial trial y_1
      Partns[0][0]=1;Partns[0][1]=0;// model 0 = {{1}}
    }
    if(M==2){//just 2 binomial trials y_1,y_2
      Partns[0][0]=1;Partns[0][1]=2;//model 0 = {{1,2}}
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=1;//model 1 = {{1},{2}}
    }
    if(M==3){//just 3 binomial trials y_1,y_2,y_3
      Partns[0][0]=1;Partns[0][1]=6;//model 0 = {{1,2,3}}
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=5;//model 1 = {{1},{2,3}}
      Partns[2][0]=2;Partns[2][1]=1;Partns[2][2]=4;//model 2 = {{2},{1,3}}
      Partns[3][0]=2;Partns[3][1]=3;Partns[3][2]=2;//model 3 = {{3},{1,2}}
      Partns[4][0]=3;Partns[4][1]=0;Partns[4][2]=1;Partns[4][3]=3;//model 4 = {{1},{2},{3}}
    }
    if(M==4){//just 4 binomial trials y_1,y_2,y_3,y_4
      //model label 0 = [{1, 2, 3, 4}] 
      // model label 1 = [{1}, {2, 3, 4}], 
      // model label 2 = [{2}, {1, 3, 4}], 
      // model label 3 = [{3}, {1, 2, 4}], 
      // model label 4 = [{4}, {1, 2, 3}], 
      // model label 5 = [{1, 2}, {3, 4}],
      // model label 6 = [{1, 3}, {2, 4}], 
      // model label 7 = [{1, 4}, {2, 3}], 
      // model label 8 = [{1}, {2}, {3, 4}], 
      // model label 9 = [{1}, {3}, {2, 4}],
      // model label 10 = [{1}, {4}, {2, 3}], 
      // model label 11 = [{2}, {3}, {1, 4}], 
      // model label 12 = [{2}, {4}, {1, 3}], 
      // model label 13 = [{3}, {4}, {1, 2}], 
      // model label 14 = [{1}, {2}, {3}, {4}], 
      //model 0 = [{1, 2, 3, 4}] 
      Partns[0][0]=1;Partns[0][1]=14;
      // model 1 = [{1}, {2, 3, 4}], 
      Partns[1][0]=2;Partns[1][1]=0;Partns[1][2]=13;
      // model 2 = [{2}, {1, 3, 4}], 
      Partns[2][0]=2;Partns[2][1]=1;Partns[2][2]=12;
      // model 3 = [{3}, {1, 2, 4}], 
      Partns[3][0]=2;Partns[3][1]=3;Partns[3][2]=11;
      // model 4 = [{4}, {1, 2, 3}], 
      Partns[4][0]=2;Partns[4][1]=7;Partns[4][2]=6;
      // model 5 = [{1, 2}, {3, 4}],
      Partns[5][0]=2;Partns[5][1]=2;Partns[5][2]=10;
      // model 6 = [{1, 3}, {2, 4}], 
      Partns[6][0]=2;Partns[6][1]=4;Partns[6][2]=9;
      // model 7 = [{1, 4}, {2, 3}], 
      Partns[7][0]=2;Partns[7][1]=8;Partns[7][2]=5;
      // model 8 = [{1}, {2}, {3, 4}], 
      Partns[8][0]=3;Partns[8][1]=0;Partns[8][2]=1;Partns[8][3]=10;
      // model 9 = [{1}, {3}, {2, 4}],
      Partns[9][0]=3;Partns[9][1]=0;Partns[9][2]=3;Partns[9][3]=9;
      // model 10 = [{1}, {4}, {2, 3}], 
      Partns[10][0]=3;Partns[10][1]=0;Partns[10][2]=7;Partns[10][3]=5;
      // model 11 = [{2}, {3}, {1, 4}], 
      Partns[11][0]=3;Partns[11][1]=1;Partns[11][2]=3;Partns[11][3]=8;
      // model 12 = [{2}, {4}, {1, 3}], 
      Partns[12][0]=3;Partns[12][1]=1;Partns[12][2]=7;Partns[12][3]=4;
      // model 13 = [{3}, {4}, {1, 2}], 
      Partns[13][0]=3;Partns[13][1]=3;Partns[13][2]=7;Partns[13][3]=2;
      // model 14 = [{1}, {2}, {3}, {4}], 
      Partns[14][0]=4;Partns[14][1]=0;Partns[14][2]=1;Partns[14][3]=3;Partns[14][4]=7;
    }
    int NumModels=MbyBMatrix[M][0];// numer of models based on the number of trials M
    //cout << NumModels << endl; getchar();
    for(int m=0; m < NumModels; m++){
      int ModelDim = Partns[m][0];//number of dimensions (or blocks in partition) of model m
      ivector domain (1, ModelDim); // vector of intervals, with indices 1...2
      LabBox Ldomain;
      for (int i = 1; i <= ModelDim; i++){
        //if working with log ikelihood then domain should be a 
        //proper subset of [0,1] to ensure local Lipschitz condition
        //let it approach [0,1] as needed by your problem
        if(LogPi) {
          //getchar(); RHS is just 1 when SmallestPositive is too small
          domain[i] = interval(SmallestPositive, (1.0-SmallestPositive));
        }
        else domain[i] = interval (0.0, 1.0);
      }
      Ldomain.Box = domain;
      Ldomain.L = m;
      //cout << "Ldomain.L: " << Ldomain.L << endl;
      LabDomainList.push_back (Ldomain);
      //model prior weights are the same
      modelprior.push_back(1.0);
      Succss[m][0]= ModelDim;
      for(int Blki=1; Blki<= Succss[m][0]; Blki++){
        int BlkIndx = Partns[m][Blki];
        int BlkiSz = Blocks[BlkIndx][0];
        Succss[m][Blki]=0;Trials[m][Blki]=0;Failrs[m][Blki]=0;
        interval ibc(1.0,1.0);
        for(int j=1; j<= BlkiSz; j++){
          Trials[m][Blki] += n[ Blocks[BlkIndx][j]];
          Succss[m][Blki] += y[ Blocks[BlkIndx][j]];
          if(!(IgnoreBinomCoeffs)) 
            ibc *= binary_coefficient(n[Blocks[BlkIndx][j]],y[Blocks[BlkIndx][j]]);
        }
        Failrs[m][Blki] = Trials[m][Blki] - Succss[m][Blki];
        Ibc[m][Blki] = ibc;
        Rbc[m][Blki] = (0.5*(Inf(ibc) + Sup(ibc)));
/*
        cout << "m = " << m << " Blki = " << Blki << '\t'
             <<  Succss[m][Blki] << '\t' << Failrs[m][Blki] << '\t' << Trials[m][Blki]  
             << '\t' << Ibc[m][Blki] << '\t' << Rbc[m][Blki] << endl;
        getchar();
*/
      }
    }
    cout << "end of FBinomialPartition constructor. \n";

}

// vector<LabBox> FBinomialPartition::get_domain(){ return LabDomainList; }

interval
FBinomialPartition::operator () (const LabBox & X)
const
{
    // L(n, y, q1) = power(q1, y) * power(1.0 - q1, n - y) * bc[0];
    n_interval_calls++;
    interval result;
    int ModelDim = Partns[X.L][0];
    // 1.0/(unit volume in any dimension, ie, uniform on unit hypercube)
    interval prior_density(1.0, 1.0); //change prior density here if needed
    ///ivector Xc(1, ModelDim); // vector of intervals, with indices 1...2
    ///for(int i=1; i<=ModelDim;i++){
    ///  Xc[i]=1.0-X.Box[i];if(Inf(Xc[i]) <= 0.0){ SetInf(Xc[i], 0.0); }
    ///}
    if(UsingLogDensity){
      interval LogLikelihood(0.0,0.0); 
      for(int i=1; i<=ModelDim; i++){
        interval q = X.Box[i];
        interval qc = 1.0 - q; if(Inf(qc) <= 0.0) SetInf(qc, SmallestPositive);
        LogLikelihood += Succss[X.L][i]*ln(q) + (Failrs[X.L][i])*ln(qc) + ln(Ibc[X.L][i]);
        ///LogLikelihood += Succss[X.L][i]*ln(X.Box[i]) + (Failrs[X.L][i])*ln(Xc[i]) + ln(Ibc[X.L][i]);
      }
      result = ln(modelprior[X.L]) + LogLikelihood + ln(prior_density);
    }
    else{
      interval likelihood(1.0,1.0);
      for(int i=1; i<=ModelDim; i++){
        interval q = X.Box[i];
        interval qc = 1.0 - q; if(Inf(qc) <= 0.0){ SetInf(qc, 0.0); }
        likelihood *= power(q,Succss[X.L][i]) * power(qc,Failrs[X.L][i]) * Ibc[X.L][i];
        //cout << "model = " << X.L << '\t' << "S F Ibc = " << Succss[X.L][i] << '\t' << Failrs[X.L][i] << '\t' << Ibc[X.L][i] << endl; getchar();
        //likelihood *= power(X.Box[i],Succss[X.L][i]) * power(Xc[i],Failrs[X.L][i]) * Ibc[X.L][i];
      }
      result = modelprior[X.L]*likelihood*prior_density;
    }
    //cout << result << '\n';getchar();
    return result;
}

real
FBinomialPartition::operator () (const LabPnt & X)
    const
{    
    n_real_calls++;
    real result;
    int ModelDim = Partns[X.L][0];
    real prior_density = 1.0;  // 1.0/(unit volume)
    if(UsingLogDensity){
      real LogLikelihood = 0.0; 
      for(int i=1; i<=ModelDim; i++){
        real q = X.Pnt[i];
/*//////////////////////-just for monsoon data as it is too small likelihood
        real qc = 1-q;
        if (q<=0) q=SmallestPositive;
        if (qc<=0) qc=SmallestPositive;
        LogLikelihood += Succss[X.L][i]*ln(q) + Failrs[X.L][i]*ln(qc) + ln(Rbc[X.L][i]);
*//////////////////////end - just for monsoon data as it is too small likelihood
///////////////////////-usual case; if likelihood is NOT too small for screen
        LogLikelihood += Succss[X.L][i]*ln(q) + Failrs[X.L][i]*ln(1.0-q) + ln(Rbc[X.L][i]);
///////////////////////-end of usual case; 
      }
      result = ln(modelprior[X.L]) + LogLikelihood + ln(prior_density);
    }
    else{
      real likelihood=1.0;
      for(int i=1; i<=ModelDim; i++){
        real q = X.Pnt[i];
        likelihood *= power(q,Succss[X.L][i]) * power(1.0-q,Failrs[X.L][i]) * Rbc[X.L][i];
      }
      result = modelprior[X.L]*likelihood*prior_density;
    }
    return result;
}

interval
FBinomialPartition::binary_coefficient(int N, int Y)
{
    int kmax = (Y <= N - Y)? Y: N - Y;
    int nf = N; 
    int df = 1;
    interval nmrtr(1,1);
    interval dnmntr(1,1);
    for(int k = 1; k <= kmax; k++){
        nmrtr *= nf; nf--;
        dnmntr *= df; df++;
    }
    return nmrtr/dnmntr;
} 

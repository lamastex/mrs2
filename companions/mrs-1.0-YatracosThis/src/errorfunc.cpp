/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
* Copyright (C) 2011 Gloria Teng
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

/*! \file errorfunc.cpp
\brief L1 error function definitions
*/

#include "errorfunc.hpp"
#include "spsnode.hpp"
#include "adaptivehistogram.hpp"
#include "AIAsubpaving.hpp" // to use Volume
#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <vector>   // to use the  stl::vector container
#include <set>      // to use the stl::multiset container
#include <algorithm>// to use stl::algorithms
#include <list>     // to use stl:: lists
#include <map>      // to use stl::maps
#include <fstream>  // for ifstream, ofstream
#include <sstream>  // to be able to manipulate strings as streams
#include <iomanip> // format manipulation on streams
#include <exception> // use exceptions
#include <math.h> // math library
#include <gsl/gsl_math.h> // to use the constant M_PI 
#include <gsl/gsl_randist.h>
//these are needed for 2D taylor integration
#include "../examples/StatsSubPav/ExactInt/Int.h"  
#include "../examples/StatsSubPav/ExactInt/dim2taylor.hpp"
// to access gsl_matrix elements
#include "gsl/gsl_matrix.h"
// to perform vector-matrix operations
#include "gsl/gsl_blas.h"

using namespace subpavings;
using namespace std;

//-------------IAE for Bivariate gaussian, Levy 2D, Rosen 2D----------------//
/*! \brief IAE for the bivariate gaussian
*/
d2t BiGOP (d2tv X, interval fhat) 
{
	//cout << "Calling BiGOP: " << endl;
   // Parameters specific to the Bivariate Gaussian target
	real rsigma_x = 1.0;
	real rsigma_y = 1.0;
	real rrho = 0;
  d2t f = taylor::init_const(X[1].order(),interval(0.0));
  real det = 1.0/(2*M_PI*rsigma_x*rsigma_y*sqrt(1-sqr(rrho)));
  f = sqr(X[1]/rsigma_x) + sqr(X[2]/rsigma_y) - 
            (2*rrho*X[1]*X[2])/(rsigma_x*rsigma_y);
  f = det * exp (-((1.0/2*(1-sqr(rrho))) * f));
  
  d2t result = taylor::init_const(X[1].order(),interval(0.0));

  //split the integrand to get positive values only (absolute values) 
  if ( (Sup(f[0][0]) < Inf(fhat)) ) {	  
	  result = fhat - f;
	  //cout << "fhat - f: " << result[0][0] << "\n" << endl;
  }  
  else if ((Sup(f[0][0]) > Inf(fhat))) { 
	   result = f - fhat; 
		//cout << "f-fhat: " << result[0][0] << "\n" << endl;
  }   
  return result;
}

/*! \brief IAE for the Levy 2D
*/
d2t LevyOP (d2tv X, interval fhat) 
{
	//cout << "Calling LevyOP: " << endl;
  // Parameters specific to the Levy target 
	real Temperature = 40.0;
	real Center1 = 1.42513; 
	real Center2 = 0.80032; 
	real GlobalMax = 176.14;
  
  d2t isum = taylor::init_const(X[1].order(),interval(0.0));
  d2t jsum = taylor::init_const(X[1].order(),interval(0.0));

  for (int i = 1; i <= 5; i++)
  {
    isum = isum + i * cos ((i - 1) * X[1] + (i));
    jsum = jsum + i * cos ((i + 1) * X[2] + (i));
  }
                    // Avoid real conversion error
  d2t hh = isum * jsum + sqr (X[1] + Center1) +
    sqr (X[2] + Center2);
  hh = hh + GlobalMax;  
  // TEMPERATURE = 1, 4, 40, 400, 4000
  d2t f = exp (-hh / Temperature);

  //integrand  
  d2t result = taylor::init_const(X[1].order(),interval(0.0));
  //split the integrand to get positive values only (absolute values) 
  if ( (Sup(f[0][0]) < Inf(fhat)) ) {	  
	  result = fhat - f;
	  //cout << "fhat - f: " << result[0][0] << "\n" << endl;
  }  
  else if ((Sup(f[0][0]) > Inf(fhat))) { 
	   result = f - fhat; 
		//cout << "f-fhat: " << result[0][0] << "\n" << endl;
  }   
  
  return result;
}

/*! \brief IAE for a Rosenbrock 2D.
*/
d2t RosenOP (d2tv X, interval fhat) 
{
	//cout << "Calling RosenOP: " << endl;
	// Parameters specific to the Rosenbrock 2D target 
	real Tinverse = 1.0;
	real Height = 100.0;
  
  d2t f = taylor::init_const(X[1].order(),interval(0.0));
  for (int i = 1; i < 2; i++) //2nd term should be size_k
    {
      f = f + (Height * sqr(X[i+1] - sqr(X[i])) +
        sqr(X[i] - 1.0));
    }
  f = exp (-(Tinverse * f));

  //integrand  
  d2t result = taylor::init_const(X[1].order(),interval(0.0));
  //split the integrand to get positive values only (absolute values) 
  if ( (Sup(f[0][0]) < Inf(fhat)) ) {	  
	  result = fhat - f;
	  //cout << "fhat - f: " << result[0][0] << "\n" << endl;
  }  
  else if ((Sup(f[0][0]) > Inf(fhat))) { 
	   result = f - fhat; 
		//cout << "f-fhat: " << result[0][0] << "\n" << endl;
  }   
  return result;
}
//-----------------End of 2D integrations----------------------------------//

//-------------Family of functions for finite mixtures---------------------//
/*! \brief Get the probability density function of a finite mixture r.v.
*/
double FinMixPDF(double x, vector<double> &W, vector<double> &M, 
                          vector<double> &S)
{
int Ncomp = W.size();
double PDF = 0;
	for (int c=0; c < Ncomp; c++){
		double z = pow((x-M[c])/S[c], 2);
		PDF += W[c]*exp(-0.5*z)/(S[c]*sqrt(2*M_PI));
	}  
return PDF;
}

/*! \brief Get the cumulative distribution function of a finite mixture r.v.
*/
double FinMixCDF(double x1, double x2, vector<double>& W, 
                          vector<double>& M, vector<double>& S)
{
double Ncomp = W.size();
double CDF = 0;
for (int c=0; c < Ncomp; c++){
CDF += 0.5*(1 + erf((x2-M[c])/S[c]/sqrt(2))) -
       0.5*(1 + erf((x1-M[c])/S[c]/sqrt(2)));
//cout << "CDF is: " << CDF << endl;
}
return CDF;
}

/*! \brief Get the absolute error of a finite mixture at x.
*/
double FinMixAbs(double x, void * params)
{
FinMix mixt = *(FinMix *) params;
double FinMixAbs = fabs(mixt.fhat - FinMixPDF(x, mixt.W, mixt.M, mixt.S));
return FinMixAbs;
}

/*! \brief Calculate the IAE at boundaries of a finite mixture.
*/
dotprecision dpFinMixIAEBoun(double x1, double x2, FinMix& mixt)
{	
	dotprecision dpFinMixIAEBoun;
	dpFinMixIAEBoun = 0.0;
	double Ncomp = (mixt.W).size();
	double cdfLeft = 0.0;
	double cdfRight = 0.0;
	int c;

	for (c=0; c < Ncomp; c++) {
		cdfLeft += mixt.W[c]*0.5*(1 + erf((x1-mixt.M[c])/mixt.S[c]/sqrt(2)));
		cdfRight += 1-mixt.W[c]*0.5*(1 + erf((x2-mixt.M[c])/mixt.S[c]/sqrt(2)));
	}

	accumulate(dpFinMixIAEBoun, cdfLeft, 1.0);
	accumulate(dpFinMixIAEBoun, cdfRight, 1.0);

	return dpFinMixIAEBoun;
} 

/*! \brief Function to find the number of generated U(0,1) data that is less 
     the weight of the components.
*/
void findComp(vector<double> & u, int& intp, const int n, 
              double* weight, int w)
{     
	// set up an array that checks whether a U(0,1) r.v. is a member of that 
	// component or not. '1' if true; '0' if false.
	//vector<int> u_one(n);
	
	int j;
	// check if the U(0,1) r.v. is a member of the component or not.
	if (w==0) {
		for (j=0; j<n; j++) {
			if (u[j] <= *(weight + w))
			{ intp++; }
			//cout<<u_one[j1]<<endl;
		}
	}
	else {
	//cout << "Checking for members between " << *(WeightCumPtr + m-1) << 
	//" and " << *(WeightCumPtr + m) << endl;
		for (j=0; j<n; j++) {
			bool a = (u[j] > (*(weight + w-1)));
			bool b = ( u[j] <= (*(weight + w)));
			if (a==1 && b==1) { intp++; }
			//cout<<u_one[j]<<endl;
		}
	}

	// count how many '1's are there in u_rv. This is the number of members 
	//in component m.
	//intp = (int) count(u_one, u_one+n, 1);
} //end find_comp()

/*! \brief Function to get cumulative sum*/
void cumsum(vector<double> weight, double* w) {
	double cum =0;
	for (size_t j = 0; j < weight.size(); j++) {
		cum += weight[j];
		*w++ = cum;
	}
} // end of function cumsum()




//-----------------End of family of functions for finite mixtures-----------//

//---------------Family of functions for regular histograms-----------------//
/*! \brief Make a regular histogram (using Dominic's histogram.m file)
*/
void makeRegularHist(RegHist& myRegHist, const RVecData& sortedData, 
								ivector theBox, int bwmethod)
{
	// vector for heights
	vector<double> heights;

	size_t n = sortedData.size();
	//cout << "there are " << n << " points." << endl;

	// determine the bandwidth
	double n3 = pow(n, -1.0/3.0); 
	double n5 = pow(n, -0.2); 
	double n7 = pow(n, -1.0/7.0);
	//cout << "n3: " << n3 << "\tn5: " << n5 << "\tn7: " << n7 << endl;

	/*
	// put the data into a list to sort the data
	RVecDataCItr rvecIt;
	list<rvector> rvecList;
	list<rvector>::iterator rvecListIt;
	for (rvecIt = rvec.begin(); rvecIt < rvec.end(); rvecIt++){
		rvector thisrv(1);
		thisrv = *rvecIt;
		rvecList.push_back(thisrv);
	}
	rvecList.sort(); //sort the data
	RVecData sortedData; // put back into vector
	for (rvecListIt = rvecList.begin(); rvecListIt != rvecList.end(); rvecListIt++)
	{ 
		rvector thisrv(1);
		thisrv = *rvecListIt;
		sortedData.push_back(thisrv);
	}
	*/
	
	// get the interquartile range
	int upperQ = ceil(0.75*n);  //upper quartile
	int lowerQ = ceil(0.25*n);  //lower quartile
	double upperQx = _double((sortedData[upperQ-1])[1]);
	double lowerQx = _double((sortedData[lowerQ-1])[1]);
	double xiq = upperQx - lowerQx; // interquartile range

	double xsd = myStd(sortedData);
	// determine which sigma to use
	double sigma;
	if (xiq == 0) { sigma = xsd; }
	else { sigma = min(xsd, (xiq/1.349)); }

	// determine which bandwidth method to use
	double bw = 0.0;
	
	cout << "determine which bandwidth method to use:" << endl;
	// Scott if bwmethod == 0
	if (bwmethod == 0) { bw = 3.4908 * sigma * n3; } 
	// Wand's one stage if bwmethod == 1 
	else if (bwmethod == 1) {
		double g11 = 1.3041 * sigma * n5;
		bw  = 1.8171 * pow(-psi(sortedData, g11, 2), (-1.0/3.0)) * n3;
	}
	// Wand's two stage if bwmethod == 2
	else if (bwmethod == 2) {
		double g22 = 1.2407 * sigma * n7;
		double g21 = 0.9558 * pow((psi(sortedData,g22,4)),(-0.2)) * n5; 
		bw  = 1.8171 * pow((-psi(sortedData,g21,2)),(-1.0/3.0)) * n3;
	}

	// Determine bin origin
	//cout << "determine the bin origin" << endl;
	rvector xmin = sortedData[0]; //the minimum value;
	rvector xmax = sortedData[n-1]; //the maximum value
	double xrange =  _double(xmax[1]) - _double(xmin[1]); // range of data	
	int nbin  = ceil(xrange/ bw); // number of bins
	real xoffset = _real((nbin * bw - xrange) / 2); // offset
	rvector xlow = Inf(theBox);
	rvector xupp = Sup(theBox);
	real bbeg = max(xlow[1], xmin[1] - xoffset); // bin origin
	real bend = min(xupp[1], xmax[1] + xoffset); // bin end
	real bwR = (bend - bbeg) / (nbin*1.0); // binwidth
	myRegHist.binwidth = bwR;

	cout << "there are " << nbin << " bins" << endl;
	int J = 0;
	for (int i = 0; i < nbin; i++) {
		// bin edges
		myRegHist.LowerBoxes.push_back(bbeg + bwR*i);
		myRegHist.UpperBoxes.push_back(bbeg + bwR*(i+1));

		//cout << "getting the counts:" << endl;
		size_t P = 0;
		for (size_t j = J; j < n; j++) {
			rvector thisrv(1);
			thisrv = sortedData[j];
			if (thisrv[1] >= (bbeg + bwR*i) && thisrv[1] < (bbeg + bwR*(i+1)) ) {
				P += 1; // Count frequencies:
			}
			else { J = j+1; break; }
		}
		myRegHist.heights.push_back((P*1.0)/(n*1.0*_double(bwR)));  //height	
	}
}  

void makeRegularHist(RegHist& myRegHist, const RVecData& sortedData, 
								ivector theBox, double bw)
{
	// vector for heights
	vector<double> heights;

	size_t n = sortedData.size();
	//cout << "there are " << n << " points." << endl;

	/*
	// put the data into a list to sort the data
	RVecDataCItr rvecIt;
	list<rvector> rvecList;
	list<rvector>::iterator rvecListIt;
	for (rvecIt = rvec.begin(); rvecIt < rvec.end(); rvecIt++){
		rvector thisrv(1);
		thisrv = *rvecIt;
		rvecList.push_back(thisrv);
	}
	rvecList.sort(); //sort the data
	RVecData sortedData; // put back into vector
	for (rvecListIt = rvecList.begin(); rvecListIt != rvecList.end(); rvecListIt++)
	{ 
		rvector thisrv(1);
		thisrv = *rvecListIt;
		sortedData.push_back(thisrv);
	}	*/

	// get the interquartile range
	int upperQ = ceil(0.75*n);  //upper quartile
	int lowerQ = ceil(0.25*n);  //lower quartile
	double upperQx = _double((sortedData[upperQ-1])[1]);
	double lowerQx = _double((sortedData[lowerQ-1])[1]);
	double xiq = upperQx - lowerQx; // interquartile range

	// determine which sigma to use
	double xsd = myStd(sortedData); // standard deviation
	double sigma;
	if (xiq == 0) { sigma = xsd; }
	else { sigma = min(xsd, (xiq/1.349)); }

	// Determine bin origin
	//cout << "determine the bin origin" << endl;
	rvector xmin = sortedData[0]; //the minimum value;
	rvector xmax = sortedData[n-1]; //the maximum value
	double xrange =  _double(xmax[1]) - _double(xmin[1]); // range of data	
	int nbin  = ceil(xrange/ bw); // number of bins
	real xoffset = _real((nbin * bw - xrange) / 2); // offset
	rvector xlow = Inf(theBox);
	rvector xupp = Sup(theBox);
	real bbeg = max(xlow[1], xmin[1] - xoffset); // bin origin
	real bend = min(xupp[1], xmax[1] + xoffset); // bin end
	real bwR = (bend - bbeg) / (nbin*1.0); // binwidth
	myRegHist.binwidth = bwR;

	int J = 0;
	for (int i = 0; i < nbin; i++) {
		// bin edges
		myRegHist.LowerBoxes.push_back(bbeg + bwR*i);
		myRegHist.UpperBoxes.push_back(bbeg + bwR*(i+1));

		//cout << "getting the counts:" << endl;
		size_t P = 0;
		for (size_t j = J; j < n; j++) {
			rvector thisrv(1);
			thisrv = sortedData[j];						
			if (thisrv[1] >= (bbeg + bwR*i) && thisrv[1] < (bbeg + bwR*(i+1)) ) {
				P += 1; // Count frequencies:
			}
			else { J = j+1; break; }
		}
		myRegHist.heights.push_back((P*1.0)/(n*1.0*_double(bwR)));  //height	
	}
}  

/*! \brief Function required for regular histogram
*/
double psi(const RVecData& rvec, double g, double r)
{
	// value to be returned by this function
	double myPsi = 0;
	
	//sample size
	int n = rvec.size();
	
	//put the sample into a vector<double>
	//cout << "data in rvector type:" << endl;
	//vector<double> theData;
	for (int i=0; i < n; i++){
		rvector thisrv = rvec[i];
	//	cout << thisrv[1] << "\t";
	}
//	cout << "\n" << endl;
	
	//data-based value
	double c = pow((double(n)),(-2)) * pow(g,(-r-1));
   //cout << "C: " << c << endl;
	
   if (n < 1000) {		   
			vector<double> XX, XX2;
			
			for (int i = 0; i < n; i++) {
				rvector rv1 = rvec[i];
				for (int j=0; j < n; j++) {
					rvector rv2 = rvec[j];		
					XX.push_back((_double(rv2[1])-_double(rv1[1]))/g);			
					XX2.push_back(pow((_double(rv2[1])-_double(rv1[1]))/g, 2));
				}			
			}
		
		// get the normal probability density values corresponding to X
		vector<double> vecNPD;
		vecNPD = gaussian(vecNPD, XX, 1);
	 
		if (r == 2) {
		   for (int i = 0; i < XX2.size(); i++) {
				myPsi += (XX2[i]-1)*vecNPD[i];
			}
			
			myPsi *= c;
		//	cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		} // end of r==2
		
		
		else if (r == 4) {
			for (int i = 0; i < XX2.size(); i++) {
			   double XX4 = XX2[i] * XX2[i];
				myPsi += (XX4 - 6*XX2[i] + 3) * vecNPD[i];
			}			
			myPsi *= c;
		//	cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		}
		
		else if (r == 6) {
			for (int i = 0; i < XX2.size(); i++) {
			   double XX4 = XX2[i] * XX2[i];
				double XX6 = XX4 * XX2[i];
				myPsi += (XX6 - 15*XX4 + 45*XX2[i] - 15) * vecNPD[i];
			}			
			myPsi *= c;
		//	cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		}
		
		else {
			cout << "Error: Input r for Function PSI must be 2, 4 or 6." << endl;
		}
	
	} // end of if n < 1000
	
	else { // n >= 1000   
		rvector xmin = rvec[0]; //the minimum value;
		rvector xmax = rvec[n-1]; //the maximum value
		int m = 500; 
		rvector d =  (xmax(1) - xmin(1)) / (m - 1);
		
		/*
		rmatrix Data(n, 1);
		rmatrix Ones(n, 1);
		for (int i = 0; i < n; i++) {
			Data[i+1][1] = (rvec[i][1] - xmin[1])/d[1];
			Ones[i+1][1] = 1;
		}*/

		cout << "get c" << endl;
		vector<rvector> C; 
		for (int j = 0; j < m; j++) {

			rvector indC(1);
			indC[1] = 0;

			/*rmatrix J(n, 1);
			for (int i = 0; i < n; i++) {
				J[i+1][1] = -j + 2;
			}*/

			//else { //cout << ((Ones - abs(Data + J)) >= 0) << endl; } 
			for (int i=0; i < n; i++) {
				rvector newX(1);
				rvector thisrv = rvec[i];
				newX[1] = thisrv[1] - xmin[1]; 
			 	indC[1] += max((1 - abs((newX[1] / d[1]) - j + 1 + 1)), 0);
			}
				C.push_back(indC);
			//}
		}

		cout << "get cc" << endl;
		vector<rvector> CC;
		for (int i = 0; i < C.size(); i++) {
			for (int j=0; j < C.size(); j++) {
				rvector thisrv(1);
				thisrv[1] = C[j]*C[i];
				CC.push_back(thisrv);
			}
		}
	  
		cout << "get jj and jj2" << endl;
		vector<double> JJ, JJ2;  		
		for (int i = 0; i < m; i++) {
			rvector rv1(1);
			rv1[1] = i;
			for (int j=0; j < m; j++) {
				rvector rv2(1);
				rv2[1] = j;		
				JJ.push_back(_double(d[1])*( _double(rv2[1]) - _double(rv1[1]))/g);			
				JJ2.push_back(pow(_double(d[1])*(_double(rv2[1]) - _double(rv1[1]))/g, 2));
				//cout << "JJ: " << _double(d[1])*( _double(rv2[1]) - _double(rv1[1]))/g << endl;
			}			
		}
    
		cout << "get gaussian" << endl;
		vector<double> vecJJNPD;
		vecJJNPD = gaussian(vecJJNPD, JJ,1);
	
		cout << "get cphi " << endl;
		vector<double> CPhi;
	   for (int i = 0; i < vecJJNPD.size(); i++) {
			rvector thisrv(1);
			thisrv = CC[i];
			CPhi.push_back(_double(thisrv[1])*vecJJNPD[i]);
		}
			
		cout << "get psi " << endl;
		if  (r == 2) {		   
			for (int i = 0; i < JJ2.size(); i++) {				
				myPsi += ( JJ2[i] - 1) * CPhi[i];
			}			
			myPsi *= c;
	//		cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		} // end of if r = 2
		 
		else if (r == 4) {		  
		  for (int i = 0; i < JJ2.size(); i++) {
				double JJ4 = JJ2[i] * JJ2[i];
				myPsi += (JJ4 - 6*JJ2[i] + 3) * CPhi[i];
			}			
			myPsi *= c;
	//		cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		} // end of r = 4
		 
		else if (r == 6) {
		  for (int i = 0; i < JJ2.size(); i++) {
				double JJ4 = JJ2[i] * JJ2[i];
				double JJ6 = JJ4 * JJ2[i];
				myPsi += (JJ6 - 15*JJ4 + 45*JJ2[i] - 15) * CPhi[i];
			}			
			myPsi *= c;
		//	cout << "r=" << r << "\tmyPsi: " << myPsi << endl; 
		} // end of r =6 
		
		else {
        cout << "Error: Input r for Function PSI must be 2, 4 or 6." << endl;
      }    
   } // end of n >= 1000 		
		
		return myPsi;
	} // end of function psi

          			
/*! Get mean of RVecData
*/
double myMean(const RVecData& rvec)
{
	//cout << "getmean" << endl;
	size_t n = rvec.size(); // number of points
	// get the sum
	double mySum = 0;
	for (int i=0; i < n; i++){
		rvector thisrv = rvec[i];
		mySum += _double(thisrv[1]);
	}
   double theMean = mySum/n;
	
	return theMean; 
}

/*! Get standard deviation of RVecData
*/
double myStd(const RVecData& rvec)
{
	//cout << "get std dev" << endl;
	size_t n = rvec.size(); // number of points
	double theMean = myMean(rvec);
	// get the sum of squares of the deviation
	double mySquaredSum = 0;
	for (int i=0; i < n; i++){
		rvector thisrv = rvec[i];
		mySquaredSum += pow(_double(thisrv[1]) - theMean, 2);
	}
	double theStd = sqrt(mySquaredSum/(n-1));
	return theStd; 
}

/*! Gaussian probability density function.
*/
std::vector<double> gaussian(std::vector<double> &vecNPD, 
                               std::vector<double> &XX, double SS)
{
	vector<double>::iterator vecIt;
	
	for (vecIt = XX.begin(); vecIt < XX.end(); vecIt++){
		double X = *vecIt;
      vecNPD.push_back(exp(-(X*X)/(2*SS*SS))/(sqrt(2*M_PI)*SS));
	}

 return vecNPD;
}

/*! Get the IAE for a finite mixture distribution for regular histograms
*/
real getRegHistFinMixIAE(size_t n, RegHist & myRegHist, FinMix & mixt)
{
//----------------get the IAE-----------------------------------------------
dotprecision dpIAE, dpIAEBoun;
dpIAE = 0.0;
int Nbin=(myRegHist.UpperBoxes).size();

gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
double result = 0.0;
double error;
gsl_function F;

F.function = &FinMixAbs;
F.params =  &mixt;

for (int j=0; j< Nbin; j++){
  mixt.fhat = myRegHist.heights[j];
  double xupp = _double(myRegHist.UpperBoxes[j]);
  double xlow = _double(myRegHist.LowerBoxes[j]);
  gsl_integration_qags(&F, xlow, xupp, 0, 1e-7, 1000, w, &result, &error);
  accumulate(dpIAE, result, 1.0);
}

/*
// Accounting for the boundaries
ivector theBoxVec;
interval boxes;
//upper bound
theBoxVec = myRegHist.theBoxes[Nbin-1];
boxes = theBoxVec[1];
double xupp1 = _double(Sup(boxes));
//lower bound
theBoxVec = myRegHist.theBoxes[0];
boxes = theBoxVec[1];
double xlow1 = _double(Inf(boxes));
dpIAEBoun = dpFinMixIAEBoun(xlow1, xupp1, mixt);
dpIAE += dpIAEBoun;
*/

// cast dot precision to real
real FinMixIAE = rnd(dpIAE);

// free the workspace
gsl_integration_workspace_free (w);

return FinMixIAE;
}

/*! Get the IAE for a finite mixture distribution for regular histograms
*/
interval getRegHistFinMixIntervalIAE(size_t n, RegHist & myRegHist, FinMix & mixt, double tol, int deg)
{
	interval totalArea(0.0);
	size_t Nbin = myRegHist.heights.size();
	
	for (size_t j=0; j< Nbin; j++){
		double fhat = myRegHist.heights[j];
		double xupp = _double(myRegHist.UpperBoxes[j]);
		double xlow = _double(myRegHist.LowerBoxes[j]);
		
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;

		//---------find the root at this domain
		// make an intervalw object using thisBox
		intervalw thisIntW(xlow, xupp);
		interval thisInt(xlow, xupp);

		// find the root
		//cout << "finding roots at this node " << thisInt << endl;
		bisect(thisIntW, tol, fhat, rootVec, mixt.W, mixt.M, mixt.S); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				//cout << "root " << i << ": " << rootVec[i] << endl;
				//first insert into uniqueRootVec
				uniqueRootVec.push_back(rootVec[i]);
				//cout << i-1 << "\t" << i << "\t" << i+1 << endl;
				//now check for uniqueness
				if (((i-1) >= 0) && (i < rootVec.size())) {
					//cout << rootVec[i] << "\t" << rootVec[i-1] << endl;
					bool uniq = (subset(abs(rootVec[i] - rootVec[i-1]), intervalw(0, 1e-10)));
					if ( uniq ) { 
						//cout << "this root has a duplicate" << endl;
						uniqueRootVec.pop_back(); }
				}
			}
			//cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;
			
			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				//cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
					//cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = getL1error(fhat, thisInt, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = getL1error(fhat, thisSubIntLeft, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = getL1error(fhat, thisSubIntRight, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				//cout << "let's have a look at all the roots:" << endl;
				//for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					//cout << uniqueRootVec[i] << endl;
				//}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt, deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast, deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = getL1error(fhat, thisSubIntFirst, deg, tol, mixt.W, mixt.M, mixt.S);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = getL1error(fhat, thisSubInt,  deg, tol, mixt.W, mixt.M, mixt.S);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = getL1error(fhat, thisSubIntLast,  deg, tol, mixt.W, mixt.M, mixt.S);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty

	} // end of iterating through the leaf nodes
	
	//cout << "IAE: " << totalArea << endl;
	return totalArea;
}

/*! Get the IAE for a laplace distribution for regular histograms
*/
interval getRegHistLaplaceIntervalIAE(size_t n, RegHist & myRegHist, double tol, int deg)
{
	interval totalArea(0.0);
	size_t Nbin = myRegHist.heights.size();
	
	for (size_t j=0; j< Nbin; j++){
		double fhat = myRegHist.heights[j];
		double xupp = _double(myRegHist.UpperBoxes[j]);
		double xlow = _double(myRegHist.LowerBoxes[j]);
		
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;

		//---------find the root at this domain
		// make an intervalw object using thisBox
		intervalw thisIntW(xlow, xupp);
		interval thisInt(xlow, xupp);

		// find the root
		//cout << "finding roots at this node " << thisInt << endl;
		LaplaceBisect(thisIntW, tol, fhat, rootVec); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = LaplaceGetL1error(fhat, thisInt, deg, tol);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				//cout << "root " << i << ": " << rootVec[i] << endl;
				//first insert into uniqueRootVec
				uniqueRootVec.push_back(rootVec[i]);
				//cout << i-1 << "\t" << i << "\t" << i+1 << endl;
				//now check for uniqueness
				if (((i-1) >= 0) && (i < rootVec.size())) {
					//cout << rootVec[i] << "\t" << rootVec[i-1] << endl;
					bool uniq = (subset(abs(rootVec[i] - rootVec[i-1]), intervalw(0, 1e-10)));
					if ( uniq ) { 
						//cout << "this root has a duplicate" << endl;
						uniqueRootVec.pop_back(); }
				}
			}
			//cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;
			
			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				//cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
					//cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = LaplaceGetL1error(fhat, thisInt, deg, tol);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntLeft, deg, tol);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = LaplaceGetL1error(fhat, thisSubIntRight, deg, tol);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				//cout << "let's have a look at all the roots:" << endl;
				//for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					//cout << uniqueRootVec[i] << endl;
				//}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = LaplaceGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LaplaceGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LaplaceGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty

	} // end of iterating through the leaf nodes
	
	//cout << "IAE: " << totalArea << endl;
	return totalArea;
}

/*! Get the IAE for a lognormal distribution for regular histograms
*/
interval getRegHistLognormalIntervalIAE(size_t n, RegHist & myRegHist, double tol, int deg)
{
	interval totalArea(0.0);
	size_t Nbin = myRegHist.heights.size();
	
	for (size_t j=0; j< Nbin; j++){
		double fhat = myRegHist.heights[j];
		double xupp = _double(myRegHist.UpperBoxes[j]);
		double xlow = _double(myRegHist.LowerBoxes[j]);
		
		//a container for the roots at this leaf node
		vector<intervalw> rootVec;

		//---------find the root at this domain
		// make an intervalw object using thisBox
		intervalw thisIntW(xlow, xupp);
		interval thisInt(xlow, xupp);

		// find the root
		//cout << "finding roots at this node " << thisInt << endl;
		LognormalBisect(thisIntW, tol, fhat, rootVec); 

		//---------find the area at this domain and take the absolute value
		//if rootVec is empty, there are no roots - so we can integrate over
		//this domain
		if ((rootVec.size() == 0)) { 
			//cout << "no roots at " << thisInt << endl;
			//get the L1 error
			interval diffArea = LognormalGetL1error(fhat, thisInt, deg, tol);
			//add to totalArea
			totalArea += diffArea;
		} //end of rootVec is empty

		else { //if rootVec is not empty
			vector<intervalw> uniqueRootVec;
			// make the elements in vector unique
			for (int i = 0; i < (rootVec.size()); i++) {
				//cout << "root " << i << ": " << rootVec[i] << endl;
				//first insert into uniqueRootVec
				uniqueRootVec.push_back(rootVec[i]);
			//cout << i-1 << "\t" << i << "\t" << i+1 << endl;
				//now check for uniqueness
				if (((i-1) >= 0) && (i < rootVec.size())) {
					//cout << rootVec[i] << "\t" << rootVec[i-1] << endl;
					bool uniq = (subset(abs(rootVec[i] - rootVec[i-1]), intervalw(0, 1e-10)));
					if ( uniq ) { 
						//cout << "this root has a duplicate" << endl;
						uniqueRootVec.pop_back(); }
				}
			}
			//cout << "==There are " << uniqueRootVec.size() << " unique root(s)==" << endl;
			
			// if there's only 1 root
			if (uniqueRootVec.size() == 1) {
				//cout << "there is only one root.." << endl;
				// is the root at the left or right boundary?
				if ( (abs(Inf(thisInt) - inf(rootVec[0])) < 1e-10) || 
					  (abs(Sup(thisInt) - inf(rootVec[0])) < 1e-10) ) {
					//cout << "there's a root at the left/right boundary:" << rootVec[0] << endl;
					interval diffArea = LognormalGetL1error(fhat, thisInt, deg, tol);
					totalArea += diffArea;
				}
				else { // the root is not at the boundaries
					//cout << "no root at the boundaries" << endl;
					//get the left sub-interval
					interval thisSubIntLeft = interval(Inf(thisInt), sup(uniqueRootVec[0]));
					//cout << "left interval: " << thisSubIntLeft << endl; 
					interval diffArea = LognormalGetL1error(fhat, thisSubIntLeft, deg, tol);
					totalArea += diffArea;
					
					//get the right sub-interval
					//get the left sub-interval
					interval thisSubIntRight = interval(inf(uniqueRootVec[0]), Sup(thisInt));
					//cout << "right interval: " << thisSubIntRight << endl; 
					diffArea = LognormalGetL1error(fhat, thisSubIntRight, deg, tol);
					totalArea += diffArea;
				}
			} // end of rootVec.size() == 1

				// if there is more than 1 root
			else {
				//cout << "let's have a look at all the roots:" << endl;
				//for (size_t i = 0; i < uniqueRootVec.size(); i++) {
					//cout << uniqueRootVec[i] << endl;
				//}

				//first check if the first root is at the boundary
				//cout << "check boundaries: " << Inf(thisInt) << "\t" << inf(rootVec[0]) << endl;
				if ( abs(Inf(thisInt) - inf(uniqueRootVec[0])) < 1e-10 ) {
					//cout << "there's a root at the leftmost boundary:" << endl;
					interval thisSubIntFirst = interval(Inf(thisInt), sup(uniqueRootVec[1]));
					//cout << "0-th interval:" << thisSubIntFirst << endl; 
					interval diffArea = LognormalGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						//cout << "the " << i+1 << "-th root is: " << rootVec[i+1] << endl;
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << i << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (excep the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea =LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of if first root is the boundary
				
				else {
					//cout << "root not at boundary" << endl;
					//if it is not the boundary, make the first sub-interval
					interval thisSubIntFirst = interval(Inf(thisInt), sup(rootVec[0]));
					//cout << "0-th interval: " << thisSubIntFirst << endl; 
					interval diffArea = LognormalGetL1error(fhat, thisSubIntFirst, deg, tol);
					totalArea += diffArea;
					
					// now iterate through each root (except the first and last) and 
					// get the sub-itnervals
					//cout << "iterating through each root" << endl;
					for (size_t i = 0; i < (uniqueRootVec.size() - 1); i++) {
						if ( (i+1) > uniqueRootVec.size() ) { // already no more roots
							//cout << inf(rootVec[i]) << "\t" << Sup(thisInt) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), Sup(thisInt));
							//cout << "the " << i << "-th interval: " << thisSubInt << endl;
							interval diffArea =LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
						else { //there are still more roots
							//cout << inf(rootVec[i]) << "\t" << sup(rootVec[i+1]) << endl;
							interval thisSubInt = interval(inf(uniqueRootVec[i]), sup(uniqueRootVec[i+1]));
							//cout << "the " << i+1 << "-th interval: " << thisSubInt << endl;
							interval diffArea = LognormalGetL1error(fhat, thisSubInt, deg, tol);
							totalArea += diffArea;
						}
					} // end of iterate through each root (except the first and last)
					
					// now check if the last root is at the boundary
					if ( abs(Sup(thisInt) - sup(uniqueRootVec[uniqueRootVec.size()-1])) < 1e-10 ) {
						//cout << "there's a root at the rightmost boundary:" << endl;
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-2]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea =LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					}
					else { //the last root is not at the boundary
						interval thisSubIntLast = interval(inf(uniqueRootVec[uniqueRootVec.size()-1]), Sup(thisInt));
						//cout << "last interval: " << thisSubIntLast << endl; 
						interval diffArea = LognormalGetL1error(fhat, thisSubIntLast, deg, tol);
						totalArea += diffArea;
					} 
				} // end of first root is not the boundary
			} // end of rootVec.size() > 1
		} // end of rootVec is not empty

	} // end of iterating through the leaf nodes
	
	//cout << "IAE: " << totalArea << endl;
	return totalArea;
}


/*! Get the IAE for a uniform mixture distribution for regular histograms 
*/
real getRegHistUnifIAE(RegHist & myRegHist, 
                       AdaptiveHistogram & myPart, size_t n, double weight, 
                       std::vector<int> holesLoc)
{
	//-------setting up containers-------------------------------
   dotprecision dpIAE;    // use type dotprecision for summation  
   dpIAE=0.0;

   // get the true height, rueF of the corresponding box in myPart
	SPSnodePtrs trueLeaves;
	SPSnodePtrsItr trueIt;
	(myPart).getSubPaving()->getLeaves(trueLeaves);

	ivector temp; //for the intersections
	
   //go through all the leaves in the regular histogram
	int Nbin=(myRegHist.UpperBoxes).size();

   
   for(int i = 0; i < Nbin; i++) {
		// get the height of this leaf
      double fhat = myRegHist.heights[i];
		
		//make this box into an ivector
		interval iBox(myRegHist.LowerBoxes[i], myRegHist.UpperBoxes[i]);
		ivector thisBox(1);
		thisBox[1] = iBox;
		
		//cout << "====checking " << thisBox << endl;
      //cout << "fhat for box " << ":" << fhat << endl;

		size_t L = 0;
		for (trueIt = trueLeaves.begin(); trueIt < trueLeaves.end(); trueIt++) {
			//cout << "----True leaf: " << (*trueIt)->getBox() << "\t" << endl;
			ivector trueBox = (*trueIt)->getBox();

			double trueF;
			if (  holesLoc[L] == 0 ) { trueF = 0; }
			else { trueF = weight/((*trueIt)->nodeVolume()); }
			//cout << "pdf: " << trueF << "------" << endl;
			
			// if this is contained in trueBox
			if ( thisBox <= trueBox || thisBox == trueBox ) {
				//use the volume of this
				real r = Volume(thisBox)*(fhat - trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
				//can move on to next leaf rather than iterating thru all trueBoxes
				//think about this later
			} //end of if this box is in trueBox
			
			// if this contains trueBox
			else if ( trueBox <= thisBox) {
				//use the volume of trueBox
				real r = ((*trueIt)->nodeVolume())*(fhat - trueF);
				//cout << "r: " << r << "\t" << abs(r) << endl;
				accumulate(dpIAE, abs(r), 1.0);
			} //end of if trueBox is in this box
			
			// if this is partially contained in trueBox 
			else if 	(Intersection(temp, thisBox, trueBox)) {
				if (Inf(temp) != Sup(temp)){
					double volume = Volume(temp);
					real r = volume*(fhat - trueF);
					//cout << "r: " << r << "\t" << abs(r) << endl;
					accumulate(dpIAE, abs(r), 1.0);
				}
			}
			L++;
		} // end of going through trueBoxes
	} // end of going through thisBoxes
	
		
		
   //cast dotprecision to real
   real unifIAE = rnd(dpIAE);
	return unifIAE;	               
}

//output regular histogram to .txt
void outputRegHistToTxt(RegHist& myRegHist, std::string& s)
{
	int prec = 5; // precision for output
	ofstream os(s.c_str()); 
	//output the nodeVolume, height, boxes
	for (size_t i = 0; i < myRegHist.UpperBoxes.size(); i++) {
		os << myRegHist.heights[i];
		streamsize oldPrec = os.precision();
		os << setprecision(prec);
		// intervals of box using Inf and Sup
		os << "\t" << myRegHist.LowerBoxes[i];
		os << "\t" << myRegHist.UpperBoxes[i];
		os << setprecision(oldPrec);
		os << endl;
	}
}
//---------------End of family of function for regular histograms------------//

//----------functions for finding roots----------------------------------//
//return the interval function?
intervalw F (const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S) 
{ return value(f(x, W, M, S)); }

//return the derivative of the function?
intervalw DF(const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S) 
{ return deriv(f(x, W, M, S)); }


//newton's routine
// need to modify this to suit fhat?
intervalw N (const intervalw &x, vector<double>& W, 
					vector<double>& M, vector<double>& S, double fhat)
{
  intervalw midX(mid(x));
  return midX - F(midX, W, M, S)/DF(x, W, M, S);
}

//find the root using interval newton method
void findRoot(const intervalw &domain, vector<double>& W, vector<double>& M, 
					vector<double>& S, double fhat, vector<intervalw>& rootVec)
{
  //cout << "finding root at " << domain << " for " << fhat << endl;
  
  intervalw newX       = domain;
  intervalw oldX       = domain + 1;
  bool     rootUnique = false;
  bool     rootExists = true;

  while( (newX != oldX) && rootExists ) {
    oldX = newX;
    if( !intersect(newX, N(oldX, W, M, S, fhat), oldX) ) 
      rootExists = false;
    if ( subset(newX, oldX) ) 
      rootUnique = true;
  }
  if ( rootExists ) {
    cout << newX;
    rootVec.push_back(newX);
    if ( rootUnique )
      cout << " contains a unique root." << endl;
    else
      cout << " may contain a simple root." << endl;
  } 
}

//bisect the domain and decide which root-finding routine to use
void bisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootVec, vector<double>& W, 
				vector<double>& M, vector<double>& S)
{
	//cout << "root finding routine: " << endl;
	//cout << x << "\t" << F(x, W, M, S)  << "\t" << DF(x, W, M, S) << endl;

	//if the function is twice differentiable?
	if ( !subset(0.0, DF(x, W, M, S)) ) {
		//cout << "Sending " << x << " to the Newton operator..." << endl;
		findRoot(x, W, M, S, fhat, rootVec);
	}
	
	else {
	 // if the function is not differentiable
		if ( subset(fhat, F(x, W, M, S)) ) {
			//cout << diam(x) << "\t" << TOL << endl;
			if ( diam(x) < TOL ) {
				//cout << x << " may contain roots. " << endl;
				rootVec.push_back(x); //keep the roots in a container
			}
			else {
				bisect(intervalw(inf(x), mid(x)), TOL, fhat, rootVec, W, M, S);
				bisect(intervalw(mid(x), sup(x)), TOL, fhat, rootVec, W, M, S);
			}
		}
	}
}

//----------functions for integration----------------------------------
interval riemannTerm(pfcn f, interval X, int Deg, vector<double>& W, 
				vector<double>& M, vector<double>& S) {
  
  interval Mid = interval(mid(X));

  // Taylor series...
  // cout << "taylor series:" << endl;
  itaylor fx  = f(itaylor::variable(Mid, Deg), W, M, S);
  interval sum = fx[0]*(diam(X))/2;
  for (int k = 2; k <= Deg; k += 2) 
   { sum += fx[k]*pow((diam(X))/2, k + 1)/(k + 1); }

  // Remainder term...
  // cout << "remainder term: " << endl;
  itaylor Fx  = f(itaylor::variable(X, Deg), W, M, S);
  //cout << "Fx: " << Fx << endl;
  interval eps = abs(Fx[Deg] - fx[Deg]); 
  //cout << "eps: " << eps << endl;
  sum += interval(-1.0*Sup(eps), Sup(eps))*pow((diam(X))/2, Deg + 1)/(Deg + 1); 

  return 2*sum;
}

interval integrate(pfcn f, interval X, int Deg, double Tol, vector<double>& W, 
				vector<double>& M, vector<double>& S) 
{
  //cout << "-------integrating over the domain -------------" << X << endl;
  //cout << "get  riemann term of " << X << endl;
  interval sum = riemannTerm(f, X, Deg, W, M, S);
  //cout << "sum: " << sum << "\t diam(sum):" << diam(sum) << "\t Tol: " << Tol <<  endl;
  
  if ( diam(sum) <= Tol ) {
    //cout << "Domain: " << X << "\t Sum: " << sum << endl;
    return sum;
	}
  else {  
    //cout << "*****diam(sum) > tol: inf(x), mid(x) + mid(x), sup(x)*****" << endl;
    return integrate(f, interval(Inf(X), mid(X)), Deg, Tol/2, W, M, S) + \
				integrate(f, interval(mid(X), Sup(X)), Deg, Tol/2, W, M, S);
  }
}

interval getL1error(double fhat, interval& thisInt, int Deg, double TOL, 
							vector<double>& W, vector<double>& M, vector<double>& S)
{
	//cout << "==========get L1 error for " << thisInt <<  endl;
	//hard-code this temporarily
	double Tol = 0.0000001;
	//cout << Tol << endl;
	
	//get the area of the histogram at this interval
	real histArea = diam(thisInt) * fhat;

	//integrate the function at this sub-interval
	interval fArea = integrate(integrand, thisInt, Deg, Tol, W, M, S);

	//cout << "get the differences " << endl;
	//get the differences of the areas
	interval diffArea = abs(abs(fArea) - histArea);
	//cout << "fArea: " << fArea << "\t" << "\t histArea: " << histArea << endl;
	//cout << diffArea << endl;

	return diffArea;
}

//----the true density for gaussian mixtures-------//
//for root finding routine
ia_ad f(const intervalw &x, vector<double>& W, vector<double>& M, vector<double>& S)
{
	ia_ad PDF(ia_ad::variable(x));

	intervalw startPDF(0,0);
	PDF = startPDF;
	
	size_t Ncomp = W.size();
	for (size_t c=0; c < Ncomp; c++){
		//cout << c << "-th component:" << W[c] << "\t" << M[c] << "\t" << S[c] << endl;
		intervalw z = pow((x-M[c])/S[c], 2);
		intervalw expPart = intervalw(exp(-0.5*sup(z)), exp(-0.5*inf(z)));
		PDF = PDF + expPart*(W[c]/(S[c]*sqrt(2*M_PI)));
	}
	//cout << "Domain: " << x << "\t PDF: " << PDF << endl;

	return PDF;
}

//for integration routine
typedef itaylor (*pfcn)(const itaylor &, vector<double>&, vector<double>&, vector<double>&);

itaylor integrand(const itaylor &x, vector<double>& W, vector<double>& M, vector<double>& S) 
{ 
	
	itaylor PDF(itaylor::variable(interval(0,0), orderOf(x)));
	//cout << "-------initializing variable PDF: ---------" << PDF << endl;

	//cout << "x: " << x[0] << endl;
	size_t Ncomp = W.size();
	for (size_t c=0; c < Ncomp; c++){
		//cout << c << "-th component: " << W[c] << "\t" << M[c] << "\t" << S[c] << endl;

		itaylor z1 = (x-M[c])/S[c]; 
		//when inf is negative and sup is positive - split into two parts 
		if (Inf(z1[0]) < 0 && Sup(z1[0]) > 0) {  
			//cout << "===========" << endl;
			//cout << z1[0] << endl;
			itaylor z = z1*z1;
			//cout << z[0] << endl;
			Inf(z[0]) = 0;
			//cout << z[0] << endl;
			PDF = PDF + exp(-0.5*z)*(W[c]/(S[c]*sqrt(2*M_PI)));
		}  
		else {
			itaylor z = z1*z1;
			//cout << "z: " << z[0] << endl;
			itaylor expPart = exp(-0.5*z);
			//cout << "exp part:" << expPart[0] << endl;
			PDF = PDF + expPart*(W[c]/(S[c]*sqrt(2*M_PI)));
		}
	}
	//cout << "Domain: " << x[0] << "\t PDF: " << PDF[0] << endl;
	
	return PDF;
}

//========functions for laplace distribution================//
//The Laplace Distribution
//for root finding routine
ia_ad LaplacePDF(const intervalw &x)
{
	ia_ad PDF(ia_ad::variable(x));

	//strictly negative cases
	if ( sup(x) < 0.0 ) {
		PDF = 0.5*intervalw(exp(inf(x)), exp(sup(x)));
	}
	//strictly positive
	else if ( inf(x) >= 0.0 ) {
		PDF = 0.5*intervalw(exp(-sup(x)), exp(-inf(x)));
	}
	// 0 is inside the interval x
	else {
		 // this is a symmetric density
		 if ( fabs(inf(x)) >= sup(x) ) {
			PDF = 0.5*intervalw(exp(inf(x)), exp(0));
		 }
		 else { 
			 PDF = 0.5*intervalw(exp(-sup(x)), exp(0));
		 }
	}
	
	//cout << "Domain: " << x << "PDF: " << PDF << endl;

	return PDF;
}

//for integration routine
typedef itaylor (*pfcnLaplace)(const itaylor &);

itaylor LaplaceIntegrand(const itaylor &x) 
{ 
	//cout << "-------Integrand: ---------" << x[0] << "\t";
	itaylor PDF(itaylor::variable(interval(0,0), orderOf(x)));

	if ( Sup(x[0]) < 0.0 ) {
		PDF = 0.5*exp(-(-x));
	}
	else if ( Inf(x[0]) >= 0.0 ) { 
		PDF = 0.5*exp(-(x));
	}
	// 0 is inside the interval x
	else { 
		//cout << gsl_ran_laplace_pdf(_double(Inf(x[0])), 1) << "\t" << gsl_ran_laplace_pdf(_double(Sup(x[0])), 1) << endl;
		if ( abs(Inf(x[0])) >= Sup(x[0]) ) {
			PDF = 0.5*exp(x);
		}
		else {
			PDF = 0.5*exp(-x);
		}
	}
	
	//cout << "PDF: " << PDF[0] << endl;
	return PDF;
}

//----------functions for finding roots: Laplace----------//
//return the interval function?
intervalw LaplaceF (const intervalw &x) 
{ return value(LaplacePDF(x)); }

//return the derivative of the function?
intervalw LaplaceDF(const intervalw &x) 
{ return deriv(LaplacePDF(x)); }


//newton's routine
// need to modify this to suit fhat?
intervalw LaplaceNewton(const intervalw &x, double fhat)
{
  intervalw midX(mid(x));
  return midX - LaplaceF(midX)/LaplaceDF(x);
}

//find the root using interval newton method
void LaplacefindRoot(const intervalw &domain, double fhat, vector<intervalw>& rootVec)
{
  //cout << "finding root at " << domain << " for " << fhat << endl;
 
  intervalw newX       = domain;
  intervalw oldX       = domain + 1;
  bool     rootUnique = false;
  bool     rootExists = true;

  while( (newX != oldX) && rootExists ) {
    oldX = newX;
    if( !intersect(newX, LaplaceNewton(oldX, fhat), oldX) ) 
      rootExists = false;
    if ( subset(newX, oldX) ) 
      rootUnique = true;
  }
  if ( rootExists ) {
    cout << newX;
    rootVec.push_back(newX);
    if ( rootUnique )
      cout << " contains a unique root." << endl;
    else
      cout << " may contain a simple root." << endl;
  } 
}

//bisect the domain and decide which root-finding routine to use
void LaplaceBisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootVec)
{
	//cout << "===========root finding routine at domain: " << x <<  endl;
	//cout << LaplacePDF(x) << endl;
	//cout << "gsl: " << gsl_ran_laplace_pdf(inf(x), 1) << "\t" <<  gsl_ran_laplace_pdf(sup(x), 1) << endl;
	//cout << x << "\t" << LaplaceF(x)  << "\t" << LaplaceDF(x) << endl;

	//if the function is twice differentiable?
	if ( !subset(0.0, LaplaceDF(x)) ) {
		cout << "Sending " << x << " to the Newton operator..." << endl;
		LaplacefindRoot(x, fhat, rootVec);
		cerr << "check this!" << endl;
		exit(0);
	}
	
	else {
	 // if the function is not differentiable
		//cout << "------------compare fhats" << endl;
		//cout.precision(10);
		//cout << intervalw(fhat) << "\t" << LaplaceF(x) << endl;
		//cout << ( fhat <= sup(LaplaceF(x)) ) << endl;
		
		if ( subset(intervalw(fhat), LaplaceF(x)) ) {
			//cout << "check tolerance: " << endl;
			//cout << diam(x) << "\t" << TOL << endl;
			if ( diam(x) < TOL ) {
				//cout << diam(x) << "\t" << TOL << endl;
				//cout << x << " may contain roots. " << endl;
				rootVec.push_back(x); //keep the roots in a container
			}
			else {
				//cout << "bisect: " << endl;
				//cout << "left" << endl;
				LaplaceBisect(intervalw(inf(x), mid(x)), TOL, fhat, rootVec);
				//cout << "right" << endl;
				LaplaceBisect(intervalw(mid(x), sup(x)), TOL, fhat, rootVec);
			}
		}
		//else { cout << "fhat not here. " << endl; }
	}
}

//----------functions for integration----------------------------------
interval LaplaceRiemannTerm(pfcnLaplace f, interval X, int Deg) {
  
  interval Mid = interval(mid(X));

  // Taylor series...
  // cout << "taylor series:" << endl;
  itaylor fx  = LaplaceIntegrand(itaylor::variable(Mid, Deg));
  interval sum = fx[0]*(diam(X))/2;
  for (int k = 2; k <= Deg; k += 2) 
   { sum += fx[k]*pow((diam(X))/2, k + 1)/(k + 1); }

  // Remainder term...
  // cout << "remainder term: " << endl;
  itaylor Fx  = LaplaceIntegrand(itaylor::variable(X, Deg));
  //cout << "Fx: " << Fx << endl;
  interval eps = abs(Fx[Deg] - fx[Deg]); 
  //cout << "eps: " << eps << endl;
  sum += interval(-1.0*Sup(eps), Sup(eps))*pow((diam(X))/2, Deg + 1)/(Deg + 1); 

  return 2*sum;
}

interval LaplaceIntegrate(pfcnLaplace f, interval X, int Deg, double Tol) 
{
  //cout << "-------integrating over the domain -------------" << X << endl;
  //cout << "get  riemann term of " << X << endl;
  interval sum = LaplaceRiemannTerm(f, X, Deg);
  //cout << "sum: " << sum << "\t diam(sum):" << diam(sum) << "\t Tol: " << Tol <<  endl;
  
  if ( diam(sum) <= Tol ) {
    //cout << "Domain: " << X << "\t Sum: " << sum << endl;
    return sum;
	}
  else {  
    //cout << "*****diam(sum) > tol: inf(x), mid(x) + mid(x), sup(x)*****" << endl;
    //cout << "sum: " << sum << "\t Tol: " << Tol << endl;
    return LaplaceIntegrate(f, interval(Inf(X), mid(X)), Deg, Tol/2) + \
				LaplaceIntegrate(f, interval(mid(X), Sup(X)), Deg, Tol/2);
  }
}

interval LaplaceGetL1error(double fhat, interval& thisInt, int Deg, double TOL)
{
//	cout << "==========get L1 error for " << thisInt <<  endl;
	//hard-code this temporarily
	double Tol = 0.0000001;
	//cout << Tol << endl;
	
	//get the area of the histogram at this interval
	real histArea = diam(thisInt) * fhat;

	//integrate the function at this sub-interval
	interval fArea = LaplaceIntegrate(LaplaceIntegrand, thisInt, Deg, Tol);
	//cout << "integrate: " << fArea << endl;

	//cout << "get the differences " << endl;
	//get the differences of the areas
	//double up = gsl_cdf_laplace_P(_double(Sup(thisInt)), 1);
	//double low = gsl_cdf_laplace_P(_double(Inf(thisInt)), 1);
	//cout << "gsl: " << up-low << endl;
	interval diffArea = abs(abs(fArea) - histArea);
	//cout << "fArea: " << fArea << "\t" << diam(thisInt) << 
	//"\t fhat:" << fhat << "\t histArea: " << histArea << endl;
	//cout << diffArea << endl;

	return diffArea;
}
//============end of functions for Laplace============================//

//========functions for lognormal=================================//
//for root finding routine
ia_ad LognormalPDF(const intervalw &x)
{
	ia_ad PDF(ia_ad::variable(x));
	
	PDF = 1/x*1/sqrt(2*M_PI)*exp(-0.5*log(x)*log(x));
	
	return PDF;
	
	//cout << "Domain: " << x << "PDF: " << PDF << endl;
}

//for integration routine
typedef itaylor (*pfcnLognormal)(const itaylor &);

itaylor LognormalIntegrand(const itaylor &x) 
{ 
	//cout << "-------Integrand: ---------" << endl;
	itaylor PDF(itaylor::variable(interval(0,0), orderOf(x)));
	
	PDF = 1/x*1/sqrt(2*M_PI)*exp(-0.5*log(x)*log(x));

	//cout << PDF[0] << endl;
	
	return PDF;
}

//----------functions for finding roots: Lognormal----------//
//return the interval function?
intervalw LognormalF (const intervalw &x) 
{ return value(LognormalPDF(x)); }

//return the derivative of the function?
intervalw LognormalDF(const intervalw &x) 
{ return deriv(LognormalPDF(x)); }

//newton's routine
// need to modify this to suit fhat?
intervalw LognormalNewton(const intervalw &x, double fhat)
{
  intervalw midX(mid(x));
  return midX - LognormalF(midX)/LognormalDF(x);
}

//find the root using interval newton method
void LognormalfindRoot(const intervalw &domain, double fhat, vector<intervalw>& rootVec)
{
  cout << "finding root at " << domain << " for " << fhat << endl;
 
  intervalw newX       = domain;
  intervalw oldX       = domain + 1;
  bool     rootUnique = false;
  bool     rootExists = true;

  while( (newX != oldX) && rootExists ) {
    oldX = newX;
    if( !intersect(newX, LognormalNewton(oldX, fhat), oldX) ) 
      rootExists = false;
    if ( subset(newX, oldX) ) 
      rootUnique = true;
  }
  if ( rootExists ) {
    cout << newX;
    rootVec.push_back(newX);
    if ( rootUnique )
      cout << " contains a unique root." << endl;
    else
      cout << " may contain a simple root." << endl;
  } 
}

//bisect the domain and decide which root-finding routine to use
void LognormalBisect(const intervalw &x, const double &TOL, double &fhat, 
				vector<intervalw>& rootVec)
{
	//cout << "===========root finding routine at domain: " << x <<  endl;
	//cout << x << "\t" << LaplaceF(x)  << "\t" << LognormalDF(x) << endl;

	//if the function is twice differentiable?
	if ( !subset(0.0, LognormalDF(x)) ) {
		cout << "Sending " << x << " to the Newton operator..." << endl;
		LognormalfindRoot(x, fhat, rootVec);
		cerr << "check this!" << endl;
		exit(1);
	}
	
	else {
	 // if the function is not differentiable
	//	cout << "------------compare fhats" << endl;
//		cout.precision(10);
	//	cout << intervalw(fhat) << "\t" << LaplaceF(x) << endl;
		//cout << ( fhat <= sup(LaplaceF(x)) ) << endl;
		
		if ( subset(intervalw(fhat), LognormalF(x)) ) {
		//	cout << "check tolerance: " << endl;
	//		cout << diam(x) << "\t" << TOL << endl;
			if ( diam(x) < TOL ) {
				//cout << diam(x) << "\t" << TOL << endl;
				//cout << x << " may contain roots. " << endl;
				rootVec.push_back(x); //keep the roots in a container
			}
			else {
		//		cout << "bisect: " << endl;
		//		cout << "left" << endl;
				LognormalBisect(intervalw(inf(x), mid(x)), TOL, fhat, rootVec);
			//	cout << "right" << endl;
				LognormalBisect(intervalw(mid(x), sup(x)), TOL, fhat, rootVec);
			}
		}
		//else { cout << "fhat not here. " << endl; }
	}
}

//----------functions for integration----------------------------------
interval LognormalRiemannTerm(pfcnLaplace f, interval X, int Deg) {
  
  interval Mid = interval(mid(X));

  // Taylor series...
  // cout << "taylor series:" << endl;
  itaylor fx  = LognormalIntegrand(itaylor::variable(Mid, Deg));
  interval sum = fx[0]*(diam(X))/2;
  for (int k = 2; k <= Deg; k += 2) 
   { sum += fx[k]*pow((diam(X))/2, k + 1)/(k + 1); }

  // Remainder term...
  // cout << "remainder term: " << endl;
  itaylor Fx  = LognormalIntegrand(itaylor::variable(X, Deg));
  //cout << "Fx: " << Fx << endl;
  interval eps = abs(Fx[Deg] - fx[Deg]); 
  //cout << "eps: " << eps << endl;
  sum += interval(-1.0*Sup(eps), Sup(eps))*pow((diam(X))/2, Deg + 1)/(Deg + 1); 

  return 2*sum;
}

interval LognormalIntegrate(pfcnLognormal f, interval X, int Deg, double Tol) 
{
  //cout << "-------integrating over the domain -------------" << X << endl;
  //cout << "get  riemann term of " << X << endl;
  interval sum = LognormalRiemannTerm(f, X, Deg);
  //cout << "sum: " << sum << "\t diam(sum):" << diam(sum) << "\t Tol: " << Tol <<  endl;
  
  if ( diam(sum) <= Tol ) {
    //cout << "Domain: " << X << "\t Sum: " << sum << endl;
    return sum;
	}
  else {  
    //cout << "*****diam(sum) > tol: inf(x), mid(x) + mid(x), sup(x)*****" << endl;
    return LognormalIntegrate(f, interval(Inf(X), mid(X)), Deg, Tol/2) + \
				LognormalIntegrate(f, interval(mid(X), Sup(X)), Deg, Tol/2);
  }
}

interval LognormalGetL1error(double fhat, interval& thisInt, int Deg, double TOL)
{
	//cout << "==========get L1 error for " << thisInt <<  endl;
	//cout.precision(10);
	
	//hard-code this temporarily
	double Tol = 0.0000001;
	//cout << Tol << endl;
	
	//get the area of the histogram at this interval
	real histArea = diam(thisInt) * fhat;

	//integrate the function at this sub-interval
	interval fArea = LognormalIntegrate(LognormalIntegrand, thisInt, Deg, Tol);
	//cout << "integrate: " << fArea << endl;

	//cout << "get the differences " << endl;
	//get the differences of the areas
	//double up = gsl_cdf_lognormal_P(_double(Sup(thisInt)), 0, 1);
	//double low = gsl_cdf_lognormal_P(_double(Inf(thisInt)), 0, 1);
	//cout << "gsl: " << up-low << endl;
	interval diffArea = abs(abs(fArea) - histArea);
	//cout << "fArea: " << fArea << "\t" << diam(thisInt) << 
	//"\t fhat:" << fhat << "\t histArea: " << histArea << endl;
	//cout << diffArea << endl;

	return diffArea;
}
//==========end of functions for lognormal==========================//

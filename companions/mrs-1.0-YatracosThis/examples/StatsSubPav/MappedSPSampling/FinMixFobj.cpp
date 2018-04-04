/*! \file FinMixFobj.cpp
\brief MappedSPnode example 1-d function object class.

This example is the pdf of a gaussian mixture.

*/

#include "FinMixFobj.hpp"
#include <cmath> //to use 2_PI
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

FinMixFobj::FinMixFobj()
        : W(1.0), M(0.0), S(1.0) {};

FinMixFobj::FinMixFobj(vector<double> WW, vector<double> MM, vector<double> SS)
        : W(WW), M(MM), S(SS) {};


interval FinMixFobj::operator()(const interval& ival) const
{
	interval PDF(0,0);
	size_t Ncomp = W.size();
	for (size_t c=0; c < Ncomp; c++){
		interval z = power((ival - M[c])/S[c],2);
		PDF += W[c]/(sqrt(2*M_PI)*S[c])*exp(-0.5*z);
	}  
	 return PDF;
}

real FinMixFobj::operator()(const real& r) const
{
	real PDF = 0;
	size_t Ncomp = W.size();
	for (size_t c=0; c < Ncomp; c++){
		real z = power((r-M[c])/S[c], 2);
		PDF += W[c]*exp(-0.5*z)/(S[c]*sqrt(2*M_PI));
	}
	return PDF;
}


/*! \file GaussianFobj1D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate gaussian.

*/

#include "GaussianFobj1D.hpp"
#include <cmath> //to use M_PI
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

interval GaussianFobj1D::operator()(const cxsc::interval& ival1) const
{
	real a = power(2*M_PI, 1.0/2.0);
   interval b = -0.5 * (power(ival1,2));
   interval IntPDF = 1.0/a * exp(b);

	return IntPDF;
}

real GaussianFobj1D::operator()(const cxsc::real& r1) const
{
    real a = power(2*M_PI, 1.0/2.0);
    real b = -0.5 * (power(r1,2));
    real RePDF = 1.0/a * exp(b);
    
    return RePDF;
}


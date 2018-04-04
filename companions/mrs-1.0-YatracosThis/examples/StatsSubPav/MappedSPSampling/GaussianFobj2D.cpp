/*! \file GaussianFobj2D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate gaussian.

*/

#include "GaussianFobj2D.hpp"
#include <cmath> //to use M_PI
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

interval GaussianFobj2D::operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const
{
	real a = power(2*M_PI, 2.0/2.0);
   interval b = -0.5 * (power(ival1,2) + power(ival2, 2));
   interval IntPDF = 1.0/a * exp(b);
 //cout << ival1 << "\t" << ival2 << "\t" << IntPDF << endl;
	return IntPDF;
}

real GaussianFobj2D::operator()(const cxsc::real& r1,
                                        const cxsc::real& r2) const
{
    real a = power(2*M_PI, 2.0/2.0);
    real b = -0.5 * (power(r1,2) + power(r2, 2));
    real RePDF = 1.0/a * exp(b);
 //   cout << r1 << "\t" << r2 << "\t" << RePDF << endl;
    return RePDF;
}


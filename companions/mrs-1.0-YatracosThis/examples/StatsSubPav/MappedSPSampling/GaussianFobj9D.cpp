/*! \file GaussianFobj10D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate gaussian.

*/

#include "GaussianFobj9D.hpp"
#include <cmath> //to use M_PI
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

interval GaussianFobj9D::operator()(
			const cxsc::interval& ival1,
			const cxsc::interval& ival2,
			const cxsc::interval& ival3,
			const cxsc::interval& ival4,
			const cxsc::interval& ival5,
			const cxsc::interval& ival6,
			const cxsc::interval& ival7,
			const cxsc::interval& ival8,
			const cxsc::interval& ival9) const
{
	real a = pow(2*M_PI, 5);
   interval b = -0.5 * 
   (
   power(ival1,2) + 
   power(ival2,2) +
   power(ival3,2) +
   power(ival4,2) +
   power(ival5,2) +
   power(ival6,2) +
   power(ival7,2) +
   power(ival8,2) +
   power(ival9,2)
   );
   
   interval IntPDF = 1.0/a * exp(b);

	return IntPDF;
}

real GaussianFobj9D::operator()(
				const cxsc::real& r1,
				const cxsc::real& r2,
				const cxsc::real& r3,
				const cxsc::real& r4,
				const cxsc::real& r5,
				const cxsc::real& r6,
				const cxsc::real& r7,
				const cxsc::real& r8,
				const cxsc::real& r9) const
{
    real a = pow(2*M_PI, 5);
    real b = -0.5 * (power(r1, 2)+ power(r2, 2)+ power(r3, 2)+
				power(r4, 2)+ power(r5, 2)+ power(r6, 2)+ power(r7, 2)+
				power(r8, 2)+ power(r9, 2));
    real RePDF = 1.0/a * exp(b);
    
  //  cout << r1 << r2 << r3 << r4 << r5 << r6 << r7 << r8 << r9 << r10 << RePDF << endl;
    
    return RePDF;
}

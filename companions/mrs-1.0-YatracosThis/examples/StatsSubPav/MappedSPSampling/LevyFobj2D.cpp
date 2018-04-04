/*! \file RosenFobj2D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate Rosen.

*/

#include "LevyFobj2D.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;


//Parameters specific to the Levy target
real Temperature = 40.0;
real Center1 = 1.42513; 
real Center2 = 0.80032; 
real GlobalMax = 176.14;
real DomainLimit = 10.0;

interval LevyFobj2D::operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const
{
    
	interval isum(0,0);
	interval jsum(0,0);
    
    for (int i = 1; i <= 5; i++)	{
		isum = isum + i * cos((i - 1) * ival1 + (i));
		jsum = jsum + i * cos((i + 1) * ival2 + (i));
	}
                    // Avoid real conversion error
	interval hh = isum * jsum + (ival1 + Center1)*(ival1 + Center1) + (ival2 + Center2)*(ival2 + Center2);
	hh = hh + GlobalMax;  
	interval result = exp (-hh / Temperature);
	return result;

}

real LevyFobj2D::operator()(const cxsc::real& r1,
                                        const cxsc::real& r2) const
{
	real isum=0;
	real jsum = 0;
    
    for (int i = 1; i <= 5; i++)	{
		isum = isum + i * cos((i - 1) * r1 + (i));
		jsum = jsum + i * cos((i + 1) * r2 + (i));
	}
                    // Avoid real conversion error
	real hh = isum * jsum + (r1 + Center1)*(r1 + Center1) + (r2 + Center2)*(r2 + Center2);
	hh = hh + GlobalMax;  
	real result = exp (-hh / Temperature);
	return result;
}


  

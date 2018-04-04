/*! \file RosenFobj2D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate Rosen.

*/

#include "RosenFobj2D.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

real Tinverse = 1.0;
real Height = 100.0;

interval RosenFobj2D::operator()(const cxsc::interval& ival1,
                              const cxsc::interval& ival2) const
{
//	cout << "=======int=========" << endl;

	ivector ival(2);
	ival[1] = ival1;
	ival[2] = ival2;

//	cout << ival << endl;

	interval result(0.0, 0.0);
	
	int a = Lb(ival), z = Ub(ival);
	
	for (size_t i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (ival[i] - sqr (ival[i - 1])) +
        sqr (ival[i - 1] - 1.0));
    }

//	cout << result << endl;

	result = exp (-(Tinverse * result));
	
//	cout << "Result: " << result << endl ;
//	cout << "=======int=========" << endl;
	
	return result;
	
}


real RosenFobj2D::operator()(const cxsc::real& r1,const cxsc::real& r2) const
{   
	real result = 0.0;
   
   rvector r(2);
   r[1] = r1;
   r[2] = r2;
   
  // cout << "=======real=========" << endl;
  // cout << r << endl;
   
    int a = Lb(r), z = Ub(r);
   
    for (size_t i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (r[i] - sqr (r[i - 1])) +
        sqr (r[i - 1] - 1.0));
    }
   
   //cout << result << endl;
   
	result = exp (-(Tinverse * result));

//	cout << "Result: " << result << endl; 
//	cout << "=======real=========" << endl;

   return result;
}

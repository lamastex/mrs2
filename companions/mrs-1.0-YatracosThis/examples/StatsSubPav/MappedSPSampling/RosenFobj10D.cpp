/*! \file RosenFobj10D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate Rosen.

*/

#include "RosenFobj10D.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

interval RosenFobj10D::operator()(
const cxsc::interval& ival1,
const cxsc::interval& ival2,
const cxsc::interval& ival3,
const cxsc::interval& ival4,
const cxsc::interval& ival5,
const cxsc::interval& ival6,
const cxsc::interval& ival7,
const cxsc::interval& ival8,
const cxsc::interval& ival9,
const cxsc::interval& ival10
) const
{

	real Tinverse = 1.0;
	real Height = 100.0;
  
	ivector RR(10);
	RR[1] = ival1;
	RR[2] = ival2;
	RR[3] = ival3;
	RR[4] = ival4;
	RR[5] = ival5;
	RR[6] = ival6;
	RR[7] = ival7;
	RR[8] = ival8;
	RR[9] = ival9;
	RR[10] = ival10;
	
	int a = Lb(RR), z = Ub(RR);

	//cout << "============\n" << RR << endl;


	interval result(0.0);
	
    for (size_t i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (RR[i] - sqr (RR[i - 1])) +
        sqr (RR[i - 1] - 1.0));
    }
 
	//cout << result << endl;
	result = exp (-(Tinverse * result));
	//cout << "Result: " << result << endl  << endl;

	return result;
}

real RosenFobj10D::operator()(
const cxsc::real& r1,
const cxsc::real& r2,
const cxsc::real& r3,
const cxsc::real& r4,
const cxsc::real& r5,
const cxsc::real& r6,
const cxsc::real& r7,
const cxsc::real& r8,
const cxsc::real& r9,
const cxsc::real& r10) const
{   
   real Tinverse = 1.0;
	real Height = 100.0;

	real result = 0.0;
	
	//mrs
	rvector RR(10);
	RR[1] = r1;
	RR[2] = r2;
	RR[3] = r3;
	RR[4] = r4;
	RR[5] = r5;
	RR[6] = r6;
	RR[7] = r7;
	RR[8] = r8;
	RR[9] = r9;
	RR[10] = r10;
	
	cout << RR << endl;
	
	int a = Lb(RR), z = Ub(RR);
	
    for (size_t i = a + 1; i <= z; i++)
    {
      result = result + (Height * sqr (RR[i] - sqr (RR[i - 1])) +
        sqr (RR[i - 1] - 1.0));
    }
 
	cout << result << endl;
	result = exp (-(Tinverse * result));
	cout << "Result: " << result << endl << endl;
	cout << "=============" << endl;
	
	return result;
}

/*! \file
\brief Definitions for Rosenbrock function..

*/

#include "RosenFobj.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
//using namespace subpavings;

//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif


RosenFobj::RosenFobj() : tInverse(1.0), height(100.0) {}

RosenFobj::RosenFobj(real ti, real h) : tInverse(ti), height(h) {}

interval RosenFobj::operator()(const cxsc::ivector& ivec) const
{
	interval result(0.0, 0.0);
	
	int a = Lb(ivec), z = Ub(ivec);
	
	if (a == z)
	{
		result = (height * sqr (1.0 - sqr (ivec[a])) + sqr (ivec[a] - 1.0));
	}
   
	else {
    
		for (size_t i = a + 1; i <= z; ++i) {
		  result = result + (height * sqr (ivec[i] - sqr (ivec[i - 1])) +
			sqr (ivec[i - 1] - 1.0));
		}
	}
	
	#ifdef MYDEBUG
		cout << "Rosen(ivector ) result is: " << result << endl << endl;
	#endif
	
	
	return result;
	
}


real RosenFobj::operator()(const cxsc::rvector& r) const
{   
	real result = 0.0;
   
    int a = Lb(r), z = Ub(r);
	
	if (a == z)
	{
		result = (height * sqr (1.0 - sqr (r[a])) + sqr (r[a] - 1.0));
	}
   
	else {
		for (size_t i = a + 1; i <= z; i++) {
		  result = result + (height * sqr (r[i] - sqr (r[i - 1])) +
			sqr (r[i - 1] - 1.0));
		}
	}
	
	#ifdef MYDEBUG
		cout << "Rosen(rvector ) result is: " << result << endl << endl; 
	#endif
	
	// turns the Rosenbrock function into a density
	#if(0)
		result = exp (-(tInverse * result));

		#ifdef MYDEBUG
			cout << "Result after exp is: " << result << endl << endl; 
		#endif
	#endif

   return result;
}

std::string RosenFobj::getName() const
{
	return std::string("RosenFunction");
}

RosenFobj::~RosenFobj(){}

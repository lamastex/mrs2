/*! \file
\brief MappedSPnode example Rosenbrock example.

*/

#include "RosenDensityFobj.hpp"



using namespace cxsc;
using namespace std;


//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif


RosenDensityFobj::RosenDensityFobj() : tInverse(1.0), height(100.0) {}

RosenDensityFobj::RosenDensityFobj(real ti, real h) : tInverse(ti), height(h) {}

interval RosenDensityFobj::operator()(const cxsc::ivector& ivec) const
{
	interval result(0.0, 0.0);
	
	int a = Lb(ivec), z = Ub(ivec);
	
	if (a == z)
	{
		result = (height * cxsc::sqr (1.0 - cxsc::sqr (ivec[a])) + cxsc::sqr (ivec[a] - 1.0));
	}
   
	else {
    
		for (size_t i = a + 1; i <= z; ++i) {
		  result = result + (height * cxsc::sqr (ivec[i] - cxsc::sqr (ivec[i - 1])) +
								cxsc::sqr (ivec[i - 1] - 1.0));
		}
	}
	
	#ifdef MYDEBUG
		cout << "Rosen(ivector ) result is: " << result << endl << endl;
	#endif
	
	result = cxsc::exp (-(tInverse * result));
		
	#ifdef MYDEBUG
		cout << "Result afer exp is: " << result << endl << endl;
	#endif

	
	return result;
	
}


real RosenDensityFobj::operator()(const cxsc::rvector& r) const
{   
	real result = 0.0;
   
    int a = Lb(r), z = Ub(r);
	
	if (a == z)
	{
		result = (height * cxsc::sqr (1.0 - cxsc::sqr (r[a])) + cxsc::sqr (r[a] - 1.0));
	}
   
	else {
		for (size_t i = a + 1; i <= z; i++) {
		  result = result + (height * cxsc::sqr (r[i] - cxsc::sqr (r[i - 1])) +
							cxsc::sqr (r[i - 1] - 1.0));
		}
	}
	
	#ifdef MYDEBUG
		cout << "Rosen(rvector ) result is: " << result << endl << endl; 
	#endif
	
	result = cxsc::exp (-(tInverse * result));

	#ifdef MYDEBUG
		cout << "Result after exp is: " << result << endl << endl; 
	#endif


   return result;
}

std::string RosenDensityFobj::getName() const
{
	return std::string("RosenDensity");
}

RosenDensityFobj::~RosenDensityFobj(){}

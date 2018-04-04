/*! \file
\brief MappedSPnode example Levy example.

*/

#include "LevyDensityFobj2D.hpp"
#include "cxsc.hpp"
#include "ivector.hpp"
#include <stdexcept>

using namespace cxsc;
using namespace std;

//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif


LevyDensityFobj2D::LevyDensityFobj2D() : 
		temperature(40.0), center1(1.42513), center2(0.80032), 
		globalMax(176.14) {}

//LevyDensityFobj2D::LevyDensityFobj2D(real ti, real h) : tInverse(ti), height(h) {}

interval LevyDensityFobj2D::operator()(const cxsc::ivector& ivec) const
{
	if(VecLen(ivec) != 2) {
		throw std::invalid_argument("LevyDensityFobj2D::operator()(const cxsc::ivector& ivec): ivec not 2 dimensions");
	}
	
	interval isum(0,0);
	interval jsum(0,0);
    
    for (int i = 1; i <= 5; i++)	{
		interval ii(i,i);
		isum = isum + ii * cxsc::cos((ii - 1) * ivec[1] + ii);
		jsum = jsum + ii * cxsc::cos((ii + 1) * ivec[2] + ii);
	}
	// Avoid real conversion error
	interval int_pow = cxsc::interval(2.0,2.0);
	interval hh = isum * jsum + cxsc::pow((ivec[1] + center1),int_pow) + pow((ivec[2] + center2),int_pow);
	hh = hh + globalMax;  
	interval result = exp (-hh / temperature);
	return result;
	
	
	
}


real LevyDensityFobj2D::operator()(const cxsc::rvector& r) const
{   
	if(VecLen(r) != 2) {
		throw std::invalid_argument("LevyDensityFobj2D::operator()(const cxsc::rvector& r): r not 2 dimensions");
	}
	
	real isum(0.0);
	real jsum(0.0);
    
    for (int i = 1; i <= 5; i++)	{
		isum = isum + i * cxsc::cos((i - 1) * r[1] + i);
		jsum = jsum + i * cxsc::cos((i + 1) * r[2] + i);
	}
	// Avoid real conversion error
	real hh = isum * jsum + cxsc::pow((r[1] + center1),2.0) + pow((r[2] + center2),2.0);
	hh = hh + globalMax;  
	real result = exp (-hh / temperature);
	return result;
	
}

std::string LevyDensityFobj2D::getName() const
{
	return std::string("LevyDensity");
}

LevyDensityFobj2D::~LevyDensityFobj2D(){}

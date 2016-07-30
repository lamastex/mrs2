
/*! \file
\brief Definitions for multivariate sphere  
* function object class.
* 
* f(x) = 
*/

#include "SphereFobj.hpp"
#include "cxsc.hpp"
#include "toolz.hpp"

#include <stdexcept>

using namespace cxsc;
using namespace std;
using namespace subpavings;


SphereFobj::SphereFobj() :name("Sphere"), cLen(0) {}

SphereFobj::SphereFobj(const cxsc::rvector& c) :name("sphere"),
				centre(c), cLen(VecLen(c)) {}

cxsc::interval SphereFobj::operator()(const cxsc::ivector& ivec) const
{
	int lb = Lb(ivec);
	int ub = Ub(ivec);
	
	if (cLen && (ub-lb+1 != cLen))
			throw std::runtime_error(
				"SphereFobj::operator() : dimensions incompatible");
	
	cxsc::interval result;
	
	if (cLen) {
	
		int clb = Lb(centre);
		result = cxsc::sqr(ivec[lb] - centre[clb]);
		
		for (int i = 1; i <= ub-lb; ++i) 
					result += cxsc::sqr(ivec[lb+i] - centre[clb+i]);
			
	}
	else {
	
		result = cxsc::sqr(ivec[lb]);
		
		for (int i = 1; i <= ub-lb; ++i) {
		
			result += cxsc::sqr(ivec[lb+i]);
				
		}
			
	}
	
	return cxsc::sqrt(result);
}

cxsc::real SphereFobj::operator()(const cxsc::rvector& r) const
{
	int lb = Lb(r);
	int ub = Ub(r);
	
	if (cLen && (ub-lb+1 != cLen))
			throw std::runtime_error(
				"SphereFobj::operator() : dimensions incompatible");

	cxsc::real result;
	
	if (cLen) {
	
		int clb = Lb(centre);
		result = cxsc::sqr(r[lb] - centre[clb]);
			
		for (int i = 1; i <= ub-lb; ++i) 
					result += cxsc::sqr(r[lb+i] - centre[clb+i]);
			
	}
	else {
	
		result = cxsc::sqr(r[lb]);
			
		for (int i = 1; i <= ub-lb; ++i) 
					result += cxsc::sqr(r[lb+i]);
			
	}
	
	return cxsc::sqrt(result);
}


std::string SphereFobj::getName() const
{
	return name;
}

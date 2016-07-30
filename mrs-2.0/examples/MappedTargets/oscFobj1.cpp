
/*! \file
\brief Definitions for multivariate oscillating  
* function object class.
* 
* f(x) defined for 0 <= x <= 1
*/

#include "oscFobj1.hpp"
#include "cxsc.hpp"
#include "toolz.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;


OscFobj::OscFobj() :name("Osc"), a(2.0), b(10.0), c(3.0) {}

cxsc::interval OscFobj::operator()(const cxsc::ivector& ivec) const
{
	int lb = Lb(ivec);
	int ub = Ub(ivec);
	
	cxsc::interval result = this->operator()(ivec[lb]);
	
		
	for (int i = lb+1; i <= ub; ++i) {
		
		
		result *= this->operator()(ivec[i]);
		
	}
	
	return result;
}

cxsc::real OscFobj::operator()(const cxsc::rvector& r) const
{
	int lb = Lb(r);
	int ub = Ub(r);
	
	cxsc::real retr = this->operator()(r[lb]);
	
	for (int i = lb+1; i <= ub && retr > 0.0; ++i) {
		retr *= (this->operator()(r[i]));
	}
    return retr;

}


cxsc::interval OscFobj::operator()(const cxsc::interval& ival) const
{
	cxsc::interval int_a(a, a);
	cxsc::interval reti = pow(ival,int_a) + (ival+1) * pow(cxsc::sin(ival*Pi_real*b),int_a) * pow(cxsc::cos(ival*Pi_real*c),int_a);
	if (Inf(ival) < 0 || Sup(ival) > 1) SetInf(reti, 0.0);
	return reti;
}

cxsc::real OscFobj::operator()(const cxsc::real& r) const
{
	
	cxsc::real retr = pow(r,a) + (r+1) * pow(cxsc::sin(r*Pi_real*b),a) * pow(cxsc::cos(r*Pi_real*c),a);
	if (r < 0.0 || r > 1.0) {
		 retr = 0.0;
	
	}
	return retr;
}

std::string OscFobj::getName() const
{
	return name;
}

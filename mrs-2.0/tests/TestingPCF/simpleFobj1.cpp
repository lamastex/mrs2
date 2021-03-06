
/*! \file
\brief Definitions for multivariate simple  
* function object class.
* 
* f(x) = x for 0 <= x <= 1
*/

#include "simpleFobj1.hpp"
#include "cxsc.hpp"
#include "toolz.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;


SimpleFobj1::SimpleFobj1() :name("simple1") {}

cxsc::interval SimpleFobj1::operator()(const cxsc::ivector& ivec) const
{
	int lb = Lb(ivec);
	int ub = Ub(ivec);
	
	cxsc::interval result( this->operator()(Inf(ivec[lb])), 
					this->operator()(Sup(ivec[lb])) );
	
	for (int i = lb+1; i <= ub; ++i) {
		result*= cxsc::interval( this->operator()(Inf(ivec[i])), 
					this->operator()(Sup(ivec[i])) );
	}
	
	return result;
}

cxsc::real SimpleFobj1::operator()(const cxsc::rvector& r) const
{
	int lb = Lb(r);
	int ub = Ub(r);
	
	cxsc::real retr = this->operator()(r[lb]);
	
	for (int i = lb+1; i <= ub && retr > 0.0; ++i) {
		retr *= (this->operator()(r[i]));
	}
    return retr;

}

cxsc::real SimpleFobj1::operator()(const cxsc::real& r) const
{
	cxsc::real retr(r);
	if (r < 0.0 || r > 1.0) {
		retr = 0.0;
	}
	return retr;
}

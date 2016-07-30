
/*! \filenction object class.

This example is f(x) = exp(-x^2)

*/

#include "ExampleFobj1D_1.hpp"
//#include "cxsc.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;




cxsc::interval ExampleMappedFobj1D_1::operator()(const cxsc::interval& ival) const
{
    cxsc::interval eval = -ival*ival;
    cxsc::interval retint(exp(Inf(eval)),exp(Sup(eval)));
    return retint;

}

cxsc::real ExampleMappedFobj1D_1::operator()(const cxsc::real& r) const
{
    cxsc::real retr = exp(-r*r);
    return retr;

}


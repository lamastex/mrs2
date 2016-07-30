
/*! \file
\brief MappedSPnode example 1-d function object class.

This example is f(x) = (0.5x)^2

*/

#include "ExampleFobj1D_2.hpp"
//#include "cxsc.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;




cxsc::interval ExampleMappedFobj1D_2::operator()(const cxsc::interval& ival) const
{
    cxsc::interval eval = (0.5*ival)*(0.5*ival);
    cxsc::interval retint(Inf(eval),Sup(eval));
    return retint;

}

cxsc::real ExampleMappedFobj1D_2::operator()(const cxsc::real& r) const
{
    cxsc::real retr = (0.5*r)*(0.5*r);
    return retr;

}


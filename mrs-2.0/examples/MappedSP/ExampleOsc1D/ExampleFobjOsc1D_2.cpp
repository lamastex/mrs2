
/*! \file
\brief MappedSPnode example 1-d oscillating function object class.

This example is \f$f(x) = (cosx)^a\f$

*/

#include "ExampleFobjOsc1D_2.hpp"

using namespace std;
using namespace subpavings;


ExampleFobjOsc1D_2::ExampleFobjOsc1D_2()
        : a(1.0) {};

ExampleFobjOsc1D_2::ExampleFobjOsc1D_2(int aa)
        : a(aa) {};


cxsc::interval ExampleFobjOsc1D_2::operator()(const cxsc::interval& ival) const
{
    return cxsc::power(cxsc::cos(ival),a);


}

cxsc::real ExampleFobjOsc1D_2::operator()(const cxsc::real& r) const
{
    return cxsc::power(cxsc::cos(r),a);
}


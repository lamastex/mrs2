
/*! \file
\brief MappedSPnode example 1-d oscillating function object class.

This example is \f$f(x) = (sinx)^b\f$

*/

#include "ExampleFobjOsc1D_3.hpp"

using namespace std;
using namespace subpavings;


ExampleFobjOsc1D_3::ExampleFobjOsc1D_3()
        : b(1.0) {};

ExampleFobjOsc1D_3::ExampleFobjOsc1D_3(int bb)
        : b(bb) {};


cxsc::interval ExampleFobjOsc1D_3::operator()(const cxsc::interval& ival) const
{
    return cxsc::power(cxsc::sin(ival),b);

}

cxsc::real ExampleFobjOsc1D_3::operator()(const cxsc::real& r) const
{
    return cxsc::power(cxsc::sin(r),b);

}


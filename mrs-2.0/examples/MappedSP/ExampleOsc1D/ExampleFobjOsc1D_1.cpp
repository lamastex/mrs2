
/*! \file
\brief MappedSPnode example 1-d oscillating function object class.

This example is \f$f(x) = exp(-ax^b)(1+csin(ax^b * tan(b\pi)))\f$

*/

#include "ExampleFobjOsc1D_1.hpp"

using namespace std;
using namespace subpavings;


ExampleFobjOsc1D_1::ExampleFobjOsc1D_1()
        : a(1.0), b(1.0), c(1.0) {};

ExampleFobjOsc1D_1::ExampleFobjOsc1D_1(cxsc::real aa, cxsc::real bb, cxsc::real cc)
        : a(aa), b(bb), c(cc) {};


cxsc::interval ExampleFobjOsc1D_1::operator()(const cxsc::interval& ival) const
{
    cxsc::interval pow_int(b,b);
    return exp(-a*cxsc::pow(ival,pow_int))*(1 + c*cxsc::sin(a*cxsc::pow(ival,pow_int)*cxsc::tan(b*PI)));

}

cxsc::real ExampleFobjOsc1D_1::operator()(const cxsc::real& r) const
{
    return exp(-a*cxsc::pow(r,b))*(1 + c*cxsc::sin(a*cxsc::pow(r,b)*cxsc::tan(b*PI)));

}


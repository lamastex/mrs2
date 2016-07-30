
/*! \file
\brief MappedSPnode example 1-d function object class.

This example is f(x) = sine(x)

*/

#include "ExampleFobjSine.hpp"
#include "ExampleFobjSinePI.hpp"

//#include "cxsc.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;



cxsc::interval ExampleMappedFobjSine::operator()(const cxsc::interval& ival) const
{
    return sin(PI*f*ival);

}

cxsc::real ExampleMappedFobjSine::operator()(const cxsc::real& r) const
{
    return sin(PI*f*r);


}



/*! \file
\brief MappedSPnode example 1-d function object class.

This example is f(x) = sum of weighted sine(x)

*/

#include "ExampleFobjSineSum.hpp"
#include "ExampleFobjSinePI.hpp"

//#include "cxsc.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;



cxsc::interval ExampleMappedFobjSineSum::operator()(const cxsc::interval& ival) const
{
    cxsc::interval sum(0.0, 0.0);
    int n = 1;
    while (n <= f)
    {
        sum += (1.0/n)*sin(2*PI*n*ival);
        n += 2;

    }
    return 4.0/PI*sum;

}

cxsc::real ExampleMappedFobjSineSum::operator()(const cxsc::real& r) const
{
    cxsc::real sum = 0;
    int n = 1;
    while (n <= f)
    {
        sum += (1.0/n)*sin(2*PI*n*r);
        n += 2;

    }
    return 4.0/PI*sum;

}


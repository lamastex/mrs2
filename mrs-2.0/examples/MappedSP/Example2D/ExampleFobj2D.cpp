
/*! \file
\brief MappedSPnode example 2-d function object class.

This example is f(x) = exp(-x*y)

*/

#include "ExampleFobj2D.hpp"
//#include "cxsc.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;




cxsc::interval ExampleMappedFobj2D::operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const
{
    cxsc::interval eval = -ival1*ival2;
    //std::cout << eval << endl;
    cxsc::interval retint(exp(Inf(eval)),exp(Sup(eval)));
    //std::cout << retint << endl;
    return retint;

}

cxsc::real ExampleMappedFobj2D::operator()(const cxsc::real& r1,
                                        const cxsc::real& r2) const
{
    cxsc::real retr = exp(-r1*r2);
    return retr;

}


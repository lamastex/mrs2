#ifndef __INT_H__
#define __INT_H__

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <time.h>
#include <stack>
#include "rmath.hpp"
#include "imath.hpp"
#include "interval.hpp"
#include "ivector.hpp"
#include "imatrix.hpp"
#include "itaylor.hpp"
#include "dim2taylor.hpp"

using namespace std;
using namespace cxsc;


interval integrateWithSplitting(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector), 
				const ivector &domain, int order, real tol);

interval integrate(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector), 
		   ivector &domain, int order);

interval integrateWithSplitting(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector, interval), 
				interval fhat, const ivector &domain, int order, real tol);

interval integrate(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector, interval), 
		   interval fhat, ivector &domain, int order);
#endif

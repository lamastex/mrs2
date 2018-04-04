#include "Int.h"
#include "dim2taylor.hpp"
#include <cmath>

using namespace std;
using namespace cxsc;

typedef taylor::dim2taylor d2t;
typedef taylor::dim2taylor_vector d2tv;


// Parameters specific to the Bivariate Gaussian target
real sigma_x = 1.0;
real sigma_y = 1.0;
real rho = 0;

d2t BiGOP (d2tv X, interval fhat) {
 // X[1].print_dim2taylor();
 // d2t result;
  d2t result = taylor::init_const(X[1].order(),interval(0.0));
  real det = 1.0/(2*M_PI*sigma_x*sigma_y*sqrt(1-sqr(rho)));  
  result = sqr(X[1]/sigma_x) + sqr(X[2]/sigma_y) - (2*rho*X[1]*X[2])/(sigma_x*sigma_y);
  result = det * exp (-((1.0/2*(1-sqr(rho))) * result)) - fhat;
  result = sqr(result);
  result = sqrt(result);
  
return result;
}

int main() {
  taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector, interval);
  ivector domain(2);
  

  interval fhat;
  fhat = interval(1.5, 1.5);
  domain[1]=interval(0, 10.0);
  domain[2]=interval(-10.0, 0);
  cout << "Box is: " << domain << endl;

  real tol=1e-6;
  int o=6;
  interval result;
  testpnt=BiGOP; 

  result=integrateWithSplitting(testpnt,fhat,domain,o,tol);


  return 1;
}


/*taylor::dim2taylor testfcn(taylor::dim2taylor_vector x) {
  //return x[1]*x[2];  //[0.25000,0.25000]
  return exp(-(sqr(x[1])+sqr(x[2]))/2.0);
}*/


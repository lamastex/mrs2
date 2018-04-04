#include "Int.h"
#include "dim2taylor.hpp"
#include <cmath>
#include <vector>
#include <iterator>
#include <fstream>
#include <sstream>
#include <time.h>


using namespace std;
using namespace cxsc;

typedef taylor::dim2taylor d2t;
typedef taylor::dim2taylor_vector d2tv;


// Parameters specific to the Rosenbrock target
real Tinverse = 1.0;
real Height = 100.0;
real DomainLimit = 5.0;

d2t RosenOP (d2tv X) {
 // X[1].print_dim2taylor();
 // d2t result;
  d2t result = taylor::init_const(X[1].order(),interval(0.0));

  for (int i = 1; i < 2; i++) //2nd term should be size_k
    {
      result = result + (Height * sqr(X[i+1] - sqr(X[i])) +
        sqr(X[i] - 1.0));
    }
  
  result = exp (-(Tinverse * result));

  return result;
}

int main() {
  
	real tol=1e-7;
	int o=16;

	clock_t start, end;
	start = clock();

	taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector);
	ivector domain(2);

	domain[1]=interval(-2.960054820867426, 3.186979924344668);
	domain[2]=interval(-0.337486934775497, 9.999699890093452);
	cout << "Box is: " << domain << endl;
	
	interval result;
	testpnt=RosenOP; 

	result=integrateWithSplitting(testpnt,domain,o,tol);

	end = clock();
	double timing;
	timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
	cout << "Computing time : " << timing << " s."<< endl;

  return 1;
}


/*taylor::dim2taylor testfcn(taylor::dim2taylor_vector x) {
  //return x[1]*x[2];  //[0.25000,0.25000]
  return exp(-(sqr(x[1])+sqr(x[2]))/2.0);
}*/


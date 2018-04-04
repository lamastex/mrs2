#include "Int.h"
#include "dim2taylor.hpp"

using namespace std;
using namespace cxsc;


typedef taylor::dim2taylor d2t;
typedef taylor::dim2taylor_vector d2tv;


d2t testfcn(d2tv x) {
  //return x[1]*x[2];  //[0.25000,0.25000]
  return exp(-(sqr(x[1])+sqr(x[2]))/2.0);
}


int main() {
  taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector);
  ivector domain(2);

  domain[1]=interval(-10,10);
  domain[2]=interval(-10,10);
 
  cout<<SetPrecision(9,8);
  interval result(-100,100);
  //real tol=0.0001;
  //for(int o=2; o<32; o+=2) {
  real tol=0.01;
  for(int o=2; o<8; o+=2) {
    //  int o=12;
    testpnt=testfcn;
    result=result&integrateWithSplitting(testpnt,domain,o,tol); 
    cout<<"order = "<<o<<". Result = "<<result<<endl;
  }
  return 1;
}



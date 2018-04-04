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

//Parameters specific to the Levy target
real Temperature = 40.0;
real Center1 = 1.42513; 
real Center2 = 0.80032; 
real GlobalMax = 176.14;
real DomainLimit = 10.0;

d2t LevyOP (d2tv X) {
  //X[1].print_dim2taylor();
  d2t isum = taylor::init_const(X[1].order(),interval(0.0));
  d2t jsum = taylor::init_const(X[1].order(),interval(0.0));

  for (int i = 1; i <= 5; i++)
  {
    isum = isum + i * cos ((i - 1) * X[1] + (i));
    jsum = jsum + i * cos ((i + 1) * X[2] + (i));
  }
                    // Avoid real conversion error
  d2t hh = isum * jsum + sqr (X[1] + Center1) +
    sqr (X[2] + Center2);
  hh = hh + GlobalMax;  
  // TEMPERATURE = 1, 4, 40, 400, 4000
  d2t result = exp (-hh / Temperature);

  return result;
}

int main() {
  taylor::dim2taylor (*testpnt)(taylor::dim2taylor_vector);
  ivector domain(2);
  
 double minSide = 0.15625;

vector<double> domain1low;
vector<double>::iterator it_domain1low;

vector<double> domain1upp;
vector<double>::iterator it_domain1upp;

vector<double> domain2low;
vector<double>::iterator it_domain2low;
vector<double> domain2upp;
vector<double>::iterator it_domain2upp;

vector<real> areaLevy;
vector<real>::iterator it_area;


  real tol=1e-7;
  int o=16;


clock_t start, end;
start = clock();

int totalSide=128;
for (int i=0; i < totalSide; i++)
  for (int j=0; j<totalSide;j++)
{

cout << "i: " << i << endl;
cout << "j: " << j << endl;

double domain1low1, domain1upp1, domain2low1, domain2upp1;

domain1low1 = -10+minSide*i;
domain1upp1 = -10+minSide*(i+1);
domain2low1 = -10+minSide*j;
domain2upp1 = -10+minSide*(j+1);

domain[1] = interval(-10+minSide*i, -10+minSide*(i+1));
domain[2] = interval(-10+minSide*j, -10+minSide*(j+1));

cout << domain[1] << domain[2] << endl;

 
  interval result;
  testpnt=LevyOP; 
  result=integrateWithSplitting(testpnt,domain,o,tol);

real area1 = Sup(result);  



domain1low.push_back(domain1low1);
domain1upp.push_back(domain1upp1);
domain2low.push_back(domain2low1);
domain2upp.push_back(domain2upp1);
areaLevy.push_back(area1);

}

end = clock();
double timing;
  timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
        cout << "Computing time : " << timing << " s."<< endl;
    



string outputFileName1, outputFileName2, outputFileName3, outputFileName4, outputFileName5;

ofstream oss;
oss << scientific;
oss.precision(5);

outputFileName1 = "Levy4096domain1low.txt";
outputFileName2 = "Levy4096domain1upp.txt";
outputFileName3 = "Levy4096domain2low.txt";
outputFileName4 = "Levy4096domain2upp.txt";
outputFileName5 = "Levy4096area.txt";

oss.open(outputFileName1.c_str());
for (it_domain1low = domain1low.begin(); it_domain1low < domain1low.end(); it_domain1low++){
oss << (*it_domain1low) << endl;
}
oss << flush;
oss.close();
cout << "Domain1low output to " << outputFileName1 << endl;


oss.open(outputFileName2.c_str());
for (it_domain1upp= domain1upp.begin(); it_domain1upp < domain1upp.end(); it_domain1upp++){
oss << (*it_domain1upp) << endl;
}
oss << flush;
oss.close();
cout << "Domain2 output to " << outputFileName2 << endl;


oss.open(outputFileName3.c_str());
for (it_domain2low = domain2low.begin(); it_domain2low < domain2low.end(); it_domain2low++){
oss << (*it_domain2low) << endl;
}
oss << flush;
oss.close();
cout << "Domain2low output to " << outputFileName3 << endl;


oss.open(outputFileName4.c_str());
for (it_domain2upp= domain2upp.begin(); it_domain2upp < domain2upp.end(); it_domain2upp++){
oss << (*it_domain2upp) << endl;
}
oss << flush;
oss.close();
cout << "Domain2upp output to " << outputFileName4 << endl;

oss.open(outputFileName5.c_str());
for (it_area = areaLevy.begin(); it_area < areaLevy.end(); it_area++){
oss << (*it_area) << endl;
}
oss << flush;
oss.close();
cout << "Area output to " << outputFileName5 << endl;








  return 1;

}


/*taylor::dim2taylor testfcn(taylor::dim2taylor_vector x) {
  //return x[1]*x[2];  //[0.25000,0.25000]
  return exp(-(sqr(x[1])+sqr(x[2]))/2.0);
}*/


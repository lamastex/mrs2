
/*! \file
\brief Test RosenFobj.


*/

#include "RosenFobj.hpp"


#include <iostream>


using namespace std;
using namespace subpavings;


int main(int argc, char* argv[])
{

    int dims = 2;
	double intside = 10.0;
	
	if (argc > 1) dims = atoi(argv[1]);
	if (argc > 2) intside = atof(argv[2]);
	
    cxsc::ivector pavingBox(dims);
    cxsc::interval pavingInterval(-intside,intside);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
	
	RosenFobj realF;
		
	cxsc::ivector testBox(dims);
	cxsc::interval testInterval(-intside/2,intside/2);
	for(int k=1; k <= dims; k++) testBox[k] = testInterval;
	cxsc::interval ival = realF(testBox);
	
	cxsc::rvector testRvec(dims);
	cxsc::real testReal(intside/2);
	for(int k=1; k <= dims; k++) testRvec[k] = testReal;
	cxsc::real r = realF(testRvec);
	
	cout << "testBox = " << testBox << endl;
	cout << "Rosen(testBox) = " << ival << endl;
	cout << "testRvec = " << testRvec << endl;
	cout << "Rosen(testRvec) = " << r << endl;


	
    return 0;

}



/*! \file
\brief Test LevyDensityFobj2D.


*/

#include "LevyDensityFobj2D.hpp"


#include <iostream>

using namespace std;
using namespace subpavings;


int main(int argc, char* argv[])
{

    cout << "\nLevy Density\n" << endl;

    int dims = 2; // Levy can only do 2D
	double intside = 10.0;
	
	if (argc > 1) dims = atoi(argv[1]);
	if (argc > 2) intside = atof(argv[2]);
	
    cxsc::ivector pavingBox(dims);
    cxsc::interval pavingInterval(-intside,intside);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
	
	LevyDensityFobj2D realF;
		
	cxsc::ivector testBox(dims);
	cxsc::interval testInterval(-intside/2,intside/2);
	for(int k=1; k <= dims; k++) testBox[k] = testInterval;
	cxsc::interval ival = realF(testBox);
	
	cxsc::rvector testRvec(dims);
	cxsc::real testReal(intside/2);
	for(int k=1; k <= dims; k++) testRvec[k] = testReal;
	cxsc::real r = realF(testRvec);
	
	cout << "testBox = " << testBox << endl;
	cout << "Levy(testBox) = " << ival << endl;
	cout << "testRvec = " << testRvec << endl;
	cout << "Levy(testRvec) = " << r << endl;

	
			
	
    return 0;

}

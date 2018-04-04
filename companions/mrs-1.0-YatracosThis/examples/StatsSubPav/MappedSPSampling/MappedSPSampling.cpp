/*! \file
\brief MappedSPnode example for Gaussian objects main.
*/

#include "GaussianFobj1D.hpp"
#include "GaussianFobj2D.hpp"
#include "GaussianFobj10D.hpp"
#include "GaussianFobj100D.hpp"
//#include "GaussianFobj1000D.hpp"
#include "RosenFobj2D.hpp"
#include "RosenFobj10D.hpp"
#include "RosenFobj100D.hpp"
//#include "RosenFobj1000D.hpp"
#include "LevyFobj2D.hpp"

#include "realmappedspnode.hpp"

#include "mappedspnodevisitor_expand.hpp"

#include <fstream>  // for ifstream, ofstream

using namespace std;
using namespace subpavings;

void output(string& filename,  const SPnode& node)
{
        // To generate a file output
        ofstream os(filename.c_str());         // Filename, c-string version
        if (os.is_open()) {

            node.leavesOutputTabs(os); // the output
            std::cout << "The output of the example"
                    << " has been written to " << filename << std::endl << std::endl;
            os.close();
        }
        else {
            std::cerr << "Error: could not open file named "
                << filename << std::endl << std::endl;
        }
}

int main(int argc, char* argv[])
{
    //=======user defined parameters================================//
    	if ( argc != 5 ) {
    cerr << "Syntax: MappedFunctions dims lb ub numLeaves" << endl;
    exit(0);
	}

    int dims = atoi(argv[1]);
    //double tolerance = atof(argv[1]);
    
	   size_t critLeaves = atof(argv[4]);
    cout << critLeaves << endl;
    
    real tolerance = 0;
   
    //========define the function object
    //RosenFobj2D realF;
    //LevyFobj2D realF;
    GaussianFobj10D realF;
    //GaussianFobj2D realF;
    //=====make the root box===========================//
    ivector pavingBox(dims);
    real lb = atof(argv[2]);
    real ub = atof(argv[3]);
    interval pavingInterval(lb, ub);
    
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
    //RealMappedSPnode nodeRoot1(pavingBox); // make a MappedSPnode object
    RealMappedSPnode nodeRoot2(pavingBox); // make a MappedSPnode object

	
   // interval ival = realF(pavingBox);

    MappedSPnodeVisitorExpand expander(realF, tolerance);
    //nodeRoot1.accept(expander);    
    nodeRoot2.priorityAccept(expander, critLeaves);
	 
	 cout.precision(20);
	 cout << "Tolerance is: " << tolerance << endl;
	 //cout << "Number of leaves is: " << nodeRoot1.getNumLeaves() << endl;
	cout << "Number of leaves is: " << nodeRoot2.getNumLeaves() << endl;

	//cout << nodeRoot.getLeafLevelsString() << endl;
   //string filename = "EstFunction1.txt";
   //output(filename, nodeRoot1);

	string filename = "EstFunction2.txt";
   output(filename, nodeRoot2);

    return 0;

}

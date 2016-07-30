
/*! \file
\brief MappedSPnode example 1-d main.


*/

#include "ExampleFobj1D_1.hpp"
#include "ExampleFobj1D_2.hpp"

#include "realmappedspnode.hpp"
#include "realexpander_estimate.hpp"

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

int main()
{

        cout << "\nEx1D\n" << endl;

        cxsc::real tolerance = 0.07;

        int dims = 1;
        ivector pavingBox(dims);
        cxsc::interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


        RealMappedSPnode nodeRoot1(pavingBox); // make a MappedSPnode object
        RealMappedSPnode nodeRoot2(pavingBox); // make a MappedSPnode object

        ExampleMappedFobj1D_1 realF1;
        ExampleMappedFobj1D_2 realF2;

        RealExpanderEstimator estimator1(realF1, tolerance);
        RealExpanderEstimator estimator2(realF2, tolerance);

        nodeRoot1.acceptSPExpandVisitor(estimator1);
        nodeRoot2.acceptSPExpandVisitor(estimator2);

        string filename = "Ex1D_pt07_1.txt";
        output(filename, nodeRoot1);

        filename = "Ex1D_pt07_2.txt";
        output(filename, nodeRoot2);


        RealMappedSPnode nodeRoot = nodeRoot1 + nodeRoot2;
        filename = "Ex1D_pt07_added.txt";
        output(filename, nodeRoot);



    return 0;

}

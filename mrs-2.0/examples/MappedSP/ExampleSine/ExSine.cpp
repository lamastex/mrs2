
/*! \file
\brief MappedSPnode example 1-d main.

This example is f(x) = exp(-x^2)

*/

#include "ExampleFobjSine.hpp"
#include "ExampleFobjSineSum.hpp"
#include "ExampleFobjSinePI.hpp"

#include "realmappedspnode.hpp"

#include "realexpander_estimate.hpp"

#include <sstream>  // to be able to manipulate strings as streams
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

    cxsc::real tolerance = 0.05;

    string filenameRoot = "ExSine";

    cout << "\nExSine\n" << endl;

    int dims = 1;
    ivector pavingBox(dims);
    cxsc::interval pavingInterval(-2,2);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

    int f = 1;

    RealMappedSPnode nodeRootSquare (pavingBox);

    for (f = 1; f < 6; f+=2) {


        RealMappedSPnode nodeRoot(pavingBox); // make a MappedSPnode object

        ExampleMappedFobjSine realF(2*f);

        RealExpanderEstimator estimator(realF, tolerance);

        nodeRoot.acceptSPExpandVisitor(estimator);

        std::ostringstream stm;
        stm << filenameRoot << "_f" << f << "_pt_" << static_cast<int>(_double(tolerance)*100) << ".txt";
        string filename = stm.str();

        output(filename, nodeRoot);

        cxsc::real mult(4.0/(f*PI));

        nodeRoot *= mult;

        nodeRootSquare += nodeRoot;

        stm.str("");
        stm << filenameRoot << "_f" << f << "_square_pt_" << static_cast<int>(_double(tolerance)*100) << ".txt";
        filename = stm.str();

        output(filename, nodeRootSquare);

    }

    for (f=f; f < 20; f+=2) {


        RealMappedSPnode nodeRoot(pavingBox); // make a MappedSPnode object

        ExampleMappedFobjSine realF(2*f);

        RealExpanderEstimator estimator(realF, tolerance);

        nodeRoot.acceptSPExpandVisitor(estimator);

        cxsc::real mult(4.0/(f*PI));

        nodeRoot *= mult;

        nodeRootSquare += nodeRoot;

    }
    f-=2;

    std::ostringstream stm;
    stm << filenameRoot << "_f" << f << "_square_pt_" << static_cast<int>(_double(tolerance)*100) << ".txt";
    string filename = stm.str();

    output(filename, nodeRootSquare);

    // do the mapped subpaving for the whole function
    RealMappedSPnode nodeRootOverall(pavingBox); // make a MappedSPnode object

    ExampleMappedFobjSineSum realFOverall(f);

    RealExpanderEstimator estimator(realFOverall, tolerance);

    nodeRootOverall.acceptSPExpandVisitor(estimator);

    stm.str("");
    stm << filenameRoot << "_f" << f << "_overall_pt_" << static_cast<int>(_double(tolerance)*100) << "_n.txt";
    filename = stm.str();

    output(filename, nodeRootOverall);

    // look at the difference between nodeRootSquare and nodeRootOverall
    RealMappedSPnode nodeRootOverallDifference = nodeRootSquare - nodeRootOverall;

    stm.str("");
    stm << filenameRoot << "_f" << f << "_difference_pt_" << static_cast<int>(_double(tolerance)*100) << "_n.txt";
    filename = stm.str();

    output(filename, nodeRootOverallDifference);

    return 0;

}

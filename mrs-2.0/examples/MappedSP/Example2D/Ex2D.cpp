
/*! \file
\brief MappedSPnode example 2-d main.


*/

#include "ExampleFobj2D.hpp"

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

    cout << "\nEx2D\n" << endl;

    int dims = 2;
    ivector pavingBox(dims);
    cxsc::interval pavingInterval(-1,1);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


    RealMappedSPnode nodeRoot(pavingBox); // make a MappedSPnode object

    ExampleMappedFobj2D realF;
    cxsc::interval ival = realF(cxsc::interval(-1.0,1.0),cxsc::interval(-1.0,1.0));
    cout << ival << endl;


    RealExpanderEstimator estimator(realF, 0.05);

    nodeRoot.acceptSPExpandVisitor(estimator);

    string filename = "Ex2D_pt05.txt";
    output(filename, nodeRoot);



    return 0;

}

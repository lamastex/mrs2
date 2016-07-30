
/*! \file
\brief MappedSPnode example 1-d oscillating main.


*/

#include "ExampleFobjOsc1D_1.hpp"
#include "ExampleFobjOsc1D_2.hpp"
#include "ExampleFobjOsc1D_3.hpp"
#include "ExampleFobjOscPI.hpp"


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

        cout << "\nExOsc1D\n" << endl;

/*
        cxsc::real tolerance = 0.0001;

        int dims = 1;
        ivector pavingBox(dims);
        cxsc::interval pavingInterval((1.0/10.0),10000.0);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


        RealMappedSPnode nodeRoot1(pavingBox); // make a MappedSPnode object
        //RealMappedSPnode nodeRoot2(pavingBox); // make a MappedSPnode object

        cxsc::real a = 1.0/8.0;
        cxsc::real b = 9.0/20.0;
        cxsc::real c = 1.0/2.0;


        ExampleFobjOsc1D_1 realF1(a,b,c);
        //ExampleMappedFobj1D_2 realF2;

        RealExpanderEstimator estimator1(realF1, tolerance);
        //RealExpanderEstimator estimator2(realF2, tolerance);

        nodeRoot1.acceptSPExpandVisitor(estimator1);
        //nodeRoot2.acceptSPExpandVisitor(estimator2);

        string filename = "Ex1DOsc_pt0001_1.txt";
        output(filename, nodeRoot1);

        //filename = "Ex1D_pt1_2.txt";
        //output(filename, nodeRoot2);



*/
        cxsc::real tolerance = 0.01;

        int dims = 1;
        ivector pavingBox(dims);
        cxsc::interval pavingInterval(-5.0,5.0);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


        RealMappedSPnode nodeRoot2(pavingBox); // make a MappedSPnode object
        RealMappedSPnode nodeRoot3(pavingBox); // make a MappedSPnode object

        int a = 3;
        ExampleFobjOsc1D_2 realF2(a);

        int b = 1;
        ExampleFobjOsc1D_3 realF3(b);

        RealExpanderEstimator estimator2(realF2, tolerance);
        RealExpanderEstimator estimator3(realF3, tolerance);

        nodeRoot2.acceptSPExpandVisitor(estimator2);
        nodeRoot3.acceptSPExpandVisitor(estimator3);

        string filename = "Ex1D_Cos_pt01_2.txt";
        output(filename, nodeRoot2);

        filename = "Ex1D_Sin_pt01_3.txt";
        output(filename, nodeRoot3);

        RealMappedSPnode nodeRootSum = nodeRoot2 + nodeRoot3;
        filename = "Ex1D_CosPlusSin_pt01_3.txt";
        output(filename, nodeRootSum);

        RealMappedSPnode nodeRootSubt = nodeRoot2 - nodeRoot3;
        filename = "Ex1D_CosMinusSin_pt01_3.txt";
        output(filename, nodeRootSubt);


    return 0;

}

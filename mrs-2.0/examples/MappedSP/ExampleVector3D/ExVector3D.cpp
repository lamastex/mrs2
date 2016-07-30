
/*! \file
\brief MappedSPnode example with vectors main.

This example is f(x) = exp(-x^2)

*/

#include "Vec.hpp"


#include "mappedspnode.hpp"


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

    cout << "\nExVector3D\n" << endl;

    int dims = 2;
    ivector pavingBox(dims);
    cxsc::interval pavingInterval(-4,4);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


    MappedSPnode<Vec> nodeOne(pavingBox); // make a MappedSPnode object
    MappedSPnode<Vec> nodeTwo(pavingBox); // make a MappedSPnode object

    // split the nodes to particular shapes
    nodeOne.splitToShape("3,3,2,1"); // split it
    nodeTwo.splitToShape("3,4,4,2,2,3,3"); // split it


    // allocate ranges for nodeOne
    Vec rangeOne1(0.0,0.0,0.0);
    Vec rangeOne2(0.0,0.0,0.0);
    Vec rangeOne3(0.0,0.0,0.0);
    Vec rangeOne4(0.3,0.1,0.3);
    Vec rangeOne5(0.2,0.1,-1.0);
    Vec rangeOne6(0.3,0.6,0.5);
    Vec rangeOne7(0.2,0.2,0.8);


    vector< Vec > rangesOne;
    rangesOne.push_back(rangeOne1);
    rangesOne.push_back(rangeOne2);
    rangesOne.push_back(rangeOne3);
    rangesOne.push_back(rangeOne4);
    rangesOne.push_back(rangeOne5);
    rangesOne.push_back(rangeOne6);
    rangesOne.push_back(rangeOne7);



    nodeOne.allocateRanges(rangesOne, 0);


    // allocate ranges for nodeTwo
    Vec  rangeTwo1(0.0,0.0,0.0);
    Vec  rangeTwo2(0.0,0.0,0.0);
    Vec  rangeTwo3(0.0,0.0,0.0);
    Vec  rangeTwo4(0.2,0.5,0.-.8);
    Vec  rangeTwo5(0.0,0.0,0.0);
    Vec  rangeTwo6(0.-1,0.-.5,0.3);
    Vec  rangeTwo7(0.0,0.0,0.6);
    Vec  rangeTwo8(-1.0,0.5,0.5);
    Vec  rangeTwo9(0.0,0.0,0.0);
    Vec  rangeTwo10(0.5,-.3,0.3);
    Vec  rangeTwo11(0.0,0.0,0.0);
    Vec  rangeTwo12(-0.5,0.4,0.3);
    Vec  rangeTwo13(0.3,-.7,0.1);

    vector< Vec > rangesTwo;
    rangesTwo.push_back(rangeTwo1);
    rangesTwo.push_back(rangeTwo2);
    rangesTwo.push_back(rangeTwo3);
    rangesTwo.push_back(rangeTwo4);
    rangesTwo.push_back(rangeTwo5);
    rangesTwo.push_back(rangeTwo6);
    rangesTwo.push_back(rangeTwo7);
    rangesTwo.push_back(rangeTwo8);
    rangesTwo.push_back(rangeTwo9);
    rangesTwo.push_back(rangeTwo10);
    rangesTwo.push_back(rangeTwo11);
    rangesTwo.push_back(rangeTwo12);
    rangesTwo.push_back(rangeTwo13);

    nodeTwo.allocateRanges(rangesTwo, 0);


    string filename = "ExVec3D_1.txt";
    output(filename, nodeOne);

    filename = "ExVec3D_2.txt";
    output(filename, nodeTwo);

   // ----------- try addition

    MappedSPnode<Vec> addition = nodeOne + nodeTwo;
    filename = "ExVec3D_Addition.txt";
    output(filename, addition);

    // ----------- try subtraction

    MappedSPnode<Vec> subtraction = nodeOne - nodeTwo;
    filename = "ExVec3D_Subtraction.txt";
    output(filename, subtraction);

    // ----------- try multiplication

    MappedSPnode<Vec> multiplication = nodeOne * nodeTwo;
    filename = "ExVec3D_Multiplication.txt";
    output(filename, multiplication);

    return 0;

}

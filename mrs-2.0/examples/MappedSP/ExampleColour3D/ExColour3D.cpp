
/*! \file
\brief MappedSPnodeSingleRange example with colours 3-d main.

*/

#include "RGBColour.hpp"


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

    cout << "\nExColour3D\n" << endl;

    int dims = 3;
    ivector pavingBox(dims);
    cxsc::interval pavingInterval(-4,4);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


    MappedSPnode<RGBColour> nodeOne(pavingBox); // make a MappedSPnode object

    MappedSPnode<RGBColour> nodeTwo(pavingBox); // make a MappedSPnode object

    // split the nodes to particular shapes
    nodeOne.splitToShape("3,3,2,1"); // split it
    nodeTwo.splitToShape("3,4,4,2,2,3,3"); // split it


    // allocate ranges for nodeOne
    RGBColour rangeOne1(1.0,1.0,1.0);
    RGBColour rangeOne2(0.8,0.8,0.8);
    RGBColour rangeOne3(0.5,0.2,0.3);
    RGBColour rangeOne4(1.0,1.0,0.0);
    RGBColour rangeOne5(0.8,0.1,0.0);
    RGBColour rangeOne6(0.3,0.6,0.5);
    RGBColour rangeOne7(1.0,0.9,0.9);


    vector< RGBColour> rangesOne;
    rangesOne.push_back(rangeOne1);
    rangesOne.push_back(rangeOne2);
    rangesOne.push_back(rangeOne3);
    rangesOne.push_back(rangeOne4);
    rangesOne.push_back(rangeOne5);
    rangesOne.push_back(rangeOne6);
    rangesOne.push_back(rangeOne7);

	nodeOne.allocateRanges(rangesOne, 0);

	// allocate ranges for nodeTwo
    RGBColour  rangeTwo1(1.0,1.0,1.0);
    RGBColour  rangeTwo2(0.6,0.6,0.6);
    RGBColour  rangeTwo3(0.3,0.5,0.6);
    RGBColour  rangeTwo4(0.2,0.5,0.0);
    RGBColour  rangeTwo5(0.1,0.0,0.6);
    RGBColour  rangeTwo6(0.0,0.5,0.5);
    RGBColour  rangeTwo7(0.0,0.0,0.6);
    RGBColour  rangeTwo8(0.8,0.9,1.0);
    RGBColour  rangeTwo9(0.4,0.4,0.4);
    RGBColour  rangeTwo10(0.1,0.0,0.3);
    RGBColour  rangeTwo11(0.3,0.4,0.1);
    RGBColour  rangeTwo12(0.0,0.4,0.0);
    RGBColour  rangeTwo13(0.3,0.0,0.1);

    vector< RGBColour > rangesTwo;
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


    string filename = "ExColour3D_1.txt";
    output(filename, nodeOne);

    filename = "ExColour3D_2.txt";
    output(filename, nodeTwo);

    // ----------- try addition

    MappedSPnode<RGBColour> addition = nodeOne + nodeTwo;
    filename = "ExColour3D_Addition.txt";
    output(filename, addition);

    // ----------- try subtraction

    MappedSPnode<RGBColour> subtraction = nodeOne - nodeTwo;
    filename = "ExColour3D_Subtraction.txt";
    output(filename, subtraction);
	
	// ----------- Try slice
	try {
		cout << "\nslices on copy of node1 on dim 1" << endl;
		
		{
			MappedSPnode<RGBColour> temp1(nodeOne);
		
			vector < int > sliceDims;
			int d = 1;
			sliceDims.push_back(d);
			
			vector < real > slicePts;
			double dr = -2.0;
			slicePts.push_back(real(dr));
			
			temp1.slice(sliceDims, slicePts);
			
			ostringstream oss;
			oss << "nodeOnesliced_d" << d << "_" << dr << ".txt";
			string s = oss.str();
			output(s, temp1);
		}
		cout << "\nslices on copy of node1 on dim 2" << endl;
		{
			MappedSPnode<RGBColour> temp1(nodeOne);
		
			vector < int > sliceDims;
			int d = 2;
			sliceDims.push_back(d);
			
			vector < real > slicePts;
			double dr = -2.0;
			slicePts.push_back(real(dr));
			
			temp1.slice(sliceDims, slicePts);
			
			ostringstream oss;
			oss << "nodeOnesliced_d" << d << "_" << dr << ".txt";
			string s = oss.str();
			output(s, temp1);
		}
		{
			MappedSPnode<RGBColour> temp1(nodeOne);
		
			vector < int > sliceDims;
			int d = 2;
			sliceDims.push_back(d);
			
			vector < real > slicePts;
			double dr = 2.0;
			slicePts.push_back(real(dr));
			
			temp1.slice(sliceDims, slicePts);
			
			ostringstream oss;
			oss << "nodeOnesliced_d" << d << "_" << dr << ".txt";
			string s = oss.str();
			output(s, temp1);
		}
		{
			MappedSPnode<RGBColour> temp1(nodeOne);
		
			vector < int > sliceDims;
			int d1 = 1;
			sliceDims.push_back(d1);
			int d2 = 3;
			sliceDims.push_back(d2);
			
			
			vector < real > slicePts;
			double dr1 = -2.0;
			slicePts.push_back(real(dr1));
			double dr2 = -2.0;
			slicePts.push_back(real(dr2));
			
			temp1.slice(sliceDims, slicePts);
			
			ostringstream oss;
			oss << "nodeOnesliced_d" << d1 << "," << d2 << "_" << dr1 << "," << dr2 << ".txt";
			string s = oss.str();
			output(s, temp1);
		}
		{
			MappedSPnode<RGBColour> temp1(nodeOne);
		
			vector < int > sliceDims;
			int d1 = 1;
			sliceDims.push_back(d1);
			int d2 = 3;
			sliceDims.push_back(d2);
			
			
			vector < real > slicePts;
			double dr1 = -2.0;
			slicePts.push_back(real(dr1));
			double dr2 = 2.0;
			slicePts.push_back(real(dr2));
			
			temp1.slice(sliceDims, slicePts);
			
			ostringstream oss;
			oss << "nodeOnesliced_d" << d1 << "," << d2 << "_" << dr1 << "," << dr2 << ".txt";
			string s = oss.str();
			output(s, temp1);
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do slice of node 1:\n" << msg << endl;
	}
	

    return 0;

}

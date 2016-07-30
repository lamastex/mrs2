
/*! \file
\brief MappedSPnode examples 2-d main.


*/

#include "RosenFobj.hpp"

#include "realmappedspnode.hpp"
#include "realexpander_estimate.hpp"
#include "intervalmappedspnode.hpp"
#include "intervalexpander_estimate.hpp"
#include "functionestimator_interval.hpp"
#include "fei_evalobj.hpp"

#include <ostream>
#include <fstream>

using namespace std;
using namespace subpavings;

void output(string& filename,  const IntervalMappedSPnode& node);
void output(string& filename,  const RealMappedSPnode& node);



int main(int argc, char* argv[])
{
	
	cout << "\nRosenbrock\n" << endl;

    int dims = 2;
	double intside = 1.0;
	
	if (argc > 1) dims = atoi(argv[1]);
	if (argc > 2) intside = atof(argv[2]);
	
    cxsc::ivector pavingBox(dims);
    cxsc::interval pavingInterval(-intside,intside);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
	
	RosenFobj realF;
	
	/* expanding a box using a RealExpanderEstimator to create an
	 * approximation of the function with a tree of RealMappedSPnodes.
	 * This is effectively doing, at the node level, what a 
	 * FunctionEstimatorReal is be managing for us, ie this is exposing
	 * the the nodes themselves.  This is really just for testing and
	 * for interest.  It is better to use the FunctionEstimatorReal
	 * (or the FunctionEstimatorInterval) */ 
	{
		cout << "\n\nReal mapping, using RealExpanderEstimator" << endl;
		
		RealMappedSPnode nodeRoot(pavingBox); // make a MappedSPnode object

		// have to specify a tolerance for the ExanderEstimator
		RealExpanderEstimator estimator(realF, 5.0);

		cout << "About to expand ..." << endl;

		nodeRoot.acceptSPExpandVisitor(estimator);

		ostringstream oss;  
		oss << "Rosen" << dims << "DCoarseReal.txt";
		std::string filename = oss.str();
		output(filename, nodeRoot);
	}
	
	/* expanding a box using an IntervalExpanderEstimator to create an
	 * approximation of the function with a tree of IntervalMappedSPnodes.
	 * This is effectively doing, at the node level, what a 
	 * FunctionEstimatorInterval is be managing for us, ie this is exposing
	 * the the nodes themselves.  This is really just for testing and
	 * for interest.  It is better to use the FunctionEstimatorInterval
	 * (or the FunctionEstimatorReal) */ 
	{
		cout << "\n\nInterval mapping" << endl;
		
		IntervalMappedSPnode nodeRoot(pavingBox); // make a MappedSPnode object

		// have to specify a tolerance for the ExanderEstimator
		IntervalExpanderEstimator estimator(realF, 17.67);

		cout << "About to expand ..." << endl;

		nodeRoot.acceptSPExpandVisitor(estimator);

		ostringstream oss;  
		oss << "Rosen" << dims << "DCoarseInterval.txt";
		std::string filename = oss.str();
		output(filename, nodeRoot);
	}
	
	/* Creating a function estimate using a priority queue, using the
	 * default (Reimann) priority measure.*/
	{
		cout << "\n\npriority queue estimator" << endl;
		
		FunctionEstimatorInterval fei(pavingBox, realF);
				
		size_t maxLeaves = 7000;
		
		LOGGING_LEVEL logging = NOLOG;
		
		fei.prioritySplit(maxLeaves, logging);
				
		ostringstream oss;  
		oss << "Rosen" << dims << "D_PQ_l" << maxLeaves << ".txt";
		std::string filename = oss.str();
		int prec = 5;
		fei.outputToTxtTabs(filename, prec, true);
	}
	return 0;

}

void output(string& filename,  const RealMappedSPnode& node)
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

void output(string& filename,  const IntervalMappedSPnode& node)
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

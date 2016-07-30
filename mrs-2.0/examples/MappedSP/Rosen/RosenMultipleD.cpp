
/*! \file
\brief Show gain from spilt and pull up for lots of different dimensions.

Run code as used for Fuzz-IEEE paper to show interval enclosures using
Reimann-Diff priority function.


*/

#include "RosenFobj.hpp"


#include "functionestimator_interval.hpp"

#include "sptools.hpp"

#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace subpavings;




int main(int argc, char* argv[])
{
	
	
    cout << "\nRosenbrock with Reimann Measurer\n" << endl;

    int minDim = 2;
    int maxDim = 3;
    
    std::vector < int > power_vec; 
	{
		int tmp[] = {3, 4, 5, 6, 7}; 
		power_vec.insert (power_vec.begin(), tmp, tmp+sizeof(tmp) / sizeof(int));
	}
    
	double intside = 1.0;
	
	cxsc::interval pavingInterval(-intside,intside);
    
    LOGGING_LEVEL logging = NOLOG;
				
    for (int dims = minDim; dims < maxDim + 1; ++dims) {
		
		cout << "\n\nD = " << dims << endl;
					
		cxsc::ivector pavingBox(dims);
		cxsc::interval pavingInterval(-intside,intside);
		for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
		
		//start file for this D
		std::string thisFilename;
		{
			ostringstream oss;
			oss << "MultRosen_D" << dims << ".txt";
			thisFilename = oss.str();
		}
		outputFileStart(thisFilename);
		// output headings to file
		ofstream os(thisFilename.c_str(), ios::app);         // append
		if (os.is_open()) {
			os << "l'=10^\tA1\tA2\tA3\tA3-A2\tA3/A2\tmd2\tmd3" << endl;
			
			os.close();
		}
		else {
			std::cerr << "Error: could not open file named "
				<< thisFilename << std::endl << std::endl;
		}
		
		for (size_t iPower = 0; iPower < power_vec.size(); ++iPower) {
			
			size_t maxLeavesPrime = 10;
			int lPower = power_vec[iPower];
			for (int i = 1; i < lPower; ++i) maxLeavesPrime*= 10;
			
			size_t maxLeaves = 2 * maxLeavesPrime;
				
	
			RosenFobj realF;
				
			
			cout << "\nmaxLeavesPrime = " << maxLeavesPrime<< endl;
			
			cxsc::real area1;
			cxsc::real area2;
			cxsc::real area3;
			cxsc::real maxDiam2;
			cxsc::real maxDiam3;
			
			{
			
				FunctionEstimatorInterval fei(pavingBox, realF);
						
				fei.prioritySplit(maxLeaves, logging);
				
				area1 = fei.getTotalAreaOfIntervalBand();
				
				cout << "Pq split to " << maxLeaves 
					<< " gives getTotalAreaOfIntervalBand() = " << area1 << endl;
					
				fei.hullPropagation();
				
				fei.priorityMerge(maxLeavesPrime, logging);
				
				area2 = fei.getTotalAreaOfIntervalBand();
				cout << "After hull propagate and pullup to " << maxLeavesPrime 
					<< " getTotalAreaOfIntervalBand() = " << area2 << endl;
					
				maxDiam2 = fei.getMaximumIntervalDiameter();
				
				cout << "and max diam is " << maxDiam2 << endl;
			}
			{
		
				FunctionEstimatorInterval feiShadow(pavingBox, realF);
						
				feiShadow.prioritySplit(maxLeavesPrime, logging);
				
				area3 = feiShadow.getTotalAreaOfIntervalBand();
				
				cout << "In contrast, Pq split to " << maxLeavesPrime
					<< " gives getTotalAreaOfIntervalBand() = " << area3 << endl;
					
				maxDiam3 = feiShadow.getMaximumIntervalDiameter();
				
				cout << "and max diam is " << maxDiam3 << endl;
					
				cout << "Gain is " << (area3 - area2) << endl;
				cout << "Gain ratio is " << (area3/area2) << endl;
				cout << endl;
			}
				
			// output to file
			ofstream os(thisFilename.c_str(), ios::app);         // append
			if (os.is_open()) {
				os << lPower << "\t"
					<< _double(area1) << "\t"
					<< _double(area2) << "\t"
					<< _double(area3) << "\t"
					<< _double(area3-area2) << "\t"
					<< _double(area3/area2) << "\t"
					<< _double(maxDiam2) << "\t"
					<< _double(maxDiam3) << endl;
				os.close();
			}
			else {
				std::cerr << "Error: could not open file named "
					<< thisFilename << std::endl << std::endl;
			}
			
		} // end max leaves loop
		
		//no more in this file
		
	} // end dims loop
	return 0;

}



/*! \file
\brief MappedSPnode example 2-d, with hull and pullup and marginals.

As used for example some papers. 


*/

#include "LevyDensityFobj2D.hpp"

#include "functionestimator_interval.hpp"
#include "functionestimator_real.hpp"

#include <ostream>
#include <fstream>

using namespace std;
using namespace subpavings;



int main(int argc, char* argv[])
{

    cout << "\nLevy Density\n" << endl;

    int dims = 2; // Levy can only do 2D
	double intside = 10.0;
	
	if (argc > 1) intside = atof(argv[1]);
	
    cxsc::ivector pavingBox(dims);
    cxsc::interval pavingInterval(-intside,intside);
    for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;
	
	LevyDensityFobj2D realF;
	

	int prec = 5; // for estimator output precision

	/* Creating a function estimate using a priority queue, hull propagate,
	 * and merge up, using the default (Reimann) priority measure.*/
	{
		cout << "\n\npriority queue and pull up with interval estimator,"
			<< " and make real estimator" << endl;
		
		FunctionEstimatorInterval fei(pavingBox, realF);
				
		size_t maxLeaves = 2000;
		//for (int i = 1; i < d; ++i) maxLeaves*=maxLeaves;
		
		LOGGING_LEVEL logging = NOLOG;
		
		fei.prioritySplit(maxLeaves, logging);
		
		cxsc::real area1 = fei.getTotalAreaOfIntervalBand();
		cout << "Pq split to " << maxLeaves 
			<< " gives getTotalAreaOfIntervalBand() = " << area1 << endl;
			
		
		maxLeaves = maxLeaves*7/10; // interval division
		
		fei.hullPropagation();
		
		fei.priorityMerge(maxLeaves, logging);
		
		cxsc::real area2 = fei.getTotalAreaOfIntervalBand();
		cout << "After hull propagate and pullup to " 
			<< maxLeaves << " getTotalAreaOfIntervalBand() = " 
			<< area2 << endl;
	
		{
			
			FunctionEstimatorInterval feiShadow(pavingBox, realF);
				
			LOGGING_LEVEL logging = NOLOG;
			
			feiShadow.prioritySplit(maxLeaves, logging);
			
			cxsc::real areaShadow = feiShadow.getTotalAreaOfIntervalBand();
			cout << "In contrast a q split to " << maxLeaves 
				<< " gives getTotalAreaOfIntervalBand() = " 
				<< areaShadow << endl;
			
		}
		
		cout << "\nmake a piecewise constant function" << endl;
		PiecewiseConstantFunction pcf = fei.makePiecewiseConstantFunction();
		
		real unnormArea = pcf.getTotalIntegral();
		cout << "area under the range for the pcf is " << unnormArea << endl;
						
		{
			ostringstream oss;
			oss << "Levy" << dims << "D_RealFromPQPullUp_l" 
				<< maxLeaves << "_db" << intside << ".txt";
			string s(oss.str());
			pcf.outputToTxtTabs(s, prec, true);
		}
		
		
		cout << "\nnormalise pcf estimate" << endl;
		PiecewiseConstantFunction pcfn = pcf.makeNormalised();
		real normArea = pcfn.getTotalIntegral();
		cout << "area under the range for the normalised pcf is " << normArea << endl;
						
		{
			ostringstream oss;
			oss << "Levy" << dims << "D_RealFromPQPullUp_l" 
				<< maxLeaves << "_db" << intside << "Norm.txt";
			string s(oss.str());
			pcfn.outputToTxtTabs(s, prec, true);
		}
		
		cout << "\nmarginalise the normalised pcf" << endl;
		
		for (int md = 1; md <= dims; ++md) {
			cout << "\nmarginalise on dimension " << md << endl;
			std::vector<int> reqDims;
			reqDims.push_back(md);
			PiecewiseConstantFunction pcfnm = pcfn.makeMarginal(reqDims);
			real pcfnMargArea = pcfnm.getTotalIntegral();
			cout << "area under the range for the marginalised normalised real estimator is " 
					<< pcfnMargArea << endl;
							
			{
				ostringstream oss;
				oss << "Levy" << dims << "D_RealFromPQPullUp_l" 
					<< maxLeaves << "_db" << intside << "Norm_m" << md << ".txt";
				string s(oss.str());
				pcfnm.outputToTxtTabs(s, prec, true);
			}
		}
		
		cout << "\ncoverage regions for the normalised pcf" << endl;
		int precCov = 10;
		{
			cxsc:: real cov = 1.0; // coverage 1 needed for making the plots
			ostringstream oss;
			oss << "Levy" << dims << "D_RealCoverageAll.txt";
			string s(oss.str());
			pcfn.outputCoverageRegion(s, cov, precCov, true);
		}
		
		{
			cxsc:: real cov = 0.5;
			cout << "\ncoverage for " << cov << endl;
			
			ostringstream oss;
			oss << "Levy" << dims << "D_RealCoverageRegion_c" 
				<<_double(cov) << ".txt";
			string s(oss.str());
			pcfn.outputCoverageRegion(s, cov, precCov, true);
			
		}
		
		{
			cxsc:: real cov = 0.9;
			cout << "\ncoverage for " << cov << endl;
			
			ostringstream oss;
			oss << "Levy" << dims 
				<< "D_RealCoverageRegion_c" <<_double(cov) << ".txt";
			string s(oss.str());
			pcfn.outputCoverageRegion(s, cov, precCov, true);
			
		}
		{
			cxsc:: real cov = 0.1;
			cout << "\ncoverage for " << cov << endl;
			
			ostringstream oss;
			oss << "Levy" << dims 
				<< "D_RealCoverageRegion_c" <<_double(cov) << ".txt";
			string s(oss.str());
			pcfn.outputCoverageRegion(s, cov, precCov, true);
			
		}
			
	}
    return 0;

}


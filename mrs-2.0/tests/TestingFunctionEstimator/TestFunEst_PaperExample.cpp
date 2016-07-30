/*
* Copyright (C) 2011, 2012 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


/*! \file
\brief Testing function estimation using oscillator function
* used for example of pq split and pull up (there and back again) 
* in various papers.
 */

#include "functionestimator_interval.hpp"
#include "functionestimator_real.hpp"
#include "piecewise_constant_function.hpp"
#include "intervalmappedspnode_measurers.hpp"
#include "oscFobj1.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>

using namespace cxsc;
using namespace std;
using namespace subpavings;



		
int main()
{
		int prec = 5; // default precision for output files
		
		interval simpleFunctionDomainInterval(0.0,1.0);
		interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
		interval overlapSimpleFunctionDomainInterval(Inf(simpleFunctionDomainInterval)-0.5*Mid(simpleFunctionDomainInterval), Mid(simpleFunctionDomainInterval));
						
		
		
		#define PQHULLPAPEREXAMPLE
		#ifdef PQHULLPAPEREXAMPLE
		try {
			cout << "\npriority queue with hull and priority merge estimators for article example" << endl;
			
			// dimensions - note fewer
			for (int d = 1; d < 2; ++d) {	
				
				ivector pavingBox(d);
				interval pavingInterval(0.5, 1.0);
				for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
				
				OscFobj fobj;

				FunctionEstimatorInterval fei(pavingBox, fobj);
				
				size_t maxLeaves = 50;
				for (int i = 1; i < d; ++i) maxLeaves*=maxLeaves;
				
				
				#if(0)
					LOGGING_LEVEL loggingSplit = NOLOG;
				#endif
				#if(1)
					LOGGING_LEVEL loggingSplit = LOGSAMPLES;
				#endif
				
				// pq on area (note - first draft used prioritySplitOnGain
				fei.prioritySplit(maxLeaves, loggingSplit);
				
				size_t leavesBefore = fei.getRootLeaves();
				
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					fei.outputToTxtTabs(s, prec, true);
				}
				
				cout << "function estimation has " << leavesBefore << " leaves" << endl;
				cxsc::real areaBefore = fei.getTotalAreaOfIntervalBand();
				cout << "getTotalAreaOfIntervalBand() = " << areaBefore << endl;
			
				
				maxLeaves = maxLeaves*9/10; // interval division
				
				cout << "\nsee what would happen if we had gone directly for " << maxLeaves << " leaves" << endl;
				 
				FunctionEstimatorInterval feiShadow(pavingBox, fobj);
				
				//note first draft used prioritySplitOnGain;
			
				LOGGING_LEVEL loggingSplitAgain = NOLOG;
				
				feiShadow.prioritySplit(maxLeaves, loggingSplitAgain);
				
				cxsc::real areaShadow = feiShadow.getTotalAreaOfIntervalBand();
				cout << "in that case, getTotalAreaOfIntervalBand() = " << areaShadow << endl;
			
				cout << "benefit from extra leaves is = " 
						<< areaShadow << " - " << areaBefore << " = "
						<< (areaShadow - areaBefore) << endl;
				
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					feiShadow.outputToTxtTabs(s, prec, true);
				}
				
				cout << "\n about to do hull propagation and priority merge" << endl;
				fei.hullPropagation();
				
				cout << "\nPriority merge to " << maxLeaves << " leaves" << endl;
				
				clock_t start = clock();
				
				#if(1)
					LOGGING_LEVEL loggingMerge = NOLOG;
				#endif
				#if(0)
					LOGGING_LEVEL loggingMerge = LOGSAMPLES;
				#endif
				// merge on area - note first draft used priorityMergeOnLoss
				fei.priorityMerge(maxLeaves, loggingMerge);
				
				// stop recording time here
				clock_t end = clock();	
				double timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
				cout << "Computing time for merge up in estimate: " << timing << " s."<< endl;
	
				
								
				size_t leavesAfter = fei.getRootLeaves();
				assert(leavesBefore >= leavesAfter);
				cxsc::real areaAfter = fei.getTotalAreaOfIntervalBand();
				cout << "After propagation and priority merge, function estimation has " << leavesAfter << " leaves" << endl;
				cout << "getTotalAreaOfIntervalBand() = " << areaAfter << endl;
				cout << "gain from going down and pulling up to " << maxLeaves << "rather than going straight there is " 
						<< areaShadow << " - " << areaAfter << " = "
						<< (areaShadow - areaAfter) << endl;
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorAfterPMD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					fei.outputToTxtTabs(s, prec, true);
				}
				
			}
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\npriority queue with hull and priority merge estimators:\n" << msg << endl;
			throw;
		}
		#endif
		
		#define PQHULLGIFEXAMPLE
		#ifdef PQHULLGIFEXAMPLE
		try {
			cout << "\npriority queue with hull and priority merge estimators for gif" << endl;
			
			// dimensions - note fewer
			for (int d = 1; d < 2; ++d) {	
			//for (int d = 1; d < maxDims; ++d) {	
				
				ivector pavingBox(d);
				interval pavingInterval(0.0, 1.0);
				for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
				
				OscFobj fobj;

				FunctionEstimatorInterval fei(pavingBox, fobj);
				
				size_t maxLeaves = 100;
				for (int i = 1; i < d; ++i) maxLeaves*=maxLeaves;
				
				#if(0)
					LOGGING_LEVEL loggingSplit = NOLOG;
				#endif
				#if(1)
					LOGGING_LEVEL loggingSplit = LOGSAMPLES;
				#endif
				
				fei.prioritySplit(maxLeaves, loggingSplit);
				
				size_t leavesBefore = fei.getRootLeaves();
				
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					fei.outputToTxtTabs(s, prec, true);
				}
				
				cout << "function estimation has " << leavesBefore << " leaves" << endl;
				cxsc::real areaBefore = fei.getTotalAreaOfIntervalBand();
				cout << "getTotalAreaOfIntervalBand() = " << areaBefore << endl;
			
				
				maxLeaves = maxLeaves*5/10; // interval division
				
				cout << "\nsee what would happen if we had gone directly for " << maxLeaves << " leaves" << endl;
				 
				FunctionEstimatorInterval feiShadow(pavingBox, fobj);
				
				//note first draft used prioritySplitOnGain;
				#if(1)
					LOGGING_LEVEL loggingSplitAgain = NOLOG;
				#endif
				#if(0)
					LOGGING_LEVEL loggingSplitAgain = LOGSAMPLES;
				#endif
				
				feiShadow.prioritySplit(maxLeaves, loggingSplitAgain);
				
				cxsc::real areaShadow = feiShadow.getTotalAreaOfIntervalBand();
				cout << "in that case, getTotalAreaOfIntervalBand() = " << areaShadow << endl;
			
				cout << "benefit from extra leaves is = " 
						<< areaShadow << " - " << areaBefore << " = "
						<< (areaShadow - areaBefore) << endl;
				
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					feiShadow.outputToTxtTabs(s, prec, true);
				}
				
				cout << "\n about to do hull propagation and priority merge" << endl;
				fei.hullPropagation();
				
				cout << "\nPriority merge to " << maxLeaves << " leaves" << endl;
				
				clock_t start = clock();
				
				#if(0)
					LOGGING_LEVEL loggingMerge = NOLOG;
				#endif
				#if(1)
					LOGGING_LEVEL loggingMerge = LOGSAMPLES;
				#endif
				// merge on area - note first draft used priorityMergeOnLoss
				fei.priorityMerge(maxLeaves, loggingMerge);
				
				// stop recording time here
				clock_t end = clock();	
				double timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
				cout << "Computing time for merge up in estimate: " << timing << " s."<< endl;
	
				
								
				size_t leavesAfter = fei.getRootLeaves();
				assert(leavesBefore >= leavesAfter);
				cxsc::real areaAfter = fei.getTotalAreaOfIntervalBand();
				cout << "After propagation and priority merge, function estimation has " << leavesAfter << " leaves" << endl;
				cout << "getTotalAreaOfIntervalBand() = " << areaAfter << endl;
				cout << "gain from going down and pulling up to " << maxLeaves << "rather than going straight there is " 
						<< areaShadow << " - " << areaAfter << " = "
						<< (areaShadow - areaAfter) << endl;
				{
					ostringstream oss;
					oss << "PQEstimatorOscillatorAfterPMD" << d << "_l" << maxLeaves << ".txt";
					string s(oss.str());
					fei.outputToTxtTabs(s, prec, true);
				}
				
			}
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\npriority queue with hull and priority merge estimators:\n" << msg << endl;
			throw;
		}
		#endif
		
		
	cout << "\nEnd test\n" << endl;

    return 0;

} // end of test program




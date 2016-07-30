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
\brief Testing function estimation using FunctionEstimatorInterval
* and FunctionEstimatorReal, including using these to make a 
* PiecewiseConstantFunction
 */

#include "functionestimator_interval.hpp"
#include "functionestimator_real.hpp"
#include "piecewise_constant_function.hpp"
#include "intervalmappedspnode_measurers.hpp"
#include "sp_check_visitor.hpp"
#include "fei_evalobj.hpp"
#include "spnode.hpp"
#include "subpaving_exception.hpp"
#include "simpleFobj1.hpp"
#include "simpleFobj2.hpp"
#include "oscFobj1.hpp"

#include "cxsc.hpp"

#include <fstream> 
#include <sstream>  
#include <ostream>  
#include <cassert>

using namespace cxsc;
using namespace std;
using namespace subpavings;


void testBasic();
void testBruteForceEstimation(); 
void testBruteForceHullPropagation(); 
void testPQEstimationGain(); 
void testIntelligentBruteForceEstimation(); 
void testComparePQEstimates(); 
void testHullPropagationReimann(); 
void testMakePCFFromFEI(); 
void testMakePCFFromFER(); 


		
int main()
{


	testBasic();
	testBruteForceEstimation();
	testBruteForceHullPropagation(); 
	testPQEstimationGain(); 
	testIntelligentBruteForceEstimation(); 
	testHullPropagationReimann(); 
	testMakePCFFromFEI(); 
	testMakePCFFromFER(); 


	cout << "\nEnd test\n" << endl;

    return 0;

} // end of test program


void testBasic() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	interval overlapSimpleFunctionDomainInterval(Inf(simpleFunctionDomainInterval)-0.5*Mid(simpleFunctionDomainInterval), Mid(simpleFunctionDomainInterval));
					
	
	try {
		cout << "\nTry the simple function" << endl;
		
		SimpleFobj2 fobj;
		
		for (int d = 1; d < 4; ++d) {	
			cout << "\n\n---- For dimension = " << d << " ----\n"<< endl;
			{
				rvector insidePt(d);
				for(int k=1; k <= d; k++) insidePt[k] = mid(simpleFunctionDomainInterval);
				cout << "Point inside domain box is " << insidePt; 
				
				cxsc::real realInsidePtImage = fobj(insidePt);
				
				cout << "Real image of point inside domain box is " << realInsidePtImage << endl;
			}
			{
				rvector partInsidePt(d);
				partInsidePt[1] = Sup(simpleFunctionDomainInterval)+ 1.0;
				for(int k=2; k <= d; k++) partInsidePt[k] = mid(simpleFunctionDomainInterval);
				cout << "\nPoint only partly inside domain box is " << partInsidePt; 
				
				cxsc::real realPartInsidePtImage = fobj(partInsidePt);
				
				cout << "Real image of point only partly inside domain box is " << realPartInsidePtImage << endl;
			}
			{
				ivector box(d);
				for(int k=1; k <= d; k++) box[k] = simpleFunctionDomainInterval;
				cout << "\nFunction domain box is " << box; 
				
				cxsc::interval intervalImage = fobj(box);
				
				cout << "Interval image of whole domain box is " << intervalImage  << endl;
				
				cxsc::real midImage = fobj.imageMid(box);
				
				cout << "Real mid-image of whole domain box is " << midImage  << endl;
			}
			{
				ivector partBox(d);
				for(int k=1; k <= d; k++) partBox[k] = halfSimpleFunctionDomainInterval;
				cout << "\nPart of domain box is " << partBox; 
				
				cxsc::interval intervalPartImage = fobj(partBox);
				
				cout << "Interval image of part of domain box is " << intervalPartImage  << endl;
	
				cxsc::real midImagePartBox = fobj.imageMid(partBox);
				
				cout << "Real mid-image of part of domain box is " << midImagePartBox  << endl;
			}
			{

				ivector overlapBox(d);
				overlapBox[1] = overlapSimpleFunctionDomainInterval;
				for(int k=2; k <= d; k++) overlapBox[k] = simpleFunctionDomainInterval;
				cout << "\nBox overlapping edge of domain is " << overlapBox; 
				
				cxsc::interval intervalOverlapImage = fobj(overlapBox);
				
				cout << "Interval image of box overlapping edge of domain is " << intervalOverlapImage  << endl;
		
				cxsc::real midImageOverlapBox = fobj.imageMid(overlapBox);
				
				cout << "Real mid-image of box overlapping edge of domain is " << midImageOverlapBox  << endl;
			}
		}
	}
	catch (std::exception& ee) {
		cout << "\nFailed to test simple function:\n" << ee.what() << endl;
		throw;
	}
	
	try {
		cout << "\nNode constructed with default box" << endl;
		
		cxsc::ivector uselessBox;

		SimpleFobj2 fobj;

		FunctionEstimatorInterval temp(uselessBox, fobj);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct FunctionEstimatorInterval with default box:\n" << ee.what() << endl;
	}
	
	try {
		cout << "\nNode constructed with default box" << endl;
		
		cxsc::ivector uselessBox;

		SimpleFobj2 fobj;

		FunctionEstimatorReal temp(uselessBox, fobj);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (subpavings::MalconstructedBox_Error& ee) {
		cout << "\nFailed to construct FunctionEstimatorReal with default box:\n" << ee.what() << endl;
	}
	
	int maxDims = 3;
	
	try {
		cout << "\nBasic estimators, using intervals" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			ostringstream oss;

			oss << "BasicEstimatorIntervalD" << d << ".txt";
			
			string s(oss.str());
			
			fei.outputToTxtTabs(s, prec, true);
			
			cout << "The diameter of the range of the root box is ";
			cout << fei.getRootRangeDiameter() << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators using intervals:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nBasic estimators, using reals" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorReal fer(pavingBox, fobj);
			
			ostringstream oss;

			oss << "BasicEstimatorRealD" << d << ".txt";
			
			string s(oss.str());
			
			fer.outputToTxtTabs(s, prec, true);
			
			cout << "The area under the range over the paving is ";
			cout << fer.getTotalIntegralOfRealEstimate() << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators using reals:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest copy constructor for interval estimators" << endl;
		
		int d = 1;
		ivector pavingBox(d);
		ivector partBox(d);
		rvector insidePt(d);
		interval pavingInterval = simpleFunctionDomainInterval;
		for(int k=1; k <= d; k++) {
			pavingBox[k] = pavingInterval;
			partBox[k] = halfSimpleFunctionDomainInterval;
			insidePt[k] = mid(simpleFunctionDomainInterval);
		}
		
		cout << "Point inside domain box is " << insidePt; 
		cout << "\nPart of domain box is " << partBox; 
		
		
		SimpleFobj1 fobj1;
		
		FunctionEstimatorInterval fei1(pavingBox, fobj1);
		
		string fei1Summary = fei1.stringSummary();
		cout << "Summary of fei1  is:\t" << fei1Summary << endl;
		ostringstream oss1;
		fei1.outputRootToStreamTabs(oss1, prec);
		string fei1Root = oss1.str();
		cout << "fei1 root is:\t" << fei1Root << endl;
		
		cout << "For fei1:" << endl;
		cxsc::real fei1RealInsidePtImage = fei1.getFobjReference()(insidePt);
		cout << "Real image of point inside domain box is " << fei1RealInsidePtImage << endl;
		cxsc::interval fei1IntervalPartImage = fei1.getFobjReference()(partBox);
		cout << "Interval image of part of domain box is " << fei1IntervalPartImage  << endl;
		cxsc::real fei1MidImagePartBox = fei1.getFobjReference().imageMid(partBox);
		cout << "Real mid-image of part of domain box is " << fei1MidImagePartBox  << endl;
	
		cout << "Construct fei2 using copy constructor from fei1:" << endl;
		
		FunctionEstimatorInterval fei2(fei1);
		string fei2SummaryBefore = fei2.stringSummary();
		
		cout << "Summary of fei2 is:\t" << fei2SummaryBefore << endl;
		ostringstream oss2;
		fei2.outputRootToStreamTabs(oss2, prec);
		string fei2Root = oss2.str();
		cout << "fei2 root is:\t" << fei2Root << endl;
		
		assert (fei1Root == fei2Root);
		
		cout << "For fei2:" << endl;
		cxsc::real fei2RealInsidePtImage = fei2.getFobjReference()(insidePt);
		cout << "Real image of point inside domain box is " << fei2RealInsidePtImage << endl;
		assert (fei1RealInsidePtImage == fei2RealInsidePtImage);
		cxsc::interval fei2IntervalPartImage = fei2.getFobjReference()(partBox);
		cout << "Interval image of part of domain box is " << fei2IntervalPartImage  << endl;
		assert (fei1IntervalPartImage == fei2IntervalPartImage);
		cxsc::real fei2MidImagePartBox = fei2.getFobjReference().imageMid(partBox);
		cout << "Real mid-image of part of domain box is " << fei2MidImagePartBox  << endl;
		assert (fei1MidImagePartBox == fei2MidImagePartBox);
		
		cout << "\nPassed all asserts" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do test copy constructor:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nTest copy constructor for real estimators" << endl;
		
		int d = 1;
		ivector pavingBox(d);
		ivector partBox(d);
		rvector insidePt(d);
		interval pavingInterval = simpleFunctionDomainInterval;
		for(int k=1; k <= d; k++) {
			pavingBox[k] = pavingInterval;
			partBox[k] = halfSimpleFunctionDomainInterval;
			insidePt[k] = mid(simpleFunctionDomainInterval);
		}
		
		SimpleFobj2 fobj1;
		
		FunctionEstimatorReal fer1(pavingBox, fobj1);
		
		string fer1Summary = fer1.stringSummary();
		cout << "Summary of fer1  is:\t" << fer1Summary << endl;
		ostringstream oss1;
		fer1.outputRootToStreamTabs(oss1, prec);
		string fer1Root = oss1.str();
		cout << "fer1 root is:\t" << fer1Root << endl;
		
		cout << "Construct fer2 using copy constructor from fer1:" << endl;
		
		FunctionEstimatorReal fer2(fer1);
		string fer2SummaryBefore = fer2.stringSummary();
		
		cout << "Summary of fer2 is:\t" << fer2SummaryBefore << endl;
		ostringstream oss2;
		fer2.outputRootToStreamTabs(oss2, prec);
		string fer2Root = oss2.str();
		cout << "fer2 root is:\t" << fer2Root << endl;
		
		assert (fer1Root == fer2Root);
		
		cout << "\nPassed all asserts" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do test copy constructor:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct interval estimator from spn subpaving with no box (this should fail)" << endl;
		
		SimpleFobj1 fobj;
		
		SPnode spn;

		FunctionEstimatorInterval fei(spn, fobj);
		
		throw std::logic_error("Should not be able to get here");
		
	}
	catch (subpavings::NoBox_Error& nbe) {
		cout << "\nException:\n" << nbe.what() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct from spn test:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nConstruct real estimator from spn subpaving with no box (this should fail)" << endl;
		
		SimpleFobj1 fobj;
		
		SPnode spn;

		FunctionEstimatorReal fer(spn, fobj);
		
		throw std::logic_error("Should not be able to get here");
		
	}
	catch (subpavings::NoBox_Error& nbe) {
		cout << "\nException:\n" << nbe.what() << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed construct from spn test:\n" << msg << endl;
		throw;
	}
	
	{ //  two blocks below use the string vectors and the shape
		std::vector < std::string > splitToShapeStringsIntervals;
		std::vector < std::string > splitToShapeStringsReals;
		string shape("2,3,3,1");
				
		try {
			cout << "\nConstruct from spn subpaving" << endl;
			
			// dimensions
			for (int d = 1; d < maxDims+1; ++d) {	
				
				ivector pavingBox(d);
				interval pavingInterval = simpleFunctionDomainInterval;
				for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
				
				SimpleFobj2 fobj;
				
				SPnode spn(pavingBox);
				spn.splitRootToShape(shape);
				
				cout << "Address of spn to be used for constructor is " << (&spn) << endl;
				
				{
					cout << "\nConstructing interval estimators using " << (&spn) << endl;
				
					FunctionEstimatorInterval fei(spn, fobj);
					
					ostringstream oss;
					fei.outputToStreamTabs(oss);
					splitToShapeStringsIntervals.push_back(oss.str());
		
					cout << "String summary is " << fei.stringSummary();
					cout << "The number of leaves is ";
					cout << fei.getRootLeaves();
					cout << " and the interval band area is ";
					cout << fei.getTotalAreaOfIntervalBand() << endl;
				}
				{
					cout << "\nConstructing real estimators using " << (&spn) << endl;
				
					FunctionEstimatorReal fer(spn, fobj);
					
					ostringstream oss;
					fer.outputToStreamTabs(oss);
					splitToShapeStringsReals.push_back(oss.str());
		
					cout << "String summary is " << fer.stringSummary();
					cout << "The number of leaves is ";
					cout << fer.getRootLeaves();
					cout << " and the area under the function estimate is ";
					cout << fer.getTotalIntegralOfRealEstimate() << endl;
				}
				
			}
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to construct from spn subpaving:\n" << msg << endl;
			throw;
		}
		
		
		try {
			cout << "\nSplit to shape" << endl;
			
			// dimensions
			for (int d = 1; d < maxDims+1; ++d) {	
				
				ivector pavingBox(d);
				interval pavingInterval = simpleFunctionDomainInterval;
				for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
				
				SimpleFobj2 fobj;

				{
					cout << "\nwith interval estimators ";
					FunctionEstimatorInterval fei(pavingBox, fobj);
					
					fei.splitToShape(shape);
		
					{
						ostringstream oss;
						fei.outputToStreamTabs(oss);
						
						assert(oss.str() == splitToShapeStringsIntervals.at(d-1));
						cout << "Passed assert that result is same as constructing with spn split to shape" << endl;
					}
		
					{
						ostringstream oss;
						oss << "SplitToShapeIntervalD" << d << ".txt";
						string s(oss.str());
						fei.outputToTxtTabs(s, prec, true);
					}
					
					cout << "The number of leaves is ";
					cout << fei.getRootLeaves();
					cout << " and the interval band area is ";
					cout << fei.getTotalAreaOfIntervalBand() << endl;
				}
				{
					cout << "\nwith real estimators ";
					FunctionEstimatorReal fer(pavingBox, fobj);
					
					fer.splitToShape(shape);
		
					{
						ostringstream oss;
						fer.outputToStreamTabs(oss);
						
						assert(oss.str() == splitToShapeStringsReals.at(d-1));
						cout << "Passed assert that result is same as constructing with spn split to shape" << endl;
					}
		
					{
						ostringstream oss;
						oss << "SplitToShapeRealD" << d << ".txt";
						string s(oss.str());
						fer.outputToTxtTabs(s, prec, true);
					}
					
					cout << "The number of leaves is ";
					cout << fer.getRootLeaves();
					cout << " and the area under the real estimate is ";
					cout << fer.getTotalIntegralOfRealEstimate() << endl;
				}
			}
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do split to shape:\n" << msg << endl;
			throw;
		}
	}
}

void testBruteForceEstimation() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	
	int maxDims = 3;
		
	std::vector<size_t> leafNumbers(maxDims,1);
	
	try {
		cout << "\nbrute force estimators, intervals" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			real tolerance = 0.5;
			fei.bruteForceEstimate(tolerance);
			
			ostringstream oss;

			oss << "BruteForceEstimatorIntervalD" << d << "_t" << _double(tolerance) << ".txt";
			
			string s(oss.str());
			
			fei.outputToTxtTabs(s, prec, true);
			
			size_t leaves = fei.getRootLeaves();
			leafNumbers[d-1] = leaves;
			cxsc::real area = fei.getTotalAreaOfIntervalBand();
			cout << "function estimation has " << leaves << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << area << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators, intervals:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nbrute force estimators, reals" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorReal fer(pavingBox, fobj);
			
			real tolerance = 0.5;
			fer.bruteForceEstimate(tolerance);
			
			ostringstream oss;

			oss << "BruteForceEstimatorRealD" << d << "_t" << _double(tolerance) << ".txt";
			
			string s(oss.str());
			
			fer.outputToTxtTabs(s, prec, true);
			
			size_t leaves = fer.getRootLeaves();
			cxsc::real area = fer.getTotalIntegralOfRealEstimate();
			cout << "function estimation has " << leaves << " leaves" << endl;
			cout << "getTotalIntegralOfRealEstimate() = " << area << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators, reals:\n" << msg << endl;
		throw;
	}
}

void testBruteForceHullPropagation() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	
	int maxDims = 3;
		
	try {
		cout << "\nHull propagation with brute force estimators on an oscillating function" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			OscFobj fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			real tolerance = 0.1;
			tolerance *= (d*d);
			fei.bruteForceEstimate(tolerance);
			size_t leavesBefore = fei.getRootLeaves();
			cxsc::real areaBefore = fei.getTotalAreaOfIntervalBand();
			cxsc::real rootRangeBefore = fei.getRootRangeDiameter();
			cout << "\nBefore propagation, function estimation has " << leavesBefore << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << areaBefore << endl;
			cout << "getRootRangeDiameter() = " << rootRangeBefore << endl;
			{
				ostringstream oss;
	
				oss << "BruteForceEstimatorBeforePropagationD" << d << "_t" << _double(tolerance) << ".txt";
				
				string s(oss.str());
				
				fei.outputToTxtTabs(s, prec, true);
			}
			
			fei.hullPropagation();
			
			size_t leavesAfter = fei.getRootLeaves();
			assert(leavesBefore == leavesAfter);
			cxsc::real areaAfter = fei.getTotalAreaOfIntervalBand();
			assert(areaBefore == areaAfter);
			cxsc::real rootRangeAfter = fei.getRootRangeDiameter();
			assert(rootRangeBefore >= rootRangeAfter);
			cout << "After propagation, function estimation has " << leavesAfter << " leaves (should be the same)" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << areaAfter << " (should be the same)" << endl;
			cout << "getRootRangeDiameter() = " << rootRangeAfter << endl;
			{
				ostringstream oss;
	
				oss << "BruteForceEstimatorAfterPropagationD" << d << "_t" << _double(tolerance) << ".txt";
				
				string s(oss.str());
				
				fei.outputToTxtTabs(s, prec, true);
			}
			
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do hull propagation:\n" << msg << endl;
		throw;
	}
	catch (...) {
		cout << "\nUnknown exception" << endl;
		throw;
	}
}	

/* Gain measure is an alternative to Reimann measure. It looks one step 
 * ahead to the 'gain' in tighter interval enclosure area from a split
 * -- we decided that this had many disadvantages */
void testPQEstimationGain() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	
	try {
		
		cout << "\npriority queue estimators using gain measure" << endl;
		
		// dimensions
		for (int d = 1; d < 2; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			size_t maxLeaves = 10;
			
			LOGGING_LEVEL logging = NOLOG;
			
			fei.prioritySplitOnGain(maxLeaves, logging);
			
			ostringstream oss;

			oss << "PQEstimatorGainD" << d << "_l" << maxLeaves << ".txt";
			
			string s(oss.str());
			
			fei.outputToTxtTabs(s, prec, true);
			
			cout << "function estimation has " << fei.getRootLeaves() << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << fei.getTotalAreaOfIntervalBand() << endl;
		
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators:\n" << msg << endl;
		throw;
	}
	
}

void testIntelligentBruteForceEstimation() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	
	int maxDims = 3;
		
	try {
		
		cout << "\n\nintelligent brute force estimate and pq compared" << endl;
		
		std::vector < size_t > maxLeavesVec;
		maxLeavesVec.push_back(10);
		maxLeavesVec.push_back(100);
		maxLeavesVec.push_back(1000);
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			OscFobj fobj;
			
			cxsc::real tol(0.0000001);
			

			for (std::vector < size_t >::iterator it = maxLeavesVec.begin()
						; it < maxLeavesVec.end()
						; ++it) {

				size_t maxLeaves = *it;
				
				
				cout << "\nd = " << d << "maxLeaves = " 
					<< maxLeaves << " tol = " << _double(tol) << endl;
				
				{
					cout << "\nbreadth first brute force estimate" << endl;
					
					FunctionEstimatorInterval fei(pavingBox, fobj);
					// visitor to check nodes using interval image tolerance requirement
					IntervalImageToleranceCheck nodeChecker(fobj, tol);
					
					fei.breadthFirstBruteForceEstimate(nodeChecker, maxLeaves, 1234);
					
					ostringstream oss;
		
					oss << "OscIntBFEst" << d << "_l" << maxLeaves << ".txt";
					
					string s(oss.str());
					
					
					cout << "getTotalAreaOfIntervalBand() = " << fei.getTotalAreaOfIntervalBand() << endl;
				}
				{
					cout << "\npq estimate using area" << endl;
					
					FunctionEstimatorInterval fei(pavingBox, fobj);
					
					LOGGING_LEVEL logging = NOLOG;
					
					// priority split using default Reimann measure
					fei.prioritySplit(maxLeaves, logging);
					
					ostringstream oss;
		
					oss << "OscPQAreaEst" << d << "_l" << maxLeaves << ".txt";
					
					string s(oss.str());
					
					
					cout << "getTotalAreaOfIntervalBand() = " << fei.getTotalAreaOfIntervalBand() << endl;
				}
			}// end for
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators:\n" << msg << endl;
		throw;
	}
	
}

/* comparing pq using the 'gain' measure and the Reimann ('area') measures */
void testComparePQEstimates() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
	
	int maxDims = 3;
	
	std::vector<size_t> leafNumbers(maxDims,1);
		
	try {
		cout << "\npriority queue estimators using gain" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			size_t maxLeaves = leafNumbers[d-1];
			
			LOGGING_LEVEL logging = NOLOG;
			
			fei.prioritySplitOnGain(maxLeaves, logging);
			
			ostringstream oss;

			oss << "PQEstimatorGainD" << d << "_l" << maxLeaves << ".txt";
			
			string s(oss.str());
			
			fei.outputToTxtTabs(s, prec, true);
			
			cout << "function estimation has " << fei.getRootLeaves() << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << fei.getTotalAreaOfIntervalBand() << endl;
		
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\npriority queue estimators using total area (Reimann)" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims+1; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			SimpleFobj2 fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			size_t maxLeaves = leafNumbers[d-1];

			LOGGING_LEVEL logging = NOLOG;
			
			// default Reimann measure
			fei.prioritySplit(maxLeaves, logging);
			
			ostringstream oss;

			oss << "PQEstimatorAreaD" << d << "_l" << maxLeaves << ".txt";
			
			string s(oss.str());
			
			fei.outputToTxtTabs(s, prec, true);
			
			cout << "function estimation has " << fei.getRootLeaves() << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << fei.getTotalAreaOfIntervalBand() << endl;
		
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do basic estimators:\n" << msg << endl;
		throw;
	}
}

/* split, hull propagate and merge, using Reimann measure */
void testHullPropagationReimann() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
		
	try {
		cout << "\npriority queue (Reiman priority function) with hull and priority merge estimators" << endl;
		
		// dimensions - note fewer
		for (int d = 1; d < 2; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			OscFobj fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			size_t maxLeaves = 20000;
			for (int i = 1; i < d; ++i) maxLeaves*=maxLeaves;
			
			LOGGING_LEVEL logging = NOLOG;
			
			// default Reimann measure
			fei.prioritySplit(maxLeaves, logging);
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
		
			
			maxLeaves = maxLeaves*2/10; // interval division
			
			cout << "\nsee what would happen if we had gone directly for " << maxLeaves << " leaves" << endl;
			 
			FunctionEstimatorInterval feiShadow(pavingBox, fobj);
			
			feiShadow.prioritySplit(maxLeaves, logging);
			cxsc::real areaShadow = feiShadow.getTotalAreaOfIntervalBand();
			cout << "in that case, getTotalAreaOfIntervalBand() = " << areaShadow << endl;
		
			
			{
				ostringstream oss;
				oss << "PQEstimatorOscillatorD" << d << "_l" << maxLeaves << ".txt";
				string s(oss.str());
				feiShadow.outputToTxtTabs(s, prec, true);
			}
			
			cout << "\n about to do hull propagation and priority merge" << endl;
			fei.hullPropagation();
			
			FEICritLeaves_LTE pmLeaves(maxLeaves);
			cout << "\nPriority merge to " << maxLeaves << " leaves" << endl;
			
			// merge (also using Reimann measure)
			clock_t start = clock();
			fei.priorityMerge(maxLeaves, logging);
			
			// stop recording time here
			clock_t end = clock();	
			double timing = ((static_cast<double>(end - start)) / CLOCKS_PER_SEC);
			cout << "Computing time for merge up in estimate: " << timing << " s."<< endl;

			
							
			size_t leavesAfter = fei.getRootLeaves();
			assert(leavesBefore >= leavesAfter);
			cxsc::real areaAfter = fei.getTotalAreaOfIntervalBand();
			cout << "After propagation and priority merge, function estimation has " << leavesAfter << " leaves" << endl;
			cout << "getTotalAreaOfIntervalBand() = " << areaAfter << endl;
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
}

void testMakePCFFromFEI() 
{
	int prec = 5; // default precision for output files
	
	interval simpleFunctionDomainInterval(0.0,1.0);
	interval halfSimpleFunctionDomainInterval(Mid(simpleFunctionDomainInterval), Sup(simpleFunctionDomainInterval));
		
	int maxDims = 3;
	
	try {
		cout << "\nMake a piecewise constant function after priority queue pull up" << endl;
		
		// dimensions
		for (int d = 1; d < maxDims; ++d) {	
			
			ivector pavingBox(d);
			interval pavingInterval = simpleFunctionDomainInterval;
			for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
			
			OscFobj fobj;

			FunctionEstimatorInterval fei(pavingBox, fobj);
			
			size_t maxLeaves = 100;
			for (int i = 1; i < d; ++i) maxLeaves*=maxLeaves;
			
			LOGGING_LEVEL logging = NOLOG;
			
			//pq on using Reimann measure
			fei.prioritySplit(maxLeaves, logging);
			
			maxLeaves = maxLeaves*9/10; // interval division
			
			cout << "\nabout to do hull propagation and priority merge" << endl;
			fei.hullPropagation();
			
			cout << "\nPriority merge to " << maxLeaves << " leaves" << endl;
							
			fei.priorityMerge(maxLeaves, logging);
			size_t feiLeaves = fei.getRootLeaves();
			cout << "number of leaves in fei is " << feiLeaves << endl;
			cout << "string summary is" << endl;
			cout << fei.stringSummary() << endl;
			
			cout << "\nabout to make real estimator fer" << endl;
			
			
			PiecewiseConstantFunction pcf = fei.makePiecewiseConstantFunction();
			size_t pcfLeaves = pcf.getRootLeaves();
			assert(feiLeaves == pcfLeaves);
			cout << "number of leaves in pcf is " << pcfLeaves << endl;
			cout << "string summary is" << endl;
			cout << pcf.stringSummary() << endl;
							
			{
				ostringstream oss;
				oss << "PCFOscillatorFromFei" << d << "_l" << pcfLeaves << ".txt";
				string s(oss.str());
				pcf.outputToTxtTabs(s, prec, true);
			}
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make a piecewise constant function after priority queue pull up:\n" << msg << endl;
		throw;
	}
}

void testMakePCFFromFER() 
{
	int prec = 5; // default precision for output files

	try {
		
		cout << "\n\nMake a piecwise constant function from a function estimator real:\n" << endl;
		
		
		int d = 2; // dimension of the box to sample data from
		ivector pavingBox(d);
		interval pavingInterval(0,1);
		for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

		SimpleFobj2 fobj;

		FunctionEstimatorReal fer1(pavingBox, fobj);
		
		std::string split = "2,3,4,4,2,3,3";
		
		try {
			fer1.splitToShape(split);
			cout << "String summmary of FunctionEstimatorReal is: " << endl;
			cout << fer1.stringSummary() << endl;
			string s1 = "splitToShapeFer1.txt";
			fer1.outputToTxtTabs(s1, prec, true);
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do splitToShape:\n" << msg << endl;
		}
		
		try {
			PiecewiseConstantFunction pcf = fer1.makePiecewiseConstantFunction();
			
			cout << "String summmary of PiecewiseConstantFunction is: " << endl;
			cout << pcf.stringSummary() << endl;
			
			string s1 = "PCFfromFer1.txt";
			pcf.outputToTxtTabs(s1, prec, true);
							
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to .makePiecewiseConstantFunction:\n" << msg << endl;
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make a piecwise constant function from a function estimator real:\n" << msg << endl;
		throw;
	}
}		




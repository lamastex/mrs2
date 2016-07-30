/*
* Copyright (C) 2011 Jennifer Harlow
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
\brief Testing AdaptiveHistogramCollators.

Run the executable and then use the shell script testing_collhists.sh
to run checks out output.

 */

#include "testing_tools.hpp"
#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data

#include <iostream>
#include <iterator>
#include <fstream>  // input and output streams
#include <stdexcept>
#include <cassert> // assert
#include <cfloat> // for DBL_EPSILON

using namespace cxsc;
using namespace std;
using namespace subpavings;

void outputADHC(const std::string& s, const AdaptiveHistogramCollator& adhc, const int prec = 5);
void swapCheckOutput(const std::string& s, const AdaptiveHistogramCollator& adhc); 
std::vector< RealVec > makeRangesZero(cxsc::real rootvol);
std::vector < RealVec > makeRanges1(cxsc::real rootvol);		
std::vector < RealVec > makeRanges2(cxsc::real rootvol);		
std::vector < RealVec > makeRanges3(cxsc::real rootvol);		
RVecData& getData1(RVecData& data);
RVecData& getData2(RVecData& data);
RVecData& getData3(RVecData& data);
RVecData& getDataExtra1(RVecData& data);
RVecData& getDataExtra2(RVecData& data);

		
int main()
{
	try {
			cout << "\nDefault constructor" << endl;
			
			AdaptiveHistogramCollator temp;
			
			string s1 = "defaultConstructorCollatorHist.txt";
			outputADHC(s1, temp, 10);
			
			assert(temp.isEmptyCollation());
			cout << "Passed assert that collator is empty" << endl;
			assert( checkFileLines(s1, 0) );
			cout << "Passed assert that output file is empty" << endl;
			
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default constructed AdaptiveHistogramCollator:\n" << msg << endl;
		throw;
	}
	
	try {
			cout << "\nCopy constructor on default" << endl;
			
			AdaptiveHistogramCollator temp1;
			AdaptiveHistogramCollator temp2(temp1);
			
			string s = "copyOfdefaultConstructor.txt";
			outputADHC(s, temp2, 10);
			
			assert(temp2.getNumberCollated() == 0);
			cout << "Passed assert that collator is empty" << endl;
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
				
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction of default constructed AdaptiveHistogramCollator:\n" << msg << endl;
			throw;
	}
	
	AdaptiveHistogramCollator nothingCollatedColl;
	
	try {
		cout << "\nAssignment copy of default constructor" << endl;
		
		AdaptiveHistogramCollator temp1;
		AdaptiveHistogramCollator temp2 = temp1;
		
		string s1 = "assignmentCopyOfdefaultConstructor.txt";
		outputADHC(s1, temp2, 10);
		
		assert(temp2.isEmptyCollation());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s1, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
		nothingCollatedColl = temp1;
		string s2 = "nothingCollatedColl.txt";
		outputADHC(s2, nothingCollatedColl, 10);
			
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do assignment copy of default constructed AdaptiveHistogramCollator:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\naddition with default constructed collator and another default default collator" << endl;
		
		AdaptiveHistogramCollator temp1;
		AdaptiveHistogramCollator temp2;
		
		AdaptiveHistogramCollator temp3 = temp1 + temp2;
		string s = "additionDefaultAndDefault.txt";
		outputADHC(s, temp3, 10);
		
		assert(temp3.isEmptyCollation());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addition with default constructed collator and another default collator:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nconstructor with default histogram (this should fail)" << endl;
		
		AdaptiveHistogram temp1;
		AdaptiveHistogramCollator temp2(temp1);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do constructor with default histogram:\n" << msg << endl;
	}
	
	int d = 2; // dimension of the box to sample data from
    ivector pavingBox(d);
    interval pavingInterval(-2,2);
    for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	
	try {
		cout << "\nconstructor with histogram with box but no data" << endl;
		
		AdaptiveHistogram temp1(pavingBox);
		AdaptiveHistogramCollator temp2(temp1);
		
		string s = "defaultConstructorWithHistogramWithBoxNoData.txt";
		outputADHC(s, temp2, 10);
		
		assert(temp2.getNumberCollated() == 1);
		cout << "Passed assert that collator has one element" << endl;
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
			
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do constructor with histogram with box but no data:\n" << msg << endl;
		throw;
	}
	
	
	try {
		cout << "\naddToCollation with default constructed collator and default histogram" << endl;
		
		AdaptiveHistogramCollator temp1;
		AdaptiveHistogram adh;
		
		temp1.addToCollation(adh);
		string s = "addToCollationDefaultAndDefaultHist.txt";
		outputADHC(s, temp1, 10);
		
		assert(temp1.isEmptyCollation());
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addToCollation with default constructed collator and default histogram:\n" << msg << endl;
		throw;
	}
	
	
	AdaptiveHistogramCollator zeroValueColl;

	try {
		cout << "\naddToCollation with default constructed collator and histogram with box but no data" << endl;
		
		AdaptiveHistogramCollator temp1;
		AdaptiveHistogram adh(pavingBox);
		
		temp1.addToCollation(adh);
		string s1 = "addToCollationDefaultAndNoDataHist.txt";
		outputADHC(s1, temp1, 10);
		
		assert(temp1.getNumberCollated() == 1);
		
		zeroValueColl = temp1;
		string s2 = "zeroValueColl.txt";
		outputADHC(s2, zeroValueColl, 10);
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addToCollation with default constructed collator and histogram with box but no data:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\naddToCollation for collator of histogram with box but no data, and another hist with wrong dimensions (this should fail)" << endl;
		
		AdaptiveHistogram adh1(pavingBox);
		
		AdaptiveHistogramCollator temp(adh1);
		
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;

		AdaptiveHistogram adh2(pavingBox2);
		
		temp.addToCollation(adh2);
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addToCollation for collator of histogram with box but no data, and another hist with wrong dimensions:\n" << msg << endl;
	}
	
	try {
		cout << "\naddToCollation for collator of histogram with box but no data, and another hist with wrong box sides (this should fail)" << endl;
		
		AdaptiveHistogram adh1(pavingBox);
		
		AdaptiveHistogramCollator temp(adh1);
		
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d);
		pavingBox2[1] = pavingInterval;
		interval pavingInterval2(-2,3);
		pavingBox2[2] = pavingInterval2;

		AdaptiveHistogram adh2(pavingBox2);
		
		temp.addToCollation(adh2);
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addToCollation for collator of histogram with box but no data, and another hist with wrong box sides:\n" << msg << endl;
	}
	

	// make some histograms and then use to make collators

	std::string split1 = "3,4,4,2,2,4,4,3";
	std::string split2 = "2,3,4,4,2,3,3";
	std::string split3 = "1,2,3,3";
	
	AdaptiveHistogram hist1(pavingBox);
    hist1.splitToShape(split1);

    AdaptiveHistogram hist2(pavingBox);
    hist2.splitToShape(split2);

    AdaptiveHistogram hist3(pavingBox);
    hist3.splitToShape(split3);

    // put in the data in a 'pulse' with no further splitting
	RVecData data1;
	data1 = getData1(data1);
	bool successfulInsertionHistFirst = hist1.insertFromRVec(data1);
    if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
	
	RVecData data2;
	data2 = getData2(data2);
	bool successfulInsertionHistSecond = hist2.insertFromRVec(data2);
    if (!successfulInsertionHistSecond) cout << "unsuccessful insertion 2" << endl;

    RVecData data3;
	data3 = getData3(data3);
	bool successfulInsertionHistThird = hist3.insertFromRVec(data3);
    if (!successfulInsertionHistThird) cout << "unsuccessful insertion 3" << endl;

	assert((successfulInsertionHistFirst && successfulInsertionHistSecond 
					&& successfulInsertionHistThird) == true);
	
	AdaptiveHistogramCollator collHist1;
	AdaptiveHistogramCollator collHist2;
	AdaptiveHistogramCollator collHist3;
					
    try {
		cout << "\nIndividual collator constructors with histograms" << endl;
		
		{
			assert(collHist1.isEmptyCollation());
			AdaptiveHistogramCollator temp(hist1);
			collHist1 = temp;
			string s = "collHist1.txt";
			outputADHC(s, collHist1, 10);
			assert(collHist1.getNumberCollated() == 1);
		}
		{
			assert(collHist2.isEmptyCollation());
			AdaptiveHistogramCollator temp(hist2);
			collHist2 = temp;
			string s = "collHist2.txt";
			outputADHC(s, collHist2, 10);
			assert(collHist2.getNumberCollated() == 1);
		}
		{
			assert(collHist3.isEmptyCollation());
			AdaptiveHistogramCollator temp(hist3);
			collHist3 = temp;
			string s = "collHist3.txt";
			outputADHC(s, collHist3, 10);
			assert(collHist3.getNumberCollated() == 1);
		}
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct collators with individual histograms:\n" << msg << endl;
		throw;
	}
	
	AdaptiveHistogramCollator twoHists;
	try {
		cout << "\nCollator constructed by constructing from hist1 and then += collator for hist2" << endl;
		
		AdaptiveHistogramCollator temp(hist1);
		temp += collHist2;
		assert(temp.getNumberCollated() == 2);
		
		string s1 = "twoHistsCollator.txt";
		outputADHC(s1, temp, 10);
		
		twoHists = temp;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collator constructed by adding all three hists to collation:\n" << msg << endl;
		throw;
	} 	
	AdaptiveHistogramCollator threeHists;
	try {
		cout << "\nCollator constructed by adding all three hists to collation" << endl;
		
		AdaptiveHistogramCollator temp;
		
		temp.addToCollation(hist1);
		temp.addToCollation(hist2);
		temp.addToCollation(hist3);
		
		assert(temp.getNumberCollated() == 3);
		
		string s1 = "threeHistsCollator.txt";
		outputADHC(s1, temp, 10);
		
		threeHists = temp;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collator constructed by adding all three hists to collation:\n" << msg << endl;
		throw;
	} 	
	
	try {
		cout << "\nCollator constructed by addition of individual collations" << endl;
		
		AdaptiveHistogramCollator temp = collHist1 + collHist2 + collHist3;
		
		assert(temp.getNumberCollated() == 3);
		
		string s1 = "threeHistsCollatorFromAddition.txt";
		outputADHC(s1, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collator constructed by addition of individual collations:\n" << msg << endl;
		throw;
	} 	
	
	try {
		cout << "\nCollator constructed by += addition of individual collations" << endl;
		
		AdaptiveHistogramCollator temp;
		AdaptiveHistogramCollator temp1(collHist1);
		AdaptiveHistogramCollator temp2(collHist2);
		AdaptiveHistogramCollator temp3(collHist3);
		
		temp += temp1 += temp2 += temp3;
		
		assert(temp.getNumberCollated() == 3);
		
		string s1 = "threeHistsCollatorFromPlusEqualsAddition.txt";
		outputADHC(s1, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collator constructed by +=addition of individual collations:\n" << msg << endl;
		throw;
	}	 	
 		
	try {
		cout << "\nCopy assignment of threeHistsCollator" << endl;
		
		AdaptiveHistogramCollator temp = threeHists;
		
		assert(temp.getNumberCollated() == 3);
		
		string s = "copyAssignmentThreeHists.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy assignment of threeHistsCollator:\n" << msg << endl;
		throw;
	}	
	try {
		cout << "\nCopy construction from threeHistsCollator" << endl;
		
		AdaptiveHistogramCollator temp(threeHists);
		
		assert(temp.getNumberCollated() == 3);
		
		string s = "copyConstructorFromThreeHists.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy construction from threeHistsCollator:\n" << msg << endl;
		throw;
	}	
	
	
	try {
		cout << "\nswap of copy of threeHistsCollator and copy of twoHistsCollator" << endl;
		
		AdaptiveHistogramCollator temp1(threeHists);
		AdaptiveHistogramCollator temp2(twoHists);
		
		string s11 = "swapCheckCopyOfThreeHistsBeforeSwap.txt";
		swapCheckOutput(s11, temp1);
		string s12 = "swapCheckCopyOfTwoHistsBeforeSwap.txt";
		swapCheckOutput(s12, temp2);
		
		std::swap(temp1, temp2);
		
		string s21 = "swapCheckShouldBeLikeCopyOfThreeHistsAfterSwap.txt";
		swapCheckOutput(s21, temp2);
		string s22 = "swapCheckShouldBeLikeCopyOfTwoHistsAfterSwap.txt";
		swapCheckOutput(s22, temp1);
		
		string s1 = "shouldBeCopyOfTwoHists_afterSwap.txt";
		outputADHC(s1, temp1, 10);
		string s2 = "shouldBeCopyOfThreeHists_afterSwap.txt";
		outputADHC(s2, temp2, 10);
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do std::swap of copy of threeHistsCollator and copy of twoHistsCollator:\n" << msg << endl;
		throw;
	}
	
	
	try {
		cout << "\ncopy of threeHistCollator += collator with nothing collated" << endl;
		
		AdaptiveHistogramCollator temp(threeHists);
		
		temp += nothingCollatedColl;
		
		assert(temp.getNumberCollated() == 3);
		
		string s = "shouldBeCopyOfThreeHists_afterPlusEqualsEmptyCollator.txt";
		outputADHC(s, temp, 10);
				
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of threeHistCollator += collator with nothing collated:\n" << msg << endl;
		throw;
	}
	
	AdaptiveHistogramCollator threeHistsPlusZeroValueColl;
	try {
		cout << "\ncopy of threeHistCollator += collator with zero value" << endl;
		
		AdaptiveHistogramCollator temp(threeHists);
		
		temp += zeroValueColl;
		
		assert(temp.getNumberCollated() == 4);
		
		threeHistsPlusZeroValueColl = temp;
		string s = "threeHistsPlusZeroValueCollator.txt";
		outputADHC(s, temp, 10);
				
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of threeHistCollator += collator with zero value:\n" << msg << endl;
		throw;
	}
		
	// check normalise with empty collators
	
	try {
		cout << "\ncall makeNormalised on an empty collator (this should fail)" << endl;
		
		AdaptiveHistogramCollator temp;
		AdaptiveHistogramCollator norm = temp.makeNormalised();
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeNormalised on an empty collator:\n" << msg << endl;
	}
	
	try {
		cout << "\nmakeNormalised with zero data collator (this should fail since no normaliser)" << endl;
		
		AdaptiveHistogramCollator temp = zeroValueColl.makeNormalised();
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeNormalised with zero data collator:\n" << msg << endl;
	}

	// check average with empty collators
	
	try {
		cout << "\ncall makeAverage on an empty collator  (this should fail)" << endl;
		
		AdaptiveHistogramCollator temp;
		AdaptiveHistogramCollator av = temp.makeAverage();
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeAverage on an empty collator:\n" << msg << endl;
	}
	
	try {
		cout << "\nmakeAverage with zero value collator (this should work)" << endl;
		
		AdaptiveHistogramCollator temp = zeroValueColl.makeAverage();
		
		string s = "averageFromZeroValueColl.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeAverage with zero data collator:\n" << msg << endl;
		throw;
	}
	
	// check marginalise with empty nodes and children
	
	std::vector<int> reqDims;
	reqDims.push_back(1);
	
	try {
		cout << "\ncall makeMarginal on nothingCollated collator (this should fail)" << endl;
		
		AdaptiveHistogramCollator temp;
		AdaptiveHistogramCollator m = temp.makeMarginal(reqDims);
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeMarginal on nothingCollated collator:\n" << msg << endl;
	}
	
	try {
		cout << "\ncall makeMarginal on zeroValueColl (this should work)" << endl;
		
		AdaptiveHistogramCollator temp = zeroValueColl.makeMarginal(reqDims);
		
		string s = "zeroValueCollMakeMarginal.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeMarginal on zeroValueColl:\n" << msg << endl;
		throw;
	}
	
	// averaging and normalisation
	
	AdaptiveHistogramCollator averageOfThreeHists;
	try {
		cout << "\nmakeAverage with threeHistsColl" << endl;
		
		AdaptiveHistogramCollator temp = threeHists.makeAverage();
				
		assert(temp.getNumberCollated() == 1);
		
		averageOfThreeHists = temp;
		
		string s = "averageFromThreeHistsColl.txt";
		outputADHC(s, averageOfThreeHists, 10);
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeAverage from threeHistsColl:\n" << msg << endl;
		throw;
	}
	
	AdaptiveHistogramCollator normalisationOfThreeHists;
	try {
		cout << "\nmakeNormalised with threeHistsColl" << endl;
		
		AdaptiveHistogramCollator temp = threeHists.makeNormalised();
		
		assert(temp.getNumberCollated() == 1);
		
		normalisationOfThreeHists = temp;
		
		string s = "normalisationFromThreeHistsColl.txt";
		outputADHC(s, normalisationOfThreeHists, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeNormalised threeHistsColl:\n" << msg << endl;
		throw;
	}

	try {
		cout << "\nmake average from copy of threeHists plus nothingCollatedColl" << endl;
		
		AdaptiveHistogramCollator temp1(threeHists);
		temp1+= nothingCollatedColl;
		assert(temp1.getNumberCollated() == 3);
		AdaptiveHistogramCollator temp2 = temp1.makeAverage();
		
		string s = "averageThreeHistsPlusNothingCollatedColl.txt";
		outputADHC(s, temp2, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make average from threeHistsPlusZeroValueColl:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nmake average from threeHistsPlusZeroValueColl" << endl;
		
		AdaptiveHistogramCollator temp = threeHistsPlusZeroValueColl.makeAverage();
		
		string s = "averageThreeHistsPlusZeroValueColl.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make average from threeHistsPlusZeroValueColl:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nmake normalisation from threeHistsPlusZeroValueColl" << endl;
		
		AdaptiveHistogramCollator temp = threeHistsPlusZeroValueColl.makeNormalised();
		
		string s = "normalisationThreeHistsPlusZeroValueColl.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make normalisation from threeHistsPlusZeroValueColl:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nmake average from collation of four of collHist1" << endl;
		
		AdaptiveHistogramCollator temp = (collHist1 + collHist1 + collHist1 + collHist1).makeAverage();
		
		string s = "averageOfCollationOfFourOfCollHist1.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make average from collation of four of collHist1:\n" << msg << endl;
		throw;
	}

	
	//marginalise dimensions
	try {
		cout << "\nmakeMarginal on collHist1 with no dims required (this should fail)" << endl;
		
		std::vector<int> r;
		//reqDims.push_back(1);
	
		AdaptiveHistogramCollator temp = collHist1.makeMarginal(r);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to makeMarginal on collHist1 with no dims required:\n" << msg << endl;
	}
	try {
		cout << "\nmakeMarginal on collHist1 with req dims not present (this should fail)" << endl;
		
		std::vector<int> r;
		r.push_back(1);
		r.push_back(3);
	
		AdaptiveHistogramCollator temp = collHist1.makeMarginal(r);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to makeMarginal on collhist1 with req dims not present:\n" << msg << endl;
	}
	
	try {
		cout << "\nmakeMarginal on collHist1 with all dims required" << endl;
		
		std::vector<int> r;
		r.push_back(1);
		r.push_back(2);
	
		AdaptiveHistogramCollator temp = collHist1.makeMarginal(r);
		string s = "makeMarginalCollHist1AllDimsReq.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to makeMarginal on collHist1 with all dims required:\n" << msg << endl;
		throw;
	}
	
	
	try {
		cout << "\nmakeMarginal on collhist1 on dim 1" << endl;
		
		std::vector<int> r;
		r.push_back(1);
		
		AdaptiveHistogramCollator temp = collHist1.makeMarginal(r);
		string s = "makeMarginalCollHist1OnDim1.txt";
		outputADHC(s, temp, 10);
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to makeMarginal on coll1 on dim 1:\n" << msg << endl;
		throw;
	}
	
	try {
		cout << "\nmakeMarginal on collhist1 on dim 2" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		
		AdaptiveHistogramCollator temp = collHist1.makeMarginal(r);
		string s = "makeMarginalCollHist1OnDim2.txt";
		outputADHC(s, temp, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to makeMarginal on coll1 on dim 2:\n" << msg << endl;
		throw;
	} 
	
	try {
		cout << "\nadd collhist1 and makeMarginal on collhist1 on dim 1 (this should fail)" << endl;
		
		std::vector<int> r;
		r.push_back(1);
		
		AdaptiveHistogramCollator temp1 = collHist1.makeMarginal(r);
		AdaptiveHistogramCollator temp2 = temp1 + collHist1;
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to add collhist1 and makeMarginal on collhist1 on dim 1:\n" << msg << endl;
	} 
	try {
		cout << "\nadd hist1 to collation made by makeMarginal on collhist1 on dim 1 (this should fail)" << endl;
		
		std::vector<int> r;
		r.push_back(1);
		
		AdaptiveHistogramCollator temp1 = collHist1.makeMarginal(r);
		temp1.addToCollation(hist1);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to add hist1 to collation made by makeMarginal on collhist1 on dim 1:\n" << msg << endl;
	} 
	AdaptiveHistogramCollator makeMarginalThreeHistsOnDim2;
	try {
		cout << "\nmakeMarginalOnThreeHists" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		
		AdaptiveHistogramCollator temp = threeHists.makeMarginal(r);
		
		string s = "makeMarginalOnThreeHistsOnDim2.txt";
		outputADHC(s, temp, 10);
		
		makeMarginalThreeHistsOnDim2 = temp;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeMarginalOnThreeHists:\n" << msg << endl;
		throw;
	}
	AdaptiveHistogramCollator makeMarginalOnAverageOfThreeHistsOnDim2;
	try {
		cout << "\nmakeMarginalOnAverageThreeHists" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		
		AdaptiveHistogramCollator temp = averageOfThreeHists.makeMarginal(r);
		
		string s = "makeMarginalOnAverageOfThreeHistsOnDim2.txt";
		outputADHC(s, temp, 10);
		
		makeMarginalOnAverageOfThreeHistsOnDim2 = temp;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeMarginalOnAverageThreeHists:\n" << msg << endl;
		throw;
	}
	AdaptiveHistogramCollator collationMarginalsOfThreeHistsOnDim2;
	try {
		cout << "\ncollation of of marginals of collHist1, collHist2, collHist3 on dim2" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		AdaptiveHistogramCollator temp1 = collHist1.makeMarginal(r);
		AdaptiveHistogramCollator temp2 = collHist2.makeMarginal(r);
		AdaptiveHistogramCollator temp3 = collHist3.makeMarginal(r);
		
		AdaptiveHistogramCollator temp = temp1 + temp2 + temp3;
		
		string s = "collationMarginalsOfThreeHistsOnDim2.txt";
		outputADHC(s, temp, 10);
		collationMarginalsOfThreeHistsOnDim2 = temp;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collation of of marginals of collHist1, collHist2, collHist3 on dim2:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nmakeAverage on collation of of marginals of collHist1, collHist2, collHist3 on dim2" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		
		AdaptiveHistogramCollator temp = collationMarginalsOfThreeHistsOnDim2.makeAverage();
		
		string s = "averageOfCollationMarginalsOfThreeHistsOnDim2.txt";
		outputADHC(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeAverage on collation of of marginals of collHist1, collHist2, collHist3 on dim2:\n" << msg << endl;
		throw;
	}
	try {
		cout << "\nAverage makeMarginalOnThreeHists" << endl;
		
		AdaptiveHistogramCollator temp = makeMarginalThreeHistsOnDim2.makeAverage();
		
		string s = "averageOfMakeMarginalThreeHistsOnDim2.txt";
		outputADHC(s, temp, 10);
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do average of makeMarginalThreeHists:\n" << msg << endl;
		throw;
	}
	
// findCoverage and findEmpiricalDensity	
	{
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		try {
			
			cout << "\nCoverage with collator with nothing collated (should fail)" << endl;
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
	
			double cov = nothingCollatedColl.findCoverage(pt);
			
			throw std::logic_error("Should not be able to do this");
		
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage with collator with nothing collated:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nDensity with collator with nothing collated (should fail)" << endl;
			
			double ed = nothingCollatedColl.findEmpiricalDensity(pt);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density with with collator with nothing collated:\n" << msg << endl;
		}
	}
	{
		cxsc::rvector pt(3);
		pt[1] = -6;
		pt[2] = 0;
		pt[2] = 1;
			
			
		try {
		
			cout << "\nCoverage with collHist1 with 3-d point (should fail)" << endl;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = collHist1.findCoverage(pt);
			
			throw std::logic_error("Should not be able to do this");
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage for collHist1 with 3d point:\n" << msg << endl;
		}
		
		try {
		
			cout << "\nDensity with hist1 with 3-d point (should fail)" << endl;
			
			double ed = collHist1.findEmpiricalDensity(pt);
			
			throw std::logic_error("Should not be able to do this");
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density for collHist1 with 3d point:\n" << msg << endl;
		}
	}
	{
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		try {
			
			cout << "\nCoverage with collator with zero value" << endl;
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
	
			double cov = zeroValueColl.findCoverage(pt);
			double covShouldBe = 0.0;
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
		
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage with collator with zero value:\n" << msg << endl;
			throw;
		}
		
		try {
		
			cout << "\nDensity with collator with zero value" << endl;
			
			double ed = zeroValueColl.findEmpiricalDensity(pt);
			double edShouldBe = 0.0;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density with collator with zero value:\n" << msg << endl;
			throw;
		}
	}
	{
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = -3;
		
		try {
			
			cout << "\nCoverage with collHist1 for point not in root collator" << endl;
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
	
			double cov = collHist1.findCoverage(pt);
			double covShouldBe = 0.0;
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
		
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do Coverage with collHist1 for point not in root collator:\n" << msg << endl;
			throw;
		}
		
		try {
		
			cout << "\nDensity with collHist1 for point not in root collator" << endl;
			
			double ed = collHist1.findEmpiricalDensity(pt);
			double edShouldBe = 0.0;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do density with collHist1 for point not in root collator:\n" << msg << endl;
			throw;
		}
	}
	
	try {
			
		cout << "\nCoverage and density for point in collHist1 (XLR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -2;
		pt[2] = -2; // in XLLL
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (0.0);
		double edShouldBe = (0.0);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XLLL):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XLR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -2;
		pt[2] = 2; // in XLR
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (0.0);
		double edShouldBe = (0.0);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XLR):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRL)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 2;
		pt[2] = -2; // in XRL
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25 - 0.125 - 0.0625*2);
		double edShouldBe = (0.03125);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRL):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRRR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 2;
		pt[2] = 2; // in XRRR
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25 - 0.125);
		double edShouldBe = (0.06250);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRRR):\n" << msg << endl;
		throw;
	}


	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRL)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 2;
		pt[2] = -2; // in XRL
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25 - 0.125 - 0.0625*2);
		double edShouldBe = (0.03125);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRL):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRRLL)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0; // in XRRLL
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0);
		double edShouldBe = (0.375);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRRLL):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XLLRR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -1;
		pt[2] = -1; // in XLLRR
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375);
		double edShouldBe = (0.25);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XLLRR):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRRLR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 1; // in XRRLR
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25);
		double edShouldBe = (0.125);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRRLR):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRRR)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 1;
		pt[2] = 0; // in XRRR
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25 - 0.125);
		double edShouldBe = (0.0625);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRRR):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XRL)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = -2; // in XRL
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
	
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (1.0 - 0.375 - 0.25 - 0.125 - 0.0625*2);
		double edShouldBe = (0.03125);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XRL):\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for point in collHist1 (XLLRL)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -1;
		pt[2] = (-1.0*(1.0+DBL_EPSILON)); // in XLLRL
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,16);
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		cout << cxsc::RestoreOpt;
			
		double cov = collHist1.findCoverage(pt);
		double ed = collHist1.findEmpiricalDensity(pt);
		
		double covShouldBe = (0.0);
		double edShouldBe = (0.0);
		
		assert(cov == covShouldBe);
		cout << "Passed assert that coverage is " << covShouldBe << endl;
		assert(ed == edShouldBe);
		cout << "Passed assert that empirical density is " << edShouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for point in collHist1 (XLLRL):\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nCoverage and density for collator hist with 3 leaves all with same density" << endl;
		
		AdaptiveHistogram temp(pavingBox);
		temp.splitToShape("2,2,1");

		// put in the data in a 'pulse' with no further splitting
		RVecData data;
		data = getDataExtra1(data);
		bool successfulInsertion = temp.insertFromRVec(data);
		if (successfulInsertion) {
			
			AdaptiveHistogramCollator coll(temp);
		
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = -2; // in XLL
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = coll.findCoverage(pt);
			double ed = coll.findEmpiricalDensity(pt);
			
			/* although pt is in XLL, which in the 
			 * ordering of the leaves by height happens
			 * to come last, coverage should use first
			 * histogram element of same height as XLL
			 * so in this case coverage is 1, ie any 
			 * chance ordering amongst leaves of the same height
			 * should not make difference to the answer*/
			
			double covShouldBe = (1.0);
			double edShouldBe = (1.0/(4*4.0));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
		}
		else cout << "data insertion unsuccessful" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for collator hist with 3 leaves all with same density:\n" << msg << endl;
		throw;
	}
	try {
	
		cout << "\nCoverage and density for collator hist with 3 leaves, smallest two with same density" << endl;
		
		AdaptiveHistogram temp(pavingBox);
		temp.splitToShape("2,2,1");

		// put in the data in a 'pulse' with no further splitting
		RVecData data;
		data = getDataExtra2(data);
		bool successfulInsertion = temp.insertFromRVec(data);
		if (successfulInsertion) {

			AdaptiveHistogramCollator coll(temp);
			
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = -2; // in XLL
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
		
			double cov = coll.findCoverage(pt);
			double ed = coll.findEmpiricalDensity(pt);
			
			/* although pt is in XLL, which in the 
			 * ordering of the leaves by height happens
			 * to come last, coverage should use first
			 * histogram element of same height as XLL,
			 * which is XR in the ordering used, 
			 * so in this case coverage only skips
			 * over the tallest node XLR (2 data points), ie any 
			 * chance ordering amongst leaves of the same height
			 * should not make difference to the answer*/
			
			
			double covShouldBe = ((5.0-2.0)/5);
			double edShouldBe = (1.0/(5*4.0));
			
			assert(cov == covShouldBe);
			cout << "Passed assert that coverage is " << covShouldBe << endl;
			assert(ed == edShouldBe);
			cout << "Passed assert that empirical density is " << edShouldBe << endl;
		}
		else cout << "data insertion unsuccessful" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do Coverage and density for collator hist with 3 leaves, smallest two with same density:\n" << msg << endl;
		throw;
	}
			
	
	
	

	// L1 distance to average of other - test tricky cases 
	
	try {
			
		cout << "\nL1 distance nothing collated coll against collHist1 (should fail)  " << endl;
		
		RealVec distances = nothingCollatedColl.getL1DistancesToAverage(collHist1);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance nothing collated coll against collHist1:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collHist1 against nothing collated coll (should fail)  " << endl;
		
		RealVec distances = collHist1.getL1DistancesToAverage(nothingCollatedColl);
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against nothing collated coll:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collHist1 against collator with box of incorrect dimensions (should fail)  " << endl;
		
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;

		AdaptiveHistogram adh(pavingBox2);
		AdaptiveHistogramCollator temp(adh);
				
		RealVec distances = collHist1.getL1DistancesToAverage(temp);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against collator with box of incorrect dimensions:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collHist1 against collator with box of incorrect side lengths (should fail)  " << endl;
		
		ivector pavingBox2 = pavingBox;
		interval pavingInterval2(-3,2);
		pavingBox2[2] = pavingInterval2;

		AdaptiveHistogram adh(pavingBox2);
		AdaptiveHistogramCollator temp(adh);
				
		RealVec distances = collHist1.getL1DistancesToAverage(temp);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against collator with box of incorrect side lengths:\n" << msg << endl;
	}
	
	
	try {
			
		cout << "\nL1 distance nothingCollatedColl against collHist1 (this should fail - no box for this)" << endl;
		
		RealVec distances (10);
		
		distances = nothingCollatedColl.getL1DistancesToAverage(distances, collHist1);
		
		throw std::logic_error("should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance nothingCollatedColl against collHist1:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collHist1 against nothingCollatedColl (this should fail)" << endl;
		
		RealVec distances (10);
		
		distances = collHist1.getL1DistancesToAverage(distances, nothingCollatedColl);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against nothingCollatedColl:\n" << msg << endl;
	}
	
	
	try {
			
		cout << "\nL1 distance zero value collator against collHist1" << endl;
		
		RealVec distances (10);
		
		distances = zeroValueColl.getL1DistancesToAverage(distances, collHist1);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(1.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance zero value collator against collHist1:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nL1 distance collHist1 against zeroDataNode  " << endl;
		
		RealVec distances;
		
		distances = collHist1.getL1DistancesToAverage(distances, zeroValueColl);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(1.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against zero value collator:\n" << msg << endl;
		throw;
	}

	// L1 distance to average of itself - test tricky cases 
	try {
			
		cout << "\nL1 distance empty collator against itself (should fail)  " << endl;
		
		const AdaptiveHistogram adh;
		const AdaptiveHistogramCollator temp (adh);
		
		RealVec distances = temp.getL1DistancesToAverage();
		
		throw std::logic_error("Should not be able to do this");
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance empty collator against itself:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance nothingCollatedColl against itself (should fail)" << endl;
		
		RealVec distances (10);
		
		distances = nothingCollatedColl.getL1DistancesToAverage(distances);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance nothingCollatedColl against itself:\n" << msg << endl;
	}
		
	try {
			
		cout << "\nL1 distance zero value against itself (should work because distance of anything to itself is 0)" << endl;
		
		RealVec distances (10);
		
		distances = zeroValueColl.getL1DistancesToAverage(distances);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(0.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance zero value coll against itself:\n" << msg << endl;
		throw;
	}
	
	// L1 distance to another node with real nodes
	
	try {
			
		cout << "\nL1 distance collHist1 against collHist1  " << endl;
		
		RealVec distances;
		
		distances = collHist1.getL1DistancesToAverage(distances, collHist1);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(0.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against collHist1:\n" << msg << endl;
		throw;
	}
	try {
			
		cout << "\nL1 distance collHist1 against itself  " << endl;
		
		RealVec distances;
		
		distances = collHist1.getL1DistancesToAverage(distances);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(0.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against itself:\n" << msg << endl;
		throw;
	}
	
	cxsc::real dis2_3(0.0);
	try {
			
		cout << "\nL1 distance collHist2 against collHist3  " << endl;
		
		RealVec distances;
		
		distances = collHist2.getL1DistancesToAverage(distances, collHist3);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(11.5/6.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
					> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		dis2_3 = distances.at(0);
		
		cout << "Passed assert that this distance dis2_3 is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist2 against collHist3:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nL1 distance collation of 2 collHist2's against collnode 3  " << endl;
		
		AdaptiveHistogramCollator temp(collHist2);
		temp += collHist2;
		
		RealVec distances;
		
		distances = temp.getL1DistancesToAverage(distances, collHist3);
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		cout << cxsc::RestoreOpt;
		
		{
			cxsc::real shouldBe(11.5/6.0);
			
			bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
			assert(!isNotEqual);
			
			cout << "Passed assert that distance of first hist is " << shouldBe << endl;
		}
		{
			cxsc::real shouldBe(11.5/6.0);
			
			bool isNotEqual = cxsc::abs(distances.at(1) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(1), shouldBe);
			assert(!isNotEqual);
			
			cout << "Passed assert that distance of second hist is " << shouldBe << endl;
		}
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collation of 2 collHist2's against collnode 3:\n" << msg << endl;
		throw;
	}
	cxsc::real dis1_2(0.0);
	
	try {
			
		cout << "\nL1 distance collHist1 against collHist2  " << endl;
		
		RealVec distances;
		
		distances = collHist1.getL1DistancesToAverage(distances, collHist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(14.0/8.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		dis1_2 = distances.at(0);
		
		cout << "Passed assert that this distance dis1_2 is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collHist1 against collHist2:\n" << msg << endl;
		throw;
	}
	RealVec distancesAgainstAverageOfThree;
	try {
			
		cout << "\nL1 distance for threeHists against its own average  " << endl;
		
		distancesAgainstAverageOfThree = threeHists.getL1DistancesToAverage(distancesAgainstAverageOfThree);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances against average are " << distancesAgainstAverageOfThree << endl;
	
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for threeHists against its own average:\n" << msg << endl;
		throw;
	}
	try {
			
		cout << "\nL1 distance for threeHists against averageOfThreeHists  " << endl;
		
		RealVec distances;
		distances = threeHists.getL1DistancesToAverage(distances, averageOfThreeHists);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances against averageOfThreeHists are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		RealVecItr ait = distancesAgainstAverageOfThree.begin();
		for (RealVecItr it = distances.begin(); it < distances.end(); ++it) {
			
			bool isNotEqual = cxsc::abs( (*it) - (*ait) ) 
				> DBL_EPSILON * cxsc::max((*it), (*ait));
			assert(!isNotEqual);
			++ait;
		}
		cout << "Passed asserts that distances are same as for threeHists against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for threeHists against averageOfThreeHists:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nL1 distance for collHist1 against threeHists  " << endl;
		
		RealVec distances = collHist1.getL1DistancesToAverage(threeHists);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs( distances.at(0) - distancesAgainstAverageOfThree.at(0) ) 
				> DBL_EPSILON * cxsc::max(distances.at(0), distancesAgainstAverageOfThree.at(0));
		assert(!isNotEqual);
		
		cout << "Passed assert that distance is same as using threeHists against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for collHist1 against threeHists:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nL1 distance for averageOfThreeHists against collHist2 " << endl;
		
		RealVec distances = averageOfThreeHists.getL1DistancesToAverage(collHist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs( distances.at(0) - distancesAgainstAverageOfThree.at(1) ) 
			> DBL_EPSILON * cxsc::max(distances.at(0), distancesAgainstAverageOfThree.at(1));
		assert(!isNotEqual);
		
		cout << "Passed assert that distance is same as using threeHists against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for averageOfThreeHists against collHist2:\n" << msg << endl;
		throw;
	}
	
	try {
			
		cout << "\nL1 distance for averageOfThree against threeHists " << endl;
		
		RealVec distances = averageOfThreeHists.getL1DistancesToAverage(threeHists);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cxsc::real shouldBe(0.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe)
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for averageOfThree against threeHists:\n" << msg << endl;
		throw;
	}
	
// test L1 distance from collator to histogram

	try {
			
		cout << "\nL1 distance from collator to histogram with null root paving (should fail) " << endl;
		
		AdaptiveHistogram adh;
		
		RealVec distances = collHist1.getL1Distances(adh);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do distance from collator to null root paving:\n" << msg << endl;
	}

	try {
			
		cout << "\nL1 distance collHist1 against hist with box of incorrect dimensions (should fail)  " << endl;
		
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval;

		AdaptiveHistogram adh(pavingBox2);
				
		RealVec distances = collHist1.getL1Distances(adh);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against hist with box of incorrect dimensions:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collHist1 against hist with box of incorrect side lengths (should fail)  " << endl;
		
		ivector pavingBox2 = pavingBox;
		interval pavingInterval2(-3,2);
		pavingBox2[2] = pavingInterval2;

		AdaptiveHistogram adh(pavingBox2);
				
		RealVec distances = collHist1.getL1Distances(adh);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against hist with box of incorrect side lengths:\n" << msg << endl;
	}
	
	
	try {
			
		cout << "\nL1 distance for threeHists against histogram with box but no data  " << endl;
		
		const AdaptiveHistogram temp(pavingBox);
		
		RealVec distances;
		distances = threeHists.getL1Distances(distances, temp);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for threeHists against histogram with box but no data are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		cxsc::real shouldBe(1.0);	
		for (RealVecItr it = distances.begin(); it < distances.end(); ++it) {
		
			bool isNotEqual = cxsc::abs( (*it) - shouldBe )
				> DBL_EPSILON * cxsc::max((*it), shouldBe);
			
			assert(!isNotEqual);
		}
		cout << "Passed asserts that distances are all " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for threeHists against histogram with box but no data:\n" << msg << endl;
		throw;
	}	

	try {
			
		cout << "\nL1 distance for collHist1 against hist1  " << endl;
		
		RealVec distances = collHist1.getL1Distances(hist1);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for collHist1 against hist1 are " << distances << endl;
		
		cxsc::real shouldBe(0.0);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is " << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for collHist1 against hist1 :\n" << msg << endl;
		throw;
	}	
	try {
			
		cout << "\nL1 distance for collHist3 against hist2  " << endl;
		
		RealVec distances = collHist3.getL1Distances(hist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for collHist1 against spsHist1 are " << distances << endl;
		
		cxsc::real shouldBe(dis2_3);
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is same as between coll nodes 2 and 3" << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for collHist3 against hist2 :\n" << msg << endl;
		throw;
	}	
	try {
			
		cout << "\nL1 distance for averageOfThreeHists against hist3  " << endl;
		
		RealVec distances = averageOfThreeHists.getL1Distances(hist3);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for averageOfThree against hist3 are " << distances << endl;
		
		cxsc::real shouldBe(distancesAgainstAverageOfThree.at(2));
		cout << "cxsc::abs(distances.at(0) - shouldBe) is " << cxsc::abs(distances.at(0) - shouldBe);
		cout << endl;
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
			> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
		assert(!isNotEqual);
		
		cout << "Passed assert that this distance is same as third element in distancesAgainstAverageOfThree" << shouldBe << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for averageOfThree against spsHist3:\n" << msg << endl;
		throw;
	}	
	try {
			
		cout << "\nL1 distance for threeHists against hist2  " << endl;
		
		RealVec distances = threeHists.getL1Distances(hist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for threeHists against hist2 are " << distances << endl;
		
		{
			cxsc::real shouldBe(dis1_2);
			bool isNotEqual = cxsc::abs(distances.at(0) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(0), shouldBe);
			assert(!isNotEqual);
		}
		{
			cxsc::real shouldBe(0.0);
			bool isNotEqual = cxsc::abs(distances.at(1) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(1), shouldBe);
			assert(!isNotEqual);
		}
		{
			cxsc::real shouldBe(dis2_3);
			bool isNotEqual = cxsc::abs(distances.at(2) - shouldBe) 
				> DBL_EPSILON * cxsc::max(distances.at(2), shouldBe);
			assert(!isNotEqual);
		}
		cout << cxsc::RestoreOpt;
		cout << "Passed assert that these distances are (dis1_2, 0.0, dis2_3)" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for threeHists against hist21:\n" << msg << endl;
		throw;
	}	

	//exporting and importing
	{
		string s("emptyCollatorExport.txt");
			
		try {
			cout << "\nexport a collator with no box" << endl;
			
			AdaptiveHistogramCollator temp;
			
			temp.exportCollator(s, 10);
			cout << "Collator exported to " << s << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to export empty collator:\n" << msg << endl;
			throw;
		}
		try {
				
			cout << "\nRe-import export of empty collator  " << endl;
			
			AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
			
			assert(import.isEmptyCollation());
			cout << "passed assert that import is empty" << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to import from export of empty collator:\n" << msg << endl;
			throw;
		}	
	}
	
		//exporting and importing
	{
		string s("zeroValueCollatorExport.txt");
		int prec = 1;
			
		try {
			cout << "\nexport a zero value collator with just root node" << endl;
			
			AdaptiveHistogramCollator temp;
			
			int d2 = 3; // dimension of the box to sample data from
			ivector pavingBox2(d2);
			for(int k=1; k <= d2; k++) {
				interval ival;
				Inf(ival) = Inf(pavingInterval)*k;
				Sup(ival) = Sup(pavingInterval)*k;
				pavingBox2[k] = ival;
			}

			AdaptiveHistogram adh2(pavingBox2);
			
			temp.addToCollation(adh2);
			
			temp.exportCollator(s, prec);
			cout << "Collator exported to " << s << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to export a zero value collator with just root node:\n" << msg << endl;
			throw;
		}
		try {
				
			cout << "\nRe-import export of zero value collator collator  " << endl;
			
			AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
			
			string reExportFile("zeroValueCollatorImportExport.txt");
			
			import.exportCollator(reExportFile, prec);
			
			cout << "Collator imported re-exported to " << reExportFile << endl;
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to import from export of zero value collator:\n" << msg << endl;
			throw;
		}	
	}
	
	try {
		
		string s("import_test_badnodelevels1.txt");
				
		cout << "\nImport from file " << s << " (this should fail)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
		throw logic_error("Should not be able to do this");
	}
	catch (logic_error& le) {
		std::string msg(le.what());
		cout << "\n*********** Failed test:\n" << msg << endl;
		cout << endl;
		throw;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nPassed test:\n" << msg << endl;
		cout << endl;
	}


	try {
		
		string s("import_test_badnodelevels2.txt");
				
		cout << "\nImport from file " << s << " (this should fail)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
		throw logic_error("Should not be able to do this");
	}
	catch (logic_error& le) {
		std::string msg(le.what());
		cout << "\n*********** Failed test:\n" << msg << endl;
		cout << endl;
		throw;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nPassed test:\n" << msg << endl;
		cout << endl;
	}


	try {
		
		string s("import_test_badbox1.txt");
				
		cout << "\nImport from file " << s << " (this should fail)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
		throw logic_error("Should not be able to do this");
	}
	catch (logic_error& le) {
		std::string msg(le.what());
		cout << "\n*********** Failed test:\n" << msg << endl;
		cout << endl;
		throw;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nPassed test:\n" << msg << endl;
		cout << endl;
	}

	try {
		
		string s("import_test_badbox2.txt");
				
		cout << "\nImport from file " << s << " (import successful, but Sup of second box is 0)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
		cout << "Box imported is " << import.getSubPaving()->getBox() << endl;
		
		assert(import.getSubPaving()->getBox() != threeHists.getSubPaving()->getBox());
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed test:\n" << msg << endl;
		cout << endl;
		throw;
	}


	try {
		
		string s("import_test_extrarange.txt");
				
		cout << "\nImport from file " << s << " (this should print a warning)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed test:\n" << msg << endl;
		cout << endl;
		throw;
	}

	try {
		
		string s("import_test_missingrange.txt");
				
		cout << "\nImport from file " << s << " (this should fail)" << endl;
		
		AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);
		
		throw logic_error("Should not be able to do this");
	}
	catch (logic_error& le) {
		std::string msg(le.what());
		cout << "\n*********** Failed test:\n" << msg << endl;
		cout << endl;
		throw;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nPassed test:\n" << msg << endl;
		cout << endl;
	}
	
	
	{
		string exportFile("threeHistsExport1.txt");
		int prec = 10;
		
		try {
				
			cout << "\nExport threeHists with same precision as original output file used " << endl;
			
			threeHists.exportCollator(exportFile, prec);
			
			cout << "Three hists exported to " << exportFile << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to export threeHists:\n" << msg << endl;
			throw;
		}	
		try {
				
			cout << "\nRe-import export of threeHists  " << endl;
			
			AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(exportFile);
			
			string s1 = "threeHistsCollatorImported1.txt";
			outputADHC(s1, import, 10);
			
			string reExportFile("threeHistsImportExport1.txt");
			
			import.exportCollator(reExportFile, prec);
			
			cout << "Three hists imported re-exported to " << reExportFile << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to import from export of threeHists:\n" << msg << endl;
			throw;
		}
		try {
	
			string s("import_test_nogap.txt");
					
			cout << "\nImport from file " << s << " (this should be okay and match three hists)" << endl;
			
			AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(s);

			string reExportFile("threeHistsImportExport3.txt");
			
			import.exportCollator(reExportFile, prec);
			
			cout << "Three hists imported re-exported to " << reExportFile << endl;

			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed test:\n" << msg << endl;
			throw;
		}	
	}
	{
		string exportFile("threeHistsExport2.txt");
		int prec = 5;
		
		try {
				
			cout << "\nExport threeHists with less precision as original output file used " << endl;
			
			threeHists.exportCollator(exportFile, prec);
			
			cout << "Three hists exported to " << exportFile << endl;
			
			
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to export threeHists:\n" << msg << endl;
			throw;
		}	
		try {
				
			cout << "\nRe-import export of threeHists  " << endl;
			
			AdaptiveHistogramCollator import = AdaptiveHistogramCollator::importCollator(exportFile);
			
			string s1 = "threeHistsCollatorImported2.txt";
			outputADHC(s1, import, 10);
			
			string reExportFile("threeHistsImportExport2.txt");
			
			import.exportCollator(reExportFile, prec);
			
			cout << "Three hists imported re-exported to " << reExportFile << endl;
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to import from export of threeHists:\n" << msg << endl;
			throw;
		}	
	}
	cout << "\nEnd of test\n" << endl;

    return 0;

} // end of test program


void outputADHC(const std::string& s, const AdaptiveHistogramCollator& adhc, const int prec)
{
	adhc.outputToTxtTabs(s, prec, false);
	cout << s << " output to file" << endl;
}

void swapCheckOutput(const std::string& s, const AdaptiveHistogramCollator& adhc) 
{
	std::ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << "This address = " << (&adhc) << endl;
		os << "address of RootCollator = " << adhc.getSubPaving() << endl;
		
		os << endl;
	
		os.close();
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


std::vector< RealVec > makeRangesZero(cxsc::real rootvol)
{
	std::vector < RealVec > result(1);
	
	result[0].push_back(cxsc::real(0.0));//X
	
	return result;
	
}

std::vector< RealVec > makeRanges1(cxsc::real rootvol)
{
	std::vector < RealVec > result(15);
	
	result[0].push_back(cxsc::real((8.0/8.0)*(1.0/rootvol)));//X
	result[1].push_back(cxsc::real((2.0/8.0)*(2.0/rootvol)));//XL
	result[2].push_back(cxsc::real((2.0/8.0)*(4.0/rootvol)));//XLL
	result[3].push_back(cxsc::real(0.0));//XLLL
	result[4].push_back(cxsc::real((2.0/8.0)*(8.0/rootvol)));//XLLR
	result[5].push_back(cxsc::real(0.0));//XLLRL
	result[6].push_back(cxsc::real((2.0/8.0)*(16.0/rootvol)));//XLLRR
	result[7].push_back(cxsc::real(0.0));//XLR
	result[8].push_back(cxsc::real((6.0/8.0)*(2.0/rootvol)));//XR
	result[9].push_back(cxsc::real((1.0/8.0)*(4.0/rootvol)));//XRL
	result[10].push_back(cxsc::real((5.0/8.0)*(4.0/rootvol)));//XRR
	result[11].push_back(cxsc::real((4.0/8.0)*(8.0/rootvol)));//XRRL
	result[12].push_back(cxsc::real((3.0/8.0)*(16.0/rootvol)));//XRRLL
	result[13].push_back(cxsc::real((1.0/8.0)*(16.0/rootvol)));//XRRLR
	result[14].push_back(cxsc::real((1.0/8.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< RealVec > makeRanges2(cxsc::real rootvol)
{
	std::vector < RealVec > result(13);
	
	result[0].push_back(cxsc::real((2.0/2.0)*(1.0/rootvol)));//X
	result[1].push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XL
	result[2].push_back(cxsc::real(0.0));//XLL
	result[3].push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XLR
	result[4].push_back(cxsc::real(0.0));//XLRL
	result[5].push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XLRR
	result[6].push_back(cxsc::real((1.0/2.0)*(16.0/rootvol)));//XLRRL
	result[7].push_back(cxsc::real(0.0));//XLRRR
	result[8].push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XR
	result[9].push_back(cxsc::real(0.0));//XRL
	result[10].push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XRR
	result[11].push_back(cxsc::real(0.0));//XRRL
	result[12].push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< RealVec > makeRanges3(cxsc::real rootvol)
{
	std::vector < RealVec > result(7);
	
	result[0].push_back(cxsc::real((3.0/3.0)*(1.0/rootvol)));//X
	result[1].push_back(cxsc::real((1.0/3.0)*(2.0/rootvol)));//XL
	result[2].push_back(cxsc::real((2.0/3.0)*(2.0/rootvol)));//XR
	result[3].push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRL
	result[4].push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRR
	result[5].push_back(cxsc::real((1.0/3.0)*(8.0/rootvol)));//XRRL
	result[6].push_back(cxsc::real(0.0));//XRRR
	
	return result;
	
}


RVecData& getData1(RVecData& data) 
{
    int d = 2;
	{ //XLLRR (-1,0), (-1,0)
		rvector thisrv(d);
		thisrv[1] = -0.5;
        thisrv[2] = -0.5;
		data.push_back(thisrv);
	}
	{ //XLLRR (-1,0), (-1,0)
		rvector thisrv(d);
		thisrv[1] = -0.75;
        thisrv[2] = -0.25;
		data.push_back(thisrv);
	}
	{ //XRL (0,2), (-2,0)
		rvector thisrv(d);
		thisrv[1] = 1.0;
        thisrv[2] = -1.0;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.5;
        thisrv[2] = 0.5;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.25;
        thisrv[2] = 0.75;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.75;
        thisrv[2] = 0.25;
		data.push_back(thisrv);
	}
	{ //XRRLR (0,1), (1,2)
		rvector thisrv(d);
		thisrv[1] = 0.5;
        thisrv[2] = 1.5;
		data.push_back(thisrv);
	}
	{ //XRRR (1,2), (0,2)
		rvector thisrv(d);
		thisrv[1] = 1.5;
        thisrv[2] = 1.0;
		data.push_back(thisrv);
	}
	return data;
}


RVecData& getData2(RVecData& data) 
{
    int d = 2;
	{ //XLRRL (-1,0), (0,1)
		rvector thisrv(d);
		thisrv[1] = -0.6;
        thisrv[2] = 0.4;
		data.push_back(thisrv);
	}
	{ //XRRR (1,2), (0,2)
		rvector thisrv(d);
		thisrv[1] = 1.3;
        thisrv[2] = 0.9;
		data.push_back(thisrv);
	}
	
	return data;
}


RVecData& getData3(RVecData& data) 
{
    int d = 2;
	{ //XL (-2,0), (-2,2)
		rvector thisrv(d);
		thisrv[1] = -1.6;
        thisrv[2] = 0.2;
		data.push_back(thisrv);
	}
	{ //XRL (0,2), (-2,0)
		rvector thisrv(d);
		thisrv[1] = 1.7;
        thisrv[2] = -0.9;
		data.push_back(thisrv);
	}
	{ //XRRL (0,1), (0,2)
		rvector thisrv(d);
		thisrv[1] = 0.8;
        thisrv[2] = 1.9;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getDataExtra1(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getDataExtra2(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -1;
        thisrv[2] = 1;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	
	return data;
}



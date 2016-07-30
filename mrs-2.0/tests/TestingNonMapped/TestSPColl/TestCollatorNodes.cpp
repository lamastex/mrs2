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
\brief Testing CollatorSPnodes

Run the executable and then use the shell script testing_colls.sh
to run checks out output.

 */

#include "testing_tools.hpp"
#include "histall.hpp"  // headers for the histograms
#include "collatorspnode.hpp"
#include "dataprep.hpp" // headers for getting data
#include "subpaving_exception.hpp"

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iterator>  
#include <cfloat> // DBL_EPSILON
#include <cassert> // assert

using namespace cxsc;
using namespace std;
using namespace subpavings;

void outputCollNode(const std::string& s, const CollatorSPnode& spn, const int prec = 5);
void outputCollAllNodes(const std::string& s, const CollatorSPnode& spn, const int prec = 5);
void swapCheckOutput(std::string& s, const CollatorSPnode& spn, int level = 0);
std::vector< RealVec > makeRangesZero(cxsc::real rootvol);
std::vector < RealVec > makeRanges1(cxsc::real rootvol);		
std::vector < RealVec > makeRanges2(cxsc::real rootvol);		
std::vector < RealVec > makeRanges3(cxsc::real rootvol);		
RVecData& getData1(RVecData& data);
RVecData& getData2(RVecData& data);
RVecData& getData3(RVecData& data);
		
int main()
{
	
	try {
			cout << "\nDefault constructor" << endl;
			
			CollatorSPnode temp;
			
			
			string s = "defaultConstructor.txt";
			outputCollNode(s, temp);
			assert(temp.isEmptyRangeCollection());
			cout << "Passed assert that collator is empty" << endl;
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
			
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default constructed CollatorSPnode:\n" << msg << endl;
	}
	
	
	try {
		cout << "\nNode constructed with default box (this should fail)" << endl;
		
		cxsc::ivector uselessBox;
		CollatorSPnode temp(uselessBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& ee) {
		cout << "\nFailed to construct NewCollatorSpnode with default box:\n" << ee.what() << endl;
	}


	try {
			cout << "\nCopy constructor on default" << endl;
			
			CollatorSPnode temp1;
			CollatorSPnode temp2(temp1);
			
			string s = "copyOfdefaultConstructor.txt";
			outputCollNode(s, temp2);
			assert(temp2.isEmptyRangeCollection());
			cout << "Passed assert that collator is empty" << endl;
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
				
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction of default constructed CollatorSPnode:\n" << msg << endl;
		}
	try {
		cout << "\nAssignment copy of default constructor" << endl;
		
		CollatorSPnode temp1;
		CollatorSPnode temp2 = temp1;
		
		string s = "assignmentCopyOfdefaultConstructor.txt";
		outputCollNode(s, temp2);
		assert(temp2.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
	
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do assignment copy of default constructed CollatorSPnode:\n" << msg << endl;
	}
	try {
		cout << "\nSplit to shape on default constructed node (this should fail)" << endl;
		
		CollatorSPnode temp;
		temp.splitRootToShape("1,1");
		
		throw std::logic_error("Should not be able to do this");

	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do splitToShape of default constructed CollatorSPnode:\n" << msg << endl;
	}
	
	try {
		cout << "\naddPavings with default constructed node and NULL" << endl;
		
		CollatorSPnode temp1;
		CollatorSPnode * temp2 = NULL;
		
		temp1.addPaving(temp2);
		string s = "addPavingDefaultAndNULL.txt";
		outputCollNode(s, temp1);
		assert(temp1.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addPavings with default constructed node and NULL:\n" << msg << endl;
	}
	
	try {
		cout << "\naddPavings with default constructed node and another default node" << endl;
		
		CollatorSPnode temp1;
		CollatorSPnode temp2;
		
		temp1.addPaving(&temp2);
		string s = "addPavingDefaultAndDefault.txt";
		outputCollNode(s, temp1);
		assert(temp1.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addPavings with default constructed node and another default node:\n" << msg << endl;
	}

	try {
		cout << "\naddition with default constructed node and another default node" << endl;
		
		CollatorSPnode temp1;
		CollatorSPnode temp2;
		
		CollatorSPnode temp3 = temp1 + temp2;
		string s = "additionDefaultAndDefault.txt";
		outputCollNode(s, temp3);
		assert(temp3.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addition with default constructed node and another default node:\n" << msg << endl;
	}
	
	int d = 2; // dimension of the box to sample data from
    ivector pavingBox(d);
    interval pavingInterval(-2,2);
    for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	try {
		cout << "\nconstructor with box" << endl;
		
		CollatorSPnode temp(pavingBox);
		string s = "constructorWithBox.txt";
		outputCollNode(s, temp, 10);
		assert(temp.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do constructor with box:\n" << msg << endl;
	}
	
	try {
		cout << "\ncopy constructor of constructor with box" << endl;
		
		CollatorSPnode temp1(pavingBox);
		CollatorSPnode temp2(temp1);
		string s = "copyConstructorOfConstructorWithBox.txt";
		outputCollNode(s, temp2, 10);
		assert(temp2.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy constructor of constructor with box:\n" << msg << endl;
	}

	CollatorSPnode nothingCollatedNode;

	try {
		cout << "\ncopy assignment of constructor with box, to make nothingCollatedNode" << endl;
		
		CollatorSPnode temp(pavingBox);
		nothingCollatedNode = temp;
		string s1 = "copyAssignmentOfConstructorWithBox.txt";
		outputCollNode(s1, nothingCollatedNode, 10);
		assert(nothingCollatedNode.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;
		string s2 = "nothingCollatedNode.txt";
		outputCollNode(s2, nothingCollatedNode, 10);
		assert( checkFileLines(s2, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy assignment of constructor with box:\n" << msg << endl;
	}
	try {
		cout << "\ntotal value x vol for nothingCollatedNode" << endl;
		
		cxsc::real valvol = nothingCollatedNode.getTotalAbsValueTimesVol();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for nothingCollatedNode is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 0.0);
		cout << "Passed assert that Val x vol == 0.0" << endl;

	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for nothingCollatedNode:\n" << msg << endl;
	}
	
	try {
		cout << "\nAdding subpavings with different number of dimensions (this should fail)" << endl;
		
		CollatorSPnode temp1(pavingBox);
		
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval2(-2,2);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval2;

		CollatorSPnode temp2(pavingBox2);
		temp1 += temp2;
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do adding subpavings with different number of dimensions:\n" << msg << endl;
	}
	
	try {
		cout << "\nAdding subpavings with different side lengths (this should fail)" << endl;
		
		CollatorSPnode temp1(pavingBox);
		
		int d2 = 2; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval2(-2,1);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval2;

		CollatorSPnode temp2(pavingBox2);
		temp1 += temp2;
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do adding subpavings with different side lengths:\n" << msg << endl;
	}

    CollatorSPnode collnode1(pavingBox);
    CollatorSPnode collnode2(pavingBox);
    CollatorSPnode collnode3(pavingBox);
    
	std::string split1 = "3,4,4,2,2,4,4,3";
	std::string split2 = "2,3,4,4,2,3,3";
	std::string split3 = "1,2,3,3";
		
	try {
		cout << "\nSplit to shape"  << endl;
		collnode1.splitRootToShape(split1);
		collnode2.splitRootToShape(split2);
		collnode3.splitRootToShape(split3);

		string s1 = "splitToShapeCollNode1.txt";
		outputCollNode(s1, collnode1);
		string s2 = "splitToShapeCollNode2.txt";
		outputCollNode(s2, collnode2);
		string s3 = "splitToShapeCollNode3.txt";
		outputCollNode(s3, collnode3);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do splitToShape:\n" << msg << endl;
	}
	
	try {
		cout << "\nAllocate ranges"  << endl;
		std::vector< RealVec > ranges = makeRanges1(collnode1.nodeRealVolume());
		collnode1.allocateRanges(ranges);
		string s = "collNode1.txt";
		outputCollNode(s, collnode1, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for collnode1:\n" << msg << endl;
	}
	try {
		std::vector< RealVec > ranges = makeRanges2(collnode2.nodeRealVolume());
		collnode2.allocateRanges(ranges);
		string s = "collNode2.txt";
		outputCollNode(s, collnode2, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for collnode2:\n" << msg << endl;
	}
	try {
		std::vector< RealVec > ranges = makeRanges3(collnode3.nodeRealVolume());
		collnode3.allocateRanges(ranges);
		string s = "collNode3.txt";
		outputCollNode(s, collnode3,10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for collnode3:\n" << msg << endl;
	}
	try {
		cout << "\nCheck total value x vol"  << endl;
		cxsc::real valvol = collnode1.getTotalAbsValueTimesVol();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for collnode1 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for collnode1:\n" << msg << endl;
	}
	try {
		cxsc::real valvol = collnode2.getTotalAbsValueTimesVol();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for collnode2 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for collnode2:\n" << msg << endl;
	}
	try {
		cxsc::real valvol = collnode3.getTotalAbsValueTimesVol();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for collnode3 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for collnode3:\n" << msg << endl;
	}
	try {
		cout << "\nCopy constructor of collnode1" << endl;
		
		CollatorSPnode cpy1(collnode1);
		
		string s = "copyConstructorCollNode1.txt";
		outputCollNode(s, cpy1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy construction of collnode1:\n" << msg << endl;
	}
	CollatorSPnode copyColl1;
	CollatorSPnode copyColl2;
	CollatorSPnode copyColl3;
	
	try {
		cout << "\nswap of copy of collnode1 and copy of collnode2" << endl;
		
		CollatorSPnode temp1(collnode1);
		CollatorSPnode temp2(collnode2);
		
		string s11 = "copyOfCollNode1BeforeSwap.txt";
		swapCheckOutput(s11, temp1);
		string s12 = "copyOfCollNode2BeforeSwap.txt";
		swapCheckOutput(s12, temp2);
		
		std::swap(temp1, temp2);
		
		string s21 = "copyOfShouldBeLikeCollNode2AfterSwap.txt";
		swapCheckOutput(s21, temp1);
		string s22 = "copyOfShouldBeLikeCollNode1AfterSwap.txt";
		swapCheckOutput(s22, temp2);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do std::swap of collnode1:\n" << msg << endl;
	}
	
	
	try {
		cout << "\nCopy assignment of collnode1" << endl;
		
		copyColl1 = collnode1;
		
		string s = "copyAssignmentCollNode1.txt";
		outputCollNode(s, copyColl1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy assignment of collnode1:\n" << msg << endl;
	}
	try {
		cout << "\nDefault copyColl2 = default copy collnode3 + collnode2" << endl;
		
		copyColl2 = copyColl3 + collnode2;
		
		string s = "shouldBeSameAsCollNode2.txt";
		outputCollNode(s, copyColl2, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default copyColl2 = default copy collnode3 + collnode2:\n" << msg << endl;
	}
	
	// make some histograms and get the SPS nodes out of them

	AdaptiveHistogram myHistFirst(pavingBox);
    myHistFirst.splitToShape(split1);

    AdaptiveHistogram myHistSecond(pavingBox);
    myHistSecond.splitToShape(split2);

    AdaptiveHistogram myHistThird(pavingBox);
    myHistThird.splitToShape(split3);

    // put in the data in a 'pulse' with no further splitting
	RVecData data1;
	data1 = getData1(data1);
	bool successfulInsertionHistFirst = myHistFirst.insertFromRVec(data1);
    if (!successfulInsertionHistFirst) cout << "unsuccessful insertion 1" << endl;
	
	RVecData data2;
	data2 = getData2(data2);
	bool successfulInsertionHistSecond = myHistSecond.insertFromRVec(data2);
    if (!successfulInsertionHistSecond) cout << "unsuccessful insertion 2" << endl;

    RVecData data3;
	data3 = getData3(data3);
	bool successfulInsertionHistThird = myHistThird.insertFromRVec(data3);
    if (!successfulInsertionHistThird) cout << "unsuccessful insertion 3" << endl;

	assert((successfulInsertionHistFirst && successfulInsertionHistSecond 
					&& successfulInsertionHistThird) == true);
					
    SPSnode spsHist1;
	SPSnode spsHist2;
	SPSnode spsHist3;
		
	try {
		cout << "\nCollator constructor with empty SPS node (this should fail)" << endl;
		
		CollatorSPnode temp(spsHist1);
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct collator with empty SPS node:\n" << msg << endl;
	}
	try {
		cout << "\nCollator constructor with SPS node with no data (this should be okay)" << endl;
		
		SPSnode temp1(pavingBox);
		
		CollatorSPnode temp2(temp1);
		assert(temp2.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;
		
		string s = "collatorFromSPSWithNoData.txt";
		outputCollNode(s, temp2, 10);
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to construct collator with SPS node with no data:\n" << msg << endl;
	} 	
 	
	try {
		cout << "\nMake SPSnodes" << endl;
		
		spsHist1 = *(myHistFirst.getSubPaving());
		spsHist2 = *(myHistSecond.getSubPaving());
		spsHist3 = *(myHistThird.getSubPaving());
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make SPSnodes:\n" << msg << endl;
	}
	
	CollatorSPnode collCopyHist1SPSnode;
	CollatorSPnode collCopyHist2SPSnode;
	CollatorSPnode collCopyHist3SPSnode;
	 
	try {
		cout << "\nMake SPSnodes with copy constructor and assignment" << endl;
		CollatorSPnode temp1(spsHist1);
		collCopyHist1SPSnode = temp1;
		CollatorSPnode temp2(spsHist2);
		collCopyHist2SPSnode = temp2;
		CollatorSPnode temp3(spsHist3);
		collCopyHist3SPSnode = temp3;
		string s1 = "sameViaSPSForCollNode1.txt";
		outputCollNode(s1, collCopyHist1SPSnode, 10);
		string s2 = "sameViaSPSForCollNode2.txt";
		outputCollNode(s2, collCopyHist2SPSnode, 10);
		string s3 = "sameViaSPSForCollNode3.txt";
		outputCollNode(s3, collCopyHist3SPSnode, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make SPSnodes:\n" << msg << endl;
	} 
	
	try {
		cout << "\nDefault += collnode3 to make copyColl3" << endl;
		
		copyColl3 += collnode3;
		
		string s = "shouldBeSameAsCollNode3.txt";
		outputCollNode(s, copyColl3, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default += collnode3:\n" << msg << endl;
	}
	
	try {
		cout << "\nnothingCollatedNode + collnode3" << endl;
		
		CollatorSPnode temp = nothingCollatedNode + collnode3;
		
		string s = "shouldBeSameAsCollNode3FromNothingCollatedPlusColl3.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do nothingCollatedNode + collnode3:\n" << msg << endl;
	}
	try {
		cout << "\ncollnode3 + nothingCollatedNode" << endl;
		
		CollatorSPnode temp = collnode3 + nothingCollatedNode;
		
		string s = "shouldBeSameAsCollNode3FromColl3PlusNothingCollated.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collnode3 + nothingCollatedNode:\n" << msg << endl;
	}
	
	try {
		cout << "\naddPavings with copy of collnode1 and NULL" << endl;
		
		CollatorSPnode * nullNode = NULL;
		copyColl1.addPaving(nullNode);
		
		string s = "addPavingsCopyOfColl1AndNULL.txt";
		outputCollNode(s, copyColl1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addPavings with copy of collnode1 and NULL:\n" << msg << endl;
	}
	CollatorSPnode add1And2;
	try {
		cout << "\naddition with collnode1 and collnode2" << endl;
		
		add1And2 = collnode1 + collnode2;
		
		string s = "additionColl1AndColl2.txt";
		outputCollNode(s, add1And2, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addition with collnode1 and collnode2:\n" << msg << endl;
	}
	CollatorSPnode anotherCopy1PlusEqualCopy2;
	try {
		cout << "\ncopy of collnode1 += collnode2" << endl;
		anotherCopy1PlusEqualCopy2 = collnode1;
		
		anotherCopy1PlusEqualCopy2 += copyColl2;
		
		assert(anotherCopy1PlusEqualCopy2.getSizeRangeCollection() == 2);
		cout << "Passed assert that collator has 2 elements" << endl;
				
		string s = "anotherCopyColl1PlusEqualColl2.txt";
		outputCollNode(s, anotherCopy1PlusEqualCopy2, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of collnode1 += collnode2:\n" << msg << endl;
	}
	
	// addition with the SPSnodes
	CollatorSPnode copyColl1ViaSPS;
	CollatorSPnode copyColl2ViaSPS;
	CollatorSPnode copyColl3ViaSPS;
	try {
		cout << "\naddPaving with default coll node and sps node equivalent of coll1" << endl;
		
		copyColl1ViaSPS.addPaving(&spsHist1);
		
		assert(copyColl1ViaSPS.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has 1 element" << endl;
		
		string s = "addPavingDefaultAndColl1ViaSPS.txt";
		outputCollNode(s, copyColl1ViaSPS, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addPaving with default coll node and sps node equivalent of coll1:\n" << msg << endl;
	}
	try {
		cout << "\naddPaving with copy of coll1 via SPS and NULL SPSnode" << endl;
		
		SPSnode* nullNode = NULL;
		copyColl1ViaSPS.addPaving(nullNode);
		
		assert(copyColl1ViaSPS.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has 1 element" << endl;
		
		string s = "addPavingCopyColl1ViaSPSAndNULL.txt";
		outputCollNode(s, copyColl1ViaSPS, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addPaving with copy of coll1 via SPS and NULL SPSnode:\n" << msg << endl;
	}
	try {
		cout << "\n+= with default coll node and sps node equivalent of coll2" << endl;
		
		copyColl2ViaSPS += spsHist2;
		
		assert(copyColl2ViaSPS.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has 1 element" << endl;

		string s = "additionDefaultAndColl2ViaSPS.txt";
		outputCollNode(s, copyColl2ViaSPS, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do += with default coll node and sps node equivalent of coll2:\n" << msg << endl;
	}
	try {
		cout << "\naddition with default coll node and sps node equivalent of coll3" << endl;
		
		CollatorSPnode temp;
		copyColl3ViaSPS = temp + spsHist3;
		
		assert(copyColl3ViaSPS.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has 1 element" << endl;

		string s = "additionDefaultAndColl3ViaSPS.txt";
		outputCollNode(s, copyColl3ViaSPS, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addition with default coll node and sps node equivalent of coll3:\n" << msg << endl;
	}
	
    
	try {
		cout << "\nanotherCopy1PlusEqualCopy2 += spsHist3" << endl;
		
		CollatorSPnode temp(anotherCopy1PlusEqualCopy2);
		temp += spsHist3;
		
		assert(temp.getSizeRangeCollection() == 3);
		cout << "Passed assert that collator has 3 elements" << endl;

		string s = "anotherCopyColl1PlusEqualColl2PlusEqualSPSHist3.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do anotherCopy1PlusEqualCopy2 += spsHist3:\n" << msg << endl;
	}
	CollatorSPnode sumOfThree;
	try {
		cout << "\naddition coll1 + coll2 + coll3" << endl;
		
		sumOfThree = collnode1 + collnode2 + collnode3;
		
		assert(sumOfThree.getSizeRangeCollection() == 3);
		cout << "Passed assert that collator has 3 elements" << endl;

		string s = "additionCollNode1CollNode2CollNode3.txt";
		outputCollNode(s, sumOfThree, 10);
		
		cout << "Check sumOfThree.getTotalAbsValueTimesVol()";
		cout <<" == collnode1.getTotalAbsValueTimesVol() ";
		cout << "+ collnode1.getTotalAbsValueTimesVol() ";
		cout << "+ collnode1.getTotalAbsValueTimesVol()" << endl;
		assert(sumOfThree.getTotalAbsValueTimesVol() 
		== (collnode1.getTotalAbsValueTimesVol() + 
			collnode2.getTotalAbsValueTimesVol() + collnode3.getTotalAbsValueTimesVol()));
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do addition coll1 + coll2 + coll3:\n" << msg << endl;
	}
	try {
		cout << "\ndefault += coll1 += coll2 += coll3" << endl;
		
		CollatorSPnode temp;
		CollatorSPnode temp1(collnode1);
		CollatorSPnode temp2(collnode2);
		CollatorSPnode temp3(collnode3);
		temp += temp1 +=temp2 += temp3;
		
		assert(temp.getSizeRangeCollection() == 3);
		cout << "Passed assert that collator has 3 elements" << endl;

		string s = "defaultPlusEqualCollNode1CollNode2CollNode3.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default += coll1 += coll2 += coll3:\n" << msg << endl;
	}
	try {
		cout << "\ncopy of sum of three + nothingCollatedNode" << endl;
		
		CollatorSPnode temp = sumOfThree + nothingCollatedNode;
		
		assert(temp.getSizeRangeCollection() == 3);
		cout << "Passed assert that collator has 3 elements" << endl;

		string s = "copySumOfThreePlusNothingCollatedNode.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of sum of three + collator with zero data:\n" << msg << endl;
	}
	CollatorSPnode noDataNode;
	CollatorSPnode sumOfThreePlusNoData;
	try {
		cout << "\ncopy of sum of three + collator with zero data" << endl;
		
		CollatorSPnode temp1(sumOfThree);
		CollatorSPnode temp2(pavingBox);
		// no split, but allocate a range
		std::vector< RealVec > zeroData = makeRangesZero(temp2.nodeRealVolume());
		temp2.allocateRanges(zeroData);
		noDataNode = temp2;
		sumOfThreePlusNoData = temp1 + noDataNode;
		
		assert(sumOfThreePlusNoData.getSizeRangeCollection() == 4);
		cout << "Passed assert that collator has 4 elements" << endl;

		string s1 = "zeroDataNode.txt";
		outputCollNode(s1, noDataNode, 10);
		assert( checkFileLines(s1, 1) );
		cout << "Passed assert that output file has one line" << endl;
		string s2 = "copySumOfThreePlusZeroData.txt";
		outputCollNode(s2, sumOfThreePlusNoData, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of sum of three + collator with zero data:\n" << msg << endl;
	}
	try {
		cout << "\nreduce coll with nothing collated" << endl;
		
		CollatorSPnode temp(pavingBox);
		temp.reduceCollation();
		assert(temp.isEmptyRangeCollection());
		cout << "Passed assert that collator is empty" << endl;

		string s = "reduceCollNothingCollated.txt";
		outputCollNode(s, temp, 10);
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce coll with nothing collated:\n" << msg << endl;
	}
	
	try {
		cout << "\nreduce copy of coll1" << endl;
		
		CollatorSPnode temp(collnode1);
		temp.reduceCollation();
		
		assert(temp.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;

		string s = "reduceCopyColl1.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of coll1:\n" << msg << endl;
	}
	try {
		cout << "\nreduce copy of sum of three" << endl;
		
		CollatorSPnode temp(sumOfThree);
		temp.reduceCollation();
		
		string s = "reduceCopySumOfThree.txt";
		outputCollNode(s, temp, 10);
		
		assert(temp.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;

		cout << "Check reduced copy of sumOfThree.getTotalAbsValueTimesVol()";
		cout <<" == sumOfThree.getTotalAbsValueTimesVol()" << endl;
		assert(temp.getTotalAbsValueTimesVol() 
		== sumOfThree.getTotalAbsValueTimesVol() );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of three:\n" << msg << endl;
	}
	try {
		cout << "\nreduce copy of sum of three + zeroData" << endl;
		
		CollatorSPnode temp(sumOfThreePlusNoData);
		temp.reduceCollation();
		
		assert(temp.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;
		
		string s = "reduceCopySumOfThreePlusNoData.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of three + zeroData:\n" << msg << endl;
	}
	
	CollatorSPnode averageOfThree;
	try {
		cout << "\nmakeAverage with sum of three" << endl;
		
		CollatorSPnode temp = sumOfThree.makeAverage();
		
		assert(temp.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;

		averageOfThree = temp;
		
		string s = "averageFromSumOfThree.txt";
		outputCollNode(s, averageOfThree, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeAverage from sum of three:\n" << msg << endl;
	}
	CollatorSPnode normalisationOfSumOfThree;
	try {
		cout << "\nmakeNormalised with sum of three" << endl;
		
		CollatorSPnode temp = sumOfThree.makeNormalised();
		
		assert(temp.getSizeRangeCollection() == 1);
		cout << "Passed assert that collator has one element" << endl;

		normalisationOfSumOfThree = temp;
		
		string s = "normalisationFromSumOfThree.txt";
		outputCollNode(s, normalisationOfSumOfThree, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeNormalised from sum of three:\n" << msg << endl;
	}
	try {
		cout << "\nreduce copy of sum of three and then average it" << endl;
		
		CollatorSPnode temp(sumOfThree);
		temp.reduceCollation();
		temp.average();
		
		string s = "averagedReduceCopySumOfThree.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of three and then average it:\n" << msg << endl;
	}
	try {
		cout << "\nreduce copy of sum of three and then normalise it" << endl;
		
		CollatorSPnode temp(sumOfThree);
		temp.reduceCollation();
		temp.normalise();
		
		string s = "normalisedReduceCopySumOfThree.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of three and then normalise it:\n" << msg << endl;
	}
	try {
		cout << "\nmake average from copy of sum of three plus zeroData" << endl;
		
		CollatorSPnode temp = sumOfThreePlusNoData.makeAverage();
		
		string s = "averagedCopySumOfThreePlusZeroData.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make average from copy of sum of three plus zeroData:\n" << msg << endl;
	}
	try {
		cout << "\nmake average (copy of sum of three + nothingCollatedNode)" << endl;
		
		CollatorSPnode temp = sumOfThree + nothingCollatedNode;
		temp.average();
		
		string s = "averagedCopySumOfThreePlusNothingCollatedNode.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make average (copy of sum of three + nothingCollatedNode):\n" << msg << endl;
	}
	
	try {
		cout << "\nnormalise copy of sum of three plus zeroData" << endl;
		
		CollatorSPnode temp = sumOfThreePlusNoData.makeNormalised();
		
		string s = "normalisedCopySumOfThreePlusZeroData.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make normalised from copy of sum of three plus zeroData:\n" << msg << endl;
	}
	
	// multiplication
	try {
		cout << "\nmultiply nothingCollatedNode by 10" << endl;
		cxsc::real mult(10.0);
		CollatorSPnode temp = nothingCollatedNode * mult;
		
		string s = "nothingCollatedNodeMultByTen.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to multiply nothingCollatedNode by 10:\n" << msg << endl;
	}
	
	try {
		cout << "\nmultiply coll1 by zero" << endl;
		cxsc::real mult(0.0);
		CollatorSPnode temp = collnode1 * mult;
		
		string s = "coll1MultByZero.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to multiply coll1 by zero:\n" << msg << endl;
	}
	
	try {
		cout << "\nmultiply coll1 by one" << endl;
		cxsc::real mult(1.0);
		CollatorSPnode temp = collnode1 * mult;
		
		string s = "coll1MultByOne.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to multiply coll1 by one:\n" << msg << endl;
	}
	CollatorSPnode collnode1Times10;
	try {
		cout << "\ncopy of coll1 multiplied by ten" << endl;
		
		CollatorSPnode temp = collnode1;
		cxsc::real mult(10.0);
		temp *= mult;
		collnode1Times10 = temp;
		
		string s = "coll1MultByTen.txt";
		outputCollNode(s, collnode1Times10, 10);
		
		cout << "Check collnode1Times10.getTotalAbsValueTimesVol()";
		cout <<" == 10 x collnode1.getTotalAbsValueTimesVol()" << endl;
		assert(temp.getTotalAbsValueTimesVol() 
		== (10 * collnode1.getTotalAbsValueTimesVol()) );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of coll1 multiplied by ten:\n" << msg << endl;
	}
	CollatorSPnode sumOfThreeTimes2;
	try {
		cout << "\ncopy of sum of Three multiplied by two" << endl;
		
		CollatorSPnode temp = sumOfThree;
		cxsc::real mult(2.0);
		temp *= mult;
		sumOfThreeTimes2 = temp;
		
		string s = "sumOfThreeMultByTwo.txt";
		outputCollNode(s, sumOfThreeTimes2, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of sum of Three multiplied by two:\n" << msg << endl;
	}
	CollatorSPnode sumOfThreePlusAddSumOfThree;
	try {
		cout << "\ncopy of sum of Three += sum of three" << endl;
		
		CollatorSPnode temp = sumOfThree;
		temp += sumOfThree;
		sumOfThreePlusAddSumOfThree = temp;
		
		string s = "sumOfThreeMultPlusAddSumOfThree.txt";
		outputCollNode(s, sumOfThreePlusAddSumOfThree, 10);
		
		cout << "Check sumOfThreePlusAddSumOfThree.getTotalAbsValueTimesVol()";
		cout <<" == sumOfThreeTimes2.getTotalAbsValueTimesVol()" << endl;
		assert(sumOfThreePlusAddSumOfThree.getTotalAbsValueTimesVol() 
		== sumOfThreeTimes2.getTotalAbsValueTimesVol() );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do sum of Three += sum of three:\n" << msg << endl;
	}
	
	try {
		cout << "\nreduce copy of sum of Three += sum of three" << endl;
		
		CollatorSPnode temp(sumOfThreePlusAddSumOfThree);
		temp.reduceCollation();
		
		string s = "reduceSumOfThreeMultPlusAddSumOfThree.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of Three += sum of three:\n" << msg << endl;
	}
	try {
		cout << "\nreduce copy of sum of Three multiplied by two" << endl;
		
		CollatorSPnode temp = sumOfThreeTimes2;
		temp.reduceCollation();
		
		string s = "reduceSumOfThreeMultByTwo.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do reduce copy of sum of Three multiplied by two:\n" << msg << endl;
	}
	

	// division
	try {
		cout << "\ndivide coll1 by zero (this should fail)" << endl;
		cxsc::real mult(0.0);
		CollatorSPnode temp = collnode1 / mult;
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to divide coll1 by zero:\n" << msg << endl;
	}
	try {
		cout << "\ndivide nothingCollatedNode by 10" << endl;
		cxsc::real div(10.0);
		CollatorSPnode temp = nothingCollatedNode / div;
		
		string s = "nothingCollatedNodeDivByTen.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to divide nothingCollatedNode by 10:\n" << msg << endl;
	}
	
	try {
		cout << "\ndivide coll1 by one" << endl;
		cxsc::real mult(1.0);
		CollatorSPnode temp = collnode1 / mult;
		
		string s = "coll1DivByOne.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to divide coll1 by one:\n" << msg << endl;
	}
	try {
		cout << "\ncopy of coll1 multiplied by ten divided by 10" << endl;
		
		CollatorSPnode temp = collnode1Times10;
		cxsc::real mult(10.0);
		temp /= mult;
		
		string s = "coll1MultByTenDividedByTen.txt";
		outputCollNode(s, temp, 10);
		
		cout << "Check collnode1Times10Div10.getTotalAbsValueTimesVol()";
		cout <<" == collnode1.getTotalAbsValueTimesVol()" << endl;
		assert(temp.getTotalAbsValueTimesVol() 
		== ( collnode1.getTotalAbsValueTimesVol()) );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of coll1 multiplied by ten divided by 10:\n" << msg << endl;
	}
	try {
		cout << "\ncopy of sum of Three multiplied by two divided by two" << endl;
		
		CollatorSPnode temp = sumOfThreeTimes2;
		cxsc::real mult(2.0);
		temp /= mult;
		
		string s = "sumOfThreeMultByTwoDividedByTwo.txt";
		outputCollNode(s, temp, 10);
		
		cout << "Check sumOfThreeMultByTwoDividedByTwo.getTotalAbsValueTimesVol()";
		cout <<" == sumOfThree.getTotalAbsValueTimesVol()" << endl;
		assert(temp.getTotalAbsValueTimesVol() 
		== ( sumOfThree.getTotalAbsValueTimesVol()) );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of sum of Three multiplied by two divided by two:\n" << msg << endl;
	}
	try {
		cout << "\ncopy of sum of Three divided by three and reduced" << endl;
		
		CollatorSPnode temp = sumOfThree;
		cxsc::real mult(3.0);
		temp /= mult;
		temp.reduceCollation();
		string s = "reducedSumOfThreeDividedByThree.txt";
		outputCollNode(s, temp, 10);
		
		cout << "Check reducedSumOfThreeDividedByThree.getTotalAbsValueTimesVol()";
		cout <<" == averageOfThree.getTotalAbsValueTimesVol()" << endl;
		assert(temp.getTotalAbsValueTimesVol() 
		== ( averageOfThree.getTotalAbsValueTimesVol()) );
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy of sum of Three divided by three and reduced:\n" << msg << endl;
	}
	
	// check normalise with empty nodes and children
	try {
		cout << "\ncall normalise on an empty node (this should fail)" << endl;
		
		CollatorSPnode temp;
		temp.normalise();
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call normalise on an empty node:\n" << msg << endl;
	}
	try {
		cout << "\ncall makeNormalised on an empty node (this should fail)" << endl;
		
		CollatorSPnode temp;
		CollatorSPnode norm = temp.makeNormalised();
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeNormalised on an empty node:\n" << msg << endl;
	}
	
	try {
		cout << "\nnormalise with node with nothing collated  (this should fail)" << endl;
		
		CollatorSPnode temp(pavingBox);
		temp.normalise();
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do normalise with node with nothing collated:\n" << msg << endl;
	}
	
	try {
		cout << "\nmakeNormalised with zero data node (this should fail since no normaliser)" << endl;
		
		CollatorSPnode temp = noDataNode.makeNormalised();
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeNormalised with zero data node:\n" << msg << endl;
	}

	try {
		cout << "\ncall normalise on a child node (this should fail)" << endl;
		CollatorSPnode temp(collnode1);
		temp.getRightChild()->normalise();
		
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call normalise on a child node:\n" << msg << endl;
	}
	
	
	try {
		cout << "\ncall makeNormalised on a child node" << endl;
		CollatorSPnode temp(collnode1);
		CollatorSPnode norm = temp.getRightChild()->makeNormalised();
		string s = "coll1RightChildMakeNormalised.txt";
		outputCollNode(s, norm, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeNormalised on a child node:\n" << msg << endl;
	}
	try {
		cout << "\ncopy of right child of normalised collator 1" << endl;
		CollatorSPnode temp(collnode1);
		temp.normalise();
		CollatorSPnode rc( *(temp.getRightChild() ) );
		string s = "coll1NormalisedRightChild.txt";
		outputCollNode(s, rc, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to get child node of normalised collator:\n" << msg << endl;
	}
	// check average with empty nodes and children
	try {
		cout << "\ncall average on an empty node  (this should fail)" << endl;
		
		CollatorSPnode temp;
		temp.average();
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call average on an empty node:\n" << msg << endl;
	}
	try {
		cout << "\ncall makeAverage on an empty node  (this should fail)" << endl;
		
		CollatorSPnode temp;
		CollatorSPnode av = temp.makeAverage();
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeAverage on an empty node:\n" << msg << endl;
	}
	try {
		cout << "\naverage with node with nothing collated  (this should fail)" << endl;
		
		CollatorSPnode temp(pavingBox);
		temp.average();
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do average with node with nothing collated:\n" << msg << endl;
	}
	try {
		cout << "\nmakeAverage with zero data node (this should work)" << endl;
		
		CollatorSPnode temp = noDataNode.makeAverage();
		
		string s = "averageFromZeroDataNode.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeAverage with zero data node:\n" << msg << endl;
	}
	try {
		cout << "\ncall average on a child node (this should fail)" << endl;
		CollatorSPnode temp(collnode1);
		temp.getRightChild()->average();

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call average on a child node:\n" << msg << endl;
	}
	try {
		cout << "\ncall makeAverage on a child node" << endl;
		CollatorSPnode temp(collnode1);
		CollatorSPnode av = temp.getRightChild()->makeAverage();
		string s = "coll1RightChildMakeAverage.txt";
		outputCollNode(s, av, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeAverage on a child node:\n" << msg << endl;
	}
	
	try {
		cout << "\ncopy right child of collator 1, divide by number collated, and reduce" << endl;
		CollatorSPnode temp(*collnode1.getRightChild());
		cxsc::real div(1.0*collnode1.getSizeRangeCollection());
		temp /= div;
		temp.reduceCollation();
		string s = "coll1RightChildHandRolledAverage.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to make a hand rolled average on a child node:\n" << msg << endl;
	}
	
	// check marginalise with empty nodes and children
	
	std::vector<int> reqDims;
	reqDims.push_back(1);
	try {
		cout << "\ncall marginalise on an empty node (this should fail)" << endl;
		
		CollatorSPnode temp;
		temp.marginalise(reqDims);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call marginalise on an empty node:\n" << msg << endl;
	}
	try {
		cout << "\ncall makeMarginal on an empty node (this should fail)" << endl;
		
		CollatorSPnode temp;
		CollatorSPnode m = temp.makeMarginalised(reqDims);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeMarginal on an empty node:\n" << msg << endl;
	}
	try {
		cout << "\ncall marginalise on nothingCollatedNode(this should work)" << endl;
		
		CollatorSPnode temp(nothingCollatedNode);
		temp.marginalise(reqDims);
		
		string s = "nothingCollatedNodeMarginalised.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call marginalise on nothingCollatedNode:\n" << msg << endl;
	}
	
	try {
		cout << "\ncall marginalise on a child node (this should fail)" << endl;
		CollatorSPnode temp(collnode1);
		temp.getRightChild()->marginalise(reqDims);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call marginalise on a child node:\n" << msg << endl;
	}
	try {
		cout << "\ncall makeMarginal on a child node" << endl;
		CollatorSPnode temp(collnode1);
		CollatorSPnode m = temp.getRightChild()->makeMarginalised(reqDims);
		string s = "coll1RightChildMakeMarginal.txt";
		outputCollNode(s, m, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to call makeMarginal on a child node:\n" << msg << endl;
	}
	
	try {
		cout << "\nmakeMarginalised with zero data node (result will have one value of 0 in range collation)" << endl;
		
		std::vector < int > reqDims(1,1);
		CollatorSPnode temp = noDataNode.makeMarginalised(reqDims);
		
		string s = "marginalisedFromZeroDataNode.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do makeMarginalised with zero data node:\n" << msg << endl;
	}
	
	//marginalise dimensions
	try {
		cout << "\nmarginalise copy of coll1 with no dims required (this should fail)" << endl;
		
		CollatorSPnode temp = collnode1;
		std::vector<int> r;
		//reqDims.push_back(1);
	
		temp.marginalise(r);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to marginalise copy of coll1 with no dims required:\n" << msg << endl;
	}
	try {
		cout << "\nmarginalise copy of coll1 with all dims required" << endl;
		
		CollatorSPnode temp = collnode1;
		assert(temp.getDimension() == 2);
		std::vector<int> r;
		r.push_back(1);
		r.push_back(2);
		assert(temp.getDimension() == 2);
	
		temp.marginalise(r);
		string s = "collnode1MarginalisedAllReq.txt";
		outputCollNode(s, temp, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to marginalise copy of coll1 with all dims required:\n" << msg << endl;
	}
	try {
		cout << "\nmarginalise copy of coll1 with req dims not present (this should fail)" << endl;
		
		CollatorSPnode temp = collnode1;
		std::vector<int> r;
		r.push_back(1);
		r.push_back(3);
	
		temp.marginalise(r);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to marginalise copy of coll1 with req dims not present:\n" << msg << endl;
	}
	
	
	try {
		cout << "\nmarginalise copy of collnode1 on dim 1" << endl;
		
		CollatorSPnode temp = collnode1;
		assert(temp.getDimension() == 2);
		std::vector<int> r;
		r.push_back(1);
		
		temp.marginalise(r);
		assert(temp.getDimension() == 1);
		cout << "Passed assert that collator has one dim" << endl;
		
		string s = "collnode1MarginalisedOn1.txt";
		outputCollNode(s, temp, 10);
		
		cout << "Check collnode1.getTotalAbsValueTimesVol() == temp.getTotalAbsValueTimesVol()";
		cout << "collnode1Marginalised.getTotalAbsValueTimesVol()" << endl;
		assert(collnode1.getTotalAbsValueTimesVol() ==
			temp.getTotalAbsValueTimesVol());
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to marginalise copy of coll1 on dim 1:\n" << msg << endl;
	}
	
	try {
		cout << "\nmarginalise copy of sumOfThree on dim 2" << endl;
		
		CollatorSPnode temp = sumOfThree;
		assert(temp.getDimension() == 2);
		//outputCollAllNodes("sumOfThreeAllOutput.txt", temp, 5);
		std::vector<int> r;
		r.push_back(2);
		
		temp.marginalise(r);
		assert(temp.getDimension() == 1);
		assert(temp.getSizeRangeCollection() == 3);
		cout << "Passed asserts that collator has one dim and 3 elements" << endl;
		string s = "sumOfThreeMarginalisedOn2.txt";
		outputCollNode(s, temp, 10);
		
		cout << "Check sumOfThree.getTotalAbsValueTimesVol() == ";
		cout << "sumOfThreeMarginalised.getTotalAbsValueTimesVol()" << endl;
		assert(sumOfThree.getTotalAbsValueTimesVol() ==
			temp.getTotalAbsValueTimesVol());
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to marginalise copy of sumOfThree on dim 2:\n" << msg << endl;
	} 
	
	try {
		cout << "\nsum of marginals of collnode1, collnode2, collnode3" << endl;
		
		std::vector<int> r;
		r.push_back(2);
		CollatorSPnode temp1 = collnode1.makeMarginalised(r);
		CollatorSPnode temp2 = collnode2.makeMarginalised(r);
		CollatorSPnode temp3 = collnode3.makeMarginalised(r);
		
		CollatorSPnode temp = temp1 + temp2 + temp3;
		assert(temp.getDimension() == 1);
		assert(temp.getSizeRangeCollection() == 3);
		cout << "Passed asserts that collator has one dim and 3 elements" << endl;
		
		string s = "sumOfThreeMarginalsOn2.txt";
		outputCollNode(s, temp, 10);
		cout << "Check sumOfThree.getTotalAbsValueTimesVol() == ";
		cout << "sumOfThreeMarginalsOn2Reduced.getTotalAbsValueTimesVol()" << endl;
		assert(sumOfThree.getTotalAbsValueTimesVol() ==
			temp.getTotalAbsValueTimesVol());
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do sum of marginals of collnode1, collnode2, collnode3:\n" << msg << endl;
	}
	//find containing node
	try {
		
		cout << "\nFind containing node with empty paving (should fail)" << endl;
		
		CollatorSPnode temp;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const CollatorSPnode* node = temp.findContainingNode(pt);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with empty paving:\n" << msg << endl;
	}
	try {
		
		cout << "\nFind containing node with zero data node" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const CollatorSPnode* node = noDataNode.findContainingNode(pt);
		
		std::string expected("X");
		assert(node->getNodeName() == expected);
		cout << "Passed assert that containing node is " << expected << endl;
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with zero data node:\n" << msg << endl;
	}
	try {
		
		cout << "\nFind containing node with nothingCollatedNode" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const CollatorSPnode* node = noDataNode.findContainingNode(pt);
		
		std::string expected("X");
		assert(node->getNodeName() == expected);
		cout << "Passed assert that containing node is " << expected << endl;
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with nothingCollatedNode:\n" << msg << endl;
	}
	
	
	try {
		
		cout << "\nFind containing node with coll1 (should not find containing node)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -3;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const CollatorSPnode* node = collnode1.findContainingNode(pt);
		
		assert(node == NULL);
		cout << "Passed assert that no containing node found" << endl;
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with coll1:\n" << msg << endl;
	}
	
	try {
		
		cout << "\nFind containing node with coll1 (should find all)" << endl;
		
		{
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = -2;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XLLL");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = 2;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XLR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = 2;
			pt[2] = -2;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XRL");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = 2;
			pt[2] = 2;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XRRR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with coll1:\n" << msg << endl;
	}
	
	try {
		
		cout << "\nMore find containing node with coll1 (should find)" << endl;
		
		{
			cxsc::rvector pt(d);
			pt[1] = 0;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XRRLL");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = 1.0;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XRRR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XLR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		{
			cxsc::rvector pt(d);
			pt[1] = -1;
			pt[2] = -1;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const CollatorSPnode* node = collnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XLLRR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with coll1:\n" << msg << endl;
	}

	// L1 distance to average of other - test tricky cases 
	try {
			
		cout << "\nL1 distance collnode1 against null pointer (should fail)  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		CollatorSPnode* temp2 = NULL;
		
		RealVec distances = temp1.getL1DistancesToAverage(temp2);
		
		throw std::logic_error("Should not be able to do this");

	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do collnode1 against null pointer:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance empty paving against collnode1 (should fail)  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		const CollatorSPnode temp2;
		
		RealVec distances = temp2.getL1DistancesToAverage(&temp1);
		
		throw std::logic_error("Should not be able to do this");

	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance empty paving against collnode1:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collnode1 against empty paving (should fail)  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		const CollatorSPnode temp2;
		
		RealVec distances = temp1.getL1DistancesToAverage(&temp2);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collnode1 against empty paving:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collnode1 against paving with box of incorrect dimensions (should fail)  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		int d2 = 3; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval2(-2,2);
		for(int k=1; k <= d2; k++) pavingBox2[k] = pavingInterval2;

		CollatorSPnode temp2(pavingBox2);
				
		RealVec distances = temp1.getL1DistancesToAverage(&temp2);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against paving with box of incorrect dimensions:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collnode1 against paving with box of incorrect side lengths (should fail)  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		int d2 = 2; // dimension of the box to sample data from
		ivector pavingBox2(d2);
		interval pavingInterval2(-3,2);
		pavingBox2[1] = pavingInterval;
		pavingBox2[2] = pavingInterval2;

		CollatorSPnode temp2(pavingBox2);
				
		RealVec distances = temp1.getL1DistancesToAverage(&temp2);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance against paving with box of incorrect side lengths:\n" << msg << endl;
	}
	
	
	try {
			
		cout << "\nL1 distance nothingCollatedNode against collnode1" << endl;
		
		const CollatorSPnode temp1(collnode1);
				
		RealVec distances (10);
		
		distances = nothingCollatedNode.getL1DistancesToAverage(distances, &temp1);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		assert(distances.empty());
		
		cout << "Passed assert that distances is empty" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance nothingCollatedNode against collnode1:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collnode1 against nothingCollatedNode (this should fail)" << endl;
		
		const CollatorSPnode temp1(collnode1);
				
		RealVec distances (10);
		
		distances = temp1.getL1DistancesToAverage(distances, &nothingCollatedNode);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance collnode1 against nothingCollatedNode:\n" << msg << endl;
	}
	
	
	try {
			
		cout << "\nL1 distance zero data node against collnode1" << endl;
		
		const CollatorSPnode temp1(collnode1);
				
		RealVec distances (10);
		
		distances = noDataNode.getL1DistancesToAverage(distances, &temp1);
		
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
		cout << "\nFailed to do L1 distance zero data node against collnode1:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance collnode1 against zeroDataNode  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances, &noDataNode);
		
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
		cout << "\nFailed to do L1 distance collnode1 against zero value collator:\n" << msg << endl;
	}
	
	// L1 distance to average of itself - test tricky cases 
	try {
			
		cout << "\nL1 distance empty paving against itself (should fail)  " << endl;
		
		const CollatorSPnode temp;
		
		RealVec distances = temp.getL1DistancesToAverage();

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance empty paving  against itself:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance nothingCollatedNode against itself (should fail)" << endl;
		
		RealVec distances (10);
		
		distances = nothingCollatedNode.getL1DistancesToAverage(distances);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance nothingCollatedNode against itself:\n" << msg << endl;
	}
		
	try {
			
		cout << "\nL1 distance zero data node against itself (should work because distance of anything to itself is 0)" << endl;
		
		RealVec distances (10);
		
		distances = noDataNode.getL1DistancesToAverage(distances);
		
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
		cout << "\nFailed to do L1 distance zero data node against itself:\n" << msg << endl;
	}
	
	// L1 distance to another node with real nodes
	
	try {
			
		cout << "\nL1 distance collnode1 against collnode1  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances, &temp1);
		
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
		cout << "\nFailed to do L1 distance collnode1 against collnode1:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance collnode1 against itself  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances);
		
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
		cout << "\nFailed to do L1 distance collnode1 against itself:\n" << msg << endl;
	}
	
	cxsc::real dis2_3(0.0);
	try {
			
		cout << "\nL1 distance collnode2 against collnode3  " << endl;
		
		const CollatorSPnode temp1(collnode2);
		const CollatorSPnode temp2(collnode3);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances, &temp2);
		
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
		cout << "\nFailed to do L1 distance collnode2 against collnode3:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance collation of 2 collnode2's against collnode 3  " << endl;
		
		CollatorSPnode temp1(collnode2);
		temp1.addPaving(&collnode2);
		const CollatorSPnode temp2(collnode3);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances, &temp2);
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
		cout << "\nFailed to do L1 distance collation of 2 collnode2's against collnode 3:\n" << msg << endl;
	}
	cxsc::real dis1_2(0.0);
	
	try {
			
		cout << "\nL1 distance collnode1 against collnode2  " << endl;
		
		const CollatorSPnode temp1(collnode1);
		const CollatorSPnode temp2(collnode2);
		
		RealVec distances;
		
		distances = temp1.getL1DistancesToAverage(distances, &temp2);
		
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
		cout << "\nFailed to do L1 distance collnode1 against collnode2:\n" << msg << endl;
	}
	RealVec distancesAgainstAverageOfThree;
	try {
			
		cout << "\nL1 distance for sumOfThree against its own average  " << endl;
		
		const CollatorSPnode temp1(sumOfThree);
		distancesAgainstAverageOfThree = temp1.getL1DistancesToAverage(distancesAgainstAverageOfThree);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances against average are " << distancesAgainstAverageOfThree << endl;
	
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for sumOfThree against its own average:\n" << msg << endl;
	}
	try {
			
		cout << "\nL1 distance for sumOfThree against averageOfThree  " << endl;
		
		const CollatorSPnode temp1(sumOfThree);
		
		RealVec distances;
		distances = temp1.getL1DistancesToAverage(distances, &averageOfThree);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances against averageOfThree are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		RealVecItr ait = distancesAgainstAverageOfThree.begin();
		for (RealVecItr it = distances.begin(); it < distances.end(); ++it) {
			
			bool isNotEqual = cxsc::abs( (*it) - (*ait) ) 
				> DBL_EPSILON * cxsc::max((*it), (*ait));
			assert(!isNotEqual);
			++ait;
		}
		cout << "Passed asserts that distances are same as for sumOfThree against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for sumOfThree against averageOfThree:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance for collnode1 against sumOfThree  " << endl;
		
		RealVec distances = collnode1.getL1DistancesToAverage(&sumOfThree);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs( distances.at(0) - distancesAgainstAverageOfThree.at(0) ) 
				> DBL_EPSILON * cxsc::max(distances.at(0), distancesAgainstAverageOfThree.at(0));
		assert(!isNotEqual);
		
		cout << "Passed assert that distance is same as using sumOfThree against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for collnode1 against sumOfThree:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance for averageOfThree against collnode2 " << endl;
		
		RealVec distances = averageOfThree.getL1DistancesToAverage(&collnode2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances are " << distances << endl;
		
		cout << cxsc::RestoreOpt;
		
		bool isNotEqual = cxsc::abs( distances.at(0) - distancesAgainstAverageOfThree.at(1) ) 
			> DBL_EPSILON * cxsc::max(distances.at(0), distancesAgainstAverageOfThree.at(1));
		assert(!isNotEqual);
		
		cout << "Passed assert that distance is same as using sumOfThree against its own average" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance for averageOfThree against collnode2:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance for averageOfThree against sumOfThree " << endl;
		
		RealVec distances = averageOfThree.getL1DistancesToAverage(&sumOfThree);
		
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
		cout << "\nFailed to do L1 distance for averageOfThree against sumOfThree:\n" << msg << endl;
	}
	
// test L1 distance from collator to spsnode

	try {
			
		cout << "\nL1 distance from collator to null sps node (should fail) " << endl;
		
		SPSnode* spn = NULL;
		
		RealVec distances = collnode1.getL1Distances(spn);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do distance from collator to null sps node:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance from collator to empty sps node (should fail) " << endl;
		
		SPSnode spn;
		
		RealVec distances = collnode1.getL1Distances(&spn);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do L1 distance from collator to empty sps node:\n" << msg << endl;
	}
	
	try {
			
		cout << "\nL1 distance for sumOfThree against spsnode with no data  " << endl;
		
		const SPSnode temp(pavingBox);
		
		RealVec distances;
		distances = sumOfThree.getL1Distances(distances, &temp);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for sumOfThree against spsnode with no data are " << distances << endl;
		
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
		cout << "\nFailed to do L1 distance for sumOfThree against spsnode with no data:\n" << msg << endl;
	}	

	try {
			
		cout << "\nL1 distance for collnode1 against spsHist1  " << endl;
		
		RealVec distances = collnode1.getL1Distances(&spsHist1);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for collnode1 against spsHist1 are " << distances << endl;
		
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
		cout << "\nFailed to do L1 distance for collnode1 against spsHist1 :\n" << msg << endl;
	}	
	try {
			
		cout << "\nL1 distance for collnode3 against spsHist2  " << endl;
		
		RealVec distances = collnode3.getL1Distances(&spsHist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for collnode1 against spsHist1 are " << distances << endl;
		
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
		cout << "\nFailed to do L1 distance for collnode3 against spsHist2 :\n" << msg << endl;
	}	
	try {
			
		cout << "\nL1 distance for averageOfThree against spsHist3  " << endl;
		
		RealVec distances = averageOfThree.getL1Distances(&spsHist3);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for averageOfThree against spsHist3 are " << distances << endl;
		
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
	}	
	try {
			
		cout << "\nL1 distance for sumOfThree against spsHist2  " << endl;
		
		RealVec distances = sumOfThree.getL1Distances(&spsHist2);
		
		cout << cxsc::SaveOpt << cxsc::SetPrecision(23,15);
		cout << "Distances for sumOfThree against spsHist2 are " << distances << endl;
		
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
		cout << "\nFailed to do L1 distance for averageOfThree against spsHist3:\n" << msg << endl;
	}	


	cout << "\nEnd of test\n" << endl;

    return 0;

} // end of test program


void outputCollNode(const std::string& s, const CollatorSPnode& spn, const int prec)
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		spn.leavesOutputTabs(os, prec);
		os.close();
		cout << s << " output to file" << endl;
	}
	else cout << "Error opening file " << s << endl;
}

void outputCollAllNodes(const std::string& s, const CollatorSPnode& spn, const int prec)
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		spn.nodesAllOutput(os, 1, prec);
		os.close();
		cout << s << " output to file" << endl;
	}
	else cout << "Error opening file " << s << endl;
}

void swapCheckOutput(std::string& s, const CollatorSPnode& spn, int level)
{
	ofstream os;
	if (spn.getParent() == NULL) os.open(s.c_str()); // don't append
	else os.open(s.c_str(), ios_base::app); //append
	if (os.is_open()) {

		// do me
		for (int i = 0; i < level; ++i) {
			os << "\t";
		}
			 
		os << spn.nodeStringSummary() << endl;
		// recurse
		
		os.close();
		if ( spn.getLeftChild() ) swapCheckOutput(s, *spn.getLeftChild(), level+1);
		if ( spn.getRightChild() ) swapCheckOutput(s, *spn.getRightChild(), level+1);
		
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


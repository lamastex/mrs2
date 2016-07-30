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
\brief Testing RealMappedSPnode
 */

#include "realmappedspnode.hpp"
#include "subpaving_exception.hpp"

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <iterator>  
#include <cfloat> // DBL_EPSILON
#include <cassert> // assert

using namespace cxsc;
using namespace std;
using namespace subpavings;

void outputNode(const std::string& s, const RealMappedSPnode& spn, const int prec = 5);
void outputAllNodes(const std::string& s, const RealMappedSPnode& spn, const int prec = 5);
void swapCheckOutput(std::string& s, const RealMappedSPnode& spn, int level = 0);
std::vector< real > makeRangesZero(cxsc::real rootvol);
std::vector < real > makeRanges1(cxsc::real rootvol);		
std::vector < real > makeRanges2(cxsc::real rootvol);		
std::vector < real > makeRanges3(cxsc::real rootvol);		
bool checkFileLines(const std::string& s, std::size_t expectedLines);
		
int main()
{
	
	try {
			cout << "\nDefault constructor" << endl;
			
			RealMappedSPnode temp;
			
			
			string s = "defaultConstructor.txt";
			outputNode(s, temp);
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
			
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do default constructed RealMappedSPnode:\n" << msg << endl;
	}
	
	
	try {
		cout << "\nNode constructed with default box (this should fail)" << endl;
		
		cxsc::ivector uselessBox;
		RealMappedSPnode temp(uselessBox);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& ee) {
		cout << "\nFailed to construct NewCollatorSpnode with default box:\n" << ee.what() << endl;
	}


	try {
			cout << "\nCopy constructor on default" << endl;
			
			RealMappedSPnode temp1;
			RealMappedSPnode temp2(temp1);
			
			string s = "copyOfdefaultConstructor.txt";
			outputNode(s, temp2);
			assert( checkFileLines(s, 0) );
			cout << "Passed assert that output file is empty" << endl;
				
		}
		catch (std::exception& ee) {
			std::string msg(ee.what());
			cout << "\nFailed to do copy construction of default constructed RealMappedSPnode:\n" << msg << endl;
		}
	try {
		cout << "\nAssignment copy of default constructor" << endl;
		
		RealMappedSPnode temp1;
		RealMappedSPnode temp2 = temp1;
		
		string s = "assignmentCopyOfdefaultConstructor.txt";
		outputNode(s, temp2);
		assert( checkFileLines(s, 0) );
		cout << "Passed assert that output file is empty" << endl;
	
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do assignment copy of default constructed RealMappedSPnode:\n" << msg << endl;
	}
	try {
		cout << "\nSplit to shape on default constructed node (this should fail)" << endl;
		
		RealMappedSPnode temp;
		temp.splitRootToShape("1,1");
		
		throw std::logic_error("Should not be able to do this");

	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do splitToShape of default constructed RealMappedSPnode:\n" << msg << endl;
	}
	
	int d = 2; // dimension of the box to sample data from
    ivector pavingBox(d);
    interval pavingInterval(-2,2);
    for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;

	try {
		cout << "\nconstructor with box" << endl;
		
		RealMappedSPnode temp(pavingBox);
		string s = "constructorWithBox.txt";
		outputNode(s, temp, 10);
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do constructor with box:\n" << msg << endl;
	}
	
	try {
		cout << "\ncopy constructor of constructor with box" << endl;
		
		RealMappedSPnode temp1(pavingBox);
		RealMappedSPnode temp2(temp1);
		string s = "copyConstructorOfConstructorWithBox.txt";
		outputNode(s, temp2, 10);
		assert( checkFileLines(s, 1) );
		cout << "Passed assert that output file has one line" << endl;
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy constructor of constructor with box:\n" << msg << endl;
	}

	

    RealMappedSPnode rmspnode1(pavingBox);
    RealMappedSPnode rmspnode2(pavingBox);
    RealMappedSPnode rmspnode3(pavingBox);
    
	std::string split1 = "3,4,4,2,2,4,4,3";
	std::string split2 = "2,3,4,4,2,3,3";
	std::string split3 = "1,2,3,3";
		
	try {
		cout << "\nSplit to shape"  << endl;
		rmspnode1.splitRootToShape(split1);
		rmspnode2.splitRootToShape(split2);
		rmspnode3.splitRootToShape(split3);

		string s1 = "splitToShapeNode1.txt";
		outputNode(s1, rmspnode1);
		string s2 = "splitToShapeNode2.txt";
		outputNode(s2, rmspnode2);
		string s3 = "splitToShapeNode3.txt";
		outputNode(s3, rmspnode3);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do splitToShape:\n" << msg << endl;
	}
	
	try {
		cout << "\nAllocate ranges"  << endl;
		std::vector< real > ranges = makeRanges1(rmspnode1.nodeRealVolume());
		rmspnode1.allocateRanges(ranges);
		string s = "RMSPNode1.txt";
		outputNode(s, rmspnode1, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for rmspnode1:\n" << msg << endl;
	}
	try {
		std::vector< real > ranges = makeRanges2(rmspnode2.nodeRealVolume());
		rmspnode2.allocateRanges(ranges);
		string s = "RMSPNode2.txt";
		outputNode(s, rmspnode2, 10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for rmspnode2:\n" << msg << endl;
	}
	try {
		std::vector< real > ranges = makeRanges3(rmspnode3.nodeRealVolume());
		rmspnode3.allocateRanges(ranges);
		string s = "RMSPNode3.txt";
		outputNode(s, rmspnode3,10);
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to allocate ranges for rmspnode3:\n" << msg << endl;
	}
	try {
		cout << "\nCheck total value x vol"  << endl;
		cxsc::real valvol = rmspnode1.getTotalLeafAreaRangeWithBox();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for rmspnode1 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for rmspnode1:\n" << msg << endl;
	}
	try {
		cxsc::real valvol = rmspnode2.getTotalLeafAreaRangeWithBox();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for rmspnode2 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for rmspnode2:\n" << msg << endl;
	}
	try {
		cxsc::real valvol = rmspnode3.getTotalLeafAreaRangeWithBox();
		cout << cxsc::SaveOpt;
		cout << cxsc::Scientific << cxsc::SetPrecision(23,15);
		cout << "total value x vol for rmspnode3 is " << valvol << endl;
		cout << cxsc::RestoreOpt;
		assert( valvol == 1.0);
		cout << "Passed assert that Val x vol == 1.0" << endl;
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to calculate val x vol for rmspnode3:\n" << msg << endl;
	}
	try {
		cout << "\nCopy constructor of rmspnode1" << endl;
		
		RealMappedSPnode cpy1(rmspnode1);
		
		string s = "copyConstructorNode1.txt";
		outputNode(s, cpy1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy construction of rmspnode1:\n" << msg << endl;
	}
	RealMappedSPnode copyColl1;
	RealMappedSPnode copyColl2;
	RealMappedSPnode copyColl3;
	
	try {
		cout << "\nswap of copy of rmspnode1 and copy of rmspnode2" << endl;
		
		RealMappedSPnode temp1(rmspnode1);
		RealMappedSPnode temp2(rmspnode2);
		
		string s11 = "copyOfNode1BeforeSwap.txt";
		swapCheckOutput(s11, temp1);
		string s12 = "copyOfNode2BeforeSwap.txt";
		swapCheckOutput(s12, temp2);
		
		std::swap(temp1, temp2);
		
		string s21 = "copyOfShouldBeLikeNode2AfterSwap.txt";
		swapCheckOutput(s21, temp1);
		string s22 = "copyOfShouldBeLikeNode1AfterSwap.txt";
		swapCheckOutput(s22, temp2);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do std::swap of rmspnode1:\n" << msg << endl;
	}
	
	
	try {
		cout << "\nCopy assignment of rmspnode1" << endl;
		
		copyColl1 = rmspnode1;
		
		string s = "copyAssignmentNode1.txt";
		outputNode(s, copyColl1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do copy assignment of rmspnode1:\n" << msg << endl;
	}
	
	//find containing node
	try {
		
		cout << "\nFind containing node with empty paving (should fail)" << endl;
		
		RealMappedSPnode temp;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const RealMappedSPnode* node = temp.findContainingNode(pt);
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with empty paving:\n" << msg << endl;
	}
	
	try {
		
		cout << "\nFind containing node just root box" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = 0;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		RealMappedSPnode tmp(pavingBox);
		
		const RealMappedSPnode* node = tmp.findContainingNode(pt);
		
		std::string expected("X");
		assert(node->getNodeName() == expected);
		cout << "Passed assert that containing node is " << expected << endl;
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with just root box node:\n" << msg << endl;
	}
	
	
	try {
		
		cout << "\nFind containing node with node1 (should not find containing node)" << endl;
		
		cxsc::rvector pt(d);
		pt[1] = -3;
		pt[2] = 0;
		
		cout << "Point ";
		prettyPrint(cout, pt);
		cout << endl;
		
		const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
		
		assert(node == NULL);
		cout << "Passed assert that no containing node found" << endl;
		
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with node1:\n" << msg << endl;
	}
	
	try {
		
		cout << "\nFind containing node with node1 (should find all)" << endl;
		
		{
			cxsc::rvector pt(d);
			pt[1] = -2;
			pt[2] = -2;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XRRR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with node1:\n" << msg << endl;
	}
	
	try {
		
		cout << "\nMore find containing node with node1 (should find)" << endl;
		
		{
			cxsc::rvector pt(d);
			pt[1] = 0;
			pt[2] = 0;
			
			cout << "Point ";
			prettyPrint(cout, pt);
			cout << endl;
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
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
			
			const RealMappedSPnode* node = rmspnode1.findContainingNode(pt);
			assert(node != NULL);
			//cout << "Containing node is " << node->getNodeName() << endl;
			std::string expected("XLLRR");
			assert(node->getNodeName() == expected);
			cout << "Passed assert that containing node is " << expected << endl;
		}
		
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do find containing node with node1:\n" << msg << endl;
	}
	
	#define SLICE
	#ifdef SLICE
	try {
		cout << "\nslice on copy of rmspnode3 on dim 1" << endl;
		
		RealMappedSPnode temp1(rmspnode3);
		
		vector < int > sliceDims;
		int d = 1;
		sliceDims.push_back(d);
		
		vector < real > slicePts;
		double dr = 0.5;
		slicePts.push_back(real(dr));
		
		temp1.slice(sliceDims, slicePts);
		ostringstream oss;
		oss << "node3sliced_d" << d << "_" << dr << ".txt";
		string s = oss.str();
		outputNode(s, temp1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do slice rmspnode3:\n" << msg << endl;
	}
	
	try {
		cout << "\nslice on copy of rmspnode3 on dim 2" << endl;
		
		RealMappedSPnode temp1(rmspnode3);
		
		vector < int > sliceDims;
		int d = 2;
		sliceDims.push_back(d);
		
		vector < real > slicePts;
		double dr = 0.0;
		slicePts.push_back(real(dr));
		
		temp1.slice(sliceDims, slicePts);
		
		ostringstream oss;
		oss << "node3sliced_d" << d << "_" << dr << ".txt";
		string s = oss.str();
		outputNode(s, temp1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do slice rmspnode3:\n" << msg << endl;
	}
	
	try {
		cout << "\nslice on copy of rmspnode1 on dim 1" << endl;
		
		RealMappedSPnode temp1(rmspnode1);
		
		vector < int > sliceDims;
		int d = 1;
		sliceDims.push_back(d);
		
		vector < real > slicePts;
		double dr = 1.0;
		slicePts.push_back(real(dr));
		
		temp1.slice(sliceDims, slicePts);
		ostringstream oss;
		oss << "node1sliced_d" << d << "_" << dr << ".txt";
		string s = oss.str();
		outputNode(s, temp1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do slice rmspnode1:\n" << msg << endl;
	}
	
	try {
		cout << "\nslice on copy of rmspnode1 on dim 2" << endl;
		
		RealMappedSPnode temp1(rmspnode1);
		
		vector < int > sliceDims;
		int d = 2;
		sliceDims.push_back(d);
		
		vector < real > slicePts;
		double dr = 1.0;
		slicePts.push_back(real(dr));
		
		temp1.slice(sliceDims, slicePts);
		
		ostringstream oss;
		oss << "node1sliced_d" << d << "_" << dr << ".txt";
		string s = oss.str();
		outputNode(s, temp1, 10);
		
	}
	catch (std::exception& ee) {
		std::string msg(ee.what());
		cout << "\nFailed to do slice rmspnode3:\n" << msg << endl;
	}
	#endif
	
	cout << "\n\nend of test program" << endl;
    return 0;

} // end of test program


void outputNode(const std::string& s, const RealMappedSPnode& spn, const int prec)
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

void outputAllNodes(const std::string& s, const RealMappedSPnode& spn, const int prec)
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

void swapCheckOutput(std::string& s, const RealMappedSPnode& spn, int level)
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

std::vector< real > makeRangesZero(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real(0.0));//X
	
	return result;
	
}

std::vector< real > makeRanges1(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((8.0/8.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((2.0/8.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/8.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real(0.0));//XLLL
	result.push_back(cxsc::real((2.0/8.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real(0.0));//XLLRL
	result.push_back(cxsc::real((2.0/8.0)*(16.0/rootvol)));//XLLRR
	result.push_back(cxsc::real(0.0));//XLR
	result.push_back(cxsc::real((6.0/8.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/8.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((5.0/8.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((4.0/8.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((3.0/8.0)*(16.0/rootvol)));//XRRLL
	result.push_back(cxsc::real((1.0/8.0)*(16.0/rootvol)));//XRRLR
	result.push_back(cxsc::real((1.0/8.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges2(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((2.0/2.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real(0.0));//XLL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real(0.0));//XLRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((1.0/2.0)*(16.0/rootvol)));//XLRRL
	result.push_back(cxsc::real(0.0));//XLRRR
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real(0.0));//XRL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real(0.0));//XRRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges3(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((3.0/3.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((1.0/3.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/3.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((1.0/3.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((1.0/3.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real(0.0));//XRRR
	
	return result;
	
}



bool checkFileLines(const std::string& s, std::size_t expectedLines)
{
	bool retValue = false;
	
	std::ifstream dataFile(s.c_str());
	
	if (dataFile.is_open())
	{
		std::size_t readNonBlankLines = 0;
		std::string line;
	
		while (dataFile.good() )
		{
			getline (dataFile,line);
			if (!line.empty()) readNonBlankLines++;
		}
		retValue = (readNonBlankLines == expectedLines);	
		
	}
	else { // dataFile not open
		std::cerr << "Unable to open file " << s << std::endl;
	}

	return retValue;
}

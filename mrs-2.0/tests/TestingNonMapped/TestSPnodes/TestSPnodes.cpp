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

/*!/ \file
\brief Testing SPnodes
*/

#include "spnode.hpp"
#include "sptypes.hpp"
#include "toolz.hpp"
#include "subpaving_exception.hpp"
#include "testing_tools.hpp"

// include fstream so as to be able to output a file
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>
#include <iterator>
#include <numeric>
#include <functional>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace subpavings;


RVecData& getData2(RVecData& data);
void checkNode(const SPnode * const spn);
void test();


int main()
{
	
	test();
	
	return 0;
}





void test()
{
	try {
		cout << "\nDefault constructor" << endl;
		
		subpavings::SPnode spn;
		
	}
	catch (std::exception const& e) {
		cout << "Exception\n" << e.what() << endl;
		
	}
	try {
		cout << "\nVolume default constructed node (this should fail)" << endl;
		
		subpavings::SPnode spn;
		spn.nodeVolume();

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& e) {
		cout << "Exception\n" << e.what() << endl;
		
	}
	
	try {
		cout << "\nTry to expand default node  (this should fail)" << endl;
		
		subpavings::SPnode spn;
		spn.nodeExpand();

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& e) {
		cout << "Exception\n" << e.what() << endl;
		
	}
	
	try {
		cout << "\nDemonstrate the perils of the default box  (this should fail)" << endl;
		
		cxsc::ivector fred;
		cout << "VecLen is " << VecLen(fred) << endl;
		cout << "Lb is " << Lb(fred) << " and Ub is " << Ub(fred) << endl;
		//referencing fred[1] will give a segmentation fault
		double vol = Volume(fred);

		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::exception const& e) {
		cout << "Exception thrown:\n" << e.what() << endl;
		
	}
	
	try {
		cout << "\nConstructor with default box  (this should fail)" << endl;
		
		cxsc::ivector fred;
		subpavings::SPnode fredowner(fred);
		
		throw std::logic_error("Should not be able to do this");
	}
	catch (std::exception const& e) {
		cout << "Exception " << e.what() << endl;
		
	}
	
	try	{
		cout << "\nTry 1-d box" << endl;
		
		int dim = 1;
		cxsc::ivector box(dim);
		cxsc::interval el(-2,2);
		for (int k = 0; k < dim; ++k) {
			box[k+1] = el;
		}
		
		subpavings::SPnode spn(box);
		
		assert(spn.getDimension() == dim);
		cout << "\nPassed assert that box has dim " << dim << endl;
		
	}
	catch (std::exception const& e) {
		cout << "Failed to construct 1-d box:\n" << e.what() << endl;
	}

	try	{
		cout << "\nTry 2-d box" << endl;
		
		int dim = 2;
		cxsc::ivector box(dim);
		cxsc::interval el(-2,2);
		for (int k = 0; k < dim; ++k) {
			box[k+1] = el;
		}
		
		subpavings::SPnode spn(box);
		
		assert(spn.getDimension() == dim);
		cout << "\nPassed assert that box has dim " << dim << endl;
		
	}
	catch (std::exception const& e) {
		cout << "Failed to construct 2-d box:\n" << e.what() << endl;
	}
	
	try	{
		cout << "\nDemonstrate what happens if we try to split a box that only has thin intervals (this should fail)" << endl;
		
		int dim = 2;
		cxsc::ivector box(dim);
		cxsc::interval el(2,2); // thin interval
		for (int k = 0; k < dim; ++k) {
			box[k+1] = el;
		}
		
		subpavings::SPnode spn(box);
		
		assert(spn.getDimension() == dim);
		
		spn.nodeExpand();
		throw std::logic_error("Should not be able to do this");

		
	}
	catch (std::exception const& e) {
		cout << "Failed to split thin interval box:\n" << e.what() << endl;
	}

	try	{
		
		
		int dim = 3;
		cxsc::ivector box(dim);
		cxsc::interval el(-2,2);
		for (int k = 0; k < dim; ++k) {
			box[k+1] = el;
		}
		
		
		subpavings::SPnode spn(box);
		
		spn.splitRootToShape("2,3,3,1");
		cout << "\n\n" << endl;

		{
			
				checkNode(&spn);
		}
		{
			SPnode* spnc= spn.getLeftChild();
			checkNode(spnc);
		}
		{
			SPnode* spnc= spn.getLeftChild()->getLeftChild();
			checkNode(spnc);
		}
		{
			SPnode* spnc= spn.getLeftChild()->getRightChild();
			checkNode(spnc);
		}
		{
			SPnode* spnc= spn.getLeftChild()->getRightChild()->getLeftChild();
			checkNode(spnc);
		}
		{
			SPnode* spnc= spn.getLeftChild()->getRightChild()->getRightChild();
			checkNode(spnc);
		}
		{
			SPnode* spnc= spn.getRightChild();
			checkNode(spnc);
		}
	}
	catch (std::exception const& e) {
		cout << e.what() << endl;
	} 
	
		// try to split to shape from vector
	{	
		int dim = 3;
		cxsc::ivector pavingBox1(dim);
		cxsc::interval el(-2,2);
		for (int k = 0; k < dim; ++k) {
			pavingBox1[k+1] = el;
		}
	
		try {
			
			cout << "\ntry splitRootAtLeastToShape, empty reqDepths" << endl;
			
			std::vector < size_t > reqDepths;
			
			SPnode spn(pavingBox1);
			
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
			std::string msg(ia.what());
			cout << "\nFailed to do splitRootAtLeastToShape, empty reqDepths:\n" << msg << endl;
		}
		
		try {
			
			cout << "\ntry splitRootAtLeastToShape, invalid reqDepths - too short" << endl;
			
			size_t depths[] = {1};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPnode spn(pavingBox1);
			
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
			cout << "\nFailed to do splitRootAtLeastToShape, invalid reqDepths:\n" << ia.what() << endl;
		}
		
		try {
			
			cout << "\ntry splitRootAtLeastToShape, invalid reqDepths " << endl;
			
			size_t depths[] = {1,2,2,3,4};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPnode spn(pavingBox1);
			
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
			cout << "\nFailed to do splitRootAtLeastToShape, invalid reqDepths:\n" << ia.what() << endl;
		}
		
		try {
			
			cout << "\ntry splitRootAtLeastToShape, invalid reqDepths " << endl;
			
			size_t depths[] = {1,2,2,1};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPnode spn(pavingBox1);
			
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
			cout << "\nFailed to do splitRootAtLeastToShape, invalid reqDepths:\n" << ia.what() << endl;
		}
		
			
		
		{
			
			cout << "\ntry splitRootAtLeastToShape, nodes too small" << endl;
			
			int dim = 1;
			cxsc::ivector box1(dim);
			cxsc::interval el(-cxsc::MinReal,cxsc::MinReal);
			for (int k = 0; k < dim; ++k) {
				box1[k+1] = el;
			}
			size_t depths[] = {1,2,3,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			SPnode spn(box1);
			
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (success) throw std::logic_error("Should have been able to do that successfully");
						
			cout << "After operation, success = " << success << " and leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;

			
		}
		
		
		{
			
			cout << "\nsplitRootAtLeastToShape valid reqDims" << endl;
			
			size_t depths[] = {1,2,3,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
	
			
			SPnode spn(pavingBox1);
						
			bool success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			cout << "After operation, leaf node levels string is " << spn.getLeafNodeLevelsString() << endl;
		}
		{
			
			cout << "\nsplitRootAtLeastToShape with less split than exists already" << endl;
			
			size_t depths1[] = {2,2,2,4,4,3};
			std::vector < size_t > reqDepths1 (depths1, depths1 + sizeof(depths1) / sizeof(size_t) );
	
			SPnode spn(pavingBox1);
						
			bool success = spn.splitRootAtLeastToShape(reqDepths1);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			std::string firstState = spn.getLeafNodeLevelsString();
			cout << "After first operation, leaf node levels string is " << firstState << endl;
			size_t depths[] = {1,2,3,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
			
			success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			std::string secondState = spn.getLeafNodeLevelsString();
			cout << "After second operation, leaf node levels string is " << secondState << endl;
			assert(firstState == secondState);
	
		}
				
			
		{
			
			cout << "\nsplitRootAtLeastToShape with exactly same split as exists already" << endl;
			
			size_t depths1[] = {2,2,2,4,4,3};
			std::vector < size_t > reqDepths1 (depths1, depths1 + sizeof(depths1) / sizeof(size_t) );
	
			SPnode spn(pavingBox1);
						
			bool success = spn.splitRootAtLeastToShape(reqDepths1);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			std::string firstState = spn.getLeafNodeLevelsString();
			cout << "After first operation, leaf node levels string is " << firstState << endl;
			size_t depths[] = {2,2,2,4,4,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
			
			success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			std::string secondState = spn.getLeafNodeLevelsString();
			cout << "After second operation, leaf node levels string is " << secondState << endl;
			assert(firstState == secondState);
	
		}
		
		{
			
			cout << "\nsplitRootAtLeastToShape with more splits in some places" << endl;
			
			size_t depths1[] = {2,2,2,4,4,3};
			std::vector < size_t > reqDepths1 (depths1, depths1 + sizeof(depths1) / sizeof(size_t) );
	
			SPnode spn(pavingBox1);
						
			bool success = spn.splitRootAtLeastToShape(reqDepths1);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			std::string firstState = spn.getLeafNodeLevelsString();
			cout << "After first operation, leaf node levels string is " << firstState << endl;
			size_t depths[] = {3,4,4,2,2,4,5,5,3};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
			
			success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			std::string secondState = spn.getLeafNodeLevelsString();
			cout << "After second operation, leaf node levels string is " << secondState << endl;
			assert(firstState.size() < secondState.size());
			
		}
		{
			
			cout << "\nsplitRootAtLeastToShape with more splits everywhere" << endl;
			
			size_t depths1[] = {2,2,2,4,4,3};
			std::vector < size_t > reqDepths1 (depths1, depths1 + sizeof(depths1) / sizeof(size_t) );
	
			SPnode spn(pavingBox1);
						
			bool success = spn.splitRootAtLeastToShape(reqDepths1);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			
			std::string firstState = spn.getLeafNodeLevelsString();
			cout << "After first operation, leaf node levels string is " << firstState << endl;
			size_t depths[] = {4,4,5,5,5,6,6,3,3,3,4,4,5,5,6,6,5,4,5,5};
			std::vector < size_t > reqDepths (depths, depths + sizeof(depths) / sizeof(size_t) );
			
			success = spn.splitRootAtLeastToShape(reqDepths);
			
			if (!success) throw std::logic_error("Should have been able to do that");
			std::string secondState = spn.getLeafNodeLevelsString();
			cout << "After second operation, leaf node levels string is " << secondState << endl;
			assert(firstState.size() < secondState.size());
			
		}
	}
	
	
	cout << "\n\nEnd of spnode tests\n" << endl;

}


void checkNode(const SPnode * const spn)
{
	int depth = spn->getTreeHeight();
	int nodedepth = spn->getNodeDepth();
	double smallestVol = spn->getSmallestLeafVol();
	double largestVol = spn->getLargestLeafVol();

	cout << spn->getNodeName() << " nodedepth is " << nodedepth << " and tree height is " << depth << endl;
	cout << "isLeaf " << (spn->isLeaf()? "Yes" : "No") << " and has leaf sibling " << (spn->hasLeafSibling()? "Yes" : "No") << endl;
	cout <<  "leafNodeLevelsStrng is " << spn->getLeafNodeLevelsString() << endl;
	cout <<  "largestLeafVol is " << largestVol << " and smallestLeafVol is " << smallestVol << endl;
	cout <<  endl;
}


RVecData& getData2(RVecData& data) 
{
    int d = 2;
	{ 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -1;
        thisrv[2] = 1;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = 5;
		data.push_back(thisrv);
	}
	{ 
		rvector thisrv(d);
		thisrv[1] = -5;
        thisrv[2] = -5;
		data.push_back(thisrv);
	}
	
	return data;
}

/*
* Copyright 2012 Jennifer Harlow
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
\brief RealMappedShapeNode and associated non-member function definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "realmappedshapenode.hpp"

#include "subpaving_exception.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <stdexcept>

#include <numeric> // for inner_product

#include <iterator> // ostream iterator 

#include <cmath>

//#include <cassert> // assertions

//#define MYDEBUG

#ifdef NDEBUG
	#undef MYDEBUG
#endif

using namespace std;
using namespace subpavings;


// ----------------------- RealMappedShapeNode class definitions ------------------

// ---------------- public methods



// constructor with backtransformation
RealMappedShapeNode::RealMappedShapeNode(const RealMappedSPnode& rmsp,
		const std::vector < std::vector < double > >& backtransform,
		const std::vector < double >& scale,
		const std::vector < double >& shift)
		:  verticesPtr(NULL), range(rmsp.getRange()),
		nodeName(rmsp.getNodeName()),
		parent(NULL), leftChild(NULL), rightChild(NULL)
{
	
	if (rmsp.isEmpty()) {
		throw NoBox_Error("RealMappedShapeNode::RealMappedShapeNode(...)");
	}
	try {
		
		try {
			ivector box = rmsp.getBox();
			extractVertices(box, backtransform);
			scaleAndShiftShape(scale, shift);
		}
		catch (std::invalid_argument& ia) {
			throw std::invalid_argument(
				string("RealMappedShapeNode::RealMappedShapeNode(...) :\n")
					+ ia.what() );
		}
		
		#ifdef MYDEBUG
			// output the name, range
			
			cout << nodeName << "\t" << range << "\t";

			// followed by vertices all on one line
			for (std::vector < std::vector < real > >::const_iterator 
						it = verticesPtr->begin();
						it < verticesPtr->end(); ++it) {
				
				 ostream_iterator< real > out_it (cout,"\t");
				 copy ( it->begin(), it->end(), out_it );
			}
			cout << endl;
		#endif
		//recursion on the children
		if (rmsp.getLeftChild()) {
			leftChild=new RealMappedShapeNode(*rmsp.getLeftChild(),
												backtransform, 
												scale, shift);
			leftChild->parent = this;
		}
		else leftChild=NULL;

		if (rmsp.getRightChild()) {
			rightChild=new RealMappedShapeNode(*rmsp.getRightChild(),
												backtransform, 
												scale, shift);
			rightChild->parent = this;
		}
		else rightChild=NULL;
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}

// constructor with no backtransformation
RealMappedShapeNode::RealMappedShapeNode(const RealMappedSPnode& rmsp,
		const std::vector < double >& scale,
		const std::vector < double >& shift)
		:  verticesPtr(NULL), range(rmsp.getRange()),
		nodeName(rmsp.getNodeName()),
		parent(NULL), leftChild(NULL), rightChild(NULL)
		
{
	
	if (rmsp.isEmpty()) {
		throw NoBox_Error("RealMappedShapeNode::RealMappedShapeNode(...)");
	}
	try {
		
		try {
			ivector box = rmsp.getBox();
			extractVertices(box);
			scaleAndShiftShape(scale, shift);
		}
		catch (std::invalid_argument& ia) {
			throw std::invalid_argument(
				string("RealMappedShapeNode::RealMappedShapeNode(...) :\n")
					+ ia.what() );
		}
		
		#ifdef MYDEBUG
			// output the name, range
			
			cout << nodeName << "\t" << range << "\t";

			// followed by vertices all on one line
			for (std::vector < std::vector < real > >::const_iterator 
						it = verticesPtr->begin();
						it < verticesPtr->end(); ++it) {
				
				 ostream_iterator< real > out_it (cout,"\t");
				 copy ( it->begin(), it->end(), out_it );
			}
			cout << endl;
		#endif
		//recursion on the children
		if (rmsp.getLeftChild()) {
			leftChild=new RealMappedShapeNode(*rmsp.getLeftChild(),
												scale, shift);
			leftChild->parent = this;
		}
		else leftChild=NULL;

		if (rmsp.getRightChild()) {
			rightChild=new RealMappedShapeNode(*rmsp.getRightChild(),
												scale, shift);
			rightChild->parent = this;
		}
		else rightChild=NULL;
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}

// constructor with box and range
RealMappedShapeNode::RealMappedShapeNode(const cxsc::ivector& box, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				real r)	:  verticesPtr(NULL), range(r),
					nodeName("X"),
					parent(NULL), leftChild(NULL), rightChild(NULL)
					
{
	
	try {
		
		try {
			extractVertices(box, backtransform);
			scaleAndShiftShape(scale, shift);
		}
		catch (std::invalid_argument& ia) {
			throw std::invalid_argument(
				string("RealMappedShapeNode::RealMappedShapeNode(...) :\n")
					+ ia.what() );
		}
		
		#ifdef MYDEBUG
			// output the name, range
			
			cout << nodeName << "\t" << range << "\t";

			// followed by vertices all on one line
			for (std::vector < std::vector < real > >::const_iterator 
						it = verticesPtr->begin();
						it < verticesPtr->end(); ++it) {
				
				 ostream_iterator< real > out_it (cout,"\t");
				 copy ( it->begin(), it->end(), out_it );
			}
			cout << endl;
		#endif
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}


// constructor with box and range
RealMappedShapeNode::RealMappedShapeNode(const cxsc::ivector& box, 
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				real r)	:  verticesPtr(NULL), range(r),
					nodeName("X"),
					parent(NULL), leftChild(NULL), rightChild(NULL)
					
{
	
	try {
		
		try {
			extractVertices(box);
			scaleAndShiftShape(scale, shift);
		}
		catch (std::invalid_argument& ia) {
			throw std::invalid_argument(
				string("RealMappedShapeNode::RealMappedShapeNode(...) :\n")
					+ ia.what() );
		}
		
		#ifdef MYDEBUG
			// output the name, range
			
			cout << nodeName << "\t" << range << "\t";

			// followed by vertices all on one line
			for (std::vector < std::vector < real > >::const_iterator 
						it = verticesPtr->begin();
						it < verticesPtr->end(); ++it) {
				
				 ostream_iterator< real > out_it (cout,"\t");
				 copy ( it->begin(), it->end(), out_it );
			}
			cout << endl;
		#endif
		
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}


//Copy constructor
//copies from the given node downwards
RealMappedShapeNode::RealMappedShapeNode(const RealMappedShapeNode& other) : 
		verticesPtr(NULL), range(other.getRange()), 
		nodeName(other.nodeName),
		parent(NULL), 
		leftChild(NULL), rightChild(NULL)
		
{
	
	try {
		// deep copy of the vertices 
		verticesPtr = 
			new std::vector < std::vector < real > >(*(other.verticesPtr));

		//recursion on the children
		if (other.leftChild) {
			leftChild=new RealMappedShapeNode(*other.leftChild);
			leftChild->parent = this;
		}
		else leftChild=NULL;

		if (other.rightChild) {
			rightChild=new RealMappedShapeNode(*other.rightChild);
			rightChild->parent = this;
		}
		else rightChild=NULL;
	}
	catch (exception const& e) {
		constructor_error_handler();
	}
}


// Destructor.
RealMappedShapeNode::~RealMappedShapeNode()
{
	
	try {
		delete verticesPtr;
		verticesPtr = NULL;
		delete leftChild;
		leftChild = NULL;
		delete rightChild;
		rightChild = NULL;
	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in RealMappedShapeNode destructor:\n" << ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}


//copy assignment operator
//copies from this node downwards
// no check for self-assignment - copy elision and copy-and-swap idiom
RealMappedShapeNode& RealMappedShapeNode::operator=(RealMappedShapeNode rhs)
{
	rhs.swap(*this);
	return *this;
}



size_t RealMappedShapeNode::getDimension() const
{
	size_t d = 1;
	
	size_t v = 2;
	
	size_t s = verticesPtr->size();
	
	while (v < s) {
		++d;
		v *= 2;
	}
	
	return d;
}


real RealMappedShapeNode::getRange() const
{
	return range;
}


std::vector < std::vector < real > >
				RealMappedShapeNode::getVertices() const
{
	return *verticesPtr;
}


// Accessor for the parent of a node.
//Returns a copy of the pointer to parent node.
RealMappedShapeNode* RealMappedShapeNode::getParent() const
{ return parent; }

// Accessor for the left child of a node.
// Returns a copy of the pointer to leftChild node.
RealMappedShapeNode* RealMappedShapeNode::getLeftChild() const
{ return leftChild; }

// Accessor for the right child of a node.
// Returns a copy of the pointer to rightChild node.
RealMappedShapeNode* RealMappedShapeNode::getRightChild() const
{ return rightChild; }



// Get the node name.
std::string RealMappedShapeNode::getNodeName() const
{
	return nodeName;
}





// Check if this RealMappedShapeNode is a leaf.
bool RealMappedShapeNode::isLeaf() const
{return ( (NULL == leftChild) && (NULL == rightChild)); }



//Output for all the  leaf boxes in this, using tab delimiters
std::ostream& RealMappedShapeNode::leavesOutputTabs(std::ostream &os) const
{
	// uses  member function leafOutputTabs to generate node output
	if ( isLeaf() ) { 
		leafOutputTabs(os);
		os << "\n";

	}

	//recurse on the children
	if (getLeftChild()!=NULL) {
		getLeftChild()->leavesOutputTabs(os);
	}

	if (getRightChild()!=NULL) {
		getRightChild()->leavesOutputTabs(os);
	}
	return os;

}

//Output for all the  leaf boxes in this, using tab delimiters
std::ostream& RealMappedShapeNode::leavesOutputTabs(std::ostream &os, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	leavesOutputTabs(os);
	os << cxsc::RestoreOpt;
	
	return os;

}

//Output for all the nodes in this, using tab delimiters
std::ostream& RealMappedShapeNode::nodesAllOutput(std::ostream &os,
									int level) const
{

	leafOutputTabs(os);
	os << "\n";
	
	if(getLeftChild()!=NULL) {
		for (int i = 0; i < level; i++) {
			os << "\t";
		}
		getLeftChild()->nodesAllOutput(os, (level+1));
	}

	if(getRightChild()!=NULL) {
		for (int i = 0; i < level; i++) {
			os << "\t";
		}
		getRightChild()->nodesAllOutput(os, (level+1));
	}

	return os;
}

//Output for all the nodes in this
std::ostream& RealMappedShapeNode::nodesAllOutput(std::ostream &os,
									int level, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	nodesAllOutput(os, level);
	os << cxsc::RestoreOpt;
	
	return os;

}



void RealMappedShapeNode::swap (RealMappedShapeNode& rmsn) //throw()
{
	std::swap(verticesPtr, rmsn.verticesPtr); 
	std::swap(range, rmsn.range);
	std::swap(nodeName, rmsn.nodeName);
	std::swap(parent, rmsn.parent);
	// can just swap child pointers
	std::swap(leftChild, rmsn.leftChild);
	std::swap(rightChild, rmsn.rightChild);
	// children have to be repointed to the new swapped parents
	if (leftChild != NULL) leftChild->parent = this;
	if (rightChild != NULL) rightChild->parent = this;
	if (rmsn.leftChild != NULL) rmsn.leftChild->parent = &rmsn;
	if (rmsn.rightChild != NULL) rmsn.rightChild->parent = &rmsn;
	
}

// -----------------private methods

void RealMappedShapeNode::scaleAndShiftShape(
			const std::vector < double >& scale,
			const std::vector < double >& shift)
{
	
	size_t d = verticesPtr->back().size();
	
	if ( scale.size() != d ) {
		throw std::invalid_argument("scale dimensions");
	}
	if ( shift.size() != d ) {
		throw std::invalid_argument("shift dimensions");
	}
	if ( !(*std::min_element(scale.begin(), scale.end()) > 0.0) ) {
		throw std::invalid_argument("scalar <= 0.0");
	}
	
	/* apply scale to the vertices and also adjust the range
	 * then apply the shift */
	#ifdef MYDEBUG
		// output the name, range
		cout << "\n before scaling and shifting" << endl;
		
		cout << nodeName << "\t" << range << "\t";

		// followed by vertices all on one line
		for (std::vector < std::vector < real > >::const_iterator 
					it = verticesPtr->begin();
					it < verticesPtr->end(); ++it) {
			
			 ostream_iterator< real > out_it (cout,"\t");
			 copy ( it->begin(), it->end(), out_it );
		}
		cout << endl;
	#endif
	
	size_t v = verticesPtr->size();
	
	/* mult ith element in every vertex by ith element in scale */
	for (size_t i = 0; i < d; ++i) {
		
		double scalar = scale[i];
		double shifter = shift[i];
		
		range /= scalar; // adjust range for scalar on this dimension
		
		for (size_t j = 0; j < v; ++j) {
			
			verticesPtr->at(j)[i] *= scalar;
			verticesPtr->at(j)[i] += shifter;
		}
	
	}

}

void RealMappedShapeNode::extractVertices(const ivector& box,
			const std::vector < std::vector < double > >& backtransform)
{
	int d = VecLen(box);
	
	if ( d < 0 ) {
		throw std::invalid_argument("backtransform dimensions");
	}
	
	if ( (backtransform.size() != static_cast<size_t>(d)) || 
			(backtransform.back().size() != static_cast<size_t>(d)) ) {
		throw std::invalid_argument("backtransform dimensions");
	}
	
	// 2^d vertices, all d long
	std::vector < std::vector < real > > tmp;
	
	size_t v = static_cast<size_t>( std::pow(2.0, d)); // number of vertices
	
	tmp.reserve(v);
	
	std::vector < real > vertex(d, 0.0);
	
	buildVertices(box, d, tmp, vertex);	
	
	// not the final content but gets vertices to the right size
	verticesPtr = new std::vector < std::vector < real > >(tmp);
	
	
	#ifdef MYDEBUG
			// output the name, range
			cout << "\n before back transformation" << endl;
			
			cout << nodeName << "\t" << range << "\t";

			// followed by vertices all on one line
			for (std::vector < std::vector < real > >::const_iterator 
						it = tmp.begin();
						it <tmp.end(); ++it) {
				
				 ostream_iterator< real > out_it (cout,"\t");
				 copy ( it->begin(), it->end(), out_it );
			}
			cout << endl;
		#endif
	
	if ( d > 1) { // backtransformation only relevant if d > 1
		for (int i = 0; i < d; ++i) {
			
			std::vector < real > row(backtransform[i].begin(), backtransform[i].end());
			/* ith row of backtransform dot product jth tmp vertex
			 * becomes ith value in jth vertex */ 
			for (size_t j = 0; j < v; ++j) {
				
				verticesPtr->at(j)[i] = 
						std::inner_product(row.begin(), row.end(),
											tmp[j].begin(),	real(0.0));
			}
		
		}
	}
	
}


void RealMappedShapeNode::extractVertices(const ivector& box)
{
	int d = VecLen(box);
	
	// 2^d vertices, all d long
	std::vector < std::vector < real > > tmp;
	
	size_t v = static_cast<size_t>( std::pow(2.0, d)); // number of vertices
	
	tmp.reserve(v);
	
	std::vector < real > vertex(d, 0.0);
	
	buildVertices(box, d, tmp, vertex);	
	
	// not the final content but gets vertices to the right size
	verticesPtr = new std::vector < std::vector < real > >(tmp);
	
}

std::vector < std::vector < real > >& 
			RealMappedShapeNode::buildVertices(const ivector& box, 
							int d,
							std::vector < std::vector < real > >& verts,
							std::vector < real > vertex)
{
	
	std::vector < real > vertexCopy(vertex);
	vertex[d-1] = cxsc::Inf(box[d]);
	vertexCopy[d-1] = cxsc::Sup(box[d]);
		
	if (d > 1) {
		
		buildVertices(box, d-1, verts, vertex);
		buildVertices(box, d-1, verts, vertexCopy);
	}
	else { // d ==1
		verts.push_back(vertex);
		verts.push_back(vertexCopy);
	}
	
	return verts;
}	





// add lChild onto this node
void RealMappedShapeNode::nodeAddLeft(RealMappedShapeNode *lChild)
{
	if (lChild != NULL) {
		
		this->leftChild = lChild;
		leftChild->parent = this;
	}
}

// add Child onto this node
void RealMappedShapeNode::nodeAddRight(RealMappedShapeNode *rChild)
{
	if (rChild != NULL) {
				
		this->rightChild = rChild;
		rightChild->parent = this;
	}
}


// Print the details of a single leaf node, using tab delimiters
std::ostream& RealMappedShapeNode::leafOutputTabs(std::ostream &os) const
{

	// output the name, range
	os << nodeName << "\t" << range << "\t";

	// followed by vertices all on one line
	for (std::vector < std::vector < real > >::const_iterator it = verticesPtr->begin();
				it < verticesPtr->end(); ++it) {
		
		 ostream_iterator< real > out_it (os,"\t");
		 copy ( it->begin(), it->end(), out_it );
	}
	return os;

}

// Print the details of a single leaf node, using tab delimiters
std::ostream& RealMappedShapeNode::leafOutputTabs(std::ostream &os, int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);

	leafOutputTabs(os);
		
	os << cxsc::RestoreOpt;
	return os;

}

// ------------------- private methods ------------------------------

// ensure theBox and children (if any) are deleted if constructed in failed constructor
void RealMappedShapeNode::constructor_error_handler() 
{
	try {
			delete verticesPtr;
			verticesPtr = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	try {
			delete leftChild;
			leftChild = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	try {
			delete rightChild;
			rightChild = NULL;
	}
	catch (std::exception& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}

// ------------------- end of RealMappedShapeNode class definitions ---------------




// -------------------- start of RealMappedShapeNode non-member functions



//Output the leaf boxes in SubPaving spn
std::ostream & subpavings::operator<<(std::ostream &os, const RealMappedShapeNode& rpn)
{
	rpn.leavesOutputTabs(os);
	
	return os;
}





// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::RealMappedShapeNode & s1, 
			subpavings::RealMappedShapeNode & s2) // throw ()
{
	s1.swap(s2);
}

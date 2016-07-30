/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief PolygonPiecewiseConstantFunction definitions
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "polygon_piecewise_function.hpp"

#include "subpaving_exception.hpp"

#include <iostream> // to use standard input and output
#include <string>   // to use the C++ string class
#include <fstream>  // for ifstream, ofstream
//#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // use exceptions
#include <cassert>



using namespace subpavings;
using namespace std;



// -------------------implementation of PolygonPiecewiseConstantFunction class --------------


// ----------- public methods

PolygonPiecewiseConstantFunction::PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift)
{
    try {
		
		RealMappedSPnode::ConstPtrs pieces;
		std::vector < cxsc::ivector > boxes;
		pcf.getPieces(pieces); 
		boxes.push_back(pcf.getRootBox());     
        
		fillNodes(pieces, boxes, backtransform, scale, shift);
				
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

PolygonPiecewiseConstantFunction::PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < double >& scale,
				const std::vector < double >& shift)
{
    try {
        
        RealMappedSPnode::ConstPtrs pieces;
		std::vector < cxsc::ivector > boxes;
		pcf.getPieces(pieces); 
		boxes.push_back(pcf.getRootBox());       
        
		fillNodes(pieces, boxes, scale, shift);
		
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}


PolygonPiecewiseConstantFunction::PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real cov)
{
    try {
        
        RealMappedSPnode::ConstPtrs pieces;
		std::vector < cxsc::ivector > boxes;
		pcf.findCoverageRegion(pieces, boxes, cov);        
        
		fillNodes(pieces, boxes, backtransform, scale, shift);
		
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}

PolygonPiecewiseConstantFunction::PolygonPiecewiseConstantFunction(
				const PiecewiseConstantFunction& pcf, 
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real cov)
{
    try {
        
        RealMappedSPnode::ConstPtrs pieces;
		std::vector < cxsc::ivector > boxes;
		pcf.findCoverageRegion(pieces, boxes, cov); 
		
		fillNodes(pieces, boxes, scale, shift);
		
	}
    catch (exception const& e) {
		constructor_error_handler();
    }
}


// copy constructor*/
PolygonPiecewiseConstantFunction::PolygonPiecewiseConstantFunction(
								const PolygonPiecewiseConstantFunction& other)
{
    try {
		
		nodes.reserve(other.nodes.size());
		
		for (RealMappedShapeNode::PtrsConstItr it = other.nodes.begin();
			it < other.nodes.end();
			++it) {
				
			nodes.push_back(new RealMappedShapeNode((**it)));	
		}	
		areaNodes.reserve(other.areaNodes.size());
		
		for (RealMappedShapeNode::PtrsConstItr it = other.areaNodes.begin();
			it < other.areaNodes.end();
			++it) {
				
			areaNodes.push_back(new RealMappedShapeNode((**it)));	
		}	
	}
    catch (exception const& e) {
		constructor_error_handler();
	}

}


//Destructor
PolygonPiecewiseConstantFunction::~PolygonPiecewiseConstantFunction()
{
	try {
		for (RealMappedShapeNode::PtrsItr it = nodes.begin();
								it < nodes.end();
								++it) {
			delete (*it);
			(*it) = NULL;
		}
		for (RealMappedShapeNode::PtrsItr it = areaNodes.begin();
								it < areaNodes.end();
								++it) {
			delete (*it);
			(*it) = NULL;
		}

	}
	catch (exception const& e) {
		try {
			constructor_error_handler();
		}
		catch(std::exception const& ee) {
			std::cerr << "Error in PolygonPiecewiseConstantFunction destructor:\n" 
										<< ee.what() << std::endl;
		}
	} // exceptions ultimately swallowed
}

// assignment operator
PolygonPiecewiseConstantFunction& PolygonPiecewiseConstantFunction::operator=(
    PolygonPiecewiseConstantFunction rhs)
{
	rhs.swap(*this);
	return *this;

}


int PolygonPiecewiseConstantFunction::getDimensions() const
{
	 return nodes.front()->getDimension();
	
}



// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void PolygonPiecewiseConstantFunction::outputToTxtTabs(const std::string& s,
                            int prec) const
{
	outputToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void PolygonPiecewiseConstantFunction::outputToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the PolygonPiecewiseConstantFunction object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		for (RealMappedShapeNode::PtrsConstItr it = nodes.begin();
			it < nodes.end();
			++it) {
				
			(*it)->leavesOutputTabs(os, prec); 
		}	
			
		if (confirm)
			std::cout << "The output of the PolygonPiecewiseConstantFunction "
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
std::ostream& PolygonPiecewiseConstantFunction::outputToStreamTabs(
									std::ostream& os, int prec) const
{

	for (RealMappedShapeNode::PtrsConstItr it = nodes.begin();
		it < nodes.end();
		++it) {
			
		(*it)->leavesOutputTabs(os, prec); 
	}	
	return os;
		
}


// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void PolygonPiecewiseConstantFunction::outputAreaToTxtTabs(const std::string& s,
                            int prec) const
{
	outputAreaToTxtTabs(s, prec, false);
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
void PolygonPiecewiseConstantFunction::outputAreaToTxtTabs(const std::string& s,
                            int prec, bool confirm) const
{

	// To generate a file output of the PolygonPiecewiseConstantFunction object
	ofstream os(s.c_str());         // Filename, c-string version
	if (os.is_open()) {
		
		for (RealMappedShapeNode::PtrsConstItr it = areaNodes.begin();
			it < areaNodes.end();
			++it) {
				
			(*it)->leavesOutputTabs(os, prec); 
		}	
			
		if (confirm)
			std::cout << "The output of the PolygonPiecewiseConstantFunction domain"
				<< "has been written to " << s << std::endl << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to output the subpaving to a txt file
// Output goes to file named according to argument s
std::ostream& PolygonPiecewiseConstantFunction::outputAreaToStreamTabs(
									std::ostream& os, int prec) const
{

	for (RealMappedShapeNode::PtrsConstItr it = areaNodes.begin();
		it < areaNodes.end();
		++it) {
			
		(*it)->leavesOutputTabs(os, prec); 
	}	
	return os;
		
}



void PolygonPiecewiseConstantFunction::swap(PolygonPiecewiseConstantFunction& ppcf) // throw()
{
	std::swap(nodes, ppcf.nodes); // just swap the collections of nodes
	std::swap(areaNodes, ppcf.areaNodes); // just swap the collections of nodes
	
}


// --------------------------- private ---------------------------------------

void PolygonPiecewiseConstantFunction::fillNodes(
				const RealMappedSPnode::ConstPtrs& pieces,
				const std::vector < cxsc::ivector >& boxes, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift)
{
    nodes.reserve(pieces.size());
		
	for (RealMappedSPnode::ConstPtrsItr it = pieces.begin();
			it < pieces.end();
			++it) {
				
		nodes.push_back( new RealMappedShapeNode((**it), 
								backtransform, scale, shift) );		
	}
	
	areaNodes.reserve(boxes.size());
		
	for (std::vector < cxsc::ivector >::const_iterator it = boxes.begin();
			it < boxes.end();
			++it) {
				
		areaNodes.push_back( new RealMappedShapeNode(*it, 
								backtransform, scale, shift) );		
	}
	
}

void PolygonPiecewiseConstantFunction::fillNodes(
				const RealMappedSPnode::ConstPtrs& pieces, 
				const std::vector < cxsc::ivector >& boxes,
				const std::vector < double >& scale,
				const std::vector < double >& shift)
{
	nodes.reserve(pieces.size());
		
	for (RealMappedSPnode::ConstPtrsItr it = pieces.begin();
			it < pieces.end();
			++it) {
				
		nodes.push_back( new RealMappedShapeNode((**it), 
								scale, shift) );		
	}
	
	areaNodes.reserve(boxes.size());
		
	for (std::vector < cxsc::ivector >::const_iterator it = boxes.begin();
			it < boxes.end();
			++it) {
				
		areaNodes.push_back( new RealMappedShapeNode(*it, 
								scale, shift) );		
	}
}


// ensure rootPaving is deleted if constructed in failed constructor
void PolygonPiecewiseConstantFunction::constructor_error_handler() 
{
	try {
		
			for (RealMappedShapeNode::PtrsItr it = nodes.begin();
								it < nodes.end();
								++it) {
				delete (*it);
				(*it) = NULL;
			}
			for (RealMappedShapeNode::PtrsItr it = areaNodes.begin();
								it < areaNodes.end();
								++it) {
				delete (*it);
				(*it) = NULL;
			}
			
	}
	catch (std::exception const& ee) {} // catch and swallow
	
	throw; // rethrow the original exception
}


// ----------------------------- non member functions

//Output all boxes in PolygonPiecewiseConstantFunction pcf
std::ostream & subpavings::operator<<(std::ostream &os, 
				const subpavings::PolygonPiecewiseConstantFunction& ppcf)
{
    ppcf.outputToStreamTabs(os);
    return os;
}


// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap (subpavings::PolygonPiecewiseConstantFunction & p1, 
		subpavings::PolygonPiecewiseConstantFunction & p2) // throw ()
{
	p1.swap(p2);
}







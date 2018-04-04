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
\brief SubpavingException definitions
*/

#include "subpaving_exception.hpp"

#include <iostream>


using namespace subpavings;
using namespace std;


// ----------------------------- subpaving exceptions definitions

IO_Error::IO_Error(std::string ss) 
	: std::runtime_error("subpavings::IO_Error: file " + ss) {}

IO_Error::~IO_Error() throw() {}

const char* IO_Error::what() const throw() 
{
	try {
		return std::runtime_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "IO_Error::Problem in what():\n" << (e.what()) << std::endl;
	}
	
}



NullSubpavingPointer_Error::NullSubpavingPointer_Error(std::string ss) 
	: std::logic_error("subpavings::NullSubpavingPointer_Error: see " + ss) {}

NullSubpavingPointer_Error::~NullSubpavingPointer_Error() throw() {}

const char* NullSubpavingPointer_Error::what() const throw() 
{
	try {
				
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in NullSubpavingPointer_Error::what():\n" << (e.what()) << std::endl;
	}
} 



NoBox_Error::NoBox_Error(std::string ss) 
	: std::logic_error("subpavings::NoBox_Error: see " + ss) {}

NoBox_Error::~NoBox_Error() throw() {}

const char* NoBox_Error::what() const throw() 
{
	try {
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in NoBox_Error::what():\n" << (e.what()) << std::endl;
	}
}


MalconstructedBox_Error::MalconstructedBox_Error(std::string ss) 
	: std::logic_error("subpavings::Malconstructed_Error: see " + ss) {}

MalconstructedBox_Error::~MalconstructedBox_Error() throw() {}

const char* MalconstructedBox_Error::what() const throw() 
{
	try {
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in MalconstructedBox_Error::what():\n" << (e.what()) << std::endl;
	}
}


IncompatibleDimensions_Error::IncompatibleDimensions_Error(std::string ss) 
	: std::logic_error("subpavings::IncompatibleDimensions_Error: see " + ss) {}

IncompatibleDimensions_Error::~IncompatibleDimensions_Error() throw() {}

const char* IncompatibleDimensions_Error::what() const throw() 
{
	try {
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in IncompatibleDimensions_Error::what(): " << (e.what()) << std::endl;
	}
}



IncompatibleLabel_Error::IncompatibleLabel_Error(std::string ss) 
	: std::logic_error("subpavings::IncompatibleLabel_Error: see " + ss) {}

IncompatibleLabel_Error::~IncompatibleLabel_Error() throw() {}

const char* IncompatibleLabel_Error::what() const throw() 
{
	try {
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in IncompatibleLabel_Error::what():\n" << (e.what()) << std::endl;
	}
}



NonRootNode_Error::NonRootNode_Error(std::string ss) 
	: std::logic_error("subpavings::NonRootNode_Error: see " + ss) {}

NonRootNode_Error::~NonRootNode_Error() throw() {}

const char* NonRootNode_Error::what() const throw() 
{
	try {
		return std::logic_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in NonRootNode_Error::what(): " << (e.what()) << std::endl;
	}
}


UnfulfillableRequest_Error::UnfulfillableRequest_Error(std::string ss) 
	: std::runtime_error("subpavings::UnfulfillableRequest_Error: see " + ss) {}

UnfulfillableRequest_Error::~UnfulfillableRequest_Error() throw() {}

const char* UnfulfillableRequest_Error::what() const throw() 
{
	try {
		return std::runtime_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "Problem in UnfulfillableRequest_Error::what(): " << (e.what()) << std::endl;
	}
}

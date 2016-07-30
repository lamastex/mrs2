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
\brief SubpavingException declarations.
*/

#ifndef __SUBPAVING_EXCEPTION_HPP__
#define __SUBPAVING_EXCEPTION_HPP__

#include <stdexcept>
#include <string>     

namespace subpavings {



	/*! A class for exceptions raised by input/output errors.
	*/
	class IO_Error : public std::runtime_error
	{
	   
	   public:
	   IO_Error(std::string ss);
	   virtual ~IO_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	/*! A class for exceptions raised because a NULL subpaving pointer 
	 * is encountered.
	*/
	class NullSubpavingPointer_Error : public std::logic_error
	{
	   
	   public:
	   NullSubpavingPointer_Error(std::string ss);
	   virtual ~NullSubpavingPointer_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	/*! A class for exceptions raised because a subpaving has
	 * no box.
	*/
	class NoBox_Error : public std::logic_error
	{
	   
	   public:
	   NoBox_Error(std::string ss);
	   virtual ~NoBox_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	/*! A class for exceptions raised because a box is malconstructed.
	 * 
	 * Examples of malconstruction include a box with negative dimensions
	 * (see the no-args constructor of the cxsc::ivector) or a box
	 * where the diameter on at least one dimension is a thin interval,
	 * ie on that dimension the upper bound of the box equals 
	 * the lower bound.
	*/
	class MalconstructedBox_Error : public std::logic_error
	{
	   
	   public:
	   MalconstructedBox_Error(std::string ss);
	   virtual ~MalconstructedBox_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	/*! A class for exceptions raised because the dimensions
	 * of a subpaving box given as, or as part of, an argument or
	 * operand are incompatible with those required or accepted.
	*/
	class IncompatibleDimensions_Error : public std::logic_error
	{
	   
	   public:
	   IncompatibleDimensions_Error(std::string ss);
	   virtual ~IncompatibleDimensions_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	/*! A class for exceptions raised because the label 
	 associated with an object given as, or as part of, an argument or
	 operand is incompatible with that required or accepted.
	*/
	class IncompatibleLabel_Error : public std::logic_error
	{
	   
	   public:
	   IncompatibleLabel_Error(std::string ss);
	   virtual ~IncompatibleLabel_Error() throw();
	   virtual const char* what() const throw();
	   
	};
	
	
	/*! A class for exceptions raised because a request that can 
	 only be performed on a root node is made on a non-root node.
	*/
	class NonRootNode_Error : public std::logic_error
	{
	   
	   public:
	   NonRootNode_Error(std::string ss);
	   virtual ~NonRootNode_Error() throw();
	   virtual const char* what() const throw();
	   
	};

	/*! A class for exceptions raised because the request cannot be 
	performed on the node given the state at runtime:  for example,
	asking for the bestMerge change on a leaf node.  
	*/
	class UnfulfillableRequest_Error : public std::runtime_error
	{
	   
	   public:
	   UnfulfillableRequest_Error(std::string ss);
	   virtual ~UnfulfillableRequest_Error() throw();
	   virtual const char* what() const throw();
	   
	};

} // end namespace subpavings


#endif



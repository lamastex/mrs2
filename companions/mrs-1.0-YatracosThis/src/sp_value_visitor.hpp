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

#ifndef __SPVALUEVISITOR_HPP__
#define __SPVALUEVISITOR_HPP__


/*! \file
\brief declarations for SPValueVisitor

An interface for a type that visits \linkSPnode SPnodes\endlink 
and gives a value as the return value of the visit.
*/


namespace subpavings {

    class SPnode;

    template <typename T>
    class SPValueVisitor {

        
        public:
	
			virtual ~SPValueVisitor(){};

			//*! \brief The visit operation.*/
            virtual T visit(SPnode * spn) const = 0;
                       

    };
    // end of SPValueVisitor class

} // end namespace subpavings

#endif

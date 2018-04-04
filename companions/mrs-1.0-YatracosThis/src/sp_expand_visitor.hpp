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

#ifndef __SP_EXPANDER_HPP__
#define __SP_EXPANDER_HPP__


/*! \file
\brief declarations for SPExpandVisitor

A interface for a type that visits \linkSPnode SPnodes\endlink 
and expands them.*/



namespace subpavings {

    class SPnode;

	template <typename T>
    class SPExpandVisitor {

        
        public:
	
			virtual ~SPExpandVisitor(){};
	
			/*! \brief The visit operation.*/
            virtual T visit(SPnode * spn) const = 0;


    };
    // end of SPExpandVisitor class

} // end namespace subpavings

#endif

/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
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
\brief definitions for SPnodeVisitor.

An interface for a type that visits SPnodes

*/

#ifndef ___SPNODEVISITOR_HPP__
#define ___SPNODEVISITOR_HPP__

#include "real.hpp"
#include <gsl/gsl_rng.h>        // to know about the gsl random number generator
#include <vector>

namespace subpavings {

    class SPnode;

    class SPnodeVisitor {

        public:

            SPnodeVisitor() {};

            virtual void visit(SPnode * spn) {};

            virtual cxsc::real tellMe(SPnode * spn) {};

				//gat41
				virtual bool priorityVisit(SPnode * spn, size_t critLeaves, gsl_rng * rgsl, std::vector<cxsc::real>& eps) {};
				virtual bool priorityVisit(SPnode * spn, size_t critLeaves, std::vector<cxsc::real>& eps) {};
				//virtual cxsc::real getSPArea(SPnode * spn) {};

    };
    // end of SPnodeVisitor class

} // end namespace subpavings

#endif

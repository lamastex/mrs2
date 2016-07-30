/*
* Copyright (C) 2010, 2011, 2012 Jennifer Harlow
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
\brief MappedFobj1D definition and declaration.
*/

#ifndef __MAPPEDFOBJ2D__
#define __MAPPEDFOBJ2D__

#include "mappedFobj.hpp"
// Include cxsc types
#include "cxsc.hpp"

//using namespace cxsc;


/*! \brief An abstract class for target function objects on 2-dimensional real space

*/

namespace subpavings {

    class MappedFobj2D : public MappedFobj{

        public:
            //! a pure virtual function for interval image of a box
            virtual cxsc::interval operator()(const cxsc::ivector&) const;

            //! a pure virtual function for interval image of an interval pair
            virtual cxsc::interval operator()(const cxsc::interval&, const cxsc::interval&) const = 0;

            //! a pure virtual function for real image of an rvector
            virtual cxsc::real operator()(const cxsc::rvector&) const;

            //! a pure virtual function for real image of a real pair
            virtual cxsc::real operator()(const cxsc::real&, const cxsc::real&) const = 0;

            //! a virtual function for real image of midpoint of a box
            virtual cxsc::real imageMid(const cxsc::ivector&) const;

            //! a virtual function for real image of midpoint of an interval pair
            virtual cxsc::real imageMid(const cxsc::interval&, const cxsc::interval&) const;

            //! Destructor
            virtual ~MappedFobj2D(){}


    };

}
#endif

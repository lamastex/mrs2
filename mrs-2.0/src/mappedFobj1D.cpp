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
\brief MappedFobj1D definition.
*/


#include "mappedFobj1D.hpp"

//using namespace cxsc;

namespace subpavings {

    cxsc::interval MappedFobj1D::operator()(const cxsc::ivector& ivec) const
    {
        return this->operator()(ivec[1]);

    }

    cxsc::real MappedFobj1D::operator()(const cxsc::rvector& rv) const
    {
        return this->operator()(rv[1]);

    }

    cxsc::real MappedFobj1D::imageMid(const cxsc::ivector& ivec) const
    {

        return this->operator()(cxsc::mid(ivec[1]));

    }

    cxsc::real MappedFobj1D::imageMid(const cxsc::interval& ival) const
    {

        return this->operator()(cxsc::mid(ival));

        /*int ivL = Lb(iv);
        int ivU = Ub(iv);
        rvector rv(ivU - ivL +1 );
        for (int i = ivL; i<=ivU; i++) {
            rv[i] = mid(iv[i]);
        }
        return this->operator()(rv);
        */
    }



}

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
\brief MappedFobj definition.
*/


#include "mappedFobj.hpp"

//using namespace cxsc;

//#define MYDEBUG

namespace subpavings {

    cxsc::real MappedFobj::imageMid(const cxsc::ivector& ivec) const
    {

		int lb = Lb (ivec);
		int ub = Ub (ivec);
		rvector rv(ub-lb+1);
		for (int i = lb; i <= ub; ++i)
		{
			rv[i] = cxsc::mid(ivec[i]);
		}
		
		#ifdef MYDEBUG
			std::cout << "this midImage rvector is " << rv << std::endl;
		#endif
		
        return this->operator()(rv);

    }
	
	std::string MappedFobj::getName() const
	{
		return std::string("");
	}
}

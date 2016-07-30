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
\brief Declarations for multivariate Gaussian example function object class.
*/

#ifndef __GAUSSIANFOBJ_HPP__
#define __GAUSSIANFOBJ_HPP__


#include "mappedFobj.hpp"

namespace subpavings {

	class GaussianFobj : public subpavings::MappedFobj {

		public:
		
			//! constructor
			GaussianFobj();

			//! declare function for interval image of a box
			cxsc::interval operator()(const cxsc::ivector& ivec) const;
				
	  
			//! declare function for real image of reals
			cxsc::real operator()(const cxsc::rvector& r) const;
			
			
			// use base class for midImage
			
			//! get the function name.
			std::string getName() const;

			//! Destructor
			~GaussianFobj();

	};
}
#endif

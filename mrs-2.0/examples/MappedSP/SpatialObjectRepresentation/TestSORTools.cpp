/*
* Copyright (C) 2011, 2012 Jennifer Harlow
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
\brief Testing PiecewiseConstantFunction tools
 */

#include "TestSORTools.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;

std::vector< bool > makeRanges1()
{
	std::vector < bool > result;
	
	result.push_back(true);//X
	result.push_back(true);//XL
	result.push_back(false);//XR
	
	return result;
	
}

std::vector< bool > makeRanges2()
{
	std::vector < bool > result;
	
	result.push_back(true);//X
	result.push_back(true);//XL
	result.push_back(true);//XLL
	result.push_back(false); //XLR
	result.push_back(true);//XR
	result.push_back(false);//XRL
	result.push_back(true);//XRR
	
	return result;
	
}

std::vector< bool > makeRanges3()
{
	std::vector < bool > result;
	
	result.push_back(true);//X
	result.push_back(true);//XL
	result.push_back(true);//XLL
	result.push_back(false);//XLLL *
	result.push_back(true);//XLLR
	result.push_back(true);//XLLRL **
	result.push_back(false);//XLLRR *
	result.push_back(false);//XLR *
	result.push_back(true);//XR
	result.push_back(false);//XRL *
	result.push_back(true);//XRR
	result.push_back(true);//XRRL
	result.push_back(false);//XRRLL *
	result.push_back(true);//XRRLR **
	result.push_back(false);//XRRR *
	
	return result;
	
}

std::vector< bool > makeRanges4()
{
	std::vector < bool > result;
	
	result.push_back(true);//X
	result.push_back(true);//XL
	result.push_back(true);//XLL
	result.push_back(true);//XLLL *
	result.push_back(true);//XLLR *
	result.push_back(true);//XLR
	result.push_back(true);//XLRL *
	result.push_back(true);//XLRR *
	result.push_back(true);//XR
	result.push_back(true);//XRL 
	result.push_back(false);//XRLL *
	result.push_back(true);//XRLR *
	result.push_back(true);//XRR
	result.push_back(false);//XRRL *
	result.push_back(true);//XRRR *
	
	return result;
	
}

std::vector< bool > makeRanges5()
{
	std::vector < bool > result;
	
	result.push_back(true);//X
	result.push_back(true);//XL
	result.push_back(false);//XLL
	result.push_back(true);//XLR
	result.push_back(false);//XLRL *
	result.push_back(true);//XLRR
	result.push_back(true);//XLRRL **
	result.push_back(false);//XLRRR *
	result.push_back(true);//XR
	result.push_back(true);//XRL **
	result.push_back(true);//XRR
	result.push_back(true);//XRRL **
	result.push_back(false);//XRRR *
	
	return result;
	
}


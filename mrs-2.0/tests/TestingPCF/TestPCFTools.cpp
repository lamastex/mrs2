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

#include "TestPCFTools.hpp"


using namespace cxsc;
using namespace std;
using namespace subpavings;

std::vector< real > makeRanges1(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((8.0/8.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((2.0/8.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((6.0/8.0)*(2.0/rootvol)));//XR
	
	return result;
	
}

std::vector< real > makeRanges2(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((10.0/10.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((3.0/10.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/10.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((1.0/10.0)*(4.0/rootvol))); //XLR
	result.push_back(cxsc::real((7.0/10.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((2.0/10.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((5.0/10.0)*(4.0/rootvol)));//XRR
	
	return result;
	
}



std::vector< real > makeRanges3(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((8.0/8.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((2.0/8.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/8.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real(0.0));//XLLL
	result.push_back(cxsc::real((2.0/8.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real(0.0));//XLLRL
	result.push_back(cxsc::real((2.0/8.0)*(16.0/rootvol)));//XLLRR
	result.push_back(cxsc::real(0.0));//XLR
	result.push_back(cxsc::real((6.0/8.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/8.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((5.0/8.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((4.0/8.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((3.0/8.0)*(16.0/rootvol)));//XRRLL
	result.push_back(cxsc::real((1.0/8.0)*(16.0/rootvol)));//XRRLR
	result.push_back(cxsc::real((1.0/8.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}


bool checkSame(cxsc::real r1, const cxsc::real r2, int n) 
{
	
	return !(cxsc::abs(r1 - r2) 
			> n* DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(r1), cxsc::abs(r2))));
	
}


std::vector< real > makeRanges4(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((2.0/2.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real(0.0));//XLL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real(0.0));//XLRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((1.0/2.0)*(16.0/rootvol)));//XLRRL
	result.push_back(cxsc::real(0.0));//XLRRR
	result.push_back(cxsc::real((1.0/2.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real(0.0));//XRL
	result.push_back(cxsc::real((1.0/2.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real(0.0));//XRRL
	result.push_back(cxsc::real((1.0/2.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges5(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((10.0/10.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((6.0/10.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((2.0/10.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((4.0/10.0)*(4.0/rootvol)));//XLR
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XLRL
	result.push_back(cxsc::real((3.0/10.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((1.0/10.0)*(16.0/rootvol)));//XLRRL
	result.push_back(cxsc::real((2.0/10.0)*(16.0/rootvol)));//XLRRR
	result.push_back(cxsc::real((4.0/10.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.0/10.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((3.0/10.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((2.0/10.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((1.0/10.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRanges6(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((10.0/10.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((3.0/10.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((3.0/10.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((0.0/10.0)*(4.0/rootvol))); //XLR
	result.push_back(cxsc::real((7.0/10.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((-2.0/10.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((9.0/10.0)*(4.0/rootvol)));//XRR
	
	return result;
	
}


std::vector< real > makeRangesNegative(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((2.0/8.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((-3.0/8.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((5.0/8.0)*(2.0/rootvol)));//XR
	
	return result;
	
}

std::vector< real > makeRangesInfinite()
{
	std::vector < real > result;
	
	result.push_back(1.0);//X
	result.push_back(cxsc::Infinity);//XL
	result.push_back(cxsc::Infinity);//XR
	
	return result;
	
}


std::vector< real > makeRangesPositiveAndNegative()
{
	std::vector < real > result;
	
	result.push_back(-1.0);//X
	result.push_back(-2.0);//XL
	result.push_back(1.0);//XR
	
	return result;
	
}

std::vector< real > makeRangesPositiveAndNegativeInfinite()
{
	std::vector < real > result;
	
	result.push_back(-cxsc::Infinity);//X
	result.push_back(-cxsc::Infinity);//XL
	result.push_back(1.0);//XR
	
	return result;
	
}


std::vector< real > makeRangesPositiveAndZero()
{
	std::vector < real > result;
	
	result.push_back(1.0);//X
	result.push_back(1.0);//XL
	result.push_back(0.0);//XR
	
	return result;
	
}


std::vector< real > makeRangesInfiniteAndNegative()
{
	std::vector < real > result;
	
	result.push_back(-1.0);//X
	result.push_back(cxsc::Infinity);//XL
	result.push_back(-1.0);//XR
	
	return result;
	
}

std::vector< real > makeRangesInfiniteAndNegativeInfinite()
{
	std::vector < real > result;
	
	result.push_back(-1.0);//X
	result.push_back(cxsc::Infinity);//XL
	result.push_back(-cxsc::Infinity);//XR
	
	return result;
	
}

std::vector< real > makeRangesInfiniteAndZero()
{
	std::vector < real > result;
	
	result.push_back(cxsc::Infinity);//X
	result.push_back(cxsc::Infinity);//XL
	result.push_back(0.0);//XR
	
	return result;
	
}

std::vector< real > makeRangesNegativeAndZero()
{
	std::vector < real > result;
	
	result.push_back(-1.0);//X
	result.push_back(0.0);//XL
	result.push_back(-1.0);//XR
	
	return result;
	
}

RVecData& getData1(RVecData& data) 
{
    int d = 2;
	{ //XLLRR (-1,0), (-1,0)
		rvector thisrv(d);
		thisrv[1] = -0.5;
        thisrv[2] = -0.5;
		data.push_back(thisrv);
	}
	{ //XLLRR (-1,0), (-1,0)
		rvector thisrv(d);
		thisrv[1] = -0.75;
        thisrv[2] = -0.25;
		data.push_back(thisrv);
	}
	{ //XRL (0,2), (-2,0)
		rvector thisrv(d);
		thisrv[1] = 1.0;
        thisrv[2] = -1.0;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.5;
        thisrv[2] = 0.5;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.25;
        thisrv[2] = 0.75;
		data.push_back(thisrv);
	}
	{ //XRRLL (0,1), (0,1)
		rvector thisrv(d);
		thisrv[1] = 0.75;
        thisrv[2] = 0.25;
		data.push_back(thisrv);
	}
	{ //XRRLR (0,1), (1,2)
		rvector thisrv(d);
		thisrv[1] = 0.5;
        thisrv[2] = 1.5;
		data.push_back(thisrv);
	}
	{ //XRRR (1,2), (0,2)
		rvector thisrv(d);
		thisrv[1] = 1.5;
        thisrv[2] = 1.0;
		data.push_back(thisrv);
	}
	return data;
}


RVecData& getData2(RVecData& data) 
{
    int d = 2;
	{ //XL 
		rvector thisrv(d);
		thisrv[1] = -0.6;
        thisrv[2] = 0.4;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 1.3;
        thisrv[2] = 0.9;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getData2LeftOnly(RVecData& data) 
{
    int d = 2;
	{ //XL 
		rvector thisrv(d);
		thisrv[1] = -0.6;
        thisrv[2] = 0.4;
		data.push_back(thisrv);
	}
		
	return data;
}

RVecData& getData2RightOnly(RVecData& data) 
{
    int d = 2;
	
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 1.3;
        thisrv[2] = 0.9;
		data.push_back(thisrv);
	}
	
	return data;
}


RVecData& getData3(RVecData& data) 
{
    int d = 2;
	{ //XL (-2,0), (-2,2)
		rvector thisrv(d);
		thisrv[1] = -1.6;
        thisrv[2] = 0.2;
		data.push_back(thisrv);
	}
	{ //XRL (0,2), (-2,0)
		rvector thisrv(d);
		thisrv[1] = 1.7;
        thisrv[2] = -0.9;
		data.push_back(thisrv);
	}
	{ //XRR (0,1), (0,2)
		rvector thisrv(d);
		thisrv[1] = 0.8;
        thisrv[2] = 1.9;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getDataExtra1(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	
	return data;
}

RVecData& getDataExtra2(RVecData& data) 
{
    int d = 2;
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	{ //XR 
		rvector thisrv(d);
		thisrv[1] = 0;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -1;
        thisrv[2] = 1;
		data.push_back(thisrv);
	}
	{ //XLR 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = 2;
		data.push_back(thisrv);
	}
	{ //XLL 
		rvector thisrv(d);
		thisrv[1] = -2;
        thisrv[2] = -2;
		data.push_back(thisrv);
	}
	
	return data;
}


std::vector< real > makeRangesArithmeticSpecial1(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((18.0/36.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((6.0/36.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((5.0/36.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((2.0/36.0)*(8.0/rootvol)));//XLLL
	result.push_back(cxsc::real((3.0/36.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real((1.0/36.0)*(4.0/rootvol))); //XLR
	result.push_back(cxsc::real((12.0/36.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((5.5/36.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((2.5/36.0)*(8.0/rootvol)));//XRLL
	result.push_back(cxsc::real((3.0/36.0)*(8.0/rootvol)));//XRLR
	result.push_back(cxsc::real((6.5/36.0)*(4.0/rootvol)));//XRR
	result.push_back(cxsc::real((2.5/36.0)*(8.0/rootvol)));//XRRL
	result.push_back(cxsc::real((4.0/36.0)*(8.0/rootvol)));//XRRR
	
	return result;
	
}

std::vector< real > makeRangesArithmeticSpecial2(cxsc::real rootvol)
{
	std::vector < real > result;
	
	result.push_back(cxsc::real((12.0/24.0)*(1.0/rootvol)));//X
	result.push_back(cxsc::real((8.5/24.0)*(2.0/rootvol)));//XL
	result.push_back(cxsc::real((3.0/24.0)*(4.0/rootvol)));//XLL
	result.push_back(cxsc::real((1.0/24.0)*(8.0/rootvol)));//XLLL
	result.push_back(cxsc::real((2.0/24.0)*(8.0/rootvol)));//XLLR
	result.push_back(cxsc::real((0.75/24.0)*(16.0/rootvol)));//XLLRL
	result.push_back(cxsc::real((1.25/24.0)*(16.0/rootvol)));//XLLRR
	result.push_back(cxsc::real((5.5/24.0)*(4.0/rootvol))); //XLR
	result.push_back(cxsc::real((2.5/24.0)*(8.0/rootvol)));//XLRL
	result.push_back(cxsc::real((3.0/24.0)*(8.0/rootvol)));//XLRR
	result.push_back(cxsc::real((3.5/24.0)*(2.0/rootvol)));//XR
	result.push_back(cxsc::real((1.7/24.0)*(4.0/rootvol)));//XRL
	result.push_back(cxsc::real((1.1/24.0)*(8.0/rootvol)));//XRLL
	result.push_back(cxsc::real((0.6/24.0)*(16.0/rootvol)));//XRLLL
	result.push_back(cxsc::real((0.5/24.0)*(16.0/rootvol)));//XRLLR
	result.push_back(cxsc::real((0.6/24.0)*(8.0/rootvol)));//XRLR
	result.push_back(cxsc::real((0.4/24.0)*(16.0/rootvol)));//XRLRL
	result.push_back(cxsc::real((0.2/24.0)*(16.0/rootvol)));//XRLRR
	result.push_back(cxsc::real((1.8/24.0)*(4.0/rootvol)));//XRR
	
	return result;
	
}

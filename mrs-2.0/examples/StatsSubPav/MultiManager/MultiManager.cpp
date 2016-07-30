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


/*! \file MultiManager.cpp
\brief Testing MultiTreeManager class
*/



#include "spnodeall.hpp"
#include "spsnode.hpp"
#include "multitreemanager.hpp"

using namespace cxsc;
using namespace subpavings;
using namespace std;

int main()
{
    MultiTreeManager multiTM; // make a MultiTreeManager object

    int splits = 10;
    // make a graph of the outcome space after splits splits
    // ie, all full binary trees with splits+1 leaf nodes
    //multiTM.makeAndGraphOutcomeSpace(splits);

    multiTM.makeOutcomeSpace(splits);
    multiTM.mapPavings();

    return 0;

} // end of multimanager test

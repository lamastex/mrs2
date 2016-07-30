/*
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
*/

/*! \file
\brief declarations for testing MappedSP nodes*/


#ifndef __MAPPEDSPTESTING_HPP__
#define __MAPPEDSPTESTING_HPP__

#include "rvector.hpp"

#include <vector>
#include <stdexcept>
//forward class declaration

namespace subpavings {
    class RealMappedSPnode;
	
}

void testingArithmetic();

void testingInts();

void testingReals();

void testingRealsMarginalise();

void testingRealsMarginaliseFail();

void testingIntervals();

void testingRvectors();


void testingBools();
#endif

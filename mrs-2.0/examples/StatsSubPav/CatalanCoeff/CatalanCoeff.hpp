/*
 * Copyright (C) 2009 Raazesh Sainudiin and Gloria Teng
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

/*! \file CatalanCoeff.hpp
  \brief Declarations for Computing Catalan Coefficients and their Frequencies.
*/

//File: CatalanCoeff.hpp

#ifndef ___CC_SUBPAVTEST_HPP__
#define __CC_SUBPAVTEST_HPP__

//#include "HistogramWrappers.hpp"
#include <time.h>
#include <fstream>
#include <fstream>
#include <sstream>  // to be able to manipulate strings as streams
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <valarray>
#include <map>
#include <math.h>
#include <algorithm>

//using namespace cxsc;
using namespace std;
#endif

//rename vector<int> as State
//make this into int or short int
typedef vector<int> State;

//rename vector< vector<int> > as States
typedef vector< vector<int> > States;

/*! templatized function object for lexicographical sorting of vec
tors whose elements have total ordering
*/
template <class T>
class LexicoSorting
{
public:
bool operator() (const T& t1, const T& t2) const
{
return std::lexicographical_compare(&t1[0], &t1[t1.size()-1], &t2[0], &t2[t2.size()-1]);
}
};

/*! templatized function object for converting to strings
*/
template <class LDInt>
inline std::string to_string (const LDInt& s)
{
  std::stringstream ss;
  ss << s;
  string finalString = ss.str() + '\t';
  return finalString;
}

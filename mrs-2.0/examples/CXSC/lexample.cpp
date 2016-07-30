/* 
 * Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
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

/*! \file lexample.cpp 
\brief C-XSC Example on interval staggered arithmetic
*/
#include "l_interval.hpp"  // interval staggered arithmetic in C-XSC
#include <iostream>
using namespace cxsc;
using namespace std;

int main() 
{
  l_interval a, b;         // Multiple-precision intervals in C-XSC
  stagprec = 2;            // global integer variable      
  cout << SetDotPrecision(16*stagprec, 16*stagprec-3) << RndNext;
  // I/O for variables of type l_interval is done using
  // the long accumulator (i.e. a dotprecision variable)   

  a = 1.0;                  // a   = [1.0,1.0]       
  "[1, 2]" >> b;          // string to interval conversion b   = [1.0,2.0]        
  cout << "a - a = " << a-a << endl;
  cout << "b - b = " << b-b << endl;
  cout << "a/b = " << a/b << endl;  
}

/* --------------------------- Output ------------------------------
$ ./lexample 
a - a = [ 0.00000000000000000000000000000, 0.00000000000000000000000000000]
b - b = [-1.00000000000000000000000000000, 1.00000000000000000000000000000]
a/b = [ 0.50000000000000000000000000000, 1.00000000000000000000000000000]
------------------------------------------------------------------*/

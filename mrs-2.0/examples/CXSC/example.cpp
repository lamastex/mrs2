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

/*! \file example.cpp 
\brief C-XSC Example on interval subtraction
*/

#include "interval.hpp"  // include interval arithmetic in C-XSC
#include <iostream>      // include standard Input Output STREAM 
using namespace cxsc;    
using namespace std;

int main()
{
  interval a, b;            // Standard intervals     
  a = 1.0;                  // a   = [1.0,1.0]       
  "[1, 2]" >> b;          // string to interval conversion b   = [1.0,2.0]        
  cout << "a - a = " << a-a << endl;
  cout << "b - b = " << b-b << endl;
}

/* --------------------------- Output ------------------------------
$ ./example 
a - a = [ -0.000000,  0.000000]
b - b = [ -1.000000,  1.000000]
------------------------------------------------------------------*/

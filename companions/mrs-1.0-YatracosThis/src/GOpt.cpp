/*
 * Copyright (C) 2005, 2006, 2007, 2008 Raazesh Sainudiin
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

/*! \file      GOpt.cpp
\brief Global optimisation definitions

A set of processes for performing global optimisation on a function using
AllGOp from the C++ Toolbox for Verified Computing
*/

// global optimisation using C++ Toolbox for Verified Computing
#include "GOpt.hpp"
#include <iostream>
#include <fstream>

#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gop.hpp>  // cxsc global optimisation
// increase stack size for some special C++ compilers
#include <stacksz.hpp>


/*! \brief declare a global pointer to an Fobj

This is horrible but it is the only way around the various problems in
 implementing global optimisation for now.
*/
Fobj* fToOtp = NULL;

/*! \brief declare a global pointer to an int for target function label

This is horrible but it is the only way around the various problems in
implementing global optimisation for now.
*/
int flabel = 0;

// implementation of GOptMin
/* f is a pointer to the Fobj we will do global optimisation on
 to look for global minimums.  The label defaults to 0.
*/
void GOptMin(Fobj* f, ivector search, real t, const int label)
{
  fToOtp = f;       // point the global pointer to f
  flabel = label;   // and set flabel = label
  ivector SearchInterval = search;
  real Tolerance = t;
  // check the dimensions match
  int search_dim = Ub(search) - Lb(search) +1;
  // get the dimensions of the domain for boxes labelled 0;
  size_t func_dim = f->getLabeledDomainDim(0);
  if(static_cast<int>(func_dim) != search_dim)
  {
    std::cerr << "Error in GOpt: dimensions of search box "
      << search_dim << " do not match dimensions of function domain "
      << func_dim << endl<< std::endl;
    exit(1);
  }
  // these are filled in by running AllGOp
  interval    Minimum;
  imatrix     Opti;
  intvector   Unique;
  int         NumberOfOptis;
  int   Error;
  // running AllGop fills in Opti, Unique, NumberOfOptis, Mimimum and Error,
  // all of which are passed by reference
  AllGOp(funcHessMin, SearchInterval, Tolerance,
    Opti, Unique, NumberOfOptis, Minimum, Error);
  // AllGOp also takes a final optional parameter for maximum number of
  // optimisation attemtps to perform
  // this has default value of MaxCount = 10000 - see gop.hpp
  // print the outcomes
  printOutcomeMin(SearchInterval, Tolerance, Opti, Unique, NumberOfOptis,
                  Minimum, Error);
}

// implementation of GOptMax
/* f is a pointer to the Fobj we will do global optimisation on
 to look for global maximums.  The label defaults to 0.
*/
void GOptMax(Fobj* f, ivector search, real t, const int label)
{
  fToOtp = f;       // point the global pointer to f
  flabel = label;   // and set flabel = label
  ivector SearchInterval = search;
  real Tolerance = t;
  // check the dimensions match
  int search_dim = Ub(search) - Lb(search) +1;
                    // get the dimensions of the domain for boxes labelled 0;
  size_t func_dim = f->getLabeledDomainDim(0);
  if(static_cast<int>(func_dim) != search_dim)
  {
    std::cerr << "Error in GOpt: dimensions of search box " << search_dim
      << " do not match dimensions of function domain " << func_dim << std::endl
      << std::endl;
    exit(1);
  }
  // these are filled in by running AllGOp
  interval    Maximum;
  imatrix     Opti;
  intvector   Unique;
  int         NumberOfOptis;
  int   Error;
  // running AllGop fills in Opti, Unique, NumberOfOptis, Mimimum and Error,
  // all of which are passed by reference
  AllGOp(funcHessMax, SearchInterval, Tolerance,
         Opti, Unique, NumberOfOptis, Maximum, Error);
  // AllGOp also takes a final optional parameter for maximum number of
  // optimisation attemtps to perform
  // this has default value of MaxCount = 10000 - see gop.hpp
  // print the outcomes
  printOutcomeMax(SearchInterval, Tolerance, Opti, Unique, NumberOfOptis,
                  Maximum, Error);
}

/* Function conforming to typedef HTscalar_FctPtr.
 Acts as 'wrapper' and uses the HessType operator()(const HTvector&) for the
 type of Fobj in FtoOpt global
*/
HessType funcHessMin(const HTvector& x)
{
  HessType hh;
  // use the operator () of the function pointed to by global fToOpt to
  // calculate hh
  hh = (*fToOtp)(x, flabel);
  return hh;
}

/* Function conforming to typedef HTscalar_FctPtr.
 Acts as 'wrapper' and uses the HessType operator()(const HTvector&) for the
 type of Fobj in FtoOpt global but uses the negative of the () operator.
*/
HessType funcHessMax(const HTvector& x)
{
  HessType hh;
  // use negative of the operator () of the function pointed to by global
  // fToOpt to calculate hh
  hh = -((*fToOtp)(x, flabel));
  return hh;
}

// Tell me all about the results from global minimums
void printOutcomeMin(ivector& SearchInterval, real& Tolerance, imatrix& Opti,
                     intvector& Unique, int NumberOfOptis, interval& Minimum,
                     int Error)
{
  // Output format
  std::cout << SetPrecision(23,15) << Scientific;

  std::cout << std::endl << std::endl << std:: endl;
  std::cout << "Global optimisation: minimums " << std::endl;
  std::cout << "Results of running AllGOp with Tolerance " << Tolerance
    << std::endl;
  std::cout << "with an initial search box of: " << std::endl;
  std::cout << SearchInterval << endl<<endl;

  std::cout << "The results for global optimisation (minimums) are "
    << std::endl;

  for (int i = 1; i <= NumberOfOptis; i++)
  {
    std::cout << Opti[i];
    if (Unique[i])
      std::cout
        << " encloses a locally unique candidate for a global minimiser!";
    else
      std::cout << " may contain a local or global minimiser!";
    std::cout << std::endl << std::endl;
  }

  if (NumberOfOptis != 0)
  {
    std::cout << Minimum << std::endl
      << "encloses the global minimum value!" << std::endl << std::endl;
  }

  std::cout << NumberOfOptis << " interval enclosure(s)" << std::endl;

  if (Error)
    std::cout << endl << AllGOpErrMsg(Error) << std::endl;
  else if (NumberOfOptis == 1 && Unique[1])
    std::cout << std::endl << "We have validated that there is "
        "a unique global optimiser!" << std::endl;

  std::cout << std::endl << "End of global optimisation for minimums"
    << std::endl << std::endl;
}

// Tell me all about the results form global maximums
void printOutcomeMax(ivector& SearchInterval, real& Tolerance, imatrix& Opti,
                     intvector& Unique, int NumberOfOptis, interval& Maximum,
                     int Error)
{
                    // Output format
  std::cout << SetPrecision(23,15) << Scientific;

  std::cout << std::endl << std::endl << std:: endl;
  std::cout << "Global optimisation: Maximums" << std::endl;
  std::cout << "Results of running AllGOp with Tolerance " << Tolerance
    << std::endl;
  std::cout << "with an initial search box of: " << std::endl;
  std::cout << SearchInterval << endl<<endl;

  std::cout << "The results for global optimisation (maximums) are "
    << std::endl;

  for (int i = 1; i <= NumberOfOptis; i++)
  {
    std::cout << Opti[i];
    if (Unique[i])
      std::cout
        << " encloses a locally unique candidate for a global maximiser!";
    else
      std::cout << " may contain a local or global maximiser!";
    std::cout << std::endl << std::endl;
  }

  if (NumberOfOptis != 0)
  {
    std::cout << Maximum << std::endl
      << "encloses the global maximum value!" << std::endl << std::endl;
  }

  std::cout << NumberOfOptis << " interval enclosure(s)" << std::endl;

  if (Error)
    std::cout << endl << AllGOpErrMsg(Error) << std::endl;
  else if (NumberOfOptis == 1 && Unique[1])
    std::cout << std::endl << "We have validated that there is "
        "a unique global optimiser!" << std::endl;

  std::cout << std::endl << "End of global optimisation for maximums"
    << std::endl << std::endl;

}

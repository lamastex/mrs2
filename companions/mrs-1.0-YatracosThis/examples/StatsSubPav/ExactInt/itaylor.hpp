/*
**  Copyright (C) 1999-2006 F. Blomquist, M. Braeuer, M. Grimmer,
**                          W. Hofschuster, W. Kraemer
**                          Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

////////////////////////////////////////////////////////////////////////
//
//     Headerfile itaylor.hpp for onedimensional Taylor-arithmetic
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//     Updated by F. Blomquist, M. Grimmer
//     Extended version 05.03.2006 by M. Grimmer
////////////////////////////////////////////////////////////////////////

/*----------------------------------------------------------------------

Definition of the class itaylor:

class itaylor for calculating all Taylor coefficients up to the maximal order p.

  Elements:
            int p .................... maximal order of the Taylor coefficients
            ivector tayl ............. interval vector; Lb=0, Ub=p;

	    static ivector faks ...... storing n! with n<=170
            static int initialized ... switcher for initialization of faks.
            static void initialize().. performs initialization of faks.

  Element functions, methodes:

            see implementation

Implementation file: itaylor.cpp

-----------------------------------------------------------------------*/

#ifndef _ITAYLOR_H
#define _ITAYLOR_H

//cxsc headers
#include <imath.hpp>
#include <interval.hpp>
#include <ivector.hpp>
#include <idot.hpp>

//C++ standard headers
#include <iostream>

using namespace cxsc;
using namespace std;

namespace taylor{


enum{
    _i_ln,

    _i_tan,
    _i_cot,

    _i_asin,
    _i_acos,
    _i_atan,
    _i_acot,

    _i_tanh,
    _i_coth,

    _i_asinh,
    _i_acosh,
    _i_atanh,
    _i_acoth

};

///////////////////////////////////////////////////////////////
//
// Class itaylor
//
///////////////////////////////////////////////////////////////

// class itaylor for calculating all Taylor-coefficients up to order p
// of functions with one independent variable.

class itaylor{

 private:
  int p;        // max. Taylor-order of the object;
  ivector tayl; // Interval vector with Taylor coefficients

  //Removed by Tomas Johnson
  //  static ivector faks;
  static int initialized;
  static void initialize();

 public:
  // Constructors and Destructors:
  itaylor();
  itaylor(const itaylor& );
  explicit itaylor(int order);
  itaylor(int order, const real& value);  // (x,1,0,...,0) Conversion error!!
  itaylor(int order, const interval& value); // (x,1,0,...,0)
  ~itaylor(){;};

  // Initialization functions for independent variables (x,1,0,...,0):
  // Caution: (x,1,0,...,0)  for real x conversion errors are possible!
  friend itaylor var_itaylor(int ord, const real& c);
  friend itaylor var_itaylor(int ord, const interval& c);

  // Initialization functions for constants (c,0,0,...,0):
  // Caution: (c,0,0,...,0)  for real c conversion errors are possible!
  friend itaylor const_itaylor(int ord, const real& c);
  friend itaylor const_itaylor(int ord, const interval& c);

  // assignment operators
  itaylor operator=(const itaylor& );
  itaylor operator=(int);
  itaylor operator=(const real& );
  itaylor operator=(const interval& );
  itaylor operator=(const ivector&); //added, mg2005-08
                                //const since C-XSC 2.1, mg2006-02
  // relational operators
  int operator==(itaylor&);     //added, mg2005-08
  int operator!=(itaylor&);     //added, mg2005-08
  int operator<=(itaylor&);     //added, mg2005-08
  int operator<(itaylor&);      //added, mg2005-08
  int operator>=(itaylor&);     //added, mg2005-08
  int operator>(itaylor&);      //added, mg2005-08
  
  // component access
  interval& operator[](int n);      //added, mg2005-08
  int order() {return p;}           //added, mg2005-08

  // additional component access functions (for compatibility reasons) //mg2005-08
  friend int get_order(const itaylor& x);
  friend ivector get_all_coef(const itaylor& x);
  friend interval get_j_coef(const itaylor& x, int j);

  //Removed by Tomas Johnson
  // access to the derivative of order j:
  // friend interval get_j_derive(const itaylor& x, int j); 
 
  // Output:
  friend void print_itaylor(const itaylor& x);
                                // kept for compatibility reasons

  friend void print_itaylor(std::ostream& os, const itaylor& x, int width=0, int digits=0);
                                // added, mg2005,2006

  friend std::ostream& operator<< (std::ostream&, itaylor&);
                                // added, mg2005-08

  // Overloading the operators for elements of the class itaylor:

  // operator - :
  friend itaylor operator-(const itaylor& x);

  // operators +,-,*,/  for (itaylor, itaylor):
  friend itaylor operator-(const itaylor& x, const itaylor& y);
  friend itaylor operator+(const itaylor& x, const itaylor& y);
  friend itaylor operator*(const itaylor& x, const itaylor& y);
  friend itaylor operator/(const itaylor& x, const itaylor& y);

  // operators +,-,*,/ for (interval, itaylor):
  friend itaylor operator-(const interval& x, const itaylor& y);
  friend itaylor operator+(const interval& x, const itaylor& y);
  friend itaylor operator*(const interval& x, const itaylor& y);
  friend itaylor operator/(const interval& x, const itaylor& y);

  // operators +,-,*,/ for (itaylor, interval):
  friend itaylor operator-(const itaylor& x, const interval& y);
  friend itaylor operator+(const itaylor& x, const interval& y);
  friend itaylor operator*(const itaylor& x, const interval& y);
  friend itaylor operator/(const itaylor& x, const interval& y);

  // operators +,-,*,/ for (real, itaylor):
  friend itaylor operator-(const real& x, const itaylor& y);
  friend itaylor operator+(const real& x, const itaylor& y);
  friend itaylor operator*(const real& x, const itaylor& y);
  friend itaylor operator/(const real& x, const itaylor& y);

  // operators +,-,*,/ for (itaylor, real):
  friend itaylor operator-(const itaylor& x, const real& y);
  friend itaylor operator+(const itaylor& x, const real& y);
  friend itaylor operator*(const itaylor& x, const real& y);
  friend itaylor operator/(const itaylor& x, const real& y);

  // operators +,-,*,/ for (int, itaylor):
  friend itaylor operator-(int x, const itaylor& y);
  friend itaylor operator+(int x, const itaylor& y);
  friend itaylor operator*(int x, const itaylor& y);
  friend itaylor operator/(int x, const itaylor& y);

  // operators +,-,*,/ for (itaylor, int):
  friend itaylor operator-(const itaylor& x, int y);
  friend itaylor operator+(const itaylor& x, int y);
  friend itaylor operator*(const itaylor& x, int y);
  friend itaylor operator/(const itaylor& x, int y);

  // Overloading the standard functions:
  friend itaylor sqrt(const itaylor& x);
  friend itaylor sqrt(const itaylor& x, int n);
  friend itaylor sqrt1px2(const itaylor& x);
  friend itaylor sqrtp1m1(const itaylor& x);
  friend itaylor sqrt1mx2(const itaylor& x);
  friend itaylor sqrtx2m1(const itaylor& x);
  friend itaylor sqr(const itaylor& x);
  friend itaylor pow(const itaylor& x, const interval& alpha);

  friend itaylor exp(const itaylor& x);
  friend itaylor expm1(const itaylor& x);

  friend itaylor ln(const itaylor& x);
  friend itaylor lnp1(const itaylor& x);

  friend itaylor sin(const itaylor& x);
  friend itaylor cos(const itaylor& x);
  friend itaylor tan(const itaylor& x);
  friend itaylor cot(const itaylor& x);

  friend itaylor sinh(const itaylor& x);
  friend itaylor cosh(const itaylor& x);
  friend itaylor tanh(const itaylor& x);
  friend itaylor coth(const itaylor& x);

  friend itaylor asin(const itaylor& x);
  friend itaylor acos(const itaylor& x);
  friend itaylor atan(const itaylor& x);
  friend itaylor acot(const itaylor& x);

  friend itaylor asinh(const itaylor& x);
  friend itaylor acosh(const itaylor& x);
  friend itaylor atanh(const itaylor& x);
  friend itaylor acoth(const itaylor& x);

  friend itaylor erf(const itaylor& x);   // added, mg2006-03
  friend itaylor erfc(const itaylor& x);  // added, mg2006-03

  // Help function
  friend void f_g_u(const itaylor& f, const itaylor& g, const itaylor& u,
                    int nb_function); 

  //ALTERED by Tomas Johnson
  friend itaylor powerAtZero(const itaylor &x, int n);
};
  //Declaration of friend functions
  itaylor var_itaylor(int ord, const interval& c);
  itaylor const_itaylor(int ord, const interval& c);

} // End of namespace taylor

#endif

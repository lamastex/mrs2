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
//     Updated by F. Blomquist, M. Grimmer
//     Extended version 05.03.2006 by M. Grimmer
////////////////////////////////////////////////////////////////////////


#ifndef _DIM2TAYLOR_H
#define _DIM2TAYLOR_H

#include <iostream>

#include "ivector.hpp"
#include "imath.hpp"
#include "idot.hpp"

using namespace cxsc;
namespace taylor{

enum{_ln,_lnp1,_tan,_cot,_asin,_acos,_atan,_acot,_tanh,_coth,
       _asinh,_acosh,_atanh,_acoth,_sqrtp1m1};

/*--------------------------------------------------------------------------

class dim2taylor: Taylor arithmetic for functions of two variables

        int p.....maximum taylor order
                   
 ivector* dat.....dynamic vector of interval vectors implemented as 
                  triangle matrix 

--------------------------------------------------------------------------*/



class  dim2taylor{  
 private:
  
  int p; // Order of Taylor expansion
  ivector* dat; // pointer to a block (array) of elements of type
		// ivector, the block is realized as a triangle matrix.
  
 public:

  dim2taylor();     // Default constructor
  dim2taylor(int);  // Constructor for special order (int)
  dim2taylor(const dim2taylor& ); // Copy constructor

  //NEW, by T.Johnson
  dim2taylor(const dim2taylor&, int); //increase order by 1

  ~dim2taylor();

  dim2taylor& operator=(const dim2taylor& );
  ivector& operator[](int n) const; 

  friend dim2taylor init_var(int, int, const interval& );
  friend dim2taylor init_const(int, const interval& );

  int get_p() const {return p;};
  int order() const {return p;}  //added, mg2005-08/2005-11

  void print_dim2taylor(); // debug
  void print_dim2taylor(std::ostream&); //added, mg2005-11

  friend std::ostream& operator<< (std::ostream&, dim2taylor&);
                                        // added, mg2005-08/2005-11

// Overloading the arithmetic operators:

  friend dim2taylor operator-(const dim2taylor& s);
  friend dim2taylor operator-(const dim2taylor&, const dim2taylor& );
  friend dim2taylor operator+(const dim2taylor&, const dim2taylor& );
  friend dim2taylor operator*(const dim2taylor&, const dim2taylor& );
  friend dim2taylor operator/(const dim2taylor&, const dim2taylor& );

  friend dim2taylor operator-(const interval&, const dim2taylor& );
  friend dim2taylor operator+(const interval&, const dim2taylor& );
  friend dim2taylor operator*(const interval&, const dim2taylor& );
  friend dim2taylor operator/(const interval&, const dim2taylor& );

  friend dim2taylor operator-(const dim2taylor&, const interval& );
  friend dim2taylor operator+(const dim2taylor&, const interval& );
  friend dim2taylor operator*(const dim2taylor&, const interval& );
  friend dim2taylor operator/(const dim2taylor&, const interval& );

  // Caution: possible conversion errors
  friend dim2taylor operator-(const real&, const dim2taylor& );
  friend dim2taylor operator+(const real&, const dim2taylor& );
  friend dim2taylor operator*(const real&, const dim2taylor& );
  friend dim2taylor operator/(const real&, const dim2taylor& );

  // Caution: possible conversion errors
  friend dim2taylor operator-(const dim2taylor&, const real& );
  friend dim2taylor operator+(const dim2taylor&, const real& );
  friend dim2taylor operator*(const dim2taylor&, const real& );
  friend dim2taylor operator/(const dim2taylor&, const real& );

  friend dim2taylor operator-(int, const dim2taylor& );
  friend dim2taylor operator+(int, const dim2taylor& );
  friend dim2taylor operator*(int, const dim2taylor& );
  friend dim2taylor operator/(int, const dim2taylor& );

  friend dim2taylor operator-(const dim2taylor&, int);
  friend dim2taylor operator+(const dim2taylor&, int);
  friend dim2taylor operator*(const dim2taylor&, int);
  friend dim2taylor operator/(const dim2taylor&, int);

  friend dim2taylor sqr(const dim2taylor& );
  friend dim2taylor sqrt(const dim2taylor& );
  friend dim2taylor sqrt1px2(const dim2taylor& );
  friend dim2taylor sqrtx2m1(const dim2taylor& );
  friend dim2taylor sqrtp1m1(const dim2taylor& );
  friend dim2taylor pow(const dim2taylor& , const interval& );
  friend dim2taylor power(const dim2taylor& , int );

  friend dim2taylor exp(const dim2taylor& );
  friend dim2taylor ln(const dim2taylor& );
  friend dim2taylor lnp1(const dim2taylor& );

  friend dim2taylor sin(const dim2taylor& );
  friend dim2taylor cos(const dim2taylor& );
  friend dim2taylor tan(const dim2taylor& );
  friend dim2taylor cot(const dim2taylor& );

  friend dim2taylor sinh(const dim2taylor& );
  friend dim2taylor cosh(const dim2taylor& );
  friend dim2taylor tanh(const dim2taylor& );
  friend dim2taylor coth(const dim2taylor& );

  friend dim2taylor asin(const dim2taylor& );
  friend dim2taylor acos(const dim2taylor& );
  friend dim2taylor atan(const dim2taylor& );
  friend dim2taylor acot(const dim2taylor& );

  friend dim2taylor asinh(const dim2taylor& );
  friend dim2taylor acosh(const dim2taylor& );
  friend dim2taylor atanh(const dim2taylor& );
  friend dim2taylor acoth(const dim2taylor& );

  friend dim2taylor erf(const dim2taylor& );    // added, mg2006-03
  friend dim2taylor erfc(const dim2taylor& );   // added, mg2006-03

  friend void f_g_u(const dim2taylor& , const dim2taylor& , const dim2taylor& ,
			  int ); //was "int&", mg2005


  //////////////////////////////////////////////////////
  //Alterations by Tomas Johnson:
  //////////////////////////////////////////////////////

  inline void rescale(int i, int j, interval factor) {dat[i][j]*=factor;}


};

/*--------------------------------------------------------------------------

Class dim2taylor_vector: For definition and initialization of objects
                         of the class dim2taylor respective to the two
                         independent variables x and y.

--------------------------------------------------------------------------*/
class dim2taylor_vector{  // Vector for independent variables
 private:
  
  int dim; // Dimension of the vector
  int lb;
  int ub;

  int p_el;
  dim2taylor* comp; // Pointer to a dyn. array of elements of type dim2taylor
  
 public:
    // Constructors:
    dim2taylor_vector(); // Default constructor
    dim2taylor_vector(int, int, int);  // order, lb, ub;
    dim2taylor_vector(const dim2taylor_vector& );
  
    ~dim2taylor_vector();
    // Member functions
    int get_dim()  const {return dim;};   // Number of elements 
                                          // of type dim2taylor
    int get_lb ()  const {return lb;};    // lb index of the array
    int get_ub ()  const {return ub;};    // ub index of the array
    int get_p_el() const {return p_el;};  // p_el: Order of the
                                          // Taylor expansion; 
    friend int Lb(const dim2taylor_vector&);  //added, mg2005
    friend int Ub(const dim2taylor_vector&);  //added, mg2005
    // Operators
    dim2taylor_vector& operator=(const dim2taylor_vector& );
    dim2taylor& operator[](int n) const;  

    // Function for initialization of objects of type dim2taylor_vector:
    friend dim2taylor_vector init_var(int order, ivector& values);

};
  //Declaration of friend functions

  dim2taylor init_var(int, int, const interval& );
  dim2taylor init_const(int, const interval& );
  dim2taylor_vector init_var(int order, ivector& values);
} // end of namespace taylor

#endif

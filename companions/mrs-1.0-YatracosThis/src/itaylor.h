/*   File: itaylor.h

     The header file for an interval  Taylor series class. The recursive
     formulas were obtained from the books

     Aberth, O. "Precise Numerical Analysis" WCB Publishers, 1988.
     Aberth, O. "Precise Numerical Methods Using C++" Academic Press, 1998

     Author: Warwick Tucker <warwick@math.uu.se>
     Latest edit: Wed Jun 16 08:50:39 CEST 2004
*/

#ifndef __ITAYLOR_H__
#define __ITAYLOR_H__

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "real.hpp"
#include "rmath.hpp"
#include "imath.hpp"
#include "interval.hpp"

using namespace std;
using namespace cxsc;

class itaylor {
 private:
  int        order;
  interval  *coeff;
 public:
  itaylor ();
  itaylor (const itaylor  &a);
  itaylor (const interval &constant);
  itaylor (const double   &constant);

  virtual ~itaylor()    { delete coeff; }
  itaylor  &operator =  (const itaylor &a);
  interval &operator [] (int index) const;
  int      getOrder     (void) const { return order; }

  friend itaylor operator + (const itaylor  &a);
  friend itaylor operator - (const itaylor  &a);

  friend itaylor operator + (const itaylor  &a, const itaylor &b);
  friend itaylor operator - (const itaylor  &a, const itaylor &b);
  friend itaylor operator * (const itaylor  &a, const itaylor &b);
  friend itaylor operator / (const itaylor  &a, const itaylor &b);

  friend itaylor exp        (const itaylor  &a);
  friend itaylor log        (const itaylor  &a);
  friend itaylor sin        (const itaylor  &a);
  friend itaylor cos        (const itaylor  &a);
  friend itaylor tan        (const itaylor  &a);
  friend itaylor pow        (const itaylor  &a, const itaylor  &b);
  friend itaylor pow        (const itaylor  &a, const interval &t);
  friend itaylor pow        (const itaylor  &a, const int      &n);

  friend itaylor derivative (const itaylor  &a);
  friend itaylor integral   (const itaylor  &a);
  friend itaylor derivative (const itaylor  &a, int k);
  friend itaylor integral   (const itaylor  &a, int k);

  static itaylor variable   (const interval &x, int degree);
  static itaylor constant   (const interval &x, int degree);

  friend void    resize     (      itaylor  &a, int degree);
  friend void    clear      (      itaylor  &a);

  // Can the compiler distinguish between these two?
  friend itaylor round      (const itaylor  &a, int degree);
  friend void    round      (      itaylor  &a, int degree);

  friend int     orderOf    (const itaylor  &a);
  friend ostream &operator << (ostream &oS, const itaylor &a);
};

#endif // __ITAYLOR_H__

/*   File: ia_ad.h

     The header for a C++ class for first-order differentiation arithmetic, 
     including some standard functions. The underlying number field is the
     IEEE floating point intervalws.

     Author: Warwick Tucker <warwick@math.uu.se>

     Compilation: None.
     Latest edit: Mon Jun 14 15:50:21 CEST 2004
*/

#ifndef IA_AD_H
#define IA_AD_H

#include <cmath>
#include <iostream>
#include "intervalw.h"
using namespace std;

class ia_ad
{
  friend ia_ad operator + (const ia_ad &x, const ia_ad &y);
  friend ia_ad operator - (const ia_ad &x, const ia_ad &y);
  friend ia_ad operator * (const ia_ad &x, const ia_ad &y);
  friend ia_ad operator / (const ia_ad &x, const ia_ad &y);

  friend ia_ad exp (const ia_ad &x);
  friend ia_ad log (const ia_ad &x);
  friend ia_ad pow (const ia_ad &x, double a);
  friend ia_ad sin (const ia_ad &x);
  friend ia_ad cos (const ia_ad &x);
  friend ia_ad tan (const ia_ad &x);
  friend ia_ad abs (const ia_ad &x);

  friend intervalw value(const ia_ad &x) { return x.val; }
  friend intervalw deriv(const ia_ad &x) { return x.der; }

  friend ostream & operator << (ostream &, const ia_ad &);
  friend istream & operator >> (istream &, ia_ad &);
public:
  ia_ad() {}
  ia_ad(const ia_ad &x)            {val = x.val; der = x.der;}
  ia_ad(int fx)                    {val = intervalw(fx); der = intervalw(0.0);}
  ia_ad(double fx)                 {val = intervalw(fx); der = intervalw(0.0);}
  ia_ad(intervalw fx)               {val = fx;    der = intervalw(0.0);}
  ia_ad(intervalw fx, intervalw dfx) {val = fx;    der = dfx;}
  static ia_ad variable(const intervalw &x);
  static ia_ad constant(const intervalw &x);
private:
  intervalw val;
  intervalw der;
};

#endif // IA_AD_H

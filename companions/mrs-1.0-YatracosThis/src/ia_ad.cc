/*   File: ia_ad.cc

     A C++ class for first-order differentiation arithmetic, including some 
     standard functions. The underlying number field is IEEE floating point 
     intervalws.

     Author: Warwick Tucker <warwick@math.uu.se>

     Compilation: g++ -Wall -c ia_ad.cc
     Latest edit: Mon Jun 14 14:03:01 CEST 2004
*/

#include "ia_ad.h"

static double sign(intervalw x) 
{ 
  if ( sup(x) < 0.0 ) 
    return -1.0;
  else if ( 0.0 < inf(x) )
    return +1.0;
  else
    return 0.0;
}

ia_ad operator + (const ia_ad &x, const ia_ad &y)
{ return ia_ad(x.val + y.val, x.der + y.der);}

ia_ad operator - (const ia_ad &x, const ia_ad &y)
{ return ia_ad(x.val - y.val, x.der - y.der);}

ia_ad operator * (const ia_ad &x, const ia_ad &y)
{ return ia_ad(x.val*y.val, x.val*y.der + x.der*y.val); }

ia_ad operator / (const ia_ad &x, const ia_ad &y)
{ 
  if ( subset(0.0, y.val) ) {
    cerr << "Error: division by zero. Bye!" << endl;
    exit(1);
  }
  return ia_ad(x.val/y.val, (x.der - (x.val/y.val)*y.der)/y.val);
}

ia_ad ia_ad::variable(const intervalw &x)
{ return ia_ad(x, 1.0); }

ia_ad ia_ad::constant(const intervalw &x)
{ return ia_ad(x, 0.0); }

ia_ad exp (const ia_ad &x)
{
  intervalw fx = exp(x.val);
  return ia_ad(fx, x.der*fx);
}

ia_ad log (const ia_ad &x)
{
  if ( inf(x.val) < 0.0 ) {
    cerr << "Error: log of negative number. Bye!" << endl;
    exit(1);
  }
  return ia_ad(log(x.val), x.der/x.val);
}

ia_ad pow (const ia_ad &x, double a)
{
  if ( subset(0.0, x.val) && (a < 0) ) {
    cerr << "Error: negative power of zero. Bye!" << endl;
    exit(1);
  }
  return ia_ad(pow(x.val, a), x.der*a*pow(x.val, a - 1));
}

ia_ad sin (const ia_ad &x)
{ return ia_ad(sin(x.val), x.der*cos(x.val)); }

ia_ad cos (const ia_ad &x)
{ return ia_ad(cos(x.val), 0.0 - x.der*sin(x.val)); }

ia_ad tan (const ia_ad &x)
{ return sin(x)/cos(x); }

ia_ad abs (const ia_ad &x)
{ 
  if ( subset(0.0, x.val) ) {
    cerr << "Error: absolute value of zero. Bye!" << endl;
    exit(1);
  }
  return ia_ad(abs(x.val), x.der*sign(x.val));
}

ostream & operator << (ostream &os, const ia_ad &x)
{  return os << '(' << x.val << ',' << x.der << ')'; }

istream & operator >> (istream &is, ia_ad &x) 
{
  intervalw fx, df;

  is >> fx;
  is >> df;
  x = ia_ad(fx, df);
  return is;
}

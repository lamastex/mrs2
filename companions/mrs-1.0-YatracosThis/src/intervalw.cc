/*   File: intervalww.cc

     A C++ class for intervalw arithmetic, including some standard functions.
     We use directed rounding, and assume that the standard functions are 
     of maximal quality. This is NOT rigorous!

     Author: Warwick Tucker <warwick@math.uu.se>

     Compilation: g++ -Wall -O0 -c intervalw.cc
     Latest edit: Fri Jul  2 10:36:49 CEST 2004
*/

#include "intervalw.h"
#include "round.h"

static double Min(double dbl1, double dbl2)
{ return (dbl1 < dbl2 ? dbl1 : dbl2); }

static double Max(double dbl1, double dbl2)
{ return (dbl1 > dbl2 ? dbl1 : dbl2); }

intervalw operator + (const intervalw &iv1, const intervalw &iv2)
{
  intervalw result;
  setRoundDown();
  result.lo = iv1.lo + iv2.lo;
  setRoundUp();
  result.hi = iv1.hi + iv2.hi;
  setRoundNear();
  return result;
}

intervalw operator - (const intervalw &iv1, const intervalw &iv2)
{
  intervalw result;
  setRoundDown();
  result.lo = iv1.lo - iv2.hi;
  setRoundUp();
  result.hi = iv1.hi - iv2.lo;
  setRoundNear();
  return result;
}

intervalw operator * (const intervalw &iv1, const intervalw &iv2)
{
  intervalw result;
  double temp1, temp2;
  setRoundDown();
  temp1 = Min(iv1.lo * iv2.lo, iv1.lo * iv2.hi);
  temp2 = Min(iv1.hi * iv2.lo, iv1.hi * iv2.hi); 
  result.lo = Min(temp1, temp2);
  setRoundUp();
  temp1 = Max(iv1.lo * iv2.lo, iv1.lo * iv2.hi);
  temp2 = Max(iv1.hi * iv2.lo, iv1.hi * iv2.hi); 
  result.hi = Max(temp1, temp2);
  setRoundNear();
  return result;
}

intervalw operator / (const intervalw &iv1, const intervalw &iv2)
{ 
  if ( (iv2.lo <= 0.0) && (0.0 <= iv2.hi) ) {
    cerr << "Error: division by zero. Bye!" << endl;
    exit(1);
  }
  intervalw inv;
  setRoundDown();
  inv.lo = 1.0 / iv2.hi;
  setRoundUp();
  inv.hi = 1.0 / iv2.lo;
  setRoundNear();
  return (iv1 * inv); 
}

bool operator == (const intervalw &iv1, const intervalw &iv2)
{
  return ( (inf(iv1) == inf(iv2)) && (sup(iv1) == sup(iv2)) );
}

bool operator != (const intervalw &iv1, const intervalw &iv2)
{
  return !(iv1 == iv2);
}

bool subset (const intervalw &iv1, const intervalw &iv2)
{ return ( (iv2.lo <= iv1.lo) && (iv1.hi <= iv2.hi) ); }

bool strict_subset (const intervalw &iv1, const intervalw &iv2)
{ return ( (iv2.lo < iv1.lo) && (iv1.hi < iv2.hi) ); }

bool intersect (intervalw &iv, const intervalw &iv1, const intervalw &iv2)
{
  if ( (sup(iv1) < inf(iv2)) || (sup(iv2) < inf(iv1)) )
    return 0;
  else {
    iv = intervalw(max(inf(iv1), inf(iv2)), min(sup(iv1), sup(iv2)));
    return true;
  }
}

double diam (const intervalw &iv) 
{ 
  setRoundUp();
  double d = iv.hi - iv.lo;
  setRoundNear();
  return d;
}

double rad (const intervalw &iv) 
{ 
  setRoundUp();
  double r = 0.5*(iv.hi - iv.lo);
  setRoundNear();
  return r;
}

double mig (const intervalw &iv) 
{ 
  if ( (iv.lo <= 0.0) && (0.0 <= iv.hi) ) 
    return 0.0;
  else
    return Min( abs(iv.lo), abs(iv.hi) );
}

double mag (const intervalw &iv) 
{ 
 return Max( abs(iv.lo), abs(iv.hi) );
}

double mid (const intervalw &iv) 
{ 
 return 0.5*(sup(iv) + inf(iv));
}

intervalw exp (const intervalw &iv)
{
  intervalw result;
  setRoundDown();
  result.lo = exp(iv.lo);
  setRoundUp();
  result.hi = exp(iv.hi);
  setRoundNear();
  return result;
}

intervalw log (const intervalw &iv)
{
  if ( iv.lo < 0.0 ) {
    cerr << "Error: log of negative number. Bye!" << endl;
    exit(1);
  }
  intervalw result;
  setRoundDown();
  result.lo = log(iv.lo);
  setRoundUp();
  result.hi = log(iv.hi);
  setRoundNear();
  return result;
}

intervalw pow (const intervalw &iv, int n)
{
  if ( (iv.lo <= 0.0) && (0.0 <= iv.hi) && (n < 0) ) {
    cerr << "Error: negative power of zero. Bye!" << endl;
    exit(1);
  }
  intervalw result;
  if ( n > 0 ) {
    if ( n % 2 == 1 ) { // n is positive and odd.
      setRoundDown();
      result.lo = pow(iv.lo, n);
      setRoundUp();
      result.hi = pow(iv.hi, n);
    }
    else {            // n is positive and even.
      setRoundDown();
      result.lo = pow(mig(iv), n);
      setRoundUp();
      result.hi = pow(mag(iv), n);
    }
  }
  else if ( n < 0 )   // n is negative.
    return pow(1.0/iv, -n);
  else                // n is zero.
    return intervalw(1.0);
  setRoundNear();
  return result;
}

intervalw pow (const intervalw &iv, double a)
{ return exp(a*log(iv)); }

const double TWO_PI  = 8*atan(1.0);
const double PI_HALF = 2*atan(1.0);

intervalw sin (const intervalw &iv)
{
 if ( TWO_PI <= diam(iv) )
   return intervalw(-1, +1);

  intervalw ivNeg = iv / TWO_PI + 0.25;
  intervalw ivPos = iv / TWO_PI - 0.25;

  int  intLo = (int) floor(iv.lo/TWO_PI);
  int  intHi = (int)  ceil(iv.hi/TWO_PI);
  bool sNeg = false, sPos = false;
  for (int i = intLo; i <= intHi; i++) {
    if ( subset(i, ivNeg) ) sNeg = true;
    if ( subset(i, ivPos) ) sPos = true;
  }

  if ( sPos && sNeg )
    return intervalw(-1.0, +1.0);
  else if ( sPos && !sNeg )
    return intervalw(Min(sin(iv.lo), sin(iv.hi)), +1.0);
  else if ( !sPos && sNeg )
    return intervalw(-1.0, Max(sin(iv.lo), sin(iv.hi)));
  else // ( !sPos && !sNeg )
    return intervalw(Min(sin(iv.lo), sin(iv.hi)), Max(sin(iv.lo), sin(iv.hi)));
}

intervalw cos (const intervalw &iv)
{ return sin(iv + PI_HALF); }

intervalw tan (const intervalw &iv)
{ return sin(iv)/cos(iv); }

intervalw abs (const intervalw &iv)
{ return intervalw(mig(iv), mag(iv)); }

ostream & operator << (ostream &os, const intervalw &iv)
{  return os << '[' << iv.lo << ',' << iv.hi << ']'; }

istream & operator >> (istream &is, intervalw &iv) 
{
  double min, max;

  is >> min;
  is >> max;
  iv = intervalw(min, max);
  return is;
}

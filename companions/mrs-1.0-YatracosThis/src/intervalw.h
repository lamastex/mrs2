/*   File: intervalw.h

     The header for a C++ class for intervalw arithmetic, including some 
     standard functions. We use directed rounding, and assume that the 
     standard functions are of maximal quality. This is NOT rigorous!

     Author: Warwick Tucker <warwick@math.uu.se>

     Compilation: None.
     Latest edit: Fri Jul  2 10:36:49 CEST 2004
*/

#ifndef INTERVALW_H
#define INTERVALW_H

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

class intervalw
{
  friend intervalw operator + (const intervalw &iv1, const intervalw &iv2);
  friend intervalw operator - (const intervalw &iv1, const intervalw &iv2);
  friend intervalw operator * (const intervalw &iv1, const intervalw &iv2);
  friend intervalw operator / (const intervalw &iv1, const intervalw &iv2);

  friend bool operator == (const intervalw &iv1, const intervalw &iv2);
  friend bool operator != (const intervalw &iv1, const intervalw &iv2);

  friend bool     subset        (const intervalw &iv1, const intervalw &iv2);
  friend bool     strict_subset (const intervalw &iv1, const intervalw &iv2);
  friend bool     intersect     (      intervalw &iv, 
			         const intervalw &iv1, const intervalw &iv2);
  friend double   diam   (const intervalw &iv);
  friend double   rad    (const intervalw &iv);
  friend double   mig    (const intervalw &iv);
  friend double   mag    (const intervalw &iv);
  friend double   mid    (const intervalw &iv);
  friend double   inf    (const intervalw &iv){ return iv.lo; }
  friend double   sup    (const intervalw &iv){ return iv.hi; }

  friend intervalw exp (const intervalw &iv);
  friend intervalw log (const intervalw &iv);
  friend intervalw pow (const intervalw &iv, int    n);
  friend intervalw pow (const intervalw &iv, double a);
  friend intervalw sin (const intervalw &iv);
  friend intervalw cos (const intervalw &iv);
  friend intervalw tan (const intervalw &iv);
  friend intervalw abs (const intervalw &iv);

  friend ostream & operator << (ostream &, const intervalw &);
  friend istream & operator >> (istream &, intervalw &);
public:
  intervalw() {}
  intervalw(double val) {lo = val; hi = val;}
  intervalw(double min, double max) {lo = min; hi = max;}
private:
  double lo;
  double hi;
};

#endif // INTERVALW_H

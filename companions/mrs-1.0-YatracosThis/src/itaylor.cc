/*   File: itaylor.cc

     A Taylor series class using recursive formulas were obtained from the 
     books

     Aberth, O. "Precise Numerical Analysis" WCB Publishers, 1988.
     Aberth, O. "Precise Numerical Methods Using C++" Academic Press, 1998

     Author: Warwick Tucker <warwick@math.uu.se>
     Latest edit: Mon Jul  5 13:28:18 CEST 2004
*/

#include "itaylor.h"

static int Min(int n1, int n2)
{ return (n1 < n2 ? n1 : n2); }

static int Max(int n1, int n2)
{ return (n1 > n2 ? n1 : n2); }

itaylor::itaylor()
{
  order = -1;
  coeff = new interval[order + 1]; // Null pointer.
}

itaylor::itaylor(const itaylor &a)
{
  order = a.order;
  coeff = new interval[order + 1];
  memcpy(coeff, a.coeff, (order + 1) * sizeof(double));
}

itaylor::itaylor(const interval &constant)
{
  order = 0;
  coeff = new interval[order + 1];
  memcpy(coeff, &constant, (order + 1) * sizeof(interval));
}

itaylor::itaylor(const double &constant)
{
  interval iv(constant);
  order = 0;
  coeff = new interval[order + 1];
  memcpy(coeff, &iv, (order + 1) * sizeof(interval));
}

itaylor &itaylor::operator = (const itaylor &a)
{
  delete coeff;
  order = a.order;
  coeff = new interval[order + 1];
  memcpy(coeff, a.coeff, (order + 1) * sizeof(interval));
  return *this;
}

interval &itaylor::operator [] (int index) const
{
  if ( (index < 0) || (order < index) ) {
    cerr << "Error: Itaylor index out of range." << endl;
    exit(EXIT_FAILURE);
  }
  return coeff[index];
}

itaylor operator + (const itaylor &a)
{
  return a;
}

itaylor operator - (const itaylor &a)
{
  itaylor c(a);
  for (int i = 0; i <= c.order; i++)
      c[i] = -a[i];
  return c;
}

itaylor operator + (const itaylor &a, const itaylor &b)
{
  itaylor c;
  if ( a.order <= b.order ) {
    c = b;
    for (int i = 0; i <= a.order; i++)
      c[i] += a[i];
  }
  else {
    c = a;
    for (int i = 0; i <= b.order; i++)
      c[i] += b[i];
  }
  return c;
}

itaylor operator - (const itaylor &a, const itaylor &b)
{
  itaylor c;
  if ( a.order <= b.order ) {
    c = -b;
    for (int i = 0; i <= a.order; i++)
      c[i] += a[i];
  }
  else {
    c = a;
    for (int i = 0; i <= b.order; i++)
      c[i] -= b[i];
  }
  return c;
}
 
itaylor operator * (const itaylor &a, const itaylor &b)
{
  int    N = a.order + b.order;
  itaylor c = itaylor::constant(interval(0), N);
  for (int i = 0; i <= a.order; i++) 
    for (int j = 0; j <= b.order; j++) {
      c[i + j] += a[i] * b[j];
    }
  return c;
}
 
itaylor operator / (const itaylor &a, const itaylor &b)
{
  int    N = Max(a.order, b.order);
  itaylor c = itaylor::constant(interval(0), N);

  c[0] = a[0] / b[0];
  for (int n = 1; n <= N; n++) {
    interval sum = interval(0);
    if ( n <= a.order )
      sum = a[n];
    for (int i = 0; i <= n - 1; i++)
      if ( n - i <= b.order )
	sum -= c[i] * b[n - i];
    c[n] = sum / b[0];
  }
  return c;
}

itaylor exp (const itaylor &a)
{
  int    N = a.order;
  itaylor c = itaylor::constant(interval(0), N);
  c[0] = exp(a[0]);
  for (int k = 1; k <= N; k++) {
    interval sum = interval(0);
    for (int i = 1; i <= k; i++)
      sum += i * a[i] * c[k - i];
    c[k] = sum / k;
  }
  return c;
}
 
itaylor log (const itaylor &a)
{
  int    N = a.order;
  itaylor c = itaylor::constant(interval(0), N);
  c[0] = ln(a[0]);
  for (int k = 1; k <= N; k++) {
    interval sum = interval(0);
    for (int i = 1; i <= k - 1; i++)
      sum += i * c[i] * a[k - i];
    sum = a[k] - sum / k;
    c[k] = sum / a[0];
  }
  return c;
}
 
itaylor sin (const itaylor &a)
{
  int    N = a.order;
  itaylor cSin = itaylor::constant(interval(0), N);
  itaylor cCos = itaylor::constant(interval(0), N);
  cSin[0] = sin(a[0]);  cCos[0] = cos(a[0]);
  for (int n = 1; n <= N; n++) {
    interval sinSum = interval(0);
    interval cosSum = interval(0);
    for (int i = 1; i <= n; i++) {
      sinSum += i * a[i] * cCos[n - i];
      cosSum -= i * a[i] * cSin[n - i];
    }
    cSin[n] = sinSum / n;
    cCos[n] = cosSum / n;
  }
  return cSin;
}
 
itaylor cos (const itaylor &a)
{
  int    N = a.order;
  itaylor cSin = itaylor::constant(interval(0), N);
  itaylor cCos = itaylor::constant(interval(0), N);
  cSin[0] = sin(a[0]);  cCos[0] = cos(a[0]);
  for (int n = 1; n <= N; n++) {
    interval sinSum = interval(0);
    interval cosSum = interval(0);
    for (int i = 1; i <= n; i++) {
      sinSum += i * a[i] * cCos[n - i];
      cosSum -= i * a[i] * cSin[n - i];
    }
    cSin[n] = sinSum / n;
    cCos[n] = cosSum / n;
  }
  return cCos;
}
 
itaylor tan (const itaylor &a)
{
  return ( sin(a) / cos(a) );
}
 
itaylor pow (const itaylor &a, const itaylor &b)
{
  return ( exp( log(a) * b ) );
}

itaylor pow (const itaylor &a, const interval &t)
{
  int     N = a.order;
  itaylor c = itaylor::constant(interval(0), N);
  c[0] = pow(a[0], t); 
  for (int n = 1; n <= N; n++) {
    interval sum = interval(0);
    for (int i = 1; i <= n; i++)
      sum += (((t + 1) * i) / n - 1) * c[n - i] * a[i];
    c[n] = sum / a[0];
  }
  return c;
}
 
itaylor pow (const itaylor &a, const int &n)
{
  int     N = a.order;
  itaylor c = itaylor::constant(interval(0), N);
  if ( n <= 5 ) {     // The bound '5' is just a guess.
    c = a;
    for (int i = 1; i < n; i++)
      c = c * a;
  }
  else
    c = pow(a, interval(n));
  return c;
}
 
void resize(itaylor &a, int degree)
{
  delete a.coeff; 
  a.order = degree; 
  a.coeff = new interval[a.order + 1];
}
 
void clear(itaylor &a)
{
  for (int i = 0; i <= a.order; i++)
    a[i] = interval(0);
}
 
void round(itaylor &a, int degree)
{
  itaylor c = itaylor::constant(interval(0), degree);
  int    N = Min(a.order, c.order);
  for (int i = 0; i <= N; i++)
    c[i] = a[i];
  a = c;
}
 
itaylor round (const itaylor &a, int degree)
{
  itaylor c(a);
  round(c, degree);
  return c;
}

itaylor itaylor::variable (const interval &value, int degree = 1)
{
  if ( degree < 1 ) {
    cerr << "Warning: (variable) Itaylor degree must be at least 1." << endl;
    degree = 1;
  }
  itaylor c;
  resize(c, degree);
  clear(c);
  c[0] = value;
  c[1] = interval(1);
  return c;
}

itaylor itaylor::constant (const interval &value, int degree = 0)
{
  if ( degree < 0 ) {
    cerr << "Warning: (constant) Itaylor degree must be at least 0." << endl;
    degree = 0;
  }
  itaylor c;
  resize(c, degree);
  clear(c);
  c[0] = value;
  return c;
}

itaylor derivative (const itaylor &a)
{
  if ( a.order == 0 ) {
    itaylor c = itaylor::constant(interval(0));
    return c;
  }
  int N = a.order - 1;
  itaylor c = itaylor::constant(interval(0), N);
  for (int i = 0; i <= N; i++) 
    c[i] = (i + 1) * a[i + 1];
  return c;
}

itaylor derivative (const itaylor &a, int k)
{
  itaylor c(a);
  for (int i = 1; i <= k; i++) {
    itaylor b = derivative(c);
    c = b;
  }
  return c;
}
 
itaylor integral (const itaylor &a)
{
  int N = a.order + 1;
  itaylor c = itaylor::constant(interval(0), N);
  for (int i = 1; i <= N; i++) 
    c[i] = a[i - 1] / interval(i);
  return c;
}

itaylor integral (const itaylor &a, int k)
{
  itaylor c(a);
  for (int i = 1; i <= k; i++) {
    itaylor b = integral(c);
    c = b;
  }
  return c;
}

int orderOf (const itaylor &a)
{
  return a.getOrder();
}

ostream &operator << (ostream &oS, const itaylor &a)
{
  for (int i = 0; i <= a.order; i++) {
    oS.width(3); oS.precision(3); 
    oS.setf(ios_base::showpos); oS.setf(ios_base::fixed);
    oS << a[i] << "*t^"; 
    oS.unsetf(ios_base::showpos);
    oS << i << " ";
  }
  oS << endl;
  return oS;
}

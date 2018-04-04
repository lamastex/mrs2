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

//////////////////////////////////////////////////////////////
//
//       Implementation of class itaylor in itaylor.cpp     
//
//////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//     Updated by F. Blomquist, M. Grimmer
//     Extended version 05.03.2006 by M. Grimmer
////////////////////////////////////////////////////////////////////////

#include "itaylor.hpp"

///////////////////////////////////////////////////////////////
//
//                     class itaylor
//
///////////////////////////////////////////////////////////////

namespace taylor {

  //ivector itaylor::faks(1);
int itaylor::initialized=0;

void itaylor::initialize()
{
  /*
 Resize(itaylor::faks,0,170);
 itaylor::faks[0]=interval(1.0);
 itaylor::faks[1]=interval(1.0);

 for(int i=2; i<=170; i++)
   itaylor::faks[i]=itaylor::faks[i-1]*interval(i);
  */
}

//-----------------------------------------------------------------------

// Constructors:

itaylor::itaylor()
{
 if(!itaylor::initialized){itaylor::initialize();itaylor::initialized=1;};
}

//-----------------------------------------------------------------------

itaylor::itaylor(int order)
{
 if(!itaylor::initialized){itaylor::initialize();itaylor::initialized=1;};
 /*if(order<0 || order>170) 
  {
    std::cerr << "itaylor::itaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 */
 p=order;
 Resize(tayl,0,p);
}

//-----------------------------------------------------------------------

itaylor::itaylor(const itaylor& s)
{
 if(!itaylor::initialized){itaylor::initialize();itaylor::initialized=1;};
 p=s.p;
 Resize(tayl,0,p);
 tayl=s.tayl;
}

//-----------------------------------------------------------------------

itaylor::itaylor(int order, const real& value)
{
 if(!itaylor::initialized){itaylor::initialize();itaylor::initialized=1;};
 /*
 if(order<0 || order>170) 
  {
    std::cerr << "itaylor::itaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 */
 p=order;
 Resize(tayl,0,p);
 
 interval interval_value=interval(value);
 tayl[0]=interval_value;
 if(p>0)
   {
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

//-----------------------------------------------------------------------

itaylor::itaylor(int order, const interval& value)
{
 if(!itaylor::initialized){itaylor::initialize();itaylor::initialized=1;};
 /* 
if(order<0 || order>170) 
  {
    std::cerr << "itaylor::itaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 */
 p=order;
 Resize(tayl,0,p);
 tayl[0]=value;
 if(p>0)
   { 
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

//-----------------------------------------------------------------------

// Functions for initialization of independent variables:
itaylor var_itaylor(int ord, const real& x)
{
    itaylor erg(ord,x);
    return erg;
}

//-----------------------------------------------------------------------

itaylor var_itaylor(int ord, const interval& x)
{
    itaylor erg(ord,x);
    return erg;
}

//-----------------------------------------------------------------------


// Functions for initialization of constants:
itaylor const_itaylor(int ord, const real& c)
{
 itaylor erg(ord);
 erg.tayl[0]=interval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

itaylor const_itaylor(int ord, const interval& c)
{
 itaylor erg(ord);
 erg.tayl[0]=c;
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}
//-----------------------------------------------------------------------



//-----------------------------------------------------------------------
// assignment operators
//-----------------------------------------------------------------------

itaylor itaylor::operator=(const itaylor& s)
{
 p=s.p;
 Resize(tayl,0,s.p);
 tayl=s.tayl;
 return *this;
}

itaylor itaylor::operator=(int n)
{
 tayl[0] = interval(n);
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

itaylor itaylor::operator=(const real& x)
{
 tayl[0] = interval(x);
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

itaylor itaylor::operator=(const interval& x)
{
 tayl[0] = x;
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

itaylor itaylor::operator=(const ivector& iv) //added,mg,2005-08
                                              //const since C-XSC 2.1,mg,2006-02

{ 
 p=VecLen(iv)-1;
 Resize(tayl,0,p);
 tayl=iv;
 return *this;
}

//-----------------------------------------------------------------------
// relational operators
//-----------------------------------------------------------------------

int itaylor::operator==(itaylor& it)     //added, mg2005-08
{
  return ((p==it.p)&&(tayl==it.tayl));
}

int itaylor::operator!=(itaylor& it)     //added, mg2005-08
{
  return (!((p==it.p)&&(tayl==it.tayl)));
}

int itaylor::operator<=(itaylor& it)     //added, mg2005-08
{
  return ((p==it.p)&&(tayl<=it.tayl));
}

int itaylor::operator<(itaylor& it)      //added, mg2005-08
{
  return ((p==it.p)&&(tayl<it.tayl));
}

int itaylor::operator>=(itaylor& it)     //added, mg2005-08
{
  return ((p==it.p)&&(tayl>=it.tayl));
}

int itaylor::operator>(itaylor& it)      //added, mg2005-08
{
  return ((p==it.p)&&(tayl>it.tayl));
}

//-----------------------------------------------------------------------
// component access
//-----------------------------------------------------------------------

interval& itaylor::operator[](int n) //added, mg2005-08
// r/w access -> not const 
{
 return tayl[n];
}


//-----------------------------------------------------------------------

// class components:

// returning the maximal order
int get_order(const itaylor& x)
{
 return x.p;
}

//-----------------------------------------------------------------------

// returning all Taylor-coefficients by an interval vector
ivector get_all_coef(const itaylor& x)
{
 return x.tayl;
}

//-----------------------------------------------------------------------

// returning Taylor-coefficient of order j
interval get_j_coef(const itaylor& x, int j)
{
 return x.tayl[j];
}

//-----------------------------------------------------------------------
/*
// returning derivative of order j,  j <= 170;
interval get_j_derive(const itaylor& x, int j)
{
  return x.tayl[j]*itaylor::faks[j];
}
*/
//-----------------------------------------------------------------------

// Output of all taylor coefficients 
void print_itaylor(const itaylor& x)
{
 std::cerr <<"Output itaylor of order " << x.p << " " << std::endl;
 for(int i=Lb(x.tayl); i<=Ub(x.tayl);i++) 
  {
   std::cerr << "i  " << i << "  component: " << x.tayl[i] << std::endl;
  };
 std::cerr << std::endl;
}

void print_itaylor(std::ostream& os, const itaylor& x, int width, int digits)  
                                             // added, mg2005,2006
{
 os <<"Ausgabe itaylor der Ordnung " << x.p << " " << std::endl;
 for(int i=Lb(x.tayl); i<=Ub(x.tayl);i++) 
  {
   os << "i  " << i << "  component: ";
   if (width>0||digits>0) os << SetPrecision(width,digits);
   os << x.tayl[i] << std::endl;
  };
 os << std::endl;
}

//-----------------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, itaylor& x)
                                             // added, mg2005-08
{
 os <<"[itaylor object, order " << x.p << ":]" << std::endl;
 for(int i=Lb(x.tayl); i<=Ub(x.tayl);i++) 
  {
   os << "[" << i << "]";
   os << x.tayl[i] << std::endl;
  };
 os << std::endl;
 return os;
}                                


//-----------------------------------------------------------------------

// Overloading of operators:

//-----------------------------------------------------------------------

// - operator:
itaylor operator-(const itaylor& x)
{
 int order=get_order(x);
 itaylor erg(order);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= -x.tayl[j];
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:  +,-,*,/  for (itaylor, itaylor):
// All operands are independent variables.
itaylor operator-(const itaylor& x, const itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in itaylor, operator - : different orders " 
             << std::endl;
   exit(1);
  };

 itaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]-y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

itaylor operator+(const itaylor& x, const itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in itaylor, operator + : different orders " 
             << std::endl;
   exit(1);
  };

 itaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]+y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

itaylor operator*(const itaylor& x, const itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in itaylor, operator * : different orders " 
             << std::endl;
   exit(1);
  };

 itaylor erg(order1);
 interval sum; 
 idotprecision sum_idot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_idot=interval(0);
  for(int i=0; i<=j; i++)
   {
    accumulate(sum_idot, x.tayl[i],y.tayl[j-i]);
   }
  rnd(sum_idot,sum);
  erg.tayl[j]= sum;
 }
 return erg; 
}

//-----------------------------------------------------------------------

itaylor operator/(const itaylor& x, const itaylor& y)
{
 int order1(get_order(x));
 int order2(get_order(y));
 if(order1 != order2) 
  {
   std::cerr << "Error in itaylor, operator / : different orders " 
             << std::endl;
   exit(1);
  };
 
 if(0 <= y.tayl[0]) 
  {
   std::cerr << "Error in itaylor, operator / : 0 in interval " << std::endl;
   exit(1);
  };

 itaylor erg(order1);
 interval sum; 
 idotprecision sum_idot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_idot=interval(0);
  for(int i=1; i<=j; i++)
   {
    accumulate(sum_idot, y.tayl[i],erg.tayl[j-i]);
   }
  rnd(sum_idot,sum);
  erg.tayl[j]= (x.tayl[j]-sum)/y.tayl[0];
 }
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (interval, itaylor):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

itaylor operator-(const interval& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(const interval& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(const interval& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(const interval& x, const itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    interval sum;
    int order(get_order(y));
    itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (itaylor, interval):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

itaylor operator-(const itaylor& x, const interval& y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(const itaylor& x, const interval& y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(const itaylor& x, const interval& y)
{
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(const itaylor& x, const interval& y)
{
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (real, itaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

itaylor operator-(const real& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(const real& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(const real& x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(const real& x, const itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    interval sum;
    int order(get_order(y));
    itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (itaylor, real):
// The operand of type real is assumed to be a constant and not
// an independent variable!

itaylor operator-(const itaylor& x, const real& y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(const itaylor& x, const real& y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(const itaylor& x, const real& y)
{
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(const itaylor& x, const real& y)
{
    if (y==0)
    {
	std::cerr << "Error in itaylor: division by 0" 
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

// Operators with two operands:   +,-,*,/  for (int, itaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

itaylor operator-(int x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = -y;
    erg.tayl[0] = interval(x) - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(int x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    erg = y;
    erg.tayl[0] = interval(x) + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(int x, const itaylor& y)
{
    int order(get_order(y));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = interval(x)*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(int x, const itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    interval sum;
    int order(get_order(y));
    itaylor w(order);
    w.tayl[0] = interval(x) / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (itaylor, int):
// The operand of type real is assumed to be a constant and not
// an independent variable!

itaylor operator-(const itaylor& x, int y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - interval(y);
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator+(const itaylor& x, int y)
{
    int order(get_order(x));
    itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + interval(y);
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator*(const itaylor& x, int y)
{
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*interval(y);
    return erg;
}

//-----------------------------------------------------------------------

itaylor operator/(const itaylor& x, int y)
{
    if (y==0)
    {
	std::cerr << "Error in itaylor: division by 0" 
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/interval(y);
    return erg;
}

//-----------------------------------------------------------------------

// Overloading the standard functions:

//-----------------------------------------------------------------------

// Help function
void f_g_u(const itaylor& f, const itaylor& g, const itaylor& u, int nb_function)
{
 int order1=get_order(f);
 int order2=get_order(g);
 int order3=get_order(u);

 // The following errors should be caught before 
 // but for security here again:
 if(order1 != order2) 
  {
   std::cerr << "Error1 in f_g_u: different orders " 
             << std::endl;
   exit(1);
  };

 if(order3 != order2) 
  {
   std::cerr << "Error2 in f_g_u: different orders " << std::endl;
   exit(1);
  };

 if(0 <= g.tayl[0]) 
  {
   std::cerr << "Error in f_g_u : wrong argument " << std::endl;
   exit(1);
  }; 

 switch(nb_function) // element No. 0
   {
    case _i_ln:f.tayl[0]=ln(u.tayl[0]); break;

    case _i_tan:f.tayl[0]=tan(u.tayl[0]); break;
    case _i_cot:f.tayl[0]=cot(u.tayl[0]); break;

    case _i_asin:f.tayl[0]=asin(u.tayl[0]); break;
    case _i_acos:f.tayl[0]=acos(u.tayl[0]); break;
    case _i_atan:f.tayl[0]=atan(u.tayl[0]); break;
    case _i_acot:f.tayl[0]=acot(u.tayl[0]); break;

    case _i_tanh:f.tayl[0]=tanh(u.tayl[0]); break;
    case _i_coth:f.tayl[0]=coth(u.tayl[0]); break; 

    case _i_asinh:f.tayl[0]=asinh(u.tayl[0]); break;
    case _i_acosh:f.tayl[0]=acosh(u.tayl[0]); break;
    case _i_atanh:f.tayl[0]=atanh(u.tayl[0]); break;

    case _i_acoth:f.tayl[0]=acoth(u.tayl[0]); break;
   }
 
 // remaining elements:
 interval sum; 
 for(int j=1; j<=Ub(f.tayl); j++) 
 {
     sum = interval(0);
     for(int i=1; i<=j-1; i++)
	 sum += interval(i)*f.tayl[i]*g.tayl[j-i];
     f.tayl[j] = (u.tayl[j]-sum/interval(j)) / g.tayl[0];
 }
}

// sqr-function
itaylor sqr(const itaylor& x)
{
    idotprecision idot;
    interval sum;
    int order(get_order(x)),m;
    itaylor erg(order);

    erg.tayl[0] = sqr(x.tayl[0]);
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	idot = 0;
	for (int j=0; j<=m-1; j++) 
	    accumulate(idot,x.tayl[j],x.tayl[k-j]);
	rnd(idot,sum);
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(x.tayl[m]); // k even 
    }
    return erg; 
}

//-----------------------------------------------------------------------

// Square-root

itaylor sqrt(const itaylor& x)
{
    idotprecision idot;
    interval sum,h;
    int order(get_order(x)),m;
    itaylor erg(order);
    if (0<=x.tayl[0]) 
    {
	std::cerr << "Error in itaylor, sqrt: 0 in interval" << std::endl;
	exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0]);
    h = erg.tayl[0];
    times2pown(h,1);
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	idot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	rnd(idot,sum);
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]) / h;
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(x,n)

itaylor sqrt(const itaylor& x, int n)
{
    int order(get_order(x));
    itaylor erg(order);
    if (0<=x.tayl[0]) 
    {
	std::cerr << "Error in itaylor, sqrt(x,n): 0 in interval" << std::endl;
	exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0],n); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (interval(k-j)/real(n)-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}
//-----------------------------------------------------------------------

// sqrt(1-x^2):
itaylor sqrt1mx2(const itaylor& x)
{
    idotprecision idot;
    interval sum,h;
    int order(get_order(x)),m;
    itaylor erg(order), g(order);

    if (Inf(x.tayl[0])<=-1 || Sup(x.tayl[0])>=1)
    {
	std::cerr << "Error in itaylor, sqrt1mx2: wrong argument" << std::endl;
	exit(1);
    };
    erg.tayl[0]=sqrt(real(1)-sqr(x.tayl[0])); // =sqrt1mx2(x.tayl[0]); Blomi 
    h = real(-1)/erg.tayl[0];
    times2pown(h,-1);
    g = sqr(x);
    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	idot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	rnd(idot,sum);
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (g.tayl[k]+erg.tayl[k])*h;
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(x^2-1):
itaylor sqrtx2m1(const itaylor& x)
{
    const real c = 30.0; 
    idotprecision idot;
    interval sum,h;
    int order(get_order(x)),m;
    itaylor erg(order), g(order);

    if (Disjoint(x.tayl[0],interval(-1,1))==0)
    {
	std::cerr << "Error in itaylor, sqrtx2m1: wrong argument" << std::endl;
	exit(1);
    };

    if (Inf(x.tayl[0])>c) erg = x*sqrt1mx2(real(1)/x);
    else if (Sup(x.tayl[0])<-c) erg = -x*sqrt1mx2(real(1)/x);
    else {
	erg.tayl[0]=sqrt(sqr(x.tayl[0])-real(1)); // =sqrtx2m1(x.tayl[0]); Blomi 
	g = sqr(x);
	h = real(1)/erg.tayl[0];
	times2pown(h,-1);
	for (int k=1; k<=order; k++)
	{
	    m = (k+1)/2;
	    idot = 0;
	    for (int j=1; j<=m-1; j++) 
		accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	    rnd(idot,sum);
	    times2pown(sum,1);  // Multiplication with 2
	    erg.tayl[k] = sum;
	    if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	    erg.tayl[k] = (g.tayl[k]-erg.tayl[k])*h;
	}
    }
    return erg;
}

// sqrt(1+x)-1
itaylor sqrtp1m1(const itaylor& x)
{
    int order(get_order(x)),m;
    itaylor erg(order);
    idotprecision idot;
    interval h,Ne;

    if (Inf(x.tayl[0])<=-1)
    {
	std::cerr << "Error in itaylor, sqrtp1m1: wrong argument" << std::endl;
	exit(1);
    };

    erg.tayl[0] = sqrtp1m1(x.tayl[0]);
    Ne = real(1.0)+erg.tayl[0];
    times2pown(Ne,1);

    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	idot = 0;
	for (int j=1; j<=m-1; j++)
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	rnd(idot,h);
	times2pown(h,1);
	erg.tayl[k] = h;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]);
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]) / Ne; 
    }
    return erg;
} // sqrtp1m1

//-----------------------------------------------------------------------

// power-function
itaylor pow(const itaylor& x, const interval& alpha)
{
    int order(get_order(x));
    itaylor erg(order);

    if (0<=x.tayl[0])
    {
	std::cerr << "Error in itaylor, pow(x,a): 0 in interval x" 
                  << std::endl;
	exit(1);
    };

    erg.tayl[0] = pow(x.tayl[0],alpha); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (interval(k-j)*alpha-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}


//-----------------------------------------------------------------------

// Exponential-function
itaylor exp(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order);

    erg.tayl[0] = exp(x.tayl[0]); // element No. 0; function value
    for(int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*erg.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] /= interval(k);
    }
    return erg; 
}

//-----------------------------------------------------------------------

// exp(x)-1;
itaylor expm1(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order);

    erg.tayl[0] = exp(x.tayl[0]); // element No. 0; function value
    for(int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*erg.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] /= interval(k);
    }
    erg.tayl[0] = exp(x.tayl[0])-real(1); // = expm1(x.tayl[0]); Blomi
    return erg; 
}

//-----------------------------------------------------------------------

// Logarithm-function
itaylor ln(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);

    f_g_u(f,x,x,_i_ln);
  
    return f;
}

//-----------------------------------------------------------------------

// ln(1+x)
itaylor lnp1(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order), g(order);

    g = interval(1) + x;
    erg.tayl[0] = lnp1(x.tayl[0]);
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=1; j<=k-1; j++)
	    erg.tayl[k] += interval(j) * erg.tayl[j] * g.tayl[k-j];
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]/interval(k)) / g.tayl[0];
    }
    return erg;
}

//-----------------------------------------------------------------------

// Sinus-function
itaylor sin(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg1(order);   // sin
    itaylor erg2(order);   // cos
    interval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= real(-1.0)/interval(j)*s2;
    }
    return erg1; 
}

//-----------------------------------------------------------------------

// Cosinus-function
itaylor cos(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg1(order);   // sin
    itaylor erg2(order);   // cos
    interval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= real(-1.0)/interval(j)*s2;
    }
 return erg2; 
}

//-----------------------------------------------------------------------

// Tangens-function
itaylor tan(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=sqr(cos(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, tan : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_tan);
  
    return f;
}

//-----------------------------------------------------------------------

// Cotangens-function
itaylor cot(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=-sqr(sin(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, cot : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_cot);
  
    return f;
}

//-----------------------------------------------------------------------

// Sinushyperbolicus-function
itaylor sinh(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg1(order);  // sinh
    itaylor erg2(order);  // cosh
    interval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= s2/interval(j);
    }
    return erg1; 
}

//-----------------------------------------------------------------------

// Cosinushyperbolicus-function
itaylor cosh(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg1(order); // sinh
    itaylor erg2(order); // cosh
    interval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remaining elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= s2/interval(j);
    }
    return erg2; 
}

//-----------------------------------------------------------------------

//Tangenshyperbolicus-function
itaylor tanh(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=sqr(cosh(x)); 
 
    f_g_u(f,g,x,_i_tanh);
  
    return f;
}

//-----------------------------------------------------------------------

// Cotangenshyperbolicus-function
itaylor coth(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=-sqr(sinh(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, coth : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_coth);
  
    return f;
}

//-----------------------------------------------------------------------

//Arcsinusfunktion
itaylor asin(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, asin : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_asin);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccosinusfunktion
itaylor acos(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g = -sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, acos : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acos);
  
    return f;
}

//-----------------------------------------------------------------------

//Arctan-function
itaylor atan(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=interval(1.0)+sqr(x);

    f_g_u(f,g,x,_i_atan);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccotan-function
itaylor acot(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=-(interval(1.0)+sqr(x));

    f_g_u(f,g,x,_i_acot);
  
    return f;
}

//-----------------------------------------------------------------------

//Areasinh-function
itaylor asinh(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=sqrt1px2(x);

    f_g_u(f,g,x,_i_asinh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacosh-function
itaylor acosh(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=sqrtx2m1(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, acosh : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acosh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areatanh-function
itaylor atanh(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, atanh : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_atanh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacotanh-function
itaylor acoth(const itaylor& x)
{
    int order=get_order(x);
    itaylor f(order);
    itaylor g(order);

    g=interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in itaylor, acoth : wrong argument " << std::endl;
	exit(1);
    };  

    f_g_u(f,g,x,_i_acoth);
  
    return f;
}

//-----------------------------------------------------------------------

//Error function "erf" //added, mg2006-03
itaylor erf(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order);
    itaylor g(order);

    g=exp(-sqr(x));
   
    erg.tayl[0] = erf(x.tayl[0]); // element No. 0; function value
    
    for(int k=1; k<=order; k++)
    {
        erg.tayl[k] = 0;
        
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*g.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] = 2*erg.tayl[k]/(sqrt(Pi())*interval(k));
    }
    return erg; 
}

//-----------------------------------------------------------------------

//Complementary Error function "erfc" //added, mg2006-03
itaylor erfc(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order);

    erg=interval(1)-erf(x); 

    return erg;
}


//-----------------------------------------------------------------------

// sqrt(1+x^2)
itaylor sqrt1px2(const itaylor& x)
{
    int order(get_order(x));
    itaylor erg(order);
    const real c = 500.0;

    if (Inf(x.tayl[0]) > c) erg = x*sqrt(real(1)+real(1)/sqr(x));
    else if (Sup(x.tayl[0]) < -c) erg = -x*sqrt(real(1)+real(1)/sqr(x)); 
    else erg = sqrt(real(1)+sqr(x));  
    return erg;
}



//-----------------------------------------------------------------------

//ALTERED by Tomas Johnson
////////////////////////////////////////////////////////////////////////////////////
itaylor powerAtZero(const itaylor &x, int n) {
  int order(get_order(x));
  itaylor erg(order);
  
  if (0<=x.tayl[1])
    {
      std::cerr << "Error in itaylor, pow(x,a): 0 in the derivative of x" 
		<< std::endl;
      exit(1);
    };
  
  for (int k=0; k<n; k++)
    erg.tayl[k]=0.0;

  erg.tayl[n] = power(x.tayl[1],n); // element No. n
  for (int k=n+1; k<=order; k++)
    {
      erg.tayl[k] = 0;
      for (int j=0; j<=k-1-n; j++)
	erg.tayl[k] += (interval(k-n-j)*n-interval(j))
		           * erg.tayl[j+n] * x.tayl[k+1-n-j];
      erg.tayl[k] /= (interval(k-n)*x.tayl[1]);
    }
  return erg;
}
////////////////////////////////////////////////////////////


} // End of namespace taylor


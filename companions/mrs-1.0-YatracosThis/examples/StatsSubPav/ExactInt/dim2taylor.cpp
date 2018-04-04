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

#include "dim2taylor.hpp"

using namespace std;
using namespace cxsc;

//------------------------------------------------------------------------

// class dim2taylor: 2-dim. Taylor arithmetic 

//------------------------------------------------------------------------


namespace taylor{
// Constructors:

dim2taylor::dim2taylor()
{
 p=0;
 dat=new ivector[p+1];
 for(int i=0; i<=p; i++) Resize(dat[i], 0, p-i);
}

//-----------------------------------------------------------------------

dim2taylor::dim2taylor(int order)
{
 p=order;
 dat=new ivector[p+1];
 for(int i=0; i<=p; i++) Resize(dat[i], 0, p-i);
}

//------------------------------------------------------------------------

dim2taylor::dim2taylor(const dim2taylor& s)
{
 p=s.p;
 dat=new ivector[p+1];
 for(int i=0; i<=p ;i++) dat[i]=s.dat[i];
}
//------------------------------------------------------------------------

dim2taylor::~dim2taylor()
{
 delete[] dat;
 dat=NULL;
}

//------------------------------------------------------------------------
// Assignment operator:

dim2taylor& dim2taylor::operator=(const dim2taylor& s)
{
 if(this != &s)
 {
  delete[] dat;
  dat=NULL;
  p=s.p;
  dat=new ivector[p+1];

  for(int i=0; i<=p ; i++) dat[i]=s.dat[i]; 

 }
 return *this;
}

//------------------------------------------------------------------------

ivector& dim2taylor::operator[](int n) const
{
    return dat[n];
} 

//------------------------------------------------------------------------
  
dim2taylor init_var(int order, int nr, const interval& value)
{

 dim2taylor t(order);

 if( (nr<1) && (nr>2) ) std::cerr << "Error in dim2taylor::init_var" 
                                  << std::endl;

 t[0][0]=value;

 if(t.p>0)
  {
   if(nr==1) t[1][0]=interval(1.0);
   else t[1][0]=interval(0.0);

   if(nr==2) t[0][1]=interval(1.0);
   else t[0][1]=interval(0.0);
   
   for(int i=2; i<=t.p; i++) {t[0][i]=interval(0.0);}   // Rest 0. line
   for(int i=1; i<=t.p-1; i++) {t[1][i]=interval(0.0);} // Rest 1. line  

   // remaining elements
   for(int j=2; j<=t.p ; j++) 
     for(int i=0; i<=t.p-j ; i++) t[j][i]=interval(0.0);       
  }

 return t;
}

//------------------------------------------------------------------------

dim2taylor init_const(int order, const interval& value)
{

 dim2taylor t(order);

 if(t.p>0)
  {
   for(int j=0; j<=t.p ; j++) 
     for(int i=0; i<=t.p-j ; i++) t[j][i]=interval(0.0);       
  }

 t[0][0]=value;

 return t;
}

//------------------------------------------------------------------------

void dim2taylor::print_dim2taylor() // debug
{
    for(int j=0; j<=p ; j++) 
    {
	for(int i=0; i<=p-j ; i++) std::cout << dat[j][i]<< " ";
	std::cout << std::endl;
    }
    std::cout << std::endl;
}

//------------------------------------------------------------------------

void dim2taylor::print_dim2taylor(std::ostream& os) 
                                 // overloaded for ostream parameter,
                                 // mg2005-08/2005-11
{
 for(int j=0; j<=p ; j++) 
  {
   for(int i=0; i<=p-j ; i++) std::cout << dat[j][i]<< " ";
   os << std::endl;
  }
 os << std::endl;
}

//------------------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, dim2taylor& d2t)
                                 // added, mg2005-08/2005-11
{
 os <<"[dim2taylor object, order " << d2t.p << ":]" << std::endl;
 for(int j=0; j<=d2t.p ; j++) 
  {
   for(int i=0; i<=d2t.p-j ; i++) std::cout << d2t.dat[j][i]<< " ";
   os << std::endl;
  }
 os << std::endl;
 return os;
}                                 


//------------------------------------------------------------------------

// Overloading of the elementary operators 

//------------------------------------------------------------------------

dim2taylor operator-(const dim2taylor& s)
{
 int order=s.p;
 dim2taylor erg(order);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=-s[j][i];
  }
 
 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const dim2taylor& s, const dim2taylor& t) 
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     std::cerr << "dim2taylor operator- : Operands with different orders";
     std::cerr << std::endl;
     exit(1);
   }

 dim2taylor erg(order1);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=s[j][i]-t[j][i];
  }
 
 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator+(const dim2taylor& s, const dim2taylor& t)
{
int order1=s.p;
int order2=t.p;

 if(order1 != order2)
   {
     cerr << "dim2taylor operator+ : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 dim2taylor erg(order1);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=s[j][i]+t[j][i];
  }

 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator*(const dim2taylor& s, const dim2taylor& t)
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     cerr << "dim2taylor operator* : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 dim2taylor erg(order1);

 idotprecision sum_idot;

 for(int k=0; k<=erg.p; k++)
  {
    for(int i=0; i<=k; i++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i; l++) // calculating erg(i,k-i)
	{
	 for(int m=0; m<=k-i; m++)
	  {
	    accumulate(sum_idot, s[l][m], t[i-l][k-i-m]);
	  } 
	}  

      rnd(sum_idot, erg[i][k-i]);      
     }
  }

 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator/(const dim2taylor& s, const dim2taylor& t)
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     cerr << "dim2taylor operator/ : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 dim2taylor erg(order1);

 idotprecision sum_idot;

 interval h, sum;

 if(0.0 <= t[0][0]) 
  {
    cerr << "dim2taylor operator/ : 0 in denominator" << endl;
    exit(1);
  }
 h=interval(1.0)/t[0][0];
 

 for(int k=0; k<=erg.p; k++) // calculating erg(i,k-i)
  {
    for(int i=0; i<=k; i++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i; l++) // Without Coeff. (l,m)=(0,0)
	{
	 for(int m=1; m<=k-i; m++) // m >= 1
	  {
	    accumulate(sum_idot, t[l][m], erg[i-l][k-i-m]);
	  }  //for m
	}    //for l

      for(int l=1; l<=i; l++) //  m=0, l!=0
	{
	 accumulate(sum_idot, t[l][0], erg[i-l][k-i]);
	}  // for l

      rnd(sum_idot, sum);  
      erg[i][k-i]=h*(s[i][k-i]-sum);
     }
  }
 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const interval& s, const dim2taylor& t)
{
 dim2taylor erg(t.p);
 dim2taylor s_ty=init_const(t.p, s);

 erg=s_ty-t;
 return erg;
}

dim2taylor operator+(const interval& s, const dim2taylor& t)
{
 dim2taylor erg(t.p);
 dim2taylor s_ty = init_const(t.p, s);

 erg = s_ty + t;
 return erg;
}

dim2taylor operator*(const interval& s, const dim2taylor& t)
{
 dim2taylor erg(t.p);
 dim2taylor s_ty=init_const(t.p, s);

 erg=s_ty*t;
 return erg;
}

dim2taylor operator/(const interval& s, const dim2taylor& t)
{
 dim2taylor erg(t.p);
 dim2taylor s_ty=init_const(t.p, s);

 erg=s_ty/t;
 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const dim2taylor& s, const interval& t)
{
 dim2taylor erg(s.p);
 dim2taylor t_ty=init_const(s.p, t);

 erg=s-t_ty;
 return erg;
}

dim2taylor operator+(const dim2taylor& s, const interval& t)
{
 dim2taylor erg(s.p);
 dim2taylor t_ty = init_const(s.p, t);

 erg=s+t_ty;
 return erg;
}

dim2taylor operator*(const dim2taylor& s, const interval& t)
{
 dim2taylor erg(s.p);
 dim2taylor t_ty=init_const(s.p, t);

 erg=s*t_ty;
 return erg;
}

dim2taylor operator/(const dim2taylor& s, const interval& t)
{
 dim2taylor erg(s.p);
 dim2taylor t_ty=init_const(s.p, t);

 erg=s/t_ty;
 return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const real& s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);
    erg = s_ty - t; 
    return erg;
}

// -------------------------------------------------------

dim2taylor operator+(const real& s, const dim2taylor& t)
{    
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);
    erg = s_ty + t;
    return erg;
}


dim2taylor operator*(const real& s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i); 
    erg = s_ty * t; 
    return erg;
}

dim2taylor operator/(const real& s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);  
    erg = s_ty / t; 
    return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const dim2taylor& s, const real& t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s - t_ty; 
 return erg;
}

dim2taylor operator+(const dim2taylor& s, const real& t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s + t_ty; 
    return erg;
}

dim2taylor operator*(const dim2taylor& s, const real& t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s * t_ty;
    return erg;
}

dim2taylor operator/(const dim2taylor& s, const real& t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s / t_ty; 
    return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(int s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);
    erg = s_ty - t; 
    return erg;
}

dim2taylor operator+(int s, const dim2taylor& t)
{    
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);
    erg = s_ty + t;
    return erg;
}


dim2taylor operator*(int s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i); 
    erg = s_ty * t; 
    return erg;
}

dim2taylor operator/(int s, const dim2taylor& t)
{   
    dim2taylor erg(t.p);
    interval s_i(s);
    dim2taylor s_ty = init_const(t.p,s_i);  
    erg = s_ty / t; 
    return erg;
}

//------------------------------------------------------------------------

dim2taylor operator-(const dim2taylor& s, int t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s - t_ty; 
 return erg;
}

dim2taylor operator+(const dim2taylor& s, int t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s + t_ty; 
    return erg;
}

dim2taylor operator*(const dim2taylor& s, int t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s * t_ty;
    return erg;
}

dim2taylor operator/(const dim2taylor& s, int t)
{   
    dim2taylor erg(s.p);
    interval t_i(t);
    dim2taylor t_ty = init_const(s.p,t_i);  
    erg = s / t_ty; 
    return erg;
}


//------------------------------------------------------------------------

// Overloading the elementary functions

//------------------------------------------------------------------------


//------------------------------------------------------------------------

dim2taylor sqr(const dim2taylor& s)
{
 dim2taylor erg(s.p);
 
 idotprecision sum_idot;
 interval sum1=interval(0.0);
 interval sum2=interval(0.0);

 erg[0][0]=sqr(s[0][0]); //Koeff. (0,0)

 if(erg.p > 0)
  {
    
    for(int k=0; k<=erg.p; k++) // calculating the remaining coefficients
      {
	for(int i=0; i<=k; i++)
	  {
	    if(i%2==1) // i: odd
	      {
		sum_idot=interval(0.0);
		for(int l=0; l<=(i-1)/2; l++) // calculating erg(i,k-i)
		  {
		    for(int m=0; m<=k-i; m++)
		      {
			accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
		      } //for m
		  }     //for l

		rnd(sum_idot, sum1);
		erg[i][k-i]=interval(2.0)*sum1;	  	
	      }
	    else //i%2==0: i: even 
	      {
		if( (k-i)%2==1 ) //k-i: odd
		  {
		    sum_idot=interval(0.0);	
		    for(int l=0; l<=(i-2)/2; l++) // calculating erg(i,k-i)
		      {
			for(int m=0; m<=k-i; m++)
			  {
			    accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
			  } //for m
		      }     //for l
		    rnd(sum_idot, sum1);
       
		    sum_idot=interval(0.0);
		    for(int m=0; m<=(k-i-1)/2; m++)
		      {
			accumulate(sum_idot, s[i/2][m], s[i/2][k-i-m]);
		      }//for m
		    rnd(sum_idot, sum2);

		    erg[i][k-i]=interval(2.0)*(sum1+sum2);  
             
		  }
		else //(k-i)%2==0: k-i even
		  {
		    sum_idot=interval(0.0);	
		    for(int l=0; l<=(i-2)/2; l++) // calculating erg(i,k-i)
		      {
			for(int m=0; m<=k-i; m++)
			  {
			    accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
			  } // for m
		      }     // for l
		    rnd(sum_idot, sum1);
       
		    sum_idot=interval(0.0);
		    for(int m=0; m<=(k-i-2)/2; m++)
		      {
			accumulate(sum_idot, s[i/2][m], s[i/2][k-i-m]);
		      } // for m
		    rnd(sum_idot, sum2);

		    erg[i][k-i]=interval(2.0)*(sum1+sum2);
		    erg[i][k-i]=erg[i][k-i]+sqr(s[i/2][(k-i)/2]);   
		  }
	      }
	  } // for i
      }     // for k
  }
 return erg;
}

//------------------------------------------------------------------------

dim2taylor sqrt(const dim2taylor& s)
// New fast version, Blomquist 17.10.2005;
// In the sum (5.25) of Braeuers thesis by twos summands are equal,
// so the evaluation can be simplified. The first and the last
// summands are equal and must not be added. After the simplification
// only the first summand must not be added, see:
// if (z==1) continue; e.g. skipping the first summand.
{ 
 dim2taylor erg(s.p);
 idotprecision sum_idot;
 interval sum(0.0);

 erg[0][0] = sqrt(s[0][0]);  // Coeff. (0,0) --> function value
  
  if(erg.p > 0)
  { 
   if(0.0 <= erg[0][0])
      {
	cout << "error here: " << erg[0][0] << endl;
	 
	cerr << "Error in dim2taylor sqrt: 0 in argument interval";
	cerr << endl;
	exit(1);
      }
    	
		
    interval h = erg[0][0]; 
    times2pown(h,1); // fast multiplication with 2;
    int ki,il,z;
    for(int k=1; k<=erg.p; k++) // calculating all other coefficients
    {   
	for(int i=0; i<=k; i++)
	{   
	    sum_idot=interval(0.0);
    // do not sum (l,m)=(0,0) and (l,m)=(i,k-i), see continue assignments
	    ki = k-i;   z = 0; // z: Numbering the following summands 
	    if (i%2==1) // i: odd
	    { 
		for(int l=0; l<=(i-1)/2; l++) 
		{
		    il = i-l;   
		    for(int m=0; m<=ki; m++)
		    {
			z++; // Numbering the summands, z = 1,2,3,...
			// (l,m) = (0,0) is the first summand (z=1),
			if (z==1) continue; // skipping accumulate(...)
			accumulate(sum_idot, erg[l][m], erg[il][ki-m]);
		    }
		}
		rnd(sum_idot,sum);
		times2pown(sum,1);
		erg[i][ki] = (s[i][ki]-sum)/h;
	    }
	    else // i: even
	    {
		if (ki%2==1) // ki=k-i: odd
		{
		    for(int l=0; l<=(i-2)/2; l++) 
		    {
			il = i-l;
			for(int m=0; m<=ki; m++)
			{
			    z++;
			    if (z==1) continue;
			    accumulate(sum_idot,erg[l][m],erg[il][ki-m]);
			}
		    }
		    for(int m=0; m<=(ki-1)/2; m++)
		    {
			z++;
			if (z==1) continue;
			accumulate(sum_idot,erg[i/2][m],erg[i/2][ki-m]);
		    }
		    rnd(sum_idot,sum);
		    times2pown(sum,1);
		    erg[i][ki] = (s[i][ki]-sum)/h;
		}
		else // ki=k-i: even
		{
		    for(int l=0; l<=(i-2)/2; l++)
		    {
			il = i-l;
			for(int m=0; m<=ki; m++)
			{
			    z++;
			    if (z==1) continue;
			    accumulate(sum_idot,erg[l][m],erg[il][ki-m]);
			}
		    }
		    for(int m=0; m<=(ki-2)/2; m++)
		    {
			z++;
			if (z==1) continue;
			accumulate(sum_idot,erg[i/2][m],erg[i/2][ki-m]);
		    }
		    rnd(sum_idot,sum);
		    times2pown(sum,1);
		    sum += sqr(erg[i/2][ki/2]);
		    erg[i][ki] = (s[i][ki]-sum)/h;
		}
	    }
	} // for i
      }   // for k
  }

 return erg;
 }

//------------------------------------------------------------------------

dim2taylor sqrt1px2(const dim2taylor& x){
    // sqrt(1+x^2);
    dim2taylor erg;
    if (Inf(x[0][0])>1) erg = x * sqrt(1+sqr(1/x));
    else 
	if (Sup(x[0][0])<-1)
	    erg = -x * sqrt(1+sqr(1/x));
	else erg = sqrt(1+sqr(x));
    return erg; 
}

//------------------------------------------------------------------------

dim2taylor sqrtx2m1(const dim2taylor& x){
    // sqrt(x^2-1);
    dim2taylor erg;
    if (Inf(x[0][0])>4) erg = x * sqrt( 1-sqr(1/x) );
    else 
	if (Sup(x[0][0])<-4)
	    erg = -x * sqrt( 1-sqr(1/x) );
	else erg = sqrt(sqr(x)-1);
    return erg; 
    }


//------------------------------------------------------------------------

dim2taylor pow(const dim2taylor& s, const interval& alpha)
{
 dim2taylor erg(s.p);

 idotprecision sum_idot;
 interval sum1, sum2, h;

 erg[0][0]=pow(s[0][0], alpha); 
 
 for(int j=1; j<=erg.p; j++)
   {
     sum_idot=interval(0.0);
     for(int i=0; i<=j-1; i++)
       {
	 h=alpha*(interval(j)-interval(i))-interval(i);
	 accumulate(sum_idot, h*erg[0][i], s[0][j-i]);
       }
     rnd(sum_idot, sum1);
     erg[0][j]=sum1/interval(j)/s[0][0];
   }

 for(int i=1; i<=erg.p; i++) // now calculating the remaining erg(i,k)
  {
    for(int k=0; k<=erg.p-i; k++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i-1; l++) // Koeff. (l,m)=(0,0) nicht summieren
	{
	 h=alpha*(interval(i)-interval(l))-interval(l);

	 for(int m=0; m<=k; m++)
	  {
	    accumulate(sum_idot, h*erg[l][m], s[i-l][k-m]);
	  }
	}
      rnd(sum_idot, sum1); 
        
      sum_idot=interval(0.0);
      for(int m=1; m<=k; m++) 
	{
	 accumulate(sum_idot, s[0][m], erg[i][k-m]);
	}
      rnd(sum_idot, sum2); 

      erg[i][k]=(sum1/interval(i)-sum2)/s[0][0];
     }
  } 
 
 return erg;
}

//------------------------------------------------------------------------

dim2taylor power(const dim2taylor& s, int n)
{
 dim2taylor erg(s.p);

 idotprecision sum_idot;
 interval sum1, sum2, h;

 erg[0][0]=power(s[0][0], n); 
 
 for(int j=1; j<=erg.p; j++)
   {
     sum_idot=interval(0.0);
     for(int i=0; i<=j-1; i++)
       {
	 h=interval(n)*(interval(j)-interval(i))-interval(i);
	 accumulate(sum_idot, h*erg[0][i], s[0][j-i]);
       }
     rnd(sum_idot, sum1);
     erg[0][j]=sum1/interval(j)/s[0][0];
   }

 for(int i=1; i<=erg.p; i++) // now calculating the remaining erg(i,k)
  {
    for(int k=0; k<=erg.p-i; k++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i-1; l++) // Koeff. (l,m)=(0,0) nicht summieren
	{
	 h=interval(n)*(interval(i)-interval(l))-interval(l);

	 for(int m=0; m<=k; m++)
	  {
	    accumulate(sum_idot, h*erg[l][m], s[i-l][k-m]);
	  }
	}
      rnd(sum_idot, sum1); 
        
      sum_idot=interval(0.0);
      for(int m=1; m<=k; m++) 
	{
	 accumulate(sum_idot, s[0][m], erg[i][k-m]);
	}
      rnd(sum_idot, sum2); 

      erg[i][k]=(sum1/interval(i)-sum2)/s[0][0];
     }
  } 
 
 return erg;
}

//------------------------------------------------------------------------

dim2taylor exp(const dim2taylor& s)
{
 dim2taylor erg(s.p);

 idotprecision sum_idot;
 interval sum;

 erg[0][0]=exp(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating erg(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot, h*erg[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     rnd(sum_idot, sum);  
	     erg[i][k-i]=sum/interval(k);
	   }
       }
   }

 return erg;
}

//------------------------------------------------------------------------

dim2taylor ln(const dim2taylor& s)
{
 dim2taylor f(s.p);

 if(0<=s[0][0])
   {
     cerr << "Error in dim2taylor ln : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, s, s, _ln);
 return f;
}

//------------------------------------------------------------------------

dim2taylor lnp1(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);
 g = interval(1) + s;

 if(0<=interval(1)+s[0][0])
   {
     cerr << "Error in dim2taylor lnp1 : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _lnp1);
 return f;
}
 
//------------------------------------------------------------------------

dim2taylor sqrtp1m1(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);
 g = 2*sqrt(1+s);;

 if(0<=interval(1)+s[0][0])
   {
     cerr << "Error in dim2taylor sqrtp1m1 : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _sqrtp1m1);
 return f;
}
 

//------------------------------------------------------------------------

dim2taylor sin(const dim2taylor& s)
{
 dim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 interval sum1, sum2;

 erg1[0][0]=sin(s[0][0]);
 erg2[0][0]=cos(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     rnd(sum_idot1, sum1);  
	     rnd(sum_idot2, sum2);

	     erg1[i][k-i]=sum1/interval(k);
	     erg2[i][k-i]=-sum2/interval(k);
	   }
       }
   }

 return erg1;
}


//------------------------------------------------------------------------

dim2taylor cos(const dim2taylor& s)
{
 dim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 interval sum1, sum2;

 erg1[0][0]=sin(s[0][0]);
 erg2[0][0]=cos(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     rnd(sum_idot1, sum1);  
	     rnd(sum_idot2, sum2);

	     erg1[i][k-i]=sum1/interval(k);
	     erg2[i][k-i]=-sum2/interval(k);
	   }
       }
   }

 return erg2;
}

//------------------------------------------------------------------------

dim2taylor tan(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=sqr(cos(s));

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor tan : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _tan);
 return f;
}


//------------------------------------------------------------------------

dim2taylor cot(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=-sqr(sin(s));

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor cot : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _cot);
 return f;
}

//------------------------------------------------------------------------

dim2taylor sinh(const dim2taylor& s)
{
 dim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 interval sum1, sum2;

 erg1[0][0]=sinh(s[0][0]);
 erg2[0][0]=cosh(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   }  // for m
	       }      // for l

	     rnd(sum_idot1, sum1);  
	     rnd(sum_idot2, sum2);

	     erg1[i][k-i]=sum1/interval(k);
	     erg2[i][k-i]=sum2/interval(k);
	   }
       }
   }

 return erg1;
}

//------------------------------------------------------------------------

dim2taylor cosh(const dim2taylor& s)
{
 dim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 interval sum1, sum2;

 erg1[0][0]=sinh(s[0][0]);
 erg2[0][0]=cosh(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     rnd(sum_idot1, sum1);  
	     rnd(sum_idot2, sum2);

	     erg1[i][k-i]=sum1/interval(k);
	     erg2[i][k-i]=sum2/interval(k);
	   }
       }
   }

 return erg2;
}

//------------------------------------------------------------------------

dim2taylor tanh(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=sqr(cosh(s));

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor tanh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _tanh);
 return f;
}

//------------------------------------------------------------------------

dim2taylor coth(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=-sqr(sinh(s));

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor coth : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _coth);
 return f;
}

//------------------------------------------------------------------------

dim2taylor asin(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=sqrt( interval(1.0)-sqr(s) );

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor asin : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _asin);
 return f;
}

//------------------------------------------------------------------------

dim2taylor acos(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=-sqrt( interval(1.0)-sqr(s) );

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor acos : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acos);
 return f;
}

//------------------------------------------------------------------------

dim2taylor atan(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=interval(1.0)+sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor atan : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _atan);
 return f;
}

//------------------------------------------------------------------------

dim2taylor acot(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=-(interval(1.0)+sqr(s));

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor acot : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acot);
 return f;
}

//------------------------------------------------------------------------

dim2taylor asinh(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

// g=sqrt(interval(1.0)+sqr(s));
 g = sqrt1px2(s); // Blomquist, 05.10.05;

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor asinh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _asinh);
 return f;
}

//------------------------------------------------------------------------

dim2taylor acosh(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

// g=sqrt(sqr(s)-interval(1.0));
 g = sqrtx2m1(s); // Blomquist, 05.10.05;

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor acosh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acosh);
 return f;
}

//------------------------------------------------------------------------

dim2taylor atanh(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=interval(1.0)-sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor atanh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _atanh);
 return f;
}

//------------------------------------------------------------------------

dim2taylor acoth(const dim2taylor& s)
{
 dim2taylor f(s.p), g(s.p);

 g=interval(1.0)-sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in dim2taylor acoth : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acoth);
 return f;
}

//------------------------------------------------------------------------

//Error function "erf" //added, mg2006-03
dim2taylor erf(const dim2taylor& s)
{
 dim2taylor erg(s.p);

 idotprecision sum_idot;
 interval sum;

 dim2taylor g(s.p);

 g=exp(-sqr(s));

 erg[0][0]=erf(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating erg(i,k-i)
                                     // equiv. to erg(k1,k2) in Braeuer Thesis
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot, h*g[l][m], s[i-l][k-i-m]);
                                             //Braeuer:s[k1-j1][k2-j2]
		   } // for m
	       }     // for l

	     rnd(sum_idot, sum);  
	     erg[i][k-i]=interval(2)*sum/(sqrt(Pi())*interval(k));
	   }
       }
   }

 return erg;

}

//------------------------------------------------------------------------

//Complementary Error function "erfc" //added, mg2006-03
dim2taylor erfc(const dim2taylor& s)
{
 dim2taylor erg(s.p);

 erg=interval(1)-erf(s);

 return erg;

}


//------------------------------------------------------------------------

void f_g_u(const dim2taylor& f, const dim2taylor& g, const dim2taylor& u,
		 int _fkt) 
{
 idotprecision sum_idot;

 interval sum1, sum2;
 
 switch(_fkt)   { 
   case _ln:       {f[0][0]=ln(u[0][0]);       break;}
   case _lnp1:     {f[0][0]=lnp1(u[0][0]);     break;}
   case _sqrtp1m1: {f[0][0]=sqrtp1m1(u[0][0]); break;}
   case _tan:      {f[0][0]=tan(u[0][0]);      break;}
   case _cot:      {f[0][0]=cot(u[0][0]);      break;}

   case _asin:  {f[0][0]=asin(u[0][0]);  break;}
   case _acos:  {f[0][0]=acos(u[0][0]);  break;}
   case _atan:  {f[0][0]=atan(u[0][0]);  break;}
   case _acot:  {f[0][0]=acot(u[0][0]);  break;}

   case _tanh:  {f[0][0]=tanh(u[0][0]);  break;}
   case _coth:  {f[0][0]=coth(u[0][0]);  break;}

   case _asinh:  {f[0][0]=asinh(u[0][0]);  break;}
   case _acosh:  {f[0][0]=acosh(u[0][0]);  break;}
   case _atanh:  {f[0][0]=atanh(u[0][0]);  break;}
   case _acoth:  {f[0][0]=acoth(u[0][0]);  break;}
   } 

 if(f.p>0)
   {
     for(int k=1; k<=f.p; k++)
       {
	 sum_idot=interval(0.0);

	 for(int j=1; j<=k-1; j++)
	   {
	     accumulate(sum_idot, interval(j)*f[0][j], g[0][k-j]);
	   }
	 rnd(sum_idot, sum1);
	 f[0][k]=(u[0][k]-sum1/interval(k))/g[0][0];
       }

     for(int i=1; i<=f.p; i++) // now calculating the remainding f(i,j)
       {
	 for(int j=0; j<=f.p-i; j++)
	   {
	     sum_idot=interval(0.0);
	     
	     for(int l=1; l<=i-1; l++) 
	       {
		 for(int m=0; m<=j; m++)
		  {
		  accumulate(sum_idot, interval(l)*f[l][m], g[i-l][j-m]);
		  }
	       }
	     rnd(sum_idot, sum1); 
	     
	     sum_idot=interval(0.0);
	     for(int m=1; m<=j; m++) 
	       {
		 accumulate(sum_idot, g[0][m], f[i][j-m]);
	       }
	     rnd(sum_idot, sum2); 
	     
	     f[i][j]=(u[i][j]-sum1/interval(i)-sum2)/g[0][0];
	   }
       } 
   }
}


//------------------------------------------------------------------------

// Class dim2taylor_vector

//------------------------------------------------------------------------
dim2taylor_vector::dim2taylor_vector() // Default constructor
{
    dim=2;
    lb=1;
    ub=2;
    p_el=1;
    comp = new dim2taylor[dim];   
    for (int i=0; i<dim; i++)
      comp[i]=dim2taylor(p_el);
    // dim constructors of the form
    // dim2taylor(int p_el)  are called to create a block (an array) of
    // dim elements of type  dim2taylor.
}


dim2taylor_vector::dim2taylor_vector(int ordnung, int Lb, int Ub)
{
    lb=Lb;
    ub=Ub;
    dim=ub-lb+1;
    p_el=ordnung;
    comp = new dim2taylor[dim];   
    for (int i=0; i<dim; i++)
      comp[i]=dim2taylor(p_el);
    // dim constructors of the form
    // dim2taylor(int p_el)  are called to create a block (an array) of
    // dim elements of type  dim2taylor.
}

//------------------------------------------------------------------------

dim2taylor_vector::dim2taylor_vector(const dim2taylor_vector& s)
{
 dim=s.dim;
 lb=s.lb;
 ub=s.ub;
 p_el=s.p_el;
 comp=new dim2taylor[dim]; 
          //Taylor order p_el and memory are handled by dim2taylor.operator=
 for(int i=0; i<dim; i++) comp[i]=s.comp[i];
}

//------------------------------------------------------------------------

dim2taylor_vector::~dim2taylor_vector()
{
 delete[] comp;
 comp=NULL;
}

//------------------------------------------------------------------------

dim2taylor_vector& dim2taylor_vector::operator=(const dim2taylor_vector& s)
{
 if(this!=&s)
 {
  delete[] comp;

  p_el=s.p_el;
  dim=s.dim;
  lb=s.lb;
  ub=s.ub;

 comp=new dim2taylor[dim]; 
          //Taylor order p_el and memory are handled by dim2taylor.operator=
  for(int i=0; i<dim ; i++) comp[i]=s.comp[i];
 }

 return *this;
}

//------------------------------------------------------------------------

dim2taylor& dim2taylor_vector::operator[](int n) const
{
    return comp[n-lb];  // allowed indices: n = lb,lb+1,...,ub;
} 

//------------------------------------------------------------------------

dim2taylor_vector init_var(int order, ivector& values)
{
 int dimension = Ub(values)-Lb(values)+1;

 if(dimension != 2) 
  {
   cerr << "Error in dim2taylor_vector::init_var" << endl;
   cerr << "! 2-dimensional Taylor arithmetic !" << endl; 
   exit(1);
  } 

 dim2taylor_vector erg(order, Lb(values), Ub(values));

 erg[Lb(values)] = init_var(order,1,values[Lb(values)]);
 erg[Ub(values)] = init_var(order,2,values[Ub(values)]);

 return erg;
}

//------------------------------------------------------------------------

int Lb(const dim2taylor_vector& d2tv)  //added, mg2005
{
  return d2tv.lb;
}

int Ub(const dim2taylor_vector& d2tv)  //added, mg2005
{
  return d2tv.ub;
}


//NEW: By T. Johnson
//------------------------------------------------------------------------

dim2taylor::dim2taylor(const dim2taylor& s, int n)
{
 p=s.p+1;
 dat=new ivector[p+1];
 for(int i=0; i<=p-1 ;i++) {
   dat[i]=s.dat[i];
   Resize(dat[i], 0, p-i);
   dat[i][p-i]=0.0;
 }
 Resize(dat[p],0,0);
 dat[p][0]=0.0;
 
}
//------------------------------------------------------------------------


} // end of namespace taylor


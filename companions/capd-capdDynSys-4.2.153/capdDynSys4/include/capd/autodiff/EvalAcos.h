/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalAcos.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ACOS_H_
#define _CAPD_AUTODIFF_EVAL_ACOS_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Acos ------------------------------------

namespace Acos
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned j=1;j<coeffNo;++j)
        temp += double(j) * result[j]*right[coeffNo-j];
      result[coeffNo] = -(left[coeffNo] + temp/(double)coeffNo)/(*right);
    }else{
      *result = acos(*left);
    }
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    // shift pointers to proper coefficients
    T* resultDer = result + order;
    const T* leftDer = left + order + coeffNo;

    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
    {
      T temp = *leftDer;
      for(unsigned j=0;j<coeffNo;++j)
        temp += resultDer[j] * right[coeffNo-j];
      resultDer[coeffNo] = -temp/(*right);
    }
  }

  template<class T>
  inline void evalCn(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    Multiindex a(dim),b(dim),c(dim);
    unsigned i,k,p;
    for(i=shift*order+coeffNo;i<end*order+coeffNo;i+=order)
      result[i] = TypeTraits<T>::zero();

    for(i=1;i<(degree+1)/2;++i){
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      const unsigned shiftB = binomial(dim+(degree-i)-1,dim);
      do{
        b.clear();
        b[0]=degree-i;
        const unsigned s1 = (shiftA + a.index(i))*order;
        do{
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          T r1 = TypeTraits<T>::zero(), r2 = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) {
            r1 += result[s1+k]*right[s2+coeffNo-k];
            r2 += result[s2+k]*right[s1+coeffNo-k];
          }
          result[(shift + c.index(degree))*order+coeffNo] += double(a[p])*r1 + double(b[p])*r2;
        }while(b.hasNext());
      }while(a.hasNext());
    }

    if(degree%2==0){
      i = degree/2;
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      do{
        b.clear();
        b[0]=i;
        const unsigned s1 = (shiftA + a.index(i))*order;
        do{
          const unsigned s2 = (shiftA + b.index(i))*order;
          if(s1<s2) continue;
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          T r1 = TypeTraits<T>::zero(), r2 = TypeTraits<T>::zero();
          if(s1!=s2){
            for(k=0;k<=coeffNo;++k) {
              r1 += result[s1+k]*right[s2+coeffNo-k];
              r2 += result[s2+k]*right[s1+coeffNo-k];
            }
            result[(shift + c.index(degree))*order+coeffNo] += double(a[p])*r1 + double(b[p])*r2;
          }
          else{
            for(k=0;k<=coeffNo;++k) r1 += result[s1+k]*right[s1+coeffNo-k];
            result[(shift + c.index(degree))*order+coeffNo] += double(a[p])*r1;
          }
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      for(k=0;k<coeffNo;++k)
        result[s+coeffNo] += c[p]*result[s+k]*right[coeffNo-k];

      result[s+coeffNo] = -(left[s+coeffNo] + result[s+coeffNo]/double(c[p]))/(*right);
    }while(c.hasNext());
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      default:
        evalC1(left,right,result,dim,order,coeffNo);
        for(unsigned i=2;i<=degree;++i)
          evalCn(i,left,right,result,dim,order,coeffNo);
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = acos(*left);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    T* resultDer = result + order;
    const T* leftDer = left + order;
    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
      *resultDer = -(*leftDer)/(*right);
  }

  template<class T>
  inline void evalCnHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    Multiindex a(dim),b(dim),c(dim);
    unsigned i,p;

    for(i=shift*order;i<end*order;i+=order)
      result[i] = TypeTraits<T>::zero();

    for(i=1;i<(degree+1)/2;++i){
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      const unsigned shiftB = binomial(dim+(degree-i)-1,dim);
      do{
        b.clear();
        b[0]=degree-i;
        const unsigned sa = (shiftA + a.index(i))*order;
        do{
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned sb = (shiftB + b.index(degree-i))*order;
          result[(shift + c.index(degree))*order] += a[p]*result[sa]*right[sb] + b[p]*result[sb]*right[sa];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    if(degree%2==0){
      i = degree/2;
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      do{
        b.clear();
        b[0]=i;
        const unsigned sa = (shiftA + a.index(i))*order;
        do{
          const unsigned sb = (shiftA + b.index(degree-i))*order;
          if(sa<sb) continue;
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          if(sa!=sb)
            result[(shift + c.index(degree))*order] += a[p]*result[sa]*right[sb] + b[p]*result[sb]*right[sa];
          else
            result[(shift + c.index(degree))*order] += a[p]*result[sa]*right[sa];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      result[s] = -(left[s] + result[s]/double(c[p]))/(*right);
    }while(c.hasNext());
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 0:
        evalC0HomogenousPolynomial(left,right,result);
        break;
      default:
        evalCnHomogenousPolynomial(degree,left,right,result,dim,order);
    }
  }
}

namespace AcosFunTime
{
  template<class T>
  inline void evalC0(const T* left, T* right, T* result, const unsigned coeffNo)
  {
    Acos::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    Acos::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = acos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = acos(*left);
  }
}

namespace AcosTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Acos::evalC0(left,right,result,coeffNo);
/*
    switch(coeffNo)
    {
      case 0:
        *result = asin(*left); break;
      case 1:
        result[1] = left[1]/(*right); break;
      default:
        const T tmp =
            double(coeffNo-1)*result[coeffNo-1]*right[1] +
            double(coeffNo-2)*result[coeffNo-2]*right[2];

        result[coeffNo] = (left[coeffNo] - tmp/double(coeffNo))/(*right);
    }
    */
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = acos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = acos(*left);
  }
}

namespace AcosConst
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      {}
    else
      *result = acos(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = acos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = acos(*left);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Acos);
CAPD_MAKE_CLASS_NODE(AcosConst);
CAPD_MAKE_CLASS_NODE(AcosTime);
CAPD_MAKE_CLASS_NODE(AcosFunTime);

}} // namespace capd::autodiff

#endif

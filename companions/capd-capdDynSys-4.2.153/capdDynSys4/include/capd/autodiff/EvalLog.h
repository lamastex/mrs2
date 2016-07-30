/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalLog.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_LOG_H_
#define _CAPD_AUTODIFF_EVAL_LOG_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Log ------------------------------------

namespace Log
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned j=1;j<coeffNo;++j)
        temp += T(j) * result[j] * left[coeffNo-j];
      temp /= (T)coeffNo;
      result[coeffNo] = (left[coeffNo] - temp)/(*left);
    }else
      *result = log(*left);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    // shift pointers to proper coefficients
    T* resultDer = result + order;
    const T* leftDer = left + order;

    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
    {
      T temp = leftDer[coeffNo];
      for(unsigned j=0;j<coeffNo;++j)
        temp -= resultDer[j] * left[coeffNo-j];
      resultDer[coeffNo] = temp/(*left);
    }
  }

  template<class T>
  inline void evalCn(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
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
        const unsigned sa = (shiftA + a.index(i))*order;
        do{
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned sb = (shiftB + b.index(degree-i))*order;
          T ra = TypeTraits<T>::zero(), rb = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) {
            ra += result[sa+k]*left[sb+coeffNo-k];
            rb += result[sb+k]*left[sa+coeffNo-k];
          }
          result[(shift + c.index(degree))*order+coeffNo] -= a[p]*ra + b[p]*rb;
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
          T ra = TypeTraits<T>::zero(), rb = TypeTraits<T>::zero();
          if(sa!=sb){
            for(k=0;k<=coeffNo;++k) {
              ra += result[sa+k]*left[sb+coeffNo-k];
              rb += result[sb+k]*left[sa+coeffNo-k];
            }
            result[(shift + c.index(degree))*order+coeffNo] -= a[p]*ra + b[p]*rb;
          }
          else{
            for(k=0;k<=coeffNo;++k) ra += result[sa+k]*left[sa+coeffNo-k];
            result[(shift + c.index(degree))*order+coeffNo] -= a[p]*ra;
          }
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      T r = TypeTraits<T>::zero();
      for(k=0;k<coeffNo;++k) {
        r += left[coeffNo-k]*result[s+k];
      }
      result[s+coeffNo] -= c[p]*r;
      result[s+coeffNo] /= typename TypeTraits<T>::Real(c[p]);
      result[s+coeffNo] += left[s+coeffNo];
      result[s+coeffNo] /= (*left);
    }while(c.hasNext());
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
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
    *result = log(*left);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    T* resultDer = result + order;
    const T* leftDer = left + order;

    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
      *resultDer = (*leftDer) / (*left);
  }

  template<class T>
  inline void evalCnHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
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
          result[(shift + c.index(degree))*order] -= a[p]*result[sa]*left[sb] + b[p]*result[sb]*left[sa];
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
            result[(shift + c.index(degree))*order] -= a[p]*result[sa]*left[sb] + b[p]*result[sb]*left[sa];
          else
            result[(shift + c.index(degree))*order] -= a[p]*result[sa]*left[sa];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      result[s] /= typename TypeTraits<T>::Real(c[p]);
      result[s] += left[s];
      result[s] /= (*left);
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

namespace LogFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Log::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    Log::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = log(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = log(*left);
  }
}

namespace LogTime
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T temp = T(coeffNo-1) * result[coeffNo-1];
      temp /= (T)coeffNo;
      result[coeffNo] = (left[coeffNo] - temp)/(*left);
    }else
      *result = log(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = log(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = log(*left);
  }
}


namespace LogConst
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = log(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = log(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = log(*left);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Log);
CAPD_MAKE_CLASS_NODE(LogConst);
CAPD_MAKE_CLASS_NODE(LogTime);
CAPD_MAKE_CLASS_NODE(LogFunTime);

}} // namespace capd::autodiff

#endif

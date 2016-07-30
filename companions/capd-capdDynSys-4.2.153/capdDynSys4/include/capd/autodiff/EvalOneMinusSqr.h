/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalOneMinusSqr.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ONE_MINUS_SQR_H_
#define _CAPD_AUTODIFF_EVAL_ONE_MINUS_SQR_H_

#include "capd/autodiff/EvalSqr.h"

namespace capd{
namespace autodiff{

// -------------------- OneMinusSqr ------------------------------------

namespace OneMinusSqr{

  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = -Sqr::sqrProduct(left,coeffNo);
    if(coeffNo==0) *result += TypeTraits<T>::one();
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const T* leftDer = left+order;
    result += order+coeffNo;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,result+=order)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned i=0;i<=coeffNo;++i)
        temp += left[i] * leftDer[coeffNo-i];
      *result = -2.*temp;
    }
  }


  template<class T>
  inline void evalC2(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left+order;
    T* resultDer = result + order + coeffNo;
    // begin of C^2
    const T* leftHess = leftDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      // case dx^2 and dx
      T temp = TypeTraits<T>::zero();
      T temp2 = TypeTraits<T>::zero();
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
      {
        temp += left[i]*leftHess[j];
        temp2 += left[i]*leftDer[j];
      }
      *resultDer = -2.*temp2;
      *resultHess = -2.*temp - Sqr::sqrProduct(leftDer,coeffNo);

      leftHess += order;
      resultHess += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,leftHess+=order,resultHess+=order)
      {
        temp = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          temp += left[i] * leftHess[j];
          temp += leftDer[i] * leftDer2[j];
        }
        *resultHess = -2.*temp;
      }
    }
  }


  template<class T>
  inline void evalC3(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC2(left,right,result,dim,order,coeffNo);
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      T temp = TypeTraits<T>::zero();
      // case dxdxdx
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
      {
        temp += left[i]* left[i111+j];
        temp += left[i1+i] * left[i11+j];
      }
      result[i111+coeffNo] = -2.*temp;

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        temp = TypeTraits<T>::zero();
        T temp2 = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          temp += left[i]*left[i112+j];    // 0,xxy
          temp += left[i1+i]*left[i12+j];  // x,xy
          temp += left[i2+i]*left[i11+j];  // y,xx

          temp2 += left[i]*left[i122+j];   // 0,xyy
          temp2 += left[i1+i]*left[i22+j]; // x,yy
          temp2 += left[i2+i]*left[i12+j]; // y,xy
        }
        result[i112+coeffNo] = -2.*temp;
        result[i122+coeffNo] = -2.*temp2;

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            temp += left[i] * left[i123+j];    // 0,xyz
            temp += left[i1+i] * left[i23+j];  // x,yz
            temp += left[i2+i] * left[i13+j];  // y,xz
            temp += left[i3+i] * left[i12+j];  // z,xy
          }
          result[i123+coeffNo] = -2.*temp;
        }
      }
    }
  }

  template<class T>
  void evalCn(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
    {
      T r = TypeTraits<T>::zero();
      for(k=0;k<=coeffNo;++k)
        r -= left[k]*left[i+coeffNo-k];
      result[i+coeffNo] = r;
    }

    Multiindex a(dim),b(dim),c(dim);
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
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          T r = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) r += left[s1+k]*left[s2+coeffNo-k];
          result[(shift + c.index(degree))*order + coeffNo] -= r;
        }while(b.hasNext());
      }while(a.hasNext());
    }

    if(degree%2==0)
    {
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
          if(s1>s2) continue;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          if(s1!=s2){
            T r = TypeTraits<T>::zero();
            for(k=0;k<=coeffNo;++k) r += left[s1+k]*left[s2+coeffNo-k];
            result[(shift + c.index(degree))*order + coeffNo] -= r;
          }
          else
            result[(shift + c.index(degree))*order + coeffNo] -= 0.5*Sqr::sqrProduct(left+s1,coeffNo);
        }while(b.hasNext());
      }while(a.hasNext());
    }

    for(i=shift*order+coeffNo;i<end*order+coeffNo;i+=order)
      result[i] *= typename TypeTraits<T>::Real(2.);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 2:
        evalC2(left,right,result,dim,order,coeffNo);
        break;
      case 1:
        evalC1(left,right,result,dim,order,coeffNo);
        break;
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      default:
        evalC3(left,right,result,dim,order,coeffNo);
        for(unsigned i=4;i<=degree;++i)
          evalCn(i, left,right, result, dim, order, coeffNo);
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    (*result) = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* /*right*/, T* result, unsigned dim, unsigned order)
  {
    const T* leftDer = left + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
      *resultDer = -2.*(*left)*(*leftDer);
  }

  template<class T>
  inline void evalC2HomogenousPolynomial(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    const  T* leftDer = left + order;
    T* resultDer = result + order;
    const unsigned s = dim*order;
    const T* leftHess = leftDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      // case dxdx
      *resultHess  = -2.*(*left)*(*leftHess) - sqr(*leftDer);

      // case dxdy
      resultHess +=order;
      leftHess += order;
      const T* leftDer2 = leftDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,resultHess+=order,leftHess+=order,leftDer2+=order)
        *resultHess  = -2.*( (*left)*(*leftHess)  +  (*leftDer)*(*leftDer2));
    }
  }

  template<class T>
  inline void evalC3HomogenousPolynomial(const T* left, const T* /*right*/, T* result, unsigned dim, unsigned order)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;

      // case dxdxdx
      result[i111]= -2.*(*left* left[i111] + left[i1]*left[i11]);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111+order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        result[i112] = -2.*(*left*left[i112] + left[i1]*left[i12] + left[i2]*left[i11]);
        result[i122] = -2.*(*left*left[i122] + left[i1]*left[i22] + left[i2]*left[i12]);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
          result[i123] = -2.*(*left*left[i123] + left[i1]*left[i23] + left[i2]*left[i13] + left[i3]*left[i12]);
      }
    }
  }  // evalC3

  template<class T>
  void evalCnHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
      result[i] = -left[0]*left[i];

    Multiindex a(dim),b(dim),c(dim);
    for(i=1;i<(degree+1)/2;++i){
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      const unsigned shiftB = binomial(dim+(degree-i)-1,dim);
      do{
        b.clear();
        b[0]=degree-i;
        do{
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          result[(shift + c.index(degree))*order] -= left[(shiftA + a.index(i))*order]*left[(shiftB + b.index(degree-i))*order];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    if(degree%2==0)
    {
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
          if(s1>s2) continue;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          if(s1!=s2)
            result[(shift + c.index(degree))*order] -= left[s1]*left[s2];
          else
            result[(shift + c.index(degree))*order] -= 0.5*sqr(left[s1]);
        }while(b.hasNext());
      }while(a.hasNext());
    }

    for(i=shift*order;i<end*order;i+=order)
      result[i] *= typename TypeTraits<T>::Real(2.);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    switch(degree)
    {
      case 3:
        evalC3HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 2:
        evalC2HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 1:
        evalC1HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 0:
        evalC0HomogenousPolynomial(left,right,result);
        break;
      default:
        evalCnHomogenousPolynomial(degree, left,right, result, dim, order);
    }
  }

}

// -------------------- OneMinusSqrFunTime ------------------------------------

namespace OneMinusSqrFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, unsigned coeffNo)
  {
    OneMinusSqr::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    OneMinusSqr::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = TypeTraits<T>::one()-sqr(*left);
  }

}

// -------------------- SqrTime ------------------------------------

namespace OneMinusSqrTime
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 2: result[2] = -1.; break;
      case 1: result[1] = -2.*(*left); break;
      case 0: *result = TypeTraits<T>::one()-sqr(*left);
    }
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = TypeTraits<T>::one()-sqr(*left);
  }
}

// -------------------- SqrConst ------------------------------------

namespace OneMinusSqrConst
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    if(degree==0 and coeffNo==0)
      *result = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = TypeTraits<T>::one()-sqr(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = TypeTraits<T>::one()-sqr(*left);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(OneMinusSqr);
CAPD_MAKE_CLASS_NODE(OneMinusSqrTime);
CAPD_MAKE_CLASS_NODE(OneMinusSqrFunTime);
CAPD_MAKE_CLASS_NODE(OneMinusSqrConst);

}} // namespace capd::autodiff

#endif

/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalSqrt.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SQRT_H_
#define _CAPD_AUTODIFF_EVAL_SQRT_H_
#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Sqrt ------------------------------------

namespace Sqrt
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      const int n = coeffNo;
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<n;++j)
        temp += (n-3*j) * result[j] * left[n-j];
      result[n] = 0.5*temp/((double)n * (*left));
    }else
      *result = sqrt(*left);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    const T* leftDer = left + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      T temp1 = result[coeffNo] * (*leftDer);
      T temp2 = capd::TypeTraits<T>::zero();
      for(unsigned j=0;j<coeffNo;++j)
      {
        temp1 += result[j] * leftDer[coeffNo-j];
        temp2 += resultDer[j] * left[coeffNo-j];
      }
      resultDer[coeffNo] = (temp1*0.5-temp2)/(*left);
    }
  } // evalC1

  template<class T>
  inline void evalC2(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,dim,order,coeffNo);

    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left+order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = leftDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      // case dx^2
      T temp1 = result[coeffNo] * (*leftHess);
      T temp2 = resultDer[coeffNo] * (*leftDer);
      T temp3 = capd::TypeTraits<T>::zero();
      for(unsigned j=0;j<coeffNo;++j)
      {
        temp1 += result[j] * leftHess[coeffNo-j];
        temp2 += resultDer[j] * leftDer[coeffNo-j];
        temp3 += resultHess[j] * left[coeffNo-j];
      }
      resultHess[coeffNo] = (0.5*temp1-0.25*temp2-temp3)/(*left);

      leftHess += order;
      resultHess += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      T* resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,resultDer2+=order,leftHess+=order,resultHess+=order)
      {
        temp1 = result[coeffNo] * (*leftHess) + resultDer[coeffNo] *(*leftDer2);
        temp2 = resultDer2[coeffNo] *(*leftDer);
        for(unsigned j=0;j<coeffNo;++j)
        {
          temp1 += result[j] * leftHess[coeffNo-j] + resultDer[j] * leftDer2[coeffNo-j];
          temp2 += resultDer2[j] * leftDer[coeffNo-j] + resultHess[j]*left[coeffNo-j];
        }
        resultHess[coeffNo] = (0.5*temp1-temp2)/(*left);
      }
    }
  } // evalC2

  template<class T>
  inline void evalC3(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC2(left,right,result,dim,order,coeffNo);

    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      T temp1 = result[coeffNo]*left[i111] - result[i11+coeffNo]*left[i1];
      T temp2 = capd::TypeTraits<T>::zero();
      // case dxdxdx
      for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
      {
        temp1 += result[i]*left[i111+j] - result[i11+i]*left[i1+j];
        temp2 += result[i111+i]*left[j];
      }
      result[i111+coeffNo] = (0.5*temp1-temp2)/(*left);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        temp1 = result[coeffNo]* left[i112] + result[i1+coeffNo]*left[i12] + result[i11+coeffNo]*left[i2];
        temp2 = result[i12+coeffNo]  * left[i1] + result[i2+coeffNo]* left[i11];
        T temp3 = result[coeffNo]* left[i122] + result[i2+coeffNo]*left[i12] + result[i22+coeffNo]*left[i1];
        T temp4 = result[i12+coeffNo]  * left[i2] + result[i1+coeffNo]* left[i22];
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        {
          temp1 += result[i]* left[i112+j] + result[i1+i]*left[i12+j] + result[i11+i]*left[i2+j];
          temp2 += result[i12+i]  * left[i1+j] + result[i2+i]* left[i11+j] + result[i112+i]*left[j];
          temp3 += result[i]* left[i122+j] + result[i2+i]*left[i12+j] + result[i22+i]*left[i1+j];
          temp4 += result[i12+i]  * left[i2+j] + result[i1+i]* left[i22+j] + result[i122+i]*left[j];
        }
        result[i112+coeffNo] = (0.5*temp1-temp2)/(*left);
        result[i122+coeffNo] = (0.5*temp3-temp4)/(*left);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          temp1 = result[i12+coeffNo]*left[i3] + result[i2+coeffNo]*left[i13] + result[i1+coeffNo]*left[i23] + result[coeffNo]*left[i123];
          temp2 = result[i13+coeffNo]*left[i2] + result[i23+coeffNo]*left[i1] + result[i3+coeffNo]*left[i12];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp1 += result[i12+i]*left[i3+j] + result[i2+i]*left[i13+j] + result[i1+i]*left[i23+j] + result[i]*left[i123+j];
            temp2 += result[i13+i]*left[i2+j] + result[i23+i]*left[i1+j] + result[i3+i]*left[i12+j] + result[i123+i]*left[j];
          }
          result[i123+coeffNo] = (0.5*temp1-temp2)/(*left);
        }
      }
    }
  } // evalC3

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
        const unsigned s1 = (shiftA + a.index(i))*order;
        do{
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          T r1 = TypeTraits<T>::zero(), r2 = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) {
            r1 += result[s1+k]*left[s2+coeffNo-k];
            r2 += result[s2+k]*left[s1+coeffNo-k];
          }
          result[(shift + c.index(degree))*order+coeffNo] += ((0.5*b[p]-a[p])*r1 + (0.5*a[p]-b[p])*r2);
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
              r1 += result[s1+k]*left[s2+coeffNo-k];
              r2 += result[s2+k]*left[s1+coeffNo-k];
            }
            result[(shift + c.index(degree))*order+coeffNo] += ((0.5*b[p]-a[p])*r1 + (0.5*a[p]-b[p])*r2);
          }
          else{
            for(k=0;k<=coeffNo;++k) r1 += result[s1+k]*left[s1+coeffNo-k];
            result[(shift + c.index(degree))*order+coeffNo] -= (0.5*a[p])*r1;
          }
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      T r1 = result[coeffNo]*left[s];
      T r2 = TypeTraits<T>::zero();
      for(k=0;k<coeffNo;++k) {
        r1 += result[k]*left[s+coeffNo-k];
        r2 += result[s+k]*left[coeffNo-k];
      }
      result[s+coeffNo] += c[p]*(0.5*r1 - r2);
      result[s+coeffNo] /= (*left * c[p]);
    }while(c.hasNext());
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
          evalCn(i,left,right,result,dim,order,coeffNo);
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left+order;
    T* resultDer = result+order;
    const T temp = 0.5*(*result)/(*left);
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
      *resultDer =  temp * (*leftDer);
  }

  template<class T>
  inline void evalC2HomogenousPolynomial(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left+order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = leftDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      // case dx^2
      *resultHess = (0.5*(*result) * (*leftHess)-0.25*(*resultDer) * (*leftDer))/(*left);

      leftHess += order;
      resultHess += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      T* resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,resultDer2+=order,leftHess+=order,resultHess+=order)
        *resultHess = (0.5*((*result) * (*leftHess) + (*resultDer) *(*leftDer2)) - (*resultDer2) *(*leftDer))/(*left);
    }
  } // evalC2

  template<class T>
  inline void evalC3HomogenousPolynomial(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      result[i111] = 0.5*((*result)*left[i111] - result[i11]*left[i1])/(*left);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        T temp1 = (*result)* left[i112] + result[i1]*left[i12] + result[i11]*left[i2];
        T temp2 = result[i12]  * left[i1] + result[i2]* left[i11];
        T temp3 = (*result)* left[i122] + result[i2]*left[i12] + result[i22]*left[i1];
        T temp4 = result[i12]  * left[i2] + result[i1]* left[i22];

        result[i112] = (0.5*temp1-temp2)/(*left);
        result[i122] = (0.5*temp3-temp4)/(*left);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          temp1 = result[i12]*left[i3] + result[i2]*left[i13] + result[i1]*left[i23] + (*result)*left[i123];
          temp2 = result[i13]*left[i2] + result[i23]*left[i1] + result[i3]*left[i12];
          result[i123] = (0.5*temp1-temp2)/(*left);
        }
      }
    }
  } // evalC3

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
        const unsigned s1 = (shiftA + a.index(i))*order;
        do{
          p = sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          result[(shift + c.index(degree))*order] += ((0.5*b[p]-a[p])*result[s1]*left[s2] + (0.5*a[p]-b[p])*result[s2]*left[s1]);
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
          //T r1 = TypeTraits<T>::zero(), r2 = TypeTraits<T>::zero();
          if(s1!=s2)
            result[(shift + c.index(degree))*order] += ((0.5*b[p]-a[p])*result[s1]*left[s2] + (0.5*a[p]-b[p])*result[s2]*left[s1]);
          else
            result[(shift + c.index(degree))*order] -= (0.5*a[p])*result[s1]*left[s1];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      result[s] += (c[p]*0.5)* (*result*left[s]);
      result[s] /= (*left * c[p]);
    }while(c.hasNext());
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
        evalCnHomogenousPolynomial(degree,left,right,result,dim,order);
    }
  }

}

// -------------------- SqrtFunTime ------------------------------------

namespace SqrtFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Sqrt::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    Sqrt::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = sqrt(*left);
  }
}

// -------------------- SqrtTime ------------------------------------

namespace SqrtTime
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = (1.5-coeffNo) * result[coeffNo-1]/((double)coeffNo * (*left));
    else
      *result = sqrt(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = sqrt(*left);
  }

}

// -------------------- SqrtConst ------------------------------------

namespace SqrtConst
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned /*coeffNo*/)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned /*coeffNo*/)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = sqrt(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = sqrt(*left);
  }

}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Sqrt);
CAPD_MAKE_CLASS_NODE(SqrtConst);
CAPD_MAKE_CLASS_NODE(SqrtTime);
CAPD_MAKE_CLASS_NODE(SqrtFunTime);

}} // namespace capd::autodiff

#endif

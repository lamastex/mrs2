/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalSinCos.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SIN_COS_H_
#define _CAPD_AUTODIFF_EVAL_SIN_COS_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Sin ------------------------------------

namespace Sin
{
  template<class T>
  inline void evalC0(const T* left, T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T tempSin = capd::TypeTraits<T>::zero();
      T tempCos = capd::TypeTraits<T>::zero();
      for(unsigned j=1;j<=coeffNo;++j)
      {
        tempSin += double(j) * right[coeffNo-j] * left[j];
        tempCos -= double(j) * result[coeffNo-j] * left[j];
      }
      result[coeffNo] = tempSin/(double)coeffNo;
      right[coeffNo] = tempCos/(double)coeffNo;
    }
    else
    {
      (*result) = sin(*left);
      (*right) = cos(*left);
    }
  }


  template<class T>
  inline void evalC1(const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    const T* leftDer = left+order;
    T* rightDer = right+order+coeffNo;
    T* resultDer = result+order+coeffNo;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      T tempSin = capd::TypeTraits<T>::zero();
      T tempCos = capd::TypeTraits<T>::zero();
      for(unsigned j=0;j<=coeffNo;++j)
      {
        tempSin += leftDer[coeffNo-j]*right[j];
        tempCos -= leftDer[coeffNo-j]*result[j];
      }
      *resultDer = tempSin;
      *rightDer = tempCos;
    }
  } // evalC1

  template<class T>
  inline void evalC2(const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,dim,order,coeffNo);

    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left + order;
    T* rightDer = right + order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = leftDer + s;
    T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      // case dx^2
      T tempSin1 = TypeTraits<T>::zero();
      T tempSin2 = TypeTraits<T>::zero();
      T tempCos1 = TypeTraits<T>::zero();
      T tempCos2 = TypeTraits<T>::zero();
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
      {
        tempSin1 += leftDer[i]*rightDer[j];
        tempSin2 += leftHess[i]*right[j];
        tempCos1 -= leftDer[i]*resultDer[j];
        tempCos2 -= leftHess[i]*result[j];
      }
      resultHess[coeffNo] = 0.5*tempSin1 + tempSin2;
      rightHess[coeffNo]  = 0.5*tempCos1 + tempCos2;

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      T* rightDer2 = rightDer + order;
      T* resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,rightDer2+=order,resultDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        tempSin1 = TypeTraits<T>::zero();
        tempCos1 = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          tempSin1 += leftDer[i] * rightDer2[j];
          tempSin1 += leftHess[i] * right[j];
          tempCos1 -= leftDer[i] * resultDer2[j];
          tempCos1 -= leftHess[i] * result[j];
        }
        resultHess[coeffNo] = tempSin1;
        rightHess[coeffNo] = tempCos1;
      }
    }
  }  // evalC2

  /// hand optimized code for third order jet propagation of sine
  template<class T>
  inline void evalC3(const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC2(left,right,result,dim,order,coeffNo);

    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      T tempSin = TypeTraits<T>::zero();
      T tempCos = TypeTraits<T>::zero();
      // case dxdxdx
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
      {
        tempSin += left[i1+i] * right[i11+j];
        tempSin += 2.*left[i11+i] * right[i1+j];
        tempSin += 3.*left[i111+i]* right[j];

        tempCos -= left[i1+i] * result[i11+j];
        tempCos -= 2.*left[i11+i] * result[i1+j];
        tempCos -= 3.*left[i111+i]* result[j];
      }
      result[i111+coeffNo] = tempSin/3.;
      right[i111+coeffNo] = tempCos/3.;

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        tempSin = TypeTraits<T>::zero();
        tempCos = TypeTraits<T>::zero();
        T tempSin2 = TypeTraits<T>::zero();
        T tempCos2 = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          tempSin += left[i1+i]*right[i22+j];  // x,yy
          tempSin += left[i12+i]*right[i2+j];  // xy,y
          tempSin += left[i122+i]*right[j];    // xyy,0

          tempSin2 += left[i2+i]*right[i11+j]; // y,xx
          tempSin2 += left[i12+i]*right[i1+j]; // xy,x
          tempSin2 += left[i112+i]*right[j];   // xxy,0

          tempCos -= left[i1+i]*result[i22+j];  // x,yy
          tempCos -= left[i12+i]*result[i2+j];  // xy,y
          tempCos -= left[i122+i]*result[j];    // xxy,0

          tempCos2 -= left[i2+i]*result[i11+j]; // y,xx
          tempCos2 -= left[i12+i]*result[i1+j]; // xy,x
          tempCos2 -= left[i112+i]*result[j];   // xxy,0
        }
        result[i112+coeffNo] = tempSin2;
        result[i122+coeffNo] = tempSin;
        right[i112+coeffNo] = tempCos2;
        right[i122+coeffNo] = tempCos;

        // case dxdydz, assume x<y<z
        unsigned i3 = i2 + order;
        unsigned i123 = i122 + order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          tempSin = TypeTraits<T>::zero();
          tempCos = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            tempSin += left[i1+i]   * right[i23+j];  // x,yz
            tempSin += left[i12+i]  * right[i3+j];  // xy,z
            tempSin += left[i13+i]  * right[i2+j];  // xz,y
            tempSin += left[i123+i] * right[j];      // xyz,0

            tempCos -= left[i1+i]   * result[i23+j];  // x,yz
            tempCos -= left[i12+i]  * result[i3+j];  // xy,z
            tempCos -= left[i13+i]  * result[i2+j];  // xz,y
            tempCos -= left[i123+i] * result[j];      // xyz,0
          }
          result[i123+coeffNo] = tempSin;
          right[i123+coeffNo] = tempCos;
        }
      }
    }
  }  // evalC3

  template<class T>
  inline void evalCn(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    Multiindex a(dim),b(dim),c(dim);
    unsigned i,k,p;
    for(i=shift*order+coeffNo;i<end*order+coeffNo;i+=order)
      right[i] = result[i] = TypeTraits<T>::zero();


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
          T sin1 = TypeTraits<T>::zero(), sin2 = TypeTraits<T>::zero();
          T cos1 = TypeTraits<T>::zero(), cos2 = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) {
            sin1 += right[s1+k]*left[s2+coeffNo-k];
            cos1 += result[s1+k]*left[s2+coeffNo-k];
            sin2 += right[s2+k]*left[s1+coeffNo-k];
            cos2 += result[s2+k]*left[s1+coeffNo-k];
          }
          const unsigned s3 = (shift + c.index(degree))*order+coeffNo;
          result[s3] += b[p]*sin1 + a[p]*sin2;
          right[s3] -= b[p]*cos1 + a[p]*cos2;
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
          const unsigned s3 = (shift + c.index(degree))*order+coeffNo;
          T sin1 = TypeTraits<T>::zero(), sin2 = TypeTraits<T>::zero();
          T cos1 = TypeTraits<T>::zero(), cos2 = TypeTraits<T>::zero();
          if(s1!=s2){
            for(k=0;k<=coeffNo;++k) {
              sin1 += right[s1+k]*left[s2+coeffNo-k];
              cos1 += result[s1+k]*left[s2+coeffNo-k];
              sin2 += right[s2+k]*left[s1+coeffNo-k];
              cos2 += result[s2+k]*left[s1+coeffNo-k];
            }
            result[s3] += b[p]*sin1 + a[p]*sin2;
            right[s3] -= b[p]*cos1 + a[p]*cos2;
          }
          else{
            for(k=0;k<=coeffNo;++k){
              sin1 += right[s1+k]*left[s1+coeffNo-k];
              cos1 += result[s1+k]*left[s1+coeffNo-k];
            }
            result[s3] += a[p]*sin1;
            right[s3] -= a[p]*cos1;
          }
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;

      T sin1 = right[coeffNo]*left[s];
      T cos1 = result[coeffNo]*left[s];
      for(k=0;k<coeffNo;++k) {
        sin1 += right[k]*left[s+coeffNo-k];
        cos1 += result[k]*left[s+coeffNo-k];
      }
      result[s+coeffNo] /= typename TypeTraits<T>::Real(c[p]);
      right[s+coeffNo]  /= typename TypeTraits<T>::Real(c[p]);

      result[s+coeffNo] += sin1;
      right[s+coeffNo]  -= cos1;
    }while(c.hasNext());
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
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
  inline void evalC0HomogenousPolynomial(const T* left, T* right, T* result)
  {
    (*result) = sin(*left);
    (*right) = cos(*left);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, T* right, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left + order;
    T* rightDer = right + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      *resultDer = (*leftDer)*(*right);
      *rightDer = -(*leftDer)*(*result);
    }
  }

  template<class T>
  inline void evalC2HomogenousPolynomial(const T* left, T* right, T* result, const unsigned dim, const unsigned order)
  {
    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left + order;
    T* rightDer = right + order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = leftDer + s;
    T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      // case dx^2
      *resultHess = 0.5*(*leftDer)*(*rightDer) + (*leftHess)*(*right);
      *rightHess  = -0.5*(*leftDer)*(*resultDer) - (*leftHess)*(*result);

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      T* rightDer2 = rightDer + order;
      T* resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,rightDer2+=order,resultDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        *resultHess = (*leftDer) * (*rightDer2) + (*leftHess) * (*right);
        *rightHess = -(*leftDer) * (*resultDer2) - (*leftHess) * (*result);
      }
    }
  }

  template<class T>
  inline void evalC3HomogenousPolynomial(const T* left, T* right, T* result, unsigned dim, unsigned order)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      // case dxdxdx
      result[i111] = (left[i1] * right[i11] + 2.*left[i11] * right[i1])/3. + left[i111]* (*right);
      right[i111] = -(left[i1] * result[i11] + 2.*left[i11] * result[i1])/3. - left[i111]* (*result);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;

        result[i112] = left[i2]*right[i11] + left[i12]*right[i1] + left[i112]*(*right);
        result[i122] = left[i1]*right[i22] + left[i12]*right[i2] + left[i122]*(*right);
        right[i112] = -(left[i2]*result[i11] + left[i12]*result[i1] + left[i112]*(*result));
        right[i122] = -(left[i1]*result[i22] + left[i12]*result[i2] + left[i122]*(*result));

        // case dxdydz, assume x<y<z
        unsigned i3 = i2 + order;
        unsigned i123 = i122 + order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          result[i123] = left[i1]*right[i23] + left[i12]*right[i3] + left[i13]*right[i2] + left[i123]*(*right);
          right[i123] = -(left[i1]*result[i23] + left[i12]*result[i3] + left[i13]*result[i2] + left[i123]*(*result));
        }
      }
    }
  }

  template<class T>
  inline void evalCnHomogenousPolynomial(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    Multiindex a(dim),b(dim),c(dim);
    unsigned i,p;//k,
    for(i=shift*order;i<end*order;i+=order)
      right[i] = result[i] = TypeTraits<T>::zero();


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
          const unsigned s3 = (shift + c.index(degree))*order;
          result[s3] += b[p]*right[s1]*left[s2] + a[p]*right[s2]*left[s1];
          right[s3] -= b[p]*result[s1]*left[s2] + a[p]*result[s2]*left[s1];
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
          const unsigned s3 = (shift + c.index(degree))*order;
          if(s1!=s2){
            result[s3] += b[p]*right[s1]*left[s2] + a[p]*right[s2]*left[s1];
            right[s3] -= b[p]*result[s1]*left[s2] + a[p]*result[s2]*left[s1];
          }
          else{
            result[s3] += b[p]*right[s1]*left[s1];
            right[s3] -= b[p]*result[s1]*left[s1];
          }
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = findMax(c.begin(),dim);
      const unsigned s = (shift + c.index(degree))*order;
      result[s] /= typename TypeTraits<T>::Real(c[p]);
      right[s]  /= typename TypeTraits<T>::Real(c[p]);

      result[s] += (*right)*left[s];
      right[s]  -= (*result)*left[s];
    }while(c.hasNext());
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order)
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

namespace SinConst{

  template<class T>
  inline void evalC0(const T* left, T* right, T* result, const unsigned coeffNo)
  {
    if(!coeffNo)
    {
      *result = sin(*left);
      *right = cos(*left);
    }
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, T* result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
    {
      *result = sin(*left);
      *right = cos(*left);
    }
  }
}

namespace SinTime{

  template<class T>
  inline void evalC0(const T* left, T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      result[coeffNo] = right[coeffNo-1]/(double)coeffNo;
      right[coeffNo] = -result[coeffNo-1]/(double)coeffNo;
    }
    else
    {
      (*result) = sin(*left);
      (*right) = cos(*left);
    }
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, T* result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order)
  {
    SinConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

namespace SinFunTime{

  template<class T>
  inline void evalC0(const T* left, T* right, T* result, const unsigned coeffNo)
  {
    Sin::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    Sin::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, T* result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, T* result, const unsigned dim, const unsigned order)
  {
    SinConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Sin);
CAPD_MAKE_CLASS_NODE(SinConst);
CAPD_MAKE_CLASS_NODE(SinTime);
CAPD_MAKE_CLASS_NODE(SinFunTime);

}} // namespace capd::autodiff

#endif

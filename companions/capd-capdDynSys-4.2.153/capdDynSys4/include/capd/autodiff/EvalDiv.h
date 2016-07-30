/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalDiv.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_DIV_H_
#define _CAPD_AUTODIFF_EVAL_DIV_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Div  -------------------------------

namespace Div
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    T temp = left[coeffNo];
    for(unsigned j=0;j<coeffNo;++j)
      temp -= result[j] * right[coeffNo-j];
    result[coeffNo] = temp/(*right);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const T* leftDer = left + order + coeffNo;
    const T* rightDer = right + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      T temp = (*leftDer) - result[coeffNo] * (*rightDer);
      for(unsigned j=0;j<coeffNo;++j)
        temp -= ( result[j]* rightDer[coeffNo-j] + resultDer[j] * right[coeffNo-j] );
      resultDer[coeffNo] = temp/(*right);
    }
  }

  /// hand optimized code for second order jet propagation of division
  template<class T>
  inline void evalC2(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,dim,order,coeffNo);
    const unsigned s = dim*order;
    // begin of C^1
    const T* rightDer = right+order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = left + order + s + coeffNo;
    const T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,rightDer+=order,resultDer+=order)
    {
      // case dxdx
      T temp = (*leftHess) - result[coeffNo]*(*rightHess) - resultDer[coeffNo]*(*rightDer);
      for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        temp -= ( result[i]* rightHess[j] + resultDer[i]*rightDer[j] + resultHess[i] * right[j] );
      resultHess[coeffNo] = temp/(*right);

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      T* resultDer2 = resultDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,resultDer2+=order,rightDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        temp = (*leftHess) - result[coeffNo]*(*rightHess) - resultDer2[coeffNo]*(*rightDer) - resultDer[coeffNo]*(*rightDer2);
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        {
          temp -= result[i] * rightHess[j];
          temp -= resultDer[i] * rightDer2[j];
          temp -= resultDer2[i] * rightDer[j];
          temp -= resultHess[i] * right[j];
        }
        resultHess[coeffNo] = temp/(*right);
      }
    }
  }

  /// hand optimized code for third order jet propagation of division
  template<class T>
  inline void evalC3(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC2(left,right,result,dim,order,coeffNo);
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      T temp = left[i111+coeffNo] - result[coeffNo]*right[i111] - result[i1+coeffNo]*right[i11] - result[i11+coeffNo]*right[i1];
      // case dxdxdx
      for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
      {
        temp -= result[i]* right[i111+j];
        temp -= result[i1+i] * right[i11+j];
        temp -= result[i11+i] * right[i1+j];
        temp -= result[i111+i]* right[j];
      }
      result[i111+coeffNo] = temp/(*right);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        temp = left[i112+coeffNo]
              - result[coeffNo]     * right[i112]
              - result[i1+coeffNo]  * right[i12]
              - result[i2+coeffNo]  * right[i11]
              - result[i12+coeffNo] * right[i1]
              - result[i11+coeffNo] * right[i2];
        T temp2 = left[i122+coeffNo]
              - result[coeffNo]     * right[i122]
              - result[i1+coeffNo]  * right[i22]
              - result[i2+coeffNo]  * right[i12]
              - result[i12+coeffNo] * right[i2]
              - result[i22+coeffNo] * right[i1];
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        {
          temp -= result[i]*right[i112+j];    // 0,xxy
          temp -= result[i1+i]*right[i12+j];  // x,xy
          temp -= result[i2+i]*right[i11+j];  // y,xx
          temp -= result[i11+i]*right[i2+j];  // xx,y
          temp -= result[i12+i]*right[i1+j];  // xy,x
          temp -= result[i112+i]*right[j];    // xxy,0

          temp2 -= result[i]*right[i122+j];   // 0,xyy
          temp2 -= result[i1+i]*right[i22+j]; // x,yy
          temp2 -= result[i2+i]*right[i12+j]; // y,xy
          temp2 -= result[i12+i]*right[i2+j]; // xy,y
          temp2 -= result[i22+i]*right[i1+j]; // yy,x
          temp2 -= result[i122+i]*right[j];   // xyy,0
        }
        result[i112+coeffNo] = temp/(*right);
        result[i122+coeffNo] = temp2/(*right);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          temp = left[i123+coeffNo]
                - result[coeffNo]     * right[i123]
                - result[i1+coeffNo]  * right[i23]
                - result[i2+coeffNo]  * right[i13]
                - result[i3+coeffNo]  * right[i12]
                - result[i12+coeffNo] * right[i3]
                - result[i13+coeffNo] * right[i2]
                - result[i23+coeffNo] * right[i1];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp -= result[i] * right[i123+j];    // 0,xyz
            temp -= result[i1+i] * right[i23+j];  // x,yz
            temp -= result[i2+i] * right[i13+j];  // y,xz
            temp -= result[i3+i] * right[i12+j];  // z,xy
            temp -= result[i12+i] * right[i3+j];  // xy,z
            temp -= result[i13+i] * right[i2+j];  // xz,y
            temp -= result[i23+i]*right[i1+j];    // yz,x
            temp -= result[i123+i]*right[j];      // xyz,0
          }
          result[i123+coeffNo] = temp/(*right);
        }
      }
    }
  }

  template<class T>
  void evalCn(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
    {
      T r = left[i+coeffNo] - result[0]*right[i+coeffNo];
      for(k=0;k<coeffNo;++k)
        r -= result[coeffNo-k]*right[i+k] + result[i+k]*right[coeffNo-k];
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
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          T r = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) r += result[s1+k]*right[s2+coeffNo-k] + result[s2+k]*right[s1+coeffNo-k];
          result[(shift + c.index(degree))*order+coeffNo] -= r;
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
        b[0]=degree-i;
        const unsigned s1 = (shiftA + a.index(i))*order;
        do{
          const unsigned s2 = (shiftA + b.index(degree-i))*order;
          if(s1<s2) continue;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          T r = TypeTraits<T>::zero();
          if(s1!=s2)
            for(k=0;k<=coeffNo;++k) r += result[s1+k]*right[s2+coeffNo-k] + result[s2+k]*right[s1+coeffNo-k];
          else
            for(k=0;k<=coeffNo;++k) r += result[s1+k]*right[s1+coeffNo-k];
          result[(shift + c.index(degree))*order+coeffNo] -= r;
        }while(b.hasNext());
      }while(a.hasNext());
    }
    for(i=shift*order+coeffNo;i<(unsigned)(end*order+coeffNo);i+=order)
      result[i] /= (*right);
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
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left)/(*right);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left + order;
    const T* rightDer = right + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
      *resultDer = ((*leftDer) - (*result) * (*rightDer))/(*right);
  }

  template<class T>
  inline void evalC2HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    const unsigned s = dim*order;
    // begin of C^1
    const T* rightDer = right+order;
    T* resultDer = result + order;
    // begin of C^2
    const T* leftHess = left + order + s;
    const T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,rightDer+=order,resultDer+=order)
    {
      // case dxdx
      *resultHess = ((*leftHess) - (*result)*(*rightHess) - (*resultDer)*(*rightDer))/(*right);

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      T* resultDer2 = resultDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,resultDer2+=order,rightDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
        *resultHess = ((*leftHess) - (*result)*(*rightHess) - (*resultDer2)*(*rightDer) - (*resultDer)*(*rightDer2))/(*right);
    }
  }

  /// hand optimized code for third order jet propagation of division
  template<class T>
  inline void evalC3HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      result[i111] = (left[i111] - (*result)*right[i111] - result[i1]*right[i11] - result[i11]*right[i1])/(*right);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        result[i112] = (left[i112] - (*result)* right[i112]
                          - result[i1]  * right[i12]
                          - result[i2]  * right[i11]
                          - result[i12] * right[i1]
                          - result[i11] * right[i2])/(*right);
        result[i122] = (left[i122] - (*result) * right[i122]
                             - result[i1]  * right[i22]
                             - result[i2]  * right[i12]
                             - result[i12] * right[i2]
                             - result[i22] * right[i1])/(*right);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          result[i123] = (left[i123] - (*result) * right[i123]
                            - result[i1]  * right[i23]
                            - result[i2]  * right[i13]
                            - result[i3]  * right[i12]
                            - result[i12] * right[i3]
                            - result[i13] * right[i2]
                            - result[i23] * right[i1])/(*right);
        }
      }
    }
  }

  template<class T>
  void evalCnHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
      result[i] = left[i] - result[0]*right[i];

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
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          result[(shift + c.index(degree))*order] -= result[s1] * right[s2] + result[s2] * right[s1];
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
          if(s1>s2) continue;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          if(s1!=s2)
            result[(shift + c.index(degree))*order] -= result[s1] * right[s2] + result[s2] * right[s1];
          else
            result[(shift + c.index(degree))*order] -= result[s1] * right[s1];
        }while(b.hasNext());
      }while(a.hasNext());
    }
    for(i=shift*order;i<(unsigned)(end*order);i+=order)
      result[i] /= (*right);
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

// -------------------- DivVarByConst  -------------------------------

namespace DivVarByConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo]/(*right);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left+coeffNo;
    T* resultDer = result+coeffNo;
    for(unsigned derNo=0;derNo<=dim;++derNo,leftDer+=order,resultDer+=order)
      (*resultDer) = (*leftDer)/(*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left)/(*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      left += s;
      result += s;
      for(unsigned derNo=0;derNo<binomial(dim-1+degree,degree);++derNo,left+=order,result+=order)
        *result = (*left) / (*right);
    }else
      *result = (*left)/(*right);
  }
}

// -------------------- DivVarByFunTime  -------------------------------

namespace DivVarByFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Div::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    for(unsigned derNo=0;derNo<binomial(dim+degree,dim);++derNo,left+=order,result+=order)
      Div::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left)/(*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    DivVarByConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- DivVarByTime  -------------------------------

namespace DivVarByTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = (left[coeffNo]-result[coeffNo-1])/(*right);
    else
      *result = (*left)/(*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    if(coeffNo){ 
      for(unsigned derNo=0;derNo<binomial(dim+degree,dim);++derNo,left+=order,result+=order)
        result[coeffNo] = (left[coeffNo]-result[coeffNo-1])/(*right);
    } else {
      for(unsigned derNo=0;derNo<binomial(dim+degree,dim);++derNo,left+=order,result+=order)
        *result = (*left)/(*right);
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left)/(*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    DivVarByConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- DivTimeByConst -------------------------------

namespace DivTimeByConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo<=1)
      result[coeffNo] = left[coeffNo] / (*right);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) / (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) / (*right);
  }
}

// -------------------- DivFunTimeByConst -------------------------------

namespace DivFunTimeByConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] / (*right);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) / (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) / (*right);
  }
}

// -------------------- DivFunTimeByTime -------------------------------

namespace DivFunTimeByTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    DivVarByTime::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) / (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) / (*right);
  }
}

// -------------------- DivFunTimeByFunTime -------------------------------

namespace DivFunTimeByFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Div::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) / (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) / (*right);
  }
}

// -------------------- DivConstByConst -------------------------------

namespace DivConstByConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = (*left) / (*right);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) / (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) / (*right);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Div);
CAPD_MAKE_CLASS_NODE(DivVarByConst);
CAPD_MAKE_CLASS_NODE(DivVarByTime);
CAPD_MAKE_CLASS_NODE(DivVarByFunTime);
CAPD_MAKE_CLASS_NODE(DivTimeByConst);
CAPD_MAKE_CLASS_NODE(DivFunTimeByConst);
CAPD_MAKE_CLASS_NODE(DivFunTimeByTime);
CAPD_MAKE_CLASS_NODE(DivFunTimeByFunTime);
CAPD_MAKE_CLASS_NODE(DivConstByConst);

}} // namespace capd::autodiff

#endif

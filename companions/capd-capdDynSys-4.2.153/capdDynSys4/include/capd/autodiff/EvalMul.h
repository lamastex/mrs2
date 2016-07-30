/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalMul.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_MUL_H_
#define _CAPD_AUTODIFF_EVAL_MUL_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Mul ------------------------------------

namespace Mul
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    T temp = TypeTraits<T>::zero();
    for(unsigned i=0;i<=coeffNo;++i)
      temp += left[i]* right[coeffNo-i];
    result[coeffNo] = temp;
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    const T* leftDer = left+order;
    const T* rightDer = right+order;
    result += (order+coeffNo);
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,result+=order)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        temp += (left[i]*rightDer[j] + leftDer[i]*right[j]);
      *result = temp;
    }
  } // evalC1

  /// hand optimized code for second order jet propagation of multiplication
  template<class T>
  inline void evalC2(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);

    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left+order;
    const T* rightDer = right+order;
    T* resultDer = result + order + coeffNo;
    // begin of C^2
    const T* leftHess = leftDer + s;
    const T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      // case dx^2 and dx
      T temp = TypeTraits<T>::zero();
      T temp2 = TypeTraits<T>::zero();
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
      {
        temp += (left[i]*rightHess[j] + leftDer[i]*rightDer[j] + leftHess[i]*right[j]);
        temp2 += (left[i]*rightDer[j] + leftDer[i]*right[j]);
      }
      *resultDer = temp2;
      *resultHess = temp;

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,rightDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        temp = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          temp += left[i] * rightHess[j];
          temp += leftDer[i] * rightDer2[j];
          temp += leftDer2[i] * rightDer[j];
          temp += leftHess[i] * right[j];
        }
        *resultHess = temp;
      }
    }
  }  // evalC2

  /// hand optimized code for third order jet propagation of multiplication
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
        temp += left[i]* right[i111+j];
        temp += left[i1+i] * right[i11+j];
        temp += left[i11+i] * right[i1+j];
        temp += left[i111+i]* right[j];
      }
      result[i111+coeffNo] = temp;

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
          temp += left[i]*right[i112+j];    // 0,xxy
          temp += left[i1+i]*right[i12+j];  // x,xy
          temp += left[i2+i]*right[i11+j];  // y,xx
          temp += left[i11+i]*right[i2+j];  // xx,y
          temp += left[i12+i]*right[i1+j];  // xy,x
          temp += left[i112+i]*right[j];    // xxy,0

          temp2 += left[i]*right[i122+j];   // 0,xyy
          temp2 += left[i1+i]*right[i22+j]; // x,yy
          temp2 += left[i2+i]*right[i12+j]; // y,xy
          temp2 += left[i12+i]*right[i2+j]; // xy,y
          temp2 += left[i22+i]*right[i1+j]; // yy,x
          temp2 += left[i122+i]*right[j];   // xyy,0
        }
        result[i112+coeffNo] = temp;
        result[i122+coeffNo] = temp2;

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          T temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            temp += left[i] * right[i123+j];    // 0,xyz
            temp += left[i1+i] * right[i23+j];  // x,yz
            temp += left[i2+i] * right[i13+j];  // y,xz
            temp += left[i3+i] * right[i12+j];  // z,xy
            temp += left[i12+i] * right[i3+j];  // xy,z
            temp += left[i13+i] * right[i2+j];  // xz,y
            temp += left[i23+i]*right[i1+j];    // yz,x
            temp += left[i123+i]*right[j];      // xyz,0
          }
          result[i123+coeffNo] = temp;
        }
      }
    }
  }  // evalC3

  template<class T>
  void evalCn(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
    {
      T r = TypeTraits<T>::zero();
      for(k=0;k<=coeffNo;++k)
        r += left[k]*right[i+coeffNo-k] + left[i+coeffNo-k]*right[k];
      result[i+coeffNo] = r;
    }

    Multiindex a(dim),b(dim),c(dim);
    for(i=1;i<(degree+1)/2;++i){
      a.clear();
      a[0]=i;
      const unsigned shiftA = binomial(dim+i-1,dim);
      const unsigned shiftB = binomial(dim+(degree-i)-1,dim);
      do{
        const unsigned s1 = (shiftA + a.index(i))*order;
        b.clear();
        b[0]=degree-i;
        do{
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          const unsigned s2 = (shiftB + b.index(degree-i))*order;
          T r = TypeTraits<T>::zero();
          for(k=0;k<=coeffNo;++k) r += left[s1+k]*right[s2+coeffNo-k] + left[s2+k]*right[s1+coeffNo-k];
          result[(shift + c.index(degree))*order + coeffNo] += r;
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
        const unsigned s1 = (shiftA + a.index(i))*order;
        b.clear();
        b[0]=i;
        do{
          const unsigned s2 = (shiftA + b.index(i))*order;
          if(s1>s2) continue;
          for(k=0;k<dim;++k) c[k]=a[k]+b[k];
          T r = TypeTraits<T>::zero();
          if(s1!=s2)
            for(k=0;k<=coeffNo;++k) r += left[s1+k]*right[s2+coeffNo-k] + left[s2+k]*right[s1+coeffNo-k];
          else
            for(k=0;k<=coeffNo;++k) r += left[s1+k]*right[s1+coeffNo-k];
          result[(shift + c.index(degree))*order + coeffNo] += r;
        }while(b.hasNext());
      }while(a.hasNext());
    }
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
    (*result) = (*left)*(*right);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left + order;
    const T* rightDer = right + order;
    T* resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
      *resultDer = (*left)*(*rightDer) + (*leftDer) *(*right);
  }

  template<class T>
  inline void evalC2HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left + order;
    const T* rightDer = right + order;
    T* resultDer = result + order;
    const unsigned s = dim*order;
    const T* leftHess = leftDer + s;
    const T* rightHess = rightDer + s;
    T* resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      // case dxdx
      *resultHess  = (*left)     * (*rightHess)
                   + (*leftDer)  * (*rightDer)
                   + (*leftHess) * (*right);

      // case dxdy
      resultHess +=order;
      leftHess += order;
      rightHess += order;
      const T* leftDer2 = leftDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,resultHess+=order,leftHess+=order,rightHess+=order,leftDer2+=order,rightDer2+=order)
      {
        *resultHess  = (*left)     * (*rightHess)
                     + (*leftDer2) * (*rightDer)
                     + (*leftDer)  * (*rightDer2)
                     + (*leftHess) * (*right);
      }
    }
  }

  template<class T>
  inline void evalC3HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      result[i111] = (*left)*right[i111] + left[i1]*right[i11] + left[i11]*right[i1] + left[i111]*(*right);

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        result[i112] = (*left)*right[i112] + left[i1]*right[i12] + left[i2]*right[i11] + left[i11]*right[i2] + left[i12]*right[i1] + left[i112]*(*right);
        result[i122] = (*left)*right[i122] + left[i1]*right[i22] + left[i2]*right[i12] + left[i12]*right[i2] + left[i22]*right[i1] + left[i122]*(*right);

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)

          result[i123]
              = (*left)    * right[i123]
              + left[i1]   * right[i23]
              + left[i2]   * right[i13]
              + left[i3]   * right[i12]
              + left[i12]  * right[i3]
              + left[i13]  * right[i2]
              + left[i23]  * right[i1]
              + left[i123] * (*right);

      }
    }
  }  // evalC3

  template<class T>
  void evalCnHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    using capd::vectalg::Multiindex;
    const unsigned shift = binomial(dim+degree-1,dim);
    const unsigned end = binomial(dim+degree,dim);
    unsigned i,k;
    for(i=shift*order;i<end*order;i+=order)
      result[i] = *(left)*right[i] + left[i]*(*right);

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
          result[(shift + c.index(degree))*order] += left[s1]*right[s2] + left[s2]*right[s1];
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
          result[(shift + c.index(degree))*order] +=
              (s1!=s2) ? left[s1]*right[s2] + left[s2]*right[s1] : left[s1]*right[s1];
        }while(b.hasNext());
      }while(a.hasNext());
    }
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

// -------------------- MulConstByVar -------------------------------

namespace MulConstByVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = (*left) * right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    right+= coeffNo;
    result+= coeffNo;
    for(unsigned derNo=0;derNo<=dim;++derNo,right+=order,result+=order)
      *result = (*left) * (*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      right += s;
      result += s;
      for(unsigned derNo=0;derNo<(unsigned)binomial(dim-1+degree,degree);++derNo,right+=order,result+=order)
        *result = (*left) * (*right);
    } else
      *result = (*left) * (*right);
  }
}

// -------------------- MulConstByConst -------------------------------

namespace MulConstByConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = (*left) * (*right);
  }


  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) * (*right);
  }
}

// -------------------- MulConstByFunTime -------------------------------

namespace MulConstByFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = (*left) * right[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) * (*right);
  }
}

// -------------------- MulConstByTime -------------------------------

namespace MulConstByTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 0:
        *result = (*left) * (*right); break;
      case 1:
        result[1] = *left;
    }
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) * (*right);
  }
}

// -------------------- MulTimeByVar -------------------------------

namespace MulTimeByVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 0:
      *result = (*left) * (*right);
      break;
    default:
      result[coeffNo] = (*left) * right[coeffNo] + right[coeffNo-1];
    }
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 0:
      for(unsigned i=0;i<=dim;++i,result+=order,right+=order)
        *result = (*left) * (*right);
      break;
    default:
      result+= coeffNo;
      right += coeffNo;
      for(unsigned i=0;i<=dim;++i,result+=order,right+=order)
        *result = (*left) * (*right) + *(right-1);
    }
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,degree)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    MulConstByVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- MulTimeByFunTime -------------------------------

namespace MulTimeByFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    MulTimeByVar::evalC0(left,right,result,coeffNo);
  }


  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    MulTimeByVar::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) * (*right);
  }
}


// -------------------- MulFunTimeByVar -------------------------------

namespace MulFunTimeByVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    result+= coeffNo;
    for(unsigned derNo=0;derNo<=dim;++derNo,result+=order,right+=order)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        temp += left[i]* right[j];
      *result = temp;
    }
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,degree)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    MulConstByVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- MulFunTimeByFunTime -------------------------------

namespace MulFunTimeByFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = (*left) * (*right);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = (*left) * (*right);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Mul);
CAPD_MAKE_CLASS_NODE(MulConstByVar);
CAPD_MAKE_CLASS_NODE(MulConstByConst);
CAPD_MAKE_CLASS_NODE(MulConstByTime);
CAPD_MAKE_CLASS_NODE(MulConstByFunTime);
CAPD_MAKE_CLASS_NODE(MulTimeByVar);
CAPD_MAKE_CLASS_NODE(MulTimeByFunTime);
CAPD_MAKE_CLASS_NODE(MulFunTimeByVar);
CAPD_MAKE_CLASS_NODE(MulFunTimeByFunTime);

}} // namespace capd::autodiff

#endif

/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalNaturalPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_NATURALPOW_H_
#define _CAPD_AUTODIFF_EVAL_NATURALPOW_H_

#include <algorithm>
#include "capd/autodiff/EvalPow.h"

namespace capd{
namespace autodiff{

/**
 * Auxiliary function.
 * Computes d^i x^c where c is integer, i>0 and x can be zero
 * @param x - array of coefficients of coefficient
 */

template<class T>
inline
T* evalC0SingularNaturalPow(const unsigned coeffNo, const T* x, T* temp1, T* temp2, const unsigned c)
{
  for(unsigned i=1;i<c;++i)
  {
    for(unsigned p=0;p<=coeffNo;++p)
    {
      temp2[p] = TypeTraits<T>::zero();
      for(unsigned j=0;j<=p;++j)
        temp2[p] += x[j]*temp1[p-j];
    }
    std::swap(temp1,temp2);
  }
  return temp1;
}

// -------------------- NaturalPow ------------------------------------

namespace NaturalPow
{

  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(coeffNo){
      if(!(isSingular(*left)))
      {
        T temp = capd::TypeTraits<T>::zero();
        for(int j=0;j<(int)coeffNo;++j)
          temp += (c*((int)coeffNo-j)-j) * result[j]* left[coeffNo-j];
        result[coeffNo] = temp/((double)coeffNo * (*left));
      }
      else
      {
        T* t = new T[2*(coeffNo+1)];
        std::copy(left,left+coeffNo+1,t);
        T* p = evalC0SingularNaturalPow(coeffNo,left,t,t+coeffNo+1,c-1);

        T temp = capd::TypeTraits<T>::zero();
        for(unsigned i=0;i<=coeffNo;++i)
          temp += left[i]*p[coeffNo-i];
        result[coeffNo] = temp;
        delete[]t;
      }
    }
    else
      *result = power(*left,c);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    const T* leftDer = left + order;
    T* resultDer = result + order;

    if(coeffNo)
    {
      if(!(isSingular(*left))) // d^i x^m where i>0 and x!=0
      {
        // C^0 part
        T temp = capd::TypeTraits<T>::zero();
        for(int j=0;j<(int)coeffNo;++j)
          temp += (c*((int)coeffNo-j)-j) * result[j] * left[coeffNo-j];
        result[coeffNo] = temp/((double)coeffNo * (*left));

        // C^1 part
        for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
        {
          T temp1 = result[coeffNo] * (*leftDer);
          T temp2 = capd::TypeTraits<T>::zero();
          for(unsigned j=0;j<coeffNo;++j)
          {
            temp1 += result[j] * leftDer[coeffNo-j];
            temp2 += left[coeffNo-j] * resultDer[j];
          }
          resultDer[coeffNo] = (temp1*(*right)-temp2)/(*left);
        }
      }
      else // d^i x^m where i>0 and x is singular. Recompute coefficients from nonrecursive formula
      {
        T* t = new T[2*(coeffNo+1)];
        std::copy(left,left+coeffNo+1,t);
        T* p = evalC0SingularNaturalPow(coeffNo,left,t,t+coeffNo+1,c-1);

        T temp = capd::TypeTraits<T>::zero();
        for(unsigned i=0;i<=coeffNo;++i)
          temp += left[i]*p[coeffNo-i];
        result[coeffNo] = temp;

        for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
        {
          temp = capd::TypeTraits<T>::zero();
          for(unsigned j=0;j<=coeffNo;++j)
            temp += p[j] * leftDer[coeffNo-j];
          resultDer[coeffNo] = ((double)c)*temp;
        }
        delete[]t;
      }
    }
    else
    {
      *result = power(*left,c);
      T temp = power(*left,c-1)*c;
      for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
        *resultDer = temp*(*leftDer);
    }
  } // evalC1

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1(left,right,result,dim,order,coeffNo);
        break;
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      default:
        throw std::logic_error("Jet propagation of NaturalPower is not implemented for degree>1");
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    const T* leftDer = left + order;
    T* resultDer = result + order;
    const int c = toInt(leftBound(*right));
    T t = c*power(*left,c-1);
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      *resultDer = t*(*leftDer);
    }
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
        throw std::logic_error("Jet propagation of NaturalPower is not implemented for degree>1");
    }
  }
}

// -------------------- NaturalPowFunTime ------------------------------------

namespace NaturalPowFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }


  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- NaturalPowTime ------------------------------------

namespace NaturalPowTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }


  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = power(*left, toInt(leftBound(*right)));
  }
}

// -------------------- NaturalPowConst ------------------------------------

namespace NaturalPowConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }


  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = power(*left, toInt(leftBound(*right)));
  }
}

// -------------------- IntegerPow ------------------------------------

namespace IntegerPow
{

  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j]* left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = power(*left,c);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 2:
        evalC0(left,right,result,coeffNo);
        Pow::evalC2(left,right,result,dim,order,coeffNo);
        break;
      case 1:
        evalC0(left,right,result,coeffNo);
        Pow::evalC1(left,right,result,dim,order,coeffNo);
        break;
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        Pow::evalC3(left,right,result,dim,order,coeffNo);
        for(unsigned i=4;i<=degree;++i)
          Pow::evalCn(i,left,right,result,dim,order,coeffNo);
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    switch(degree)
    {
      case 3:
        Pow::evalC3HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 2:
        Pow::evalC2HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 1:
        Pow::evalC1HomogenousPolynomial(left,right,result,dim,order);
        break;
      case 0:
        evalC0HomogenousPolynomial(left,right,result);
        break;
      default:
        Pow::evalCnHomogenousPolynomial(degree,left,right,result,dim,order);
    }
  }
}

// -------------------- IntegerPowFunTime ------------------------------------

namespace IntegerPowFunTime
{

  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    IntegerPow::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- IntegerPowTime ------------------------------------

namespace IntegerPowTime
{

  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(coeffNo)
      result[coeffNo] = (c-((int)coeffNo-1)) * result[coeffNo-1]/((double)coeffNo * (*left));
    else
      *result = power(*left,c);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- IntegerPowConst ------------------------------------

namespace IntegerPowConst
{

  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo!=0)
      *result = power(*left,toInt(leftBound(*right)));
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- HalfIntegerPow ------------------------------------


namespace HalfIntegerPow
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    double c = toDouble(leftBound(*right));
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = power(sqrt(*left),(unsigned)(2*c));
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
        temp2 += left[coeffNo-j] * resultDer[j];
      }
      resultDer[coeffNo] = (temp1*(*right)-temp2)/(*left);
    }
  } // evalC1

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1(left,right,result,dim,order,coeffNo);
        break;
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      default:
        throw std::logic_error("Jet propagation of HalfIntegerPow is not implemented for degree>1");
    }
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    const unsigned c = toInt(2.*leftBound(*right));
    *result = power(sqrt(*left),c);
  }

  template<class T>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    Pow::evalC1HomogenousPolynomial(left,right,result,dim,order);
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
        throw std::logic_error("Jet propagation of HalfIntegerPow is not implemented for degree>1");
    }
  }
}

// -------------------- HalfIntegerPowFunTime ------------------------------------


namespace HalfIntegerPowFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    const double c = toDouble(leftBound(*right));
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = power(sqrt(*left),(unsigned)(2*c));
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    const double c = toDouble(leftBound(*right));
    *result = power(sqrt(*left),(unsigned)(2*c));
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- HalfIntegerPowTime ------------------------------------


namespace HalfIntegerPowTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    const double c = toDouble(leftBound(*right));
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = power(sqrt(*left),(unsigned)(2*c));
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    unsigned c = toInt(2.*leftBound(*right));
    *result = power(sqrt(*left),c);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

// -------------------- HalfIntegerPowConst ------------------------------------


namespace HalfIntegerPowConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==0){
      const int c = toInt(2.*leftBound(*right));
      *result = power(sqrt(*left),c);
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
    const unsigned c = toInt(2.*leftBound(*right));
    *result = power(sqrt(*left),c);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      evalC0HomogenousPolynomial(left,right,result);
  }
}

CAPD_MAKE_CLASS_NODE(NaturalPow);
CAPD_MAKE_CLASS_NODE(NaturalPowConst);
CAPD_MAKE_CLASS_NODE(NaturalPowTime);
CAPD_MAKE_CLASS_NODE(NaturalPowFunTime);

CAPD_MAKE_CLASS_NODE(IntegerPow);
CAPD_MAKE_CLASS_NODE(IntegerPowConst);
CAPD_MAKE_CLASS_NODE(IntegerPowTime);
CAPD_MAKE_CLASS_NODE(IntegerPowFunTime);

CAPD_MAKE_CLASS_NODE(HalfIntegerPow);
CAPD_MAKE_CLASS_NODE(HalfIntegerPowConst);
CAPD_MAKE_CLASS_NODE(HalfIntegerPowTime);
CAPD_MAKE_CLASS_NODE(HalfIntegerPowFunTime);

}} // namespace capd::autodiff

#endif

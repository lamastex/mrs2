/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalSub.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SUB_H_
#define _CAPD_AUTODIFF_EVAL_SUB_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- Sub ------------------------------------

namespace Sub
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    left+=coeffNo;
    right+=coeffNo;
    result+=coeffNo;
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<=dim;++derNo,left+=order, right+=order,result+=order)
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      left += s;
      right += s;
      result +=s;
      for(unsigned derNo=0;derNo<binomial(dim-1+degree,degree);++derNo,left+=order,right+=order,result+=order)
        *result = *left - *right;
    } else
      *result = *left - *right;
  }
}

// -------------------- ConstMinusVar  --------------------------

namespace ConstMinusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = -right[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = - (*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      right += s;
      result +=s;
      for(unsigned derNo=0;derNo<binomial(dim-1+degree,degree);++derNo,right+=order,result+=order)
        *result = -(*right);
    } else
      *result = *left - *right;
  }
}

// -------------------- ConstMinusFunTime  --------------------------

namespace ConstMinusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = -right[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- ConstMinusTime  --------------------------

namespace ConstMinusTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==1)
      result[coeffNo] = -TypeTraits<T>::one();
    else if (coeffNo==0)
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- ConstMinusConst  --------------------------

namespace ConstMinusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- TimeMinusConst  --------------------------

namespace TimeMinusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==1)
      result[coeffNo] = TypeTraits<T>::one();
    else if (coeffNo==0)
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- TimeMinusFunTime  --------------------------

namespace TimeMinusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo>1)
      result[coeffNo] = - right[coeffNo];
    else
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- TimeMinusVar  --------------------------

namespace TimeMinusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    TimeMinusFunTime::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = - (*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    ConstMinusVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- FunTimeMinusVar  --------------------------

namespace FunTimeMinusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = - (*right);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    ConstMinusVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- FunTimeMinusFunTime  --------------------------

namespace FunTimeMinusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- FunTimeMinusTime  --------------------------

namespace FunTimeMinusTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo>1)
      result[coeffNo] = left[coeffNo];
    else
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- FunTimeMinusConst  --------------------------

namespace FunTimeMinusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = left[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left - *right;
  }
}

// -------------------- VarMinusConst --------------------------

namespace VarMinusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = left[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    left+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,left+=order,result+=order)
      *result = *left;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      left += s;
      result +=s;
      for(unsigned derNo=0;derNo<binomial(dim-1+degree,degree);++derNo,left+=order,result+=order)
        *result = *left;
    } else
      *result = *left - *right;
  }
}

// -------------------- VarMinusFunTime  --------------------------

namespace VarMinusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    left+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,left+=order,result+=order)
      *result =  (*left);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    VarMinusConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- VarMinusTime  --------------------------

namespace VarMinusTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    FunTimeMinusTime::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    left+=(coeffNo+order);
    result+=(coeffNo+order);
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<dim;++derNo,left+=order,result+=order)
      *result = (*left);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left - *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    VarMinusConst::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Sub);
CAPD_MAKE_CLASS_NODE(ConstMinusVar);
CAPD_MAKE_CLASS_NODE(ConstMinusFunTime);
CAPD_MAKE_CLASS_NODE(ConstMinusTime);
CAPD_MAKE_CLASS_NODE(ConstMinusConst);
CAPD_MAKE_CLASS_NODE(TimeMinusConst);
CAPD_MAKE_CLASS_NODE(TimeMinusFunTime);
CAPD_MAKE_CLASS_NODE(TimeMinusVar);
CAPD_MAKE_CLASS_NODE(FunTimeMinusConst);
CAPD_MAKE_CLASS_NODE(FunTimeMinusTime);
CAPD_MAKE_CLASS_NODE(FunTimeMinusFunTime);
CAPD_MAKE_CLASS_NODE(FunTimeMinusVar);
CAPD_MAKE_CLASS_NODE(VarMinusConst);
CAPD_MAKE_CLASS_NODE(VarMinusTime);
CAPD_MAKE_CLASS_NODE(VarMinusFunTime);

}} // namespace capd::autodiff

#endif

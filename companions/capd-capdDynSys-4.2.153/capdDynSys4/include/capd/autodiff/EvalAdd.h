/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalAdd.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ADD_H_
#define _CAPD_AUTODIFF_EVAL_ADD_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// ---------------------------- Add  ------------------------------------
namespace Add
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    // shift pointers to proper coefficients
    left += coeffNo;
    right += coeffNo;
    result += coeffNo;
    //computing C^0 and C^1 parts
    for(unsigned derNo=0;derNo<=dim;++derNo,left+=order, right+=order,result+=order)
      *result = *left + *right;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
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
        *result = *left + *right;
    } else
      *result = *left + *right;
  }
}

// -------------------- ConstPlusVar --------------------------
namespace ConstPlusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = right[coeffNo];
    else
      *result = *left + *right;
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+= (coeffNo+order);
    result+= (coeffNo+order);
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = *right;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
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
        *result = *right;
    } else
      *result = *left + *right;
  }
}

// -------------------- ConstPlusConst --------------------------
namespace ConstPlusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = *left + *right;
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left + *right;
  }
}

// -------------------- ConstPlusTime --------------------------
namespace ConstPlusTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 1:
      result[1] = TypeTraits<T>::one();
      break;
    case 0:
      *result = (*left) + (*right);
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
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left + *right;
  }
}

// -------------------- ConstPlusFunTime --------------------------
namespace ConstPlusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = right[coeffNo];
    else
      *result = (*left) + (*right);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left + *right;
  }
}

// -------------------- TimePlusVar --------------------------
namespace TimePlusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 0:
      *result = (*left) + (*right);
      break;
    case 1:
      result[1] = TypeTraits<T>::one() + right[1];
      break;
    default:
      result[coeffNo] = right[coeffNo];
    }
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+= (coeffNo+order);
    result+= (coeffNo+order);
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = *right;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    ConstPlusVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- TimePlusFunTime --------------------------
namespace TimePlusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    TimePlusVar::evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left + *right;
  }
}

// -------------------- FunTimePlusVar --------------------------
namespace FunTimePlusVar
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    right+= (coeffNo+order);
    result+= (coeffNo+order);
    for(unsigned derNo=0;derNo<dim;++derNo,right+=order,result+=order)
      *result = *right;
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(T* left, T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order)
  {
    ConstPlusVar::evalHomogenousPolynomial(degree,left,right,result,dim,order);
  }
}

// -------------------- FunTimePlusFunTime --------------------------

namespace FunTimePlusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* right, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, T* result)
  {
    *result = *left + *right;
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = *left + *right;
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(Add);
CAPD_MAKE_CLASS_NODE(ConstPlusVar);
CAPD_MAKE_CLASS_NODE(ConstPlusConst);
CAPD_MAKE_CLASS_NODE(ConstPlusTime);
CAPD_MAKE_CLASS_NODE(ConstPlusFunTime);
CAPD_MAKE_CLASS_NODE(TimePlusVar);
CAPD_MAKE_CLASS_NODE(TimePlusFunTime);
CAPD_MAKE_CLASS_NODE(FunTimePlusVar);
CAPD_MAKE_CLASS_NODE(FunTimePlusFunTime);

}} // namespace capd::autodiff

#endif

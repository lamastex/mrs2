/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file EvalUnaryMinus.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_UNARYMINUS_H_
#define _CAPD_AUTODIFF_EVAL_UNARYMINUS_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{

// -------------------- UnaryMinus -----------------------------

namespace UnaryMinus
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = -left[coeffNo];
  }

  template<class T>
  inline void evalC1(const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    left+=coeffNo;
    result+=coeffNo;
    // compute both C^0 and C^1
    for(unsigned derNo=0;derNo<=dim;++derNo,left+=order, result+=order)
      *result = -(*left);
  }

  template<class T>
  inline void eval(const unsigned degree, const T* left, const T* right, T* result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    evalC1(left,right,result,binomial(dim+degree,dim)-1,order,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = -(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned dim, const unsigned order)
  {
    if(degree)
    {
      const unsigned s = binomial(dim+degree-1,dim)*order;
      left += s;
      result +=s;
      for(unsigned derNo=0;derNo<binomial(dim-1+degree,degree);++derNo,left+=order,result+=order)
        *result = -(*left);
    } else
      *result = -(*left);
  }
}

// -------------------- UnaryMinusFunTime -----------------------------

namespace UnaryMinusFunTime
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    result[coeffNo] = -left[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    result[coeffNo] = -left[coeffNo];
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = -(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = -(*left);
  }
}

// -------------------- UnaryMinusTime -----------------------------

namespace UnaryMinusTime
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo<2)
      result[coeffNo] = -left[coeffNo];
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = -(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = -(*left);
  }
}

// -------------------- UnaryMinusConst -----------------------------

namespace UnaryMinusConst
{
  template<class T>
  inline void evalC0(const T* left, const T* /*right*/, T* result, const unsigned coeffNo)
  {
    if(coeffNo)
      *result = -(*left);
  }

  template<class T>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, T* result, const unsigned /*dim*/, const unsigned /*order*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, T* result)
  {
    *result = -(*left);
  }

  template<class T>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, T* result, const unsigned /*dim*/, const unsigned /*order*/)
  {
    if(degree==0)
      *result = -(*left);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_CLASS_NODE(UnaryMinus);
CAPD_MAKE_CLASS_NODE(UnaryMinusConst);
CAPD_MAKE_CLASS_NODE(UnaryMinusTime);
CAPD_MAKE_CLASS_NODE(UnaryMinusFunTime);

}} // namespace capd::autodiff

#endif

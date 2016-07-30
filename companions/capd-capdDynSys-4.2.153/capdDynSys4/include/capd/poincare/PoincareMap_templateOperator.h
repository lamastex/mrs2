/// @addtogroup poincare
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file  PoincareMap_templateOperator.h
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_
#define _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_

#include <cassert>
#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/BasicPoincareMap.hpp"

namespace capd{
namespace poincare{

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::computePoincareMap(T& originalSet, int n)
{
  // We move the set to be very close to the section
  T setAfterTheSection = this->reachSection(originalSet, n);
  VectorType bound = (VectorType)originalSet;
  ScalarType oneStepReturnTime;

  // try to cross section in one step and use Newton method
  if(this->crossSection<T>(originalSet,oneStepReturnTime,bound)){
    *(this->returnTime) = originalSet.getCurrentTime() + oneStepReturnTime;
    originalSet = setAfterTheSection;
    if(this->derivativeOfFlow!=0)
      *(this->derivativeOfFlow) = this->m_solver.getCurve().derivative(oneStepReturnTime);
    if(this->hessianOfFlow!=0)
      *(this->hessianOfFlow) = this->m_solver.getCurve().hessian(oneStepReturnTime);
    if(this->jet!=0){
      this->m_solver.getCurve().eval(oneStepReturnTime,*jet);
      (*jet)() = (typename JetType::VectorType)bound;
    }

    return bound;
  }
  // if not possible, try to cross section by the old method.
  return this->crossSection<T>(originalSet, setAfterTheSection);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, const VectorType& c, const MatrixType& A, ScalarType& out_returnTime, int n)
{
  // The originalSet after this function contains a set just after section.
  this->derivativeOfFlow = NULL;
  this->hessianOfFlow = NULL;
  this->jet = NULL;
  this->returnTime = &out_returnTime;

  // We move the set to be very close to the section
  const size_type D = A.numberOfRows();
  T setAfterTheSection = this->reachSection(originalSet, n);
  VectorType bound = (VectorType)setAfterTheSection;
  ScalarType oneStepReturnTime;
  // try to cross section in one step and use Newton method to resolve for return time
  if(this->crossSection<T>(originalSet,oneStepReturnTime,bound)){
    SaveStepControl<Solver> ssc(this->m_solver);
/*  Here we use the following estimates
    t - bound for return time
    A - const matrix
    c - const vector
    coordinates on section are given by A*(x(t)-c)

    t0 := mid(t)
    dt := t-t0

    x0 := mid(x(t0))
    dx := x(t0) - x0

    Using Taylor expansion of order 2 wrt. time and then mean value form we get:

    A*(x(t)-c) = A*(x(t0)-c)  +  A*f(x(t0))*dt  +  0.5*A*Df(x(t))*f(x(t))*dt^2
               = A*(x(t0)-c)  +  A*f(x0)*dt  +  A*Df(x(t0))*dt*deltaX  +  0.5*A*Df(x(t))*f(x(t))*dt^2

   VERY important:
   - first multiply "thin objects", like [ A*f(x0) ] * dt
   - evaluation of A*(x(t0)-c) must take into account representation of x(t0).
     This is hidden in the class that represents the set - usually this is doubleton or tripleton.
*/
    MatrixType M(D,D);
    *(this->returnTime) = originalSet.getCurrentTime() + oneStepReturnTime;
    ScalarType t0, dt;
    oneStepReturnTime.split(t0,dt);

    // estimate [ A*Df(x(t)) ] * [ f(x(t)) * (0.5*dt^2) ]
    VectorType fx = this->getVectorField()(*(this->returnTime),bound,M);
    VectorType result = (A*M)*(fx*(typename capd::TypeTraits<ScalarType>::Real(0.5)*sqr(dt)));
    this->m_solver.setStep(t0);
    originalSet.move(this->m_solver);
    VectorType X = (VectorType)(originalSet);
    //VectorType x0(D), dx(D);
    VectorType x0 = X, dx=X;
    split(X,x0,dx);

    // add [ A*Df(x(t0))] *  [ dt*dx ]
    M = A*this->getVectorField().derivative(originalSet.getCurrentTime(),X);
    VectorType second = M*(dx*dt);

    // add [ A*f(x0) ] * dt. VERY important: first multiply "thin objects", i.e. A and f(x0).
    VectorType first = (A*this->getVectorField()(originalSet.getCurrentTime(),x0))*dt;

    // intersect bounds and add to result
    VectorType Y = intersection(originalSet.affineTransformation(M,x0)*dt,first+second);
    result += Y;
    
    // eventually add A*(x(t0)-c) taking into account representation of the set. This is hidden in the class that represents the set.
    VectorType Y0 = originalSet.affineTransformation(A,c);
    result += originalSet.affineTransformation(A,c);

    // reset the set
    originalSet = setAfterTheSection;

    // intersect with naive estimation
    return intersection(A*(bound-c),result);
  }

  // if not possible, try to cross section by the old method.
  return A*(this->crossSection<T>(originalSet, setAfterTheSection)-c);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, ScalarType& out_returnTime, int n)
{
  // The originalSet after this function contains a set just after section.
  this->derivativeOfFlow = NULL;
  this->hessianOfFlow = NULL;
  this->jet = NULL;
  this->returnTime = &out_returnTime;

  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, ScalarType& out_returnTime, int n)
{
  // the originalSet after this function contains a set just after section
  // we move the set close to section
  this->derivativeOfFlow = &der;
  this->hessianOfFlow = NULL;
  this->jet = NULL;
  this->returnTime = &out_returnTime;

  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,der,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, HessianType& hessian, ScalarType& out_returnTime, int n)
{
  // the originalSet after this function contains a set just after section
  // we move the set close to section
  this->derivativeOfFlow = &der;
  this->hessianOfFlow = &hessian;
  this->jet = NULL;
  this->returnTime = &out_returnTime;

  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, HessianType& hessian, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,der,hessian,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, typename T::JetType& result, ScalarType& out_returnTime, int n)
{
  // the originalSet after this function contains a set just after section
  // we move the set close to section
  this->derivativeOfFlow = NULL;
  this->hessianOfFlow = NULL;
  this->jet = &result;
  this->returnTime = &out_returnTime;
  VectorType r = this->computePoincareMap(originalSet,n);
  result() = r;
  return r;
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, typename T::JetType& result, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,result,returnTime,n);
}

}} // namespace capd::poincare

#endif // _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_

/// @}

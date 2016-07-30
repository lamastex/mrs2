/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICCNSOLVER_HPP_
#define _CAPD_DYNSYS_BASICCNSOLVER_HPP_

#include "capd/basicalg/power.h"
#include "capd/dynsys/BasicCnSolver.h"
#include "capd/dynsys/BasicSolver.hpp"
#include "capd/diffAlgebra/BasicCnCurve.hpp"
#include "capd/diffAlgebra/CnCurve.hpp"

namespace capd{
namespace dynsys{

template <typename MapT, typename StepControlT,typename CurveT>
BasicCnSolver<MapT,StepControlT,CurveT>::~BasicCnSolver(){}

// ---------------------------- CONSTRUCTORS ---------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
BasicCnSolver<MapT,StepControlT,CurveT>::BasicCnSolver(
      MapType& vectorField,
      size_type order,
      const StepControlType& stepControl
    )   : capd::dynsys::StepControlInterface<StepControlType,ScalarType>(stepControl),
          SolutionCurve(0.0,0.0,vectorField.dimension(),order,vectorField.degree()),
          m_vField(&vectorField),
          m_fixedTimeStep(TypeTraits<ScalarType>::zero()),
          m_step(TypeTraits<ScalarType>::zero())
{
  this->m_vField->setOrder(order + 1);
  this->m_vField->differentiateTime();
  Multiindex::generateList(this->dimension(),this->degree(),this->m_listIndices);
  if(order<=0) throw std::logic_error("BasicCnTaylor constructor: order must be a positive integer");
}

// -----------------------------------------------------------------------------


template <typename MapType, typename StepControlType, typename CurveT>
typename BasicCnSolver<MapType, StepControlType, CurveT>::VectorType
BasicCnSolver<MapType, StepControlType, CurveT>::operator()(VectorType v)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = v;
  this->evalAndSum(v);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnSolver<MapT,StepControlT,CurveT>::operator()(VectorType v,MatrixType& der)
{
  // set initial condition
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = v[i];
    for(size_type j=0;j<this->dimension();++j)
      this->coefficient(i,j,0) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();
  }
  this->evalAndSum(v,der);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnSolver<MapT,StepControlT,CurveT>::operator()(
      VectorType v, const MatrixType& D, MatrixType& out_der
    )
{
  // set initial condition
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = v[i];
    for(size_type j=0;j<this->dimension();++j)
      this->coefficient(i,j,0) = D(i+1,j+1);
  }
  this->evalAndSum(v,out_der);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnSolver<MapT,StepControlT,CurveT>::operator()(VectorType v,MatrixType& der, HessianType& hessian)
{
  // set initial condition
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = v[i];
    for(size_type j=0;j<this->dimension();++j)
    {
      this->coefficient(i,j,0) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();
      for(size_type c=j;c<this->dimension();++c)
        this->coefficient(i,j,c,0) = TypeTraits<ScalarType>::zero();
    }
  }
  this->evalAndSum(v,der,hessian);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnSolver<MapT,StepControlT,CurveT>::operator()(
      VectorType v, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  // set initial condition
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = v[i];
    for(size_type j=0;j<this->dimension();++j)
    {
      this->coefficient(i,j,0) = D(i+1,j+1);
      for(size_type c=j;c<this->dimension();++c)
        this->coefficient(i,j,c,0) = H(i,j,c);
    }
  }
  this->evalAndSum(v,out_der,out_hessian);
  return v;
}


// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnSolver<MapT,StepControlT,CurveT>::operator()(ScalarType& t, JetType& coeff)
{
  this->setCurrentTime(t);
  this->setInitialCondition(coeff);
  this->evalAndSum(coeff);
  t += this->m_step;
  return VectorType(coeff);
}


// -----------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
void BasicCnSolver<MapType, StepControlType, CurveT>::setOrder(size_type order)
{
  if(order > this->getAllocatedOrder())
  {
    this->m_vField->setOrder(order + 1);
    this->m_vField->differentiateTime();
  }
  SolutionCurve::setOrder(order);
}

// ---------------------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
typename BasicCnSolver<MapType, StepControlType, CurveT>::ScalarType
BasicCnSolver<MapType, StepControlType, CurveT>::getCoeffNorm(size_type r, size_type degree) const
{
  ScalarType result = TypeTraits<ScalarType>::zero();
  size_type i;
  for(i=0;i<this->dimension();++i)
  {
    result = capd::max(result,capd::abs(this->centerCoefficient(i,r)));
    result = capd::max(result,capd::abs(this->coefficient(i,r)));
  }
  if(degree)
  {
    const typename CurveT::JetType& c = this->getCoefficients()[r];

    for(i=0;i<this->dimension();++i)
    {
      typename CurveT::JetType::const_iterator b = c.begin(i,1), e=c.end(i,degree);
      while(b!=e)
        result = capd::max(result,capd::abs(*b++));
    }
  }

  return result;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v,MatrixType& der, HessianType& hessian)
{
  size_type i,j,c;

  int r=this->getOrder();
  this->m_vField->computeODECoefficients(this->getCoefficients(),2,r);
  this->computeTimeStep(v);

  for(i=0;i<this->dimension();++i)
  {
    v[i] = this->centerCoefficient(i,r) = this->coefficient(i,r);
    for(j=0;j<this->dimension();++j)
    {
      der(i+1,j+1) = this->coefficient(i,j,r);
      for(c=j;c<this->dimension();++c)
        hessian(i,j,c) = this->coefficient(i,j,c,r);
    }
  }
  for(r=this->getOrder()-1;r>=0;--r)
  {
    for(i=0;i<this->dimension();++i)
    {
      this->centerCoefficient(i,r) = this->coefficient(i,r);
      v[i] = v[i]*this->m_step + this->coefficient(i,r);
      for(j=0;j<this->dimension();++j)
      {
        der(i+1,j+1) = der(i+1,j+1)*this->m_step + this->coefficient(i,j,r);
        for(c=j;c<this->dimension();++c)
          hessian(i,j,c) = hessian(i,j,c)*this->m_step + this->coefficient(i,j,c,r);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v,MatrixType& der)
{
  size_type i,j;

  int r=this->getOrder();
  this->m_vField->computeODECoefficients(this->getCoefficients(),1,r);
  this->computeTimeStep(v);

  for(i=0;i<this->dimension();++i)
  {
    v[i] = this->centerCoefficient(i,r) = this->coefficient(i,r);
    for(j=0;j<this->dimension();++j)
      der(i+1,j+1) = this->coefficient(i,j,r);
  }
  for(r=this->getOrder()-1;r>=0;--r)
  {
    for(i=0;i<this->dimension();++i)
    {
      this->centerCoefficient(i,r) = this->coefficient(i,r);
      v[i] = v[i]*this->m_step + this->coefficient(i,r);
      for(j=0;j<this->dimension();++j)
        der(i+1,j+1) = der(i+1,j+1)*this->m_step + this->coefficient(i,j,r);
    }
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v)
{
  int r=this->getOrder();
  this->m_vField->computeODECoefficients(this->getCoefficientsAtCenter(),r);
  this->computeTimeStep(v);
  v = this->getCoefficientsAtCenter()[r];
  for(--r; r >= 0; --r)
    capd::vectalg::multiplyAssignObjectScalarAddObject(v,this->m_step,this->getCoefficientsAtCenter()[r]);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnSolver<MapT,StepControlT,CurveT>::evalAndSum(JetType& v)
{
  this->m_vField->computeODECoefficients(this->getCoefficients(),v.degree(),this->getOrder());
  this->computeTimeStep(v());

  for(size_type i=0;i<this->dimension();++i)
  {
    ScalarType* p = &this->coefficient(i,this->getOrder());
    this->centerCoefficient(i,this->getOrder()) = *p;
    typename JetType::iterator b = v.begin(i), e = v.end(i);
    for(;b!=e;++b,++p)
      *b = *p;
    for(int r = this->getOrder()-1;r>=0;--r)
    {
      ScalarType* p = &this->coefficient(i,r);
      this->centerCoefficient(i,r) = *p;
      typename JetType::iterator b = v.begin(i), e = v.end(i);
      for(;b!=e;++b,++p)
        *b = (*b)*this->m_step + (*p);
    }
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnSolver<MapT,StepControlT,CurveT>::setInitialCondition(const JetType& v)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    typename JetType::const_iterator b = v.begin(i), e=v.end(i);
    ScalarType* p = &this->coefficient(i,0);
    while(b!=e)
    {
      *p = *b;
      b++;
      p++;
    }
  }
}

//###########################################################//

template <typename MapType, typename StepControlType, typename CurveT>
const typename BasicCnSolver<MapType, StepControlType, CurveT>::SolutionCurve&
BasicCnSolver<MapType, StepControlType, CurveT>::getCurve()
{
  this->setDomain(0.,rightBound(this->m_step));
  return *this;
}

}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICCNSOLVER_HPP_

/// @}

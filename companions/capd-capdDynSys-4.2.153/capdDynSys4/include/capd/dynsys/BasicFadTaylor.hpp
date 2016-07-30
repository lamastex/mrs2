/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicFadTaylor.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICFADTAYLOR_HPP_
#define _CAPD_DYNSYS_BASICFADTAYLOR_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/basicalg/factrial.h"

#include "capd/dynsys/BasicFadTaylor.h"
#include "capd/diffAlgebra/FadCurve.hpp"
#include "capd/diffAlgebra/Curve.hpp"

namespace capd {
namespace dynsys {

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::recordDags()
{
  m_centerOut = m_vectorField(m_time,this->m_center);
  m_out = m_vectorField(m_ftime,this->m_in);
  m_remOut = m_vectorField(m_time,this->m_rem);
  m_jacRemOut = m_vectorField(m_ftime,this->m_jacRem);
  m_time[1] = TypeTraits<ScalarType>::one();
  m_ftime[1].x() = TypeTraits<ScalarType>::one();
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::reset()
{
  m_time.reset();
  m_ftime.reset();
  for(size_type i=0;i<this->dimension();++i)
  {
    m_centerOut[i].reset();
    m_out[i].reset();
    m_remOut[i].reset();
    m_jacRemOut[i].reset();
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
BasicFadTaylor<FadMapT,StepControlT>::BasicFadTaylor(MapType& f, size_type _order, const StepControlT& stepControl)
  : capd::dynsys::StepControlInterface<StepControlT,ScalarType>(stepControl),
    SolutionCurve(0.0, 0.0,f.dimension(),_order,1/*degree*/),
    m_vectorField(f),
    m_centerOut(f.dimension(),true),
    m_out(f.dimension(),true),
    m_remOut(f.dimension(),true),
    m_jacRemOut(f.dimension(),true),
    m_fixedTimeStep(TypeTraits<ScalarType>::zero()),
    m_step(TypeTraits<ScalarType>::zero())
{
  this->recordDags();
  if(_order<=0) throw std::logic_error("BasicFadTaylor constructor: order must be a positive integer");
}


//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::setInitialCondition(const VectorType& u, TVector& in)
{
  for(size_type i=0;i<this->dimension();++i)
    in[i][0] = u[i];
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::setInitialCondition(const VectorType& u, TFVector& in)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    in[i][0] = u[i];
    differentiate(in[i][0],i,this->dimension());
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::setInitialCondition(const VectorType& u, const MatrixType& M, TFVector& in)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    in[i][0] = u[i];
    differentiate(in[i][0],i,this->dimension());
    for(size_type j=0;j<this->dimension();++j)
      in[i][0].d(j) = M(i+1,j+1);
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::sumTaylorSeries(VectorType& u, TVector& in, size_type order)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    u[i] = in[i][order];
    for(int r=order-1;r>=0;--r)
      u[i] = u[i]*m_step + in[i][r];
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::sumTaylorSeries(MatrixType& M,  TFVector& in, size_type order)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    for(size_type j=0;j<this->dimension();++j)
    {
      M(i+1,j+1) = in[i][order].d(j);
      for(int r=order-1;r>=0;--r)
        M(i+1,j+1) = M(i+1,j+1)*m_step + in[i][r].d(j);
    }
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void BasicFadTaylor<FadMapT,StepControlT>::sumTaylorSeries(VectorType& u, MatrixType& M,  TFVector& in, size_type order)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    u[i] = in[i][order].x();
    for(int r=order-1;r>=0;--r)
      u[i] = u[i]*m_step + in[i][r].x();

    for(size_type j=0;j<this->dimension();++j)
    {
      M(i+1,j+1) = in[i][order].d(j);
      for(int r=order-1;r>=0;--r)
        M(i+1,j+1) = M(i+1,j+1)*m_step + in[i][r].d(j);
    }
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename BasicFadTaylor<FadMapT,StepControlT>::VectorType
BasicFadTaylor<FadMapT,StepControlT>::operator()(VectorType u)
{
  this->setInitialCondition(u,this->m_center);
  this->computeCoeff(this->m_center,this->m_centerOut,this->getOrder());
  this->computeTimeStep(u);
  this->sumTaylorSeries(u,this->m_center,this->getOrder());
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename BasicFadTaylor<FadMapT,StepControlT>::VectorType
BasicFadTaylor<FadMapT,StepControlT>::operator()(VectorType u, const MatrixType& derivative, MatrixType& resultDerivative)
{
  this->setInitialCondition(u,derivative,this->m_in);
  this->computeCoeff(this->m_in,this->m_out,this->getOrder());
  this->computeTimeStep(u);
  this->sumTaylorSeries(u,resultDerivative,this->m_in,this->getOrder());
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename BasicFadTaylor<FadMapT,StepControlT>::VectorType
BasicFadTaylor<FadMapT,StepControlT>::operator()(VectorType u, MatrixType& resultDerivative)
{
  this->setInitialCondition(u,this->m_in);
  this->computeCoeff(this->m_in,this->m_out,this->getOrder());
  this->computeTimeStep(u);
  this->sumTaylorSeries(u,resultDerivative,this->m_in,this->getOrder());
  return u;
}

//###########################################################//

template <typename FadMapT, typename StepControlType>
typename BasicFadTaylor<FadMapT, StepControlType>::ScalarType
BasicFadTaylor<FadMapT, StepControlType>::getCoeffNorm(size_type r,size_type degree) const
{
  ScalarType result = TypeTraits<ScalarType>::zero();
  for(size_type i=0;i<this->dimension();++i)
  {
    result = capd::max(result,capd::abs(this->centerCoefficient(i,r)));
    result = capd::max(result,capd::abs(this->coefficient(i,r)));
  }
  if(degree){
    for(size_type i=0;i<this->dimension();++i)
      for(size_type j=0;j<this->dimension();++j)
        result = capd::max(result,capd::abs(this->coefficient(i,j,r)));
  }
  return result;
}

//###########################################################//

template <typename FadMapT, typename StepControlType>
const typename BasicFadTaylor<FadMapT, StepControlType>::SolutionCurve& BasicFadTaylor<FadMapT, StepControlType>::getCurve()
{
  this->setDomain(0.,rightBound(m_step));
  return *this;
}

}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICFADTAYLOR_HPP_
/// @}

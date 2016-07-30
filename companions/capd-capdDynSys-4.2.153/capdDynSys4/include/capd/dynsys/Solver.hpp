/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Solver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_SOLVER_HPP_
#define _CAPD_DYNSYS_SOLVER_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/basicalg/factrial.h"

#include "capd/dynsys/Solver.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/DynSys.hpp"
#include "capd/dynsys/BasicSolver.hpp"
#include "capd/dynsys/enclosure.hpp"
#include "capd/dynsys/approveRemainder.h"

namespace capd{
namespace dynsys{

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
Solver<MapType, StepControlType,CurveType>::Solver(MapType& vectorField, size_type a_order, const StepControlType& stepControl)
  : BaseTaylor(vectorField,a_order,stepControl), implicitCurve(1.0,1.0,vectorField.dimension(),a_order,vectorField.degree())
{}

//###########################################################//

template <typename MapType, typename StepControlType, typename CurveT>
void Solver<MapType, StepControlType, CurveT>::setOrder(size_type order)
{
  BaseTaylor::setOrder(order);
  this->implicitCurve.setOrder(order);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType,StepControlType,CurveType>::VectorType
Solver<MapType,StepControlType,CurveType>::Phi(const ScalarType& t, const VectorType& v)
{
  VectorType result(v.dimension());
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(v,this->getOrder());
  BaseTaylor::sumTaylorSeries(result,this->getCoefficientsAtCenter(),this->getOrder());
  return result;
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::MatrixType
Solver<MapType, StepControlType,CurveType>::JacPhi(const ScalarType& t, const VectorType &iv)
{
  BaseTaylor::computeCoefficients(t,iv,this->getOrder());

  // the summation of the Taylor series
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  MatrixType result = matrixCoeff[this->getOrder()];
  for(int r=this->getOrder()-1;r>=0;--r)
    capd::vectalg::multiplyAssignObjectScalarAddObject(result,this->m_step,matrixCoeff[r]);

  return result;
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::VectorType
Solver<MapType, StepControlType,CurveType>::Remainder(const ScalarType& t, const VectorType &iv, VectorType& o_enc)
{
  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());
  // we compute an enclosure for \varphi([0,timestep],iv)
  o_enc = this->enclosure(t,iv);
  this->computeRemainderCoefficients(t + I*this->m_step,o_enc);
  return this->getRemainderCoefficients()[this->getOrder()+1]*power(this->m_step,this->getOrder()+1);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::VectorType
Solver<MapType, StepControlType,CurveType>::enclosure(const ScalarType& t, const VectorType &x)
///< the function finds an enclosure for \varphi([0,step],x)
{
  return capd::dynsys::enclosure(*this->m_vField,t,x, this->m_step);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::MatrixType
Solver<MapType, StepControlType,CurveType>::jacEnclosure(const ScalarType& t, const VectorType &enc, ScalarType* logNormOfDerivative)
/**< the function finds enclosure for Jacobian matrix (variational part)
     source- "C^1-Lohner algorithm" by P. Zgliczynski */
{
  return capd::dynsys::jacEnclosure(*this->m_vField,t,this->m_step,enc,capd::vectalg::EuclLNorm<VectorType,MatrixType>(),logNormOfDerivative);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::MatrixType
Solver<MapType, StepControlType,CurveType>::jacEnclosure(const ScalarType& t, const VectorType &enc)
/**< the function finds enclosure for Jacobian matrix (variational part)
     source- "C^1-Lohner algorithm" by P. Zgliczynski */
{
  return capd::dynsys::jacEnclosure(*this->m_vField,t,this->m_step,enc,capd::vectalg::EuclLNorm<VectorType,MatrixType>());
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void Solver<MapType, StepControlType,CurveType>::JacRemainder(
      const ScalarType& t,
      const VectorType &vecEnc,
      const MatrixType &jacEnc,
      VectorType &Rem,
      MatrixType &jacRem
  )
/// the remainder for varaiational part is computed
/// vecEnc - enclosure for \varphi([0,step],x)
/// jacEnc -  enclosure for \partial/\partial x \varphi([0,step],x)
{
  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());
  this->computeRemainderCoefficients(t+I*this->m_step,vecEnc,jacEnc);
  ScalarType fac = power(this->m_step,this->getOrder()+1);
  Rem = this->getRemainderCoefficients()[this->getOrder()+1] * fac;
  jacRem = this->getMatrixRemainderCoefficients()[this->getOrder()+1] * fac;
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void Solver<MapType,StepControlType,CurveType>::encloseC0Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi
  )
{
  this->computeTaylorCoefficients(t,x,xx);
  capd::dynsys::computeAndApproveRemainder(*this,t,xx,o_rem,o_enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,this->getCoefficientsAtCenter(),this->getMatrixCoefficients(),this->getOrder());
}


//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void Solver<MapType,StepControlType,CurveType>::encloseC1Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc
  )
{
  capd::diffAlgebra::C1TimeJet<MatrixType> rem(&o_rem,&o_jacRem);
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);

  this->computeTaylorCoefficients(t,x,xx);
  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,this->getCoefficientsAtCenter(),this->getMatrixCoefficients(),this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename Solver<MapType, StepControlType,CurveType>::ScalarType
Solver<MapType, StepControlType,CurveType>::getCoeffNorm(size_type r, size_type degree) const
{
  typename TypeTraits<ScalarType>::Real result = 0;
  for(size_type i=0;i<this->dimension();++i){
    result = capd::max(result,rightBound(abs(this->remainderCoefficient(i,r))));
    result = capd::max(result,rightBound(abs(this->coefficient(i,r))));
  }
  if(degree)
  {
    for(size_type i=0;i<this->dimension();++i){
      for(size_type j=0;j<this->dimension();++j){
        result = capd::max(result,rightBound(abs(this->remainderCoefficient(i,j,r))));
        result = capd::max(result,rightBound(abs(this->coefficient(i,j,r))));
      }
    }
  }
  return ScalarType(result);
}

}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_SOLVER_HPP_

/// @}

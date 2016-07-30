/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file CnSolver_phi.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNSOLVER_PHI_HPP_
#define _CAPD_DYNSYS_CNSOLVER_PHI_HPP_

#include "capd/dynsys/CnSolver.h"

namespace capd{
namespace dynsys{

template<typename MapT, typename StepControlT, typename CurveT>
typename CnSolver<MapT,StepControlT,CurveT>::VectorType
CnSolver<MapT,StepControlT,CurveT>::Phi(const ScalarType&t, const VectorType &iv)
{
  this->setCurrentTime(t);
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = iv;
  this->m_vField->computeODECoefficients(coeff,this->getOrder());

  // summation of the Taylor series
  VectorType v = coeff[this->getOrder()];
  for(int r = this->getOrder() - 1; r >= 0; --r)
    capd::vectalg::multiplyAssignObjectScalarAddObject(v,this->m_step,coeff[r]);
  return v;
}

//###########################################################//

template<typename MapT, typename StepControlT, typename CurveT>
typename CnSolver<MapT,StepControlT,CurveT>::MatrixType
CnSolver<MapT,StepControlT,CurveT>::JacPhi(const ScalarType& t, const VectorType &iv)
{
  MatrixType result(dimension(),dimension());
  this->setCurrentTime(t);
  BaseTaylor::operator()(iv,result);
  return result;
}

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_CNSOLVER_PHI_HPP_

/// @}

/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicC2Solver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICC2SOLVER_HPP_
#define _CAPD_DYNSYS_BASICC2SOLVER_HPP_

#include <string>
#include <stdexcept>

#include "capd/dynsys/BasicC2Solver.h"
#include "capd/dynsys/BasicSolver.hpp"
#include "capd/diffAlgebra/BasicC2Curve.hpp"
#include "capd/diffAlgebra/C2Curve.hpp"

namespace capd{
namespace dynsys{

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
BasicC2Solver<MapType,StepControlType,CurveT>::BasicC2Solver(
      MapType& vectorField,
      size_type order,
      const StepControlType& stepControl
  ) : BaseTaylor(vectorField,order,stepControl)
{
  if(this->getVectorField().degree()<2)
    this->getVectorField().setDegree(2);
}


//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
typename BasicC2Solver<MapType,StepControlType,CurveT>::VectorType
BasicC2Solver<MapType,StepControlType,CurveT>::operator()(VectorType v,MatrixType& der, HessianType& hessian)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();
  coeff[0] = v;
  matrixCoeff[0].setToIdentity();
  hessianCoeff[0].clear();

  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,this->getOrder());
  this->computeTimeStep(v);

  this->sumTaylorSeries(v,der,hessian,coeff,matrixCoeff,hessianCoeff,this->getOrder());
  return v;
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
typename BasicC2Solver<MapType,StepControlType,CurveT>::VectorType
BasicC2Solver<MapType,StepControlType,CurveT>::operator()(
      VectorType x, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();
  coeff[0] = x;
  matrixCoeff[0] = D;
  hessianCoeff[0] = H;

  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,this->getOrder());
  this->computeTimeStep(x);

  this->sumTaylorSeries(x,out_der,out_hessian,coeff,matrixCoeff,hessianCoeff,this->getOrder());
  return x;
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
void BasicC2Solver<MapType,StepControlType,CurveT>::sumTaylorSeries(
      VectorType& v,
      MatrixType& der,
      HessianType& hessian,
      VectorType* coeff,
      MatrixType* matrixCoeff,
      HessianType* hessianCoeff,
      size_type order
  )
{
  // summation of the Taylor series
  v = coeff[order];
  der = matrixCoeff[order];
  hessian = hessianCoeff[order];
  for(int r = order - 1; r >= 0; --r)
  {
    capd::vectalg::multiplyAssignObjectScalarAddObject(v,this->m_step,coeff[r]);
    capd::vectalg::multiplyAssignObjectScalarAddObject(der,this->m_step,matrixCoeff[r]);
    capd::vectalg::multiplyAssignObjectScalarAddObject(hessian,this->m_step,hessianCoeff[r]);
  }
}

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICC2SOLVER_HPP_

/// @}

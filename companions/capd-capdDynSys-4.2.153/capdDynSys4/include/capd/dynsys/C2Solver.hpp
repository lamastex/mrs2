/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file C2Solver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C2SOLVER_HPP_
#define _CAPD_DYNSYS_C2SOLVER_HPP_

#include <string>
#include <stdexcept>

#include "capd/dynsys/C2Solver.h"
#include "capd/dynsys/Solver.hpp"
#include "capd/dynsys/BasicC2Solver.hpp"
#include "capd/vectalg/Norm.hpp"
#include "approveRemainder.h"

namespace capd{
namespace dynsys{

//###########################################################//

template<typename MapType,typename StepControlT, typename CurveT>
C2Solver<MapType,StepControlT,CurveT>::C2Solver(MapType& vectorField,size_type order)
  : BasicSolver<MapType,StepControlT,CurveT>(vectorField,order),
    Solver<MapType,StepControlT,CurveT>(vectorField,order),
    BasicC2Solver<MapType,StepControlT,CurveT>(vectorField,order)
{}

//###########################################################//

template<typename MapType,typename StepControlT, typename CurveT>
void C2Solver<MapType,StepControlT,CurveT>::c2Enclosure(const VectorType& enc, MatrixType& o_jacEnc, HessianType& o_hessianEnc)
{
  capd::dynsys::c2Enclosure(*(this->m_vField),this->m_step,enc,o_jacEnc,o_hessianEnc);
}

// ####################################################################

template<typename MapType,typename StepControlT, typename CurveT>
void C2Solver<MapType,StepControlT,CurveT>::c2Remainder(
      const VectorType& Enc,
      const MatrixType& jacEnc,
      const HessianType& hessianEnc,
      VectorType& o_Rem,
      MatrixType& o_jacRem,
      HessianType& o_hessianRem
  )
{
  VectorType* remCoeff = this->getRemainderCoefficients();
  MatrixType* matrixRemCoeff = this->getMatrixRemainderCoefficients();
  HessianType* hessianRemCoeff = this->getHessianRemainderCoefficients();
  remCoeff[0] = Enc;
  matrixRemCoeff[0] = jacEnc;
  hessianRemCoeff[0] = hessianEnc;
  this->m_vField->computeODECoefficients(remCoeff,matrixRemCoeff,hessianRemCoeff,this->getOrder()+1);
  ScalarType fac = power(this->m_step,this->getOrder()+1);
  o_Rem = remCoeff[this->getOrder()+1] * fac;
  o_jacRem = matrixRemCoeff[this->getOrder()+1] * fac;

  typename HessianType::iterator i1=o_hessianRem.begin(), i2=o_hessianRem.end(),
                                 j=hessianRemCoeff[this->getOrder()+1].begin();
  while(i1!=i2)
  {
    *i1 = *j * fac;
    ++i1;
    ++j;
  }
}

// ####################################################################

template<typename MapType,typename StepControlT, typename CurveT>
void C2Solver<MapType,StepControlT,CurveT>::encloseC2Map(
    const ScalarType& t,
    const VectorType& x,
    const VectorType& xx,
    VectorType& o_phi,
    VectorType& o_rem,
    VectorType& o_enc,
    MatrixType& o_jacPhi,
    MatrixType& o_jacRem,
    MatrixType& o_jacEnc,
    HessianType& o_hessianPhi,
    HessianType& o_hessianRem,
    HessianType& o_hessianEnc
  )
{
  VectorType* coeff = this->getCoefficients();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();

  const int order=this->getOrder();
  coeff[0] = xx;
  matrixCoeff[0].setToIdentity();
  hessianCoeff[0].clear();
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,order);
  this->computeCoefficientsAtCenter(x,order);
  
  capd::diffAlgebra::C2TimeJet<MatrixType> rem(&o_rem,&o_jacRem,&o_hessianRem);
  capd::diffAlgebra::C2TimeJet<MatrixType> enc(&o_enc,&o_jacEnc,&o_hessianEnc);

  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,o_hessianPhi,this->getCoefficientsAtCenter(),matrixCoeff,hessianCoeff,order);
}

// ####################################################################

template<typename MapType,typename StepControlT, typename CurveT>
void C2Solver<MapType,StepControlT,CurveT>::computeRemainder(ScalarType t, const VectorType& xx, C2TimeJetType& o_enc, C2TimeJetType& o_rem)
{
  o_enc.vector() = this->enclosure(t,xx);
  this->c2Enclosure(o_enc.vector(),o_enc.matrix(),o_enc.hessian());
  this->c2Remainder(o_enc.vector(),o_enc.matrix(),o_enc.hessian(),o_rem.vector(),o_rem.matrix(),o_rem.hessian());
}

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_C2SOLVER_HPP_

/// @}

/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file FadTaylorHOE.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_FADTAYLORHOE_HPP_
#define _CAPD_DYNSYS_FADTAYLORHOE_HPP_

#include "capd/dynsys/FadTaylorHOE.h"
#include "capd/dynsys/FadTaylor.hpp"
#include "capd/dynsys/highOrderEnclosure.h"

namespace capd{
namespace dynsys{

template<class FadMapT, typename StepControlT>
FadTaylorHOE<FadMapT,StepControlT>::FadTaylorHOE(MapType& f, size_type _order, const StepControlT& _stepControl)
  : BaseTaylor(f,_order,_stepControl)
{}

//###########################################################//

template<typename MapType, typename StepControlType>
typename FadTaylorHOE<MapType,StepControlType>::VectorType
FadTaylorHOE<MapType,StepControlType>::enclosure(const ScalarType& t, const VectorType& x, VectorType& outRemainder)
{
  bool stepChangeAllowance = this->isStepChangeAllowed();
  this->turnOffStepControl();
  this->setCurrentTime(t);
  this->computeCoefficients(x,this->getOrder());
  VectorType result(x.dimension());
  capd::dynsys::highOrderEnclosure(t,*this,outRemainder,result);
  this->onOffStepControl(stepChangeAllowance);
  return result;
}

//###########################################################//

template<typename MapType, typename StepControlType>
void FadTaylorHOE<MapType,StepControlType>::c1Enclosure(const ScalarType& t, const VectorType& x, VectorType& o_enc, MatrixType& o_jacEnc)
{
  bool stepChangeAllowance = this->isStepChangeAllowed();
  this->turnOffStepControl();
  this->setInitialCondition(x,this->m_in);
  this->setCurrentTime(t);
  this->computeCoeff(this->m_in,this->m_out,this->getOrder());

  capd::diffAlgebra::C1TimeJet<MatrixType> rem(x.dimension());
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);

  capd::dynsys::highOrderEnclosure(t,*this,rem,enc);

  this->onOffStepControl(stepChangeAllowance);
}

//###########################################################//

template<typename MapType, typename StepControlType>
typename FadTaylorHOE<MapType,StepControlType>::VectorType
FadTaylorHOE<MapType,StepControlType>::enclosure(const ScalarType& t, const VectorType& x)
{
  VectorType remainder(x.dimension());
  return this->enclosure(t,x,remainder);
}

//###########################################################//

template<typename MapType, typename StepControlType>
typename FadTaylorHOE<MapType,StepControlType>::VectorType
FadTaylorHOE<MapType,StepControlType>::Remainder(const ScalarType& t, const VectorType &x, VectorType& out_enc)
{
  VectorType remainder(x.dimension());
  out_enc = this->enclosure(t,x,remainder);
  return remainder;
}

//###########################################################//

template<typename MapType, typename StepControlType>
void FadTaylorHOE<MapType,StepControlType>::encloseC0Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi
  )
{
  // here we compute all the coefficients for phi(t) and DPhi(t)
  this->computeTaylorCoefficients(t,x,xx);

  this->computeTimeStep(t,xx);

  capd::dynsys::highOrderEnclosure(t,*this,o_rem,o_enc);

  this->sumTaylorSeries(o_phi,this->m_center,this->getOrder());
  this->sumTaylorSeries(o_jacPhi,this->m_in,this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlType>
void FadTaylorHOE<MapType,StepControlType>::encloseC1Map(
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
  // here we compute all the coefficients for phi(t) and DPhi(t)
  this->computeTaylorCoefficients(t,x,xx);

  this->computeTimeStep(t,xx);

  capd::diffAlgebra::C1TimeJet<MatrixType> rem(&o_rem,&o_jacRem);
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);
  capd::dynsys::highOrderEnclosure(t,*this,rem,enc);

  this->sumTaylorSeries(o_phi,this->m_center,this->getOrder());
  this->sumTaylorSeries(o_jacPhi,this->m_in,this->getOrder());
}

}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_FADTAYLORHOE_HPP_

/// @}

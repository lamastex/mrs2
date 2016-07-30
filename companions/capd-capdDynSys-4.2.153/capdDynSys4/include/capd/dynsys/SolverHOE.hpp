/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file SolverHOE.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_SOLVERHOE_HPP_
#define _CAPD_DYNSYS_SOLVERHOE_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/dynsys/highOrderEnclosure.h"
#include "capd/dynsys/Solver.hpp"
#include "capd/dynsys/SolverHOE.h"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace dynsys{

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
SolverHOE<MapType,StepControlType,CurveType>::SolverHOE(MapType& vectorField, size_type order, const StepControlType& stepControl)
  : BasicSolver<MapType,StepControlType,CurveType>(vectorField,order,stepControl),
    BaseTaylor(vectorField,order,stepControl)
{}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename SolverHOE<MapType,StepControlType,CurveType>::VectorType
SolverHOE<MapType,StepControlType,CurveType>::enclosure(const ScalarType& t, const VectorType& x, VectorType& out_Remainder)
///<the function finds an enclosure for \varphi([0,step],x)
{
  capd::poincare::SaveStepControl<SolverHOE> saveStepControl(*this);
  this->turnOffStepControl();
  this->setCurrentTime(t);
  this->computeCoefficients(x,this->getOrder());
  VectorType result(x.dimension());
  capd::dynsys::highOrderEnclosure(t,*this,out_Remainder,result);
  return result;
}


//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void SolverHOE<MapType,StepControlType,CurveType>::c1Enclosure(const ScalarType& t, const VectorType& x, VectorType& o_enc, MatrixType& o_jacEnc)
///<the function finds an enclosure for \varphi([0,step],x)
{
  capd::poincare::SaveStepControl<SolverHOE> saveStepControl(*this);
  this->turnOffStepControl();
  this->setCurrentTime(t);
  this->computeCoefficients(x,this->getOrder());
  capd::diffAlgebra::C1TimeJet<MatrixType> rem(this->dimension());
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);
  capd::dynsys::highOrderEnclosure(t,*this,rem,enc);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename SolverHOE<MapType,StepControlType,CurveType>::VectorType
SolverHOE<MapType,StepControlType,CurveType>::enclosure(const ScalarType& t, const VectorType& x)
///<the function finds an enclosure for \varphi([0,step],x)
{
  VectorType remainder(x.dimension());
  return this->enclosure(t,x,remainder);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
typename SolverHOE<MapType,StepControlType,CurveType>::VectorType
SolverHOE<MapType,StepControlType,CurveType>::Remainder(const ScalarType& t, const VectorType &x, VectorType& out_enc)
///<the function finds an enclosure for \varphi([0,step],x)
{
  VectorType remainder(x.dimension());
  out_enc = this->enclosure(t,x,remainder);
  return remainder;
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void SolverHOE<MapType,StepControlType,CurveType>::computeRemainder(
      ScalarType t,
      const VectorType& /*xx*/,
      VectorType& o_enc,
      VectorType& o_rem 
  ){
  capd::dynsys::highOrderEnclosure(t,*this,o_rem,o_enc);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
void SolverHOE<MapType,StepControlType,CurveType>::computeRemainder(
      ScalarType t, 
      const VectorType& /*xx*/, 
      capd::diffAlgebra::C1TimeJet<MatrixType>& o_enc, 
      capd::diffAlgebra::C1TimeJet<MatrixType>& o_rem  
  )
{
  capd::dynsys::highOrderEnclosure(t,*this,o_rem,o_enc);
}

}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_SOLVERHOE_HPP_

/// @}

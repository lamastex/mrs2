/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file SolverHOE.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_SOLVERHOE_H_
#define _CAPD_DYNSYS_SOLVERHOE_H_

#include <string>
#include <stdexcept>
#include "capd/dynset/C0Set.h"
#include "capd/dynset/C1Set.h"
#include "capd/dynsys/Solver.h"

namespace capd{
namespace dynsys{

template <
  typename MapT,
  typename StepControlT = capd::dynsys::ILastTermsStepControl,
  typename CurveT = capd::diffAlgebra::Curve< capd::diffAlgebra::BasicCurve<typename MapT::MatrixType> >
>
class SolverHOE: public Solver<MapT,StepControlT,CurveT>
{
public:
  typedef MapT MapType;
  typedef StepControlT StepControlType;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename ScalarType::BoundType BoundType;
  typedef Solver<MapT,StepControlT,CurveT> BaseTaylor;
  typedef typename MatrixType::size_type size_type;

  SolverHOE(MapType& vField, size_type order, const StepControlT& stepControl = StepControlT());

  /// Main method for simultaneous validation of the existence of solutions and computation of Lagrange remainder.
  VectorType enclosure(const ScalarType& t, const VectorType& x, VectorType& out_Remainder);

  // Overridden methods from DynSys
  VectorType Remainder(const ScalarType& t, const VectorType& x, VectorType& out_enc);
  VectorType enclosure(const ScalarType& t, const VectorType& x) ;
  void computeRemainder(ScalarType t, const VectorType& xx, VectorType& o_enc, VectorType& o_rem);
  void computeRemainder(ScalarType t, const VectorType& xx, capd::diffAlgebra::C1TimeJet<MatrixType>& o_enc, capd::diffAlgebra::C1TimeJet<MatrixType>& o_rem);

  // Overridden methods from C1DynSys
  /// computes simultaneously an enclosure for the solutions to ODE and associated variational equations
  void c1Enclosure(const ScalarType& t, const VectorType& x, VectorType& o_enc, MatrixType& o_jacEnc);

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1SetMove is used meaning that this solver can integrate C^0 and C^1 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
	  C1SetMove<SolverHOE,SetType>::move(set,*this);
  }
  
  template <typename SetType>
  void operator()(SetType& set, SetType & result){
    this->saveCurrentSet(set);
	  C1SetMove<SolverHOE,SetType>::move(set, result, *this);
  }


protected:
  void operator=(const SolverHOE& ){}
  SolverHOE(const SolverHOE& t) : BasicSolver<MapT,StepControlT>(t), BaseTaylor(t){}
};

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_SOLVERHOE_H_

/// @}

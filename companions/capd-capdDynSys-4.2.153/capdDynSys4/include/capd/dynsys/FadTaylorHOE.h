/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file FadTaylor.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_FADTAYLORHOE_H_
#define _CAPD_DYNSYS_FADTAYLORHOE_H_

#include "capd/dynsys/FadTaylor.h"

namespace capd{
namespace dynsys{

template<class FadMapT, typename StepControlT = capd::dynsys::ILastTermsStepControl >
class FadTaylorHOE : public FadTaylor<FadMapT,StepControlT>
{
public:
  typedef FadMapT MapType;
  typedef FadTaylor<FadMapT,StepControlT> BaseTaylor;
  typedef StepControlT StepControlType;
  typedef typename BaseTaylor::ScalarType ScalarType;
  typedef typename BaseTaylor::FScalar FScalar;
  typedef typename BaseTaylor::TFScalar TFScalar;
  typedef typename BaseTaylor::MatrixType MatrixType;
  typedef typename BaseTaylor::VectorType VectorType;
  typedef typename BaseTaylor::FVector FVector;
  typedef typename BaseTaylor::TFVector TFVector;
  typedef typename BaseTaylor::FunctionType FunctionType;
  typedef typename MatrixType::size_type size_type;

  FadTaylorHOE(MapType& f, size_type _order, const StepControlT& _stepControl=StepControlT());

  /// Main method for simultaneous validation of the existence of solutions and computation of Lagrange remainder.
  VectorType enclosure(const ScalarType& t, const VectorType& x, VectorType& out_Remainder);

  // Overridden methods from DynSys
  VectorType Remainder(const ScalarType& t, const VectorType &iv, VectorType &out_enc);
  VectorType enclosure(const ScalarType& t, const VectorType& x) ;
  void encloseC0Map(
      const ScalarType& t,  //< @param[in] current time of ODE
      const VectorType& x0, //< @param[in] an internal point of x, usually center of x
      const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
      VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
      VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi  //< @param[out] bound for derivative Dphi(x)
  );

  // Overridden methods from C1DynSys

  /// computes simultaneously an enclosure for the solutions to ODE and associated variational equations
  void c1Enclosure(const ScalarType& t, const VectorType& x, VectorType& o_enc, MatrixType& o_jacEnc);
  void encloseC1Map(
      const ScalarType& t,  //< @param[in] current time of ODE
      const VectorType& x0, //< @param[in] an internal point of x, usually center of x
      const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
      VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
      VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi, //< @param[out] bound for derivative Dphi(x)
      MatrixType& o_jacRem, //< @param[out] bound for the error of numerical method over the time step for variational equation
      MatrixType& o_jacEnc  //< @param[out] enclosure of all trajectories of variational equations with initial condition set to Identity over the time interval (time step of numerical method)
  );

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1SetMove is used meaning that this solver can integrate C^0 and C^1 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
	  C1SetMove<FadTaylorHOE,SetType>::move(set,*this);
  }
  
  /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0, C^1 set representing initial conditions
   /// @param[out] result    on return contains image of the set   
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
    this->saveCurrentSet(set);
	  C1SetMove<FadTaylorHOE,SetType>::move(set, result, *this);
  }

}; // the end of class FadTaylorHOE

}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_FADTAYLORHOE_H_

/// @}

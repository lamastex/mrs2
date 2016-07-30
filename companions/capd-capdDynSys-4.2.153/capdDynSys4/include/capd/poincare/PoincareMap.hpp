/// @addtogroup poincare
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file  PoincareMap.hpp
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_HPP_
#define _CAPD_POINCARE_POINCARE_MAP_HPP_

#include <cassert>
#include <exception>
#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/BasicPoincareMap.hpp"

namespace capd{
namespace poincare{



/**
 *   Constructs PoincareMap for given dynamical system and section
 *
 *  @param ds        dynamical system
 *  @param section   Poincare section
 *  @param crossing  section crossing direction
 *  @param factor    time step correction factor during section crossing (should be in interval (0, 1))
**/
template <typename SolverT, typename SectionT>
PoincareMap<SolverT,SectionT>::PoincareMap(
       Solver  & ds,
       SectionType & section,
       typename BasicPoincareMap<SolverT,SectionT>::CrossingDirection crossing,
       const BoundType & errorTolerance
  ) : BasicPoincareMap<SolverT,SectionT>(ds, section, crossing, errorTolerance),
  timeStepCorrectionFactor(0.9),
  maxCorrectionAttempts(10),
  minCrossingTimeStep( capd::TypeTraits< ScalarType >::epsilon()*4),
  sectionCoordinateSystem(MatrixType::Identity(ds.dimension()))
{}

/*__________________________________________________________________________*/

template<typename SolverT, typename FunT>
void PoincareMap<SolverT,FunT>::updateJacEnclosure(CnSet& S)
{
  if(jet!=NULL)
  {
    const size_type dim = this->getVectorField().dimension();
    JetType enc(dim, S.degree());

    substitutionPowerSeries(S.getLastJetEnclosure(), S.currentSet(),enc, false);

    for (size_type p = 0; p < this->getVectorField().dimension(); ++p)
    {
      typename JetType::iterator
          b = jet->begin(p, 1),
          e = jet->end(p, S.degree());
      typename JetType::iterator i = enc.begin(p, 1);
      while (b != e)
      {
        *b = intervalHull(*b, *i);
        ++b;
        ++i;
      }
    }
  }

  this->updateTime(S);
}

/*__________________________________________________________________________*/

template<typename SolverT, typename FunT>
void PoincareMap<SolverT,FunT>::saveJacEnclosure(CnSet& S)
{
  if(jet!=NULL)
    *jet = S.currentSet();
  this->saveTime(S);
}

/*__________________________________________________________________________*/

template <typename SolverT, typename SectionT>
void PoincareMap<SolverT,SectionT>::updateJacEnclosure(C1Set& S)
{
  if(derivativeOfFlow!=NULL)
  {
    MatrixType oneStepBound = S.getLastMatrixEnclosure()*MatrixType(S);
    intervalHull(*derivativeOfFlow,oneStepBound,*derivativeOfFlow);
  }
  this->updateTime(S);
}

/*__________________________________________________________________________*/

template <typename SolverT, typename SectionT>
void PoincareMap<SolverT,SectionT>::updateJacEnclosure(C2Set& S)
{
  if(this->derivativeOfFlow!=NULL)
  {
    intervalHull(*derivativeOfFlow,S.getLastMatrixEnclosure()*MatrixType(S),*derivativeOfFlow);

    // computation of the hessian of Poincare map
    if(hessianOfFlow!=NULL)
      intervalHull(*hessianOfFlow,S.getLastMatrixEnclosure()*HessianType(S) + S.getLastHessianEnclosure()*MatrixType(S),*hessianOfFlow);
  } 
  this->updateTime(S);
}

/*__________________________________________________________________________*/

template <typename SolverT, typename SectionT>
void PoincareMap<SolverT,SectionT>::saveJacEnclosure(C2Set& S)
{
  this->saveTime(S);
  if(derivativeOfFlow!=NULL)
    *derivativeOfFlow = MatrixType(S);
  if(hessianOfFlow!=NULL)
    *hessianOfFlow = (HessianType)S;
}

/*__________________________________________________________________________*/

template <typename SolverT, typename SectionT>
void PoincareMap<SolverT,SectionT>::saveJacEnclosure(C1Set& S)
{
  this->saveTime(S);
  if(derivativeOfFlow!=NULL)
    *derivativeOfFlow = MatrixType(S);
}

}} // namespace capd::poincare

#endif // _CAPD_POINCARE_POINCARE_MAP_HPP_

/// @}

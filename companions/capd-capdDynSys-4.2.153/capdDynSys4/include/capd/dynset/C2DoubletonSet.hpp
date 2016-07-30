/////////////////////////////////////////////////////////////////////////////
/// @file C2DoubletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C2DOUBLETONSET_HPP_
#define _CAPD_DYNSET_C2DOUBLETONSET_HPP_

#include <sstream>
#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/geomset/MatrixDoubletonSet.hpp"
#include "capd/dynset/C2DoubletonSet.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/diffAlgebra/Hessian.hpp"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace dynset{

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, ScalarType t)
  : SetType(
      x,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x),
    C1BaseSet(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : SetType(
      x+r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, r0),
    C1BaseSet(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(
      x+C*r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0),
    C1BaseSet(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(
      x+C*r0+r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0, r),
    C1BaseSet(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(
      const VectorType& x,
      const MatrixType& C, const VectorType& r0,
      const MatrixType& B, const VectorType& r,
      ScalarType t
   ) : SetType(
      x+C*r0+B*r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0, B, r),
    C1BaseSet(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t)
    : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       HessianType(c0part.dimension()),
       HessianType(c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part),
   m_Chess(c0part.dimension(),c0part.dimension()),
   m_Bhess(c0part.dimension(),c0part.dimension()),
   m_invBhess(c0part.dimension(),c0part.dimension()),
   m_alpha(c0part.dimension()),
   m_HR(c0part.dimension()),
   m_HR0(c0part.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(
      const C0BaseSet & c0part,
      const C1BaseSet& c1part,
      const HessianType& H,
      const  ScalarType t
    )  : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       HessianType(H),
       HessianType(c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part),
   m_Chess(c0part.dimension(),c0part.dimension()),
   m_Bhess(c0part.dimension(),c0part.dimension()),
   m_invBhess(c0part.dimension(),c0part.dimension()),
   m_alpha(H),
   m_HR(c0part.dimension()),
   m_HR0(c0part.dimension())
{
  split(m_alpha,m_HR0);
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

// ------------------------------------------------------------------------

template<typename MatrixType, class Policies>
void C2DoubletonSet<MatrixType,Policies>::move(DynSysType& c2dynsys, C2DoubletonSet& result)
{
  const size_type dim = this->m_currentSet.dimension();
  VectorType y(dim), rem(dim), enc(dim);
  MatrixType jacPhi(dim,dim), jacEnc(dim,dim), jacRem(dim,dim);
  MatrixType B(dim,dim), deltaC(dim,dim);

  HessianType EH(dim), LH(dim), RH(dim);

  VectorType xx = this->m_currentSet;
  MatrixType V = this->m_currentMatrix;
  c2dynsys.encloseC2Map(
      this->getCurrentTime(),
      this->m_x,xx,             // input arguments
      y, rem, enc,              // output for C^0 part
      jacPhi, jacRem, jacEnc,   // output for C^1 part
      LH, RH, EH                // output for C^2 part
  );

  result.m_x = y + rem;
  result.m_C = jacPhi * m_C;
  B = result.m_B = jacPhi * m_B;

  // ---------- C^0 part -----------------

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentSet = result.m_x + result.m_C * this->m_r0 + result.m_B*this->m_r;

  // here we compute representation for the new set
  // xx is unnecessary now
  split(result.m_x, xx);
  split(result.m_C, deltaC);
  xx += deltaC * m_r0;

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB
  this->Policies::computeBinvB(result.m_B,result.m_invB,this->m_r);

  // eventually we compute new representation of r
  result.m_r = (result.m_invB * B) * m_r + result.m_invB * xx;

  // ---------- C^1 part -----------------

  MatrixType J = jacPhi + jacRem;
  result.m_D = J*this->m_D;
  B = result.m_Bjac = J*this->m_Bjac;
  result.m_Cjac = J*this->m_Cjac;

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentMatrix = J*this->m_currentMatrix;
  if(!intersection(result.m_currentMatrix,result.m_D+result.m_Cjac*m_R0 + B*m_R,result.m_currentMatrix))
     throw std::logic_error("C2Doubleton::move error: empty intersection of two enclosures of monodromy matrix.");
  // jacRem is unnecessary now
  split(result.m_D, jacRem);
  split(result.m_Cjac, deltaC);
  jacRem += deltaC * m_R0;

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB
  this->Policies::computeBinvB(result.m_Bjac,result.m_invBjac,this->m_R);

  // eventually we compute new representation of r
  result.m_R = (result.m_invBjac * B) * m_R + result.m_invBjac * jacRem;

  // ---------- C^2 part: enclosure -----------------

  // compute LH += RH;
  capd::vectalg::addAssignObjectObject(LH,RH);

  // first enclosure - store it in deltaAplha
  // newHessian = JacPhi*oldHessian + alpha
  HessianType alpha = LH*V;
  HessianType deltaAlpha = J*this->m_currentHessian + alpha;

  // second possible enclosure
  // newHessian = JacPhi*m_alpha + (JacPhi*Bhess)*HR + (JacPhi*Chess)*HR0 + alpha

  result.m_Chess = J*this->m_Chess;
  B = result.m_Bhess = J*this->m_Bhess;

  result.m_currentHessian = alpha;
  result.m_currentHessian += J*this->m_alpha;
  result.m_currentHessian += result.m_Chess*this->m_HR0;
  result.m_currentHessian += B*this->m_HR;

  if(! capd::vectalg::intersection(deltaAlpha,result.m_currentHessian,result.m_currentHessian))
    throw std::logic_error("C2Doubleton::move error: empty intersection of two enclosures of Hessian.");

  // computation of new representation
  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB
  this->Policies::computeBinvB(result.m_Bhess,result.m_invBhess,this->m_HR);

  result.m_alpha = J*this->m_alpha;
  result.m_alpha += alpha;
  capd::vectalg::split(result.m_alpha,deltaAlpha);
  capd::vectalg::split(result.m_Chess,deltaC);

  deltaAlpha += deltaC*m_HR0;
  result.m_HR = (result.m_invBhess*B)*this->m_HR;
  result.m_HR += result.m_invBhess * deltaAlpha;

  if(&result != this){
    result.m_r0 = m_r0;
    result.m_R0 = m_R0;
    result.m_HR0 = m_HR0;
  }

  // save enclosures and new time
  result.setCurrentTime(this->getCurrentTime()+c2dynsys.getStep());
  result.setLastEnclosure(enc);
  result.setLastMatrixEnclosure(jacEnc);
  result.setLastHessianEnclosure(EH);

  this->Policies::reorganizeIfNeeded(result.m_B,result.m_invB,result.m_r,result.m_C,result.m_r0);
  this->Policies::reorganizeC1IfNeeded(result.m_Bjac,result.m_invBjac,result.m_R,result.m_Cjac,result.m_R0);
  this->Policies::reorganizeC2IfNeeded(result.m_Bhess,result.m_invBhess,result.m_HR,result.m_Chess,result.m_HR0);
}

// -------------------------------------------------------------

template<typename MatrixType, class Policies>
std::string C2DoubletonSet<MatrixType,Policies>::show(void) const
{
  std::ostringstream descriptor;
  descriptor << name()
             << C0BaseSet::toString()
             << C1BaseSet::toString();
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C2RECT2SET_HPP_

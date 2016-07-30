/////////////////////////////////////////////////////////////////////////////
/// @file C2DoubletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_C2DOUBLETONSET_H_
#define _CAPD_DYNSET_C2DOUBLETONSET_H_

#include "capd/dynset/C2Set.h"
#include "capd/dynsys/C2DynSys.h"
#include "capd/geomset/CenteredDoubletonSet.h"
#include "capd/geomset/MatrixDoubletonSet.h"


namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * C2 set in doubleton form.
 *
 *  C^2-Lohner algorithm.
 */

template<typename MatrixT, class Policies>
class C2DoubletonSet: public Policies, public C2Set<MatrixT>,
                      protected capd::geomset::CenteredDoubletonSet<MatrixT>,
                      protected capd::geomset::MatrixDoubletonSet<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C2Set<MatrixT> SetType;
  typedef typename SetType::HessianType HessianType;
  typedef capd::dynsys::C2DynSys<MatrixType> DynSysType;
  typedef capd::geomset::CenteredDoubletonSet<MatrixT> C0BaseSet;
  typedef capd::geomset::MatrixDoubletonSet<MatrixT> C1BaseSet;

  C2DoubletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, const HessianType& H, ScalarType t = TypeTraits<ScalarType>::zero());

  using SetType::operator VectorType;
  using SetType::operator MatrixType;
  using SetType::operator HessianType;

  void move(DynSysType& c2dynsys);
  void move(DynSysType& c2dynsys, C2DoubletonSet& result);

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r;
    VectorType x0 = midVector(this->m_currentSet);
    VectorType gradient = f.gradient(this->m_currentSet);
    if(!intersection(f(x0) + C0BaseSet::evalAffineFunctional(gradient,x0), f(this->m_currentSet), r)){
      throw std::logic_error("C2DoubletonSet::evalAt - empty intersection!. Report this error to CAPD developers!");
    }
    return r;
  }

  virtual std::string name() const { return "C2DoubletonSet"; }
  std::string show(void) const;

  using C0BaseSet::get_r0;
  using C0BaseSet::getElement_r0;
  using C0BaseSet::get_x;
  using C0BaseSet::getElement_x;
  using C0BaseSet::get_B;
  using C0BaseSet::get_invB;
  using C0BaseSet::getElement_B;
  using C0BaseSet::getRow_B;
  using C0BaseSet::getColumn_B;
  using C0BaseSet::get_C;
  using C0BaseSet::getElement_C;
  using C0BaseSet::getRow_C;
  using C0BaseSet::getColumn_C;

  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_Bjac;
  using C1BaseSet::m_R0;
  using C1BaseSet::m_Cjac;
  using C1BaseSet::get_D;
  using C1BaseSet::getElement_D;
  using C1BaseSet::getRow_D;
  using C1BaseSet::getColumn_D;
  using C1BaseSet::get_R;
  using C1BaseSet::getElement_R;
  using C1BaseSet::getRow_R;
  using C1BaseSet::getColumn_R;
  using C1BaseSet::get_invBjac;
  using C1BaseSet::get_Bjac;
  using C1BaseSet::getElement_Bjac;
  using C1BaseSet::getRow_Bjac;
  using C1BaseSet::getColumn_Bjac;
  using C1BaseSet::get_R0;
  using C1BaseSet::getElement_R0;
  using C1BaseSet::getRow_R0;
  using C1BaseSet::getColumn_R0;
  using C1BaseSet::get_Cjac;
  using C1BaseSet::getElement_Cjac;
  using C1BaseSet::getRow_Cjac;
  using C1BaseSet::getColumn_Cjac;

protected:
  using C0BaseSet::m_x;
  using C0BaseSet::m_r;
  using C0BaseSet::m_r0;
  using C0BaseSet::m_B;
  using C0BaseSet::m_C;

  using C1BaseSet::m_invBjac;
  // why does it not compile????
/*
  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_R0;
  using C1BaseSet::m_Bjac;
  using C1BaseSet::m_Cjac;
*/

// C^2 part represented as alpha + Bhess*HR + Chess*HR0

  MatrixType m_Chess, m_Bhess, m_invBhess;
  HessianType m_alpha, m_HR, m_HR0;
};
/// @}

// -----------------------------------------------------------------------------

template<typename MatrixType, class Policies>
inline void C2DoubletonSet<MatrixType,Policies>::move(DynSysType& c2dynsys)
{
  move(c2dynsys,*this);
}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C2DOUBLETONSET_H_

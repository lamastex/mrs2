/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file C1DynSys.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C1DYNSYS_H_
#define _CAPD_DYNSYS_C1DYNSYS_H_

#include <string>
#include "capd/dynsys/DynSys.h"

namespace capd{
namespace dynsys{

template<typename MatrixType>
class C1DynSys : public capd::dynsys::DynSys<MatrixType>{
public:
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;

  /// given an enclosure for the trajectories, computes enclosure for variational equations
  virtual MatrixType jacEnclosure(const ScalarType& t, const VectorType& enc) = 0;

  /// computes Lagrange remainder for variational equations
  /// given enclosures for the trajectories and for variational equations over the time step
  virtual void JacRemainder(
            const ScalarType& t,
            const VectorType &vecEnclosure,
            const MatrixType &jacEnclosure,
            VectorType &Remainder,
            MatrixType &jRemainder
      ) = 0;

  using DynSys<MatrixType>::enclosure;

  /// computes simultaneously an enclosure for the solutions to ODE and associated variational equations
  virtual void c1Enclosure(const ScalarType& t, const VectorType& x, VectorType& o_enc, MatrixType& o_jacEnc){
    o_enc = this->enclosure(t,x);
    o_jacEnc = this->jacEnclosure(t,o_enc);
  }

  virtual void encloseC1Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc
  ) = 0;
};

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_C1DYNSYS_H_

/// @}

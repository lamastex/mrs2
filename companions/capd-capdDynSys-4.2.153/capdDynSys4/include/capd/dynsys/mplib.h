//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynsys/mplib.h
///
/// @author Tomasz Kapela   @date 2010-01-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_MPLIB_H_
#define _CAPD_DYNSYS_MPLIB_H_


#include "capd/basicalg/factrial.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/C2Solver.h"
#include "capd/dynsys/CnSolver.h"
#include "capd/dynsys/SolverHOE.h"
#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"

#ifdef __HAVE_MPFR__

namespace capd{
  typedef capd::dynsys::DynSys<MpIMatrix> MpIDynSys;

  typedef capd::dynsys::C2Solver<MpIMap> MpIC2OdeSolver;
  typedef capd::dynsys::CnSolver<MpIMap> MpICnOdeSolver;
  typedef capd::dynsys::SolverHOE<capd::MpIMap> MpIOdeSolver;

  typedef capd::dynsys::SolverException<MpIVector> MpISolverException;

  // classes for nonrigorous computations
  typedef capd::dynsys::BasicSolver<MpMap, capd::dynsys::DLastTermsStepControl> MpOdeSolver;
  typedef capd::dynsys::BasicC2Solver<MpMap, capd::dynsys::DLastTermsStepControl> MpC2OdeSolver;
  typedef capd::dynsys::BasicCnSolver<MpMap, capd::dynsys::DLastTermsStepControl> MpCnOdeSolver;

///@deprecated
  typedef capd::dynsys::C2Solver<MpIMap> MpIC2Taylor;
  typedef capd::dynsys::CnSolver<MpIMap> MpICnTaylor;
  typedef capd::dynsys::SolverHOE<capd::MpIMap> MpITaylor;

  typedef capd::dynsys::SolverException<MpIVector> MpITaylorException;

  // classes for nonrigorous computations
  typedef capd::dynsys::BasicSolver<MpMap, capd::dynsys::DLastTermsStepControl> MpTaylor;
  typedef capd::dynsys::BasicC2Solver<MpMap, capd::dynsys::DLastTermsStepControl> MpC2Taylor;
  typedef capd::dynsys::BasicCnSolver<MpMap, capd::dynsys::DLastTermsStepControl> MpCnTaylor;

} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // _CAPD_DYNSYS_MPLIB_H_

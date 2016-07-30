
/////////////////////////////////////////////////////////////////////////////
/// @file Solver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"
#include "capd/dynsys/BasicSolver.hpp"
#include "capd/dynsys/Solver.hpp"
#include "capd/dynsys/SolverHOE.hpp"
#include "capd/vectalg/Matrix.hpp"

#ifdef __HAVE_MPFR__

//  template class capd::dynsys::Solver<capd::MpIMap, capd::dynsys::IMpLastTermsStepControl<capd::MpInterval> >;
  template class capd::dynsys::Solver<capd::MpIMap, capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::Solver<capd::MpIMap, capd::dynsys::IEncFoundStepControl>;

//  template class capd::dynsys::BasicSolver<capd::MpIMap, capd::dynsys::IMpLastTermsStepControl<capd::MpInterval> >;
  template class capd::dynsys::BasicSolver<capd::MpIMap, capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::BasicSolver<capd::MpIMap, capd::dynsys::IEncFoundStepControl>;
  template class capd::dynsys::BasicSolver<capd::MpIMap, capd::dynsys::DLastTermsStepControl>;

//  template class capd::dynsys::BasicSolver<capd::MpMap, capd::dynsys::MpLastTermsStepControl<capd::MpFloat> >;
  template class capd::dynsys::BasicSolver<capd::MpMap, capd::dynsys::DLastTermsStepControl>;

//  template class capd::dynsys::SolverHOE<capd::MpIMap, capd::dynsys::IMpLastTermsStepControl<capd::MpInterval> >;
  template class capd::dynsys::SolverHOE<capd::MpIMap, capd::dynsys::ILastTermsStepControl>;


#endif

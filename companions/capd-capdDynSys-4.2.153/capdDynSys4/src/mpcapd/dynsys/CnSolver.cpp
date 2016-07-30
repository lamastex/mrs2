
/////////////////////////////////////////////////////////////////////////////
/// @file CnSolver.cpp
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
#include "capd/dynsys/BasicCnSolver.hpp"
#include "capd/dynsys/CnSolver.hpp"
#include "capd/map/Map.hpp"

#ifdef __HAVE_MPFR__
/*
  template class capd::dynsys::BasicCnSolver<capd::MpMap, capd::dynsys::MpLastTermsStepControl<capd::MpFloat> >;

  template class capd::dynsys::BasicCnSolver<capd::MpIMap,capd::dynsys::IMpLastTermsStepControl<capd::MpInterval> >;
  template class capd::dynsys::CnSolver<capd::MpIMap,capd::dynsys::IMpLastTermsStepControl<capd::MpInterval>  >;
*/
  template class capd::dynsys::BasicCnSolver<capd::MpMap, capd::dynsys::DLastTermsStepControl>;
  template class capd::dynsys::BasicCnSolver<capd::MpIMap,capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::CnSolver<capd::MpIMap,capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::BasicCnSolver<capd::MpIMap,capd::dynsys::IEncFoundStepControl>;
  template class capd::dynsys::CnSolver<capd::MpIMap,capd::dynsys::IEncFoundStepControl>;
#endif



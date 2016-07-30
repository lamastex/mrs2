
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


#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/dynsys/Solver.hpp"

template class capd::dynsys::Solver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::Solver<capd::IMap,capd::dynsys::ILastTermsStepControl>;

template class capd::dynsys::BasicSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::BasicSolver<capd::IMap,capd::dynsys::ILastTermsStepControl>;

template class capd::dynsys::BasicSolver<capd::DMap>;
template class capd::dynsys::BasicSolver<capd::LDMap>;

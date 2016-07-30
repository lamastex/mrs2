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

#include "capd/map/lib.h"

#include "capd/dynsys/BasicCnSolver.hpp"
#include "capd/dynsys/CnSolver.hpp"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/diffAlgebra/CnTimeJet.h"
#include "capd/diffAlgebra/CnCurve.hpp"

template class capd::dynsys::BasicCnSolver<capd::DMap>;
template class capd::dynsys::BasicCnSolver<capd::LDMap>;
template class capd::dynsys::BasicCnSolver<capd::IMap,capd::dynsys::ILastTermsStepControl>;
template class capd::dynsys::CnSolver<capd::IMap,capd::dynsys::ILastTermsStepControl>;

template class capd::dynsys::BasicCnSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::CnSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;

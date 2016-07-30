//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynsys/typedefs.h
///
/// @author Daniel Wilczak   @date 2013-01-09
//
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

namespace CAPD_USER_NAMESPACE{

typedef capd::dynsys::DynSys<CAPD_USER_NAMESPACE::IMatrix> IDynSys;
typedef capd::dynsys::SolverHOE<CAPD_USER_NAMESPACE::IMap> IOdeSolver;
typedef capd::dynsys::C2Solver<CAPD_USER_NAMESPACE::IMap> IC2OdeSolver;
typedef capd::dynsys::CnSolver<CAPD_USER_NAMESPACE::IMap> ICnOdeSolver;
typedef capd::dynsys::SolverException<CAPD_USER_NAMESPACE::IVector> ISolverException;

// classes for nonrigorous computations

typedef capd::dynsys::BasicSolver<CAPD_USER_NAMESPACE::DMap> DOdeSolver;
typedef capd::dynsys::BasicC2Solver<CAPD_USER_NAMESPACE::DMap> DC2OdeSolver;
typedef capd::dynsys::BasicCnSolver<CAPD_USER_NAMESPACE::DMap> DCnOdeSolver;

typedef capd::dynsys::BasicSolver<CAPD_USER_NAMESPACE::LDMap> LDOdeSolver;
typedef capd::dynsys::BasicC2Solver<CAPD_USER_NAMESPACE::LDMap> LDC2OdeSolver;
typedef capd::dynsys::BasicCnSolver<CAPD_USER_NAMESPACE::LDMap> LDCnOdeSolver;


///@deprecated
typedef capd::dynsys::SolverHOE<CAPD_USER_NAMESPACE::IMap> ITaylor;
typedef capd::dynsys::C2Solver<CAPD_USER_NAMESPACE::IMap> IC2Taylor;
typedef capd::dynsys::CnSolver<CAPD_USER_NAMESPACE::IMap> ICnTaylor;
typedef capd::dynsys::SolverException<CAPD_USER_NAMESPACE::IVector> ITaylorException;

// classes for nonrigorous computations

typedef capd::dynsys::BasicSolver<CAPD_USER_NAMESPACE::DMap> DTaylor;
typedef capd::dynsys::BasicC2Solver<CAPD_USER_NAMESPACE::DMap> DC2Taylor;
typedef capd::dynsys::BasicCnSolver<CAPD_USER_NAMESPACE::DMap> DCnTaylor;

typedef capd::dynsys::BasicSolver<CAPD_USER_NAMESPACE::LDMap> LDTaylor;
typedef capd::dynsys::BasicC2Solver<CAPD_USER_NAMESPACE::LDMap> LDC2Taylor;
typedef capd::dynsys::BasicCnSolver<CAPD_USER_NAMESPACE::LDMap> LDCnTaylor;

} // end of namespace capd


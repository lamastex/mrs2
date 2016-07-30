
/////////////////////////////////////////////////////////////////////////////
/// @file Node.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#include "capd/intervals/lib.h"
#include "capd/diffAlgebra/Node.hpp"

namespace capd{
namespace diffAlgebra{

template class Node<double>;
template class ConsNode<double>;
template class VarNode<double>;
template class UnaryNode<double>;
template class BinaryNode<double>;

template class SumNode<double>;
template class DifNode<double>;
template class MulNode<double>;
template class DivNode<double>;
template class SinNode<double>;
template class CosNode<double>;
template class ExpNode<double>;
template class LogNode<double>;
template class SqrtNode<double>;
template class PowNode<double>;

template class MulConsNode<double>;
template class MulParamNode<double>;
template class MulParamParamNode<double>;
template class DivByParamNode<double>;
template class DivConsByParamNode<double>;

template class Node<long double>;
template class ConsNode<long double>;
template class VarNode<long double>;
template class UnaryNode<long double>;
template class BinaryNode<long double>;

template class SumNode<long double>;
template class DifNode<long double>;
template class MulNode<long double>;
template class DivNode<long double>;
template class SinNode<long double>;
template class CosNode<long double>;
template class ExpNode<long double>;
template class LogNode<long double>;
template class SqrtNode<long double>;
template class PowNode<long double>;

template class MulConsNode<long double>;
template class MulParamNode<long double>;
template class MulParamParamNode<long double>;
template class DivByParamNode<long double>;
template class DivConsByParamNode<long double>;

template class Node<capd::DInterval>;
template class ConsNode<capd::DInterval>;
template class UnaryNode<capd::DInterval>;
template class BinaryNode<capd::DInterval>;
template class VarNode<capd::DInterval>;
template class SumNode<capd::DInterval>;
template class DifNode<capd::DInterval>;
template class MulNode<capd::DInterval>;
template class DivNode<capd::DInterval>;
template class SinNode<capd::DInterval>;
template class CosNode<capd::DInterval>;
template class ExpNode<capd::DInterval>;
template class LogNode<capd::DInterval>;
template class SqrtNode<capd::DInterval>;
template class PowNode<capd::DInterval>;

template class MulConsNode<capd::DInterval>;
template class MulParamNode<capd::DInterval>;
template class MulParamParamNode<capd::DInterval>;
template class DivByParamNode<capd::DInterval>;
template class DivConsByParamNode<capd::DInterval>;

}}


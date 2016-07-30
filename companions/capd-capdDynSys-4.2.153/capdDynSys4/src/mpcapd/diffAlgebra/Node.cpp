
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

#include "capd/multiPrec/mplib.h"
#include "capd/intervals/mplib.h"
#include "capd/diffAlgebra/Node.hpp"

namespace capd{
namespace diffAlgebra{

template class Node<capd::MpFloat>;
template class ConsNode<capd::MpFloat>;
template class VarNode<capd::MpFloat>;
template class UnaryNode<capd::MpFloat>;
template class BinaryNode<capd::MpFloat>;

template class SumNode<capd::MpFloat>;
template class DifNode<capd::MpFloat>;
template class MulNode<capd::MpFloat>;
template class DivNode<capd::MpFloat>;
template class SinNode<capd::MpFloat>;
template class CosNode<capd::MpFloat>;
template class ExpNode<capd::MpFloat>;
template class LogNode<capd::MpFloat>;
template class SqrtNode<capd::MpFloat>;
template class PowNode<capd::MpFloat>;

template class MulConsNode<capd::MpFloat>;
template class MulParamNode<capd::MpFloat>;
template class MulParamParamNode<capd::MpFloat>;
template class DivByParamNode<capd::MpFloat>;
template class DivConsByParamNode<capd::MpFloat>;

template class Node<capd::MpInterval>;
template class ConsNode<capd::MpInterval>;
template class UnaryNode<capd::MpInterval>;
template class BinaryNode<capd::MpInterval>;
template class VarNode<capd::MpInterval>;
template class SumNode<capd::MpInterval>;
template class DifNode<capd::MpInterval>;
template class MulNode<capd::MpInterval>;
template class DivNode<capd::MpInterval>;
template class SinNode<capd::MpInterval>;
template class CosNode<capd::MpInterval>;
template class ExpNode<capd::MpInterval>;
template class LogNode<capd::MpInterval>;
template class SqrtNode<capd::MpInterval>;
template class PowNode<capd::MpInterval>;

template class MulConsNode<capd::MpInterval>;
template class MulParamNode<capd::MpInterval>;
template class MulParamParamNode<capd::MpInterval>;
template class DivByParamNode<capd::MpInterval>;
template class DivConsByParamNode<capd::MpInterval>;

}}


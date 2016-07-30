/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file eval.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_HPP_
#define _CAPD_AUTODIFF_EVAL_HPP_

#define CAPD_MAKE_NODE(NodeName,ClassName,left,right,result) case NodeName : { p = new ClassName##Node<T>; p->left = left; p->right = right; p->result = result; break; }

#include <stdexcept>
#include <sstream>

#include "capd/autodiff/DagIndexer.h"
#include "capd/autodiff/EvalAcos.h"
#include "capd/autodiff/EvalAdd.h"
#include "capd/autodiff/EvalAsin.h"
#include "capd/autodiff/EvalAtan.h"
#include "capd/autodiff/EvalDiv.h"
#include "capd/autodiff/EvalExp.h"
#include "capd/autodiff/EvalLog.h"
#include "capd/autodiff/EvalMul.h"
#include "capd/autodiff/EvalNaturalPow.h"
#include "capd/autodiff/EvalOneMinusSqr.h"
#include "capd/autodiff/EvalPow.h"
#include "capd/autodiff/EvalSinCos.h"
#include "capd/autodiff/EvalSqr.h"
#include "capd/autodiff/EvalSqrt.h"
#include "capd/autodiff/EvalSub.h"
#include "capd/autodiff/EvalUnaryMinus.h"

namespace capd{
namespace autodiff{

template<class T>
void Int4ToAbstractNode(const std::vector<MyNode>& node, std::vector<AbstractNode<T>* >& out, const JetSize jetSize, T* data)
{
  out.resize(node.size());
  for(unsigned i=0;i<node.size();++i)
  {
    AbstractNode<T>* p;
    T* left = data + jetSize*node[i].left;
    T* right = data + jetSize*node[i].right;
    T* result = data + jetSize*node[i].result;

    switch(node[i].op)
    {
    CAPD_MAKE_NODE(NODE_ADD,Add,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_VAR,ConstPlusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_CONST,ConstPlusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_TIME,ConstPlusTime,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_FUNTIME,ConstPlusFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_TIME_PLUS_VAR,TimePlusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_TIME_PLUS_FUNTIME,TimePlusFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_PLUS_VAR,FunTimePlusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_PLUS_FUNTIME,FunTimePlusFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_SUB,Sub,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_CONST,ConstMinusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_VAR,ConstMinusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_TIME,ConstMinusTime,left,right,result);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_FUNTIME,ConstMinusFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_CONST,TimeMinusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_FUNTIME,TimeMinusFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_VAR,TimeMinusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_CONST,FunTimeMinusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_TIME,FunTimeMinusTime,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_FUNTIME,FunTimeMinusFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_VAR,FunTimeMinusVar,left,right,result);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_CONST,VarMinusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_TIME,VarMinusTime,left,right,result);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_FUNTIME,VarMinusFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_UNARY_MINUS,UnaryMinus,left,right,result);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_CONST,UnaryMinusConst,left,right,result);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_TIME,UnaryMinusTime,left,right,result);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_FUNTIME,UnaryMinusFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_MUL,Mul,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_VAR,MulConstByVar,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_CONST,MulConstByConst,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_TIME,MulConstByTime,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_FUNTIME,MulConstByFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_TIME_BY_VAR,MulTimeByVar,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_TIME_BY_FUNTIME,MulTimeByFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_FUNTIME_BY_VAR,MulFunTimeByVar,left,right,result);
    CAPD_MAKE_NODE(NODE_MUL_FUNTIME_BY_FUNTIME,MulFunTimeByFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_DIV,Div,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_CONST,DivVarByConst,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_TIME,DivVarByTime,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_FUNTIME,DivVarByFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_TIME_BY_CONST,DivTimeByConst,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_CONST,DivFunTimeByConst,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_TIME,DivFunTimeByTime,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_FUNTIME,DivFunTimeByFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_DIV_CONST_BY_CONST,DivConstByConst,left,right,result);

    CAPD_MAKE_NODE(NODE_SQR,Sqr,left,right,result);
    CAPD_MAKE_NODE(NODE_SQR_CONST,SqrConst,left,right,result);
    CAPD_MAKE_NODE(NODE_SQR_TIME,SqrTime,left,right,result);
    CAPD_MAKE_NODE(NODE_SQR_FUNTIME,SqrFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_SQRT,Sqrt,left,right,result);
    CAPD_MAKE_NODE(NODE_SQRT_CONST,SqrtConst,left,right,result);
    CAPD_MAKE_NODE(NODE_SQRT_TIME,SqrtTime,left,right,result);
    CAPD_MAKE_NODE(NODE_SQRT_FUNTIME,SqrtFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_POW,Pow,left,right,result);
    CAPD_MAKE_NODE(NODE_POW_CONST,PowConst,left,right,result);
    CAPD_MAKE_NODE(NODE_POW_TIME,PowTime,left,right,result);
    CAPD_MAKE_NODE(NODE_POW_FUNTIME,PowFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_NATURAL_POW,NaturalPow,left,right,result);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_CONST,NaturalPowConst,left,right,result);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_TIME,NaturalPowTime,left,right,result);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_FUNTIME,NaturalPowFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_INTEGER_POW,IntegerPow,left,right,result);
    CAPD_MAKE_NODE(NODE_INTEGER_POW_CONST,IntegerPowConst,left,right,result);
    CAPD_MAKE_NODE(NODE_INTEGER_POW_TIME,IntegerPowTime,left,right,result);
    CAPD_MAKE_NODE(NODE_INTEGER_POW_FUNTIME,IntegerPowFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_HALF_INTEGER_POW,HalfIntegerPow,left,right,result);
    CAPD_MAKE_NODE(NODE_HALF_INTEGER_POW_CONST,HalfIntegerPowConst,left,right,result);
    CAPD_MAKE_NODE(NODE_HALF_INTEGER_POW_TIME,HalfIntegerPowTime,left,right,result);
    CAPD_MAKE_NODE(NODE_HALF_INTEGER_POW_FUNTIME,HalfIntegerPowFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_EXP,Exp,left,right,result);
    CAPD_MAKE_NODE(NODE_EXP_CONST,ExpConst,left,right,result);
    CAPD_MAKE_NODE(NODE_EXP_TIME,ExpTime,left,right,result);
    CAPD_MAKE_NODE(NODE_EXP_FUNTIME,ExpFunTime,left,right,result);
    CAPD_MAKE_NODE(NODE_LOG,Log,left,right,result);
    CAPD_MAKE_NODE(NODE_LOG_CONST,LogConst,left,right,result);
    CAPD_MAKE_NODE(NODE_LOG_TIME,LogTime,left,right,result);
    CAPD_MAKE_NODE(NODE_LOG_FUNTIME,LogFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_SIN,Sin,left,right,result);
    CAPD_MAKE_NODE(NODE_SIN_CONST,SinConst,left,right,result);
    CAPD_MAKE_NODE(NODE_SIN_TIME,SinTime,left,right,result);
    CAPD_MAKE_NODE(NODE_SIN_FUNTIME,SinFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR,OneMinusSqr,left,right,result);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_CONST,OneMinusSqrConst,left,right,result);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_TIME,OneMinusSqrTime,left,right,result);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_FUNTIME,OneMinusSqrFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_ATAN,Atan,left,right,result);
    CAPD_MAKE_NODE(NODE_ATAN_CONST,AtanConst,left,right,result);
    CAPD_MAKE_NODE(NODE_ATAN_TIME,AtanTime,left,right,result);
    CAPD_MAKE_NODE(NODE_ATAN_FUNTIME,AtanFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_ASIN,Asin,left,right,result);
    CAPD_MAKE_NODE(NODE_ASIN_CONST,AsinConst,left,right,result);
    CAPD_MAKE_NODE(NODE_ASIN_TIME,AsinTime,left,right,result);
    CAPD_MAKE_NODE(NODE_ASIN_FUNTIME,AsinFunTime,left,right,result);

    CAPD_MAKE_NODE(NODE_ACOS,Acos,left,right,result);
    CAPD_MAKE_NODE(NODE_ACOS_CONST,AcosConst,left,right,result);
    CAPD_MAKE_NODE(NODE_ACOS_TIME,AcosTime,left,right,result);
    CAPD_MAKE_NODE(NODE_ACOS_FUNTIME,AcosFunTime,left,right,result);

    default:
      std::ostringstream out;
      out << "Implementation error! Node no. " << node[i].op << " is missed in switch-case block in function evaluating DAG.\nPlease report this bug to developers!\n ";
      throw std::logic_error(out.str());
    }
    out[i] = p;
  }
}

}} // namespace capd::autodiff

#endif

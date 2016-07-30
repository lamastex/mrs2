/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Node_explog.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-20012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_NODE_EXPLOG_HPP_
#define _CAPD_DIFFALGEBRA_NODE_EXPLOG_HPP_

namespace capd{
namespace diffAlgebra{

template<typename ScalarType>
ScalarType& ExpNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  (*(this->left))(i);

  if(i==0)
    return this->value[i] = exp(this->left->value[0]);

  ScalarType result(0);
  for(int j=0;j<i;++j)
    result += ScalarType(i-j) * this->value[j] * this->left->value[i-j];

  return this->value[i] = result/ScalarType(i);
}

template<typename ScalarType>
ScalarType& LogNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  (*(this->left))(i);

  if(i==0) 
    return this->value[i] = log(this->left->value[0]);

  ScalarType result(0.);
  for(int j=0;j<i;++j)
    result+= ScalarType(i-j) * this->left->value[j] * this->value[i-j];

  return this->value[i] = (this->left->value[i] - result/ScalarType(i)) / this->left->value[0];
}

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_NODE_EXPLOG_HPP_

/// @}

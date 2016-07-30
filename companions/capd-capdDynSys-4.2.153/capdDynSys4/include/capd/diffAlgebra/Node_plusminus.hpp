/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Node_plusminus.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_NODE_PLUSMINUS_HPP_
#define _CAPD_DIFFALGEBRA_NODE_PLUSMINUS_HPP_

namespace capd{
namespace diffAlgebra{

template<typename ScalarType>
ScalarType& SumNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  return this->value[i] = (*(this->left))(i) + (*(this->right))(i);
}

template<typename ScalarType>
ScalarType& DifNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  return this->value[i] = (*(this->left))(i) - (*(this->right))(i);
}

// --------------------------------- global Operators -------------------------

template<typename ScalarType>
capd::diffAlgebra::Node<ScalarType>& operator+(capd::diffAlgebra::Node<ScalarType>& x, capd::diffAlgebra::Node<ScalarType>& y)
{
  if(x.getOrder()!=y.getOrder())
    throw std::runtime_error("operator+(Node&, Node&) - incompatible dimensions");
  return *(new capd::diffAlgebra::SumNode<ScalarType>(x.getOrder(), &x,&y));
}

template<typename ScalarType>
capd::diffAlgebra::Node<ScalarType>& operator-(capd::diffAlgebra::Node<ScalarType>& x, capd::diffAlgebra::Node<ScalarType>& y)
{
  if(x.getOrder()!=y.getOrder())
    throw std::runtime_error("operator-(Node&, Node&) - incompatible dimensions");
  return *(new capd::diffAlgebra::DifNode<ScalarType>(x.getOrder(),&x,&y));
}

}} // namespace capd::diffAlgebra


#endif // _CAPD_DIFFALGEBRA_NODE_PLUSMINUS_HPP_

/// @}

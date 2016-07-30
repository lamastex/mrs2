/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Node_muldiv.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_NODE_MULDIV_HPP_
#define _CAPD_DIFFALGEBRA_NODE_MULDIV_HPP_

namespace capd{
namespace diffAlgebra{


template<typename ScalarType>
ScalarType& MulNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  (*(this->left))(i);
  (*(this->right))(i);

  ScalarType result(0);
  for(int j=0;j<=i;++j)
    result += (this->left->value[j]) * (this->right->value[i-j]) ;

  return this->value[i] = result;
}


template<typename ScalarType>
ScalarType& MulConsNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  return this->value[i] = (*(this->left))(i) * m_constant;
}

template<typename ScalarType>
ScalarType& MulParamNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  return this->value[i] = (*(this->left))(i) * (*(this->right->value));
}

template<typename ScalarType>
ScalarType& MulParamParamNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;
  if(i>0)
    return this->value[i] = ScalarType(0.);
  
  return *(this->value) = (*(this->left->value)) * (*(this->right->value));
}

template<typename ScalarType>
ScalarType& MulConsParamNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;
  if(i>0)
    return this->value[i] = ScalarType(0.);
  
  return *(this->value) = *(this->left->value) * m_constant;
}

// -----------------------------------------------------------------------


template<typename ScalarType>
ScalarType& DivNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  (*(this->left))(i);
  (*(this->right))(i);
  if (! (this->right->value[0]!=ScalarType(0.) ) )
  {
    throw(std::runtime_error("function error: division by zero"));
  }

  ScalarType result(0.);
  for(int j=0;j<i;++j)
    result+=this->value[j] * (this->right->value[i-j]);

  return this->value[i] = (this->left->value[i] - result) / (this->right->value[0]);
}

template<typename ScalarType>
ScalarType& DivByParamNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  return this->value[i] = (*(this->left))(i) / (*(this->right->value));
}

template<typename ScalarType>
ScalarType& DivConsByParamNode<ScalarType>::operator()(int i)
{
  if(i<=this->m_maxComputedDerivative)
    return this->value[i];
  this->m_maxComputedDerivative = i;

  if(i>0)
    return this->value[i] = ScalarType(0.);

  return *(this->value) = m_constant / (*(this->left->value));
}
// -----------------------------------------------------------------------


template<typename ScalarType>
capd::diffAlgebra::Node<ScalarType>& operator*(capd::diffAlgebra::Node<ScalarType>& x, capd::diffAlgebra::Node<ScalarType>& y)
{
  if(x.getOrder()!=y.getOrder())
    throw std::runtime_error("operator*(Node&, Node&) - incompatible dimensions");
  return *(new capd::diffAlgebra::MulNode<ScalarType>(x.getOrder(),&x,&y));
}

template<typename ScalarType>
capd::diffAlgebra::Node<ScalarType>& operator/(capd::diffAlgebra::Node<ScalarType>& x, capd::diffAlgebra::Node<ScalarType>& y)
{
  if(x.getOrder()!=y.getOrder())
    throw std::runtime_error("operator/(Node&, Node&) - incompatible dimensions");
  return *(new capd::diffAlgebra::DivNode<ScalarType>(x.getOrder(),&x,&y));
}

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_NODE_MULDIV_HPP_

/// @}

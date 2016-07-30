/////////////////////////////////////////////////////////////////////////////
/// @file DagIndexer.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_DAGINDEXER_HPP_
#define _CAPD_AUTODIFF_DAGINDEXER_HPP_

#include <algorithm>
#include "capd/autodiff/DagIndexer.h"

namespace capd{
namespace autodiff{

template<class Scalar>
DagIndexer<Scalar>::DagIndexer(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
  : m_order(order)
{
  this->allocate(domain,image,degree,nodes,order);
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>::DagIndexer(const DagIndexer& dag)
  : m_order(dag.m_order)
{
  this->allocateAndCopy(dag);
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>& DagIndexer<Scalar>::operator=(const DagIndexer& dag)
{
  if(this == &dag) return *this;
  delete[] m_coefficients;
  this->allocateAndCopy(dag);
  return *this;
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>::~DagIndexer()
{
  delete []m_coefficients;
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::setOrder(Order order)
{
  size_type newTimeJetSize = (order+1)*this->jetSize();
  size_type blockSize = newTimeJetSize*m_numberOfNodes;
  ScalarType* coeff = new ScalarType[blockSize];
  std::fill(coeff,coeff+blockSize,TypeTraits<ScalarType>::zero());

  for(size_type i=0;i<this->m_numberOfNodes;++i)
    getC0Coeff(coeff,VarNo(i),JetSize(newTimeJetSize),CoeffNo(0)) = getC0Coeff(m_coefficients,VarNo(i),JetSize(m_timeJetSize),CoeffNo(0));
  delete[] m_coefficients;
  m_coefficients = coeff;
  this->m_order = Order(order);
  this->m_timeJetSize = newTimeJetSize;
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::resize(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
{
  delete[] m_coefficients;
  this->allocate(domain,image,degree,nodes,order);
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::allocate(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
{
  this->m_domainDimension = domain;
  this->m_imageDimension = image;
  this->m_degree = degree;
  this->m_numberOfNodes = nodes;
  this->m_order = Order(order);
  this->m_timeJetSize = (this->m_order+1)*this->jetSize();
  size_type blockSize = this->m_numberOfNodes*this->m_timeJetSize;
  this->m_coefficients = new ScalarType[blockSize];
  std::fill(this->m_coefficients,this->m_coefficients+blockSize,TypeTraits<ScalarType>::zero());
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::allocateAndCopy(const DagIndexer& dag)
{
  this->m_domainDimension = dag.m_domainDimension;
  this->m_imageDimension = dag.m_imageDimension;
  this->m_degree = dag.m_degree;
  this->m_numberOfNodes = dag.m_numberOfNodes;
  this->m_order = dag.m_order;
  this->m_timeJetSize = dag.m_timeJetSize;

  size_type blockSize = this->m_numberOfNodes*this->m_timeJetSize;
  this->m_coefficients = new ScalarType[blockSize];
  std::copy(dag.coefficients(),dag.coefficients()+blockSize,this->m_coefficients);
}

}}

#endif

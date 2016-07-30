/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Map.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_MAP_HPP_ 
#define _CAPD_MAP_MAP_HPP_ 

#include <algorithm>
#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/autodiff/DagIndexer.hpp"
#include "capd/map/BasicFunction.hpp"
#include "capd/autodiff/eval.hpp"
#include "capd/map/Map.h"
#include "capd/basicalg/factrial.h"
#include "capd/diffAlgebra/Hessian.hpp"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/map/Map_codeTranslation.hpp"

namespace capd{
namespace map{

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>::Map(void) : m_degree(1)
{
  this->realloc(this->m_degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>::Map(const std::string &f, size_type degree) : BaseFunction(f), m_degree(degree)
{
  this->realloc(this->m_degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>::Map(const char *s, size_type degree) : BaseFunction(s), m_degree(degree)
{
  this->realloc(this->m_degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>::Map(const Map& f) : BaseFunction(f), m_degree(f.degree())
{
  this->realloc(this->m_degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>& Map<MatrixT>::operator=(const std::string &s)
{
  try{
    BaseFunction::operator=(s);
    this->realloc(this->m_degree);
  } catch(std::runtime_error &r)
  {
    BaseFunction::createDefault();
    this->createDefault();
    std::string re = "exception in Map &Map::operator=(const string &)\n";
    re += r.what();
    throw std::runtime_error(re);
  }
  return *this;
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>& Map<MatrixT>::operator=(const char *s)
{
  return this->operator=(std::string(s));
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
Map<MatrixT>& Map<MatrixT>::operator=(const Map& f)
{
  if(&f==this) return *this;
  BaseFunction::operator=(f);
  this->realloc(this->m_degree);
  return *this;
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(const VectorType& v) const
{
  using namespace capd::autodiff;
  const size_type d = v.dimension();
  if(d > (size_type)this->m_var.size())
  {
    std::ostringstream message;
    message << "Exception in Map::operator(const VectorType &,int)\n";
    message << "VectorType has " << d << " coordinates\n";
    message << "Map object has " << this->m_var.size()<< " variables";
    throw std::runtime_error(message.str());
  }

  size_type i;
  for(i=0;i<d;++i)
    this->m_dag(VarNo(i),CoeffNo(0)) = v[i];

  this->evalHomogenousPolynomial();

  ImageVectorType result(this->imageDimension());
  for(i=0;i<this->imageDimension();++i)
    result[i] = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(0));
  return result;
}

// -----------------------------------------------------------------------

template<typename MatrixT>
typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(const VectorType& u, MatrixType& der) const
{
  this->checkDegree(1);
  using namespace capd::autodiff;
  size_type i,j;
  for(i=0;i<this->dimension();++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = u[i];
    for(j=0;j<this->dimension();++j)
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();
  }

  this->evalHomogenousPolynomial();
  this->evalHomogenousPolynomial(1);

  ImageVectorType result(this->imageDimension());
  for(i=0;i<this->imageDimension();++i)
  {
    result[i] = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(0));
    for(j=0;j<this->dimension();++j)
      der(i+1,j+1) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),CoeffNo(0));
  }
  return result;
}

// -----------------------------------------------------------------------

template<typename MatrixT>
typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(const VectorType& u, MatrixType& der, HessianType& h) const
{
  this->checkDegree(2);
  using namespace capd::autodiff;
  size_type i,j,c;
  const size_type dim = this->dimension();
  for(i=0;i<dim;++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = u[i];
    for(j=0;j<dim;++j)
    {
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();
      for(c=j;c<dim;++c)
        this->m_dag(VarNo(i),DerNo(j),DerNo(c),CoeffNo(0)) = TypeTraits<ScalarType>::zero();
    }
  }

  this->evalHomogenousPolynomial();
  this->evalHomogenousPolynomial(1);
  this->evalHomogenousPolynomial(2);

  ImageVectorType result(this->imageDimension());
  for(i=0;i<this->imageDimension();++i)
  {
    result[i] = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(0));
    for(j=0;j<dim;++j)
    {
      der(i+1,j+1) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),CoeffNo(0));
      for(c=j;c<dim;++c)
        h(i,j,c) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),DerNo(c),CoeffNo(0));
    }
  }
  return result;
}

// -----------------------------------------------------------------------

template<typename MatrixT>
typename Map<MatrixT>::JetType
Map<MatrixT>::operator()(const JetType& c) const{
  this->checkDegree(c.degree());
  JetType result(this->imageDimension(),this->dimension(),c.degree());

  // first we set initial condition
  for(size_type i=0;i<c.dimension();++i)
  {
    typename JetType::const_iterator b = c.begin(i), e=c.end(i);
    ScalarType* p = this->m_dag.coefficients()+i*(this->m_dag.timeJetSize());
    while(b!=e)
    {
      *p = *b;
      b++;
      p += this->getOrder()+1;
    }
  }

  this->evalAndCopyResult(result);

  return result;
}


// -----------------------------------------------------------------------

template<typename MatrixT>
typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(const VectorType& x, JetType& c) const
{
  using namespace capd::autodiff;
  this->checkDegree(c.degree());

  // first we set initial condition

  size_type i,j;
  for(i=0;i<this->dimension();++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = x[i];
    for(j=0;j<this->dimension();++j)
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();

    typename JetType::const_iterator b = c.begin(i,2), e=c.end(i);
    ScalarType* p = this->m_dag.coefficients()+i*(this->m_dag.timeJetSize()) + (1+this->dimension())*(this->getOrder()+1);
    for(;b!=e;++b,p+=this->getOrder()+1)
      *p = TypeTraits<ScalarType>::zero();
  }

  this->evalAndCopyResult(c);

  return c();
}

// -----------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::homogenousPolynomial(MatrixType& der) const
{
  this->checkDegree(1);
  using namespace capd::autodiff;
  size_type i,j;
  const size_type dim = this->dimension();
  const size_type outdim = this->imageDimension();
  for(i=0;i<dim;++i)
    for(j=0;j<dim;++j)
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();

  this->evalHomogenousPolynomial(1);

  for(i=0;i<outdim;++i)
    for(j=0;j<dim;++j)
      der(i+1,j+1) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),CoeffNo(0));
}

// -----------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::homogenousPolynomial(const MatrixType& der, HessianType& h) const
{
  this->checkDegree(2);
  using namespace capd::autodiff;
  size_type i,j,c;
  const size_type dim = this->dimension();
  const size_type outdim = this->imageDimension();
  for(i=0;i<dim;++i)
    for(j=0;j<dim;++j)
    {
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) = der(i+1,j+1);
      for(c=j;c<dim;++c)
        this->m_dag(VarNo(i),DerNo(j),DerNo(c),CoeffNo(0)) = TypeTraits<ScalarType>::zero();
    }

  this->evalHomogenousPolynomial(1);
  this->evalHomogenousPolynomial(2);

  for(i=0;i<outdim;++i)
    for(j=0;j<dim;++j)
      for(c=j;c<dim;++c)
        h(i,j,c) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),DerNo(c),CoeffNo(0));
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::computeODECoefficients(VectorType iv[], size_type order) const
{
  this->checkOrder(order);

  using namespace capd::autodiff;
  size_type i,k;
  for(i=0;i<iv->dimension();++i)
    this->m_dag(VarNo(i),CoeffNo(0)) = iv[0][i];

  for(k=0;k<order;++k)
  {
    this->eval(CoeffNo(k));
    for(i=0;i<(size_type)iv->dimension();++i)
      iv[k+1][i] = this->m_dag(VarNo(i),CoeffNo(k+1))
                 = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(k))/(k+1);
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::computeODECoefficients(VectorType iv[], MatrixType im[], size_type order) const
{
  this->checkOrder(order);
  this->checkDegree(1);

  using namespace capd::autodiff;
  size_type i,j,k;
  const size_type dim = this->dimension();
  for(i=0;i<dim;++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = iv[0][i];
    for(j=0;j<dim;++j)
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) = im[0](i+1,j+1);
  }

  for(k=0;k<order;++k)
  {
    this->eval(1,CoeffNo(k));
    for(i=0;i<dim;++i)
    {
      iv[k+1][i] = this->m_dag(VarNo(i),CoeffNo(k+1)) = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(k))/(k+1);
      for(j=0;j<dim;++j)
        im[k+1](i+1,j+1) = this->m_dag(VarNo(i),DerNo(j),CoeffNo(k+1)) = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),CoeffNo(k))/(k+1);
    }
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::computeODECoefficients(VectorType iv[], MatrixType im[], HessianType h[], size_type order) const
{
  this->checkOrder(order);
  this->checkDegree(2);
  const size_type dim = this->dimension();
  using namespace capd::autodiff;
  size_type i,j,c,k;
  for(i=0;i<dim;++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = iv[0][i];
    for(j=0;j<dim;++j){
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) = im[0](i+1,j+1);
      for(c=j;c<dim;++c)
        this->m_dag(VarNo(i),DerNo(j),DerNo(c),CoeffNo(0)) = h[0](i,j,c);
    }
  }

  for(k=0;k<order;++k)
  {
    this->eval(2,CoeffNo(k));
    for(i=0;i<dim;++i)
    {
      iv[k+1][i] = this->m_dag(VarNo(i),CoeffNo(k+1))
                 = this->m_dag(VarNo(this->m_pos[i]),CoeffNo(k))/(k+1);

      for(j=0;j<dim;++j)
      {
        im[k+1](i+1,j+1) = this->m_dag(VarNo(i),DerNo(j),CoeffNo(k+1))
                         = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),CoeffNo(k))/(k+1);

        for(c=j;c<dim;++c)
        {
          h[k+1](i,j,c) = this->m_dag(VarNo(i),DerNo(j),DerNo(c),CoeffNo(k+1))
                        = this->m_dag(VarNo(this->m_pos[i]),DerNo(j),DerNo(c),CoeffNo(k))/(k+1);
        }
      }
    }
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::computeODECoefficients(JetType x[], size_type degree, size_type order) const
{
  this->checkOrder(order);
  this->checkDegree(degree);

  using namespace capd::autodiff;
  size_type i;

  for(i=0;i<this->dimension();++i)
  {
    typename JetType::const_iterator b = x[0].begin(i), e=x[0].end(i,degree);
    ScalarType* p = this->m_dag.coefficients()+i*(this->m_dag.timeJetSize());
    while(b!=e)
    {
      *p = *b;
      b++;
      p += this->getOrder()+1;
    }
  }

  for(size_type k=0;k<order;++k)
  {
    this->eval(degree,CoeffNo(k));
    for(i=0;i<this->dimension();++i)
    {
      typename JetType::iterator b = x[k+1].begin(i), e = x[k+1].end(i,degree);
      ScalarType* p = this->m_dag.coefficients() + this->m_pos[i]*(this->m_dag.timeJetSize())+k;
      ScalarType* q = this->m_dag.coefficients() + i*(this->m_dag.timeJetSize())+k+1;
      while(b!=e)
      {
        *b = *q = (*p)/(k+1);
        b++;
        p += this->getOrder()+1;
        q += this->getOrder()+1;
      }
    }
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::evalAndCopyResult(JetType& c) const
{
  // propagate jet
  for(size_type d=0;d<=c.degree();++d)
    this->evalHomogenousPolynomial(d);
  for(size_type i=0;i<this->imageDimension();++i)
  {
    typename JetType::iterator b = c.begin(i), e=c.end(i);
    ScalarType* p = this->m_dag.coefficients() + this->m_pos[i]*(this->m_dag.timeJetSize());
    while(b!=e)
    {
      *b = *p;
      b++;
      p += this->getOrder()+1;
    }
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::setParameters(const VectorType& values)
{
  BaseFunction::setParameters(values.begin(),values.dimension());
}

// -----------------------------------------------------------------------

template<typename MatrixT>
void Map<MatrixT>::setDegree(size_type degree)
{
  this->m_degree = degree;
  this->realloc(this->m_degree);
}

}} // namespace capd::map

#endif // _CAPD_MAP_MAP_HPP_ 

/// @}

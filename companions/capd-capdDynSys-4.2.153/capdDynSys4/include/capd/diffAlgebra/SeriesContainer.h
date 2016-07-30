/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file SeriesContainerz.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_SERIESCONTAINER_H_
#define _CAPD_DIFFALGEBRA_SERIESCONTAINER_H_

#include "capd/diffAlgebra/CnContainer.h"
#include "capd/diffAlgebra/Node.h"

namespace capd {
namespace diffAlgebra {

template <typename ScalarT>
class SeriesContainer : public capd::diffAlgebra::CnContainer<Node<ScalarT>*,0,0,0> {
public:
  typedef ScalarT ScalarType;
  typedef Node<ScalarType> NodeType;
  typedef capd::diffAlgebra::CnContainer<NodeType*,0,0,0> ContainerType;
  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;
  typedef typename ContainerType::Multipointer Multipointer;
  typedef typename ContainerType::Multiindex Multiindex;

  using ContainerType::begin;
  using ContainerType::end;

  SeriesContainer(int dim, int degree, int order) :
    ContainerType(dim, dim, degree), m_order(order) {
    createNodes(order);
  }

  SeriesContainer(int dim, int degree, bool) :
    ContainerType(dim, dim, degree), m_order(-1) {
    iterator b = begin(), e = end();
    while(b != e) {
      (*b) = 0;
      ++b;
    }
  }

  ~SeriesContainer() {
    removeNodes();
  }

  void setOrder(int a_new) {
    m_order = a_new;
    iterator b = begin(), e = end();
    while(b != e) {
      if(*b)
        (*b)->setOrder(a_new + 1, NULL);
      ++b;
    }
  }
  template <typename VectorType>
  void takeVector(const Multipointer& mp, int r, VectorType& result) const {
    typename VectorType::iterator b = result.begin(), e = result.end();
    int ind = this->index(mp);
    int d = 0;
    while(b != e && d < this->dimension()) {
      const_iterator n = begin(d, mp.dimension()) + ind;
      *b = (*n)->value[r];
      ++b;
      ++d;
    }
  }
protected:
  int m_order;  // order of taylor expansion in time variable
  void createNodes(int order) {
    iterator b = begin(), e = end();
    while(b != e) {
      (*b) = new ConsNode<ScalarType> (order + 1, ScalarType(0.));
      ++((**b).m_links);
      ++b;
    }
  }

  void removeNodes() {
    iterator b = begin(), e = end();
    while(b != e) {
      if(*b)
        if(!(--(**b).m_links))
          delete (*b);
      ++b;
    }
  }
private:
  /// we do not allow copying of these Containers.
  SeriesContainer(const SeriesContainer & s) {}

  /// we do not allow copying of these Containers.
  SeriesContainer& operator=(const SeriesContainer & s) {}
};

}} // the end of the namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_SERIESCONTAINER_H_
/// @}

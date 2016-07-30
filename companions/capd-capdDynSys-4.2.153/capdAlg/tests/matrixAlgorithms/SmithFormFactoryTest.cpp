/////////////////////////////////////////////////////////////////////////////
/// @file SmithFormFactoryTest.cpp
///
/// @author Mateusz Juda <mateusz.juda@{ii.uj.edu.pl,gmail.com}>
///
/// @date 2014-06-19
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include <capd/matrixAlgorithms/SmithFormFactory.h>
#include <capd/matrixAlgorithms/PARIInterface.h>
#include <capd/matrixAlgorithms/PARISmithForm.h>
#include <capd/auxil/Logger.h>
#include <capd/multiPrec/MpInt.h>

#include "SmithFixtures.h"


#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <sstream>
#include <memory>
#include <typeinfo>

namespace bl = boost::lambda;
namespace b = boost;

using namespace boost::assign;
using namespace boost::lambda;

using namespace capd::matrixAlgorithms;
using namespace capd::test;

BOOST_AUTO_TEST_SUITE(SmithFormFactoryTestSuite)

BOOST_AUTO_TEST_CASE_TEMPLATE(defaultSmithForm, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType Matrix;
  Fixture fixture;
  SmithFormFactory factory(false);

  for (size_t i = 0; i < fixture.matrices_matrix.size(); ++i) {
    Matrix B = fixture.matrices_matrix[i];
    std::auto_ptr<SmithForm<Matrix> > smithForm(factory(B, true, true, true, true));

    (*smithForm)();

    BOOST_CHECK_EQUAL(fixture.matrices_smithMatrix[i], B);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PARISmithFormTest, Fixture, SmithFixtures)
{
  if (!PARIInterface::enabled()) {
    return;
  }

  typedef typename Fixture::MatrixType Matrix;
  Fixture fixture;
  SmithFormFactory factory(true);

  for (size_t i = 0; i < fixture.matrices_matrix.size(); ++i) {
    Matrix B = fixture.matrices_matrix[i];
    std::auto_ptr<SmithForm<Matrix> > smithForm(factory(B, true, true, true, true));

    BOOST_ASSERT(dynamic_cast<PARISmithForm<typename Matrix::ScalarType>*>(smithForm.get()) != NULL);
    try{
      (*smithForm)();
      BOOST_CHECK_EQUAL(fixture.matrices_smithMatrix[i], B);
    } catch (std::overflow_error& ex) {
      if (!fixture.bigInt)
	throw ex;
    }
  }
}


#ifdef HAVE_MPCAPD_ALG
typedef capd::vectalg::Matrix<long long, 0, 0> ZMatrix;
typedef capd::vectalg::Matrix<capd::multiPrec::MpInt, 0, 0> MpMatrix;

MpMatrix promote(const ZMatrix& orgData)
{
  const size_t rows = orgData.numberOfRows();
  const size_t cols = orgData.numberOfColumns();

  MpMatrix promoted(rows, cols);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      ZMatrix::ScalarType scalar = orgData[i][j];
      promoted[i][j] = (long) scalar; // for now it is enough, there is no support for long long.
      //      mpz_import(promoted[i][j].get_mpz_t(), 1, -1, sizeof(scalar), 0, 0, &scalar);
    }
  }

  return promoted;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(defaultSmithFormMp, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType Matrix;
  Fixture fixture;
  SmithFormFactory factory(false);

  for (size_t i = 0; i < fixture.matrices_matrix.size(); ++i) {
    MpMatrix B = promote(fixture.matrices_matrix[i]);
    std::auto_ptr<SmithForm<MpMatrix> > smithForm(factory(B, true, true, true, true));

    (*smithForm)();

    BOOST_CHECK_EQUAL(promote(fixture.matrices_smithMatrix[i]), B);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PARISmithFormMp, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType Matrix;
  Fixture fixture;
  SmithFormFactory factory(true);

  for (size_t i = 0; i < fixture.matrices_matrix.size(); ++i) {
    MpMatrix B = promote(fixture.matrices_matrix[i]);
    std::auto_ptr<SmithForm<MpMatrix> > smithForm(factory(B, true, true, true, true));

    (*smithForm)();

    BOOST_CHECK_EQUAL(promote(fixture.matrices_smithMatrix[i]), B);
  }
}

#endif // HAVE_MPCAPD_ALG


BOOST_AUTO_TEST_SUITE_END()

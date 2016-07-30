//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file atanTest.cpp
///
/// @author Tomasz Kapela @date 2010-02-27
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

//#define BOOST_TEST_MODULE atanTest
#include "compare.h"
BOOST_AUTO_TEST_SUITE(atanSuite)

std::vector<double> computeAtanDer(MapType::VectorType & u){
  double x = u[0].leftBound();
  double y = u[1].leftBound();

  // code generated by the following Mathematica code
  // W[n_,m_]:=D[ArcTan[x*y],{x,n},{y,m}]/(n!m!)//FullSimplify
  // Table[Table[W[m-n,n]//CForm,{n,0,m}],{m,0,5}]//Flatten

  capd::rounding::DoubleRounding::roundNearest();
  double r[] = {ArcTan(x*y),y/(1 + Power(x,2)*Power(y,2)),x/(1 + Power(x,2)*Power(y,2)),-((x*Power(y,3))/Power(1 + Power(x,2)*Power(y,2),2)),(1 - Power(x,2)*Power(y,2))/Power(1 + Power(x,2)*Power(y,2),2),-((Power(x,3)*y)/Power(1 + Power(x,2)*Power(y,2),2)),(Power(y,3)*(-1 + 3*Power(x,2)*Power(y,2)))/(3.*Power(1 + Power(x,2)*Power(y,2),3)),(x*Power(y,2)*(-3 + Power(x,2)*Power(y,2)))/Power(1 + Power(x,2)*Power(y,2),3),(Power(x,2)*y*(-3 + Power(x,2)*Power(y,2)))/Power(1 + Power(x,2)*Power(y,2),3),(Power(x,3)*(-1 + 3*Power(x,2)*Power(y,2)))/(3.*Power(1 + Power(x,2)*Power(y,2),3)),(x*Power(y,5) - Power(x,3)*Power(y,7))/Power(1 + Power(x,2)*Power(y,2),4),-((Power(y,2)*(1 - 6*Power(x,2)*Power(y,2) + Power(x,4)*Power(y,4)))/Power(1 + Power(x,2)*Power(y,2),4)),-((x*y*(3 - 8*Power(x,2)*Power(y,2) + Power(x,4)*Power(y,4)))/Power(1 + Power(x,2)*Power(y,2),4)),-((Power(x,2)*(1 - 6*Power(x,2)*Power(y,2) + Power(x,4)*Power(y,4)))/Power(1 + Power(x,2)*Power(y,2),4)),(Power(x,5)*y - Power(x,7)*Power(y,3))/Power(1 + Power(x,2)*Power(y,2),4),(Power(y,5) - 10*Power(x,2)*Power(y,7) + 5*Power(x,4)*Power(y,9))/(5.*Power(1 + Power(x,2)*Power(y,2),5)),(x*Power(y,4)*(5 - 10*Power(x,2)*Power(y,2) + Power(x,4)*Power(y,4)))/Power(1 + Power(x,2)*Power(y,2),5),(y*(-1 + 15*Power(x,2)*Power(y,2) - 15*Power(x,4)*Power(y,4) + Power(x,6)*Power(y,6)))/Power(1 + Power(x,2)*Power(y,2),5),(x*(-1 + 15*Power(x,2)*Power(y,2) - 15*Power(x,4)*Power(y,4) + Power(x,6)*Power(y,6)))/Power(1 + Power(x,2)*Power(y,2),5),(Power(x,4)*y*(5 - 10*Power(x,2)*Power(y,2) + Power(x,4)*Power(y,4)))/Power(1 + Power(x,2)*Power(y,2),5),(Power(x,5) - 10*Power(x,7)*Power(y,2) + 5*Power(x,9)*Power(y,4))/(5.*Power(1 + Power(x,2)*Power(y,2),5))};
  return std::vector<double> (r,r+sizeof(r)/sizeof(double));
}


BOOST_AUTO_TEST_CASE(xatan)
{
  std::string txt = "var:x,y;fun:atan(x*y);",
              msg = "Function \"" + txt + "\"  x = " ;
  MapType f(txt,5);
  VectorType x(2);
  JetType df(1,2, 5);

  x[0] = 1.0; x[1]=1;
  std::vector<double> expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,1.0)");

  MapType g("var:x,y;fun:-atan(-x*y);",5);
  g(x,df);
  compareResults(expected, df, msg+"(1.0,1.0)");

  x[0] = 3.0; x[1]=-5.0;
  expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(3.0,-5.0)");

  x[0] = 0.0; x[1] = 0.0;
  expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");
}


using capd::autodiff::Node;

void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params*/[], int /*noParams*/)
{
  out[0] = 2*atan(in[0]*in[1]/(1+sqrt(1+sqr(in[0]*in[1]))));
}

BOOST_AUTO_TEST_CASE(xatannode)
{
  std::string msg = "Function \"2*atan(xy/(1+sqrt(1+sqr(xy))))\"  u = " ;
  MapType f(_f,2,1,0,5);
  VectorType x(2);
  JetType df(1,2,5);

  x[0] = 1.0; x[1] = -1;
  std::vector<double> expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(1.0,-1.0)");

  x[0] = 0.5; x[1] = 3;
  expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(.5,3.0)");

  x[0] = 0.0; x[1] = 0.0;
  expected = computeAtanDer(x);
  f(x,df);
  compareResults(expected, df, msg+"(0.0,0.0)");
}

BOOST_AUTO_TEST_SUITE_END()
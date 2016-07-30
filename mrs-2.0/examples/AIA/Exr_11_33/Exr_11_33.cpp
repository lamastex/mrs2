/* 
 * Copyright (C) 2008, 2009 Raazesh Sainudiin and Jennifer Harlow
 *
 * This file is part of mrs, a C++ class library for statistical set processing.
 *
 * mrs is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*! \file Exr_11_33.cpp
\brief Exercise 11.33 from Jaulin et al, Springer, 2001, p. 342.

Implementation of exercise 11.33 from Jaulin, Kieffer, Didrit and Walter, 
Applied Interval Analysis, Springer, p. 342.

*/
#include "Exr_11_33.hpp"

// These AIASubPavings are declared as global
AIASubPaving Sc, Sc1, Sc2;

// specifications of example interval boolean tests
// The boolean interval test can return BI_TRUE, BI_FALSE, or BI_INDET

AIA_BOOL_INTERVAL IBTAnnular(const ivector& x)
{
  // here we test a 2-d box for inclusion in the area between circles centred 
  // on the origin with radii 1 and 2
  interval ToInvert(1.0,2.0),Temp;
  interval Img = sqr(x[1]) + sqr(x[2]);

  if (!Intersection(Temp,Img,ToInvert)) return BI_FALSE;
  if ( Img<=ToInvert ) return BI_TRUE;

  return BI_INDET;
}

AIA_BOOL_INTERVAL IBTFdirect(const ivector& x)
{
  // A boolean interval test to illustrate SIVIA being used to evaluate the 
  // direct image of a set by a function provided the function is invertible
  // ie SIVIA will invert inverse_f

  interval Temp;
  ivector Img(2);

  // taking the function f : (x1,x2) -> (2x1-x2,-x1+2x2)
  // f is invertible and the inverse is 
  // inverse_f : (x1,x2) -> (2x1 + x2,x1+2x2)/3
  // this is an inclusion function for inverse_f

  Img[1] = (2.0*x[1]+x[2]) / 3.0;
  Img[2] = (x[1]+2.0*x[2]) / 3.0;

  return (Img<=Sc); // Sc is an AIASubPaving set up by the first example
}

AIA_BOOL_INTERVAL IBTFinverse(const ivector& x)
{
  interval Temp;
  ivector Img(2);

  // the function f is f : (x1,x2) -> (2x1-x2,-x1+2x2)

  Img[1] = 2.0*x[1]-x[2];
  Img[2] = -x[1]+2.0*x[2];

  return (Img<=Sc1);
}

// end specification of example boolean interval tests

using namespace cxsc;
using namespace std;

int main()
{
  double prec;
  clock_t start, end;

  ivector x(2);
  x[1] = interval(-5.0,5.0);
  x[2] = interval(-5.0,5.0);

  AIASubPaving A;
  A = new AIASPnode(x);

  // Using SIVIA for set inversion
  //find an AIASubPaving characterisation Sc containing the area between 
  // circles centred on the origin
  // with radii 1 and 2 (in 2 dimensional space)
  cout << "Characterization of the set Sc={(x1,x2) | 1 <= x1^2+x2^2 <= 2 }" 
       << endl;
  cout << "Enter a precision (between 1 and 0.001): ";
  cin >> prec;

  start = clock();

  // when we start we give A a box big enough to guarantee to contain 
  // the characterisation of Sc

  Sc = Sivia(IBTAnnular,A,prec);
  end = clock();

  cout << "Computing time : " 
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "Volume: " << Volume(Sc) << endl;
  cout << "Number of leaves: " << NbLeaves(Sc) << endl;

  // To realize a file output of the AIASubPaving Sc
                    // Filename
  ofstream os("AIAannular.txt");
  os << 2 << endl;  // Dimension of the AIASubPaving
                    // Root box
  os << interval(-5.0,5.0) << " "
    << interval(-5.0,5.0) << " " << endl;
                    // Precision used
  os << "Precision is " << prec << endl;
  os << Sc << endl; // AIASubPaving itself
  cout << "The output AIASubPaving has been written to AIAannular.txt" 
       << endl << endl;

  // end of testing to find reciprocal image

  // the AIASubPaving Sc that we have created is the regular AIASubPaving 
  // that covers the set
  // X = {(x1,x2) in R2 | sqr(x1) + sqr(x2) is in [1,2]} 
  // (remember that it contains this area rather than being this area, 
  // and that eps has determined how small we go in the AIASubPaving 
  // characterisation)

  // make a new AIASubPaving to provide an initial source box for the test next
  delete A;

  x[1] = interval(-5.0,5.0);
  x[2] = interval(-5.0,5.0);

  A = new AIASPnode(x);

  // testing using SIVIA to find the direct image of an invertible function
  // remember that we are only finding some upper enclosure of the direct 
  // image really
  // Note that this example will use the AIASubPaving Sc we created 
  // above - see IBTFdirect

  // ie create an AIASubPaving Sc1 containing f(Sc), where Sc was found above

  cout << "Characterization of the set Sc1=f(Sc)" << endl
    << "with f1(x) = 2*x1-x2," << endl
    << "      f2(x) = -x1+2*x2," << endl;
  cout << "by realizing the inversion of f-1 by Sivia" << endl;
  cout << "Enter a precision (between 1 and 0.01): ";
  cin >> prec;

  start = clock();
                    // Sc1 will be used by the following example
  Sc1 = Sivia(IBTFdirect,A,prec);
  end = clock();

  cout << "Computing time : " 
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "Volume: " << Volume(Sc1) << endl;
  cout << "Number of leaves: " << NbLeaves(Sc1) << endl;

  // To realize a file output of the AIASubPaving Sc1
                    // Filename
  ofstream os1("AIAdirect.txt");
  os1 << 2 << endl; // Dimension of the AIASubPaving
                    // Root box
  os1 << interval(-5.0,5.0) << " "
    << interval(-5.0,5.0) << " " << endl;
                    // Precision used
  os1 << "Precision is " << prec << endl;
                    // AIASubPaving itself
  os1 << Sc1 << endl;

  cout << "The output AIASubPaving has been written to AIAdirect.txt" 
       << endl << endl;

  // end of example for finding direct image

  // get a new AIASubPaving A to provide initial source box for next example
  delete A;

  x[1] = interval(-5.0,5.0);
  x[2] = interval(-5.0,5.0);

  A = new AIASPnode(x);

  // Image evaluation using set inversion
  // this uses the AIASubPaving Sc1 created by the above example
  // create an AIASubPaving Sc2 which contains inverse_f(Sc1)

  cout << "Characterization of the set Sc2=f-1(Sc1)" << endl
    << "with f^-1_1(x) = (2*x1+x2)/3," << endl
    << "     f^-1_2(x) = (x1+2*x2)/3," << endl;
  cout << "by realizing the inversion of f by Sivia" << endl;
  cout << "Enter a precision (between 1 and 0.01): ";
  cin >> prec;

  start = clock();
  Sc2 = Sivia(IBTFinverse,A,prec);
  end = clock();

  cout << "Computing time : " 
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "Volume: " << Volume(Sc2) << endl;
  cout << "Number of leaves: " << NbLeaves(Sc2) << endl;

  // To realize a file output of the AIASubPaving Sc
                    // Filename
  ofstream os2("AIAinverse.txt");
  os2 << 2 << endl; // Dimension of the AIASubPaving
                    // Root box
  os2 << interval(-5.0,5.0) << " "
    << interval(-5.0,5.0) << " " << endl;
                    // Precision used
  os2 << "Precision is " << prec << endl;
                    // AIASubPaving itself
  os2 << Sc2 << endl;

  cout << "The output AIASubPaving has been written to AIAinverse.txt" 
       << endl << endl;

  // end of testing SIVIA for set inversion
  // we should compare Sc2 to Sc in terms of volume and look at the effects of 
  // the precision variable

  delete A;         // delete subpavings newed in dyamic memory
  delete Sc;
  delete Sc1;
  delete Sc2;

  return 0;
}

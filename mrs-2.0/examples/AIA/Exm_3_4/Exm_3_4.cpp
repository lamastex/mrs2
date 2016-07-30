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

/*! \file Exm_3_4.cpp
\brief Example 3.4 from Jaulin et al, Springer, 2001, pp. 61-63.

Implementation of example 3.4 from Jaulin, Kieffer, Didrit and Walter, 
Applied Interval Analysis, Springer, 2001, pp. 61-63.

*/

#include "Exm_3_4.hpp"

// These AIASubPavings are declared as global
AIASubPaving Sc5, Sc6, Sc7;

// specifications of example interval boolean tests
// The boolean interval test can return BI_TRUE, BI_FALSE, or BI_INDET

AIA_BOOL_INTERVAL IBT_ex3_4(const ivector& x)
{
  // here we test a 2-d box for inclusion in the area 
  // such that x1^4 - x1^2 + 4x2^2 is in the interval [-0.1, 0.1]

  interval ToInvert(-0.1, 0.1),Temp;
  interval Img = power(x[1],4) - sqr(x[1]) + 4*sqr(x[2]);

  if (!Intersection(Temp,Img,ToInvert)) return BI_FALSE;
  if ( Img<=ToInvert ) return BI_TRUE;

  return BI_INDET;
}

AIA_BOOL_INTERVAL IBTFinverse_ex3_4(const ivector& x)
{
  ivector Img(2);
  // example in 2-d space R2
  // for f: R2 -> R2
  //f1(x1, x2) = (x1-1)^2 - 1+ x2
  //f2(x1, x2) = -(x1^2) + (x2-1)^2

  Img[1] = sqr(x[1]) - 2*x[1] + x[2];
  Img[2] = -sqr(x[1]) + sqr(x[2]) - 2*x[2] + 1;

  return (Img<=Sc6);
}

// end specification of example boolean interval tests

// specification of example interval vector function 
// (ie, vector inclusion function)

ivector IVF_ex3_4(const ivector& x)
{
  // example in 2-d space R2
  // for f: R2 -> R2
  //f1(x1, x2) = (x1-1)^2 - 1+ x2
  //f2(x1, x2) = -(x1^2) + (x2-1)^2

  ivector Img(2);

  Img[1] = sqr(x[1]) - 2*x[1] + x[2];
  Img[2] = -sqr(x[1]) + sqr(x[2]) - 2*x[2] + 1;

  return (Img);
}

// end specification of example interval vector functions

using namespace cxsc;
using namespace std;

//int main(int argc, char* argv[])
int main()
{
  double prec;
  clock_t start, end;

  ivector x(2);
  x[1] = interval(-3.0,3.0);
  x[2] = interval(-3.0,3.0);

  AIASubPaving A;
  A = new AIASPnode(x);

  // Using SIVIA for set inversion

  // Find an AIASubPaving characterisation Sc5 as in example 3.4

  cout << "Characterization of the set Sc5={(x1,x2) | -0.1 " 
       << "<= x1^4-x1^2+4x2^2 <= 0.1 }" << endl;
  cout << "Enter a precision (between 1 and 0.001): ";
  cin >> prec;

  start = clock();

  // when we start we give A a box big enough to guarantee to contain 
  // the characterisation of Sc

  Sc5 = Sivia(IBT_ex3_4,A,prec);
  end = clock();

  cout << "Computing time : " 
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "Volume: " << Volume(Sc5) << endl;
  cout << "Number of leaves: " << NbLeaves(Sc5) << endl;

  // To realize a file output of the AIASubPaving Sc
                    // Filename
  ofstream os5("AIA3_4a.txt");
  os5 << 2 << endl; // Dimension of the AIASubPaving
                    // Root box
  os5 << interval(-3.0,3.0) << " "
    << interval(-3.0,3.0) << " " << endl;
                    // Precision used
  os5 << "Precision is " << prec << endl;
                    // AIASubPaving itself
  os5 << Sc5 << endl;
  cout << "The output AIASubPaving has been written to AIA3_4a.txt" 
       << endl << endl;

  // Using ImageSp to find the image of Sc5 using
  // 	f1(x) = (x1-1)^2 - 1+ x2
  // "    f2(x) = -(x1^2) + (x2-1)^2

  cout << "Characterization of the set Sc6=f(Sc5)" << endl
    << " with Sc5 from our previous example and "<< endl
    << "with f1(x) = (x1-1)^2 - 1 +x2," << endl
    << "     f2(x) = -(x1^2) + (x2-1)^2" << endl;
  cout << "by realizing the image of f by ImageSp" << endl;
  cout << "Enter a precision (between 1 and 0.01): ";
  cin >> prec;

  start = clock();

  // use Image SP to find a characterisation of the 
  // image of Sc5 using the function in IVF
  Sc6 = ImageSp(IVF_ex3_4, Sc5, prec);

  end = clock();

  cout << "Computing time : " 
      << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "The volume is " << Volume(Sc6) << endl;
  cout << "The number of leaves is " << NbLeaves(Sc6) << endl;

  // To realize a file output of the AIASubPaving Sc
                    // Filename
  ofstream os6("AIA3_4b.txt");
  os6 << 2 << endl; // Dimension of the AIASubPaving
                    // Domain AIASubPaving
  os6 << interval(-3.0,3.0) << " "
    << interval(-3.0,3.0) << " " << endl;
                    // Precision used
  os6 << "Precision is " << prec << endl;
                    // Image AIASubPaving itself
  os6 << Sc6 << endl;

  cout << "The output AIASubPaving has been written to AIA3_4b.txt" 
       << endl << endl;

  //remake A
  delete A;

  x[1] = interval(-5.0,5.0);
  x[2] = interval(-5.0,5.0);

  A = new AIASPnode(x);

  // set inversion using SIVIA
  // this uses the AIASubPaving Sc6 created by the above example
  // create an AIASubPaving Sc7 which contains inverse_f(Sc6)

  cout << "Characterization of the set Sc7=f-1(Sc6)" << endl
    << "with f as above  and Sc6 as above " << endl;
  cout << "by realizing the inversion of f by Sivia" << endl;
  cout << "Enter a precision (between 1 and 0.01): ";
  cin >> prec;

  start = clock();
  Sc7 = Sivia(IBTFinverse_ex3_4,A,prec);
  end = clock();

  cout << "Computing time : " 
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
  cout << "Volume: " << Volume(Sc7) << endl;
  cout << "Number of leaves: " << NbLeaves(Sc7) << endl;

  // To realize a file output of the AIASubPaving Sc
                    // Filename
  ofstream os7("AIA3_4c.txt");
  os7 << 2 << endl; // Dimension of the AIASubPaving
                    // Root box
  os7 << interval(-5.0,5.0) << " "
    << interval(-5.0,5.0) << " " << endl;
                    // Precision used
  os7 << "Precision is " << prec << endl;
                    // AIASubPaving itself
  os7 << Sc7 << endl;

  cout << "The output AIASubPaving has been written to AIA3_4c.txt" 
       << endl << endl;

  delete A;         // delete all Subpavings newed in dynamic memory
  delete Sc5;
  delete Sc6;
  delete Sc7;

  return 0;
}

/*
 * Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
 * Copyright (C) 2009 Jennifer Harlow
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

/*! \file NewExr_11_33.cpp
\brief Exercise 11.33 from Jaulin et al, Springer, 2001 using class SPMinimalnode
*/

// Exercise 11.33 from Jaulin et al, Springer, 2001 using class SPMinimalnode
#include "NewExr_11_33.hpp"

// specifications of example interval boolean tests
// The boolean interval test can return BI_TRUE, BI_FALSE, or BI_INDET


BOOL_INTERVAL IBTAnnular(const ivector& x, const SPMinimalnode * const spn)
{
    ivector Img(1);
  // the inclusion function image of the given interval x
    Img[1] = sqr(x[1]) + sqr(x[2]);

    return Img<=spn;
}

// A boolean interval test to illustrate SIVIA being used to evaluate the
// direct image of a set by a function provided the function is invertible
// Tests if the image of the ivector x in question is in the supplied subpaving
// spn, ie, as used by SIVIA, SIVIA will invert inverse_f: find if the inverse
// image of the supplied ivector x (a proposed part of the image set) is in the
// domain suppaving spn

BOOL_INTERVAL IBTFdirect(const ivector& x, const SPMinimalnode * const spn)
{

    ivector Img(2);

    // the function f is f : (x1,x2) -> (2x1-x2,-x1+2x2)
    // f is invertible and the inverse is
  // inverse_f : (x1,x2) -> (2x1 + x2,x1+2x2)/3

    Img[1] = (2.0*x[1]+x[2]) / 3.0;
    Img[2] = (x[1]+2.0*x[2]) / 3.0;

    return (Img<=spn);
}

// A boolean interval test to illustrate SIVIA being used to
// evaluate the domain corresponding to an
// image subpaving spn and a function f
BOOL_INTERVAL IBTFinverse(const ivector& x, const SPMinimalnode * const spn)
{

    ivector Img(2);

    Img[1] = 2.0*x[1]-x[2];
    Img[2] = -x[1]+2.0*x[2];

    return (Img<=spn);
}

// end specification of example boolean interval tests


int main()
{

    double prec;
    clock_t start, end;

    ivector x(2);
    x[1] = interval(-5.0,5.0);
    x[2] = interval(-5.0,5.0);

    MinimalSubPaving A;
    A = new SPMinimalnode(x);

    ivector t(1);
    t[1] = interval(1.0,2.0);   // a restriction on the range of the
                            // function - we want to find the x
                    // for which f(x) is in this interval
                    // ie acts like a simple MinimalSubPaving for SIVIA

    MinimalSubPaving ToInvert = new SPMinimalnode(t);

    // Using SIVIA for set inversion
    //find an MinimalSubPaving characterisation Sc containing the
  // area between circles centred on the origin
        // with radii 1 and 2 (in 2 dimensional space)
    cout << "Characterization of the set Sc={(x1,x2) | 1 <= x1^2+x2^2 <= 2 }"
       << endl;
    cout << "Enter a precision (between 1 and 0.001): ";
    cin >> prec;

    start = clock();

    // when we start we give A a box big enough to guarantee to contain the
  // characterisation of Sc

    MinimalSubPaving Sc = Sivia(IBTAnnular, ToInvert, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc) << endl;
    cout << "The number of leaves is " << spLeaves(Sc) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os("new_annular.txt");         // Filename
        os << 2 << endl;                    // Dimension of the MinimalSubPaving
        os << interval(-5.0,5.0) << " "         // Root box
            << interval(-5.0,5.0) << " " << endl;
        os << "Precision is " << prec << endl; // Precision used
    os << Sc << endl;                   // MinimalSubPaving itself
    cout << "The output MinimalSubPaving has been written to new_annular.txt"
      << endl << endl;

    // end of testing to find reciprocal image

    // the MinimalSubPaving Sc that we have created is the regular MinimalSubPaving that
  // covers the set
    // X = {(x1,x2) in R2 | sqr(x1) + sqr(x2) is in [1,2]}
  // (remember that it contains this area rather than being this area, and
  // that eps has determined how small we go in the MinimalSubPaving characterisation)


    // make a new MinimalSubPaving to provide an initial source box for the test next
    delete A;

    x[1] = interval(-5.0,5.0);
    x[2] = interval(-5.0,5.0);

    A = new SPMinimalnode(x);

    // testing using SIVIA to find the direct image of an invertible function
    // remember that we are only finding some upper enclosure of the direct
  // image really
    // Note that this example will use the MinimalSubPaving Sc we created
  // above - see IBTFdirect

    // ie create an MinimalSubPaving Sc1 containing f(Sc), where Sc was found above

    cout << "Characterization of the set Sc1=f(Sc)" << endl
    << "with f1(x) = 2*x1-x2," << endl
        << "      f2(x) = -x1+2*x2," << endl;
    cout << "by realizing the inversion of f-1 by Sivia" << endl;
    cout << "Enter a precision (between 1 and 0.01): ";
    cin >> prec;

    start = clock();
  // Sc1 will be used by the following example
    MinimalSubPaving Sc1 = Sivia(IBTFdirect, Sc, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc1) << endl;
    cout << "The number of leaves is " << spLeaves(Sc1) << endl;

    // To realize a file output of the MinimalSubPaving Sc1
    ofstream os1("new_direct.txt");         // Filename
    os1 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os1 << interval(-5.0,5.0) << " "         // Root box
    << interval(-5.0,5.0) << " " << endl;
        os1 << "Precision is " << prec << endl; // Precision used
    os1 << Sc1 << endl;                   // MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to new_direct.txt"
       << endl << endl;

    // end of example for finding direct image


    // get a new MinimalSubPaving A to provide initial source box for next example
    delete A;

    x[1] = interval(-5.0,5.0);
    x[2] = interval(-5.0,5.0);

    A = new SPMinimalnode(x);

    // IMAGE EVALUATION USING SET-INVERSION - PART II
    // this uses the MinimalSubPaving Sc1 created by the above example
    // create an MinimalSubPaving Sc2 which contains inverse_f(Sc1)

    cout << "Characterization of the set Sc2=f-1(Sc1)" << endl
    << "with f^-1_1(x) = (2*x1+x2)/3," << endl
        << "     f^-1_2(x) = (x1+2*x2)/3," << endl;
    cout << "by realizing the inversion of f by Sivia" << endl;
    cout << "Enter a precision (between 1 and 0.01): ";
    cin >> prec;

    start = clock();
    MinimalSubPaving Sc2 = Sivia(IBTFinverse, Sc1, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc2) << endl;
    cout << "The number of leaves is " << spLeaves(Sc2) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os2("new_inverse.txt");            // Filename
    os2 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os2 << interval(-5.0,5.0) << " "         // Root box
    << interval(-5.0,5.0) << " " << endl;
        os2 << "Precision is " << prec << endl; // Precision used
        os2 << Sc2 << endl;                   // MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to new_inverse.txt"
       << endl << endl;


    // end of testing SIVIA for set inversion
    // we should compare Sc2 to Sc in terms of volume and look at the
  // effects of the precision variable

    delete A;
    delete Sc;
    delete Sc1;
    delete Sc2;

    return 0;
}

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

/*! \file NewExm_3_3.cpp
\brief Example 3.3 from Jaulin et al, Springer, 2001 using class SPMinimalnode
*/

#include "NewExm_3_3.hpp"

// specifications of example interval boolean tests
// The boolean interval test can return BI_TRUE, BI_FALSE, or BI_INDET

BOOL_INTERVAL IBTAnnular(const ivector& x, const SPMinimalnode * const spn)
{
    ivector Img(1);
  // the inclusion function image of the given interval x
    Img[1] = sqr(x[1]) + sqr(x[2]);

    return Img<=spn;
}


// end specification of example boolean interval tests

// specification of example interval vector function
// (ie, vector inclusion function)

// specification of example interval vector function
// (ie, vector inclusion function)
// as in Example 3.3 of AIA
ivector IVF_ex3_3(const ivector& x)
{
    // example in 2-d space R2
    // for f: R2 -> R2
    //f1(x1, x2) = x1*x2
    //f2(x1, x2) = x1 + x2

    ivector Img(2);

    Img[1] = x[1] * x[2];
    Img[2] = x[1] + x[2];

    return (Img);
}

// end specification of example interval vector funciions



int main()
{

    double prec;
    clock_t start, end;

    ivector x(2);
    x[1] = interval(-3.0,3.0);
    x[2] = interval(-3.0,3.0);

    MinimalSubPaving A;
    A = new SPMinimalnode(x);

    ivector t(1);
  // a restriction on the range of the function - we want to find the x
    t[1] = interval(1.0,2.0);
                    // for which f(x) is in this interval
                    // ie acts like a simple MinimalSubPaving for SIVIA

    MinimalSubPaving ToInvert = new SPMinimalnode(t);

    // Using SIVIA for set inversion
    //find an MinimalSubPaving characterisation Sc containing
  // the area between circles centred on the origin
        // with radii 1 and 2 (in 2 dimensional space)
    cout << "Characterization of the set Sc={(x1,x2) | 1 <= x1^2+x2^2 <= 2 }"
       << endl;
    cout << "Enter a precision (between 1 and 0.001): ";
    cin >> prec;

    start = clock();

    // when we start we give A a box big enough to guarantee to
  // contain the characterisation of Sc

    MinimalSubPaving Sc = Sivia(IBTAnnular, ToInvert, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc) << endl;
    cout << "The number of leaves is " << spLeaves(Sc) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os("New3_3a.txt");         // Filename
        os << 2 << endl;                    // Dimension of the MinimalSubPaving
        os << interval(-5.0,5.0) << " "         // Root box
            << interval(-5.0,5.0) << " " << endl;
        os << "Precision is " << prec << endl; // Precision used
    os << Sc << endl;                   // MinimalSubPaving itself
    cout << "The output MinimalSubPaving has been written to New3_3a.txt"
       << endl << endl;

    // the MinimalSubPaving Sc that we have created is the regular MinimalSubPaving
  // that covers the set
    // X = {(x1,x2) in R2 | sqr(x1) + sqr(x2) is in [1,2]}
  // (remember that it contains this area rather than being this area,
  // and that eps has determined how small we go in the MinimalSubPaving
  // characterisation)

    cout << "Characterization of the set Sc4=f(Sc)" << endl
        << " with Sc from our first example and "<< endl
        << "with f1(x) = x1*x2," << endl
        << "     f2(x) = x1 + x2," << endl;
    cout << "by realizing the image of f by ImageSp" << endl;
    cout << "Enter a precision (between 1 and 0.01): ";
    cin >> prec;

    start = clock();

    //use Image SP to find a characterisation of the image of Sc using the
  // function in IVF
    MinimalSubPaving Sc4 = ImageSp(IVF_ex3_3, Sc, prec);

    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

    cout << "The volume is " << spVolume(Sc4) << endl;
    cout << "The number of leaves is " << spLeaves(Sc4) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os4("New3_3d.txt");            // Filename
    os4 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os4 << interval(-3.0,3.0) << " "         // Domain MinimalSubPaving
        << interval(-3.0,3.0) << " " << endl;
    os4 << "Precision is " << prec << endl; // Precision used
    os4 << Sc4 << endl;                   // Image MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to New3_3d.txt"
       << endl << endl;

    delete A;
    delete Sc;
    delete Sc4;

    return 0;
}

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

/*! \file NewExr_11_35.cpp
\brief Exercise 11.35 from Jaulin et al, Springer, 2001 using class SPMinimalnode
*/

// implementation of example 11.35 from Jaulin et al, 2001

#include "NewExr_11_35.hpp"

// specification of example interval vector function
// (ie, vector inclusion function)
// as in Exercise 11.34 and 11.35 of AIA
ivector IVF_ex11_35(const ivector& x)
{
    // example in 2-d space R2
    // for f: R2 -> R2
    //f1(x1, x2) = x1^2 + x2 - x2^2
    //f2(x1, x2) = x1^2 + x2^2

    ivector Img(2);

    Img[1] = sqr(x[1]) + x[2] - sqr(x[2]);
    Img[2] = sqr(x[1]) + sqr(x[2]);

    return (Img);
}

// end specification of example interval vector funciions


int main()
{
    double prec;
    clock_t start, end;

    ivector x(2);
    x[1] = interval(-2.0,2.0);
    x[2] = interval(-2.0,2.0);

    MinimalSubPaving A;
    A = new SPMinimalnode(x);

    // find the image of A with the function
    // f1 = x1^2 + x2 - x2^2
    // f1 = x1^2 + x2^2

    cout << "Characterization of the set Sc3=f(A)" << endl
        << "with f1(x) = x1^2 + x2 - x2^2," << endl
        << "     f2(x) = x1^2 + x2^2," << endl;
    cout << "by realizing the image of f by ImageSp" << endl;
    cout << "Enter a precision (between 1 and 0.001): ";
    cin >> prec;

    start = clock();

    // use Image SP to find a characterisation of the image
  // of Sc using the function in IVF
    MinimalSubPaving Sc3 = ImageSp(IVF_ex11_35, A, prec);

    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

    cout << "The volume is " << spVolume(Sc3) << endl;
    cout << "The number of leaves is " << spLeaves(Sc3) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os4("new_image.txt");            // Filename
    os4 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os4 << interval(-2.0,2.0) << " "         // Domain MinimalSubPaving
        << interval(-2.0,2.0) << " " << endl;
    os4 << "Precision is " << prec << endl; // Precision used
    os4 << Sc3 << endl;                   // Image MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to new_image.txt"
       << endl << endl;

    delete A;
    delete Sc3;

    return 0;
}

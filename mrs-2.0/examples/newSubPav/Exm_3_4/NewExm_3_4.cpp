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


/*! \file NewExm_3_4.cpp
\brief Example 3.4 from Jaulin et al, Springer, 2001 using class SPMinimalnode
*/

// implementation of AIA Example 3.4 using class SPMinimalnode
#include "NewExm_3_4.hpp"

// specifications of example interval boolean tests
// The boolean interval test can return BI_TRUE, BI_FALSE, or BI_INDET
BOOL_INTERVAL IBT_ex3_4(const ivector& x, const SPMinimalnode * const spn)
{
    ivector Img(1);
  // the inclusion function image of the given interval x
    Img[1] = sqr(sqr(x[1])) - sqr(x[1]) + 4*sqr(x[2]);

    return Img<=spn;

}


// A boolean interval test to illustrate SIVIA being used
// to evaluate the domain corresponding to an
// image subpaving spn and a function f
BOOL_INTERVAL IBTFinverse_ex3_4(const ivector& x, const SPMinimalnode * const spn)
{

    //f1(x1, x2) = (x1-1)^2 - 1+ x2
    //f2(x1, x2) = -(x1^2) + (x2-1)^2

    ivector Img(2);

    Img[1] = sqr(x[1]) - 2*x[1] + x[2];
    Img[2] = -sqr(x[1]) + sqr(x[2]) - 2*x[2] + 1;

    return (Img<=spn);
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
    t[1] = interval(-0.1,0.1);
                    // for which f(x) is in this interval
                    // ie acts like a simple MinimalSubPaving for SIVIA

    MinimalSubPaving ToInvert = new SPMinimalnode(t);

    // Using SIVIA for set inversion

    //find an MinimalSubPaving characterisation Sc5 as in example 3.4

    cout << "Characterization of the set Sc5={(x1,x2) | -0.1 "
       << "<= x1^4-x1^2+4x2^2 <= 0.1 }" << endl;
    cout << "Enter a precision (between 1 and 0.001): ";
    cin >> prec;

    start = clock();

    // when we start we give A a box big enough to guarantee to contain
    // the characterisation of Sc

    MinimalSubPaving Sc5 = Sivia(IBT_ex3_4, ToInvert, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc5) << endl;
    cout << "The number of leaves is " << spLeaves(Sc5) << endl;

    // To realize a file output of the MinimalSubPaving Sc5
    ofstream os5("New3_4a.txt");         // Filename
    os5 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os5 << interval(-3.0,3.0) << " "         // Root box
        << interval(-3.0,3.0) << " " << endl;
    os5 << "Precision is " << prec << endl; // Precision used
    os5 << Sc5 << endl;                   // MinimalSubPaving itself
    cout << "The output MinimalSubPaving has been written to New3_4a.txt"
       << endl << endl;


    // Using ImageSp to find the image of Sc5 using
    //  f1(x) = (x1-1)^2 - 1+ x2
    // "    f2(x) = -(x1^2) + (x2-1)^2


    cout << "Characterization of the set Sc6=f(Sc5)" << endl
        << " with Sc5 from our previous example and "<< endl
        << "with f1(x) = (x1-1)^2 - 1 +x2," << endl
        << "     f2(x) = -(x1^2) + (x2-1)^2" << endl;
    cout << "by realizing the image of f by ImageSp" << endl;
    cout << "Enter a precision (between 1 and 0.01): ";
    cin >> prec;

    start = clock();

    // use Image SP to find a characterisation of the image of Sc5
  // using the function in IVF
    MinimalSubPaving Sc6 = ImageSp(IVF_ex3_4, Sc5, prec);

    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;

    cout << "The volume is " << spVolume(Sc6) << endl;
    cout << "The number of leaves is " << spLeaves(Sc6) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os6("New3_4b.txt");            // Filename
    os6 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os6 << interval(-3.0,3.0) << " "         // Domain MinimalSubPaving
        << interval(-3.0,3.0) << " " << endl;
    os6 << "Precision is " << prec << endl; // Precision used
    os6 << Sc6 << endl;                   // Image MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to New3_4b.txt"
       << endl << endl;

    //remake A
    delete A;

    x[1] = interval(-5.0,5.0);
    x[2] = interval(-5.0,5.0);

    A = new SPMinimalnode(x);


    // this uses the MinimalSubPaving Sc6 created by the above example
    // create an MinimalSubPaving Sc7 which contains inverse_f(Sc6)

    cout << "Characterization of the set Sc7=f-1(Sc6)" << endl
        << "with f as above  and Sc6 as above " << endl;
    cout << "by realizing the inversion of f by Sivia" << endl;
    cout << "Enter a precision (between 1 and 0.01): ";
    cin >> prec;

    start = clock();
    MinimalSubPaving Sc7 = Sivia(IBTFinverse_ex3_4, Sc6, A, prec);
    end = clock();

    cout << "Computing time : "
       << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
    cout << "Volume: " << spVolume(Sc7) << endl;
    cout << "The number of leaves is " << spLeaves(Sc7) << endl;

    // To realize a file output of the MinimalSubPaving Sc
    ofstream os7("New3_4c.txt");            // Filename
    os7 << 2 << endl;                    // Dimension of the MinimalSubPaving
    os7 << interval(-5.0,5.0) << " "         // Root box
        << interval(-5.0,5.0) << " " << endl;
    os7 << "Precision is " << prec << endl; // Precision used
    os7 << Sc7 << endl;                   // MinimalSubPaving itself

    cout << "The output MinimalSubPaving has been written to New3_4c.txt"
       << endl << endl;

    delete A;
    delete Sc5;
    delete Sc6;
    delete Sc7;

    return 0;
}

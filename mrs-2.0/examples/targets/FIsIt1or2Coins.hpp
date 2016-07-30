/* 
 * Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
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
/*! \file FIsIt1or2Coins.hpp
\brief Trans-dimensional posterior distribution over prob. of heads for possibly two coins
model declarations.
*/

#ifndef __FISIT1OR2COINS__
#define __FISIT1OR2COINS__

/*! \brief trans-dimensional posterior distribution given coin flipping data for possibly two coins (N1, n1), (N2, n2)
  as a function object class.

   Mixture of 1 parameter model (prob. of heads q for both coins
   and 2 parameter model (probs. of heads, q1, q2)
   priors are uniform for both. 
*/
class FIsIt1or2Coins: public Fobj{
    int n1; //!< n1 = number of Heads of first coin 
    int n2; //!< n2 = number of Heads of second coin 
    int N1; //!< N1 = number of tosses of first coin 
    int N2; //!< N2 = number of tosses of second coin 
    interval bc1, bc2;
    real rbc1, rbc2;
    //! this is needed sometimes to ensure local lipschitz condition is met
    real SmallestPositive;
    //! this is needed when the binamial coefficients leave the screen and we want posterior probs
    bool IgnoreBinomCoeffs;
    // vector<LabBox> LabDomainList;
    interval binary_coefficient(int N, int n);
    //! Prior prob. of model 1
    real w1; 
    //! Prior prob. of model 2
    real w2; 
    //! Volume of each labelled domain box
    vector<real> vol_labelled_domain;
    //! Track number of interval function calls
    mutable int n_interval_calls;
    //! Track number of real function calls
    mutable int n_real_calls;
public:
    //! A constructor.
    FIsIt1or2Coins(int N1, int n1, int N2, int n2, bool LogPi);
    // vector<LabBox> get_domain();
    //! interval function object
    interval operator()(const LabBox& X) const;
    //! real function object
    real operator()(const LabPnt& X) const ;   
    //! get volume of a labeled box
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    //! Get number of interval function calls
    int get_interval_calls(){ return n_interval_calls; }
    //! Get number of real function calls
    int get_real_calls(){ return n_real_calls; }
};

#endif


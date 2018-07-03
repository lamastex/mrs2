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
/*! \file FBinomialPartition.hpp
\brief Trans-dimensional posterior distribution over m binomial partition 
model declarations.
*/

#ifndef __FBINOMIALPARTITION__
#define __FBINOMIALPARTITION__

/*! \brief Trans-dimensional posterior distribution over m binomial partition
  model as a function object class.

   priors are uniform for both. 
*/
class FBinomialPartition: public Fobj{
    std::vector<int> y; //!< y = [y_1,y_2,...,y_k] number of 1s
    std::vector<int> n; //!< n = [n_1,n_2,...,n_k] number of trials
    //! array of values for the Stirling numbers of the second kind (sequence A008277 in OEIS)
    vector< vector<int> > MbyBMatrix;//this is not really used but there for giggles
    //! blocks or nonempty subsets of {1,2,3,4,5}
    vector< vector<int> > Blocks;
    //! partititon of experiment made up of blocks
    vector< vector<int> > Partns;
    vector< vector<int> > Succss;
    vector< vector<int> > Trials;
    vector< vector<int> > Failrs;
    vector< vector<interval> > Ibc; //! model-specific interval binomial coefficients
    vector< vector<real> > Rbc; //! model-specific real binomial coefficients
    // vector<LabBox> LabDomainList;
    interval binary_coefficient(int N, int Y);
    //! Prior prob. of model i
    vector<real> modelprior; 
    //! Volume of each labelled domain box
    vector<real> vol_labelled_domain;
    //! this is needed sometimes to ensure local lipschitz condition is met
    real SmallestPositive;
    //! this is needed when the binamial coefficients leave the screen and we want posterior probs
    bool IgnoreBinomCoeffs;
    //! Track number of interval function calls
    mutable int n_interval_calls;
    //! Track number of real function calls
    mutable int n_real_calls;
public:
    //! A default hard-coded constructor.
    FBinomialPartition(bool LogPi);
    //! An initialized constructor.
    FBinomialPartition(vector<int> ns, vector<int> ys, bool LogPi);
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


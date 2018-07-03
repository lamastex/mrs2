/* 
 * Copyright (C) 2005, 2006, 2007, 2008 Raazesh Sainudiin and Thomas York
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

#ifndef __SMALLCLASSES_HPP__
#define __SMALLCLASSES_HPP__

#include "cxsc.hpp"

#include <algorithm>
#include <set>
#include <vector>
#include <ostream>

//! A labeled point class
class LabPnt
{
  public:
    //! specifies the point as cxsc::rvector Pnt of the labeled point LabPnt
    rvector Pnt;
    //! specifies the label L of the labeled point LabPnt
    int L;
    //! print Pnt of the labeled point LabPnt to an output file stream out
    void Print(std::ostream& out)
    {
      out << L; for(int i=Lb(Pnt); i<=Ub(Pnt); i++)
      {
        out << " " << Pnt[i];
      }
      out << std::endl; return;
    }
};

/*! \brief A weighted labeled point class for a labeled point-valued particle.

  Used in trans-dimensional importance sampling and sequential importance 
sampling.  Ideal for labeled, point-valued, trans-dimensional particle systems.
*/
class WLabPnt: public LabPnt
{
  public:
    double qPnt;    //!< Proposal density of labeled point
    double fPnt;    //!< Target shape at labeled point -- f^*(Pnt)
    double Wt;      //!< Weight of the labeled point -- fPnt/qPnt
    
    //! Print the weight, label and the point of the particle
    void Print(std::ostream& out)
    {
      out << Wt << " " << L; for(int i=Lb(Pnt); i<=Ub(Pnt); i++)
      {
        out << " " << Pnt[i];
      }
      out << std::endl; return;
    }
};

//! A labeled box class
class LabBox
{
  public:
    //! specifies the box as cxsc::ivector Box of the labeled box LabBox
    ivector Box;
    //! specifies the label L of the labeled box LabBox
    int L;
    //! number of samples to be drawn from each labeled box LabBox
    int SamplesToDo;

    //! Output formatted with brackets etc, suitable for human reading
    void Print(std::ostream& out)
    {
      out << L;
      for(int i=Lb(Box); i<=Ub(Box); i++)
      {
        out << " " << Box[i];
      }
      out << std::endl;
      return;
    }

    /*! \brief Tab-delimited output format for a labelled box, numeric only. 
      
      Suitable to be read with MATLAB's dlmread etc.
    */
    void Output(std::ostream& out)
    {
      out << L;
      for(int i=Lb(Box); i<=Ub(Box); i++)
      {
        out << "\t" << Inf(Box[i]) << "\t" << Sup(Box[i]);
      }
      out << std::endl;
      return;
    }

};

//! A parameters class for an isotropic normal distribution
class NormalParam
{
  public:
    NormalParam(){}
    NormalParam(const NormalParam& np){ Mean = np.Mean; Sigma = np.Sigma; }
    NormalParam(rvector M, real S):Mean(M), Sigma(S){}
    rvector Mean;
    real Sigma;
};

//! A class for `labeled box with range' of a real-valued function over it.
class RangedLabBox
{
  public:
    //! The labeled box
    LabBox LBox;
    //! The interval range enclosure of the function over the labeled box
    interval BoxRE;
    //! Product of diams of intervals in labeled box, ie volume of labeled box.
    real BoxVol;
    //! Integral of the prior density over the box
    real BoxPriorIntegral;
    
    /*! \brief interval enclosure of the integral of the function over the 
      labeled box with scaling.
      
      We use a f_scale value to stay in the number screen.  Necessary for 
      product likelihod targets. 
    */
    interval BoxIntegral(bool logpi, real f_scale) const
    {
      return (logpi)? 
        exp(BoxRE-f_scale)*BoxPriorIntegral: BoxRE/f_scale*BoxPriorIntegral;
    }

    //! interval enclosure of the integral of the function over the labeled box.
    interval BoxIntegral(bool logpi) const
    {
      return (logpi)? exp(BoxRE)*BoxPriorIntegral: BoxRE*BoxPriorIntegral;
    }
    
    //! print the label, components of the box, range enclosure and integral. 
    void Print(std::ostream& out) const
    {
      out << LBox.L;
      for(int i=Lb(LBox.Box); i<=Ub(LBox.Box); i++) out << " " << LBox.Box[i];
      out << " RE: " << BoxRE <<  " BoxIntegral: " << BoxIntegral(0) << std::endl;
      return;
    }
};

/*! \brief A Function Object class for a sorting criterion between one 
  RangedLabBox and another.

  Determines a total order between two RangedLabBox es on the basis of the 
  magnitude of the diameter of the interval enclosure of the integral of the 
  function over the labeled box of each RangedLabBox.  Thus, SortBox is a type 
  that we use as a template argument for an associative container 
  RangedLabBoxSet, a type of STL-bsed set.
  Read for e.g. Josuttis1999 p. 178, p. 294 for details.

    \todo
    Read Josuttis1999 Sec 8.1.2 to extend and intervene the automatic sorting 
    criterion via Function Objects with more than one internal state at the 
    same time [perhaps also try defining different sorting criterion of the 
    same data type as in p. 178 'As a constructor parameter' with e.g. in P. 
    191].  Perhaps incorporate alternative sorting criterion into SortBox.  
    Older ideas for total order between RangedLabBoxes
    include:
    class SortByBoxREDiam { // sort by diameter of range enclosure
    public:
      bool operator() ( const RangedLabBox& P1, const RangedLabBox& P2 ) const {
        return ( diam(P1.BoxRE) <= diam(P2.BoxRE) ); }
    };
    and
    class SortByBoxVol { // sort by volume of subdomain box
    public:
      bool operator() ( const RangedLabBox& P1, const RangedLabBox& P2 ) const {
           return ( (P1.BoxVol) <= (P2.BoxVol) ); }
    };
*/
class SortBox
{
  /*! \brief Boolean flag for the natural logarithmic scale for the function. 
  
    Works when target function is a positive density (up to a normalizing 
    constant).
  */
  bool UsingLogDensity;
  public:
    //! Comparison operator
    SortBox(bool LogPi = false):UsingLogDensity(LogPi){}
    bool operator() ( const RangedLabBox& P1, const RangedLabBox& P2 ) const
    {
      //  cout << diam(P1.BoxRE) << " " << P1.BoxVol << " " 
      //       << diam(P2.BoxRE) << " " << P2.BoxVol << std::endl;
      return (UsingLogDensity)?
      (ln(P1.BoxPriorIntegral)+Sup(P1.BoxRE)+ln(1.0 - exp(-diam(P1.BoxRE))) >=
       ln(P2.BoxPriorIntegral)+Sup(P2.BoxRE)+ln(1.0 - exp(-diam(P2.BoxRE))) ):
      (diam(P1.BoxRE)*P1.BoxPriorIntegral>=diam(P2.BoxRE)*P2.BoxPriorIntegral);
    }
};

/*! \brief Sorted associative STL set container for RangedLabBoxes

  STL-based set that uses the SortBox sorting criterion.
*/
typedef  std::set<RangedLabBox, SortBox> RangedLabBoxSet;

struct MeansLabelsProportions{
  std::vector<rvector> means;
  std::vector<int> labels;
  std::vector<real> proportions;
};

//! A class for the status of a Rejection Sampler
class RSSample
{
  public:
  
    //! Number of draws from proposal distribution.
    long Nprop;
    
    //! Number of labels or topologies.
    long Ntopologies;
    
    //! The set of unique integer labels in LabDomainList.
    std::vector<int> LabelSet;
    
    //! The maximum dimension of the labeled boxes in LabDomainList
    int n_Dim_Max; 
    
    //! The envelope integral as a cxsc::real
    real EnvelopeIntegral;
    
    //! An STL vector container to store accepted samples of labeled points
    std::vector<LabPnt> Samples;
    
    //! A real estimate of the integral of the function over the domain
    /*! \todo Needs LabelSet specifics like Mean()*/
    real IntegralEstimate()
    { 
      return EnvelopeIntegral*real(int(Samples.size()))/real(Nprop);
    }
    /*
    rvector Mean(){
      vector<LabPnt>::const_iterator it = Samples.begin();
      rvector mean = it->Pnt;
      for(it++; it != Samples.end(); it++){ mean += it->Pnt; }
      mean /= (real)(int)Samples.size();
      return mean;
    }
    */
    
    //! Arithmetic mean of the sampled labeled points in a label-specific way 
    std::vector<rvector> Mean()
    {
      std::cout << "   Number of labels or topologies = " << Ntopologies << std::endl;
      std::vector<int>::const_iterator itINTV = 
        max_element(LabelSet.begin(),LabelSet.end());
      int MaxLabelNum = *itINTV;
      
      //! Flag for IF a label has been encountered in the samples
      std::vector<bool> first(MaxLabelNum+1, true);
      
      //! number of distinct labels in the samples
      std::vector<long> L_sums(MaxLabelNum+1, 0);
                    
      //! sum of distinct labels in the samples
      std::vector<rvector> sums(MaxLabelNum+1);
      /*! \todo Either replace with gsl_mean like computations due to their 
        diff eqns form or Kahan Summations Function Obj using std::transform 
        or DotAccum in c-xsc 
      */
      std::vector<LabPnt>::const_iterator it = Samples.begin();
      for(; it != Samples.end(); it++)
      {
        int label = it->L;
        if(first[label])
        {
          sums[label] = rvector(it->Pnt);
          L_sums[label] = 1;
          first[label] = false;
        }
        else
        {
          sums[label] += it->Pnt;
          L_sums[label] += 1;
        }
      }
      for(int i=0; i<Ntopologies; i++)
      {
        if(L_sums[LabelSet[i]] > 0)
        { 
          sums[LabelSet[i]] /= (real)L_sums[LabelSet[i]];
        }
        std::cout << "label: " << LabelSet[i] << "  proportion: " 
             << (real)L_sums[LabelSet[i]]/(real)(long)Samples.size() 
             << std::endl << "Labelled Mean:\n" << sums[LabelSet[i]] << std::endl;
      }
      return sums;
    }
    
    //! Arithmetic mean of the sampled labeled points, proportions and label: label-specific way 
    MeansLabelsProportions MeanLabelProportion()
    {
      bool printOut = false;
      MeansLabelsProportions mlp;
      if (printOut) std::cout << "   Number of labels or topologies = " << Ntopologies << std::endl;
      std::vector<int>::const_iterator itINTV = 
        max_element(LabelSet.begin(),LabelSet.end());
      int MaxLabelNum = *itINTV;
      
      //! Flag for IF a label has been encountered in the samples
      std::vector<bool> first(MaxLabelNum+1, true);
      
      //! number of distinct labels in the samples
      std::vector<long> L_sums(MaxLabelNum+1, 0);
                    
      //! sum of distinct labels in the samples
      std::vector<rvector> sums(MaxLabelNum+1);

      //! model labels in the samples
      std::vector<int> ModelLabels(Ntopologies);

      //! model proportions in the samples
      std::vector<real> ModelProportions(Ntopologies,0.0);

      /*! \todo Either replace with gsl_mean like computations due to their 
        diff eqns form or Kahan Summations Function Obj using std::transform 
        or DotAccum in c-xsc 
      */
      std::vector<LabPnt>::const_iterator it = Samples.begin();
      for(; it != Samples.end(); it++)
      {
        int label = it->L;
        if(first[label])
        {
          sums[label] = rvector(it->Pnt);
          L_sums[label] = 1;
          first[label] = false;
        }
        else
        {
          sums[label] += it->Pnt;
          L_sums[label] += 1;
        }
      }
      for(int i=0; i<Ntopologies; i++)
      {
        if(L_sums[LabelSet[i]] > 0)
        { 
          sums[LabelSet[i]] /= (real)L_sums[LabelSet[i]];
        }
        real prop_i =  (real)L_sums[LabelSet[i]]/(real)(long)Samples.size(); 
        ModelProportions[i]= prop_i;
        ModelLabels[i]=LabelSet[i];
        if (printOut) {
          std::cout << "label: " << LabelSet[i] << "  proportion: " 
             << prop_i
             << std::endl << "Labelled Mean:\n" << sums[LabelSet[i]] << std::endl;
        }
      }
      mlp.means = sums; mlp.labels=ModelLabels; mlp.proportions=ModelProportions;
      return mlp;
    }
    
    //! Print label-specific sample means from Mean().
    void PrintMeans(std::ostream& out)
    {
      std::vector<rvector> means = Mean();
      std::vector<rvector>::const_iterator it;
      for(size_t i=0; i<means.size(); i++)
      {
        out << "label: " << LabelSet[i] << "   mean: " 
            << means[LabelSet[i]] << std::endl;
      }
    }
    
    //! Print sampled labeled points in Samples as a matrix with TAB padding.
    void Print(std::ostream& out)
    {
      std::vector<LabPnt>::const_iterator it = Samples.begin();
      for(; it!=Samples.end(); it++)
      {
        out << it->L;
        int Current_Dim=Ub(it->Pnt)-Lb(it->Pnt)+1;
        for(int i=Lb(it->Pnt); i<=Ub(it->Pnt); i++)
        {
          out << '\t' << (it->Pnt[i]);
        }
        if(n_Dim_Max > Current_Dim) 
        {
          for(int i=Current_Dim+1; i<=n_Dim_Max; i++) out << '\t';
        }
        out << std::endl;
      }
      return;
    }

    //! Return a copy of the Samples of sampled labeled points.
    std::vector<LabPnt> getSample()
    {
      return Samples;
    }

};                  // end of RSample class

//! A class for the status of an Importance Sampler
class ISSample
{
  public:
    //! Number of topologies or model labels
    long Ntopologies;
    
    //! The set of unique integer labels in LabDomainList
    std::vector<int> LabelSet;
    
    //! The integral of the envelope function
    real EnvelopeIntegral;
    
    //! vector of weighted labeled points -- our labeled point-valued particles
    std::vector<WLabPnt> Samples;
    
    /*! \todo Needs LabelSet specifics like Mean()*/
    real IntegralEstimate(std::vector<real>& IntegralEsts)
    {
      /*! \todo replace with gsl_mean like computations due to their diff 
        eqns form or Kahan Summations 
      */
      std::vector<WLabPnt>::const_iterator it = Samples.begin();
      real wsum = it->Wt;
      for(; it != Samples.end(); it++)
      {
        int label = it->L;
        real weight = it->Wt;
        IntegralEsts[label] += weight;
        wsum += weight;
      }
      for(int i=0; i<Ntopologies; i++)
      {
        IntegralEsts[i] /= (real)(int)Samples.size(); 
      }
      return wsum/(real)(int)Samples.size();
    }
    
    //! Sample mean of lableled points.
    std::vector<rvector> Mean()
    {
      std::cout << "   Number of labels or topologies = " << Ntopologies << std::endl;
      std::vector<int>::const_iterator itINTV = 
        max_element(LabelSet.begin(),LabelSet.end());
      int MaxLabelNum = *itINTV;
      
      //! Flag for IF a label has been encountered in the samples
      std::vector<bool> first(MaxLabelNum+1, true);
                    
      //! number of distinct labels in the samples
      std::vector<real> w_sums(MaxLabelNum+1, 0);
                    
      //! sum of distinct labels in the samples
      std::vector<rvector> sums(MaxLabelNum+1);
      /*! \todo Either replace with gsl_mean like computations due to their 
        diff eqns form or Kahan Summations Function Obj using std::transform 
        or DotAccum in c-xsc */
      
      std::vector<WLabPnt>::const_iterator it = Samples.begin();
      for(; it != Samples.end(); it++)
      {
        int label = it->L;
        if(first[label])
        {
          sums[label] = rvector(it->Pnt);
          sums[label] *= it->Wt;
          w_sums[label] = it->Wt;
          first[label] = false;
        }
        else
        {
          sums[label] += it->Wt*it->Pnt;
          w_sums[label] += it->Wt;
        }
      }
      for(int i=0; i<Ntopologies; i++)
      {
        if(w_sums[LabelSet[i]] > 0){ sums[LabelSet[i]] /= w_sums[LabelSet[i]]; }
        std::cout << "label: " << LabelSet[i] << "  weight: " 
             << w_sums[LabelSet[i]] << std::endl << "mean: " << sums[LabelSet[i]] 
             << std::endl;
      }
      return sums;
    }
    
    //! Print labeled means.
    void PrintMeans(std::ostream& out)
    {
      std::vector<rvector> means = Mean();
      std::vector<rvector>::const_iterator it;
      for(int i=0; i<Ntopologies; i++)
      {
        out << "label: " << LabelSet[i] << "   mean: " 
            << means[LabelSet[i]] << std::endl;
      }
    }
    
    //Print weighted labeled points.
    void Print(std::ostream& out)
    {
      std::vector<WLabPnt>::const_iterator it = Samples.begin();
                    // it->Print(out); }
      for(; it!=Samples.end(); it++)
      {
        out << it->Wt << "     " << it->L; 
        for(int i=Lb(it->Pnt); i<=Ub(it->Pnt); i++){ out << " " << it->Pnt[i]; }
        out << std::endl;
      }
      return;
    }
};
#endif

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

/*! \file      Fobj.hpp
\brief Fobj definition and declaration.
*/

#ifndef __FOBJ__
#define __FOBJ__

#include <gop.hpp>
#include "SmallClasses.hpp"

/*! \brief An abstract class for target function objects

For use with MRS, for sampling from a unnormalized density
which is typically a likelihood times a prior density.  We know how to normalize
the prior density over any compact box.  
Currently: either unif (PriorType == 0) or exp(-x)
The function object operator () returns the likelihood, or
if UsingLogDensity is true, the log of the LogLikelihood. 

Also for use with C-XSC Global optimization routines GOptMin() and GOptMax(), 
in the context of rigorous maximum lkelihood estimation over labeled domains.
*/

class Fobj
{
  protected:
  
    //! a flag for working on the log(target shape) scale                
    bool UsingLogDensity;
  
    //! The initial collection of labelled domain boxes -- prior support                
    vector<LabBox> LabDomainList;
    
    //! The prior integral over the list of labeled domain boxes
    vector<real> LabDomainPriorIntegralList;

    //! To specify a type of prior density: uniform, exponential, user_defined.
    int PriorType;

  public:
    // ******************* pure virtual functions
    //! Get number of interval function calls
    virtual int get_interval_calls() = 0;

    //! Get number of real function calls
    virtual int get_real_calls() = 0;

    // Get number of GradType function calls
    //virtual int get_gradtype_calls() = 0;

    /*! \brief a defined pure virtual function for default Lebesgue measure 
      (volume) of a labelled box
      
      derived classes inherit interface and default behavior 
      see [Item 34,P.166, Meyers, Effective C++]
    */
    virtual real LabBoxVolume(const LabBox& LB) = 0;

    /*! \brief a pure virtual function for interval image of boxes under Fobj
    
      Derived classes inherit interface only 
      see [Item 34,P.162, Meyers, Effective C++] 
    */
    virtual interval operator()(const LabBox& LB) const = 0;

    //! a pure virtual function for real image of real points under Fobj
    virtual real operator()(const LabPnt& LP) const = 0;

    // *******************  virtual functions but not pure virtual
    //! Destructor
    virtual ~Fobj(){}

    /*! \brief a virtual function for HessType image of HTvector under Fobj
      
      Matches signature for HTscalar_FctPtr type for Global optimization 
      using AllGOp and model label if a label is specified
    */
    virtual HessType operator()(const HTvector& x, const int label = 0) const;

    //! Integral over box of the prior
    virtual real LabBoxPriorIntegral(const LabBox& LB);

    //! set up the prior integral over the list of labeled domain boxes
    virtual void SetupLabDomainPriorIntegralList()
    {
      LabDomainPriorIntegralList.resize(LabDomainList.size());
      for(size_t i=0; i<LabDomainList.size(); i++)
      {
        LabDomainPriorIntegralList[i] 
          = Fobj::LabBoxPriorIntegral(LabDomainList[i]);
      }
    }

    //! Draw a real point in labeled box from density proportional to prior
    virtual rvector DrawFromBoxPrior(const LabBox& LB, const rvector& randvec)
    {
      int lo = Lb(LB.Box);
      int hi = Ub(LB.Box);
      rvector the_point(lo, hi);
      if((hi - lo) >  (Ub(randvec) - Lb(randvec)))
      {
        cerr 
          << "in Fobj::DrawFromBoxPrior. Too few random numbers in randvec. \n";
        exit(1);
      }

      for(int i=lo; i<=hi; i++)
      {
        the_point[i] = inv_cdf(LB.Box[i], randvec[i]);
      }

      return the_point;
    }

    //! prior density of a labeled real point
    virtual real PriorDensity(const LabPnt& LP)
    {
      rvector Pnt = LP.Pnt;

      if(PriorType == 0)
      {
        return 1.0;
      }
      else
      {
        // -1*log of prior density
        real mlpd = 0.0;

        for(int i=Lb(Pnt); i<=Ub(Pnt); i++)
        {
          mlpd -= Pnt[i];
        }

        return exp(mlpd);
      }
    }

    /*! \brief return the point at which to bisect a labeled box along 
      dimension k
      
      Bisect labeled box along k-th dimension so as to get equal prior 
      probabilities in the resulting two halves.
    */
    virtual real BisectPt(LabBox& LB, int k)
    { 
      // two new boxes should each get half the prior prob.
      real x2 = inv_cdf(LB.Box[k], 0.5);
      //if(x1-x2 != 0.0){
      //cout << "bisect pt: x1, x2: " << x1 
      //  << "    " << x2 << "   " << (x1-x2) << endl;
      //}
      return x2;
    }

    /*! \brief return real x s.t. that fraction R of probability lies below x

      for prior density (either unif (PriorType == 0) or exp(-x)) normalized to 
      1 in box
    */
    virtual real inv_cdf(interval I, real R)
    {
      // cout << "In inv_cdf. PriorType: " << PriorType << endl;
      return (PriorType == 0)?
      // Inf(I) + R*(diam(I))
      // better than the line above? agrees with mid(I) when R=0.5
        (1.0 - R)*Inf(I) + R*Sup(I)
        : Inf(I) - ln(1.0 - R*(1.0 - exp(-diam(I))));
    }

    // ***************** non-virtual functions

    //! set the target scale to natural logarithm
    void setUsingLogDensity(bool LogPi)
    {
      UsingLogDensity = LogPi;
    }

    //! get the target scale being used
    bool getUsingLogDensity()
    {
      return UsingLogDensity;
    }

    //! Get the dimensions of the list of labeled domains
    int getLabeledDomainDim(int label) const
    {
      int retValue = 0;
      bool foundLabel = 0;

      vector<LabBox>::const_iterator it;

      it = LabDomainList.begin();

      while(!foundLabel && it < LabDomainList.end())
      {
        if (it->L == label)
        {
          foundLabel = 1;
          retValue = Ub(it->Box) - Lb(it->Box) + 1;
        }
        it++;
      }
      return retValue;
    }

    //! Get the list of labeled domains
    vector<LabBox> get_domain() const
    {
      return LabDomainList;
    }

    //! Get the set of unique integer labels in LabDomainList
    vector<int> get_labelset()
    {
      set<int> LS;

      for(vector<LabBox>::const_iterator it = LabDomainList.begin();
          it != LabDomainList.end(); it++) LS.insert(it->L);
      vector<int> LV;
      for(set<int>::const_iterator it = LS.begin(); 
          it != LS.end(); it++) LV.push_back(*it);
      return LV;
    }

};                  // end of Fobj declarations and definitions

// ******************************************* Fobj definitions ****************

/*! \brief get the default volume (Lebesgue) of a labeled box 

  an inline definition of the pure virtual function for default Lebesgue 
  measure (volume) of a labeled box
*/
inline real Fobj::LabBoxVolume(const LabBox& LB)
{
  real volume = 1.0;
  for (int i = 1; i <= VecLen(LB.Box); i++) volume *= diam(LB.Box[i]);
  return volume;
}

/*! \brief get the default integral of prior over a labeled box

  an inline definition of the virtual function for default integral of prior 
  over a labeled box -- not yet pure virtual
*/
inline real Fobj::LabBoxPriorIntegral(const LabBox& LB)
{
  real integral = 1.0;
  for (int i = 1; i <= VecLen(LB.Box); i++)
  {
    integral *= (PriorType == 0)? 
      diam(LB.Box[i]): exp(-Inf(LB.Box[i])) - exp(-Sup(LB.Box[i]));
  }
  return integral;
}

//! an inline definition for HessType image of HTvector under Fobj
/* get rid of warnings from -Wextra by changing signature to 
 * inline HessType Fobj::operator()(const HTvector&, const int) const */
inline HessType Fobj::operator()(const HTvector& x, const int label) const
{
  HessType hh(1);
  hh = 0.0;
  return hh;
}

// ************************************ Fobj1D definition and declaration ******

/*! \brief abstract class for one-dimensional function objects

  This is mainly used to speedup 1D targets and not supported well.
*/
class Fobj1D
{
  protected:

    bool UsingLogDensity;

    vector<LabBox> LabDomainList;

    vector<real> LabDomainPriorIntegralList;

  public:
    // Destructor
    virtual ~Fobj1D(){}

    void setUsingLogDensity(bool LogPi)
    {
      UsingLogDensity = LogPi;
    }

    virtual interval operator()(const interval& LB) const = 0;

    virtual real operator()(const real& LP) const = 0;

    bool getUsingLogDensity()
    {
      return UsingLogDensity;
    }

    vector<LabBox> get_domain() const
    {
      return LabDomainList;
    }
};                  // end of FOBJ1D declarations and definitions
#endif

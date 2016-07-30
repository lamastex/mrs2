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
/*! \file FCFN3.hpp
\brief Trans-dimensional three-taxa Cavender-Farris-Neyman phylogenetic
model declarations.
*/

#ifndef __FCFN3_HPP__
#define __FCFN3_HPP__

/*! \brief 1-dimensional Cavender-Farris-Neyman (CFN) model likelihood
  as a function object class.

  This is the ultrametric "star tree", a trifurcating tree with all three 
  branches having equal length.
*/
class FCFN3Star: public Fobj
{
  int Cid;    //!< Cid = count of invariant sites (minimal suffcient stats)
  int Cnid;   //!< Cnid = count of variant sites (minimal suffcient stats)
  real f0;    //!< Fraction of invariant sites, f0 := Cid/(Cid+Cnid)

  /*! \brief Dimensions of the labelled Domain,

  ie, number of distinct branch lengths that have this topology label.
  */
  int n_dimensions;
  real TotSites;    //!< Total number of sites
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  //   vector<LabBox> LabDomainList;
  public:
    //! A constructor.
    FCFN3Star(int Cid, int Cnid, interval Domain, bool LogPi, int Prior);
    //   vector<LabBox> get_domain();
    //! interval function object
    interval operator()(const LabBox& X) const;
    //! real function object
    real operator()(const LabPnt& X) const;
    /*! Volume of rooted tree boxes is implemented here 
      and NOT inherited from Fobj
    */
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};

/*! \brief 3-dimensional Cavender-Farris-Neyman (CFN) model likelihood 
  as a function object class
*/
class FCFN3UnRooted: public Fobj
{
  int Cxxx;      //!< Count of sites with pattern xxx (minimal suffcient stats)
  int Cxxy;      //!< Count of sites with pattern xxy (minimal suffcient stats)
  int Cyxx;      //!< Count of sites with pattern yxx (minimal suffcient stats)
  int Cxyx;      //!< Count of sites with pattern xyx (minimal suffcient stats)
  real fxxx;     //!< Fraction of sites with pattern xxx, Cxxx/TotSites
  real fxxy;     //!< Fraction of sites with pattern xxy, Cxxy/TotSites
  real fyxx;     //!< Fraction of sites with pattern yxx, Cyxx/TotSites
  real fxyx;     //!< Fraction of sites with pattern xyx, Cxyx/TotSites
  
  /*! \brief Dimensions of the labelled Domain,
    
    ie, number of distinct branch lengths that have this topology label.
  */  
  int n_dimensions; 
  
  real TotSites;    //!< Total number of sites
  interval PositiveProbInterval;
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  //   vector<LabBox> LabDomainList;
  //! Track number of grad_type function calls
  mutable int n_hesstype_calls;
  public:
    //! A constructor.
    FCFN3UnRooted(int Cxxx, int Cxxy, int Cyxx, int Cxyx, 
                  interval Domain, bool LogPi, int Prior);
    //   vector<LabBox> get_domain();
    //! interval function object
    interval operator()(const LabBox& X) const;
    //! real function object
    real operator()(const LabPnt& X) const;

    ///! Gradfunction object
    //GradType operator()(const GTvector& x, const int label = 0) const;

    //! HessType function object
    HessType operator()(const HTvector& x, const int label = 0) const;

    /*! Volume of rooted tree boxes is implemented here 
      and NOT inherited from Fobj
    */
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};

/*! \brief 2-dimensional Cavender-Farris-Neyman (CFN) model likelihood as a 
function object class
*/
class FCFN3Rooted: public Fobj
{
  int Cxxx;      //!< Count of sites with pattern xxx (minimal suffcient stats)
  int Cxxy;      //!< Count of sites with pattern xxy (minimal suffcient stats)
  int Cyxx;      //!< Count of sites with pattern yxx (minimal suffcient stats)
  int Cxyx;      //!< Count of sites with pattern xyx (minimal suffcient stats)
  real fxxx;     //!< Fraction of sites with pattern xxx, Cxxx/TotSites
  real fxxy;     //!< Fraction of sites with pattern xxy, Cxxy/TotSites
  real fyxx;     //!< Fraction of sites with pattern yxx, Cyxx/TotSites
  real fxyx;      //!< Fraction of sites with pattern xyx, Cxyx/TotSites

  
  /*! \brief Dimensions of the labelled Domain,
    
    ie, number of distinct branch lengths that have this topology label.
  */  
  int n_dimensions; 
  real TotSites;    //!< Total number of sites
  interval PositiveProbInterval;
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  //   vector<LabBox> LabDomainList;
  public:
    //! A constructor.
    FCFN3Rooted(int Cxxx, int Cxxy, int Cyxx, int Cxyx, 
                interval Domain, bool LogPi, int Prior);
    //   vector<LabBox> get_domain();
    //! interval function object
    interval operator()(const LabBox& X) const;
    //! real function object
    real operator()(const LabPnt& X) const;
    /*! Volume of rooted tree boxes specialised in our embedding 
      (implemented and NOT inherited)
      */
    virtual real LabBoxVolume(const LabBox& LB);
                    //!< Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
                    //!< Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};

/*! \brief 1,2,3-trans-dimensional Cavender-Farris-Neyman (CFN) model 
  likelihood as a function object class
*/
class FCFN3: public Fobj
{
  int Cxxx;    //!< Count of sites with pattern xxx (minimal suffcient stats)
  int Cxxy;    //!< Count of sites with pattern xxy (minimal suffcient stats)
  int Cyxx;    //!< Count of sites with pattern yxx (minimal suffcient stats)
  int Cxyx;    //!< Count of sites with pattern xyx (minimal suffcient stats)
  //! Cid = count of invariant sites (minimal suffcient stats for star tree)
  int Cid;          
  //! Cnid = count of variant sites (minimal suffcient stats for star tree)
  int Cnid;         
  real f0;     //!< Fraction of invariant sites, f0 := Cid/(Cid+Cnid)
  real fxxx;   //!< Fraction of sites with pattern xxx, Cxxx/TotSites
  real fxxy;   //!< Fraction of sites with pattern xxy, Cxxy/TotSites
  real fyxx;   //!< Fraction of sites with pattern yxx, Cyxx/TotSites
  real fxyx;   //!< Fraction of sites with pattern xyx, Cxyx/TotSites
  real TotSites; //!< Total number of sites
  //! Representable interval without complete interval arithmetic
  interval PositiveProbInterval;
  //! Volume of each labelled domain box
  vector<real> vol_labelled_domain;
  //! Track number of interval function calls
  mutable int n_interval_calls;
  //! Track number of real function calls
  mutable int n_real_calls;
  //   vector<LabBox> LabDomainList;
  public:
    //! A constructor.
    FCFN3(int Cxxx, int Cxxy, int Cyxx, int Cxyx, 
          interval Domain, bool LogPi, int Prior);
    //   vector<LabBox> get_domain();
    //! interval function object
    interval operator()(const LabBox& X) const;
    //! real function object
    real operator()(const LabPnt& X) const;
    /*! Volume of rooted tree boxes specialised in our embedding 
      (implemented and NOT inherited)
      */
    virtual real LabBoxVolume(const LabBox& LB);
    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }
    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }
};
#endif

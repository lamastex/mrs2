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

#ifndef __FWITCHESHAT_HPP__
#define __FWITCHESHAT_HPP__

//! n-dimensional witch's hat density as a function object class
class FWitchesHat: public Fobj
{
  real PeakMean, PeakRadius, PeakWeight, Rtod;
  real UniformWeight;
  real DomainVolume;
  mutable int n_interval_calls;
  mutable int n_real_calls;
  real h[11];

  //   vector<LabBox> LabDomainList;
  public:
    FWitchesHat(int n_dimensions, real PMean, real PRadius, real PWeight, 
                real UWeight, real DomainLimit, bool LogPi,int Prior);
    //   vector<LabBox> get_domain();
    interval operator()(const LabBox& X) const;
    real operator()(const LabPnt& X) const;
    virtual real LabBoxVolume(const LabBox& LB){return Fobj::LabBoxVolume(LB);}
    int get_interval_calls(){ return n_interval_calls; }
    int get_real_calls(){ return n_real_calls; }
};
#endif

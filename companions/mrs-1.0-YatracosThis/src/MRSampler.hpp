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

/*! \file:     MRSampler.hpp
\brief Moore Rejection Sampler (MRSampler) declarations

*/

#ifndef __MRSAMPLER_HPP__
#define __MRSAMPLER_HPP__

/*! \brief The Moore rejection sampler class for trans-dimensional targets 
  over labeled metric spaces.

\todo In MRSampler everything is inline for now -- we need to make this 
  truly object oriented, add a default constructor, add mrs namespace, 
  get rid of @#defines, etc. Needs about 20 hours for complete documentation...
*/
class MRSampler
{ 
  public:
    //! Initialised constructor
    MRSampler (Fobj & f, int max_n_boxes, double Alb, 
               unsigned seed = gsl_rng_default_seed, bool use_f_scale = true);
  
    //! Destructor
    ~MRSampler ();
    
    //! Return lower bound on the acceptance prob.
    double getPALB ();
    
    /*! \brief Do further partitioning until 
        acc. prob. lower bound > Alb, then setup pdf.
      */
    void Refine (double Alb);
    
    /*! \brief Refine partition by bisections until Desired_N_boxes many 
      boxes is reached, then setup pdf.
    */
    void RefineUntil (unsigned int Desired_N_boxes);
    
    //! Refine partition by doing Nbisect many bisections.
    void Refine (int Nbisect);
    
    //! Return the number of boxes in RangeDomainSet.
    int get_n_boxes ();
    
    /*! \brief Print labeled boxes in domain partition 
      DomainParts [C-XSC output format].
    */
    void Print_Domain_Partition (std::ostream& out);
    
    /*! \brief Print labeled boxes in domain partition 
      DomainParts [naive TAB-delimited numeric only format].
    */    
    void Output_Domain_Partition (std::ostream& out);

    //! Print the RangeDomainSet in tab-delimited numeric only format.
    std::ostream& MRSoutput(std::ostream &os, const double eps = 0) const;

    //! Return one sample of labeled point via rejection sampling, if possible.
    LabPnt RejectionSampleOnce (int &tries);
    
    /*! \brief Draw nRS many samples of labeled points via rejection sampling, 
      if possible.
      
      RejectionSampleMany stores the samples in the RSSample object theSample.
    */
    void RejectionSampleMany (int nRS, RSSample & theSample);
    
    //! Importance sampling with Quasi random numbers -- [Ignore: experimental]
    void ImportanceSampleManyQuasi (int NSamples, bool residual,
      ISSample & theSample);
    
    //! Importance sampling with Pseudo/Quasi random numbers
    void ImportanceSampleMany (int NSamples, bool residual, bool pseudoRNG,
      ISSample & theSample);
    
    //! Importance sampling with Pseudo random numbers
    void ImportanceSampleManyPseudo (int NSamples, bool residual,
      ISSample & theSample)
    {
      ImportanceSampleMany (NSamples, residual, true, theSample);
    }
    
    /*! Print boxes with MATLAB. 
      
      \todo Needs standardization of rendering format(s) 
      for ease of making low-dimensional pictures -- MATLAB/POVRAY/MATPLOTLIB.
    */
    void PrintBoxes (int Nprint);
    
    
    double getIU ();
    double getIL ();
    double get_unscaled_IU ()
    {
      return _double (Sup (Integral));
    }
    double getIUminusL ();
    double getUmax ();
    double getPAest ();
    double get_wsum ();
    int get_nsum ();
    double get_wmax ();
    double get_wmin ();
    // double get_f_scale(){ return _double(f_scale); }
    void updateIntegral ();
    void updateUmax ();
    int get_n_topologies(){ return static_cast<int>(DomainLabelSet.size()); }

    int get_nonresidual_samples ()
    {
      return nonresidual_samples;
    }
  private:
    // private member functions.
    
    RangedLabBox getBoxREInfo (LabBox LBox);
    
    //! Initialize pseudo and Quasi Random Number Generators in GSL.
    void InitRNG (unsigned seed);
    
    void FirstBox ();
    
    /*! \brief Bisect the box at top of priority queue. 
    */
    void Bisect ();
    
    /*! \brief Adaptively partition domain by bisecting the most prioritised 
      labeled box.
      
      The labeled box in the current partition of the domain with largest 
      (prior_integral*diam(range enclosure)) is bisected
    */
    void AdaptPartition (double Alb);
    
    void SetupPDF ();
    // added SetupWalker to separate the gsl_discrete_struct independent 
    // of SetupPDF to accomodate
    // the number_of_samples-dependent residual sampling scheme
    void SetupWalker (unsigned int number_of_samples);
    // private data members
    Fobj & F;

    vector < LabBox > Domain;
    //! The set of unique integer labels in LabDomainList
    vector<int> DomainLabelSet;
    // n_dim_max is the max dimensionality of all the labelled domains.
    unsigned int n_dim_max;
    // vector < unsigned int > n_dimensions;
    //  map< int, unsigned int > n_dimensions;
    unsigned int topologies;
    unsigned int Max_n_boxes;
    // unsigned int Number_Of_Samples;
    int nonresidual_samples;
    interval Integral;

    vector < LabBox > DomainParts;
    rvector UBox;   // Sup of range enclosure of f  (pi or ln(pi)) in box
    rvector LoBox;  // Inf of range enclosure of f  (pi or ln(pi)) in box

    real wmin, wmax, wsum;
    int nsum;
    // Umax is best upper bound to pi or ln(pi) at present; 
    // Lmax is max of lower bound function (updated every bisection)
    real Umax, Lmax, fMaxLB;

    real f_scale;
                    // set to true when Umax - fmid_max small enough;
    bool f_scaleDone;
    bool UseFScale;
    bool UsingLogDensity;

    bool own_rng;
    gsl_rng *rgsl;
    gsl_qrng *qrgsl;
    gsl_ran_discrete_t *gslpdfstruct;
    double *proposalpmf;
    double *residual_proposalpmf;
    double *proposalpdf;
    // int* nonresiduals;
    RangedLabBoxSet RangeDomainSet;
};

// non member function for the MRSampler class
// output operatore for MRSampler
// uses MRSoutput public member function
std::ostream & operator<<(std::ostream &os, const MRSampler& mrs);
#endif

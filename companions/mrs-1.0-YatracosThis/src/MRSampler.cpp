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

/*! \file:     MRSampler.cpp
\brief Moore Rejection Sampler (MRSampler) definitions

*/

// Including STD headers
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <string>
#include <math.h>
#include <getopt.h>
#include <time.h>

// Including interval arithmetic package C-XSC headers                    
#include "interval.hpp"
#include "imath.hpp"
#include "rmath.hpp"
#include "ivector.hpp"
#include "rvector.hpp"
#include "imatrix.hpp"

// Including GSL headers
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <functional>
#include <numeric>
#include <assert.h>

using namespace std;
using namespace cxsc;

// Including MRS headers
#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "toolz.hpp"
#include "MRSampler.hpp"

#define USEFSCALE 1
#define MAXNTRIES 100000

// when max of f is known to lie in a interval x 
// with log(diam(x)) < this, set scale and stop updating Umax
#define LOGDIAMFMAX (25.0)
// scale so Umax has this value
#define UmaxMAX (1.0e200)
#define RS_SQUEEZE 1
// for doing integration, get integral of lower bound 
// as sum(L_i*V_i), plus IS estimate of integral of (U - L)
#define WEIGHTED_SQUEEZE 0
// could put in neg infinity here in extended interval arithmetic
// at the expense of speed 
#define BIGNEGATIVE (-1.0e-307)
//#define DO_F_SCALING 1

// implementation of class MRSampler
// public member functions
MRSampler::MRSampler (Fobj & f, int max_n_boxes, double Alb, 
                      unsigned seed, bool use_f_scale)
:
// target density or shape
F (f),
// should be vector of LabBox
Domain (f.get_domain ()),
// should be corresponding labels
DomainLabelSet (f.get_labelset ()),
n_dim_max(0),
Max_n_boxes (max_n_boxes),
f_scale(0),
f_scaleDone(false),
UseFScale(use_f_scale),
own_rng (false),
gslpdfstruct (NULL),
proposalpmf (NULL),
residual_proposalpmf (NULL),
proposalpdf (NULL),
RangeDomainSet (SortBox (F.getUsingLogDensity ()))
{                   // constructor

  Umax = BIGNEGATIVE;
  wsum = 0.0;
  wmax = -1.0;
  wmin = 1.0e100;
  nsum = 0;

  vector < LabBox >::const_iterator it = Domain.begin();
  //n_dimensions = 1 + Ub (it->Box) - Lb (it->Box);

  for(it = Domain.begin(); it != Domain.end(); it++)
  {
    //int label = it->L;
    unsigned int n_dimensions = 1 + Ub (it->Box) - Lb (it->Box);
    //   cout << "label, n_dimensions[label]: " << label << " " 
    //          <<  n_dimensions[label] << endl;
    if(n_dimensions > n_dim_max){ n_dim_max = n_dimensions; }
  }

                    // Domain.size()
  topologies = static_cast<int>(DomainLabelSet.size());
  //cout << "number of topologies = " << topologies << endl;
  UsingLogDensity = F.getUsingLogDensity ();
  nonresidual_samples = 0;

  //  cout << "#x, f(x) " << endl;
  //   ofstream outx ("fout");
  //   for(int i=1; i<1000; i++){
  //       LabPnt P;
  //       rvector x (1, n_dimensions);
  //       x[1] = i*0.001;
  //       P.Pnt = x;
  //       P.L = 0;
  //       outx << P.Pnt[1] << "  " << F(P) << "   " << exp(F(P)) << endl;
  //   }
  //   cout << "after printing fout" << endl;
  InitRNG (seed);   // initialize the random number generator
  FirstBox ();
  cout << "after FirstBox, before Refine " << endl;
  Refine (Alb);
  cout << "after Refine " << endl;
  //cout << "making pdf structure" << endl;
  gslpdfstruct =
    gsl_ran_discrete_preproc ((size_t) RangeDomainSet.size (), proposalpmf);
}

// destructor
MRSampler::~MRSampler ()
{                   
  if (own_rng)
    gsl_rng_free (rgsl);
  gsl_qrng_free (qrgsl);
  gsl_ran_discrete_free (gslpdfstruct);
  free (proposalpmf);
  if (residual_proposalpmf != NULL)
  {
    free (residual_proposalpmf);
  }
  free (proposalpdf);
}

// Return lower bound on the acceptance prob.
double
MRSampler::getPALB ()
{
  return _double (Inf (Integral) / Sup (Integral));
}

// Do further partitioning until acc. prob. lower bound > Alb, then setup pdf
void
MRSampler::Refine (double Alb)
{
  AdaptPartition (Alb);
  SetupPDF ();
}

// Refine partition by bisections until Desired_N_boxes many boxes is reached, 
// then setup pdf
void
MRSampler::RefineUntil (unsigned int Desired_N_boxes)
{
  while (RangeDomainSet.size () < Desired_N_boxes)
  {
    Bisect ();
  }
  SetupPDF ();
}

// Refine partition by doing Nbisect many bisections
void
MRSampler::Refine (int Nbisect)
{
  for (int j = 0; j < Nbisect; j++)
  {
    Bisect ();
  }
  SetupPDF ();
}

// Return the number of boxes in RangeDomainSet, ie RangeDomainSet.size ()
int
MRSampler::get_n_boxes ()
{
  return static_cast<int>(RangeDomainSet.size ());
}

// Print the boxes in the Domain Partition to ostream out
void
MRSampler::Print_Domain_Partition (std::ostream& out)
{
  for (unsigned int ui = 0; ui < RangeDomainSet.size (); ui++)
    out << DomainParts[ui].L << endl << DomainParts[ui].Box << endl;
}

// Output the boxes in the Domain Partition to ostream out
void
MRSampler::Output_Domain_Partition (std::ostream& out)
{
  for (unsigned int ui = 0; ui < RangeDomainSet.size (); ui++)
  {
    out << DomainParts[ui].L;
    for(int i=Lb((DomainParts[ui]).Box); i<=Ub((DomainParts[ui]).Box); i++)
    {
      out << "\t" << Inf((DomainParts[ui]).Box[i]) 
          << "\t" << Sup((DomainParts[ui]).Box[i]);
    }
    out << endl;
  }
}

// Print the RangeDomainSet in tab-delimited numeric only format
/*! \todo may want the output to be padded with TABS for easy dlmread in MATLAB
  for the trans-diminsional case: same for Output_Domain_Partition
*/
std::ostream& MRSampler::MRSoutput(std::ostream &os, const double eps) const
{
  // do nothing if there is nothing in the set
  if(!RangeDomainSet.empty())
  {
    
    RangedLabBoxSet::const_iterator it;
    
    it = RangeDomainSet.begin();
    
    double vol = Volume((it->LBox).Box);
    
    while (vol>eps  && it != RangeDomainSet.end())
    {               // pull em out from top of pq
      RangedLabBox theBox = *it;
      ivector x = it->LBox.Box;
      
      // label
      os << (it->LBox.L);
      
      // range enclosure
      os << "\t" << Inf(it->BoxRE) << "\t" << Sup(it->BoxRE);
      
      //then the box
      for (int i = Lb(x); i <= Ub(x) ; i++)
      {
        os << "\t" << Inf(x[i]) << "\t" << Sup(x[i]);
      }
      
      os<<endl;
      
      vol = Volume((it->LBox).Box);
      
      it++;
      
    }               // end iteration through the set
    
  }                 // end if set not empty
  
  return os;
}

// Return one sample via rejection sampling, if possible
LabPnt
MRSampler::RejectionSampleOnce (int& tries)
{

  LabPnt proposed_LPnt;
  for (tries = 1; tries <= MAXNTRIES; tries++)
  {                 // try a bunch of times
    int proposed_index = 
                        static_cast<int>(gsl_ran_discrete (rgsl, gslpdfstruct));
    //  int ndim = VecLen(DomainParts[proposed_index].Box);
    rvector proposed_point;
    if(0)           // old way in mrs-0.1; uniform prior -- left for recall!!!
    {
      //Different boxes can be different dimensionalities.
      proposed_point = DrawUnifBox (rgsl, DomainParts[proposed_index].Box);
    }
    else
    {
      LabBox LBox = DomainParts[proposed_index];
      int BoxDim = Ub(LBox.Box) - Lb(LBox.Box) + 1;
      rvector rand_vector(1,BoxDim);
      for(int k=1; k<=BoxDim; k++)
      {
        rand_vector[k] = gsl_rng_uniform(rgsl);
      }

      proposed_point = F.DrawFromBoxPrior(DomainParts[proposed_index], 
                                          rand_vector);
    }
    proposed_LPnt.Pnt = proposed_point;
    proposed_LPnt.L = DomainParts[proposed_index].L;
    proposed_LPnt.Print(cout);

    real rand = gsl_rng_uniform (rgsl);
    real Fprop;
    if (rand > 1.0)
    {
      printf
        ("#proposed_index, UBox[proposed_index], height: %i %g %g \n",
        proposed_index, _double (UBox[proposed_index]), _double (rand));
      getchar ();
    }

    if (UsingLogDensity)
    {
      // Ah but Ubox is now the log
      //   cout << "RSonce: rany: " << rany << "  Ubox: " 
      //        << UBox[proposed_index] <<endl;
      if (((RS_SQUEEZE) && 
           (ln(rand) <= LoBox[proposed_index] - UBox[proposed_index]))
        || (ln(rand) <= F(proposed_LPnt) - UBox[proposed_index]))
      {
        return proposed_LPnt;
      }
    }
    else
    {
      real rany = rand * UBox[proposed_index];
                    // < lower bound, don't need to eval. function
      if (((RS_SQUEEZE) && (rany <= LoBox[proposed_index]))
        || (rany <= (Fprop = F (proposed_LPnt))))
      {
        return proposed_LPnt;
      }
    }
  }                 // loop over tries
  cerr << "In MRSample.SampleOnce. After " << MAXNTRIES
    << " proposals, none accepted. Acceptance prob. very low. " << endl;
  cout << "estimated acceptance prob: " << getPAest () << endl;
  getchar ();
  return proposed_LPnt;
}

//! Draw nRS many rejection samples, if possible, 
//! and store in the RSSample object theSample
void
MRSampler::RejectionSampleMany (int nRS, RSSample& theSample)
{
  LabPnt proposed_LPnt;
  int total_tries = 0;
  ////ofstream psampout("MRS_proposed.samples");// if you are healthily paranoid to see proposals!
  for (; theSample.Samples.size () < unsigned(nRS);)
  {
    int proposed_index = static_cast<int>(gsl_ran_discrete 
                                          (rgsl, gslpdfstruct));
    LabBox LBox = DomainParts[proposed_index];
    rvector proposed_point = 
              F.DrawFromBoxPrior(LBox, DrawUnifUnitBox(rgsl, VecLen(LBox.Box)));

    proposed_LPnt.Pnt = proposed_point;
    proposed_LPnt.L = DomainParts[proposed_index].L;
    ////proposed_LPnt.Print(psampout);

    real rand = gsl_rng_uniform (rgsl);
    real Fprop;
    if (rand > 1.0)
    {
      printf
        ("#proposed_index, UBox[proposed_index], height: %i %g %g \n",
        proposed_index, _double (UBox[proposed_index]), _double (rand));
      getchar ();
    }

    if (UsingLogDensity)
    {
      // real rany = rand * UBox[proposed_index];
      // cout << "RSmany. rany: " << rany << "  Ubox: " 
      //      << UBox[proposed_index] << endl;
      if (((RS_SQUEEZE) && 
           (ln(rand) <= LoBox[proposed_index] - UBox[proposed_index]))
        || ln (rand) <= F(proposed_LPnt) - UBox[proposed_index])
      {
        //  proposed_LPnt.Print(cout);
        theSample.Samples.push_back (proposed_LPnt);
      }
    }
    else
    {
      real rany = rand * UBox[proposed_index];
                    // less than lower bound, don't need to eval. function
      if (((RS_SQUEEZE) && (rany <= LoBox[proposed_index]))
        || (rany <= (Fprop = F (proposed_LPnt))))
      {
        // proposed_LPnt.Print(cout);
        theSample.Samples.push_back (proposed_LPnt);
      }
    }
    // } // loop over tries
    // if(tries == MAXNTRIES){ cout << "No RS accepted in " 
    //                    << MAXNTRIES << " tries. Refine partition?" << endl; }
    total_tries++;  // += tries;
    // cout << "RS accepted, tries: " << theSample.Samples.size() 
    //      << "  " << total_tries << endl;
  }                 // loop over samples
  theSample.Nprop = total_tries;
  theSample.Ntopologies = DomainLabelSet.size ();
  theSample.LabelSet = DomainLabelSet;
  theSample.n_Dim_Max = n_dim_max;
  theSample.EnvelopeIntegral = getIU ();
  // cerr << "In MRSample.SampleOnce. After " << MAXNTRIES
  //      << " proposals, none accepted. Acceptance prob. very low. " << endl;
  //  return proposed_LPnt;
}

void
MRSampler::ImportanceSampleMany (int NSamples, bool residual, 
                                 bool pseudoRNG, ISSample& theSample)
{
  WLabPnt Sample;
                    // n_dimensions); 
  // can we do this, and have dimensionality automatically increased 
  // as necessary by assignment?
  rvector proposed_point (1, 1);
  //double v[n_dimensions];
  cout << "n_dim_max: " << n_dim_max << endl;
                    // n_dimensions[0]);
  vector<double> v(n_dim_max);
  cout << "after alloc v " << endl;
  int NBoxes = int (RangeDomainSet.size ());

  if (residual)
  {
    SetupWalker (NSamples);
  }                 // will this work ok here?
  else
  {
    for (int i = 0; i < NBoxes; i++)
    {
      DomainParts[i].SamplesToDo = 0;
    }
  }                 //initialize

                    // N is number of residual samples
  int Nresidual = (NSamples - nonresidual_samples);
  for (int i = 0; i < Nresidual; i++)
  {                 // randomly choose the boxes for the residual samples
    DomainParts[(gsl_ran_discrete (rgsl, gslpdfstruct))].SamplesToDo++;
  }
  // now the number of samples to be drawn from each box is set, 
  // and stored in DomainParts[index].SamplesToDo
  size_t sample_num = 0;
  for (int i = 0; i < NBoxes; i++)
  {
    while (DomainParts[i].SamplesToDo > 0)
    {               //until all samples from Box i are exhausted
                    //decrement the sample that's about to be taken care of
      DomainParts[i].SamplesToDo--;
      sample_num++; //increment sample_num
      if(1)
      {
        if (pseudoRNG)
        {
          proposed_point = DrawUnifBox (rgsl, DomainParts[i].Box);
        }
        else
        {           // quasi RNG
                    // the quasi-random state is v at current sweep
          gsl_qrng_get (qrgsl, (&v[0]));
          // Note: QR vector v may have more dimensions than box; 
          // in this case don't use the extra elements.
          proposed_point = DrawQZUnifBoxV ((&v[0]), DomainParts[i].Box);
        }
      }
      else          //
      {
        LabBox LBox = DomainParts[i];
        //  rvector proposed_point;
        if (pseudoRNG)
        {
          proposed_point = 
            F.DrawFromBoxPrior(LBox, DrawUnifUnitBox(rgsl, VecLen(LBox.Box)));
          //	proposed_point = DrawUnifBox (rgsl, DomainParts[i].Box);
        }
        else
        {           // quasi RNG
                    // the quasi-random state is v at current sweep
          gsl_qrng_get (qrgsl, (&v[0]));
          // Note: QR vector v may have more dimensions than box; 
          // in this case don't use the extra elements.
          // get a rvector from a double* ???
          int ndim = VecLen(LBox.Box);
          rvector rand_vector(1, ndim);
          for(int ii=1; ii<=ndim; ii++){ rand_vector[ii] = v[ii-1]; }

          proposed_point = F.DrawFromBoxPrior(LBox, rand_vector);
          // proposed_point = DrawQZUnifBoxV ((&v[0]), DomainParts[i].Box);
        }
      }
      Sample.Pnt = proposed_point;
      Sample.L = DomainParts[i].L;
      if (UsingLogDensity)
      {
        cerr <<
          "UsingLogDensity=true is not defined for MRSampler::ImpPDSampleOnce()"
          << endl;
        exit(1);
      }
      else
      {
        Sample.fPnt = _double (F (Sample));
        if (WEIGHTED_SQUEEZE)
        {
          Sample.fPnt -= _double (LoBox[i]);
        }           // FIXME scaling??
                    //proposalpdf has been MODIFIED after pdfstruct was made
        Sample.qPnt = _double (proposalpdf[i] * F.PriorDensity(Sample));
        Sample.Wt = Sample.fPnt / Sample.qPnt;
      }
      theSample.Samples.push_back (Sample);
    }               //end while(DomainParts[i].SamplesToDo > 0)
  }                 // end i-loop
  //  Offset = (WEIGHTED_SQUEEZE)? _double(Inf(Integral)): 0.0;
  assert (sample_num == theSample.Samples.size ());
  theSample.Ntopologies = DomainLabelSet.size ();
  theSample.LabelSet = DomainLabelSet;
  theSample.EnvelopeIntegral = getIU ();
}                   // end ImportanceSampleMany

void
MRSampler::ImportanceSampleManyQuasi (int NSamples, bool residual, 
                                      ISSample & theSample)
{
  WLabPnt Sample;
  rvector proposed_point (1, 1);
  // the double vector to store the quasi-random state vector
  vector<double> v(n_dim_max);

  int NBoxes = int (RangeDomainSet.size ());

  if (residual)
  {
    SetupWalker (NSamples);
  }                 // will this work ok here?
  else
  {
    for (int i = 0; i < NBoxes; i++)
    {
      DomainParts[i].SamplesToDo = 0;
    }               //initialize
  }

                    // N is number of residual samples
  int Nresidual = (NSamples - nonresidual_samples);
  for (int i = 0; i < Nresidual; i++)
  {                 // randomly choose the boxes for the residual samples
    DomainParts[(gsl_ran_discrete (rgsl, gslpdfstruct))].SamplesToDo++;
  }
  // now the number of samples to be drawn from each box is set, 
  // and stored in DomainParts[index].SamplesToDo
  size_t sample_num = 0;
  for (int i = 0; i < NBoxes; i++)
  {
    while (DomainParts[i].SamplesToDo > 0)
    {               //until all samples from Box i are exhausted
                    //decrement the sample that's about to be taken care of
      DomainParts[i].SamplesToDo--;
      sample_num++; //increment sample_num
      if (0)
      {
        proposed_point = DrawUnifBox (rgsl, DomainParts[i].Box);
      }
      else
      {
                    // the quasi-random state is v at current sweep
        gsl_qrng_get (qrgsl, (&v[0]));
        proposed_point = DrawQZUnifBoxV ((&v[0]), DomainParts[i].Box);
      }

      Sample.Pnt = proposed_point;
      Sample.L = DomainParts[i].L;
      if (UsingLogDensity)
      {
        cerr <<
          "UsingLogDensity=true is not defined for MRSampler::ImpPDSampleOnce()"
          << endl;
        exit (1);
      }
      else
      {
        Sample.fPnt = _double (F (Sample));
        if (WEIGHTED_SQUEEZE)
        {
          Sample.fPnt -= _double (LoBox[i]);
        }           // FIXME scaling??
                    //proposalpdf has been MODIFIED after pdfstruct was made
        Sample.qPnt = _double (proposalpdf[i] * F.PriorDensity(Sample));
        Sample.Wt = Sample.fPnt / Sample.qPnt;
      }
      theSample.Samples.push_back (Sample);
    }               //end while(DomainParts[i].SamplesToDo > 0)
  }                 // end i-loop
  //  Offset = (WEIGHTED_SQUEEZE)? _double(Inf(Integral)): 0.0;
  //  theSample.N = sample_num;
  assert (sample_num == theSample.Samples.size ());
  theSample.Ntopologies = Domain.size ();
  theSample.EnvelopeIntegral = getIU ();
}

void
MRSampler::PrintBoxes (int Nprint)
{
  //this is quick and dirty Matlab code generator for visualising range enclosure in 1D cases
  Nprint = 0;       //output through cout
  int NBoxesToPrint=1000;
  RangedLabBox theBox;
  RangedLabBoxSet::const_iterator it = RangeDomainSet.begin ();
  for (int i = 0; it != RangeDomainSet.end () && i < NBoxesToPrint; it++, i++)
  {                 // pull em out from top of pq
    theBox = *it;
    ivector x = it->LBox.Box;
    //cout << "i: " << i << "  Box: " << it->LBox.Box << " RE: :" 
    //     << it->BoxRE << "  Vol: " << it->BoxVol;
    //cout << " IntDiam: " << it->BoxIntegral (UsingLogDensity) << endl;
    cout << "rectangle('Position',["
      << Inf(x[1]) << "," << Inf(it->BoxRE) << "," << Sup(x[1])-Inf(x[1]) 
      << "," << Sup(it->BoxRE)-Inf(it->BoxRE)
      << "], 'FaceColor','b')" << endl;
  } cout << endl;
}


double
MRSampler::getIU ()
{
  return _double (Sup (Integral));
}

double
MRSampler::getIL ()
{
  return _double (Inf (Integral));
}

double
MRSampler::getIUminusL ()
{
  return _double (diam (Integral));
}

double
MRSampler::getUmax ()
{
  return _double (Umax);
}

double
MRSampler::getPAest ()
{
  return _double (wsum / _double (nsum));
}

double
MRSampler::get_wsum ()
{
  return _double (wsum);
}

int
MRSampler::get_nsum ()
{
  return nsum;
}

double
MRSampler::get_wmax ()
{
  return _double (wmax);
}

double
MRSampler::get_wmin ()
{
  return _double (wmin);
}

//! update Integral using present partition
void
MRSampler::updateIntegral ()
{
  RangedLabBox theBox;
  interval localIntegral (0.0, 0.0);
  RangedLabBoxSet::const_iterator it = RangeDomainSet.begin ();
  for (; it != RangeDomainSet.end (); it++)
  {
    theBox = *it;
    //    localIntegral += theBox.BoxIntegral (UsingLogDensity);
    localIntegral += theBox.BoxIntegral (UsingLogDensity, f_scale);
  }
  cout << "in updateIntegral. IL, IU: " << Inf (Integral) << " " <<
    Sup (Integral) << endl;
  Integral = localIntegral;
}

//! Updating the upper and lower bounds for the global maximum of target F
void
MRSampler::updateUmax ()
{
  //! For each box, first get upper bound U = Sup(range enclosure), 
  //! and then Umax is the max over boxes of this
  //! and already have found during bisection Lmax, 
  //! i.e. max of over boxes of Inf(range enclosure)
  //! Lmax, Umax are rigorous lower and upper bounds for global maximum of f
  //! each time U is eval'd for a box and found to be > current Umax,
  //! eval f at midpoint of box. Then max of these is fmid_max
  RangedLabBox theBox;
  real fmid_max = BIGNEGATIVE;
  real f_scale_local;
  bool f_scaleDone_local = false;

  RangedLabBoxSet::const_iterator it = RangeDomainSet.begin ();
  bool first = true;
  for (; it != RangeDomainSet.end (); it++)
  {
    theBox = *it;
    real U = Sup (theBox.BoxRE);
    if (U > Umax || first)
    {
      Umax = U;
      cout << "Umax: " << Umax << endl;
      LabPnt lab_midpnt;
      lab_midpnt.Pnt = mid (theBox.LBox.Box);
      lab_midpnt.L = theBox.LBox.L;
      real fmid = F (lab_midpnt);
      if (fmid > fmid_max)
                    // f eval'd in middle of box with greatest U.
        fmid_max = fmid;
      first = false;
    }
  }
  //   Umax rigorous upper bound on f
                    // fMaxLB is lower bound on maximum of f
  fMaxLB = (Lmax > fmid_max) ? Lmax : fmid_max;

  if (UsingLogDensity)
  {
    // could use F - scale

    f_scale_local = Umax - ln (UmaxMAX);
    cout << "UmaxMAX, Umax, f_scale_local: " << UmaxMAX << " " 
         << Umax << " " << f_scale_local << endl;
    if (Umax - fMaxLB < LOGDIAMFMAX)
    {
      f_scaleDone_local = true;
    }
  }
  else
  {
    // could use F/f_scale
    f_scale_local = Umax / UmaxMAX;
    if (Umax < exp (LOGDIAMFMAX) * fMaxLB)
    {
      f_scaleDone_local = true;
    }
  }

  f_scale = (UseFScale)? f_scale_local: (UsingLogDensity)? 0: 1.0;
  cout << "f_scale: " << f_scale << "  " << f_scale_local << endl;
  f_scaleDone = f_scaleDone_local;
  cout << "bottom of updateUmax \n";
}

// private member functions
RangedLabBox
MRSampler::getBoxREInfo (LabBox LBox)
{ // return the RangedLabBox with information about Labelled Box LBox
  RangedLabBox result;
  ivector Box = LBox.Box;

  //  cout << "in getBoxREInfo. before F(LBox) " << endl;
  interval RangeEnclosure = F (LBox);
  //  cout << "in getBoxREInfo. after F(LBox). RE: " << RangeEnclosure << endl;
  //cout << "RangeEnclosure = " << RangeEnclosure << '\n'; getchar();

  // this is the RE of the function supplied, could be either pi or log(pi)
  result.BoxRE = RangeEnclosure;
  result.LBox = LBox;
  result.BoxVol = F.LabBoxVolume(LBox);
  result.BoxPriorIntegral = F.LabBoxPriorIntegral(LBox);

  return result;
}

// initialize rng
void
MRSampler::InitRNG (unsigned seed)
{                   
  const gsl_rng_type *Tgsl;
  gsl_rng_env_setup ();
  Tgsl = gsl_rng_default;
  rgsl = gsl_rng_alloc (Tgsl);
  gsl_rng_set (rgsl, seed);
  /* Note Regarding Quasi RNG: 
    For now we will just use one quasi rng, with dimension equal to dmax, 
    the max. dimension of all the labeled subdomains. Then when getting a 
    point in a d dimension box (d<dmax), just use the first d elements, 
    ignoring the rest. Not really sure whether this is legit.
  */
  qrgsl = gsl_qrng_alloc (gsl_qrng_sobol, n_dim_max);
  own_rng = true;
}

void
MRSampler::FirstBox ()
{ // setup the initial box (whole domain)
  RangedLabBox theBox;
  for (unsigned k = 0; k < topologies; k++)
  {
    cout << "in FirstBox, before getBoxREInfo. k: " << k << endl;
    theBox = getBoxREInfo (Domain[k]);
    theBox.Print(cout);
    cout << "in FirstBox, after getBoxREInfo " << endl;
    theBox.LBox.Print(cout);
    if (!RangeDomainSet.insert (theBox).second)
    {
      cout << "in FirstBox, failed to insert theBox. " << endl;
    }
  }
  Lmax = Inf (theBox.BoxRE);
  updateUmax ();
  cout << "in FirstBox, after updateUmax \n";
  Integral = theBox.BoxIntegral (UsingLogDensity, f_scale);
  cout << "bottom of FirstBox. \n";
  // exit(1);
}

/*! Bisect the box at the beginning of the RangeDomainSet and push two 
  resulting boxes onto set. Always bisect along the first dimension with the 
  longest side.  RangeDomainSet is STL set -- sorted associative container.
*/
void
MRSampler::Bisect ()
{
  //cout << "bisecting\n"; getchar();
  RangedLabBox theBox (*(RangeDomainSet.begin ()));
  //   cout << "top two boxes: " << endl;
  //   theBox.Print(cout);
                    // remove from the set the box to be bisected.
  RangeDomainSet.erase (RangeDomainSet.begin ());

  //   theBox = (*(RangeDomainSet.begin ()));
  //   theBox.Print(cout);
  //  (RangeDomainSet.begin())->Print(cout);

  ivector B (theBox.LBox.Box);
  int BoxDimensions = Ub(B) - Lb(B) + 1;
                    // midpoint of box
  rvector c (mid (B));
                    // Subboxes of Y
  imatrix U (0, 1, 1, BoxDimensions);

  real IBisectedBoxIntegral =  
    Inf (theBox.BoxIntegral (UsingLogDensity, f_scale));
  
  real SBisectedBoxIntegral = 
    Sup (theBox.BoxIntegral (UsingLogDensity, f_scale));
  // index of longest dimension of box
  int k = MaxDiamComp (B);
  real cc = (0)? c[k]: F.BisectPt(theBox.LBox, k);
  LabBox LBox;
  LBox.L = theBox.LBox.L;
  U[0] = U[1] = B;
  SetSup (U[0][k], cc);
  SetInf (U[1][k], cc);
  //  cout << "the 2 new boxes " << endl;
  for (int j = 0; j < 2; j++)
  {
    LBox.Box = U[j];
    //cout <<"getting rangeinfo for\n" << LBox.L << '\n' << LBox.Box; getchar();
    theBox = getBoxREInfo (LBox);
    //  cout << "j: " << j << endl; theBox.Print(cout);
    if (!RangeDomainSet.insert (theBox).second)
    {
      cout << "in Bisect. failed to insert theBox into RangeDomainSet " <<
        endl;
    }
    real newBoxL = Inf (theBox.BoxRE);
    //cout << "newBoxL = " << newBoxL; getchar();
    if (newBoxL > Lmax)
    {
      Lmax = newBoxL;
    }               // update Lmax
    Integral += theBox.BoxIntegral (UsingLogDensity, f_scale);
  }

  if (RangeDomainSet.size () % 100000 == 0)
    cout << "# boxes: " << get_n_boxes () << "  Acceptance Prob. >= " <<
      getPALB () << "\r";
  // cout << "end of bisect\n"; //getchar();

  // RangedLabBoxSet::const_iterator rdsit = RangeDomainSet.begin ();
  //  for( ; rdsit!= RangeDomainSet.end ();rdsit++){
  //      rdsit->Print(cout);
  //  }
}

// adaptively partition by bisection along longest dimension of box with 
// largest (vol*diam(range enclosure))
void
MRSampler::AdaptPartition (double Alb)
{                  
  unsigned int next_updateUmax = 32;
  //  cout << "Integral: " << Integral << endl;
  while ((RangeDomainSet.size () < Max_n_boxes) && (getPALB () < Alb))
  {
    Bisect ();
    //cout << "Integral: " << Integral << endl;
    //cout << "in AdaptPartition. "  << RangeDomainSet.size() << endl;
    if (!f_scaleDone && RangeDomainSet.size () == next_updateUmax)
    {
      cout << "in AdaptPartition before updateUmax. \n";
      updateUmax ();
      cout << "in AdaptPartition after updateUmax. \n";
      updateIntegral ();
      next_updateUmax *= 2;
    }
  }

  updateUmax ();
  cout << "in AdaptPartition after updateUmax2 \n";
  updateIntegral ();
  //   cout << "Integral: " << Integral << endl;
  cout << "# Adaptive partitioning complete. Boxes: " << RangeDomainSet.
    size () << "  Lower bound on Acceptance Prob.: " << getPALB () <<
    " IL, IU: " << getIL () << "   " << getIU () << endl;
  //  getchar();
}

void
MRSampler::SetupPDF ()
{
  //cout << "inside SetupPDF \n";
  unsigned int NBoxes = static_cast<unsigned int>(RangeDomainSet.size ());
  rvector VolBox(0, NBoxes-1);
  rvector PriorIntegralBox(0, NBoxes-1);
  Resize (UBox, 0, NBoxes - 1);
  Resize (LoBox, 0, NBoxes - 1);

  if (proposalpmf != NULL)
    free (proposalpmf);
  proposalpmf = (double *) malloc (NBoxes * sizeof (double));
  if (proposalpdf != NULL)
    free (proposalpdf);
  proposalpdf = (double *) malloc (NBoxes * sizeof (double));

  RangedLabBox theBox;
  Umax = Sup (RangeDomainSet.begin ()->BoxRE);
  RangedLabBoxSet::const_iterator it = RangeDomainSet.begin ();
  for (unsigned int ui = 0; it != RangeDomainSet.end (); it++, ui++)
  {                 // pull em out from top of pq
    theBox = *it;
    DomainParts.push_back (theBox.LBox);
    UBox[ui] = Sup (theBox.BoxRE);
    LoBox[ui] = Inf (theBox.BoxRE);
    VolBox[ui] = theBox.BoxVol;
    PriorIntegralBox[ui] = theBox.BoxPriorIntegral;
    if (UBox[ui] > Umax)
      Umax = UBox[ui];
  }

                    //getchar();
  cout << "#Using log(pi)? " << UsingLogDensity << endl;
  double pTotal = 0.0;

                    // reset proposalpmf
  for (unsigned int ui = 0; ui < NBoxes; ui++)
  {
    real BoxI = (0)? VolBox[ui]: PriorIntegralBox[ui];
    if (UsingLogDensity)
    {
      //  UBox[ui] = exp (UBox[ui] - f_scale); 
      // leave Ubox as just Sup RE of logdensity
      //  LoBox[ui] = exp (LoBox[ui] - f_scale);
      // could we use BoxIntegral here?
      if(WEIGHTED_SQUEEZE)
      {
        proposalpmf[ui] = _double((exp (UBox[ui] - f_scale) - 
                                   exp (LoBox[ui] - f_scale)) * BoxI);
      }
      else
      {
        proposalpmf[ui] = _double(exp(UBox[ui] - f_scale) * BoxI);
      }
    }
    else
    {
                    //
      if (WEIGHTED_SQUEEZE)
      {
        proposalpmf[ui] = _double ((UBox[ui] - LoBox[ui]) * BoxI);
      }
      else
      {
        proposalpmf[ui] = _double (UBox[ui] * BoxI);
      }
    }

                    //initialize
    DomainParts[ui].SamplesToDo = 0;
  }
  //make pmf a probability mass function -- normalize
  pTotal =
    std::accumulate (proposalpmf, proposalpmf + NBoxes, 0.0,
    kahan_sum < double >());

  //   UBIntegral = pTotal;
  std::transform (proposalpmf, proposalpmf + NBoxes, proposalpmf,
    bind2nd (divides < double >(), pTotal));

  //Now we will use proposalpdf to keep track of 
  // the DENSITY (normalized simple fnctn) of a proposed point
  int count0probBoxes16 = 0, count0probBoxes10 = 0, count0probBoxes6 =
    0, count0probBoxes3 = 0;
  for (unsigned int ui = 0; ui < NBoxes; ui++)
  {
    if(0)
    {
                    // (BoxI); // ????
      proposalpdf[ui] = proposalpmf[ui] / _double (VolBox[ui]);
    }
    else
    {
      proposalpdf[ui] = proposalpmf[ui] / _double (PriorIntegralBox[ui]);
    }
    // do we need proposalpdf??? 
    // it will not be const over box, if prior is not const.
    if (proposalpmf[ui] >= 1e-3)
      count0probBoxes3++;
    if (proposalpmf[ui] >= 1e-6)
      count0probBoxes6++;
    if (proposalpmf[ui] <= 1e-16)
      count0probBoxes16++;
    if (proposalpmf[ui] <= 1e-10)
      count0probBoxes10++;
  }
  // sometimes 20% of the boxes have mass <= 1e-16 !!!  this unnecessarily 
  // enlarges our priority queue and the SetupWalker
  // need regular sub-paving to make this efficient 
  // perhaps have > 1 list of boxes
  cout << "#No. of Boxes with proposal mass function <= 1e-16 " <<
    count0probBoxes16 << endl;
  cout << "#No. of Boxes with proposal mass function <= 1e-10 " <<
    count0probBoxes10 << endl;
  cout << "#No. of Boxes with proposal mass function >= 1e-6 " <<
    count0probBoxes6 << endl;
  cout << "#No. of Boxes with proposal mass function >= 1e-3 " <<
    count0probBoxes3 << endl;
}

void
MRSampler::SetupWalker (unsigned int samples)
{                   // Setup residual_proposalpmf
  unsigned int NBoxes = static_cast<unsigned int>(RangeDomainSet.size ());

  if (residual_proposalpmf != NULL)
                    // start over again
    free (residual_proposalpmf);
  residual_proposalpmf = (double *) malloc (NBoxes * sizeof (double));

  nonresidual_samples = 0;
  // preprocess proposalpmf to get gslpdfstruct, used by gsl_ran_discrete
  if (gslpdfstruct != NULL)
    gsl_ran_discrete_free (gslpdfstruct);
  if (samples != 0)
  {                 // modify pmf for residual sampling to minimizing variance
    double blowup, intpart;
    double dsamples = double (samples);
    double count_residuals = 0;
    for (unsigned int ui = 0; ui < NBoxes; ui++)
    {
      blowup = proposalpmf[ui] * dsamples;
                    // save the residue in the revided residual pmf
      count_residuals += (residual_proposalpmf[ui] = modf (blowup, &intpart));
      nonresidual_samples += (DomainParts[ui].SamplesToDo =
        int (intpart));
    }
    std::transform (residual_proposalpmf, residual_proposalpmf + NBoxes,
      residual_proposalpmf, bind2nd (divides < double >(),
      count_residuals));
  }
                    // this is residual part now
  gslpdfstruct = gsl_ran_discrete_preproc ((size_t) NBoxes, 
                                           residual_proposalpmf);
}

//Output the RangeDomain set
std::ostream & operator<<(std::ostream &os, const MRSampler& mrs)
{
  // uses public member function MRSoutput to generate output
  mrs.MRSoutput(os);

  return os;
}

/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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

/*! \file
\brief Definitions for classes measuring spsnodes.
*/

#include "spsnode_measure_obj.hpp"


namespace subpavings {

    

    /* A virtual class providing a way to measure spsnodes.*/
	SPSNodeMeasure::~SPSNodeMeasure(){}

    
    /* Class measuring count of data points associated with a node.    */
    cxsc::real SPSNodeMeasureCount::
		operator() (const SPSnode * const spn) const
	{ return (1.0*spn->getCounter()); }
    
    /* Class measuring volume of box of node.    */
    cxsc::real SPSNodeMeasureVol::operator()
		(const SPSnode * const spn) const
    { return (spn->nodeVolume()); }
    
    /*Class measuring count/volume (ie unormalised histogram 
	 * 			height) of box of node.*/
    cxsc::real SPSNodeMeasureHeight::operator()
		(const SPSnode * const spn) const
    { return (spn->getCountOverVolume()); }
   
    /* Class measuring change in EMP under COPERR from splitting 2 nodes.

    Under COPERR, EMP is -1/n^2 x sum over leaves of
    (counts in leaf squared / volume of leaf)
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing  change in the sum over leaves of
    (counts in leaf squared over volume of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPCOPERR(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering.
	* 
	So we return the negated change.
    */
	
	SPSNodeMeasureEMPSumChangeCOPERR::
				SPSNodeMeasureEMPSumChangeCOPERR(size_t n)
				 : bigN(n) {}
    
	
	cxsc::real SPSNodeMeasureEMPSumChangeCOPERR::
		SPSNodeMeasureEMPSumChangeCOPERR::operator()
					(const SPSnode * const spn) const
	{
		// note negated
	   return (-rnd(spn->getSplitChangeEMPCOPERR(bigN)));
	}



    /* Class measuring change in EMP under AIC from splitting 2 nodes.

    Under AIC, EMP is -1 x sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)))
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing the change in -1 x the sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPAIC(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering.
	
	So we return the negated value.    */
    cxsc::real SPSNodeMeasureEMPSumChangeAIC::
		operator() (const SPSnode * const spn) const
	{
		//note negated
		return (-rnd(spn->getSplitChangeEMPAIC()));
	}
   

    /* Class measuring change in EMP under COPERR from merging 2 nodes.

    Under COPERR, EMP is -1/n^2 x sum over leaves of
    (counts in leaf squared / volume of leaf)
    where n is the total number of data points in the histogram

    For two subleaf nodes we are comparing  change in the sum over leaves of
    (counts in leaf squared over volume of leaf)
    which would result if each node were to be the one to be merged.

    Merges take from the left of the queue first ("smallest")

    The smaller (more negative) the value returned by getMergeChangeEMPCOPERR(),
    the more a node will reduce or least increase the overall EMP by being
    merged, so it should be lower, ie more to left, in the ordering.
    */
	SPSNodeMeasureEMPSumChangeMergeCOPERR::
				SPSNodeMeasureEMPSumChangeMergeCOPERR(size_t n)
				 : bigN(n) {}
				
				
    cxsc::real SPSNodeMeasureEMPSumChangeMergeCOPERR::
		operator() (const SPSnode * const spn) const
	{	
		return (rnd(spn->getMergeChangeEMPCOPERR(bigN)));
	}


    /* Class measuring change in EMP under AIC from merging 2 nodes.

    Under AIC, EMP is -1 x sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)))
    where n is the total number of data points in the histogram

    For two subleaf nodes we are comparing the change in -1 x the sum over leaves
    of (counts in leaf x (ln(count in leaf /(n x vol of leaf)
    which would result if each node were to be the one to be merged.

    Merges take from the left of the queue first ("smallest")

    The smaller (more negative) the value returned by getMergeChangeEMPAIC(),
    the more a node will reduce or least increase the overall EMP by being
    merged, so it should be lower, ie more to left, in the ordering.
    */
    cxsc::real SPSNodeMeasureEMPSumChangeMergeAIC::
		operator() (const SPSnode * const spn) const
	{
		return (rnd(spn->getMergeChangeEMPAIC()));
	}
   
	
	 //From Gloria, adapted
	/* Class measuring volume * (1- empirical mass). */
    SPSNodeMeasureVolMassMinus::
				SPSNodeMeasureVolMassMinus(size_t n) : bigN(n) {}
		
	cxsc::real SPSNodeMeasureVolMassMinus::
		operator() (const SPSnode * const spn) const
	{   
		cxsc::real empMass(1.0 - (spn->getCounter()/(bigN*1.0)));
		
		cxsc::real massVol = empMass * spn->nodeRealVolume();
		
		return massVol;
	}
	
	 /*biggest absolute difference in 
	 * prospective child counts if split
	 * RC = C - LC, so |RC - LC| = |C - 2LC| */
	cxsc::real SPSNodeMeasureChildDifferential::
		operator() (const SPSnode * const spn) const
	{   
		size_t twiceLeftCount = 2*(spn->getLeftCountIfSplit());
		size_t thisCounter = spn->getCounter();
		
		return (twiceLeftCount > thisCounter ? 1.0*(twiceLeftCount - thisCounter)
							: 1.0*(thisCounter - twiceLeftCount) );
	}
	
	
     /* Class giving no measure.    */
    cxsc::real SPSNodeMeasureNothing::
		operator() (const SPSnode * const) const
	{ return cxsc::real(1.0); }
}



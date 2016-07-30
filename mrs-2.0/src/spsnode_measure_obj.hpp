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
\brief Interface and concrete classes for ways of measuring spsnodes.
*/


#ifndef __SPSNODE_MEASURE_OBJ_HPP__
#define __SPSNODE_MEASURE_OBJ_HPP__

#include "spsnode.hpp"
#include "real.hpp"


namespace subpavings {
	
	/*! \brief A virtual class providing a way to measure spsnodes.

    The measure should be calculated so that if a priority splitting queue
	was ordered by node measure, the nodes with the largest measures should be 
	split first.
    */
    class SPSNodeMeasure {

        public:
		
		virtual ~SPSNodeMeasure();

        /*! \brief Return the measure of the node. */
        virtual cxsc::real operator() (const SPSnode * const spn) const = 0;
    };
	
	
    /*! \brief Class measuring count of data points associated with a node.
    */
    class SPSNodeMeasureCount : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
    };

    /*! \brief Class comparing on volume of box of node.
    */
    class SPSNodeMeasureVol : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
    };
    
    /*! \brief Class measuring count/volume (ie histogram height) of box of node.
    */
    class SPSNodeMeasureHeight : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
    };

    
	/*! \brief Class measuring change in EMP under COPERR from splitting 2 nodes.

    Under COPERR, EMP is -1/n^2 x sum over leaves of
    (counts in leaf squared / volume of leaf)
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing  change in the sum over leaves of
    (counts in leaf squared over volume of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPCOPERR(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering, 
	so we meausure using the negated value.     */
    class SPSNodeMeasureEMPSumChangeCOPERR : public SPSNodeMeasure
    {
		public:
			explicit SPSNodeMeasureEMPSumChangeCOPERR(size_t n);
			
			cxsc::real operator() (const SPSnode * const spn) const;
        		
		private:
			SPSNodeMeasureEMPSumChangeCOPERR();
			size_t bigN;
    };


    /*! \brief Class measuring change in EMP under AIC from splitting 2 nodes.

    Under AIC, EMP is -1 x sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)))
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing the change in -1 x the sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPAIC(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering, 
	so we meausure using the negated value.     */
    class SPSNodeMeasureEMPSumChangeAIC : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
    };

    /*! \brief Class measuring change in EMP under COPERR from merging 2 nodes.

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
    class SPSNodeMeasureEMPSumChangeMergeCOPERR : public SPSNodeMeasure
   {
		public:
			explicit SPSNodeMeasureEMPSumChangeMergeCOPERR(size_t n);
			
			cxsc::real operator() (const SPSnode * const spn) const;
        		
		private:
			SPSNodeMeasureEMPSumChangeMergeCOPERR();
			size_t bigN;
    };


    /*! \brief Class measuring change in EMP under AIC from merging 2 nodes.

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
    class SPSNodeMeasureEMPSumChangeMergeAIC : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
        
    };
	
	 //From Gloria, adapted
	/*! \brief Class comparing volume * (1- empirical mass).
	*/
    class SPSNodeMeasureVolMassMinus: public SPSNodeMeasure
    {
		public:
			explicit SPSNodeMeasureVolMassMinus(size_t n);
		
			cxsc::real operator() (const SPSnode * const spn) const;
			
		
		private:
			SPSNodeMeasureVolMassMinus(); // no default constructor
			const size_t bigN;
	    
	};
	
	/*! \brief Class measuring biggest absolute difference in 
	 * prospective child counts if split.
	 * 
	This is proportional to the L1 error between histogram estimate
	\a spn and histogram estimate child nodes of \a spn.  
	
	L1 error is 1/N * { | n- 2n_l | + |n - 2n_r | }
	where N is total data points in the entire histogram,
	n is data points in the node, n_l is data points going to left
	child and n_r is data points going to right child.  
	
	n_r = n - n_l so n = 2n_r = n - 2(n - n_l) = 2n_l - n
	So
	1/N * { | n- 2n_l | + |n - 2n_r | } = 2/N * | n- 2n_l |
	and n - 2n_l = n- n_l - n_l = n_r - n_l is the difference between
	prospective child counts if split, so this measure can be
	treated as proportional to the L1 error between the histogram estimate
	on \a spn and the histogram estimate on the children of \a spn. */
   class SPSNodeMeasureChildDifferential : public SPSNodeMeasure
    {
        cxsc::real operator() (const SPSnode * const spn) const;
    };
	
    /*! \brief Class measuring nodes to give no change in ordering.

    */
    class SPSNodeMeasureNothing : public SPSNodeMeasure
    {
		cxsc::real operator() (const SPSnode * const) const;
		
    };

}



#endif



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
\brief Classes for comparing spsnodes.
*/

#ifndef ___SPSNODECOMP_HPP__
#define ___SPSNODECOMP_HPP__

namespace subpavings {

    //! Forward class declarations
    class SPSnode;


    /*! \brief A Virtual class providing a way to compare spsnodes.

    These classes create an ordering with the 'largest' at the right, or 'end'.
    This suits the implementation of priority queues for the HistogramWrapper,
    which takes nodes from the end of a multiset.
    */
    class NodeCompObj {

        public:

        /*! return true if lhs is 'smaller' (lower in priority) than rhs. */
        virtual bool operator() (const SPSnode * const lhs,
                                    const SPSnode * const rhs) const = 0;
    };


    /** @name Concrete classes derived from NodeCompObj

    These classes provide ways of comparing spsnodes.

    */
    //@{



    /*! \brief Class comparing on count of data points associated with a node.
    */
    class CompCount : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        { return (lhs->getCounter() < rhs->getCounter()); }
    };

    /*! \brief Class comparing on volume of box of node.
    */
    class CompVol : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        { return (lhs->nodeVolume() < rhs->nodeVolume()); }
    };
    
    /*! \brief Class comparing on count/volume (ie histogram height) of box of node.
    */
    class CompHeight : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        { return (lhs->getCountOverVolume() 
					< rhs->getCountOverVolume()); }
    };

    /*! \brief Class comparing change in EMP under COPERR from splitting 2 nodes.

    Under COPERR, EMP is -1/n^2 x sum over leaves of
    (counts in leaf squared / volume of leaf)
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing  change in the sum over leaves of
    (counts in leaf squared over volume of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPCOPERR(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering.
    */
    class CompEMPSumChangeCOPERR : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        {
            size_t nLhs = lhs->getRootCounter();
            size_t nRhs = rhs->getRootCounter();

            return (rnd(lhs->getSplitChangeEMPCOPERR(nLhs)) >
                            rnd(rhs->getSplitChangeEMPCOPERR(nRhs)));
        }
    };


    /*! \brief Class comparing change in EMP under AIC from splitting 2 nodes.

    Under AIC, EMP is -1 x sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)))
    where n is the total number of data points in the histogram

    For two leaf nodes we are comparing the change in -1 x the sum over leaves of
    (counts in leaf x (ln(count in leaf /(n x vol of leaf)
    which would result if each node were to be the one to be split.

    The smaller (more negative) the value returned by getSplitChangeEMPAIC(),
    the more a node will reduce or least increase the overall EMP by being
    split, so it should be higher, ie more to right, in the ordering.
    */
    class CompEMPSumChangeAIC : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        {
            return (rnd(lhs->getSplitChangeEMPAIC()) >
                        rnd(rhs->getSplitChangeEMPAIC()));
        }
    };

    /*! \brief Class comparing change in EMP under COPERR from merging 2 nodes.

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
    class CompEMPSumChangeMergeCOPERR : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        {
            size_t nLhs = lhs->getRootCounter();
            size_t nRhs = rhs->getRootCounter();

            return (rnd(lhs->getMergeChangeEMPCOPERR(nLhs)) <
                            rnd(rhs->getMergeChangeEMPCOPERR(nRhs)));
        }
    };


    /*! \brief Class comparing change in EMP under AIC from merging 2 nodes.

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
    class CompEMPSumChangeMergeAIC : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        {
            return (rnd(lhs->getMergeChangeEMPAIC()) <
                        rnd(rhs->getMergeChangeEMPAIC()));
        }
    };
	
	 //From Gloria, adapted
	/*! \brief Class comparing volume * (1- empirical mass).
	*/
    class CompVolMassMinus: public NodeCompObj
    {
		public:
			CompVolMassMinus(size_t bigN) : n(bigN) {};
		
			bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
			{   
				cxsc::real lEmpMass(1.0 - (lhs->getCounter()/(n*1.0)));
				cxsc::real rEmpMass (1.0 - (rhs->getCounter()/(n*1.0)));
				
				cxsc::real lMassVol = lEmpMass * lhs->nodeRealVolume();
				cxsc::real rMassVol = rEmpMass * rhs->nodeRealVolume();

				return ( lMassVol <rMassVol );
			}
		
		private:
			CompVolMassMinus(); // no default constructor
			const size_t n;
	    
	};
	class CompVolMassMinusAdj: public NodeCompObj
    {
		public:
			CompVolMassMinusAdj(size_t bigN) : n(bigN) {};
		
			bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
			{   
				size_t lCount = lhs->getCounter();
				size_t rCount = rhs->getCounter();
				
				if ((lCount == rCount) && 
					(lhs->getNodeDepth() == rhs->getNodeDepth())) {
						return (lhs->getNodeName() < rhs->getNodeName()) ;
				}
								
				cxsc::real lEmpMass(1.0 - (lCount/(n*1.0)));
				cxsc::real rEmpMass (1.0 - (rCount/(n*1.0)));
				
				cxsc::real lMassVol = lEmpMass * lhs->nodeRealVolume();
				cxsc::real rMassVol = rEmpMass * rhs->nodeRealVolume();

				return ( lMassVol <rMassVol );
			}
		
		private:
			CompVolMassMinusAdj(); // no default constructor
			const size_t n;
	    
	};

    /*! \brief Class comparing nodes to give no change in ordering.

    */
    class CompNothing : public NodeCompObj
    {
		bool operator()   (const SPSnode * const /*lhs*/,
								const SPSnode * const /*rhs*/) const
		{
			return false;
		}
		
    };

    //@}
}

#endif



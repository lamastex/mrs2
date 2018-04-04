/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
* 
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

/*! \file      nodecompobj.hpp
\brief Classes for comparing spsnodes.
*/

#ifndef ___SPSNODECOMP_HPP__
#define ___SPSNODECOMP_HPP__

#include "toolz.hpp" // to use MaxDiam

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
            size_t nRhs = lhs->getRootCounter();

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
            size_t nRhs = lhs->getRootCounter();

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

    /*! \brief Class comparing nodes to give no change in ordering.

    */
    class CompNothing : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        {
            return false;
        }
    };
	 
	  /*! \brief Class comparing on count/volume (ie histogram height) of box of node.
    */
    class CompHeight : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
        { return (lhs->getCounter()/lhs->nodeVolume() 
					< rhs->getCounter()/rhs->nodeVolume()); }
    };
    
    //gat41
   /*! \brief Class comparing the area of a histogram.
   */
    class CompArea : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			  size_t n = lhs->getRootCounter();
			  
			  cxsc::interval lCount = interval(lhs->getCounter()*1.0/(n*1.0));
			  cxsc::interval rCount = interval(rhs->getCounter()*1.0/(n*1.0));
			  
			  cxsc::interval lVol = interval(lhs->nodeVolume());
			  cxsc::interval rVol = interval(rhs->nodeVolume());
			  
			  cxsc::interval lMassVol = lCount * lVol;
			  cxsc::interval rMassVol = rCount * rVol;
			 
			  cxsc::real lMid = mid(lMassVol);
			  cxsc::real rMid = mid(rMassVol);
			  
			  /*
			  std::cout << n << std::endl;
			  std::cout << lhs->getNodeName() << "\t" << rhs->getNodeName() << std::endl;
			  std::cout << lCount << "\t" << rCount << std::endl;
			  std::cout << lVol << "\t" << rVol << std::endl;
			  std::cout << lMassVol << "\t" << rMassVol << std::endl;
			  std::cout << lMid << "\t" << rMid << std::endl;
			  */
			  
			  return ( (lMid < rMid) );
			  //return ( (lhs->getCounter()*1.0) * lhs->nodeVolume() 
				//			< (rhs->getCounter()*1.0) * rhs->nodeVolume()); 
		}
	};
	 
	//gat41
   /*! \brief Class comparing the "extended"area of a histogram.
   */
    class CompExtArea : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			  size_t n = lhs->getRootCounter();
			  
			  cxsc::interval lCount = interval(lhs->getCounter()*1.0/(n*1.0));
			  cxsc::interval rCount = interval(rhs->getCounter()*1.0/(n*1.0));
			  
			  cxsc::interval lVol = interval(lhs->nodeVolume());
			  cxsc::interval rVol = interval(rhs->nodeVolume());
			  
			  cxsc::interval lMassVol = lCount * lVol *lVol;
			  cxsc::interval rMassVol = rCount * rVol * rVol;
			 
			  cxsc::real lMid = mid(lMassVol);
			  cxsc::real rMid = mid(rMassVol);

			  return ( (lMid < rMid) );
			  //return ( (lhs->getCounter()*1.0) * lhs->nodeVolume() 
				//			< (rhs->getCounter()*1.0) * rhs->nodeVolume()); 
		}
	};
	 
		//gat41
   /*! \brief Class comparing the Chebyshev distance between the mean and 
    * uniform mean.
   */
    class CompMean : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			//std::cout << lhs->getNodeName() << "\t" << lhs->getChebDistMean() << std::endl;
			//std::cout << rhs->getNodeName() << "\t" << rhs->getChebDistMean() << std::endl;
			//std::cout << "=================" << std::endl;
			
			return ( (lhs->getChebDistMean()) < 
						(rhs->getChebDistMean()) );
		}
	}; 
	 
	 	//gat41
   /*! \brief Class comparing the Chebyshev distance between the mean and 
    * uniform mean multiplied with the emprical mass.
   */
    class CompMeanMass : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			//std::cout << lhs->getNodeName() << "\t" << lhs->getChebDistMean() << std::endl;
			//std::cout << rhs->getNodeName() << "\t" << rhs->getChebDistMean() << std::endl;
			//std::cout << "=================" << std::endl;

			return ( (lhs->getChebDistMean()*lhs->getEmpMass()) < 
						(rhs->getChebDistMean()*rhs->getEmpMass()) );
		}
	}; 
	 
	 	  	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompCovar : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			//std::cout << lhs->getNodeName() << "\t" << lhs->getChebDistMean() << std::endl;
			//std::cout << rhs->getNodeName() << "\t" << rhs->getChebDistMean() << std::endl;
			//std::cout << "=================" << std::endl;

			return ( (lhs->getChebDistCovar()) < 
						(rhs->getChebDistCovar()) );
		}
	}; 
	 
	  	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompCovarMass : public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			//std::cout << lhs->getNodeName() << "\t" << lhs->getChebDistMean() << std::endl;
			//std::cout << rhs->getNodeName() << "\t" << rhs->getChebDistMean() << std::endl;
			//std::cout << "=================" << std::endl;

			return ( (lhs->getChebDistCovar()*lhs->getEmpMass()) < 
						(rhs->getChebDistCovar()*rhs->getEmpMass()) );
		}
	}; 
	 
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1D: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			/*std::cout << lhs->getNodeName() << "\t" << lhs->getHellingerDist1D() << std::endl;
			std::cout << rhs->getNodeName() << "\t" << rhs->getHellingerDist1D() << std::endl;
			std::cout << "=================" << std::endl;
			*/
			return ( (lhs->getHellingerDist1D()) < 
						(rhs->getHellingerDist1D()) );
		}
	}; 
	 
	 //gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1DMass: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			/*
			std::cout << lhs->getNodeName() << "\t" << lhs->getHellingerDist1D()*lhs->getEmpMass() << std::endl;
			std::cout << rhs->getNodeName() << "\t" << rhs->getHellingerDist1D()*rhs->getEmpMass() << std::endl;
			std::cout << "=================" << std::endl;
			*/
			return ( (lhs->getHellingerDist1D()*lhs->getEmpMass()) < 
						(rhs->getHellingerDist1D()*rhs->getEmpMass()) );
		}
	}; 
	
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1DMassDiam: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			/*
			std::cout << lhs->getNodeName() << "\t" << lhs->getHellingerDist1D()*lhs->getEmpMass() << std::endl;
			std::cout << rhs->getNodeName() << "\t" << rhs->getHellingerDist1D()*rhs->getEmpMass() << std::endl;
			std::cout << "=================" << std::endl;
			*/
			 int maxdiamcomp = 0;  // to take value calculated from MaxDiam
          // find the maximum diameter, put the max dimension into maxdiamcomp
          double maxDiamL = ::MaxDiam(lhs->getBox(), maxdiamcomp);
          double maxDiamR = ::MaxDiam(rhs->getBox(), maxdiamcomp);
            
			return ( (lhs->getHellingerDist1D()*lhs->getEmpMass()*maxDiamL) < 
						(rhs->getHellingerDist1D()*rhs->getEmpMass()*maxDiamR) );
		}
	}; 
	
		 //gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDistMass: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			/*
			std::cout << lhs->getNodeName() << "\t" << lhs->getHellingerDist1D()*lhs->getEmpMass() << std::endl;
			std::cout << rhs->getNodeName() << "\t" << rhs->getHellingerDist1D()*rhs->getEmpMass() << std::endl;
			std::cout << "=================" << std::endl;
			*/
			return ( (lhs->getHellingerDist()*lhs->getEmpMass()) < 
						(rhs->getHellingerDist()*rhs->getEmpMass()) );
		}
	}; 
	
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDistMassDiam: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      { 
			/*
			std::cout << lhs->getNodeName() << "\t" << lhs->getHellingerDist1D()*lhs->getEmpMass() << std::endl;
			std::cout << rhs->getNodeName() << "\t" << rhs->getHellingerDist1D()*rhs->getEmpMass() << std::endl;
			std::cout << "=================" << std::endl;
			*/
			 int maxdiamcomp = 0;  // to take value calculated from MaxDiam
          // find the maximum diameter, put the max dimension into maxdiamcomp
          double maxDiamL = ::MaxDiam(lhs->getBox(), maxdiamcomp);
          double maxDiamR = ::MaxDiam(rhs->getBox(), maxdiamcomp);
            
			return ( (lhs->getHellingerDist()*lhs->getEmpMass()*maxDiamL) < 
						(rhs->getHellingerDist()*rhs->getEmpMass()*maxDiamR) );
		}
	}; 
	
	//gat41
   /*! \brief Class comparing volume multiplied with the inverse of the empirical mass.
   */
    class CompVolInv: public NodeCompObj
    {
        bool operator()   (const SPSnode * const lhs,
                            const SPSnode * const rhs) const
      {   
			  cxsc::interval lCount = interval(1.0/(lhs->getEmpMass()));
			  cxsc::interval rCount = interval(1.0/(rhs->getEmpMass()));
			  
			  cxsc::interval lVol = interval(lhs->nodeVolume());
			  cxsc::interval rVol = interval(rhs->nodeVolume());
			  
			  cxsc::interval lMassVol = lCount * lVol;
			  cxsc::interval rMassVol = rCount * rVol;
			 
			  cxsc::real lMid = mid(lMassVol);
			  cxsc::real rMid = mid(rMassVol);
		}
	 }; 
	
	
	
	
    //@}
}

#endif



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
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*! \file      nodecompobjval.hpp
\brief Classes for comparing SPSVnodes.
*/

#ifndef ___SPSVNODECOMPVAL_HPP__
#define ___SPSVNODECOMPVAL_HPP__

#include "toolz.hpp" //to use MaxDiam

namespace subpavings {

    //! Forward class declarations
    class SPSVnode;
    
    /*! \brief A Virtual class providing a way to compare SPSVnodes.

    These classes create an ordering with the 'largest' at the right, or 'end'.
    This suits the implementation of priority queues for the HistogramWrapper,
    which takes nodes from the end of a multiset.
    */
    class NodeCompObjVal {

        public:

        /*! return true if lhs is 'smaller' (lower in priority) than rhs. */
        virtual bool operator() (const SPSVnode * const lhs,
                                    const SPSVnode * const rhs) const = 0;
    };


    //@{



    /*! \brief Class comparing on count of data points associated with a node.
    */
    class CompCountVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
        { return (lhs->getCounter() < rhs->getCounter()); }
    };

    /*! \brief Class comparing on volume of box of node.
    */
    class CompVolVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
        { return (lhs->nodeVolume() < rhs->nodeVolume()); }
    };

    
    /*! \brief Class comparing nodes to give no change in ordering.

    */
    class CompNothingVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
        {
            return false;
        }
    };
	 
	 /*! \brief Class comparing on count/volume (ie histogram height) of box of node.
    */
    class CompHeightVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
        { return (lhs->getCounter()/lhs->nodeVolume() 
					< rhs->getCounter()/rhs->nodeVolume()); }
    }; 

  
    //gat41
   /*! \brief Class comparing the area of a histogram.
   */
    class CompAreaVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
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
			  
			  return ( (lMid < rMid) );
		}
	};
	 
	 //gat41
   /*! \brief Class comparing the Chebyshev distance between the mean and 
    * uniform mean.
   */
    class CompMeanVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 
			
			return ( (lhs->getChebDistMean()) < 
						(rhs->getChebDistMean()) );
		}
	}; 
	 
	 	//gat41
   /*! \brief Class comparing the Chebyshev distance between the mean and 
    * uniform mean multiplied with the emprical mass.
   */
    class CompMeanMassVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 

			return ( (lhs->getChebDistMean()*lhs->getEmpMass()) < 
						(rhs->getChebDistMean()*rhs->getEmpMass()) );
		}
	}; 
	 
	 	  	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompCovarVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 

			return ( (lhs->getChebDistCovar()) < 
						(rhs->getChebDistCovar()) );
		}
	}; 
	 
	  	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompCovarMassVal : public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 
			return ( (lhs->getChebDistCovar()*lhs->getEmpMass()) < 
						(rhs->getChebDistCovar()*rhs->getEmpMass()) );
		}
	}; 
	 
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1DVal: public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 
			return ( (lhs->getHellingerDist1D()) < 
						(rhs->getHellingerDist1D()) );
		}
	}; 
	 
	 //gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1DMassVal: public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 
			return ( (lhs->getHellingerDist1D()*lhs->getEmpMass()) < 
						(rhs->getHellingerDist1D()*rhs->getEmpMass()) );
		}
	}; 
	
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDist1DMassDiamVal: public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
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
    class CompHellingerDistMassVal: public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
      { 
			return ( (lhs->getHellingerDist()*lhs->getEmpMass()) < 
						(rhs->getHellingerDist()*rhs->getEmpMass()) );
		}
	}; 
	
	//gat41
   /*! \brief Class comparing the Chebyshev distance between the var-covar and 
    * uniform var-covar multiplied with the empirical mass.
   */
    class CompHellingerDistMassDiamVal: public NodeCompObjVal
    {
        bool operator()   (const SPSVnode * const lhs,
                            const SPSVnode * const rhs) const
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
	
    //@}
}

#endif



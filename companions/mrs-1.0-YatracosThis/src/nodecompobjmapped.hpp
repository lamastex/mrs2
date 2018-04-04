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

/*! \file      nodecompobjmapped.hpp
\brief Classes for comparing spnodes visited by mappedspnode.
*/

#ifndef ___NODECOMPOBJMAPPED_HPP__
#define ___NODECOMPOBJMAPPED_HPP__


using namespace cxsc;
using namespace std;


namespace subpavings {

    //! Forward class declarations
    class SPnode;
    
    /*! \brief A Virtual class providing a way to compare spnodes visited by mappedspnodes.

    These classes create an ordering with the 'largest' at the right, or 'end'.
    This suits the implementation of priority queues for the HistogramWrapper,
    which takes nodes from the end of a multiset.
    */
    class NodeCompObjMapped {

		public:

		/*! return true if lhs is 'smaller' (lower in priority) than rhs. */
		virtual bool operator() (const SPnode * const lhs,
									const SPnode * const rhs) const = 0;
	};


    //@{

   /*! \brief Class comparing the approximate area of a box by multiplying the 
    * box volume with the diameter of its range enclosure.
   */
    class CompSPArea: public NodeCompObjMapped
    {
		MappedFobj& fobj;

		CompSPArea(); // private default constructor

		public:
		CompSPArea(MappedFobj& f) : fobj(f) {}
		 
      bool operator()   (const SPnode * const lhs,
                            const SPnode * const rhs) const
      { 
			//std::cout << "----doing comparisons" << endl;
			//std::cout << lhs->getNodeName() << "\t" << lhs->getBox() << "\t";
			ivector lbox = lhs->getBox();
			interval lRange = fobj(lbox);
			interval lArea = (lhs->nodeVolume()) * (lRange);
			//riemann sum
			//cout << "left area: " << lArea << endl;
			
			//std::cout << rhs->getNodeName() << "\t" << rhs->getBox() << "\t";
			ivector rbox = rhs->getBox();
			interval rRange = fobj(rbox);
			interval rArea = (rhs->nodeVolume()) * (rRange);
			//cout << "right area: " << rArea << endl;
			
			//cout << (diam(lArea)) << "\t" << (diam(rArea)) << endl;
			//the diameters are the uncertainty in the approximate to the integral error
			return (diam(lArea) < diam(rArea));
		}
	}; 

    //@}
}

#endif

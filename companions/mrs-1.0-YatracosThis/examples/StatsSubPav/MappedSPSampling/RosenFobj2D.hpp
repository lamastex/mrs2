/*! \file RosenFobj2D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __RosenFOBJ2D_HPP__
#define __RosenFOBJ2D_HPP__


#include "mappedFobj2D.hpp"


class RosenFobj2D : public subpavings::MappedFobj2D {

    public:

        //! no argument constructor
        RosenFobj2D() {};

        //! declare function for interval image of a box
        
        //ivector
        
       // virtual cxsc::interval operator()(const cxsc::ivector& ival) const;
			
        
			virtual cxsc::interval operator()(const cxsc::interval& ival1,
                              const cxsc::interval& ival2) const;
			
  
        
        
			//! declare function for real image of reals
			//rvector
			//virtual cxsc::real operator()(const cxsc::rvector& r) const;
        
			virtual cxsc::real operator()(const cxsc::real& r1,
                                                const cxsc::real& r2) const;

        // use base class for midImage

        //! Destructor
        ~RosenFobj2D(){}


};
#endif


/*! \file
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __EXAMPLEMAPPEDSPFOBJ2D_HPP__
#define __EXAMPLEMAPPEDSPFOBJ2D_HPP__


#include "mappedFobj2D.hpp"


class ExampleMappedFobj2D : public subpavings::MappedFobj2D {

    public:

        //! no argument constructor
        ExampleMappedFobj2D() {};

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r1,
                                                const cxsc::real& r2) const;

        // use base class for midImage

        //! Destructor
        ~ExampleMappedFobj2D(){}


};
#endif

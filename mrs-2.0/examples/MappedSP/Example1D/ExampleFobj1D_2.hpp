
/*! \file
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __EXAMPLEMAPPEDSPFOBJ1D_2_HPP__
#define __EXAMPLEMAPPEDSPFOBJ1D_2_HPP__


#include "mappedFobj1D.hpp"


class ExampleMappedFobj1D_2 : public subpavings::MappedFobj1D {

    public:

        //! no argument constructor
        ExampleMappedFobj1D_2() {};

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~ExampleMappedFobj1D_2(){}


};
#endif

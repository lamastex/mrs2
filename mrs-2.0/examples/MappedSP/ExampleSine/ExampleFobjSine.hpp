
/*! \file
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __EXAMPLEMAPPEDSPFOBJSINE_HPP__
#define __EXAMPLEMAPPEDSPFOBJSINE_HPP__


#include "mappedFobj1D.hpp"



class ExampleMappedFobjSine : public subpavings::MappedFobj1D {

    private:

    //! frequency
    int f;

    public:

        //! no argument constructor
        ExampleMappedFobjSine() :f(1) {};

        //! constructor with frequency
        ExampleMappedFobjSine(int fr) : f(fr) {};

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~ExampleMappedFobjSine(){}


};
#endif

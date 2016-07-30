
/*! \file
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __EXAMPLEMAPPEDSPFOBJSINESUM_HPP__
#define __EXAMPLEMAPPEDSPFOBJSINESUM_HPP__


#include "mappedFobj1D.hpp"

class ExampleMappedFobjSineSum : public subpavings::MappedFobj1D {

    private:

    //! frequency
    int f;

    public:

        //! no argument constructor
        ExampleMappedFobjSineSum() :f(1) {};

        //! constructor with frequency
        ExampleMappedFobjSineSum(int fr) : f(fr) {};

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~ExampleMappedFobjSineSum(){}


};
#endif


/*! \file
\brief Declarations for MappedSPnode 1-d example soscillating function object class

This example is \f$f(x) = (cosx)^a\f$

*/

#ifndef __EXAMPLEMAPPEDSPFOBJOSC1D_2_HPP__
#define __EXAMPLEMAPPEDSPFOBJOSC1D_2_HPP__


#include "mappedFobj1D.hpp"
#include "ExampleFobjOscPI.hpp"


class ExampleFobjOsc1D_2 : public subpavings::MappedFobj1D {

    private:
        int a;

    public:

        //! no argument constructor
        ExampleFobjOsc1D_2();

        //! parameterised constructor
        ExampleFobjOsc1D_2(int aa);

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~ExampleFobjOsc1D_2(){}


};
#endif

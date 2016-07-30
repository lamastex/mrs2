
/*! \file
\brief Declarations for MappedSPnode 1-d example soscillating function object class
*/

#ifndef __EXAMPLEMAPPEDSPFOBJOSC1D_1_HPP__
#define __EXAMPLEMAPPEDSPFOBJOSC1D_1_HPP__


#include "mappedFobj1D.hpp"
#include "ExampleFobjOscPI.hpp"


class ExampleFobjOsc1D_1 : public subpavings::MappedFobj1D {

    private:
        cxsc::real a;
        cxsc::real b;
        cxsc::real c;

    public:

        //! no argument constructor
        ExampleFobjOsc1D_1();

        //! parameterised constructor
        ExampleFobjOsc1D_1(cxsc::real aa, cxsc::real bb, cxsc::real cc);

        //! declare function for cxsc::interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for cxsc::real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~ExampleFobjOsc1D_1(){}


};
#endif

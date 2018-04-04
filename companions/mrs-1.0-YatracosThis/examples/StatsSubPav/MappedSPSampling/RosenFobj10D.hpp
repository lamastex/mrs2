/*! \file RosenFobj10D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __RosenFOBJ10D_HPP__
#define __RosenFOBJ10D_HPP__


#include "mappedFobj10D.hpp"


class RosenFobj10D : public subpavings::MappedFobj10D {

    public:

        //! no argument constructor
        RosenFobj10D() {};

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(
const cxsc::interval& ival1,
const cxsc::interval& ival2,
const cxsc::interval& ival3,
const cxsc::interval& ival4,
const cxsc::interval& ival5,
const cxsc::interval& ival6,
const cxsc::interval& ival7,
const cxsc::interval& ival8,
const cxsc::interval& ival9,
const cxsc::interval& ival10
) const;

        //! declare function for real image of reals
        virtual cxsc::real operator()(
        const cxsc::real& r1,
const cxsc::real& r2,
const cxsc::real& r3,
const cxsc::real& r4,
const cxsc::real& r5,
const cxsc::real& r6,
const cxsc::real& r7,
const cxsc::real& r8,
const cxsc::real& r9,
const cxsc::real& r10
) const;

        // use base class for midImage

        //! Destructor
        ~RosenFobj10D(){}


};
#endif

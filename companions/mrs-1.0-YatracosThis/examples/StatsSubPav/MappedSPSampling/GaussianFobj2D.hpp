/*! \file GaussianFobj2D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __GAUSSIANFOBJ2D_HPP__
#define __GAUSSIANFOBJ2D_HPP__


#include "mappedFobj2D.hpp"


class GaussianFobj2D : public subpavings::MappedFobj2D {

    public:

        //! no argument constructor
        GaussianFobj2D() {};

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const;

        //! declare function for real image of reals
        virtual cxsc::real operator()(const cxsc::real& r1,
                                                const cxsc::real& r2) const;

        // use base class for midImage

        //! Destructor
        ~GaussianFobj2D(){}


};
#endif

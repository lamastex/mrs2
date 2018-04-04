/*! \file GaussianFobj1D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __GAUSSIANFOBJ1D_HPP__
#define __GAUSSIANFOBJ1D_HPP__


#include "mappedFobj1D.hpp"


class GaussianFobj1D : public subpavings::MappedFobj1D {

    public:

        //! no argument constructor
        GaussianFobj1D() {};

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival1) const;

        //! declare function for real image of reals
        virtual cxsc::real operator()(const cxsc::real& r1) const;

        // use base class for midImage

        //! Destructor
        ~GaussianFobj1D(){}


};
#endif

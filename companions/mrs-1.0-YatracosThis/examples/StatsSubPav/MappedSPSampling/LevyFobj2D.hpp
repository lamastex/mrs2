/*! \file LevyFobj2D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __LevyFOBJ2D_HPP__
#define __LevyFOBJ2D_HPP__


#include "mappedFobj2D.hpp"


class LevyFobj2D : public subpavings::MappedFobj2D {

    public:

        //! no argument constructor
        LevyFobj2D() {};

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival1,
                                const cxsc::interval& ival2) const;

        //! declare function for real image of reals
        virtual cxsc::real operator()(const cxsc::real& r1,
                                                const cxsc::real& r2) const;

        // use base class for midImage

        //! Destructor
        ~LevyFobj2D(){}


};
#endif

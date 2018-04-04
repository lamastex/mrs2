/*! \file GaussianFobj9D.hpp
\brief Declarations for MappedSPnode 1-d example function object class
*/

#ifndef __GAUSSIANFOBJ9D_HPP__
#define __GAUSSIANFOBJ9D_HPP__


#include "mappedFobj9D.hpp"


class GaussianFobj9D : public subpavings::MappedFobj9D {

    public:

        //! no argument constructor
        GaussianFobj9D() {};

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
			const cxsc::interval& ival9) const;

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
				const cxsc::real& r9) const;

        // use base class for midImage

        //! Destructor
        ~GaussianFobj9D(){}


};
#endif

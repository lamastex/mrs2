/*
* Copyright (C) 2010 Jennifer Harlow
*
*/

/*! \file 
\brief MappedFobj10D declarations.
*/

#ifndef __MAPPEDFOBJ10D__
#define __MAPPEDFOBJ10D__

#include "mappedFobj.hpp"
#include "cxsc.hpp"


namespace subpavings {
	
	/*! \brief An abstract class for target function objects on 10-dimensional real space.
	 
	 These target functions can be used by MappedSP types.

	*/

    class MappedFobj10D : public MappedFobj{

        public:
            //! a virtual function for interval image of a box
            virtual cxsc::interval operator()(const cxsc::ivector&) const;

            //! a pure virtual function for interval image of an interval pair
            virtual cxsc::interval operator()(
					const cxsc::interval&, const cxsc::interval&,
					const cxsc::interval&, const cxsc::interval&, 
					const cxsc::interval&, const cxsc::interval&,
					const cxsc::interval&, const cxsc::interval&, 
					const cxsc::interval&, const cxsc::interval&) const = 0;

            //! a virtual function for real image of an rvector
            virtual cxsc::real operator()(const cxsc::rvector&) const;

            //! a pure virtual function for real image of a real pair
            virtual cxsc::real operator()(
					const cxsc::real&, const cxsc::real&,
					const cxsc::real&, const cxsc::real&,
					const cxsc::real&, const cxsc::real&,
					const cxsc::real&, const cxsc::real&,
					const cxsc::real&, const cxsc::real&) const = 0;

            //! a virtual function for real image of midpoint of an interval
            virtual cxsc::real imageMid(const cxsc::ivector&) const;

            //! a virtual function for real image of midpoint of an interval pair
            virtual cxsc::real imageMid( 
					const cxsc::interval&, const cxsc::interval&,
					const cxsc::interval&, const cxsc::interval&, 
					const cxsc::interval&, const cxsc::interval&,
					const cxsc::interval&, const cxsc::interval&, 
					const cxsc::interval&, const cxsc::interval&) const;

            //! Destructor
            virtual ~MappedFobj10D(){}


    };

}
#endif

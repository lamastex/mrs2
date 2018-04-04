/*! \file
\brief MappedSPnode example Levy Density 2D example.
*/

#ifndef __LEVYFOBJ_HPP__
#define __LEVYFOBJ_HPP__


#include "mappedFobj.hpp"

/*! The Levy density is the Levy function made into a density.
 * We get the density from energy by exponentiating its negative.*/

class LevyDensityFobj2D : public subpavings::MappedFobj {

    public:

        //! constructor
        LevyDensityFobj2D();

        //! declare function for interval image of a box
        cxsc::interval operator()(const cxsc::ivector& ivec) const;
			
  
        //! declare function for real image of reals
		cxsc::real operator()(const cxsc::rvector& r) const;
        
		std::string getName() const;
		
        // use base class for midImage

        //! Destructor
        ~LevyDensityFobj2D();
	
	private:
	
		const real temperature;
		const real center1; 
		const real center2; 
		const real globalMax;
		



};
#endif

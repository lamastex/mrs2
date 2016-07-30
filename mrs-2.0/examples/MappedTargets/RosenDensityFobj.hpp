/*! \file
\brief Rosenbrock density target.
*/

#ifndef __ROSENDENSITYFOBJ_HPP__
#define __ROSENDENSITYFOBJ_HPP__


#include "mappedFobj.hpp"


class RosenDensityFobj : public subpavings::MappedFobj {

    public:

        //! no argument constructor
        RosenDensityFobj();
		
		 //! constructor
        RosenDensityFobj(real ti, real h);

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::ivector& ivec) const;
			
  
        //! declare function for real image of reals
		virtual cxsc::real operator()(const cxsc::rvector& r) const;
        
		std::string getName() const;
		
        // use base class for midImage

        //! Destructor
        ~RosenDensityFobj();
	
	private:
	
		const real tInverse;
		const real height;



};
#endif

/*! \file
\brief Declarations for Rosenbrock function.
* 
* \note Rosenbrock funtion, not density.
*/

#ifndef __RosenFOBJ_HPP__
#define __RosenFOBJ_HPP__


#include "mappedFobj.hpp"


class RosenFobj : public subpavings::MappedFobj {

    public:

        //! no argument constructor
        RosenFobj();
		
		 //! constructor
        RosenFobj(real ti, real h);

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::ivector& ivec) const;
			
  
        //! declare function for real image of reals
		virtual cxsc::real operator()(const cxsc::rvector& r) const;
        
        std::string getName() const;
		
        // use base class for midImage

        //! Destructor
        ~RosenFobj();
	
	private:
	
		const real tInverse;
		const real height;



};
#endif

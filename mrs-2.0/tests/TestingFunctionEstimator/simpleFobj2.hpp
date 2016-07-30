
/*! \file
\brief Declarations for  multivariate simple  
* function object class.
* 
* f(x) = 2x for 0 <= x <= 1
*/

#ifndef __TESTFUNCEST_SIMPLEFOBJ_HPP__
#define __TESTFUNCEST_SIMPLEFOBJ_HPP__


#include "mappedFobj.hpp"


class SimpleFobj2 : public subpavings::MappedFobj {

    public:

        //! Constructor
        SimpleFobj2();

        /*! declare function for cxsc::interval image of a box
		 * \return product over dimensions [d] interval ( operator(inf(ivec[d])), 
		 * 						 operator(sup(ivec[d])) ) 
		 * */
        virtual cxsc::interval operator()(const cxsc::ivector& ivec) const;

        /*! declare function for cxsc::real image of rvector
		 * \return product over dimensions[d] (2xr[d]) */
        virtual cxsc::real operator()(const cxsc::rvector& r) const;

        // use base class for midImage

        //! Destructor
        ~SimpleFobj2(){}
		
		private:
		
		/*! declare function for cxsc::real image of real
		* \return 2r */
        cxsc::real operator()(const cxsc::real& r) const;
		
		std::string name;
		

};
#endif

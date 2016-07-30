
/*! \file
\brief Declarations for multivariate simple  
* function object class.
* 
* f(x) = x for 0 <= x <= 1
*/

#ifndef __TESTFUNCEST_UNIFORMFOBJ_HPP__
#define __TESTFUNCEST_UNIFORMFOBJ_HPP__


#include "mappedFobj.hpp"


class SimpleFobj1 : public subpavings::MappedFobj {

    public:

        //! Constructor
        SimpleFobj1();

        /*! declare function for cxsc::interval image of a box
		 * \return product over dimensions [d] interval ( operator(inf(ivec[d])), 
		 * 						 operator(sup(ivec[d])) ) 
		 * */
        virtual cxsc::interval operator()(const cxsc::ivector& ivec) const;

        /*! declare function for cxsc::real image of rvector
		 * \return product over dimensions[d] (r[d]) */
        virtual cxsc::real operator()(const cxsc::rvector& r) const;

        // use base class for midImage

        //! Destructor
        ~SimpleFobj1(){}
		
		private:
		
		/*! declare function for cxsc::real image of real
		* \return r */
        cxsc::real operator()(const cxsc::real& r) const;
		
		std::string name;
		

};
#endif

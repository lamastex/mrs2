
/*! \file
\brief Declarations for multivariate oscillating
* function object class.
* 
*/

#ifndef __TESTFUNCEST_OSCFOBJ_HPP__
#define __TESTFUNCEST_OSCFOBJ_HPP__


#include "mappedFobj.hpp"


class OscFobj : public subpavings::MappedFobj {

    public:

        //! Constructor
        OscFobj();

        /*! declare function for cxsc::interval image of a box
		 * */
        cxsc::interval operator()(const cxsc::ivector& ivec) const;

        /*! declare function for cxsc::real image of rvector*/
		cxsc::real operator()(const cxsc::rvector& r) const;

		std::string getName() const;

        // use base class for midImage

        //! Destructor
        ~OscFobj(){}
		
		private:
		
		cxsc::interval operator()(const cxsc::interval& ival) const;
		
		cxsc::real operator()(const cxsc::real& r) const;
		
		std::string name;
		
		cxsc::real a;
        cxsc::real b;
        cxsc::real c;


};
#endif

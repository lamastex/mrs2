
/*! \file
\brief Declarations for multivariate sphere
* function object class.
* 
*/

#ifndef __SPHEREFOBJ_HPP__
#define __SPHEREFOBJ_HPP__


#include "mappedFobj.hpp"


class SphereFobj : public subpavings::MappedFobj {

    public:

        //! No-args constructor
        SphereFobj();
		
		//! Constructor
		SphereFobj(const cxsc::rvector& c);

        /*! declare function for cxsc::interval image of a box
		 * */
        cxsc::interval operator()(const cxsc::ivector& ivec) const;

        /*! declare function for cxsc::real image of rvector*/
		cxsc::real operator()(const cxsc::rvector& r) const;

		std::string getName() const;

        // use base class for midImage

        //! Destructor
        ~SphereFobj(){}
		
		private:
		
		const cxsc::rvector centre;
		
		const int cLen;
		
		std::string name;

};
#endif

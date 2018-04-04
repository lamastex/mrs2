/*! \fileFinMixFobj.hpp
\brief Declarations for MappedSPnode 1D Gaussian Mixtures function object class
*/

#ifndef __FINMIXFOBJ_HPP__
#define __FINMIXFOBJ_HPP__

#include "mappedFobj1D.hpp"
#include <vector>


class FinMixFobj : public subpavings::MappedFobj1D {
    private:
		std::vector<double> W;
		std::vector<double> M;
		std::vector<double> S;

	 public:

        //! no argument constructor
        FinMixFobj();

        //! parameterised constructor
		  FinMixFobj(std::vector<double> WW, std::vector<double> MM, 
		  					std::vector<double> SS);

        //! declare function for interval image of a box
        virtual cxsc::interval operator()(const cxsc::interval& ival) const;

        //! declare function for real image of reals
        virtual cxsc::real operator()(const cxsc::real& r) const;

        // use base class for midImage

        //! Destructor
        ~FinMixFobj(){}


};
#endif

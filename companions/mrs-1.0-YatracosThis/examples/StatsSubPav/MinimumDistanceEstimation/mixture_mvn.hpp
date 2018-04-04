/* **************************************************************************************
 *  Class to find density for a mixture of Multivariate Normal density functions
 *  and also to generate random values from the mixture distribution
 * 
 * 
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Based on the original by Ralph dos Santos Silva
 *  see http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html
 *  Copyright (C) 2006  Ralph dos Santos Silva, modified
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation
***************************************************************************************/

#ifndef __MIXTURE_MVN_HPP__
#define __MIXTURE_MVN_HPP__

#include <vector>
#include <cstddef>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "cxsc.hpp"

namespace subpavings { 
	
	namespace kde {
		
		class MixtureMVN {
			
			public:
		
				MixtureMVN(const std::vector < std::vector < double > >& m,
							const std::vector < std::vector < std::vector < double > > >& s,
							const std::vector < double >& mx,
							unsigned long int rsd);
							
				~MixtureMVN();
				
				void resetPRNG(unsigned long int rsd);
				
				void resetPRNG();
				
				std::vector < std::vector < double > >&
					prn(std::vector < std::vector < double > >& rvs,
							size_t n) const;
							
				std::vector < std::vector < real > >&
					prn(std::vector < std::vector < real > >& rvs,
							size_t n) const;
				
				std::vector < cxsc::rvector >&
					prn(std::vector < cxsc::rvector >& rvs,
							size_t n) const;
				
				std::vector < double >& prn(std::vector < double >& rv) const;
				
				std::vector < real >& prn(std::vector < real >& rv) const;
				
				cxsc::rvector& prn(cxsc::rvector& rv) const;
				
				real f(const std::vector < double >& dv) const;
				
				real f(const std::vector < real >& rv) const;
				
				real f(const cxsc::rvector& rv) const;
			
			private:
			
				MixtureMVN();
				
				MixtureMVN(const MixtureMVN& other );
				
				MixtureMVN& operator=(const MixtureMVN& rhs);
			
				size_t getDistIndex() const;
				
				void clean();
				
				
				void prn_gsl_vec(size_t i, 	gsl_vector *result) const;

				double d_gsl_vec(size_t i, const gsl_vector *x) const;
				
				size_t mixn;
				size_t dim;
				std::vector < gsl_vector* > means;
				std::vector < gsl_matrix* > scales;
				std::vector < double > mixes;
				unsigned long int seed; 
				gsl_rng* r;
			
			

		};
	
	}
}

#endif

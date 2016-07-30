/*
* Copyright (C) 2012 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


/*! \file
\brief Declarations for a class for making a kernel density 
* estimate using a data-adaptive method with a MCMC.

*/

#ifndef __MCMCKDE_HPP__
#define __MCMCKDE_HPP__

	#include "type_kde.hpp"
	#include "cxsc.hpp"

	#include <vector>
	#include <string>
	#include <cstddef> // size_t

namespace subpavings {
	namespace kde {
		
		/*! \brief A class for making a kernel density 
			* estimator using a data-adaptive method with a MCMC to get
			* the bandwidth estimates.  See
			* (Zhang, X., King, Maxwell, Hyndman, Rob. J. (2006), 'A Bayesian approach
			* to bandwidth selection for multivariate kernel density estimation',
			* Computational Statistics and Data Analysis, vol. 50, pp. 3009--3031).
			* 
			* The author of this class gratefully acknowledges that 
			* this is based on C code supplied by Profs. Zhang and King, 
			* who also granted permission to the author to use and adapt
			* this code for the author's Masters Thesis, which incorporated
			* the development of this library.  
			* 
			* \note At present the tuning parameter used is a constant
			* and cannot be set by the user. 

			*/
		class MCMCKDE : public TypeKDE {
			
			public :
			
				/*! \brief Constructor using reals.
				 * 
				 * \param dx the collection of sample data from which
				 *  to create the kernel density estimator */
				MCMCKDE(const std::vector < std::vector < real > >& dx);
				
				/*! \brief Constructor using doubles.
				 * 
				 * \param dx the collection of sample data from which
				 *  to create the kernel density estimator */
				MCMCKDE(const std::vector < std::vector < double > >& dx);
				
				/*! \brief Destructor. */
				~MCMCKDE(){}
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of cxsc::real values. 
				 * 
				 * \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const std::vector < cxsc::real >& x) const;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a vector of double values. 
				 * 
				  \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const std::vector < double >& x) const;
				
				/*! \brief Estimate kde(x), the kernel density estimator
				 * function evaluated at \a x. where \a x is given 
				 * as a cxsc::rvector value.
				 * 
				 * \note df method should be run to estimate bandwidths
				 * before kde can be evaluated. */
				cxsc::real kde(const cxsc::rvector& x) const;
				
				/*! \brief Estimate bandwidths using given parameters.
				 * 
				 * \param size_batch the number of iterations in 
				 * each batch to use 
				 * (the number of batches is used for the diagnostic
				 * information output to the file \a logFilename).
				 * \param num_batch the number of batches to use 
				 * (the number of batches is used for the diagnostic
				 * information output to the file \a logFilename).
				 * \param warm the length of the warm-up period to use, 
				 * in terms of the number of warm-up iterations.  
				 * \param size_t step controls the step used for 
				 * logging trace results to the file \a logFilename
				 * \param logFilename the name of a file to which to send
				 * logging output used to record bandwidths and diagnostic
				 * information.
				 * \return the collection of bandwidths used.  These 
				 * are the entries on the diagonal of the diagonal
				 * bandwidth matrix created by this routine.*/
				std::vector < cxsc::real > df(	size_t size_batch, 
									size_t num_batch, 
									size_t warm,
									size_t step, // just controls step for results logging, ie traces
									const std::string& logFilename,
									unsigned int seed) const;
									
									
				/*! \brief Estimate bandwidths using given parameters.
				 * 
				 * \param size_batch the number of iterations in 
				 * each batch to use 
				 * (the number of batches is used for the diagnostic
				 * information output to the file \a logFilename).
				 * \param num_batch the number of batches to use 
				 * (the number of batches is used for the diagnostic
				 * information output to the file \a logFilename).
				 * \param warm the length of the warm-up period to use, 
				 * in terms of the number of warm-up iterations.  
				 * \param size_t step controls the step used for 
				 * logging trace results to the file \a logFilename
				 * \param logFilename the name of a file to which to send
				 * logging output used to record bandwidths and diagnostic
				 * information.
				* \param resultsFilename the name of a file to which to send
				 * trace information.
				 * \return the collection of bandwidths used.  These 
				 * are the entries on the diagonal of the diagonal
				 * bandwidth matrix created by this routine.*/
				std::vector < cxsc::real > df(size_t size_batch,
						size_t num_batch,
						size_t warm,
						size_t step, // just controls step for results logging, ie traces
						const std::string& logFilename,
						const std::string& resultsFilename,
						unsigned int seed) const;
				
					

			private :

				cxsc::real _kde(const std::vector < cxsc::real >& rx) const;
								
				cxsc::real cost(const std::vector < cxsc::real >& h) const;

				cxsc::real kn_gibbs(
						cxsc::real xCost,
						std::vector < cxsc::real >& x)const;
				
				real gasdev() const;
				
				void rescaleData() const;
				
				void logResults(const std::string& resultsFilename,
					const std::vector < cxsc::real >& costs,
					const std::vector <std::vector < cxsc::real > >& results,
					size_t warm, size_t step) const;
		
				
				/* "tuning parameter mentioned in paper? */
				const static cxsc::real mutsizp;


				mutable size_t accept_h;
				mutable size_t total_h;
				mutable int iset;
				mutable cxsc::real gset;
				std::vector < std::vector < real > > data_x;
				size_t dim;
				size_t data_num;
				cxsc::real cont; 
				mutable std::vector < cxsc::real > bandwidths;
				mutable cxsc::real productBandwidths;
				mutable std::vector < std::vector < real > > rescaled_x;
				
		}; // end class
	}

} // end namespace


#endif

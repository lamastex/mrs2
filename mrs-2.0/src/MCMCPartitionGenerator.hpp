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
\brief MCMCPartitionGenerator declarations.

*/

#ifndef ___MCMC_PARTITION_HPP__
#define ___MCMC_PARTITION_HPP__

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

#include "real.hpp"

namespace subpavings {
	
	/*! \brief A class to be able to generate a random partition
	of a binary tree into a given number of leaves.

   */
    class MCMCPartitionGenerator {

    public:
        /*! \brief Default constructor.
		 * 
		 * \param seed a seed for a random number generator.        */
        explicit MCMCPartitionGenerator(unsigned long int seed);

        
        /*! \brief Destructor.    */
        ~MCMCPartitionGenerator();

        
		/*! \brief Generate a partition of \a numLeaves to return
		 * a value in {1, 2, ..., numLeaves-1}, each with
		 * equal probability.
		 * 
		 * The value returned can be thought of as the number going to 
		 * 'one side' of a partition of \a numLeaves, where the 
		 * allowable partitions, as pairs (one side, other side) are 
		 * (1, numLeaves-1), (2, numLeaves-2), ... (numLeaves-1, 1).
		 * 
		 * \internal This method uses approximations to the Catalan
		 * numbers to be able to get equally probable partitions.  The
		 * approximation may mean that function has to have several 
		 * attempts to find a return value.  If a non-empty 
		
		\param numLeaves the number of leaves to partition.
		\return A value in {1, 2, ..., numLeaves-1}.*/
		
		
		unsigned long int generateStatePartition(
			unsigned long int numLeaves) const;
		
		//calc sum(ln(1/prob))
		unsigned long int generateNaturalStatePartition(
				unsigned long int numLeaves) const;
		
		int generateKnuthDecision(unsigned long int p,
									unsigned long int q,
									bool across = false) const;
		
		cxsc::real getLnCatalanRatio(
			unsigned long int k, unsigned long int kPrime) const;
			
		void initialiseInstructions(
				size_t numLeaves ) const;
	
	
		std:: vector < unsigned long int > getInstructions() const;
		
		void clearInstructions() const;
		
		
	private:
	
		MCMCPartitionGenerator();
		
		MCMCPartitionGenerator(const MCMCPartitionGenerator& other);

		MCMCPartitionGenerator& operator=(MCMCPartitionGenerator tmp);
	
		void makeCatalans();
		
		
		
		gsl_rng* rgsl;
		
		std:: vector < unsigned long int > catalans;
		
		mutable std:: vector < unsigned long int >* instructionsPtr;
		mutable std:: vector < unsigned long int >* workspacePtr;	
			
        
    };
    // end of MCMCPartitionGenerator class
	

} // end namespace subpavings




#endif

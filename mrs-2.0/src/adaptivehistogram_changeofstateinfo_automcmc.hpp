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
\brief Declarations for type collecting info for automcmc from adaptive histograms.
*/

#ifndef __CHANGEOFSTATEINFO_AUTOMCMC_HPP__
#define __CHANGEOFSTATEINFO_AUTOMCMC_HPP__

#include "adaptivehistogram.hpp"

namespace subpavings {

	/*! \brief Type for collecting change of state information
	 * for use with automated MCMC processes.
	 * 
	 * Information collected and made available
	 * is that usable by automated
	 * MCMC processes, including Gelman-Rubin diagnostic. */ 
	class ChangeOfStateInformationAutoMCMC 
			: public AdaptiveHistogram::ChangeOfStateInformation {
			
		public:
			ChangeOfStateInformationAutoMCMC(real lp,
											size_t cl, size_t cc,
											unsigned long int tld);
			
			ChangeOfStateInformationAutoMCMC(real lp,
					const AdaptiveHistogram& adh);
			
			ChangeOfStateInformationAutoMCMC(real lp,
										size_t cl, size_t cc, 
										const AdaptiveHistogram& adh);
			
			/*! \brief Get change in log posterior from last change. */			
			real getDeltaPi() const;
			/*! \brief Get current log posterior from last change. */			
			real getCurrentLogPosterior() const;
			/*! \brief Get current number of leaves. */
			size_t getCurrentLeaves() const;
			/*! \brief Get current number of cherries. */
			size_t getCurrentCherries() const;
			/*! \brief Get current average leaf depth. */
			real getAverageLeafDepth() const;
		
		private:	
		
			ChangeOfStateInformationAutoMCMC();
			
			/*! \brief Notify this of change in log posterior from last change. */			
			void notifyDeltaPi(real dp);
			/*! \brief Notify this of a split. */
			void notifySplit(const SPnode * const spn);
			/*! \brief Notify this of a merge. */
			void notifyMerge(const SPnode * const spn);
			
			/*! \brief Current log posterior. */
			real logPosterior;
			/*! \brief Change in log posterior from last change. */
			real deltaPi;
			/*! \brief Current number of leaves. */
			size_t currentLeaves;
			/*! \brief Current number of cherries. */
			size_t currentCherries;
			/*! \brief Current total leaf depth. */
			unsigned long int totalLeafDepth;
		
	};
}

#endif

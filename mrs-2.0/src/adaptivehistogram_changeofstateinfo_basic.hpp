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
\brief Declarations for type collecting basic info from adaptive histograms.
*/

#ifndef __CHANGEOFSTATEINFO_BASIC_HPP__
#define __CHANGEOFSTATEINFO_BASIC_HPP__

#include "adaptivehistogram.hpp"

namespace subpavings {

	/*! \brief Type for collecting change of state information. */ 
	class ChangeOfStateInformationBasic 
			: public AdaptiveHistogram::ChangeOfStateInformation {
		
		public:
			ChangeOfStateInformationBasic(size_t cl, size_t cc);
			
			ChangeOfStateInformationBasic(
					const AdaptiveHistogram& adh);
			
			/*! \brief Change in log posterior from last change. */			
			real getDeltaPi() const;
			/*! \brief Get current number of leaves. */
			size_t getCurrentLeaves() const;
			/*! \brief Get current number of cherries. */
			size_t getCurrentCherries() const;
			
		private:
		
			ChangeOfStateInformationBasic();
			
			/*! \brief Notify this of change in log posterior from last change. */			
			void notifyDeltaPi(real dp);
			/*! \brief Notify this of a split. */
			void notifySplit(const SPnode * const spn);
			/*! \brief Notify this of a merge. */
			void notifyMerge(const SPnode * const spn);
			
			/*! \brief Change in log posterior from last change. */
			real deltaPi;
			/*! \brief Current number of leaves. */
			size_t currentLeaves;
			/*! \brief Current number of cherries. */
			size_t currentCherries;
			
		
	};
}

#endif

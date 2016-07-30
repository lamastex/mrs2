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
\brief Declarations for a type to do a combination of 
* carving and SEB priority queue to get starting points for MCMC
 */

#ifndef __CARVER_SEB_HPP__
#define __CARVER_SEB_HPP__

#include "adaptivehistogram.hpp"  // headers for the histograms
#include "histmcmcobjs.hpp"


#include <vector>
#include <string>

namespace subpavings {
	
	/*! \brief A static class for doing carving-SEB RPQ routines. 
	 * 
	 * \note Any pointers in \a hists at the end of the routines shown
	 * here are pointers to objects in dynamic memory (ie, newed).  
	 * These objects will need to be deleted at the end of the routine
	 * using them.  
	 * 
	 * \todo Would like to clean up this aspect of the design.  
	 * If we used boost::shared_pointer a lot of horrible things
	 * like this could be dealt with, but it adds yet another library
	 * dependency to the whole thing. */
	class CarverSEB {
		
		public :

			/*! @name Find some of the best starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as the 
			 * \a keepBest best from \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Best is defined in terms of 
			 * highest log-posterior mass.  
			 * 
			 * \internal 
			 * The original method.  Equivalent to findStartingPointsBest.
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using evaluatorCarving) and 
			 * \a carvingStarts are regarded as initial values only:
			 * if no maximums are found using these values they 
			 * may be adjusted within the method and the search proces
			 * repeated.
			 * 
			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  
	 
			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  Note that this is changed during the 
			 * operation.
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with a SEB
			 * queue will be made.
			 * \param keepBest is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to 
			 * the maximum posterior states
			 * found.  It is not guaranteed that there will be 
			 * \a keepBest of these, or indeed that there will be any.  */
			
			//@{
			
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPoints(
					AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior, size_t minPoints,
					int carvingStarts, size_t keepBest,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			
			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPoints(
					AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					int carvingStarts, size_t keepBest,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			//@}
			
			/*! @name Find some of the best starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as the 
			 * \a keep best from \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Best is defined in terms of 
			 * highest log-posterior mass.  
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using evaluatorCarving) and
			 * \a carvingStarts are regarded as initial values only:
			 * may be adjusted within the method and the search proces
			 * repeated.

			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  
			 * 
			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with a SEB
			 * queue will be made.
			 * \param keep is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to 
			 * the maximum posterior states
			 * found.  It is not guaranteed that there will be 
			 * \a keep of these, or indeed that there will be any.  */
			
			//@{
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsBest(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsBest(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			
			//@}
					
			/*! @name Find the most spread-out starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as the 
			 * \a keep most spread out maximums from 
			 * \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Spread is defined in terms
			 * of the number of leaves in the histograms.  
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using evaluatorCarving) and 
			 * \a carvingStarts are regarded as initial values only:
			 * if no maximums are found using these values they 
			 * may be adjusted within the method and the search proces
			 * repeated.

			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  
			 * 
			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with a SEB
			 * queue will be made.
			 * \param keep is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to 
			 * the well-spread-out states
			 * found.  It is not guaranteed that there will be 
			 * \a keep of these, or indeed that there will be any.  */
			
			//@{
			
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsMaxLeafSpread(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsMaxLeafSpread(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			//@}

			/*! @name Find some well spread-out starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as the 
			 * \a keep most spread out maximums from 
			 * \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Spread is defined in terms
			 * of the number of leaves in the histograms.  
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using evaluatorCarving) and 
			 * \a carvingStarts are regarded as initial values only:
			 * if no maximums are found using these values they 
			 * may be adjusted within the method and the search proces
			 * repeated.
			 * 
			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  

			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with an SEB
			 * queue will be made.
			 * \param keep is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param outOfClosest is the size of the subset
			 * of the maximum posterior
			 * points found from which to select the \a keep most
			 * spread-out points.  This using this can prevent 
			 * the starting points being so far apart that
			 * chains take too long to converge .
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to 
			 * the well-spread-out states
			 * found.  It is not guaranteed that there will be 
			 * \a keep of these, or indeed that there will be any.  */
			
			//@{
			
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsMaxLeafSpread(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					int carvingStarts,
					size_t keep,
					size_t outOfClosest,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			
			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsMaxLeafSpread(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					int carvingStarts,
					size_t keep,
					size_t outOfClosest,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			//@}
			
			/*! @name Find some over-dispersed starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as
			 * \a keep well diversified points from the
			 * \a sequence of states created by the carving-seb
			 * RPQ that gives the overall maximum log-posterior 
			 * point in  
			 * \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Diversification is defined in terms
			 * of the number of leaves in the histograms.  
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using evaluatorCarving) and 
			 * \a carvingStarts are regarded as initial values only:
			 * if no maximums are found using these values they 
			 * may be adjusted within the method and the search proces
			 * repeated.
			 * 
			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  

			 * The method identifies the the carving-seb
			 * RPQ that gives the overall maximum log-posterior 
			 * point in \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  It then identifies the 
			 * maximum log posterior state in this sequence and selects
			 * \a keep well diversified states from a sub-sequence of
			 * states in the sequence that contains the maximum
			 * log-posterior state and starts at the first state
			 * that achieves \a \f$ \alpha \f$ of that maximum.    
			 * \f$ \alpha \f$ is set within the method.
			 * 
			 * \internal \f$ \alpha \f$ = 0.95. 
			 * 
			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with a SEB
			 * queue will be made.
			 * \param keep is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to \a keep
			 * the over-diversified states found.  */
			
			//@{
				
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsOverdispersed(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsOverdispersed(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					int carvingStarts,
					size_t keep,
					bool stopOnMaxPosterior,
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);
			
			//@}


			/*! \brief Find some over-dispersed starting points.
			 * 
			 * Fills \a hists with pointers to histograms created as
			 * \a keep well diversified points from the
			 * \a sequence of states created by the carving-seb
			 * RPQ that gives the overall maximum log-posterior 
			 * point in  
			 * \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  Diversification is defined in terms
			 * of the number of leaves in the histograms.  
			 * 
			 * \note The maximum leaves in the carving queue 
			 * (specified using \a evaluatorCarving) and 
			 * \a carvingStarts are regarded as initial values only:
			 * they 
			 * may be adjusted within the method and the search proces
			 * repeated if no satisfactory maximums are found using
			 * the original values.  Any such adjustment process
			 * must be finite, so that the method will eventually
			 * always either progress or fail with an exception (not
			 * continue indefinitely).
			 * 
			 * \note Any pointers in \a hists at the end of the routine
			 * are pointers to objects in dynamic memory (ie, newed).  
			 * These objects will need to be deleted at the end of the routine
			 * using them.  

			 * \internal In the implementation used here the maximum
			 * leaves specified in \a evaluatorCarving may be
			 * increased and another attempt to find a maximum made
			 * if the maximum posterior point is found either
			 * by an SEB queue straight from the root or at the maximum
			 * carving value tried.  The maximum value to which the 
			 * maximum leaves in the carving queue may be so set 
			 * within the method is half of the maximum leaves in the
			 * SEB queue, and a maximum of 2 attempts to find a better
			 * maximum posterior point with adjusted values may be 
			 * made (i.e., 3 attempts to find a maximum log posterior
			 * point altogether).  
			 *  
			 * The method identifies the the carving-seb
			 * RPQ that gives the overall maximum log-posterior 
			 * point in \a carvingStarts+1 attempts at a
			 * combined carving-SEB rpq.  It then identifies the 
			 * maximum log posterior state in this sequence and selects
			 * \a keep well diversified states from a sub-sequence of
			 * states in the sequence that contains the maximum
			 * log-posterior state and starts at the first state
			 * that achieves \a percentSpread of that maximum. 
			 * 
			 * \param adh is the histogram to use to find the best 
			 * staritiong points.  
			 * \param hists is a container to fill with pointers to
			 * AdaptiveHistogram found during the operation.  Note that
			 * these are on the heap and will need to be deleted.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue and also the SEB-related stopping criteria, 
			 * the maximum number of points in any leaf below which 
			 * no further splitting will take place.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume of a node 
			 * to control which nodes are considered to be splittable: a node
			 * must have volume >= 2* minVol to be split, i.e when split
			 * each child will have volume >= minVol.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of any %AdaptiveHistogram in \a hists after the operation. 
			 * \param carvingStarts is the number of different 
			 * lengths of the carving queue to use before the SEB
			 * queue. In addition to \a carvingStarts attempts, one
			 * attempt with no carving at all is always made, i.e, 
			 * setting \a carvingStarts = 0 will mean that just 
			 * one attempt straight from the root with a SEB
			 * queue will be made.
			 * \param keep is the number of histogram pointers 
			 * to try to store in \a hists.
			 * \param percentSpread the fraction of the maximum 
			 * log-posterior to use to identify the first state
			 * in the sub-sequence from which the \a keep 
			 * states to be returned are selected.  Note, suitable
			 * value are about 0.8-0.95.  
			 * \param stopOnMaxPosterior is an indicator for whether
			 * the search for maximum posterior points should continue
			 * until \a evaluatorSEB's maximum leaves is reached 
			 * or terminate once it seems clear that a maximum has been
			 * found.  If false, the search will continue until 
			 * a evaluatorSEB's maximum leaves is reached.  If true
			 * the search may continue for some relatively short 
			 * number of states past a local maximum to ensure that it 
			 * it appears to be a global maximum.  
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param checkPostFileNameBase is the base for the 
			 * filenames to which to output
			 * the log from the creation of each histogram found
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain pointers to \a keep
			 * the over-diversified states found.  */
			
			//@{
			
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsOverdispersed(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					int carvingStarts,
					size_t keep,
					double percentSpread,
					bool stopOnMaxPosterior, //advise false
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& findStartingPointsOverdispersed(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					int carvingStarts,
					size_t keep,
					double percentSpread,
					bool stopOnMaxPosterior, //advise false
					const std::string& postFileName,
					const std::string& checkPostFileNameBase,	
					int prec,
					unsigned long int seed = 1234);

			//@}
			
			/*! @name Add a starting point to a container
			 * of starting points.
			 * 
			 * Adds one more pointer to \a hists, created using the 
			 * carving and SEB queue parameters specified in 
			 * \a evaluatorCarving and \a evaluatorSEB.  
			 * 
			 * \note the pointer added to \a hists will be a pointer
			 * to an object in dynamic memory and will need to 
			 * deleted at the end of the routine using it. 
			 * 
			 * \param adhBase is the histogram to use to find the best 
			 * staritiong points.  Note that this is changed during the 
			 * operation.
			 * \param hists is a container to which to add a 
			 * pointer to the new start point created in this method.  
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param evaluatorSEB is the PrioritySplitQueueEvaluator
			 * to use to control the SEB queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * SEB queue.  Note that if it also specifies a positive
			 * SEB-related stopping criteria, the maximum number of 
			 * leaves specified may not be reached.  
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of the %AdaptiveHistogram added to \a hists. 
			 * \param checkPostFileNameBase is the base for the 
			 * filename to which to output
			 * the log from the creation of the new histogram
			 * to be pointed to in \a hists carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a hists will contain one more pointer, to a state
			 * created using \a evaluatorCarving and \a evaluatorSEB.  */
			
			//@{
				
			/*! \brief Version without \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& createStartingPoint(
					const AdaptiveHistogram& adhBase,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, // correct carvedPt
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB,  // correct maxPostPt
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					const std::string& checkPostFileNameBase,	
					unsigned long int seed = 1234);
			
			/*! \brief Version with \a minVol argument. */
			static std::vector< AdaptiveHistogram* >& createStartingPoint(
					const AdaptiveHistogram& adhBase,
					std::vector< AdaptiveHistogram* >& hists,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, // correct carvedPt
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB,  // correct maxPostPt
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					const std::string& checkPostFileNameBase,	
					unsigned long int seed = 1234);
			
			//@}
			
			/*! @name Carve the given %AdaptiveHistogram.
			 * 
			 * \param adh is the histogram to carve. Note that this
			 * is changed during the 
			 * operation.
			 * \param evaluatorCarving is the PrioritySplitQueueEvaluator
			 * to use to control the carving queue.  It should specify
			 * the maximum number of leaves to which to attempt the 
			 * carving queue.
			 * \param logPrior the prior to uses to evaluate the 
			 * log-posterior.
			 * \param minPoints the minimum number of points in a node 
			 * to control which nodes are considered to be splittable.
			 * \param minVol the minimum volume that should be in any leaf node
			 * of \a adh after the operation. 
			 * \param postFileName is the filename to which to output
			 * the log from the carving queue.  If this is the empty
			 * string ("") no log will be output.
			 * \param prec is the precision to use when outputting logs.
			 * \param seed is the seed to use for the carver-SEB RPQ 
			 * process. Defaults to 1234.
			 * \post \a adh will be carved as specified by 
			 * \a evaluatorCarving.  */
			
			//@{
				
			/*! \brief Version without \a minVol argument. */
			static AdaptiveHistogram& carve(
								AdaptiveHistogram& adh,
								AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
								LogMCMCPrior& logPrior, size_t minPoints,
								const std::string& postFileName,
								int prec,
								unsigned long int seed = 1234);
			
			/*! \brief Version with \a minVol argument. */
			static AdaptiveHistogram& carve(
								AdaptiveHistogram& adh,
								AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
								LogMCMCPrior& logPrior, 
								size_t minPoints,
								double minVol,
								const std::string& postFileName,
								int prec,
								unsigned long int seed = 1234);
			
			//@}

			
		private :
		
			CarverSEB(); // private and not implemented
			
			static size_t findMax(const std::vector<real>& maxPosteriors);
			
			static std::vector< AdaptiveHistogram* >& addToHists(
					const AdaptiveHistogram& adh,
					std::vector< AdaptiveHistogram* >& histsToAddTo,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
					LogMCMCPrior& logPrior,
					size_t minPoints,
					double minVol,
					size_t toFind,
					double dblGap,
					size_t baseIndex,
					const std::vector < real >& posteriorVec,
					size_t startLeaves,
					size_t carvedPt,
					const std::string& checkPostFileNameBase,	
					const gsl_rng * rgsl_base);

			static size_t findMaxPosteriorPoints(
					const AdaptiveHistogram& adh,
					std::vector<size_t>& carvedLaunchPoints,
					std::vector<size_t>& maxPosteriorPoints,
					std::vector<real>& maxPosteriors,
					AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					int carvingStarts,
					bool stopOnMaxPosterior,
					const gsl_rng * rgsl_base,
					const std::string& postFileName,
					int prec	);
					
					
			static void outputMaxPosteriorPoints(
					std::ostream& os,
					const std::vector<size_t>& carvedLaunchPoints,
					const std::vector<size_t>& maxPosteriorPoints,
					const std::vector<cxsc::real>& maxPosteriors);

			static bool essentiallyEqual(cxsc::real r1, cxsc::real r2);
					
			static void cleanMaxPoints(
					std::vector<size_t>& carvedLaunchPoints,
					std::vector<size_t>& maxPosteriorPoints,
					std::vector<real>& maxPosteriors);


			static std::vector < size_t > keepFurthestApart(
							const std::vector<size_t>& carvedLaunchPoints,
							const std::vector<size_t>& maxPosteriorPoints,
							const std::vector<cxsc::real>& maxPosteriors,
							size_t keep,
							size_t outOfClosest);

			static std::vector < size_t > keepBest(
							const std::vector<size_t>& carvedLaunchPoints,
							const std::vector<size_t>& maxPosteriorPoints,
							const std::vector<real>& maxPosteriors,
							size_t keep);

			static std::vector< AdaptiveHistogram* >& recreateStartingPoints(
					const AdaptiveHistogram& adhBase,
					std::vector< AdaptiveHistogram* >& hists,
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorCarving, 
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					const std::string& checkPostFileNameBase,	
					const std::vector<size_t>& carvedLaunchPoints,
					const std::vector<size_t>& maxPosteriorPoints,
					const std::vector<cxsc::real>& maxPosteriors,
					const std::vector < size_t >& keepTheseOnes,
					const gsl_rng * rgsl_base);
					

			static bool recreateStartingPoint(
					AdaptiveHistogram& adh, // changed
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorCarving, 
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					size_t carvedPt,
					size_t maxPostPoint,
					std::vector<cxsc::real>& posteriorVec, 
					std::vector<cxsc::real>& loglikVec, 
					const gsl_rng * rgsl_base,
					LOGGING_LEVEL logPQ);


			static bool recreateStartingPoint(
					AdaptiveHistogram& adh, // changed
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorCarving, 
					const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
					LogMCMCPrior& logPrior, 
					size_t minPoints,
					double minVol,
					std::vector<cxsc::real>& posteriorVec, 
					std::vector<cxsc::real>& loglikVec, 
					const gsl_rng * rgsl_base,
					LOGGING_LEVEL logPQ);
			
			
					
	}; // end class SEBCarver
} // end namespace subpavings
#endif



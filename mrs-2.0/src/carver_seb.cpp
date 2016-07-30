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
\brief My version of Gloria's carving PQ to get starting points for MCMC
 */

#include "carver_seb.hpp"


#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <cassert>
#include <iterator>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <utility>
#include <stdexcept>

#include <gsl/gsl_math.h> //gsl_isnan(), isinf

//#define DOSAMPLES // get samples for each stage in the carver queue
#define DEBUG_MAXPOINTS
//#define DEBUG_MAXPOINTS_OVERDISPERSED

using namespace std;

namespace subpavings {

	

	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPoints(
			AdaptiveHistogram& adh,
			std::vector< AdaptiveHistogram* >& hists,
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			int carvingStarts, size_t keepBest,
			bool stopOnMaxPosterior,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,	
			int prec,
			unsigned long int seed)
	{
		double minVol = 0.0;
		
		return findStartingPoints(
			adh,
			hists,
			evaluatorCarving, 
			evaluatorSEB, 
			logPrior, 
			minPoints,
			minVol,
			carvingStarts,
            keepBest,
			stopOnMaxPosterior,
			postFileName,
			checkPostFileNameBase,	
			prec,
			seed);
	}
			
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPoints(
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
			unsigned long int seed)
	{
		if (carvingStarts < static_cast<int>(keepBest)-1) 
			throw std::invalid_argument(
						"findStartingPoints:\nchooseStarts < keepBest-1");

		gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
		gsl_rng_set (rgsl, seed); // change the seed
		
		// clone it
		gsl_rng * rgsl_base = gsl_rng_clone(rgsl);

		// parameters needed to start the PQ
		LOGGING_LEVEL logPQ = NOLOG; // logging level  
		
		size_t ml = evaluatorCarving.getMaxLeaves();
		int checkMaxStep = 0;
		if (carvingStarts > 0)
			checkMaxStep = static_cast<int>(ml/carvingStarts); // integer division
		
		#ifdef DEBUG_MAXPOINTS
			cout << "checkMaxStep = " << checkMaxStep << endl;
		#endif
		
		//containers to store needed outputs
		vector<real> posteriorSupportVec;
		posteriorSupportVec.reserve(ml);
		vector<real> emptyBoxVec;
		emptyBoxVec.reserve(ml);
		
		//take a copy of the base histogram
		const AdaptiveHistogram adhBase(adh);
		
		size_t startLeaves = adh.getRootLeaves();
		
		{
		
			std::vector<size_t> carvedLaunchPoints;
			std::vector<size_t> maxPosteriorPoints;
			std::vector<real> maxPosteriors;
			
			/* the maxleaves for the carving queue is a guess, so if 
			 * we can't get our starts with this, try with the max leaves
			 * that we did manage to get to, and reset everything - but
			 * we do have to make sure that our steps are at least 1
			 * long or we won't find any start points at all. */	
			int tries = 2;
			bool successfulHist = false;
			while ( !successfulHist && (tries > 0) ) {
		
				bool shiftCatalan = false;
				successfulHist = adh.prioritySplitWithSupportPosteriorMaxLik(
						evaluatorCarving, 
						evaluatorSEB, 
						logPQ, 
						minPoints,
						minVol,
						logPrior, 
						emptyBoxVec,
						posteriorSupportVec,
						checkMaxStep,
						stopOnMaxPosterior,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors,
						rgsl, shiftCatalan);
				
				tries--;

				if (!successfulHist && tries > 0 && (carvingStarts > 0)) { // not the last try
				
					ml = posteriorSupportVec.size() + startLeaves - 1;
					evaluatorCarving.setMaxLeaves(ml);
					
					checkMaxStep = static_cast<int>(ml/carvingStarts); // integer division
			
					std::cout << "\n\nResetting carving queue: " 
						<< "PQ failed at " << ml << " leaves" << endl;
					
					if (checkMaxStep) {
						cout << "changed carving max leaves to "  << evaluatorCarving.getMaxLeaves() << " leaves" << endl;
						cout << "checkMaxStep = " << checkMaxStep << "\n" << endl;
					 
						adh = adhBase;
						gsl_rng_memcpy (rgsl, rgsl_base);
						emptyBoxVec.clear();
						posteriorSupportVec.clear();
						carvedLaunchPoints.clear();
						maxPosteriorPoints.clear();
						maxPosteriors.clear();
					}
					else { // checkMaxStep == 0
						tries = 0; // force failure
						cout << "Cannot reset: checkMaxStep = " 
						<< ml << "/" << carvingStarts << " = " 
						<< checkMaxStep << "\n" << endl;
					 
					}
				}
			}// end while 

			if (!successfulHist) {
				
				std::cout << "PQ failed: number of leaves is " 
						<< (posteriorSupportVec.size() + startLeaves - 1) << endl;
				gsl_rng_free (rgsl);
				gsl_rng_free (rgsl_base);
				throw std::runtime_error("Unsuccessful posterior support pq");
			}
			
			{
				cout << "Carved launch points" << endl;
				ostream_iterator<size_t> out_it (cout,"\t");
				copy ( carvedLaunchPoints.begin(), carvedLaunchPoints.end(), out_it );
				cout <<  endl;
			}
			{
				cout << "Max posterior points" << endl;
				ostream_iterator<size_t> out_it (cout,"\t");
				copy ( maxPosteriorPoints.begin(), maxPosteriorPoints.end(), out_it );
				cout <<  endl;
			}
			{
				cout << "maxPosteriors" << endl;
				ostream_iterator<real> out_it (cout,"\t");
				copy ( maxPosteriors.begin(),maxPosteriors.end(), out_it );
				cout <<  endl;
			}

			assert(carvedLaunchPoints.size() == maxPosteriorPoints.size());
			assert(carvedLaunchPoints.size() == maxPosteriors.size());
			cout << "\nnumber of max points is " << carvedLaunchPoints.size() << endl;
			
			/*chuck out the points that did not work */
			{
				std::vector < size_t> duplicateMaxPosteriorIndices;
				for (size_t i = 0; i < maxPosteriorPoints.size() ; ++i) {
					for (size_t j = i+1; j < maxPosteriorPoints.size() ; ++j) { 
						if (maxPosteriorPoints[i] == maxPosteriorPoints[j]) {
							// does not matter if same index gets in here more than once
							duplicateMaxPosteriorIndices.push_back(j);
						}
					}
					
				}
				
				std::vector<size_t> tmpCarvedLaunchPoints;
				std::vector<size_t> tmpMaxPosteriorPoints;
				std::vector<real> tmpMaxPosteriors;
				
				for (size_t i = 0; i < maxPosteriors.size() ; ++i) {
					
					real thisPost = maxPosteriors[i];
					
					if ( gsl_isinf(_double(thisPost)) ) {
						cout << "Throwing out maxPosterior at index " << i << endl;
					}
					else {
						
						// check i is not one of the duplicate indices
						std::vector < size_t>::iterator dit 
							= std::find(duplicateMaxPosteriorIndices.begin(),
								duplicateMaxPosteriorIndices.end(), i);
						
						if (!(dit < duplicateMaxPosteriorIndices.end()) ) {		
						
							tmpCarvedLaunchPoints.push_back(carvedLaunchPoints[i]);
							tmpMaxPosteriorPoints.push_back(maxPosteriorPoints[i]);
							tmpMaxPosteriors.push_back(maxPosteriors[i]);
						}// end if not duplicate
						else {// duplicate
							cout << "Ignoring max posterior " << thisPost
								<< " because it is a duplicate" << endl;
						}
					}
				}
				tmpCarvedLaunchPoints.swap(carvedLaunchPoints);
				tmpMaxPosteriorPoints.swap(maxPosteriorPoints);
				tmpMaxPosteriors.swap(maxPosteriors);
			}
			
			/* and save the best*/
			std::vector < size_t > bestOnes;  
			{
				std::vector < real > bestPosts;  
				
				for (size_t i = 0; i < maxPosteriors.size() ; ++i) {
					
					real thisPost = maxPosteriors[i];
					
					if (bestOnes.size() < keepBest) { // will fit in somewhere
							
							bestOnes.push_back(i);
							bestPosts.push_back(thisPost);
					}
					else if (keepBest) {
						
						std::vector<real>::iterator minBestPostIt 
							= min_element(bestPosts.begin(), bestPosts.end());
						
						if (*minBestPostIt < thisPost) {
							
							*minBestPostIt = thisPost;
							bestOnes.at(
								distance(bestPosts.begin(), minBestPostIt)) = i;
						}
					}
					
				}
				
			}
			
			cout << "After removing the ones with infinite posterior and the duplicates:" << endl;
			{
				cout << "Carved launch points" << endl;
				ostream_iterator<size_t> out_it (cout,"\t");
				copy ( carvedLaunchPoints.begin(), carvedLaunchPoints.end(), out_it );
				cout <<  endl;
			}
			{
				cout << "Max posterior points" << endl;
				ostream_iterator<size_t> out_it (cout,"\t");
				copy ( maxPosteriorPoints.begin(), maxPosteriorPoints.end(), out_it );
				cout <<  endl;
			}
			{
				cout << "maxPosteriors" << endl;
				ostream_iterator<real> out_it (cout,"\t");
				copy ( maxPosteriors.begin(),maxPosteriors.end(), out_it );
				cout <<  endl;
			}
			if (!maxPosteriors.empty()) {
				std::vector < real >::iterator maxMaxPosteriorIt = 
						max_element(maxPosteriors.begin(),maxPosteriors.end());
				
				real maxMaxPosterior = *maxMaxPosteriorIt;
				
				size_t maxMaxPosteriorIndex = static_cast<size_t> (distance(maxPosteriors.begin(), 
													maxMaxPosteriorIt));
				
				size_t maxMaxLaunchPoint = carvedLaunchPoints.at(maxMaxPosteriorIndex);
				size_t maxMaxPostPoint = maxPosteriorPoints.at(maxMaxPosteriorIndex);
				
				cout << "\nmaximium is the one at index " << maxMaxPosteriorIndex << endl;
				cout << "maxMaxLaunchPoint = " << maxMaxLaunchPoint << endl;
				cout << "maxMaxPostPoint = " << maxMaxPostPoint << endl;
				cout << "maxMaxPosterior = " << maxMaxPosterior << endl;
				
				for (size_t i = 0; i < carvedLaunchPoints.size() ; ++i) {
					
					// is this index in our best ones?
					bool keep = false;
					int keepIndex;
					vector<size_t>::iterator fit 
						= find (bestOnes.begin(), bestOnes.end(), i);
					if (fit != bestOnes.end()) {
						keep = true;
						keepIndex = static_cast<int>(
								distance(bestOnes.begin(), fit)) + 1;
					

						size_t carvedPt = carvedLaunchPoints[i];
						size_t maxPostPoint = maxPosteriorPoints[i];
						real maxPost = maxPosteriors[i];
						
						AdaptiveHistogram adhmax(adhBase);
						
						// prng in pristine state
						gsl_rng * rgsl_max = gsl_rng_clone(rgsl_base);
						
						bool successPQ = false;
						std::vector<real> thisPosteriorVec;
						std::vector<real> thisLoglikVec; 
						
						AdaptiveHistogram::PrioritySplitQueueEvaluator
								thisEvaluatorCarving(evaluatorCarving.getMeasurer(), carvedPt);
						
						AdaptiveHistogram::PrioritySplitQueueEvaluator
								thisEvaluatorSEB(evaluatorSEB.getMeasurer(), maxPostPoint);
						
						#ifdef DOSAMPLES
							logPQ = LOGSAMPLES;
						#endif
						
						bool shiftCatalan = false;
												
						successPQ = adhmax.prioritySplitMaxLik(
									thisEvaluatorCarving, 
									thisEvaluatorSEB, 
									logPQ, 
									minPoints,
									minVol,
									logPrior,
									thisPosteriorVec,
									thisLoglikVec,
									rgsl_max, shiftCatalan);
						
						if (successPQ) {
							
							size_t leaves = adhmax.getRootLeaves();
							assert(leaves == maxPostPoint);
							
							real maxLik= adhmax.getLogLikelihood();
							real lnPrior = logPrior(leaves-1);
							
							real thisPost = maxLik + lnPrior;
							
							cout << "\nAt max point with leaves " << maxPostPoint
								<< " maxPosterior from earlier pq was " << maxPost
								<< " and posterior calculated here is " << thisPost << endl;
							
							// if we want to keep this one
							if (keep) {
								hists.push_back(new AdaptiveHistogram(adhmax));
							}
							
							// output log-posterior
							if (!checkPostFileNameBase.empty()) {
								ostringstream oss;
								oss << checkPostFileNameBase << "_" << carvedPt << ".txt";
								string postFileName = oss.str();
								std::vector < const subpavings::RealVec* > dataPtrs;
								dataPtrs.push_back(&thisPosteriorVec);
								dataPtrs.push_back(&thisLoglikVec);
								std::vector < std::string > colNames;
								colNames.push_back("Post");
								colNames.push_back("lnLik");
								
								outputToFileVertical(dataPtrs, 
												colNames,
												postFileName);
							}	
							
						}
						else cout << "Could not recreate state for max point with leaves " 
							<< maxPostPoint << endl;
						
						gsl_rng_free(rgsl_max);
					}
					
				}
			}
			else cout << "No achievable maximum posterior points in that run" << endl;	
			// output log-support posteriors
			{
				std::vector < const subpavings::RealVec* > dataPtrs;
				dataPtrs.push_back(&emptyBoxVec);
				dataPtrs.push_back(&posteriorSupportVec);
				std::vector < std::string > colNames;
				colNames.push_back("Empty");
				colNames.push_back("PostSupp");
				
				outputToFileVertical(dataPtrs, 
								colNames,
								postFileName,
								prec);
			}

		}
		
		try {
			gsl_rng_free (rgsl);
			rgsl = NULL;
			gsl_rng_free (rgsl_base);
			rgsl_base = NULL;
		}
		catch(...) {}// catch and swallow
		

		return hists;

	} 


	//find start points furthest apart
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsMaxLeafSpread(
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
			unsigned long int seed)
	{
		double minVol = 0.0;
		size_t outOfClosest = static_cast<size_t>(carvingStarts)+1;
		return findStartingPointsMaxLeafSpread(
			adh,
			hists,
			evaluatorCarving, 
			evaluatorSEB, 
			logPrior,
			minPoints,
			minVol,
			carvingStarts,
			keep,
			outOfClosest,
			stopOnMaxPosterior,
			postFileName,
			checkPostFileNameBase,	
			prec,
			seed);
	}

	//find start points furthest apart
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsMaxLeafSpread(
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
			unsigned long int seed)
	{
		
		size_t outOfClosest = static_cast<size_t>(carvingStarts)+1;
		
		return findStartingPointsMaxLeafSpread(
			adh,
			hists,
			evaluatorCarving, 
			evaluatorSEB, 
			logPrior,
			minPoints,
			minVol,
			carvingStarts,
			keep,
			outOfClosest,
			stopOnMaxPosterior,
			postFileName,
			checkPostFileNameBase,	
			prec,
			seed);
	}

	//find start points furthest apart
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsMaxLeafSpread(
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
			unsigned long int seed)
	{
		double minVol = 0.0;
		return findStartingPointsMaxLeafSpread(
			adh,
			hists,
			evaluatorCarving, 
			evaluatorSEB, 
			logPrior,
			minPoints,
			minVol,
			carvingStarts,
			keep,
			outOfClosest,
			stopOnMaxPosterior,
			postFileName,
			checkPostFileNameBase,	
			prec,
			seed);
	}
	
	//find start points furthest apart
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsMaxLeafSpread(
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
			unsigned long int seed)
	{
		if (outOfClosest < keep) throw std::invalid_argument(
				"findStartingPointsMaxLeafSpread(...): outOfClosest < keep");	
		if (static_cast<int>(outOfClosest) > carvingStarts+1) throw std::invalid_argument(
				"findStartingPointsMaxLeafSpread(...): outOfClosest > carvingStarts+1");	
				
		gsl_rng * rgsl = NULL;
		gsl_rng * rgsl_base = NULL;
		
		try {
		
			rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
			gsl_rng_set (rgsl, seed); // change the seed
			
			// clone it
			rgsl_base = gsl_rng_clone(rgsl);

			std::vector<size_t> carvedLaunchPoints;
			std::vector<size_t> maxPosteriorPoints;
			std::vector<real> maxPosteriors;
			
			int tries = 3;
			
			while(tries && (carvedLaunchPoints.size() < keep)) {
				
				if (!carvedLaunchPoints.empty()) {
					std::vector<size_t>().swap(carvedLaunchPoints);
					std::vector<size_t>().swap(maxPosteriorPoints);
					std::vector<real>().swap(maxPosteriors);
				}
			
				findMaxPosteriorPoints(
													adh,
													carvedLaunchPoints,
													maxPosteriorPoints,
													maxPosteriors,
													evaluatorCarving, 
													evaluatorSEB, 
													logPrior,
													minPoints,
													minVol,
													carvingStarts,
													stopOnMaxPosterior,
													rgsl,
													postFileName,
													prec);
				--tries;
				
				#ifdef DEBUG_MAXPOINTS
					outputMaxPosteriorPoints(
						cout,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors);
				#endif
				
				/*chuck out the points that did not work */
				cleanMaxPoints(carvedLaunchPoints,
								maxPosteriorPoints,
								maxPosteriors);
				
				if (tries && (carvedLaunchPoints.size() < keep) ) {
					++carvingStarts;
					#ifdef DEBUG_MAXPOINTS
						cout << "Failed to find enough points: increasing choose starts" << endl;
					#endif
				}
			}
			try {
				gsl_rng_free (rgsl);
				rgsl = NULL;
			}
			catch(...) {}// catch and swallow


			if (carvedLaunchPoints.size() >= keep) {
				

				#ifdef DEBUG_MAXPOINTS
					cout << "After removing the ones with infinite posterior and the duplicates:" << endl;
					outputMaxPosteriorPoints(
						cout,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors);
				#endif

				std::vector < size_t > keepTheseOnes = keepFurthestApart(
							carvedLaunchPoints,
							maxPosteriorPoints,
							maxPosteriors,
							keep,
							outOfClosest);
				
				recreateStartingPoints(
						adh,
						hists, // filled in process
						evaluatorCarving, 
						evaluatorSEB, 
						logPrior, 
						minPoints,
						minVol,
						checkPostFileNameBase,	
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors,
						keepTheseOnes,
						rgsl_base);

			}
			
			else throw std::runtime_error(
				"findStartingPointsMaxLeafSpread(...): Not enough achievable starts");	
			try {
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			return hists;

		}
		catch(...) {
			try {
				if (rgsl != NULL) gsl_rng_free (rgsl);
				rgsl = NULL;
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			throw;
		}
	}
	 // end of finding start points
	 
	 
	 

	//find the best starting points
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsBest(
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
			unsigned long int seed)
	{
		double minVol = 0.0;
		
		return findStartingPointsBest(
					adh,
					hists,
					evaluatorCarving, 
					evaluatorSEB, 
					logPrior,
					minPoints,
					minVol,
					carvingStarts,
					keep,
					stopOnMaxPosterior,
					postFileName,
					checkPostFileNameBase,	
					prec,
					seed);
	
	}
	
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsBest(
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
			unsigned long int seed)
	{
		
		if (carvingStarts + 1< static_cast<int>(keep)) 
			throw std::invalid_argument(
						"findStartingPointsBest:\ncarvingStarts + 1 < keep");
			
		gsl_rng * rgsl = NULL;
		gsl_rng * rgsl_base = NULL;
		
		try {
		
			rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
			gsl_rng_set (rgsl, seed); // change the seed
			
			// clone it
			rgsl_base = gsl_rng_clone(rgsl);

			
			std::vector<size_t> carvedLaunchPoints;
			std::vector<size_t> maxPosteriorPoints;
			std::vector<real> maxPosteriors;
			
			int tries = 3;
			
			while(tries && (carvedLaunchPoints.size() < keep)) {
				
				if (!carvedLaunchPoints.empty()) {
					std::vector<size_t>().swap(carvedLaunchPoints);
					std::vector<size_t>().swap(maxPosteriorPoints);
					std::vector<real>().swap(maxPosteriors);
				}
			
				findMaxPosteriorPoints(
												adh,
												carvedLaunchPoints,
												maxPosteriorPoints,
												maxPosteriors,
												evaluatorCarving, 
												evaluatorSEB, 
												logPrior,
												minPoints,
												minVol,
												carvingStarts,
												stopOnMaxPosterior,
												rgsl,
												postFileName,
												prec);
				--tries;
				
				#ifdef DEBUG_MAXPOINTS
					outputMaxPosteriorPoints(
						cout,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors);
				#endif
				
				/*chuck out the points that did not work */
				cleanMaxPoints(carvedLaunchPoints,
								maxPosteriorPoints,
								maxPosteriors);
				
				if (tries && (carvedLaunchPoints.size() < keep) ) {
					++carvingStarts;
					#ifdef DEBUG_MAXPOINTS
						cout << "Failed to find enough points: increasing choose starts" << endl;
					#endif
				}
			}
			try {
				gsl_rng_free (rgsl);
				rgsl = NULL;
			}
			catch(...) {}// catch and swallow

			if (carvedLaunchPoints.size() >= keep) {
				

				#ifdef DEBUG_MAXPOINTS
					cout << "After removing the ones with infinite posterior and the duplicates:" << endl;
					outputMaxPosteriorPoints(
						cout,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors);
				#endif

				std::vector < size_t > keepTheseOnes = keepBest(
							carvedLaunchPoints,
							maxPosteriorPoints,
							maxPosteriors,
							keep);

				recreateStartingPoints(
						adh,
						hists, // filled in process
						evaluatorCarving, 
						evaluatorSEB, 
						logPrior, 
						minPoints,
						minVol,
						checkPostFileNameBase,	
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors,
						keepTheseOnes,
						rgsl_base);

			}
			
			else throw std::runtime_error(
				"findStartingPointsMaxLeafSpread(...): Not enough achievable starts");	
			try {
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			return hists;

		}
		catch(...) {
			try {
				if (rgsl != NULL) gsl_rng_free (rgsl);
				rgsl = NULL;
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			throw;
		}
	}
	 // end of finding start points

	//find start points overdispersed
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsOverdispersed(
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
			unsigned long int seed)
	{
		double minVol = 0.0;
		double percentSpread = 0.95;
		
		
		return findStartingPointsOverdispersed(
				adh,
				hists,
				evaluatorCarving, 
				evaluatorSEB, 
				logPrior,
				minPoints,
				minVol,
				carvingStarts,
				keep,
				percentSpread,
				stopOnMaxPosterior, //always
				postFileName,
				checkPostFileNameBase,	
				prec,
				seed);
	}


	//find start points overdispersed
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsOverdispersed(
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
			unsigned long int seed)
	{
		double percentSpread = 0.95;
		
		
		return findStartingPointsOverdispersed(
				adh,
				hists,
				evaluatorCarving, 
				evaluatorSEB, 
				logPrior,
				minPoints,
				minVol,
				carvingStarts,
				keep,
				percentSpread,
				stopOnMaxPosterior, //always
				postFileName,
				checkPostFileNameBase,	
				prec,
				seed);
	}


	//find start points overdispersed
	/* Tries to increase starts if it cannot find any maxes (this is
	 * redundant with current implementation in Adaptive Hists, where it 
	 * puts in the point reached at end of SEB if no max is found) and
	 * also increases carving leaves, if it can (restricted by max leaves on seb)
	 * if the max found is at root and it is doing  carving starts*/
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsOverdispersed(
			const AdaptiveHistogram& adh,
			std::vector< AdaptiveHistogram* >& hists,
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
			LogMCMCPrior& logPrior,
			size_t minPoints,
			int carvingStarts,
			size_t keep,
			double percentSpread,
			bool stopOnMaxPosterior, //always
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,	
			int prec,
			unsigned long int seed)
	{
		
		double minVol = 0.0;
		return findStartingPointsOverdispersed(
							adh,
							hists,
							evaluatorCarving, 
							evaluatorSEB, 
							logPrior,
							minPoints,
							minVol,
							carvingStarts,
							keep,
							percentSpread,
							stopOnMaxPosterior, //always
							postFileName,
							checkPostFileNameBase,	
							prec,
							seed);
	}

	
	std::vector< AdaptiveHistogram* >& CarverSEB::findStartingPointsOverdispersed(
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
			bool stopOnMaxPosterior, //always
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,	
			int prec,
			unsigned long int seed)
	{
		#ifdef DEBUG_MAXPOINTS
			cout << "\nUsing findStartingPointsOverdispersed, keep = " << keep 
				<< " and percentSpread = " << percentSpread << endl;
		#endif
				
		if (keep < 1) {
			return hists;
		}
			
		if (!(percentSpread < 1.0))
			throw std::invalid_argument(
				"findStartingPointsOverdispersed(...) : percentSpread >= 1.0");
		if (!(percentSpread > 0.0))
			throw std::invalid_argument(
				"findStartingPointsOverdispersed(...) : percentSpread <= 0.0");
			
		gsl_rng * rgsl = NULL;
		gsl_rng * rgsl_base = NULL;
		
		
		try {
		
			rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
			gsl_rng_set (rgsl, seed); // change the seed
			
			// clone it
			rgsl_base = gsl_rng_clone(rgsl);

			std::vector<size_t> carvedLaunchPoints;
			std::vector<size_t> maxPosteriorPoints;
			std::vector<real> maxPosteriors;
			
			int tries = 3;
			
			int indexMax = 0;
			
			size_t carvingLeaves = evaluatorCarving.getMaxLeaves();
			size_t sebLeaves = evaluatorSEB.getMaxLeaves();
				
			
			while(tries && (carvedLaunchPoints.empty() || 
				(( !indexMax || (indexMax == carvingStarts) ) && (carvingStarts > 0)) ) ) {
					
				
				if (!carvedLaunchPoints.empty()) {
					std::vector<size_t>().swap(carvedLaunchPoints);
					std::vector<size_t>().swap(maxPosteriorPoints);
					std::vector<real>().swap(maxPosteriors);
				}
			
				carvingLeaves = findMaxPosteriorPoints(
											adh,
											carvedLaunchPoints,
											maxPosteriorPoints,
											maxPosteriors,
											evaluatorCarving, 
											evaluatorSEB, 
											logPrior,
											minPoints,
											minVol,
											carvingStarts,
											stopOnMaxPosterior,
											rgsl,
											postFileName,
											prec);
				--tries;
				
				#ifdef DEBUG_MAXPOINTS
					outputMaxPosteriorPoints(
						cout,
						carvedLaunchPoints,
						maxPosteriorPoints,
						maxPosteriors);
				#endif
				
				/*chuck out the points that did not work */
				cleanMaxPoints(carvedLaunchPoints,
								maxPosteriorPoints,
								maxPosteriors);
				
				/* check if we have no maxes left at all -
				 * under present implementation this should never be a problem
				 * 
				 * could alter implementation in adh back to putting in infinity
				 * and then increase the seb queue length but I am inclined not 
				 * to do this 
				 */ 
				if (tries && carvedLaunchPoints.empty() ) {
					++carvingStarts;
					#ifdef DEBUG_MAXPOINTS
						cout << "\nFailed to find a maximum: increasing choose starts" << endl;
					#endif
				}
				else {
					indexMax = findMax(maxPosteriors);
					
					/* if carving starts is 0 this is okay
					 * but otherwise if max is first or last one 
					 * and we have tries left 
					 * then we will try to keep looping */
					if ( tries && (carvingStarts > 0) 
								&& ( !indexMax || (indexMax == carvingStarts) ) 
								&& (carvingLeaves*2 < sebLeaves) ) {
						carvingLeaves *=2;
						evaluatorCarving.setMaxLeaves(carvingLeaves);
						#ifdef DEBUG_MAXPOINTS
						
							cout << "\nMax is launched from root or last start point:"
								<< " increasing carving leaves to "
								<< carvingLeaves << endl;
						#endif
					} /* but set tries to 0 if we can't increase carving leaves */
					else if (!(carvingLeaves*2 < sebLeaves)) tries = 0;
					
				}
			} // end while loop
			try {
				gsl_rng_free (rgsl);
				rgsl = NULL;
			}
			catch(...) {}// catch and swallow

			/* need at least one to find a max */
			if (carvedLaunchPoints.empty()) 
				throw std::runtime_error(
					"findStartingPointsOverdispersed(...) : could not get max");
				

			#ifdef DEBUG_MAXPOINTS
				cout << "\nAfter removing the ones with infinite posterior and the duplicates:" << endl;
				outputMaxPosteriorPoints(
					cout,
					carvedLaunchPoints,
					maxPosteriorPoints,
					maxPosteriors);
			#endif

			
			LOGGING_LEVEL logPQ = NOLOG; // logging level  
			#ifdef DOSAMPLES
				logPQ = LOGSAMPLES;
			#endif

			size_t carvedPt = carvedLaunchPoints[indexMax];
			size_t maxPostPoint = maxPosteriorPoints[indexMax];
			real maxPost = maxPosteriors[indexMax];
			
			#ifdef DEBUG_MAXPOINTS
				cout << "\nMax is at index " << indexMax 
				<< ", carving point " << carvedPt
				<< ", maxPostPoint " << maxPostPoint
				<< ", maxPost " << maxPost << endl;
				
			#endif
			
			std::vector<real> thisPosteriorVec; 
			std::vector<real> thisLoglikVec; 
			
			AdaptiveHistogram* adhmaxPtr = NULL;
			
			try {
				
				adhmaxPtr = new AdaptiveHistogram(adh);
				
				bool successPQ = recreateStartingPoint(
						*adhmaxPtr,
						evaluatorCarving, 
						evaluatorSEB, 
						logPrior, 
						minPoints,
						minVol,
						carvedPt,
						maxPostPoint,
						thisPosteriorVec, 
						thisLoglikVec, 
						rgsl_base,
						logPQ);
				
				if (!successPQ) {
					cerr << "Could not recreate state for max point with leaves "; 
					throw std::runtime_error(
						"findStartingPointsOverdispersed(...) : could not recreate max");
				}
				if (thisPosteriorVec.empty()) // should never happen
					throw std::runtime_error(
						"findStartingPointsOverdispersed(...) : no posteriors from max");
			}
			catch (...) {
				if (adhmaxPtr != NULL)  delete adhmaxPtr;
				adhmaxPtr = NULL;
				
				throw;
			}	
				
			#ifdef DEBUG_MAXPOINTS
				size_t leaves = adhmaxPtr->getRootLeaves();
				assert(leaves == maxPostPoint);
				
				real maxLik= adhmaxPtr->getLogLikelihood();
				real lnPrior = logPrior(leaves-1);
				
				real thisPost = maxLik + lnPrior;
				
				cout << "\nOverall max point with leaves " << maxPostPoint
						<< " comparison maxPosterior is " << maxPost
						<< " and posterior calculated here is " << thisPost << endl;
				cout << "(logLik here is " << maxLik << " and logPrior is " 
											<< lnPrior << ")" << endl;
			#endif
			
			// output log-posterior for max
			if (!checkPostFileNameBase.empty()) {
				ostringstream oss;
				oss << checkPostFileNameBase << "_" << carvedPt 
									<< "_" << maxPostPoint << ".txt";
				string checkPostFileName = oss.str();
				std::vector < const subpavings::RealVec* > dataPtrs;
				dataPtrs.push_back(&thisPosteriorVec);
				dataPtrs.push_back(&thisLoglikVec);
				std::vector < std::string > colNames;
				colNames.push_back("Post");
				colNames.push_back("lnLik");
				
				outputToFileVertical(dataPtrs, 
								colNames,
								checkPostFileName);
			}
			
			
			/*want to fit in ceil (keep-1)/2 points before the max*/
			size_t toFindBefore = ceil((keep-1)/2.0);
			size_t numPosts = thisPosteriorVec.size();
				
			size_t startLeaves = adh.getRootLeaves();
			size_t indexPercentPost = 0;	
			double dblGap = 1.0;
			if (keep > 1) {
				
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "\nsorting out start points from before" << endl;
					
				#endif
			
				real startPosterior = thisPosteriorVec.at(0);
				real endPosterior = thisPosteriorVec.back();
				
				real gainPost = endPosterior - startPosterior;
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "numPosts = " << numPosts << endl;
					cout << "startPosterior = " << startPosterior << endl;
					cout << "endPosterior = " << endPosterior << endl;
					cout << "gainPost = " << gainPost << endl;
					
				#endif
				
				if (gainPost > 0.0) {
				
					real gainPercentPost = startPosterior + percentSpread*gainPost;
					
					// find 95% gain point
					for (size_t i = 0; i < numPosts; ++i) {
						if (thisPosteriorVec[i] > gainPercentPost) {
							indexPercentPost = i;
							break; // break out of for loop
						}
					}
					#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
						cout << "gainPercentPost = " << gainPercentPost << endl;
						
					#endif
					
				}
				//else indexPercentPost = 0;
				
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "indexPercentPost = " << indexPercentPost << endl;
					cout << "startLeaves = " << startLeaves << endl;
					
				#endif
				
				/* check if we can get what we need before */
				if(numPosts - 1 < toFindBefore) toFindBefore = numPosts-1;
				
				if (indexPercentPost > numPosts - toFindBefore - 1)
							indexPercentPost = numPosts - toFindBefore - 1;
		
				size_t toChooseFrom = numPosts - indexPercentPost - 1;
				/* so we know toChooseFrom >= toFindBefore */
				
				if (toFindBefore > 0) 
					dblGap = toChooseFrom/toFindBefore;
				
				//if we have to do everything after
				else if (!toFindBefore && numPosts > keep) {
					dblGap = (numPosts-1.0)/(keep-1.0);
				}
				/* if !toFindBefore && numPosts <= keep dblGap = 1.0 by default */
					
				assert(!(dblGap < 1.0));
				
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "toFindBefore = " << toFindBefore << " and toChooseFrom = " << toChooseFrom << endl;
					cout << "dblGap = " << dblGap << endl;
						
				#endif
			}
			if (toFindBefore) {
				std::vector< AdaptiveHistogram* > histsBefore;
				
				addToHists(
					adh,
					histsBefore,
					evaluatorCarving, 
					evaluatorSEB, 
					logPrior,
					minPoints,
					minVol,
					toFindBefore,
					dblGap,
					indexPercentPost,
					thisPosteriorVec,
					startLeaves,
					carvedPt,
					checkPostFileNameBase,	
					rgsl_base);
					
				assert(histsBefore.size() == toFindBefore);
					
				hists.insert(hists.end(), histsBefore.begin(), histsBefore.end());
			} // end adding from before max
			
			// add the max point
			hists.push_back(adhmaxPtr);
			
			size_t toFindAfter = keep - 1 -toFindBefore;
				
			if (toFindAfter > 0) {
				
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "\nsorting out start points from after" << endl;
					
				#endif
			
				size_t evaluatorMaxLeaves = evaluatorSEB.getMaxLeaves();
				bool evaluatorUsingCritStop = evaluatorSEB.getUsingCritStop();
				
				size_t toAdd = static_cast<size_t>(toFindAfter*dblGap);
				if (toAdd) {
					/* we should reset the evaluatorSEB*/
					evaluatorSEB.setMaxLeaves(
						maxPostPoint + static_cast<size_t>(toFindAfter*dblGap));
					evaluatorSEB.setUsingCritStop(false);
					
					// and stuff the end of thisPosterorVec
					std::vector<real> tmp(toAdd, thisPosteriorVec.back());
					thisPosteriorVec.insert(thisPosteriorVec.end(), tmp.begin(), tmp.end());
				}
				
				size_t baseIndex = numPosts-1 + static_cast<size_t>(dblGap);
				
				#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
					cout << "indexMax = " << indexMax << " and toFindAfter = " << toFindAfter 
					<< " and numPosts = " << numPosts
					<< " and toAdd = " << toAdd << endl;
						
				#endif
					
				
				std::vector< AdaptiveHistogram* > histsAfter;
				addToHists(
					adh,
					histsAfter,
					evaluatorCarving, 
					evaluatorSEB, 
					logPrior,
					minPoints,
					minVol,
					toFindAfter,
					dblGap,
					baseIndex,
					thisPosteriorVec,
					startLeaves,
					carvedPt,
					checkPostFileNameBase,	
					rgsl_base);
					
				
				assert(histsAfter.size() == toFindAfter);
				
				hists.insert(hists.end(), histsAfter.begin(), histsAfter.end());
				
				if (toAdd > 0) {
					/* we should reset the evaluatorSEB*/
					evaluatorSEB.setMaxLeaves(evaluatorMaxLeaves);
					evaluatorSEB.setUsingCritStop(evaluatorUsingCritStop);
					
				}
			} // end adding from after max
			
			
			try {
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			assert(hists.size() == keep);
			
			#ifdef DEBUG_MAXPOINTS
				cout << "\nMade hist starting points, size " << hists.size() << "\n" << endl;
				cout << "Final carving starts is " << carvingStarts 
						<< " and carving leaves is " 
						<< evaluatorCarving.getMaxLeaves() << endl; 
			#endif
			
			return hists;

		}
		catch(...) {
			try {
				if (rgsl != NULL) gsl_rng_free (rgsl);
				rgsl = NULL;
				if (rgsl_base != NULL) gsl_rng_free (rgsl_base);
				rgsl_base = NULL;
			}
			catch(...) {}// catch and swallow
			
			throw;
		}
	}
	 // end of finding start points
	 
	/* add a hist ptr to hists.
	Note: the evaluator SEB is not changed in this routine so - the
			 * routine will use whatever criteria it contains - so if you
			 * definitely want to get to the max leaves specified by
			 * carver evaluator and seb evaluator, make sure that the
			 * critstops for these are 0 (ie only max leaves will stop the process) */ 
	std::vector< AdaptiveHistogram* >& CarverSEB::createStartingPoint(
			const AdaptiveHistogram& adhBase,
			std::vector< AdaptiveHistogram* >& hists,
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			const std::string& checkPostFileNameBase,	
			unsigned long int seed)
	{
		double minVol = 0.0;
		
		return createStartingPoint(
				adhBase,
				hists,
				evaluatorCarving, 
				evaluatorSEB, 
				logPrior, 
				minPoints,
				minVol,
				checkPostFileNameBase,	
				seed);
	
	}
	
	std::vector< AdaptiveHistogram* >& CarverSEB::createStartingPoint(
			const AdaptiveHistogram& adhBase,
			std::vector< AdaptiveHistogram* >& hists,
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			double minVol,
			const std::string& checkPostFileNameBase,	
			unsigned long int seed)
	{
		gsl_rng * rgsl = NULL;
			
		try {
		
			rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
			gsl_rng_set (rgsl, seed); // change the seed
			
			LOGGING_LEVEL logPQ = NOLOG; // logging level  
			#ifdef DOSAMPLES
				logPQ = LOGSAMPLES;
			#endif
			
			std::vector<real> thisPosteriorVec; 
			std::vector<real> thisLoglikVec; 
			
			bool successPQ = false;
			AdaptiveHistogram* adhPtr;
			
			try {
				
				adhPtr = new AdaptiveHistogram(adhBase);

				
				successPQ = recreateStartingPoint(
						*adhPtr,
						evaluatorCarving, 
						evaluatorSEB, 
						logPrior, 
						minPoints,
						minVol,
						thisPosteriorVec, 
						thisLoglikVec, 
						rgsl,
						logPQ);
			}
			catch (...) {
				if (adhPtr != NULL) delete adhPtr;
				adhPtr = NULL;
				
				throw;
			}	
			if (successPQ) {
				
				hists.push_back(adhPtr);
				
				// output log-posterior
				if (!checkPostFileNameBase.empty()) {
					ostringstream oss;
					oss << checkPostFileNameBase << "_" << evaluatorCarving.getMaxLeaves() 
								<< "_" << evaluatorSEB.getMaxLeaves()  << ".txt";
					string checkPostFileName = oss.str();
					
					std::vector < const subpavings::RealVec* > dataPtrs;
					dataPtrs.push_back(&thisPosteriorVec);
					dataPtrs.push_back(&thisLoglikVec);
					std::vector < std::string > colNames;
					colNames.push_back("Post");
					colNames.push_back("lnLik");
					
					outputToFileVertical(dataPtrs, 
									colNames,
									checkPostFileName);
				}	
				
			}
			else {
				cerr << "Could not create state for max point with leaves " 
				<< evaluatorSEB.getMaxLeaves() << endl;
				if (adhPtr != NULL) delete adhPtr;
				adhPtr = NULL;
			}
			
			return hists;
		}
		catch(...) {
			try {
				if (rgsl != NULL) gsl_rng_free (rgsl);
				rgsl = NULL;
			}
			catch(...) {}// catch and swallow
				
			throw;
		}

		
	}

	AdaptiveHistogram& CarverSEB::carve(
						AdaptiveHistogram& adh,
						AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving,  
						LogMCMCPrior& logPrior, 
						size_t minPoints,
						const std::string& postFileName,
						int prec,
						unsigned long int seed)
	{
		double minVol = 0.0;
		return carve(
						adh,
						evaluatorCarving,  
						logPrior,
						minPoints,
						minVol,
						postFileName,
						prec,
						seed);
	
	}
	
	AdaptiveHistogram& CarverSEB::carve(
						AdaptiveHistogram& adh,
						AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving,  
						LogMCMCPrior& logPrior, 
						size_t minPoints,
						double minVol,
						const std::string& postFileName,
						int prec,
						unsigned long int seed)
	{
		
		gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937); // set up with default seed
		gsl_rng_set (rgsl, seed); // change the seed
		
		// parameters needed to start the PQ
		LOGGING_LEVEL logPQ = NOLOG; // logging level  
		
		//containers to store needed outputs
		vector<real> posteriorVec;
		vector<real> posteriorSupportVec;
		vector<real> emptyBoxVec;
		
		bool shiftCatalan = false;
		
		bool successfulHist = adh.prioritySplitWithSupportPosterior(
				evaluatorCarving, 
				logPQ, 
				minPoints,
				minVol,
				logPrior, 
				emptyBoxVec,
				posteriorSupportVec,
				posteriorVec,
				rgsl, shiftCatalan);

		if (!successfulHist) {
			gsl_rng_free (rgsl);
			throw std::runtime_error("Unsuccessful posterior support pq");
		}
		
			
		// output
		if (!postFileName.empty() )	{
			std::vector < const subpavings::RealVec* > dataPtrs;
			dataPtrs.push_back(&emptyBoxVec);
			dataPtrs.push_back(&posteriorSupportVec);
			dataPtrs.push_back(&posteriorVec);
			std::vector < std::string > colNames;
			colNames.push_back("Empty");
			colNames.push_back("PostSupp");
			colNames.push_back("Posterior");
			
			outputToFileVertical(dataPtrs, 
							colNames,
							postFileName,
							prec);
		}

		try {
			gsl_rng_free (rgsl);
			rgsl = NULL;
		}
		catch(...) {}// catch and swallow

		return adh; // return the original histogram

	} 
	 

	// -------------------- private -------------------------
	 
	 
	/*find the maximum posterior point */ 
	size_t CarverSEB::findMax(const std::vector<real>& maxPosteriors)
	{
		
		if (maxPosteriors.empty())
			throw std::invalid_argument("findMax(const std::vector<real>&) : maxes empty");
		
		size_t numMaxPts = maxPosteriors.size();
		
		/* find the first max point */
		size_t indexMax = 0;
		for (size_t i = 0; i < numMaxPts; ++i) {
			if (maxPosteriors[i] > maxPosteriors[indexMax]) {
				indexMax = i;
			}
		}

		return indexMax;
	}
	
	/*find start the max posterior points
	* return the final max leaves in the carving queue */
	size_t CarverSEB::findMaxPosteriorPoints(
			const AdaptiveHistogram& adh,
			std::vector<size_t>& carvedLaunchPoints,
			std::vector<size_t>& maxPosteriorPoints,
			std::vector<real>& maxPosteriors,
			AdaptiveHistogram::PrioritySplitQueueEvaluator evaluatorCarving, 
			const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			double minVol,
			int carvingStarts, // how many carving queue starts to do, ie in addition to root
			bool stopOnMaxPosterior,
			const gsl_rng * rgsl_base,
			const std::string& postFileName,
			int prec	)
	{
		
		// clone rgsl
		gsl_rng * rgsl = gsl_rng_clone(rgsl_base);
		
		size_t ml = evaluatorCarving.getMaxLeaves();
		int checkMaxStep = 0;
		if (carvingStarts > 0)
			checkMaxStep = static_cast<int>(ml/carvingStarts); // integer division
		
		#ifdef DEBUG_MAXPOINTS
			cout << "\n checkMaxStep = " << checkMaxStep << endl;
		#endif
		
		// parameters needed to start the PQ
		LOGGING_LEVEL logPQ = NOLOG; // logging level  
		
		//containers to store needed outputs
		vector<real> posteriorSupportVec;
		posteriorSupportVec.reserve(ml);
		vector<real> emptyBoxVec;
		emptyBoxVec.reserve(ml);
		
		// hist to use use
		AdaptiveHistogram adhUse(adh);
		
		size_t startLeaves = adhUse.getRootLeaves();
		
		/* the maxleaves for the carving queue is a guess, so if 
		 * we can't get our starts with this, try with the max leaves
		 * that we did manage to get to, and reset everything - but
		 * we do have to make sure that our steps are at least 1
		 * long or we won't find any start points at all. */	
		int tries = 2;
		bool successfulHist = false;
		while ( !successfulHist && (tries > 0) ) {
			
			bool shiftCatalan = false;
			
			successfulHist = adhUse.prioritySplitWithSupportPosteriorMaxLik(
					evaluatorCarving, 
					evaluatorSEB, 
					logPQ, 
					minPoints,
					minVol,
					logPrior, 
					emptyBoxVec,
					posteriorSupportVec,
					checkMaxStep,
					stopOnMaxPosterior,
					carvedLaunchPoints,
					maxPosteriorPoints,
					maxPosteriors,
					rgsl, shiftCatalan);
			
			tries--;

			if (!successfulHist && (tries > 0) && (carvingStarts > 0)) { // not the last try
			
				ml = posteriorSupportVec.size() + startLeaves - 1;
				evaluatorCarving.setMaxLeaves(ml);
				
				checkMaxStep = static_cast<int>(ml/carvingStarts); // integer division
		
				std::cout << "\n\nResetting carving queue: " 
					<< "PQ failed at " << ml << " leaves" << endl;
				
				if (checkMaxStep) {
					cout << "changed carving max leaves to "  << evaluatorCarving.getMaxLeaves() << " leaves" << endl;
					cout << "checkMaxStep = " << checkMaxStep << "\n" << endl;
				 
					adhUse = adh;
					gsl_rng_memcpy (rgsl, rgsl_base);
					emptyBoxVec.clear();
					posteriorSupportVec.clear();
					carvedLaunchPoints.clear();
					maxPosteriorPoints.clear();
					maxPosteriors.clear();
				}
				else { // checkMaxStep == 0
					tries = 0; // force failure
					cout << "Cannot reset: checkMaxStep = " 
					<< ml << "/" << carvingStarts << " = " 
					<< checkMaxStep << "\n" << endl;
				 
				}
			}
		}// end while 
		
		gsl_rng_free(rgsl);
		
		if (!successfulHist) {
		
			std::cout << "PQ failed: number of leaves is " 
					<< (posteriorSupportVec.size() + startLeaves - 1) << endl;
			throw std::runtime_error("Unsuccessful posterior support pq");
		}
		
		assert(carvedLaunchPoints.size() == maxPosteriorPoints.size());
		assert(carvedLaunchPoints.size() == maxPosteriors.size());
		
		// output log-support posteriors
		if (!postFileName.empty() )	{
			std::vector < const subpavings::RealVec* > dataPtrs;
			dataPtrs.push_back(&emptyBoxVec);
			dataPtrs.push_back(&posteriorSupportVec);
			std::vector < std::string > colNames;
			colNames.push_back("Empty");
			colNames.push_back("PostSupp");
			
			outputToFileVertical(dataPtrs, 
							colNames,
							postFileName,
							prec);
		}
		
		return ml;

		
	}
			
			
	//output the max points
	void CarverSEB::outputMaxPosteriorPoints(
			std::ostream& os,
			const std::vector<size_t>& carvedLaunchPoints,
			const std::vector<size_t>& maxPosteriorPoints,
			const std::vector<real>& maxPosteriors)
	{
		cout << "\nnumber of max points is " << carvedLaunchPoints.size() << endl;
		
		{
			os << "Carved launch points" << endl;
			ostream_iterator<size_t> out_it (os,"\t");
			copy ( carvedLaunchPoints.begin(), carvedLaunchPoints.end(), out_it );
			os <<  endl;
		}
		{
			os << "Max posterior points" << endl;
			ostream_iterator<size_t> out_it (os,"\t");
			copy ( maxPosteriorPoints.begin(), maxPosteriorPoints.end(), out_it );
			os <<  endl;
		}
		{
			os << "maxPosteriors" << endl;
			ostream_iterator<real> out_it (os,"\t");
			copy ( maxPosteriors.begin(),maxPosteriors.end(), out_it );
			os <<  endl;
		}
	}

	bool CarverSEB::essentiallyEqual(cxsc::real r1, cxsc::real r2)
	{
		if (gsl_isinf(_double(r1)) && gsl_isinf(_double(r2)) ) return true;
		else if (gsl_isinf(_double(r1)) || gsl_isinf(_double(r2))) return false;
		else {
			real largest = r1;
			real smallest = r2;
			
			if (r2 > r1) {
				largest = r2;
				smallest = r1;
			}
			
			return (largest < succ(succ(smallest)));
		}
	}

	/*chuck out the points that did not work */
	void CarverSEB::cleanMaxPoints(
			std::vector<size_t>& carvedLaunchPoints,
			std::vector<size_t>& maxPosteriorPoints,
			std::vector<real>& maxPosteriors)
	{
		size_t numMaxPts = carvedLaunchPoints.size();
		
		std::vector < size_t> duplicateMaxPosteriorIndices;
		for (size_t i = 0; i < numMaxPts ; ++i) {
			for (size_t j = i+1; j < numMaxPts ; ++j) { 
				if ( (maxPosteriorPoints[i] == maxPosteriorPoints[j])
				&&
				(essentiallyEqual(maxPosteriors[i], maxPosteriors[j]))
				) {
					// does not matter if same index gets in here more than once
					duplicateMaxPosteriorIndices.push_back(j);
				}
			}
			
		}
		std::vector<size_t> tmpCarvedLaunchPoints;
		std::vector<size_t> tmpMaxPosteriorPoints;
		std::vector<real> tmpMaxPosteriors;
		
		for (size_t i = 0; i < numMaxPts ; ++i) {
			
			real thisPost = maxPosteriors[i];
			
			if ( gsl_isinf(_double(thisPost)) ) {
				#ifdef DEBUG_MAXPOINTS
					cout << "Throwing out maxPosterior at index " << i << endl;
				#endif
			}
			else {
				
				// check i is not one of the duplicate indices
				std::vector < size_t>::iterator dit 
					= std::find(duplicateMaxPosteriorIndices.begin(),
						duplicateMaxPosteriorIndices.end(), i);
				
				if (!(dit < duplicateMaxPosteriorIndices.end()) ) {		
				
					tmpCarvedLaunchPoints.push_back(carvedLaunchPoints[i]);
					tmpMaxPosteriorPoints.push_back(maxPosteriorPoints[i]);
					tmpMaxPosteriors.push_back(maxPosteriors[i]);
				}// end if not duplicate
				else {// duplicate
					#ifdef DEBUG_MAXPOINTS
						cout << "Ignoring max posterior " << thisPost
							<< " because it is a duplicate" << endl;
					#endif
				}
			}
		}
		tmpCarvedLaunchPoints.swap(carvedLaunchPoints);
		tmpMaxPosteriorPoints.swap(maxPosteriorPoints);
		tmpMaxPosteriors.swap(maxPosteriors);

	}

	std::vector < size_t > CarverSEB::keepBest(
					const std::vector<size_t>& carvedLaunchPoints,
					const std::vector<size_t>& maxPosteriorPoints,
					const std::vector<real>& maxPosteriors,
					size_t keep)
	{
		
		size_t numMaxPts = carvedLaunchPoints.size();
		
		if (numMaxPts < keep)
			throw std::runtime_error(
				"keepBest(...): not enough points to choose from");
		
		std::vector < size_t > bestOnes;  
		
		if (numMaxPts) {
			/* find the max point */
			size_t indexMax = 0;
			for (size_t i = 0; i < numMaxPts; ++i) {
				if (maxPosteriors[i] > maxPosteriors[indexMax]) {
					indexMax = i;
				}
			}
			
			size_t maxMaxPostPoint = maxPosteriorPoints[indexMax];
			
			#ifdef DEBUG_MAXPOINTS
				real maxMaxPosterior = maxPosteriors[indexMax];
					
				size_t maxMaxLaunchPoint = carvedLaunchPoints[indexMax];
				
				cout << "\nmaximium is the one at index " << indexMax << endl;
				cout << "maxMaxLaunchPoint = " << maxMaxLaunchPoint << endl;
				cout << "maxMaxPostPoint = " << maxMaxPostPoint << endl;
				cout << "maxMaxPosterior = " << maxMaxPosterior << endl;
			#endif
		
		
			std::vector < real > bestPosts;  
			
			for (size_t i = 0; i < numMaxPts ; ++i) {
				
				real thisPost = maxPosteriors[i];
				
				if (bestOnes.size() < keep) { // will fit in somewhere
						
						bestOnes.push_back(i);
						bestPosts.push_back(thisPost);
				}
				else if (keep) {
					
					std::vector<real>::iterator minBestPostIt 
						= min_element(bestPosts.begin(), bestPosts.end());
					
					if (*minBestPostIt < thisPost) {
						
						*minBestPostIt = thisPost;
						bestOnes.at(
							distance(bestPosts.begin(), minBestPostIt)) = i;
					}
				}
				
			}
			
		}
		return bestOnes;
		
	}


	std::vector < size_t > CarverSEB::keepFurthestApart(
					const std::vector<size_t>& carvedLaunchPoints,
					const std::vector<size_t>& maxPosteriorPoints,
					const std::vector<real>& maxPosteriors,
					size_t keep,
					size_t outOfClosest)
	{
		size_t numMaxPts = carvedLaunchPoints.size();
		
		if (numMaxPts < outOfClosest)
			outOfClosest = numMaxPts;
		
		if (outOfClosest < keep)
			throw std::runtime_error(
				"keepFurthestApart(...): not enough points to choose from");
		
		std::vector < size_t > keepTheseOnes;  
		
		if (numMaxPts) {
			/* find the max point */
			size_t indexMax = 0;
			for (size_t i = 0; i < numMaxPts; ++i) {
				if (maxPosteriors[i] > maxPosteriors[indexMax]) {
					indexMax = i;
				}
			}
			
			size_t maxMaxPostPoint = maxPosteriorPoints[indexMax];
			
			#ifdef DEBUG_MAXPOINTS
				real maxMaxPosterior = maxPosteriors[indexMax];
					
				size_t maxMaxLaunchPoint = carvedLaunchPoints[indexMax];
				
				cout << "\nmaximium is the one at index " << indexMax << endl;
				cout << "maxMaxLaunchPoint = " << maxMaxLaunchPoint << endl;
				cout << "maxMaxPostPoint = " << maxMaxPostPoint << endl;
				cout << "maxMaxPosterior = " << maxMaxPosterior << endl;
			#endif
			
			// make a map of abs(# leaves diff to max post point) -> index in the maxes
			typedef std::map< int, std::vector < size_t > > IndexMap;
			typedef IndexMap::iterator IndexMapIt;
			
			IndexMap leafDiffIndexMap;
			
			for (size_t i = 0; i < numMaxPts; ++i) {
				
				int inVal = maxPosteriorPoints[i] > maxMaxPostPoint ? 
								maxPosteriorPoints[i]- maxMaxPostPoint 
								: maxMaxPostPoint - maxPosteriorPoints[i];
				
				#ifdef DEBUG_MAXPOINTS_EXTRA
					cout << "i = " << i << " and inVal = " << inVal << endl;
				#endif
						
				std::pair < IndexMapIt, bool > inPair = 
				leafDiffIndexMap.insert(
								pair< int, std::vector < size_t > >(inVal, 
									std::vector<size_t>(1,i) ));
				if (!inPair.second) inPair.first->second.push_back(i);

			}
			
				
			
			std::vector<size_t> reorderedMaxPosteriorPtIndices;
			for (IndexMapIt it = leafDiffIndexMap.begin();
				it !=  leafDiffIndexMap.end();
				++it) {
					
					#ifdef DEBUG_MAXPOINTS_EXTRA
						cout << "inserting ";
						for (size_t j = 0; j < it->second.size(); ++j) cout << "\t" << it->second[j];
						cout << " into the reordered vector " << endl;
					#endif
				
					//leafDiffIndexMap will be ordered lowest diff to highest
					reorderedMaxPosteriorPtIndices.insert(reorderedMaxPosteriorPtIndices.end(),
															it->second.begin(),
															it->second.end());
					//reorderedMaxPosteriorPtIndices.push_back(it->second);
			}
			reorderedMaxPosteriorPtIndices.resize(outOfClosest);
			
			assert(reorderedMaxPosteriorPtIndices[0] == indexMax);
			
			#ifdef DEBUG_MAXPOINTS_EXTRA
				cout << "The indices of the max points, in order of distance to the max, are: " << endl;
				for (size_t i = 0; i < reorderedMaxPosteriorPtIndices.size(); ++i) 
												cout << "\t" << reorderedMaxPosteriorPtIndices[i];
				cout << endl;
			#endif
			
			/* furthest away is at index numMaxPts-1 in reorderedMaxPosteriorPtIndices
			 * and we wanted to keep keep points, and we have 1 alread (the max)*/
			
			
			assert(keep > 1);
			
			std::vector < size_t > kept;
			
			size_t trialGap = static_cast<size_t>( std::pow(2.0, keep - 2) );
			if (trialGap <= outOfClosest - 1 ) {
				
				size_t gap = trialGap;
				#ifdef DEBUG_MAXPOINTS_EXTRA
					cout << "gap = std::pow(2.0, keep - 2) = " << gap << endl;
				
				#endif
				for (size_t i = 0; i < keep; ++i) {
					size_t keepIndex = ((outOfClosest - 1) * i)/gap; 
					assert(keepIndex < outOfClosest);
					
					#ifdef DEBUG_MAXPOINTS_EXTRA
						cout << " i= " << i 
									<< " keepIndex = ceil(((outOfClosest - 1) * i)/gap) = " 
									<< keepIndex << endl;
					
					#endif
					
					kept.push_back(keepIndex);
					
					
				}
				
			}
			else {
			
				double gap = (outOfClosest - 1.0)/(keep - 1.0);
				if (gap < 1.0) gap = 1.0;
				#ifdef DEBUG_MAXPOINTS_EXTRA
					cout << "gap = " << gap << endl;
				
				#endif
				
				for (size_t i = 0; i < keep; ++i) {
					size_t keepIndex = 0;
					if (i) {
						
						keepIndex = static_cast<size_t>(ceil(gap * i));
						
						#ifdef DEBUG_MAXPOINTS_EXTRA
							cout << "i = " << i 
								<< " and keepIndex = static_cast<size_t>(ceil(gap * i)) = " 
								<< keepIndex << endl;
						
						#endif
						if (keepIndex > outOfClosest - 1) {
							keepIndex = outOfClosest - 1;
							#ifdef DEBUG_MAXPOINTS_EXTRA
								cout << "reducing keepIndex to = " << keepIndex << endl;
							
							#endif
						}
					}
					kept.push_back(keepIndex);
					size_t j = i;
					
					while ((kept[j-1] == kept[j]) && j > 1) {
						#ifdef DEBUG_MAXPOINTS_EXTRA
							cout << "but kept[j-1] = " << kept[j-1] << " == kept[j] = " << kept[j] << endl;
						
						#endif
						assert(kept[j-1] > 0);
						kept[j-1] = kept[j-1] - 1;
						#ifdef DEBUG_MAXPOINTS_EXTRA
							cout << "reducing kept[j-1] to " << kept[j-1] << endl;
						
						#endif
						--j;
					}
					if (kept[0] == kept[1]) 
						throw std::runtime_error(
							"keepFurthestApart(...) : cannot find enough separate points");
					
				}
			}
			assert(kept.size() == keep);
			
			#ifdef DEBUG_MAXPOINTS_EXTRA
				cout << "The final locations of the indices to keep are " << endl;
				for (size_t i = 0; i < keep; ++i) cout << "\t" << kept[i];
				cout << endl << endl;
			#endif
			
			for (size_t i = 0; i < keep; ++i) {
				keepTheseOnes.push_back(reorderedMaxPosteriorPtIndices[kept[i]]);
			}
			
			/* first point pushed back will be 
			 * reorderedMaxPosteriorPtIndices[0], when i = 0, 
			 * ie the maxPost point index*/
			
			assert(keepTheseOnes.size() == keep);
			
			#ifdef DEBUG_MAXPOINTS
				cout << "The final indices to keep are " << endl;
				for (size_t i = 0; i < keep; ++i) cout << "\t" << keepTheseOnes[i];
				cout << endl << endl;
			#endif
		}
		
		return keepTheseOnes;
	}
		
	std::vector< AdaptiveHistogram* >& CarverSEB::recreateStartingPoints(
			const AdaptiveHistogram& adhBase,
			std::vector< AdaptiveHistogram* >& hists,
			const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorCarving, 
			const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			double minVol,
			const std::string& checkPostFileNameBase,	
			//int prec,
			const std::vector<size_t>& carvedLaunchPoints,
			const std::vector<size_t>& maxPosteriorPoints,
			const std::vector<real>& maxPosteriors,
			const std::vector < size_t >& keepTheseOnes,
			const gsl_rng * rgsl_base)
	{
		LOGGING_LEVEL logPQ = NOLOG; // logging level  
		#ifdef DOSAMPLES
			logPQ = LOGSAMPLES;
		#endif
		
		for (size_t i = 0; i < carvedLaunchPoints.size() ; ++i) {
			
			// is this index in our ones to keep?
			vector<size_t>::const_iterator fit 
				= find (keepTheseOnes.begin(), keepTheseOnes.end(), i);
			if (fit != keepTheseOnes.end()) {
				
				size_t carvedPt = carvedLaunchPoints[i];
				size_t maxPostPoint = maxPosteriorPoints[i];
				real maxPost = maxPosteriors[i];
				
				std::vector<real> thisPosteriorVec; 
				std::vector<real> thisLoglikVec; 
				
				bool successPQ = false;
				AdaptiveHistogram* adhPtr;
				
				try {
					
					adhPtr = new AdaptiveHistogram(adhBase);

					
					successPQ = recreateStartingPoint(
							*adhPtr,
							evaluatorCarving, 
							evaluatorSEB, 
							logPrior, 
							minPoints,
							minVol,
							carvedPt,
							maxPostPoint,
							thisPosteriorVec, 
							thisLoglikVec, 
							rgsl_base,
							logPQ);
				}
				catch (...) {
					if (adhPtr != NULL) delete adhPtr;
					adhPtr = NULL;
					
					throw;
				}	
				if (successPQ) {
					
					#ifdef DEBUG_MAXPOINTS
					
						cout << "\nAt recreate point with leaves " << maxPostPoint
						<< "\nfrom carving point " << carvedPt << endl;
					
						size_t leaves = adhPtr->getRootLeaves();
						
						if (leaves != maxPostPoint) 
							cout << "\nLeaves in recreated hist is " << leaves << endl;
					
						assert(leaves == maxPostPoint);
						
						real maxLik= adhPtr->getLogLikelihood();
						real lnPrior = logPrior(leaves-1);
						
						real thisPost = maxLik + lnPrior;
						
						cout << "maxPosterior from earlier pq was " << maxPost
								<< " and posterior calculated here is " << thisPost << endl;
						cout << "(logLik here is " << maxLik << " and logPrior is " 
													<< lnPrior << ")" << endl;
					#endif
					
					hists.push_back(adhPtr);
					
					// output log-posterior
					if (!checkPostFileNameBase.empty()) {
						ostringstream oss;
						oss << checkPostFileNameBase << "_" << carvedPt 
									<< "_" << maxPostPoint << ".txt";
						string checkPostFileName = oss.str();
						
						std::vector < const subpavings::RealVec* > dataPtrs;
						dataPtrs.push_back(&thisPosteriorVec);
						dataPtrs.push_back(&thisLoglikVec);
						std::vector < std::string > colNames;
						colNames.push_back("Post");
						colNames.push_back("lnLik");
						
						outputToFileVertical(dataPtrs, 
										colNames,
										checkPostFileName);
					}	
					
				}
				else {
					cerr << "Could not recreate state for max point with leaves " 
					<< maxPostPoint << endl;
					if (adhPtr != NULL) delete adhPtr;
					adhPtr = NULL;
				}
			}
			
		}
		return hists;
	}

		/* I tried seeing if I could do recreation more efficiently - idea was:
	 * Can reuse if old hist has same carving point and max posterior larger
	 * But first, we don't necessarily have points to recreate with same
	 * carving point, and secondly hist we find to base next one on will
	 * also have the SEB splits and so more leaves than the carving stop
	 * point, so the adh priority split method will throw an exception, and
	 * it is not worth trying to keep a carved 'base' because there
	 * won't be enough using this same base to make it worth while I think*/


	bool CarverSEB::recreateStartingPoint(
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
			LOGGING_LEVEL logPQ)
	{
		
		AdaptiveHistogram::PrioritySplitQueueEvaluator
					thisEvaluatorCarving(evaluatorCarving.getMeasurer(), carvedPt);
		thisEvaluatorCarving.setUsingCritStop(false); // should be default anyway
			
		AdaptiveHistogram::PrioritySplitQueueEvaluator
					thisEvaluatorSEB(evaluatorSEB.getMeasurer(), maxPostPoint);
		thisEvaluatorSEB.setUsingCritStop(false); // should be default anyway
			
		return recreateStartingPoint(
			adh, // changed
			thisEvaluatorCarving, 
			thisEvaluatorSEB, 
			logPrior, 
			minPoints,
			minVol,
			posteriorVec, 
			loglikVec, 
			rgsl_base,
			logPQ);

	}


	bool CarverSEB::recreateStartingPoint(
			AdaptiveHistogram& adh, // changed
			const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorCarving, 
			const AdaptiveHistogram::PrioritySplitQueueEvaluator& evaluatorSEB, 
			LogMCMCPrior& logPrior, 
			size_t minPoints,
			double minVol,
			std::vector<cxsc::real>& posteriorVec, 
			std::vector<cxsc::real>& loglikVec, 
			const gsl_rng * rgsl_base,
			LOGGING_LEVEL logPQ)
	{
				
		gsl_rng * rgsl_max = NULL;
		
		
		try {
			
			// prng in pristine state
			rgsl_max = gsl_rng_clone(rgsl_base);
			
			#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
				cout << "In recreate hists,\n carved max leaves = " 
						<< evaluatorCarving.getMaxLeaves() << endl;
				cout << "In recreate hists,\n seb max leaves = " 
						<< evaluatorSEB.getMaxLeaves() << endl;
				cout << "current hist has leaves = " << adh.getRootLeaves() << endl;
			#endif
			
			bool shiftCatalan = false;
							
			bool successPQ = adh.prioritySplitMaxLik(
						evaluatorCarving, 
						evaluatorSEB, 
						logPQ, 
						minPoints,
						minVol,
						logPrior,
						posteriorVec,
						loglikVec,
						rgsl_max,
						shiftCatalan);
			
			#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
				cout << "\nhist now has leaves = " << adh.getRootLeaves() << endl;
			#endif
			
			gsl_rng_free(rgsl_max);
			rgsl_max = NULL;
			
			return successPQ;
		}
		catch (...) {
			if (rgsl_max != NULL) gsl_rng_free(rgsl_max);
			rgsl_max = NULL;
			
			throw;
		}	

	}


	 //add to some overdispersed starting points
	std::vector< AdaptiveHistogram* >& CarverSEB::addToHists(
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
			const gsl_rng * rgsl_base)
	{
		#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
            size_t numPosts = posteriorVec.size();
			cout << "\n In addToHists: toFind = " << toFind << endl;
			cout << " and baseIndex = " << baseIndex
			<< " and numPosts = " << numPosts << endl;
				
		#endif
			
		std::vector<size_t> keepTheseOnes(toFind);
		std::vector<size_t> finalCarvedLaunchPoints(toFind);
		std::vector<size_t> finalMaxPosteriorPoints(toFind);
		std::vector<real> finalMaxPosteriors(toFind);
		
		/* find the points to recreate and hold them*/
		for (size_t i = 0; i < toFind; ++i) {
			size_t thisIndex = baseIndex + static_cast<size_t>(i*dblGap);
			#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
				cout << "i = " << i << endl;
				cout << "thisIndex = " << thisIndex << endl;
				
			#endif	
			assert(thisIndex < posteriorVec.size());
			size_t thisPostPt = startLeaves+thisIndex;
			finalCarvedLaunchPoints[i] 
					= (carvedPt < thisPostPt ? carvedPt : thisPostPt);
			
			finalMaxPosteriorPoints[i] = (thisPostPt);
			#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
				cout << "finalMaxPosteriorPoints[" << i 
					<< "] = " << finalMaxPosteriorPoints[i] << endl;
				
			#endif		
			finalMaxPosteriors[i] = posteriorVec[thisIndex];
			keepTheseOnes[i] = i; 
			
			#ifdef DEBUG_MAXPOINTS_OVERDISPERSED
				cout << "finalMaxPosteriors[" << (i) 
					<< "] = " << finalMaxPosteriors[i] << endl; // but it will not always
				
			#endif				
		}
		recreateStartingPoints(
				adh,
				histsToAddTo, // filled in process
				evaluatorCarving, 
				evaluatorSEB, 
				logPrior, 
				minPoints,
				minVol,
				checkPostFileNameBase,	
				finalCarvedLaunchPoints,
				finalMaxPosteriorPoints,
				finalMaxPosteriors,
				keepTheseOnes,
				rgsl_base);
		
		assert(histsToAddTo.size() >= toFind);

		return histsToAddTo;
	}
}

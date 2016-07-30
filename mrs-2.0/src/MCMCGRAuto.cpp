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
\brief Definitions for new version of a class to do MCMC with Gelman-Rubin 
* heuristic automatic sampling rule.
 */

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCGRAuto.hpp"

#include "realmappedspnode.hpp"
#include "sptools.hpp"
#include "sptypes.hpp"

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions

#include <gsl/gsl_randist.h>




using namespace cxsc;
using namespace subpavings;
using namespace std;

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	 
							int runID, 
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType) const 

{
	const std::string& scalarsFileName("scalars.txt"); 
	
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	return doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName,
							logPrior, //gat41 
							samplingType); 	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	 
							int runID, 
							LogMCMCPrior& logPrior,  //gat41
							SAMPLING_TYPE samplingType) const

{
	const std::string& scalarsFileName("scalars.txt"); 
	
	return doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							logPrior, //gat41
							samplingType); //gat41	
}
//-----------------------------
		
PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
				std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
				int maxLoops, 
				int samplesNeeded,
				int thinout, 
				size_t minPoints,
				int rhatFlagCounterThreshold, 	 
				int runID, 
				SAMPLING_TYPE samplingType,
				LOGGING_LEVEL logging,
				LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	real minVol(0.0);
	
	const std::string& scalarsFileName("scalars.txt"); 
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
				std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
				int maxLoops, 
				int samplesNeeded,
				int thinout, 
				int samplesMonitored,
				size_t minPoints,
				int rhatFlagCounterThreshold, 	 
				int runID, 
				SAMPLING_TYPE samplingType,
				LOGGING_LEVEL logging,
				LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	const std::string& scalarsFileName("scalars.txt"); 
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

//----------------------------
			

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	 
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	SAMPLING_TYPE samplingType = NORMAL;
	LOGGING_LEVEL logging = NOLOG; 
	LOGGING_LEVEL loggingInChangeStates = NOLOG;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							real minVol,
							int rhatFlagCounterThreshold, 	 
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	SAMPLING_TYPE samplingType = NORMAL;
	LOGGING_LEVEL logging = NOLOG; 
	LOGGING_LEVEL loggingInChangeStates = NOLOG;
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior);	//gat41
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout,
							int samplesMonitored, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	 
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior,   //gat41
							SAMPLING_TYPE samplingType) const 
{
	LOGGING_LEVEL logging = NOLOG; 
	LOGGING_LEVEL loggingInChangeStates = NOLOG;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored, 
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	 
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

//------------------------------


PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(
							histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	  
							runID,
							scalarsFileName,  
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,
							size_t minPoints,
							int rhatFlagCounterThreshold, 	  
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	real minVol(0.0);
	
	return _doMCMCGRAuto(
							histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							rhatFlagCounterThreshold, 	  
							runID,
							scalarsFileName,  
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41
}

//-------------------------------

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID, 
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType) const 

{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	const std::string& scalarsFileName("scalars.txt"); 
	return doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							runID,
							scalarsFileName,
							logPrior, //gat41 
							samplingType); 	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,
							size_t minPoints,
							int runID, 
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType) const 

{
	const std::string& scalarsFileName("scalars.txt"); 
	
	return doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout,
							samplesMonitored, 
							minPoints,
							runID,
							scalarsFileName,
							logPrior, 
							samplingType); //gat41	
}

//------------------------------
		
PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
				std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
				int maxLoops, 
				int samplesNeeded,
				int thinout, 
				size_t minPoints,
				int runID, 
				SAMPLING_TYPE samplingType,
				LOGGING_LEVEL logging,
				LOGGING_LEVEL loggingInChangeStates,
				LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	real minVol(0.0);
	
	const std::string& scalarsFileName("scalars.txt"); 
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
				std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
				int maxLoops, 
				int samplesNeeded,
				int thinout, 
				int samplesMonitored,
				size_t minPoints,
				int runID, 
				SAMPLING_TYPE samplingType,
				LOGGING_LEVEL logging,
				LOGGING_LEVEL loggingInChangeStates,
				LogMCMCPrior& logPrior) const //gat41
{
	const std::string& scalarsFileName("scalars.txt"); 
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

//----------------------------

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	SAMPLING_TYPE samplingType = NORMAL;
	LOGGING_LEVEL logging = NOLOG; 
	LOGGING_LEVEL loggingInChangeStates = NOLOG;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,   
							LogMCMCPrior& logPrior, //gat41
							SAMPLING_TYPE samplingType) const 
{
	LOGGING_LEVEL logging = NOLOG; 
	LOGGING_LEVEL loggingInChangeStates = NOLOG;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName, 
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41	
}

//---------------------

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	int samplesMonitored = getDefaultSamplesMonitored();
	if (thinout > 1) samplesMonitored = (histPtrsVec.size()*samplesMonitored/thinout);
	if (!samplesMonitored) samplesMonitored = histPtrsVec.size()*2;
	
	real minVol(0.0);
	
	return _doMCMCGRAuto(
							histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName,  
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41
}

PiecewiseConstantFunction MCMCGRAuto::doMCMCGRAuto(
							std::vector< subpavings::AdaptiveHistogram* >& histPtrsVec,
							int maxLoops, 
							int samplesNeeded,
							int thinout, 
							int samplesMonitored,
							size_t minPoints,
							int runID,
							const std::string& scalarsFileName,  
							SAMPLING_TYPE samplingType,
							LOGGING_LEVEL logging,
							LOGGING_LEVEL loggingInChangeStates,
							LogMCMCPrior& logPrior) const //gat41
{
	real minVol(0.0);
	
	return _doMCMCGRAuto(
							histPtrsVec,
							maxLoops, 
							samplesNeeded,
							thinout, 
							samplesMonitored,
							minPoints,
							minVol,
							runID,
							scalarsFileName,  
							samplingType,
							logging,
							loggingInChangeStates,
							logPrior); //gat41
}

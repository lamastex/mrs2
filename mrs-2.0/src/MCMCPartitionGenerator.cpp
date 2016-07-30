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
\brief MCMCPartitionGenerator definitions.
 */
/*header for controlling debugging*/
#include "debug_control.hpp"

#include "MCMCPartitionGenerator.hpp"

#include "sptools.hpp" 

#include <iostream>  // input and output streams
#include <iterator>
#include <sstream>
#include <ctime>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <climits> // for ULONG_MAX

#include <gsl/gsl_randist.h>

#include "dot.hpp"

#include "rmath.hpp"

//#define CHECK_MIN_PROB // checking min probabilities - does not really matter
//#define CHECK_MIN // checking mins in working out probabilities
//#define MYDEBUG
//#define MYDEBUG_SPECIAL
//#define MYDEBUG_MISSED_PROBS

//#define PROBS

#ifdef NDEBUG
	#undef CHECK_MIN_PROB
	#undef CHECK_MIN
	#undef MYDEBUG
	#undef MYDEBUG_SPECIAL
	#undef MYDEBUG_MISSED_PROBS

	
#endif


#define MAX_TRIES 3 // max tries to find correct probability

//using namespace cxsc;
using namespace std;

namespace subpavings {
	
	MCMCPartitionGenerator::MCMCPartitionGenerator(unsigned long int seed)
	 : rgsl(NULL), instructionsPtr(NULL), workspacePtr(NULL)
	{

		try {
			
			makeCatalans();
			rgsl = gsl_rng_alloc(gsl_rng_mt19937);
			
			gsl_rng_set(rgsl, seed);
			
		}
		catch (...) {
			try {
				if (rgsl != NULL) gsl_rng_free(rgsl);
			}
			catch (...) {}
			clearInstructions();
			throw;
		}
	}
			
	MCMCPartitionGenerator::~MCMCPartitionGenerator()
	{
		try {
			gsl_rng_free(rgsl);
		}
		catch (...) {} // catch and swallow
		clearInstructions();
	}
	
	unsigned long int MCMCPartitionGenerator::generateStatePartition(
			unsigned long int numLeaves) const
	{
	
		#ifdef MYDEBUG
			cout << "\nIn generateState, numLeaves = " << numLeaves << endl;
			cout << "\n" << endl;
		#endif
		
		if (!numLeaves) throw std::logic_error(
			"MCMCPartitionGenerator::generateStatePartition(...) : numLeaves = 0");
		
		unsigned long int left = 0;
		
		if (numLeaves <= 2) left = numLeaves - 1; // 1 if numLeaves is 2
		
		else { // numLeaves > 2 
		
			unsigned long int k = numLeaves-1; // k > 1
		
			/* pick some partition of numLeaves at random
			 * anywhere between 1 and numLeaves-1 */
			
			// rand in (0,1]
			double rand = 1.0-gsl_rng_uniform(rgsl);
			//double rand = 0.500005; fix missed probability for 9 leaves '(also 0.499995)
			
			/* If we can use the precomputed Catalan numbers it is much quicker */
			if (numLeaves <= catalans.size()) { // catalans included C0
			
				//pick any of the numbers in 0..Ck-1
				unsigned long int catK = catalans[k];
				#ifdef MYDEBUG
					cout << "using catalans table, catK = " << catK << endl;
					cout << "1/catK = " << (cxsc::real(1.0)/(1.0*catK)) << endl;
					cout << "1/catK = " << (cxsc::interval(1.0)/(1.0*catK)) << endl;
					cout << "gsl_rng_max = " << gsl_rng_max(rgsl) << endl;
				#endif
				
				//unsigned long int choose = gsl_rng_uniform_int(rgsl, catK);
				unsigned long int sum = 0;
				
				#ifdef MYDEBUG
					cout << "rand = " << rand << endl;
				#endif
			
				
				for (unsigned long int i = 1; i <= numLeaves/2; ++i) {
					
					long unsigned catComp = catalans[i-1]*catalans[k-i];
					sum += catComp;
					
					#ifdef MYDEBUG
						cout << "adding " << catComp 
								<< " to sum, sum = " << sum << endl;
					#endif
					if (rand <= (sum*1.0)/catK) {
						left = i;
						#ifdef MYDEBUG
							cout << "picked " << left << endl;
						#endif
					
						break;
					}
					if (rand > (catK - (sum*1.0))/catK) {
						left = k - i + 1;
						#ifdef MYDEBUG
							cout << "picked " << (k - i + 1) << endl;
						#endif
						break;
					}
				}
			
			}
			
			else { // can't use catalans table
				
						
				/* prob of split 1|numLeaves-2 is Ck-1/Ck where k = numLeaves-1
				Ck-1/Ck = 0.5*k+1/(2k-1) (which = 1 if k = 1 ie numLeaves = 1*/
				double firstEdgeBand = 0.5*(k+1.0)/(2*k-1.0);
				
				// where we choose to go
				unsigned long int chosenS = 1;
				
				#ifdef MYDEBUG
					cout << "k = " << k << endl;
					cout << "rand = " << rand << endl;
					cout << "choose " << chosenS << " if rand <= " << firstEdgeBand << endl;
					cout << "choose k if rand > " << (1.0-firstEdgeBand) << endl;
				#endif
				
				if ( rand <= firstEdgeBand ) {
					left = chosenS;
					#ifdef MYDEBUG
						cout << "choose " << chosenS << endl;
					#endif
				}
				/* if k == 1, we should have gone to 1 by now */
				else if ( rand > 1.0 - firstEdgeBand ) {
					left = k - chosenS + 1;  // = k = numLeaves - 1
					#ifdef MYDEBUG
						cout << "choose k" << endl;
					#endif
				}
				/* if k == 2, we should have gone to either 1 or k = 2 by now */
				else { // not going down the edge
				
					assert( k > 2);
				
					++chosenS; // = 2;
				
					/* the next probability is firstEdgeBand*1*0.5*(k/(2k-3)), so 
					 * adding that to firstEdgeBand we get 
					 * cut-off firstEdgeBand*(1+0.5*(k/(2k-3))) */
					double secondEdgeBand = firstEdgeBand * (1+ 0.5*k/(2*k-3.0));
					
					#ifdef MYDEBUG
						cout << "not on the edge, try 1 in from edge " << k << endl;
						cout << "choose " << chosenS << " if rand <= " << secondEdgeBand << endl;
						cout << "choose k-1 = " << (k - chosenS + 1) << " if rand > " << (1.0-secondEdgeBand) << endl;
					#endif
					
					
					if (rand <= secondEdgeBand) {
						left = chosenS;
						#ifdef MYDEBUG
							cout << "choose " << chosenS << endl;
						#endif
					}
					/* if k == 3, we should have gone to either 1 or 2 or 3 by now */
					else if ( rand > 1.0 - secondEdgeBand) {
						left = k - chosenS + 1; // k-1
						#ifdef MYDEBUG
							cout << "choose k-1 = " << (k - chosenS + 1) << endl;
						#endif
					}
					/* if k == 4, we should have gone to either 1 or 2 or 3 or 4 by now */
					else {
						
						assert( k > 4);
						
						++ chosenS; //= 3;
						
						/* the next probability is 
						 * 	2 * (secondEdgeBand - firstEdgeBand) * 0.5 * (k-1.0)/(2*k-5.0)), 
						 * which we add to secondEdgeBand to get the cut-off point */
						double thirdEdgeBand = secondEdgeBand + 
							(secondEdgeBand - firstEdgeBand) * (k-1.0)/(2*k-5.0);
					
						#ifdef MYDEBUG
							cout << "not on the edge, try 2 in from edge " << k << endl;
							cout << "choose " << chosenS << " if rand <= " << thirdEdgeBand << endl;
							cout << "choose k-2 = " << (k - chosenS + 1) << " if rand > " << (1.0-thirdEdgeBand) << endl;
						#endif
						
						
						if (rand <= thirdEdgeBand) {
							left = chosenS;
							#ifdef MYDEBUG
								cout << "choose " << chosenS << endl;
							#endif
						}
						/* if k == 5, we should have gone to either 1 or 2 or 3 or 4 or 5 by now */
						else if ( rand > 1.0 - thirdEdgeBand) {
							left = k - chosenS + 1; // k-2
							#ifdef MYDEBUG
								cout << "choose k-2 = " << (k - chosenS + 1) << endl;
							#endif
						}
						/* if k == 6, we should have gone to either 1 or 2 or 3 or 4 or 5 or 6 by now */
						
						else {
							// use our approximated probabilities
								
							assert(k > 6);
							
							++chosenS; // = 4
							unsigned long int startChosenS = chosenS;
								
							cxsc::real realK(1.0*k);
							
							/* rescale rand to go from 0.0 up and in units without the big divisor */
							cxsc::real rescaledRand( (rand-thirdEdgeBand)
										* 4.0*cxsc::SqrtPi_real/(cxsc::sqrt(realK)*(realK+1)) );
							/* top of range similarly */
							cxsc::real upperBound( (1.0 - 2*thirdEdgeBand) * 
									4.0*cxsc::SqrtPi_real/(cxsc::sqrt(realK)*(realK+1)) );
							
							cxsc::real fudge(1.0);
							
							#ifdef MYDEBUG
								cout << "not on the first or second or third edge" << endl;
								cout << "rescaledRand = " << _double(rescaledRand) << endl;
								cout << "upperBound = " << _double(upperBound) << endl;
							#endif
							
							assert (rescaledRand <= upperBound);
							
							cxsc::dotprecision finalAccProbsRescaled(0.0);
							
							int tries = 0;
							while (!left && (tries < MAX_TRIES)) {
								
								cxsc::dotprecision accProbsRescaled(0.0);
								
								#ifdef MYDEBUG
									cout << "try number " << (tries + 1) << endl;
									cout << "fudge = " << fudge << endl;
								#endif
								
								#ifdef MYDEBUG_SPECIAL
									if (tries) {
										cout << "try number " << (tries + 1) << endl;
										cout << "fudge = " << fudge << endl;
									}
								#endif
							
								// step up from here, looking to send to left = S or left = k - S + 1
								for (chosenS = startChosenS ; chosenS <= numLeaves/2; ++chosenS) {
									
									cxsc::real realS(1.0*chosenS); 
									
									// what's the new top of the lower band, ie to send to left = S
									cxsc::accumulate( accProbsRescaled, 
												cxsc::exp(-1.0/8.0*(realS/(realK*(realK - realS))+1/(realS-1)))
											/cxsc::sqrt((realS - 1.0)*(realK - realS)),
												fudge/(realS*(realK - realS + 1)) );
									
									
									#ifdef CHECK_MIN
										if (cxsc::rnd(accProbsRescaled) < cxsc::MinReal) {
											cout << "\nAlert: cxsc::rnd(accProbsRescaled) = " << cxsc::rnd(accProbsRescaled)
												<< " < cxsc::MinReal = " << cxsc::MinReal << endl;
											throw std::runtime_error("MinReal exceeded: accProbsRescaled");
										}
										if (cxsc::rnd(upperBound - accProbsRescaled) < cxsc::MinReal) {
											cout << "\nAlert: cxsc::rnd(upperBound - accProbsRescaled) = " << cxsc::rnd(upperBound - accProbsRescaled)
												<< " < cxsc::MinReal = " << cxsc::MinReal << endl;
											throw std::runtime_error("MinReal exceeded: upperBound - accProbsRescaled");
										}
									#endif
									#ifdef MYDEBUG
									{
										cout << "S = " << chosenS << endl;
										cout << "fudge = " << fudge << endl;
										cout << "accProbsRescaled = " << _double(cxsc::rnd(accProbsRescaled)) << endl;
										cout << "choose " << chosenS << " if rescaledRand <= accProbsRescaled" << endl;
										cout << "choose " << (k - chosenS + 1) 
												<< " if rescaledRand > (upperBound - accProbsRescaled = " 
												<< _double(cxsc::rnd(upperBound - accProbsRescaled))  << endl;
									}
									#endif
									
									
									if (rescaledRand <= accProbsRescaled) {
										left = chosenS;
										#ifdef MYDEBUG
											cout << "choose left = " << left << endl;
										#endif
										
										break;
									}
									else if (rescaledRand > (upperBound - accProbsRescaled)) {
										left = k - chosenS + 1;
										#ifdef MYDEBUG
											cout << "choose left = " << left << endl;
										#endif
										
										break;
									}
									#ifdef MYDEBUG_SPECIAL
									{
										if (tries && left) {
											cout << "S = " << chosenS << endl;
											cout << "fudge = " << fudge << endl;
											cout << "accProbsRescaled = " << _double(cxsc::rnd(accProbsRescaled)) << endl;
											cout << "choose " << chosenS << " if rescaledRand <= accProbsRescaled" << endl;
											cout << "choose " << (k - chosenS + 1) 
													<< " if rescaledRand > (upperBound - accProbsRescaled = " 
													<< _double(cxsc::rnd(upperBound - accProbsRescaled))  << endl;
											
											cout << "\nchoose left = " << left << endl;
										}
									}
									#endif
								} // end of for loop
								
								finalAccProbsRescaled = accProbsRescaled;
								
								if (!left) {
									/* possibility we could get to the end and still have not found left,
									 * because our probabilities in this for loop are approximated
									 * and are slightly too low ... so try again with a fudge factor
									 * and work out from the middle - state we want will be towards the middle*/ 
																	
									/* gap = upperBound - 2*accProbsRescaled
									   want fudge = upperBound/(upperBound-gap) 
													= upperBound/(2*accProbsRescaled) */
									fudge = upperBound/(2*cxsc::rnd(accProbsRescaled));
									#ifdef MYDEBUG_SPECIAL
										cout << "first not found left" << endl;
										cout << "rand = " << rand << endl;
										cout << "rescaledRand = " << rescaledRand << endl;
										cout << "upperBound = " << upperBound << endl;
										cout << "2*cxsc::rnd(accProbsRescaled) = " << (2*cxsc::rnd(accProbsRescaled)) << endl;
									#endif
									
									accProbsRescaled = 0.0;
									accumulate(accProbsRescaled, 0.5, upperBound);
									
									#ifdef MYDEBUG_SPECIAL
										cout << "adjusted accProbsRescaled) = " << (cxsc::rnd(accProbsRescaled)) << endl;
									#endif
									
									// step down from here, looking to send to left = S or left = k - S + 1
									for (chosenS = numLeaves/2; chosenS >= startChosenS; --chosenS) {
										
										cxsc::real realS(1.0*chosenS); 
										
										// what's the new top of the lower band, ie to send to left = S
										cxsc::accumulate( accProbsRescaled, 
													-cxsc::exp(-1.0/8.0*(realS/(realK*(realK - realS))+1/(realS-1)))
												/cxsc::sqrt((realS - 1.0)*(realK - realS)),
													fudge/(realS*(realK - realS + 1)) );
										
										#ifdef MYDEBUG_SPECIAL
										{
											cout << "S = " << chosenS << endl;
											cout << "fudge = " << fudge << endl;
											cout << "accProbsRescaled = " << _double(cxsc::rnd(accProbsRescaled)) << endl;
											cout << "choose " << chosenS 
													<< " if rescaledRand <= 0.5*upperBound && rescaledRand > accProbsRescaled" << endl;
											cout << "else choose " << (k - chosenS + 1) 
													<< " if rescaledRand > 0.5*upperBound && rescaledRand <= (upperBound - accProbsRescaled = " 
													<< _double(cxsc::rnd(upperBound - accProbsRescaled))  << endl;
											
										}
										#endif
										
																			
										if ((rescaledRand <= 0.5*upperBound) && (rescaledRand > accProbsRescaled)) {
											left = chosenS;
											#ifdef MYDEBUG
												cout << "choose left = " << left << endl;
											#endif
											
											break;
										}
										else if ((rescaledRand > 0.5*upperBound) && 
													(rescaledRand <= (upperBound - accProbsRescaled)) ){
											left = k - chosenS + 1;
											#ifdef MYDEBUG
												cout << "choose left = " << left << endl;
											#endif
											
											break;
										}
										#ifdef MYDEBUG_SPECIAL
										{
											if (left) {
												cout << "S = " << chosenS << endl;
												cout << "fudge = " << fudge << endl;
												cout << "accProbsRescaled = " << _double(cxsc::rnd(accProbsRescaled)) << endl;
												cout << "choose " << chosenS << " if rescaledRand > accProbsRescaled" << endl;
												cout << "choose " << (k - chosenS + 1) 
														<< " if rescaledRand <= (upperBound - accProbsRescaled = " 
														<< _double(cxsc::rnd(upperBound - accProbsRescaled))  << endl;
												cout << "\nchoose left = " << left << endl;
											}
										}
										#endif
									} // end of for second loop
									
									
									
								} // end of if !left
								
								/* still possibility we could get to the end and still have not found left?
								 * because our probabilities in this for loop are approximated
								 * and are slightly too low ... so try again with a fudge factor*/ 
								if (!left) {
									
									/* gap = upperBound - 2*accProbsRescaled
									   want fudge = upperBound/(upperBound-gap) 
													= upperBound/(2*accProbsRescaled) */
									fudge = upperBound/(2*cxsc::rnd(finalAccProbsRescaled));
									#ifdef MYDEBUG_SPECIAL
										cout << "** STILL not found left, upperBound = " << upperBound << endl;
										cout << "rand = " << rand << endl;
										cout << "2*cxsc::rnd(finalAccProbsRescaled) = " << (2*cxsc::rnd(finalAccProbsRescaled)) << endl;
									#endif
								}
								
								++tries;
							} // end of tries while loop
							
							/* remote possibility we could get to the end and still have not found left?
							 * I think that if that's the case we should pick the middle*/
							if (!left) {
								chosenS = numLeaves/2; // integer division
								left = chosenS;
								if (rand > 0.5) left = k - chosenS + 1; // same if numLeaves is even
								#if defined (MYDEBUG_SPECIAL) || defined (MYDEBUG_MISSED_PROBS)
								{
									cout << "\n********* missed probability **************\n" << endl;
								
									cout << "S = " << chosenS << endl;
									cout << "rescaledRand = " << _double(rescaledRand) << endl;
									cout << "accProbsRescaled = " << _double(cxsc::rnd(finalAccProbsRescaled)) << endl;
									cout << "upperBound - accProbsRescaled = " 
											<< _double(cxsc::rnd(upperBound - finalAccProbsRescaled))  << endl;
																	
									cout << "\n*******************************************\n" << endl;
								}
								#endif
								
							}  
						} // end the else to deal with cases where we are not in first or second or third 'edge'		 
					}// end the else to deal with cases where we are not in first or second
				}// end the else to deal with cases where we are not in first 'edge'
			
			} // end else we cannot use catalans table 	
			assert (left != 0);
		} // end if numLeaves > 2
			
		#ifdef MYDEBUG
			cout << "left = " << left << endl;
		#endif
		
		if (NULL != instructionsPtr) instructionsPtr->push_back(left);
		
		return left;
		
	}
	
	//calc sum(ln(1/prob))
	unsigned long int MCMCPartitionGenerator::generateNaturalStatePartition(
		unsigned long int numLeaves) const
	{
		if (!numLeaves) throw std::logic_error(
			"MCMCPartitionGenerator::generateStatePartition(...) : numLeaves = 0");
		
		
		cxsc::dotprecision lnProb(0.0); // ln(1.0)
				
		unsigned long int left = 0;
		
		if (numLeaves <= 2) left = numLeaves - 1; // 1 if numLeaves is 2
		
		else { // numLeaves > 2 
		
			cxsc::accumulate(lnProb, 1.0, cxsc::ln(1.0*(numLeaves-1)));
			// prob= 1.0/(numLeaves-1);
			
			/* pick some partition of numLeaves at random
			 * anywhere between 1 and numLeaves-1 */
			left = gsl_rng_uniform_int(rgsl, numLeaves-1) + 1;
			// will be 1 if numLeaves = 2;
			
			//lnProb += generateStatePartitioned(state, startIndex, currentLevel+1, left, rgsl);
			
			
			//lnProb += generateStatePartitioned(state, startIndex+left, currentLevel+1, numLeaves-left, rgsl);
			
			
		}
		
		return left;
		
	}
	
	/* dec = 1 is equivalent to '(', 0 to ')'
	 * */
	int MCMCPartitionGenerator::generateKnuthDecision(
				unsigned long int p,
				unsigned long int q,
				bool across) const
	{
		int dec = 0;
		
		/* if p == 0 we know we'll generate 0 */
		if (p > 0) {
			
			dec = 1;
		
			assert (!(q < p)); // q >= p
				
			/* generate random integer in range [0 , (q+p)(q-p+1)) 
			 * gsl_rng_uniform_int(r, n) generates 0 <= u < n 
			 * (ie u in [0, n-1]) from r*/ 
			
            /* check lnMax = std::log(q+p) + std::log(q - p + 1) */
			assert ( !( (std::log(q+p) + std::log(q - p + 1)) 
                                    > std::log(gsl_rng_max(rgsl))) );
			unsigned long int x = gsl_rng_uniform_int(rgsl, (q+p)*(q-p+1) );
			if ( x < (q+1)*(q-p) ) {
				dec = 0;
			}
		}
		/* if p = 0, dec = 0 */
		
		#if(0)
			if (NULL != instructionsPtr) instructionsPtr->push_back(dec);
		#else
		if (NULL != workspacePtr && NULL != instructionsPtr) {
			size_t n = instructionsPtr->size();
			assert ( workspacePtr->size() == n );
			/* if 1, add workspace to instructions,
			 * then add a new 1 slot to workspace and instructions*/
			if (dec) {
				
				std::transform ( instructionsPtr->begin(), instructionsPtr->end(),
                             workspacePtr->begin(), instructionsPtr->begin(),
                             std::plus<unsigned long int>() );
				
				workspacePtr->push_back(1);
				instructionsPtr->push_back(1);
			}
			/* else (0), close the last open slot on the workspace*/
			else {
				for (vector<unsigned long int>::reverse_iterator 
					rit = workspacePtr->rbegin() ;
					rit < workspacePtr->rend();
					++rit ) {
				
					if (*rit) {
						*rit = 0;
						break;
					}
				}
				/* closing a cherry*/
				if (across) {
					workspacePtr->resize(n-1);
					instructionsPtr->resize(n-1);
				}
			}
		}
		#endif
		return dec;
	}
	
	cxsc::real MCMCPartitionGenerator::getLnCatalanRatio(
			unsigned long int k, unsigned long int kPrime) const
	{
		size_t catalanMaxK = catalans.size()-1; // catalans included C0
	
		if (k <= catalanMaxK && kPrime <= catalanMaxK) { 
				return (cxsc::ln(1.0*catalans[k]) 
					- cxsc::ln(1.0*catalans[kPrime]));
		}
		else if (k <= catalanMaxK && kPrime > catalanMaxK) { 
						
				return ( cxsc::ln(1.0*catalans[k]) 
					- (2.0*kPrime*cxsc::ln(2.0) 
						- 0.5*cxsc::ln(cxsc::Pi_real*(1.0*kPrime)) 
						- cxsc::ln(kPrime+1.0)- 1.0/(8*kPrime)) );
		}
		else if (k > catalanMaxK && kPrime <= catalanMaxK) { 
			
				return ( (2.0*k*cxsc::ln(2.0) 
					- 0.5*cxsc::ln(cxsc::Pi_real*(1.0*k)) 
					- cxsc::ln(k+1.0)- 1.0/(8*k)) 
						
						- cxsc::ln(1.0*catalans[kPrime]) );
		}
		else  { 
			
			/* don't do (k - kPrime because result should be an unsigned long int
			* and this will go horribly wrong if kPrime > k */
			return ( 2.0*k*cxsc::ln(2.0) - 2.0*kPrime*cxsc::ln(2.0) 
				- 0.5*(cxsc::ln(1.0*k) - cxsc::ln(1.0*kPrime)) 
				- (cxsc::ln(k + 1.0) - cxsc::ln(kPrime + 1.0))
				- 1.0/8*(1.0/k - 1.0/kPrime));
		}
	}
	
	void MCMCPartitionGenerator::initialiseInstructions(
				size_t numLeaves ) const
	{
		clearInstructions();
		instructionsPtr = new std:: vector < unsigned long int >();
		instructionsPtr->reserve(numLeaves);
		workspacePtr = new std:: vector < unsigned long int >();
	}

	std:: vector < unsigned long int > MCMCPartitionGenerator::getInstructions() const
	{
		return *instructionsPtr;
	}
	
	void MCMCPartitionGenerator::clearInstructions() const
	{
		try {
			if (instructionsPtr != NULL) {
				delete instructionsPtr;
				instructionsPtr = NULL;
			}
			if (workspacePtr != NULL) {
				delete workspacePtr;
				workspacePtr = NULL;
			}
		}
		catch (...) {}
	}
	
	void MCMCPartitionGenerator::makeCatalans()
	{
		catalans.push_back(1);
		catalans.push_back(1);
		
		size_t i = 2;
		while (catalans.back() < ULONG_MAX/4) {
			unsigned long int cat = 0;
			for (size_t j = 1; j <= i; ++j ) {
				cat+=catalans[j-1]*catalans[i-j];
				
			}
			catalans.push_back(cat);
			i++;
			
		}
	}
	
	
	
	
	
}


/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
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
\brief Definitions for function objects for MCMC with adaptive histograms.
*/


#include "histmcmcobjs.hpp"
#include "toolz.hpp"
//#include "sptypes.hpp"
//#include "cxsc.hpp"

#include <stdexcept>
#include <climits>

namespace subpavings {

	/*Abstract class for MCMC priors.
	*/
	
	LogMCMCPrior::LogMCMCPrior() : priorName("") {}

	LogMCMCPrior::LogMCMCPrior(std::string pn) : priorName(pn) {}

	LogMCMCPrior::~LogMCMCPrior() {}
	
	// give the name of the prior
	std::string LogMCMCPrior::getName() const
	{ return priorName; }



	/*A class for a log prior based on a Catalan number prior.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2)
	*/
	
	// default constructor
	LogCatalanPrior::LogCatalanPrior() : LogMCMCPrior("CatalanPrior") 
	{}

	// number of splits k
	cxsc::real LogCatalanPrior::operator()(const size_t k) const
	{
		return -2*lCk(k);

	}
	// number of splits k+1 compared to k
	// (C_k/C_(k+1))^2 = ((k+2)/(2(2k+1)))^2
	cxsc::real LogCatalanPrior::changeOnSplitOne(const size_t k) const
	{
		return 2.0*(std::log(k+2) - std::log(2.0) - std::log(2*k+1));

	}
	// number of splits k-1 compared to k
	// (C_k/C_(k-1))^2 = ((2(2k-1)/(k+1))^2
	cxsc::real LogCatalanPrior::changeOnMergeOne(const size_t k) const
	{
		return 2.0*(std::log(2.0) + std::log(2*k-1) - std::log(k+1));

	}
	



//#if(0)
/*A class for a log prior based on a Catalan number prior with temp term.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2) but also temperature term.
TODO: need to turn this into a temperature as non-negative integer translated prior:
1/(C(k+t)^2), where t \in {0,1,2,...}. So when t=0 we get natural Catalan prior 1/(C(k+t)^2)=1/(Ck^2)
	*/
	
	// default constructor
	LogCatalanTempPrior::LogCatalanTempPrior(double t) : 
    LogMCMCPrior("CatalanPrior"), temp(t) {}

	// number of splits k
	cxsc::real LogCatalanTempPrior::operator()(const size_t k) const
	{
		return -(2+temp)*lCk(k);

	}
	// number of splits k+1 compared to k at "temp" t
	// (C_k/C_(k+1))^(2+t) = ((k+2)/(2(2k+1)))^(2+t)
	cxsc::real LogCatalanTempPrior::changeOnSplitOne(const size_t k) const
	{
		return (2.0+temp)*(std::log(k+2) - std::log(2.0) - std::log(2*k+1));

	}
	// number of splits k-1 compared to k at "temp" t
	// (C_k/C_(k-1))^(2+t) = ((2(2k-1)/(k+1))^(2+t)
	cxsc::real LogCatalanTempPrior::changeOnMergeOne(const size_t k) const
	{
		return (2.0+temp)*(std::log(2.0) + std::log(2*k-1) - std::log(k+1));

	}
//#endif	






	LogTemperaturePrior::LogTemperaturePrior(double t)
				: LogMCMCPrior("TemperaturePrior"), temp(t) {}

	// number of leaves l =- splits k + 1
	cxsc::real LogTemperaturePrior::operator()(const size_t k) const
	{
		//return -1.0/temp*(k+1);//unnormalized prior
		//adding normalizing constant 092015 maths@Stockholm
		return (log(exp(1.0/temp)-1.0))-1.0/temp*(k+1);

	}
	cxsc::real LogTemperaturePrior::changeOnSplitOne(const size_t k) const
	{
		/* e^(-(1/t)*(k+2))/e^(-(1/t)*(k+1))=e^(-1/t)
		 * log(e^(-1/t))=-1/t */
		return -1.0/(temp); //unchanged for normalized version

	}
	cxsc::real LogTemperaturePrior::changeOnMergeOne(const size_t k) const
	{
		/* e^(-(1/t)*(k+1))/e^(-(1/t)*(k+2))=e^(1/t)
		 * log(e^(1/t))=1/t */
		return 1.0/(temp); //unchanged for normalized version

	}


	/*Abstract class for MCMC IMH priors.
	*/
	
	LogMCMCIMHPrior::LogMCMCIMHPrior() : priorName("") {}

	LogMCMCIMHPrior::LogMCMCIMHPrior(std::string pn) : priorName(pn) {}

	LogMCMCIMHPrior::~LogMCMCIMHPrior() {}
	
	// give the name of the prior
	std::string LogMCMCIMHPrior::getName() const
	{ return priorName; }



	/*A class for a log prior for IMH MCMC based on a Catalan number prior.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2)
	*/
	
	// default constructor
	LogCatalanIMHPrior::LogCatalanIMHPrior() : LogMCMCIMHPrior("CatalanPrior") 
	{
		makeCatalans();
	}

	// number of splits k
	cxsc::real LogCatalanIMHPrior::operator()(const size_t k) const
	{
		return -2*lCk(k);

	}
	
		
	/* pair (change, change adjusted for Catalan raio)
	 * 
	 * where change = ln ( (C_kPrime/C_k)^2) 
	 * 				= 2*ln(C_k/C_kPrime) 
	 * 				= 2(ln(C_k) - ln(C_kPrime))
	 * 
	 * and change adjusted for Catalan ratio
	 * 				= ln ( (C_kPrime/C_k)^2 * (C_k/C_kPrime)) 
	 * 				= ln(C_k/C_kPrime) 
	 * 				= ln(C_k) - ln(C_kPrime)
	 * 
	 * ie pair is (2*change, change)
	 * */
	std::pair < cxsc::real, cxsc::real> LogCatalanIMHPrior::logChange(
			unsigned long int k, unsigned long int kPrime) const
	{
		size_t catalanMaxK = catalans.size()-1; // catalans included C0
		
		cxsc::real chg(0.0);
	
		if (k <= catalanMaxK && kPrime <= catalanMaxK) { 
				chg = (cxsc::ln(1.0*catalans[k]) 
					- cxsc::ln(1.0*catalans[kPrime]));
		}
		else if (k <= catalanMaxK && kPrime > catalanMaxK) { 
						
				chg = ( cxsc::ln(1.0*catalans[k]) 
					- (2.0*kPrime*cxsc::ln(2.0) 
						- 0.5*cxsc::ln(cxsc::Pi_real*(1.0*kPrime)) 
						- cxsc::ln(kPrime+1.0)- 1.0/(8*kPrime)) );
		}
		else if (k > catalanMaxK && kPrime <= catalanMaxK) { 
			
				chg = ( (2.0*k*cxsc::ln(2.0) 
					- 0.5*cxsc::ln(cxsc::Pi_real*(1.0*k)) 
					- cxsc::ln(k+1.0)- 1.0/(8*k)) 
						
						- cxsc::ln(1.0*catalans[kPrime]) );
		}
		else  { 
			
			/* don't do (k - kPrime because result should be an unsigned long int
			* and this will go horribly wrong if kPrime > k */
			chg = ( 2.0*k*cxsc::ln(2.0) - 2.0*kPrime*cxsc::ln(2.0) 
				- 0.5*(cxsc::ln(1.0*k) - cxsc::ln(1.0*kPrime)) 
				- (cxsc::ln(k + 1.0) - cxsc::ln(kPrime + 1.0))
				- 1.0/8*(1.0/k - 1.0/kPrime));
		}
		
		return std::pair < cxsc::real, cxsc::real >(2*chg, chg); 
	}


	void LogCatalanIMHPrior::makeCatalans()
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


	

	/* Abstract class for MCMC proposal distributions.*/
	
	MCMCProposal::MCMCProposal() : proposalName("") {}

	MCMCProposal::MCMCProposal(std::string pn) : proposalName(pn) {}
	
	MCMCProposal::~MCMCProposal() {}

	std::string MCMCProposal::getName() const
	{ return proposalName; }







	/* The 'stay-split-merge' base chain proposal class, where 
	probability of staying in same state \f$ \sigma \f$ is fixed,
	the probability of a split is \f$ \frac{1-\sigma}{2} \f$,
	the probability of a merge is \f$ \frac{1-\sigma}{2} \f$ and, given
	a split,
	the probabilities of each leaf being chosen are equal, and given a merge, the
	probabilities of each cherry being chosen are equal.
	*/
	UniformSSMProposal::UniformSSMProposal(double s) 
		: MCMCProposal("UniformSSMProposal"), probSplitMerge((1-s)/2) 
	{
		if (!(s > 0)) throw std::invalid_argument("UniformSSMProposal(double)");
		if (!(s < 1)) throw std::invalid_argument("UniformSSMProposal(double)");
	}

	
	// fill a vector with probabilities as reals
	cxsc::real UniformSSMProposal::fillNodeProposalProbs(
			const size_t nLeaf, const size_t nCherry,
			RealVec& probs) const
	{
		cxsc::real retSum = 0.0;
		probs.reserve(nLeaf + nCherry);

		if (nLeaf > 0) {
			cxsc::real pLeaf = probSplitMerge/nLeaf;
			retSum += (1.0*nLeaf * pLeaf);
			probs.assign(nLeaf, pLeaf);
		}
		if (nCherry > 0) {
			cxsc::real pCherry = probSplitMerge/nCherry;
			retSum += (1.0*nCherry * pCherry);
			probs.insert(probs.end(), nCherry, pCherry);
		}

		return retSum;
	}

	// logQ(m | m') - logQ(m' | m)
	// where m' is proposal state after a split on m
	// what matters is the current number of leaves
	// and the number of cherries under the proposal
	cxsc::real UniformSSMProposal::getLogQRatioSplitProposal(
									const size_t leavesNow,
									const size_t /*cherriesNow*/,
									const size_t /*leavesProspective*/,
									const size_t cherriesProspective) const
	{
		cxsc::real retValue = 0;
		if (leavesNow > 0 && cherriesProspective > 0) {
			 retValue = log(1.0*leavesNow) - log(1.0*cherriesProspective);
		}
		else throw std::logic_error("UniformSSMProposal::getLogQRatioSplitProposal");
		// else retValue = 0 - this should never occur
		return retValue;
	}

	// logQ(m | m') - logQ(m' | m)
	// where m' is proposal state after a split on m
	// what matters is the prospective number of leaves
	// and the current number of cherries
	cxsc::real UniformSSMProposal::getLogQRatioMergeProposal(
									const size_t /*leavesNow*/,
									const size_t cherriesNow,
									const size_t leavesProspective,
									const size_t /*cherriesProspective*/) const
	{
		return log(1.0*cherriesNow) - log(1.0*leavesProspective);
	}





	/* Class where probability of split or merge is fixed and, given a split,
	the probabilities of each leaf being chosen are equal, and given a merge, the
	probabilities of each cherry being chosen are equal.
	*/
	
	// default constructor
	UniformProposal::UniformProposal() : MCMCProposal("UniformSSMProposal"),
							probSplit(0.5) {}

	UniformProposal::UniformProposal(double p) 
		: MCMCProposal("UniformSSMProposal"), probSplit(p) {}

	
	// fill a vector with probabilities as reals
	cxsc::real UniformProposal::fillNodeProposalProbs(
			const size_t nLeaf, const size_t nCherry,
			RealVec& probs) const
	{

		cxsc::real retSum = 0.0;
		probs.reserve(nLeaf + nCherry);

		if (nLeaf > 0) {
			cxsc::real pLeaf = probSplit*1.0/nLeaf;
			retSum += (1.0*nLeaf * pLeaf);
			probs.assign(nLeaf, pLeaf);
		}
		if (nCherry > 0) {
			cxsc::real pCherry = (1.0-probSplit)/nCherry;
			retSum += (1.0*nCherry * pCherry);
			probs.insert(probs.end(), nCherry, pCherry);
		}

		return retSum;
	}

	// logQ(m | m') - logQ(m' | m)
	// where m' is proposal state after a split on m
	// what matters is the current number of leaves
	// and the number of cherries under the proposal
	cxsc::real UniformProposal::getLogQRatioSplitProposal(
									const size_t leavesNow,
									const size_t /*cherriesNow*/,
									const size_t /*leavesProspective*/,
									const size_t cherriesProspective) const
	{
		cxsc::real retValue = 0;
		if (leavesNow == 0 && cherriesProspective > 0) {
			retValue = log(1.0-probSplit) - log(1.0*cherriesProspective);
		}
		else if (leavesNow > 0 && cherriesProspective == 0) {
			retValue = log(probSplit) - log(1.0*leavesNow);
		}
		else if (leavesNow > 0 && cherriesProspective > 0) {
			 retValue = log(1.0-probSplit) - log(probSplit)
						+ log(1.0*leavesNow) - log(1.0*cherriesProspective);
		}
		// else retValue = 0 - this should never occur
		return retValue;
	}

	// logQ(m | m') - logQ(m' | m)
	// where m' is proposal state after a split on m
	// what matters is the prospective number of leaves
	// and the current number of cherries
	cxsc::real UniformProposal::getLogQRatioMergeProposal(
									const size_t /*leavesNow*/,
									const size_t cherriesNow,
									const size_t leavesProspective,
									const size_t /*cherriesProspective*/) const
	{
		return log(probSplit) - log(1.0-probSplit)
					+ log(1.0*cherriesNow) - log(1.0*leavesProspective);
	}




	/* Class where probabilities of any splittable leaf or mergeable cherry being
	proposed are equal, i.e. if there is just one leaf (eg root) then it is certain
	to be proposed.
	*/
	EquiProbProposal::EquiProbProposal() : MCMCProposal("EquiprobableProposal") {}


	// fill a vector with probabilities as reals
	cxsc::real EquiProbProposal::EquiProbProposal::fillNodeProposalProbs(
									const size_t nLeaf,
									const size_t nCherry,
									RealVec& probs) const
	{
		cxsc::real retSum = 0.0;
		probs.reserve(nLeaf + nCherry);

		if (nLeaf + nCherry > 0) {
			cxsc::real pNode = 1.0/(nLeaf + nCherry);
			retSum += pNode*(1.0*nLeaf + 1.0*nCherry);

			probs.assign(nLeaf+nCherry, pNode);
		}

		return retSum;
	}

	// logQ(m | m') - logQ(m' | m)
	// where m' is proposal state after a split on m
	// the current and prospective numbers of leaves and cherries all matter
	cxsc::real EquiProbProposal::getLogQRatioSplitProposal(
									const size_t leavesNow,
									const size_t cherriesNow,
									const size_t leavesProspective,
									const size_t cherriesProspective) const
	{
		cxsc::real retValue = 0;
		if ((leavesNow + cherriesNow > 0)
				&& (leavesProspective + cherriesProspective > 0)) {
			retValue = log(1.0*(leavesNow + cherriesNow))
				- log(1.0*(leavesProspective + cherriesProspective));
		}
		// else retValue = 0 - this should never occur
		return retValue;
	}

	// same as for split
	cxsc::real EquiProbProposal::getLogQRatioMergeProposal(
									const size_t leavesNow,
									const size_t cherriesNow,
									const size_t leavesProspective,
									const size_t cherriesProspective) const
	{
		return getLogQRatioSplitProposal(leavesNow, cherriesNow,
								leavesProspective, cherriesProspective);
	}
} // end namespace subpavings




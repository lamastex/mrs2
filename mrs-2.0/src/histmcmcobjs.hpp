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
\brief Declarations for function objects for MCMC with adaptive histograms.
*/

#ifndef ___MCMCFOBJS_HPP__
#define ___MCMCFOBJS_HPP__


//#include "toolz.hpp"
#include "sptypes.hpp"
#include "cxsc.hpp"

#include <string>
#include <vector>
#include <utility>

namespace subpavings {

	/*! \brief Abstract class for MCMC priors.
	*/
	class LogMCMCPrior
	{
		
		public:

		/*! \brief No-args constructor. */
		LogMCMCPrior();

		/*! \brief Constructor. */
		LogMCMCPrior(std::string pn);
		
		virtual ~LogMCMCPrior();

		/*! \brief Get the name of the prior type */
		std::string getName() const;

		/*! \brief Get log prior value for state with \a k splits.*/
		virtual cxsc::real operator()(const size_t k) const = 0;
		
		/*! \brief Get the change in log prior value when a 
		 * state with \a k splits is split once more.*/
		virtual cxsc::real changeOnSplitOne(const size_t k) const = 0;
		
		/*! \brief Get the change in log prior value when a 
		 * state with \a k splits has one cherry merged.*/
		virtual cxsc::real changeOnMergeOne(const size_t k) const = 0;
		
		
		protected:

		std::string priorName;


	};

	/*! \brief A class for a log prior based on a Catalan number prior.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2)
	*/
	class LogCatalanPrior : public LogMCMCPrior
	{
		public:

		// default constructor
		LogCatalanPrior();
		
		~LogCatalanPrior() {};

		// number of splits k
		cxsc::real operator()(const size_t k) const;
		
		// number of splits k+1 compared to k
		// (C_k/C_(k+1))^2 = ((k+2)/(2(2k+1)))^2
		cxsc::real changeOnSplitOne(const size_t k) const;
		
		// number of splits k-1 compared to k
		// (C_k/C_(k-1))^2 = ((2(2k-1)/(k+1))^2
		cxsc::real changeOnMergeOne(const size_t k) const;
		
		
		private:
		
		std:: vector < unsigned long int > catalans;
		

	};

//#if(0)    
    /*! \brief A class for a log prior based on a Catalan number prior
     * with a temperature term.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2) but also temperature term.
	*/
	class LogCatalanTempPrior : public LogMCMCPrior
	{
		public:

		// constructor
		LogCatalanTempPrior(double t);
		
		~LogCatalanTempPrior() {};

		// number of splits k
		cxsc::real operator()(const size_t k) const;
		
		// number of splits k+1 compared to k
		// (C_k/C_(k+1))^2 = ((k+2)/(2(2k+1)))^2
		cxsc::real changeOnSplitOne(const size_t k) const;
		
		// number of splits k-1 compared to k
		// (C_k/C_(k-1))^2 = ((2(2k-1)/(k+1))^2
		cxsc::real changeOnMergeOne(const size_t k) const;
		
		
		private:
		
        double temp; // the temperature coefficient
        
		std:: vector < unsigned long int > catalans;
		
        LogCatalanTempPrior();
	

	};

//#endif	

	class LogTemperaturePrior : public LogMCMCPrior
	{
		
		public:

		explicit LogTemperaturePrior(double t);
		
		~LogTemperaturePrior() {};

		// number of leaves l = splits + 1
		cxsc::real operator()(const size_t k) const;
		
		cxsc::real changeOnSplitOne(const size_t k) const;
		
		cxsc::real changeOnMergeOne(const size_t k) const;
		
				
		private:

		double temp; // the temperature coefficient

		// default constructor is private and cannot be used outside the class
		LogTemperaturePrior() {}


	};

	/*! \brief Abstract class for IMH MCMC priors.
	*/
	class LogMCMCIMHPrior
	{
		
		public:

		/*! \brief No-args constructor. */
		LogMCMCIMHPrior();

		/*! \brief Constructor. */
		LogMCMCIMHPrior(std::string pn);
		
		virtual ~LogMCMCIMHPrior();

		/*! \brief Get the name of the prior type */
		std::string getName() const;

		/*! \brief Get log prior value for state with \a k splits.*/
		virtual cxsc::real operator()(const size_t k) const = 0;
		
		/*! \brief Get the change in log prior value when a 
		 * state with \a k splits changes to a state with \a kPrime
		 * splits, and the same adjusted for catalan ratio C_k/C_kPrime.*/
		virtual std::pair < cxsc::real, cxsc::real> logChange(
			unsigned long int k, unsigned long int kPrime) const = 0;
		
		
		protected:

		std::string priorName;


	};

	/*! \brief A class for a log prior for IMH MCMC based on a Catalan number prior.

	The prior is related to the Catalan number of the number of bisections of
	the root box k, Ck.  After k splits there are Ck distinct possible full
	binary trees and so, for each k, total probability of all the states
	resulting from k splits is proportional to 1/Ck.  If each of these
	k-split states is equally likely then then the probability of any one
	k-split state is 1/(Ck^2)
	*/
	class LogCatalanIMHPrior : public LogMCMCIMHPrior
	{
		public:

		// default constructor
		LogCatalanIMHPrior();
		
		~LogCatalanIMHPrior() {};

		// number of splits k
		cxsc::real operator()(const size_t k) const;
		
		std::pair < cxsc::real, cxsc::real> logChange(
			unsigned long int k, unsigned long int kPrime) const;

		private:
		
		void makeCatalans();
		
		std:: vector < unsigned long int > catalans;
		

	};



	/*! \brief Abstract class for MCMC proposal distributions.

	The proposal distribution function object cannot work directly or solely with
	the state of the node tree because the MCMC definition allows some leaf nodes
	to be excluded from those proposable for splitting.  Hence the proposal
	distributions should be based on the splittable leaf and cherry nodes, not all
	the leaf and cherry nodes in the tree.
	*/
	class MCMCProposal
	{
		public:

		MCMCProposal();

		explicit MCMCProposal(std::string pn);
		
		virtual ~MCMCProposal();

		std::string getName() const;

		// fill a vector with probabilities as reals
		// and return the sum of the probabilities
		virtual cxsc::real fillNodeProposalProbs(const size_t nLeaf,
											const size_t nCherry,
											RealVec& probs) const = 0;

		// logQ(m | m') - logQ(m' | m) for split proposal m'
		virtual cxsc::real getLogQRatioSplitProposal(const size_t leavesNow,
												const size_t cherriesNow,
												const size_t leavesProspective,
												const size_t cherriesProspective
												) const = 0;
		// logQ(m | m') - logQ(m' | m) for merge proposal m'
		virtual cxsc::real getLogQRatioMergeProposal(const size_t leavesNow,
												const size_t cherriesNow,
												const size_t leavesProspective,
												const size_t cherriesProspective
												) const = 0;
		
		protected:

			std::string proposalName;

		
	};

	/*! Class where probability of split or merge is fixed and, given a split,
	the probabilities of each leaf being chosen are equal, and given a merge, the
	probabilities of each cherry being chosen are equal.
	*/
	class UniformProposal : public MCMCProposal
	{

		public:

			// default constructor
			UniformProposal();
			
			explicit UniformProposal(double p);

			~UniformProposal(){};
			
			// fill a vector with probabilities as reals
			cxsc::real fillNodeProposalProbs(const size_t nLeaf, const size_t nCherry,
											RealVec& probs) const;

			// logQ(m | m') - logQ(m' | m)
			// where m' is proposal state after a split on m
			// what matters is the current number of leaves
			// and the number of cherries under the proposal
			cxsc::real getLogQRatioSplitProposal(const size_t leavesNow,
											const size_t /*cherriesNow*/,
											const size_t /*leavesProspective*/,
											const size_t cherriesProspective) const;
											
			// logQ(m | m') - logQ(m' | m)
			// where m' is proposal state after a split on m
			// what matters is the prospective number of leaves
			// and the current number of cherries
			cxsc::real getLogQRatioMergeProposal(const size_t /*leavesNow*/,
											const size_t cherriesNow,
											const size_t leavesProspective,
											const size_t /*cherriesProspective*/) const;
			
			private:

			double probSplit;

		
	};
	
	/*! \brief The 'stay-split-merge' base chain proposal class, where 
	probability of staying in same state \f$ \sigma \f$ is fixed,
	the probability of a split is \f$ \frac{1-\sigma}{2} \f$,
	the probability of a merge is \f$ \frac{1-\sigma}{2} \f$ and, given
	a split,
	the probabilities of each leaf being chosen are equal, and given a merge, the
	probabilities of each cherry being chosen are equal.
	*/
	class UniformSSMProposal : public MCMCProposal
	{

		public:

			
			explicit UniformSSMProposal(double s);

			~UniformSSMProposal(){};
			
			// fill a vector with probabilities as reals
			cxsc::real fillNodeProposalProbs(const size_t nLeaf, const size_t nCherry,
											RealVec& probs) const;

			// logQ(m | m') - logQ(m' | m)
			// where m' is proposal state after a split on m
			// what matters is the current number of leaves
			// and the number of cherries under the proposal
			cxsc::real getLogQRatioSplitProposal(const size_t leavesNow,
											const size_t /*cherriesNow*/,
											const size_t /*leavesProspective*/,
											const size_t cherriesProspective) const;
											
			// logQ(m | m') - logQ(m' | m)
			// where m' is proposal state after a split on m
			// what matters is the prospective number of leaves
			// and the current number of cherries
			cxsc::real getLogQRatioMergeProposal(const size_t /*leavesNow*/,
											const size_t cherriesNow,
											const size_t leavesProspective,
											const size_t /*cherriesProspective*/) const;
			
			private:

			UniformSSMProposal();
			
			double probSplitMerge;

		
	};

	/*! Class where probabilities of any splittable leaf or mergeable cherry being
	proposed are equal, i.e. if there is just one leaf (eg root) then it is certain
	to be proposed.
	*/
	class EquiProbProposal : public MCMCProposal
	{

		public:

			EquiProbProposal();
			
			~EquiProbProposal(){};

			// fill a vector with probabilities as reals
			cxsc::real fillNodeProposalProbs(const size_t nLeaf, const size_t nCherry,
											RealVec& probs) const;

			// logQ(m | m') - logQ(m' | m)
			// where m' is proposal state after a split on m
			// the current and prospective numbers of leaves and cherries all matter
			cxsc::real getLogQRatioSplitProposal(const size_t leavesNow,
											const size_t cherriesNow,
											const size_t leavesProspective,
											const size_t cherriesProspective) const;
			// same as for split
			cxsc::real getLogQRatioMergeProposal(const size_t leavesNow,
											const size_t cherriesNow,
											const size_t leavesProspective,
											const size_t cherriesProspective) const;
	};
} // end namespace subpavings

#endif



/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009 Jennifer Harlow
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

/*! \file      histmcmcobjs.hpp
\brief Function objects for MCMC with adaptive histograms.
*/

#ifndef ___MCMCFOBJS_HPP__
#define ___MCMCFOBJS_HPP__

namespace subpavings {

/*! \brief Abstract class for MCMC priors.
*/
class LogMCMCPrior
{
    protected:

    string priorName;

    public:

    LogMCMCPrior() : priorName("") {}

    LogMCMCPrior(string pn) : priorName(pn) {}

    virtual string getName() const = 0;

    virtual real operator()(const size_t k) const = 0;

};

/*! \brief A class for a log prior based on a Catalan number prior.

The prior is related to the Catalan number of the number of bisections of
the root box k, Ck.  After k splits there are Ck distinct possible full
binary trees and so, for each k, total probability of all the states
resulting from k splits is proportional to 1/Ck.  If if each of these
k-split states is equally likely then then the probability of any one
k-split state is 1/(Ck^2)
*/
class LogCatalanPrior : public LogMCMCPrior
{
    private:
    
    double temp;
   
    
    public:

    // default constructor
    LogCatalanPrior() {}; 
     
    //LogCatalanPrior(double t) : LogMCMCPrior("CatalanPrior"), temp(t) {}

    // give the name of the prior
    string getName() const
    { return priorName; }

    // number of splits k
    // but i think input is number of leaves - so need to subtract 1?
    real operator()(const size_t k) const
    {
			return -2*lCk(k)/1;
			//double a = 2  + 4*M_PI/pow(3, 2.5); 
			//return (-1*log(a) - 2*lCk(k-1))/temp;
    }

};

class LogTemperaturePrior : public LogMCMCPrior
{
    private:

    double temp; // the temperature coefficient

    // default constructor is private and cannot be used outside the class
    LogTemperaturePrior() {}

    public:

    LogTemperaturePrior(double t) : LogMCMCPrior("CatalanPrior"), temp(t) {}

        // give the name of the prior
    string getName() const
    { return priorName; }

    // number of leaves l
    real operator()(const size_t l) const
    {
        return -1.0/temp*l;

    }

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
    protected:

        string proposalName;

    public:

    MCMCProposal() : proposalName("") {}

    MCMCProposal(string pn) : proposalName(pn) {}

    virtual string getName() const = 0;

    // fill a vector with probabilities as reals
    // and return the sum of the probabilities
    virtual real fillNodeProposalProbs(const size_t nLeaf,
                                        const size_t nCherry,
                                        RealVec& probs) const = 0;

    // logQ(m | m') - logQ(m' | m) for split proposal m'
    virtual real getLogQRatioSplitProposal(const size_t leavesNow,
                                            const size_t cherriesNow,
                                            const size_t leavesProspective,
                                            const size_t cherriesProspective
                                            ) const = 0;
    // logQ(m | m') - logQ(m' | m) for merge proposal m'
    virtual real getLogQRatioMergeProposal(const size_t leavesNow,
                                            const size_t cherriesNow,
                                            const size_t leavesProspective,
                                            const size_t cherriesProspective
                                            ) const = 0;
};

/*! Class where probability of split or merge is fixed and, given a split,
the probabilities of each leaf being chosen are equal, and given a merge, the
probabilities of each cherry being chosen are equal.
*/
class UniformProposal : public MCMCProposal
{

    private:

        double probSplit;

    public:

        // default constructor
        UniformProposal() : MCMCProposal("UniformProposal"),
                                probSplit(0.5) {}

        UniformProposal(double p) : MCMCProposal("UniformProposal"),
                                        probSplit(p) {}

        // give the name of the proposal
        string getName() const
        { return proposalName; }

        // fill a vector with probabilities as reals
        real fillNodeProposalProbs(const size_t nLeaf, const size_t nCherry,
                                        RealVec& probs) const
        {

            real retSum = 0.0;
            probs.reserve(nLeaf + nCherry);

            if (nLeaf > 0) {
                real pLeaf = probSplit*1.0/nLeaf;
                retSum += (1.0*nLeaf * pLeaf);
                probs.assign(nLeaf, pLeaf);
            }
            if (nCherry > 0) {
                real pCherry = (1.0-probSplit)/nCherry;
                retSum += (1.0*nCherry * pCherry);
                probs.insert(probs.end(), nCherry, pCherry);
            }

            return retSum;
        }

        // logQ(m | m') - logQ(m' | m)
        // where m' is proposal state after a split on m
        // what matters is the current number of leaves
        // and the number of cherries under the proposal
        real getLogQRatioSplitProposal(const size_t leavesNow,
                                        const size_t cherriesNow,
                                        const size_t leavesProspective,
                                        const size_t cherriesProspective) const
        {
            real retValue = 0;
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
        real getLogQRatioMergeProposal(const size_t leavesNow,
                                        const size_t cherriesNow,
                                        const size_t leavesProspective,
                                        const size_t cherriesProspective) const
        {
            return log(probSplit) - log(1.0-probSplit)
                        + log(1.0*cherriesNow) - log(1.0*leavesProspective);
        }
};

/*! Class where probabilities of any splittable leaf or mergeable cherry being
proposed are equal, i.e. if there is just one leaf (eg root) then it is certain
to be proposed.
*/
class EquiProbProposal : public MCMCProposal
{

    public:

        EquiProbProposal() : MCMCProposal("EquiprobableProposal") {}

        // give the name of the proposal
        string getName() const
        { return proposalName; }

        // fill a vector with probabilities as reals
        real fillNodeProposalProbs(const size_t nLeaf, const size_t nCherry,
                                        RealVec& probs) const
        {
            real retSum = 0.0;
            probs.reserve(nLeaf + nCherry);

            if (nLeaf + nCherry > 0) {
                real pNode = 1.0/(nLeaf + nCherry);
                retSum += pNode*(1.0*nLeaf + 1.0*nCherry);

                probs.assign(nLeaf+nCherry, pNode);
            }

            return retSum;
        }

        // logQ(m | m') - logQ(m' | m)
        // where m' is proposal state after a split on m
        // the current and prospective numbers of leaves and cherries all matter
        real getLogQRatioSplitProposal(const size_t leavesNow,
                                        const size_t cherriesNow,
                                        const size_t leavesProspective,
                                        const size_t cherriesProspective) const
        {
            real retValue = 0;
            if ((leavesNow + cherriesNow > 0)
                    && (leavesProspective + cherriesProspective > 0)) {
                retValue = log(1.0*(leavesNow + cherriesNow))
                    - log(1.0*(leavesProspective + cherriesProspective));
            }
            // else retValue = 0 - this should never occur
            return retValue;
        }

        // same as for split
        real getLogQRatioMergeProposal(const size_t leavesNow,
                                        const size_t cherriesNow,
                                        const size_t leavesProspective,
                                        const size_t cherriesProspective) const
        {
            return getLogQRatioSplitProposal(leavesNow, cherriesNow,
                                    leavesProspective, cherriesProspective);
        }
};

} // end of namespace subpavings

#endif



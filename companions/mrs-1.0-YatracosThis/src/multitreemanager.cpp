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

/*! \file multitreemanager.cpp
\brief MultiTreeManager definitions
*/

#include "multitreemanager.hpp"

#include <gsl/gsl_rng.h>    // to use the gsl random number generator
#include <map>              // to use maps
#include <sstream>  // to be able to manipulate strings as streams

#include "sptools.hpp"      // to use the subpaving tools

/*! templatized function object for lexicographical sorting of vectors whose elements have total ordering
*/
template <class T>
class LexicoSorting
{
  public:
    bool operator() (const T& t1, const T& t2) const {
      return std::lexicographical_compare(&t1[0], &t1[t1.size()-1], &t2[0], &t2[t2.size()-1]);
      //return lexicographical_compare(t1.begin(), t1.end(), t2.begin(), t2.end());
    }
};


// a map for counting visits to states of Histograms encoded
    // as OrderedLeafDepths
    typedef map< IntVec, int, LexicoSorting<IntVec > > OrdLeafDepthsMap;


// ---------- implementation of MultiTreeManager class -------------

// ---------------- private methods


// recursively accumulate splitting outcomes
bool MultiTreeManager::addToOutcomeSpaceAndGraph(const string s, int toLevel,
            int thisLevel, SPSnode * tree, set<string>& lines)
{
    bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = spLeaves(tree);

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space
        string parent = "\"" + getLeafLevelsString(tree) + "\"";

        // make one copy of the tree for each leaf
        // put each copy in the temporary container
        for (size_t i = 0; i < n; i++) {
            SPSnode* copyTree = new SPSnode(*tree);
            container.push_back(copyTree);
            // take the jth copy of tree and find its leaves
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();
            // add the parent child connection to the dot graph
            string segment = "\"" + getLeafLevelsString(copyTree) + "\"";
            string line = "\t " + parent + " -> " + segment + ";";

            pair<set<string>::iterator, bool> ret;

            ret = lines.insert(line); // try and see if the line is a new one
            if (ret.second==true) outputFile(s, line);  // output if new line
        }

        // copy the contents of the temporary container to pavings
        pavings.insert(pavings.end(), container.begin(),container.end());

        SPSnodePtrsItr it;
        // for each copy in the temporary container

        for (it = container.begin(); it < container.end(); it++) {
            // recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
            done = addToOutcomeSpaceAndGraph(s, toLevel,
                                            thisLevel + 1, *it, lines);
        }

    }
    else {
        done = true;
    }
    return done;

}

// recursively accumulate splitting outcomes
bool MultiTreeManager::addToOutcomeSpace(int toLevel,
            int thisLevel, SPSnode * tree)
{
    bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = spLeaves(tree);

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space

        // make one copy of the tree for each leaf
        // put each copy in the temporary container
        for (size_t i = 0; i < n; i++) {

            SPSnode* copyTree = new SPSnode(*tree);
            container.push_back(copyTree);
            // take the jth copy of tree and find its leaves
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();

        }

        // copy the contents of the temporary container to pavings
        pavings.insert(pavings.end(), container.begin(),container.end());

        SPSnodePtrsItr it;

        // for each copy in the temporary container
        for (it = container.begin(); it < container.end(); it++) {
            // recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
            done = addToOutcomeSpace(toLevel, thisLevel + 1, *it);
        }

    }
    else {
        done = true;
    }
    return done;

}

// ---------------- public methods


// Destructor
MultiTreeManager::~MultiTreeManager()
{
    SPSnodePtrsItr it;
    for (it = pavings.begin(); it < pavings.end(); it++) {
        delete *it;
        *it = NULL;
    }
    pavings.clear();
}


// recursively accumulate splitting outcomes
void MultiTreeManager::mapPavings()
{
    // don't clear the current pavings! - that's what we use here

    // go through the pavings and record into the vector of shape count maps
    if (!pavings.empty()) {

         size_t numberPavings = pavings.size();

        cout << " mapping outcomes... there are " << numberPavings
                        << " outcomes " << endl;

        std::pair< OrdLeafDepthsMap::iterator, bool > mapBool;

        // what is the maximum number of leaves we have?
        set<size_t> leafSet;
        size_t count = 0;
        for (SPSnodePtrsItr sit = pavings.begin(); sit < pavings.end(); sit++) {
            leafSet.insert(spLeaves(*sit));
            count ++;
        }

        size_t uniqueLeaves = leafSet.size();

        vector< OrdLeafDepthsMap > vecMaps(uniqueLeaves); // vector of maps

        for (SPSnodePtrsItr it = pavings.begin(); it < pavings.end(); it++) {

            size_t thisLeaves = spLeaves((*it));

            IntVec thisLevels;
            thisLevels  = (*it)->getLeafNodeLevels(thisLevels);

            int indexer = 0; // need to find index for maps with this no. leaves
            for (int j = 1; j < thisLeaves; j++) {
                if (leafSet.count(j)) indexer++;
            }

            mapBool = vecMaps[indexer].insert(pair<IntVec, int> (thisLevels,1));

            if(!(mapBool.second)) // if its a new one, add
                (mapBool.first)->second +=1; // else increment count
        }

        string basefilename = "multimanager";
        string filename = getUniqueFilename (basefilename);

        outputFile(filename,
            "tree shape frequencies and relative frequencies in outcomes");
        std::ostringstream stm;
        stm << numberPavings << " outcomes altogether";
        outputFile(filename, stm.str());

        // print out the results
        std::cout << "tree shape frequencies (relative frequencies) in outcomes"
                << "\n";

        vector< OrdLeafDepthsMap >::iterator vecMapsIt = vecMaps.begin();

        for (set<size_t>::iterator lit = leafSet.begin(); lit != leafSet.end();
                                        lit++) {

            OrdLeafDepthsMap::iterator mapIt;

            cout << "Leaves : " << "\t" << *lit << "\t("
                            << (*vecMapsIt).size() << " unique outcomes)\n";

            std::ostringstream stm1;
            stm1 << "Leaves : \t" << (*lit) << "\t" << (*vecMapsIt).size()
                            << "\tunique outcomes";
            outputFile(filename, stm1.str());

            for(mapIt = (*vecMapsIt).begin(); mapIt != (*vecMapsIt).end();
                                        mapIt++) {

                IntVec levels = mapIt->first;

                size_t freq = (mapIt->second);
                double relfreq = (1.0*(mapIt->second))/numberPavings;

                cout << freq << " : " << "(\t"
                    << relfreq << ")\t";

                std::ostringstream stm2;
                stm2 << "\t\t\t\t" << freq << "\t"
                    << relfreq << "\t";
                string thisline =  stm2.str();

                outputFile(filename, thisline, levels);

                copy (levels.begin(), levels.end(),
                      ostream_iterator<int>(cout, ";"));
                cout << '\n';
            }
            vecMapsIt++; // move the iterator to the vec maps in step
        }
        std::cout << "Output file in " << filename << std::endl;
    }
    else std::cout << "There are no pavings in the outcome space to map"
                << std::endl;

}


// recursively accumulate splitting outcomes
void MultiTreeManager::makeAndMapOutcomeSpace(int toLevel)
{
    makeOutcomeSpace(toLevel);
    mapPavings();

}


// get the outcome space from continual splitting to level toLevel
// and do the graph
void MultiTreeManager::makeAndGraphOutcomeSpace(int toLevel)
{

    // clear the pavings
    pavings.clear();

    std::cout << "Making and graphing the outcome space" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    int i = 0;
    string baseFileName = "outputGraph";
    string suffix = ".dot";
    string s = getUniqueFilename(baseFileName, suffix);
    outputFile(s, "digraph G {"); // opening line

    // add it to the pavings data member collection
    pavings.push_back(root);

    // Send it into addtoOutComeSpace
    bool done = false;
    set<string> lines; // to check on uniqueness of lines for graph
    done = addToOutcomeSpaceAndGraph(s, toLevel, 1, root, lines);
    while (!done) {}

    outputFile(s, "}"); // closing line

    // make the image of the graph
    makeDotImage(s);

}

// get the outcome space from continual splitting to level toLevel
void MultiTreeManager::makeOutcomeSpace(int toLevel)
{

    // clear the pavings
    pavings.clear();

    std::cout << "Making the outcome space" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    int i = 0;
    // add it to the pavings data member collection
    pavings.push_back(root);

    // Send it into addtoOutComeSpace
    bool done = false;
    done = addToOutcomeSpace(toLevel, 1, root);
    while (!done) {}

}

// ----------- end of implementation of MultiTreeManager class ---------------



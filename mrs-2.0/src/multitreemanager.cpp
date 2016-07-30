/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief MultiTreeManager definitions
*/

#include "multitreemanager.hpp"

#include "sptools.hpp"      // to use the subpaving tools
#include "sptypes.hpp"      // to use the subpaving tools

#include <map>              // to use maps
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <sstream>
#include <iterator>
#include <cassert>
#include <stdexcept>

//#define MYDEBUG



namespace subpavings {
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
    typedef std::map< IntVec, int, LexicoSorting<IntVec > > OrdLeafDepthsMap;


	
	
}

using namespace subpavings;
using namespace std;

MultiTreeManager::VertexInfo::VertexInfo(const IntVec& _sl, 
										const IntVec& _cl,
										int _s,
										const std::string& _l)
		: sortedLevels(_sl), codedLevels(_cl), sum(_s), label(_l), 
			constraint_line("")
{}

bool MultiTreeManager::VertexInfo::operator< (
					const MultiTreeManager::VertexInfo& rhs) const
{
	return lexicographical_compare(codedLevels.begin(),
									codedLevels.end(),
									rhs.codedLevels.begin(),
									rhs.codedLevels.end());
}
// ---------- implementation of MultiTreeManager class -------------


// ---------------- public methods


// Destructor
MultiTreeManager::~MultiTreeManager()
{
    SPSnodePtrsItr it;
    for (it = pavings.begin(); it < pavings.end(); ++it) {
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
        for (SPSnodePtrsItr sit = pavings.begin(); sit < pavings.end(); ++sit) {
            leafSet.insert((*sit)->getNumberLeaves());
            count ++;
        }

        size_t uniqueLeaves = leafSet.size();

        vector< OrdLeafDepthsMap > vecMaps(uniqueLeaves); // vector of maps

        for (SPSnodePtrsItr it = pavings.begin(); it < pavings.end(); ++it) {

            size_t thisLeaves = (*it)->getNumberLeaves();

            IntVec thisLevels;
            thisLevels  = (*it)->getLeafNodeLevels(thisLevels);

            int indexer = 0; // need to find index for maps with this no. leaves
            for (size_t j = 1; j < thisLeaves; ++j) {
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
                                        ++lit) {

            OrdLeafDepthsMap::iterator mapIt;

            cout << "Leaves : " << "\t" << *lit << "\t("
                            << (*vecMapsIt).size() << " unique outcomes)\n";

            std::ostringstream stm1;
            stm1 << "Leaves : \t" << (*lit) << "\t" << (*vecMapsIt).size()
                            << "\tunique outcomes";
            outputFile(filename, stm1.str());

            for(mapIt = (*vecMapsIt).begin(); mapIt != (*vecMapsIt).end();
                                        ++mapIt) {

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
void MultiTreeManager::makeAndMapOutcomeSpace(int toLevel, int maxDepth)
{
	if (maxDepth == 0) maxDepth = toLevel;
	
    makeOutcomeSpace(toLevel, maxDepth);
    mapPavings();

}


// get the outcome space from continual splitting to level toLevel
// and do the graph
void MultiTreeManager::makeAndGraphOutcomeSpace(int toLevel, int maxDepth)
{

	if (maxDepth == 0) maxDepth = toLevel;
	
    // clear the pavings
    pavings.clear();

    std::cout << "Making and graphing the outcome space" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    string baseFileName = "outputGraph";
    string suffix = ".dot";
    string s = getUniqueFilename(baseFileName, suffix);
    outputFile(s, "digraph G {"); // opening line

    // add it to the pavings data member collection
    pavings.push_back(root);

    // Send it into addtoOutComeSpace
    bool done = false;
    set<string> lines; // to check on uniqueness of lines for graph
    done = addToOutcomeSpaceAndGraph(s, toLevel, 1, root, maxDepth, lines);
    while (!done) {}

    outputFile(s, "}"); // closing line

    // make the image of the graph
    makeDotImage(s);

}

// get the outcome space from continual splitting to level toLevel
void MultiTreeManager::makeOutcomeSpace(int toLevel, int maxDepth)
{

	if (maxDepth == 0) maxDepth = toLevel;
	
    // clear the pavings
    pavings.clear();

    std::cout << "Making the outcome space" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    // add it to the pavings data member collection
    pavings.push_back(root);

    // Send it into addtoOutComeSpace
    bool done = false;
    done = addToOutcomeSpace(toLevel, 1, root, maxDepth);
    

}

// get the outcome space from continual splitting to level toLevel
void MultiTreeManager::makeFinalOutcomeSpace(int toLevel, int maxDepth)
{

	if (maxDepth == 0) maxDepth = toLevel;
	
    // clear the pavings
    pavings.clear();

    std::cout << "Making the outcome space" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    // add it to the pavings data member collection
    pavings.push_back(root);

    // Send it into addtoOutComeSpace
    bool done = false;
    done = getFinalLevelOutcomeSpace(toLevel, 1, root, maxDepth);
    

}


MultiTreeManager::level_vertex_info_map 
						MultiTreeManager::makeLevelStringSumCountsMap(int toLevel)
{
	MultiTreeManager::level_vertex_info_map mp;
	
	// clear the pavings
    pavings.clear();

    std::cout << "Starting to propagate trees" << std::endl;

    //Make a node with a dummy box
    ivector pavingBox(1);
    interval pavingInterval(0,1);
    pavingBox[1] = pavingInterval;
    SPSnode* root = new SPSnode(pavingBox);

    // add it to the pavings data member collection
    pavings.push_back(root);
	
	std::pair < MultiTreeManager::vertex_info_set::iterator, bool >
		addPair = addTreeToLevelStringSumCountsMap(mp, root);

	if (!addPair.second) throw std::logic_error ("Failed to put root into map");
	
	// Send it into addtoOutComeSpace
    addToLevelStringSumCountsMap(mp, toLevel, 1, root, addPair.first);
    
	std::cout << "done map" << std::endl;
	
	return mp;
	
}


// ---------------- private methods


// recursively accumulate splitting outcomes
bool MultiTreeManager::addToOutcomeSpaceAndGraph(const string& s, int toLevel,
            int thisLevel, SPSnode * tree, int maxDepth, set<string>& lines)
{
    bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = tree->getNumberLeaves();

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space
        string parent = "\"" + getLeafLevelsString(tree) + "\"";

        // make one copy of the tree for each leaf
        // put each copy in the temporary container
        for (size_t i = 0; i < n; ++i) {
            SPSnode* copyTree = new SPSnode(*tree);
            // take the jth copy of tree and find its leaves
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();
            
			if ( copyTree->getTreeHeight() <= maxDepth ) {
				container.push_back(copyTree);
				
				// add the parent child connection to the dot graph
				string segment = "\"" + getLeafLevelsString(copyTree) + "\"";
				string line = "\t " + parent + " -> " + segment + ";";

				pair<set<string>::iterator, bool> ret;

				ret = lines.insert(line); // try and see if the line is a new one
				if (ret.second==true) outputFile(s, line);  // output if new line

			}
			
        }

        // copy the contents of the temporary container to pavings
        pavings.insert(pavings.end(), container.begin(),container.end());

        SPSnodePtrsItr it;
        // for each copy in the temporary container

        for (it = container.begin(); it < container.end(); ++it) {
            // recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
            done = addToOutcomeSpaceAndGraph(s, toLevel,
                                            thisLevel + 1, *it, maxDepth, lines);
        }

    }
    else {
        done = true;
    }
    return done;

}

// recursively accumulate splitting outcomes
bool MultiTreeManager::addToOutcomeSpace(int toLevel,
            int thisLevel, SPSnode * tree, int maxDepth)
{
    bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = tree->getNumberLeaves();

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space

        // make one copy of the tree for each leaf
        // put each copy in the temporary container
        for (size_t i = 0; i < n; ++i) {

		    SPSnode* copyTree = new SPSnode(*tree);
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();
			if ( copyTree->getTreeHeight() <= maxDepth ) {
				container.push_back(copyTree);

			}
        }

        // copy the contents of the temporary container to pavings
        pavings.insert(pavings.end(), container.begin(),container.end());

        SPSnodePtrsItr it;

        // for each copy in the temporary container
        for (it = container.begin(); it < container.end(); ++it) {
            // recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
            done = addToOutcomeSpace(toLevel, thisLevel + 1, *it, maxDepth);
        }

    }
    else {
        done = true;
    }
    return done;

}


// recursively accumulate splitting outcomes
bool MultiTreeManager::getFinalLevelOutcomeSpace(int toLevel,
            int thisLevel, SPSnode * tree, int maxDepth)
{
    bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = tree->getNumberLeaves();

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space

        // make one copy of the tree for each leaf
        // put each copy in the temporary container
        for (size_t i = 0; i < n; ++i) {

            SPSnode* copyTree = new SPSnode(*tree);
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();
			if ( copyTree->getTreeHeight() <= maxDepth ) {
				container.push_back(copyTree);

			}

        }

        // copy the contents of the temporary container to pavings
        if (thisLevel == toLevel) 
			pavings.insert(pavings.end(), container.begin(),container.end());

		else {
			SPSnodePtrsItr it;

			// for each copy in the temporary container
			for (it = container.begin(); it < container.end(); ++it) {
				// recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
				done = getFinalLevelOutcomeSpace(toLevel, thisLevel + 1, *it, maxDepth);
			}
		}

    }
    else {
        done = true;
    }
    return done;

}

bool MultiTreeManager::addToLevelStringSumCountsMap(
			MultiTreeManager::level_vertex_info_map& mp,
			int toLevel, int thisLevel, SPSnode * tree,
			MultiTreeManager::vertex_info_set::iterator parent)
{
	bool done = false;

    if (thisLevel <= toLevel) {

        // get how many leaves n there are on the tree
        size_t n = tree->getNumberLeaves();
		
		std::string parentLabel = parent->label;
		std::string parentSegment = "\"" + parentLabel + "\"";
				
		#ifdef MYDEBUG
	
			std::cout << "\nparent is\t" << parentLabel << std::endl;
	
		#endif

        // make a temporary collection container
        SPSnodePtrs container;
        container.reserve(n); // reserve space
		
		std::vector < MultiTreeManager::vertex_info_set::iterator > parentIts;
		parentIts.reserve(n);
		
		// make one copy of the tree for each leaf
        // put copies in a temporary container
        for (size_t i = 0; i < n; ++i) {

            SPSnode* copyTree = new SPSnode(*tree);
            // take the jth copy of tree and find its leaves
            SPSnodePtrs leaves;
            copyTree->getLeaves(leaves);
            // split the ith leaf of this copy of tree
            leaves[i]->nodeExpand();
			
			//find if we want to pursue this one any further
			std::pair < MultiTreeManager::vertex_info_set::iterator, bool > addPair
					= addTreeToLevelStringSumCountsMap(	mp, copyTree);
			
			// put in the graph connection no matter what
			std::string childLabel = addPair.first->label;
				
			// add the parent child connection to the dot graph
			string childSegment = "\"" + childLabel + "\"";
			string line = "\t " + parentSegment + " -> " + childSegment + ";";


			pair<set<string>::iterator, bool> ret;

			ret = (addPair.first->lines).insert(line);
			
			#ifdef MYDEBUG

				std::cout << "\ndot line is\t" << line << std::endl;
			
			#endif
			
			// but only add to the ones to keep going with if this a new place 
			if (addPair.second) {
				
				container.push_back(copyTree);
				
				parentIts.push_back(addPair.first);
				
			}
     
        }
		
		assert (container.size() == parentIts.size());

        // copy the contents of the temporary container to pavings
        pavings.insert(pavings.end(), container.begin(),container.end());

        SPSnodePtrsItr it;

        // for each copy in the temporary container, and its string
		std::vector < MultiTreeManager::vertex_info_set::iterator >::iterator
								lit = parentIts.begin();
        for (it = container.begin(); it < container.end(); ++it, ++lit) {
            // recurse addToOutComeSpace(toLevel, thisLevel+1, copy)
            done = addToLevelStringSumCountsMap(mp ,toLevel, thisLevel + 1,
													*it, *lit);
        }

    }
    else {
		
		// put in invisible constraints
		invisibleSidewaysLinks(mp);
		
        done = true;
    }
    return done;
}


std::pair < MultiTreeManager::vertex_info_set::iterator, bool >
		MultiTreeManager::addTreeToLevelStringSumCountsMap(
					MultiTreeManager::level_vertex_info_map& mp,
					SPSnode * tree)
{
	IntVec thisLevels;
	thisLevels  = tree->getLeafNodeLevels(thisLevels);
	
	size_t leaves = thisLevels.size();
	
	//sort it, in reverse order
	sort( thisLevels.begin(), thisLevels.end(), std::greater_equal<int>());
	
	#ifdef MYDEBUG
	{
		std::cout << "\nLeaves size is\t" << leaves << "\tlevels are ";
		std::ostream_iterator<int> out_it1 (std::cout, " ");
			std::copy ( thisLevels.begin(), thisLevels.end(), out_it1 );
		std::cout << std::endl;
	}
	#endif

	IntVec thisCode;
	int sum = 0;
	if (leaves > 1) { 
		
		int leavesMinusOne = leaves - 1;
	
		// code it, and might as well sum it as the same time
		thisCode = IntVec(leaves-1, 0);
		IntVec::iterator this_it = thisLevels.begin();
		for (int i = 0; i < leavesMinusOne; ++i) {
			int count = 0;
			while ( this_it < thisLevels.end() && (*this_it) == (leavesMinusOne - i) ) {
				count++;
				sum += (*this_it);
				this_it++;
				
			}
			thisCode[i] = count;
		}
	}
	#ifdef MYDEBUG
	{
		std::cout << "sum is " << sum << "\tcodes are ";
		std::ostream_iterator<int> out_it2 (std::cout, " ");
			std::copy ( thisCode.begin(), thisCode.end(), out_it2 );
		std::cout << std::endl;
	}
	#endif
	
	//find if there is a set for this number of leaves
	MultiTreeManager::vertex_info_set this_vertex_info_set;
	
	std::pair < MultiTreeManager::level_vertex_info_map::iterator, bool > mp_ins;

	mp_ins = mp.insert
		( pair < size_t, MultiTreeManager::vertex_info_set > 
			(leaves, this_vertex_info_set) ); 

	MultiTreeManager::level_vertex_info_map::iterator mp_it = mp_ins.first;
	
	//mp_it.second should point to the set for this level 
	
	std::pair < MultiTreeManager::vertex_info_set::iterator, bool > set_ins;
	
	std::string label = getCodedNodeLabel(thisLevels,
											thisCode,
											sum);
		
	set_ins = mp_it->second.insert( 
			MultiTreeManager::VertexInfo(thisLevels, thisCode, sum, label) );
	
	return set_ins; 
}

void MultiTreeManager::invisibleSidewaysLinks(
				MultiTreeManager::level_vertex_info_map& mp)
{
	// each map iterator is (level, vertex info set) 
	for (MultiTreeManager::level_vertex_info_map::const_iterator m_it = mp.begin();
			m_it != mp.end();
			++m_it) {
		// only need to do if more than 1 in set
		if ((m_it->second).size() > 1) {
			// each set iterator is iterator to a vertex info
			MultiTreeManager::vertex_info_set::const_iterator s_it = (m_it->second).begin();
			
			std::string leftLabel = s_it->label;
			std::string leftSegment = "\"" + leftLabel + "\"";
					
			s_it++; // start from here
			
			for (;
					s_it != (m_it->second).end();
					++s_it) {
					
				std::string rightLabel = s_it->label;
				std::string rightSegment = "\"" + rightLabel + "\"";
		
				// add the invisible left-right connection to the dot graph
				string line = "\t " + leftSegment + " -> " + rightSegment 
							+ " [style=invis];";

				s_it->constraint_line = line;
				
				#ifdef MYDEBUG

					std::cout << "\ninvisible dot line is\t" << line << std::endl;
				
				#endif
				
				leftSegment = rightSegment;
			}
		}
	}

}


std::string MultiTreeManager::getCodedNodeLabel(
								const IntVec& sortedLevels,
								const IntVec& codedLevels,
								size_t sum)
{
	std::ostringstream stm;
	std::string delim(" ");
	stm << "<";
	std::ostream_iterator<int> out_it2 (stm, delim.c_str());
	std::copy ( codedLevels.begin(), codedLevels.end(), out_it2 );
	long pos=stm.tellp();
	if (!codedLevels.empty()) stm.seekp (pos - 1);
	//stm << "> -> " << sum << " (" << (sum/(1.0*sortedLevels.size())) << ")";
	stm << ">";
	return stm.str(); 
}



// ----------- end of implementation of MultiTreeManager class ---------------



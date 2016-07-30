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
\brief MultiTreeManager declarations.
*/

#ifndef ___MULTITREEMANAGER_HPP__
#define ___MULTITREEMANAGER_HPP__

#include <string>   // to use the C++ string class
#include <set>      // to use the stl::set container
#include <map>
#include <utility> // for pairs


#include "spsnode.hpp"

namespace subpavings {

	/*! \brief A class which can look into the state space of SPSnode trees.

	Primary method graphOutcomeSpace creates a collection of all the different
	possible SPSnode trees down to a specified number of splits.

	The number of distinct full binary SPStrees after k splits is the Catalan
	number Ck.
	*/

	class MultiTreeManager {
	
	public:
	
		class VertexInfo {
			
			public:
			
				VertexInfo(const IntVec& _sl, const IntVec& _cl, int _s,
							const std::string& _l);
				
				const IntVec sortedLevels;
				
				const IntVec codedLevels;
				
				const int sum;
				
				const std::string label;
				
				mutable std::set < string > lines;

				mutable std::string constraint_line;
				
				bool operator< (const VertexInfo& rhs) const;

			
		};
	
		//typedef std::pair< std::string, int > string_sum;
		
		typedef std::set < VertexInfo > vertex_info_set;
		
		typedef std::map < size_t, vertex_info_set > level_vertex_info_map;
		
		/*! \brief default constructor
		*/
		MultiTreeManager() {}

		/*! \brief Destructor.
		*/
		~MultiTreeManager();

		/*! Map the leaf levels of trees currently in the pavings collection

		Creates a collection of ordered maps, one map for each number of leaves of
		any tree in the pavings collection.  A map for leaves L has keys the unique
		leaf level pattern summaries with L leaves and values the number of trees
		in the pavings collection that has that leaf level pattern summary.

		Console output gives a summary of the frequencies and relative frequencies
		of each leaf level pattern, grouped by number of leaves.
		*/
		void mapPavings();


		/*! Get and map the outcome space resulting from splitting to level toLevel.

		The outcome space is all the possible results of up to and including toLevel
		splits starting from a single root node, ie the number of trees in the
		outcome space for toLevel = k is the Catalan number Ck.

		The map creates a summary of the frequencies of each leaf level pattern
		in the outcome space as a collection of ordered maps, one map for each
		number of leaves of any tree in the outcome space collection.  A map for
		leaves L has keys the unique leaf level pattern summaries with L leaves
		and values the number of trees in the outcome space that has that
		leaf level pattern summary.

		Console output and txt file output both give a summary of the frequencies
		and relative frequencies of each leaf level pattern, grouped by number
		of leaves.

		\param toLevel is the number of splits to go up to.
		\post output file in tab delimited txt format of summary of the
		frequencies and relative frequencies of each leaf level pattern,
		grouped by number of leaves.
		*/
		void makeAndMapOutcomeSpace(int toLevel, int maxDepth = 0);


		/*! Make the outcome space from continual splitting to level \a toLevel.

		The outcome space is all the possible results of up to and including \a toLevel
		splits starting from a single root node. The number of trees in the
		outcome space for \a toLevel = k is the Catalan number Ck but this method 
		will produce as many copies of a tree as there are paths to that tree. 

		Pointers to the trees which result from the splitting are put into the data
		member pavings.

		\param toLevel is the number of splits to go up to.
		\post the data member pavings will contain pointers to all the trees in
		the outcome space.
		*/
		void makeOutcomeSpace(int toLevel, int maxDepth = 0);

		/*! Make the final outcome space from continual splitting to level \a toLevel.

		The outcome space is all the results with exactly \a toLevel
		splits starting from a single root node. The number of trees in the
		outcome space for \a toLevel = k is the Catalan number Ck but this method 
		will produce as many copies of a tree as there are paths to that tree. 

		Pointers to the trees which result from the splitting are put into the data
		member pavings.

		\param toLevel is the number of splits to go up to.
		\post the data member pavings will contain pointers to all the trees in
		the outcome space with \a toLevel splits, with multiple copies of each
		tree where there are multiple paths to it.
		*/
		void makeFinalOutcomeSpace(int toLevel, int maxDepth = 0);

		/*! Make and graph outcome space from continual splitting to level toLevel

		The outcome space is all the possible results of up to and including toLevel
		splits starting from a single root node, ie the number of trees in the
		outcome space for toLevel = k is the Catalan number Ck.

		Pointers to the trees which result from the splitting are put into the data
		member pavings.

		The method makes a dot graph showing all the unique leaf level patterns
		in the outcome space and the relationship between them, ie. which pattern
		can lead to which other patterns.

		\param toLevel is the number of splits to go up to.
		\post the data member pavings will contain pointers to all the trees in
		the outcome space.
		*/
		void makeAndGraphOutcomeSpace(int toLevel, int maxDepth = 0);
		
		level_vertex_info_map makeLevelStringSumCountsMap(int toLevel);
		
		
	private:

		/*! \brief vector of pointers to subpavings managed by this object
		*/
		SPSnodePtrs pavings;

		/*! \brief Copy constructor not implemented.
		*/
		explicit MultiTreeManager
			(const MultiTreeManager& other);

		/*! \brief Copy assignment operator not implemented.
		*/
		MultiTreeManager&
			operator=(const MultiTreeManager& rhs);

		/*! \brief Accumulate splitting outcomes in pavings and add to graph.

		Pointers to the trees which result from the splitting are put into the data
		member pavings.  For each tree made a line is added to the file for a dot
		graph showing the link from parent to child provided the link indicated by
		the line has not already been registered in the graph.

		\param s is the name of the graph to add to.
		\param toLevel is the number of splits to get to in total.
		\param thisLevel is the number of splits on this loop through the method.
		\param tree is the pointer to the tree we are currently working on.
		\param lines is a set of graph lines to make sure we don't have duplicates.
		\return true when done;
		*/
		bool addToOutcomeSpaceAndGraph(const std::string& s, int toLevel, int thisLevel,
											SPSnode * tree, int maxDepth,
											std::set<std::string>& lines);



		/*! \brief Accumulate splitting outcomes in pavings.

		Pointers to the trees which result from the splitting are put into the data
		member pavings.

		\param toLevel is the number of splits to get to in total.
		\param thisLevel is the number of splits on this loop through the method.
		\param tree is the pointer to the tree we are currently working on.
		\return true when done;
		*/
		bool addToOutcomeSpace(int toLevel,
				int thisLevel, SPSnode * tree,
				int maxDepth);

		/*! \brief Accumulate splitting outcomes in pavings but only for 
		final target level \a toLevel.

		Pointers to the trees which are the final result 
		of the splitting are put into the data
		member pavings.

		\param toLevel is the number of splits to get to in total.
		\param thisLevel is the number of splits on this loop through the method.
		\param tree is the pointer to the tree we are currently working on.
		\return true when done;		*/
		bool getFinalLevelOutcomeSpace(int toLevel,
            int thisLevel, SPSnode * tree, int maxDepth);

		bool addToLevelStringSumCountsMap(
						MultiTreeManager::level_vertex_info_map& mp,
						int toLevel, int thisLevel, SPSnode * tree, 
						MultiTreeManager::vertex_info_set::iterator parent);
		
		std::pair < MultiTreeManager::vertex_info_set::iterator, bool >
			addTreeToLevelStringSumCountsMap(
					MultiTreeManager::level_vertex_info_map& mp,
					SPSnode * tree);
		void invisibleSidewaysLinks(
				MultiTreeManager::level_vertex_info_map& mp);
		
		/*! \brief get a string of the leaf levels of given subpaving

		Left to right, 0 is root
		*/
		static std::string getLeafLevelsString(const SPSnode * const spn)
		{
			return spn->getLeafNodeLevelsString();
		}

		static std::string getCodedNodeLabel(
								const IntVec& sortedLevels,
								const IntVec& codedLevels,
								size_t sum);


	};


}

#endif



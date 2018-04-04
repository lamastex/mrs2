/* 
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin
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

/*! \file      FPhyloPOT.hpp
  \brief Declarations for example function class FPhyloPOT 
  (Phylogenetic tree by post order traversal).
*/

#ifndef __FPHYLOPOT_HPP__
#define __FPHYLOPOT_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "interval.hpp"
#include "imath.hpp"
#include "rmath.hpp"
#include "ivector.hpp"
#include <functional>
#include <exception>

#include "intvector.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <gop.hpp>  // cxsc global optimisation

using namespace std;
using namespace cxsc;

#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "PhyloTree.hpp"

//forward class declaration
// A class for a phylogenetic tree
class PhyloTree;    

struct DATA
{
  // set in ReadSequence
  vector<string> seqNames;

  vector< vector<char> > rawSequences;

  int No_seq;

  int Seq_length;

  vector<int> BPatternCounts;

  // storage for the base sequences we find in chars and in coded form
  vector<string> BPatterns;
  vector< vector<int> > baseSequences;

  int BNo_pattern;
};
// struct DATA SeqData;

/*! \brief A class to create a function object for trans-dimensional 
  phylogenetic likelihood by post-order traversal.

  Uses Jukes Cantor model, ie 4 nucleotides.

  This class creates a representation of the unrooted tree(s) or topology 
  appropriate for the number of taxa using the class PhyloNode.

  Note that at present the class is only set up to deal with unrooted trees, 
  ie the number of dimensions (branches) is determined solely by the number of 
  taxa (tree space).

  The likelihood as a real or interval or HessType, under a particular 
  topology, is calculated by post-order traversal of the tree for that topology.

  Tree structures are read in from a file.  Likelihoods are based on DNA 
  sequence data read into the function object constructor.
 */
class FPhyloPOT: public Fobj
{

  private:
    // ---------------------  data members ------------------------

    //! The maximum allowable characters in a sequence name in the sequence file
    // the enum hack, makes maxCharInSeqName a symbolic name for 20
    enum{maxCharInSeqName = 20};

    //! How many topologies (unrooted) there are for 2, 3, 4, and 5 taxa.
    enum{TwoTaxaTop = 1, ThreeTaxaTop = 1, FourTaxaTop = 3, FiveTaxaTop = 15};

    //! How many dimensions there are for 2, 3, 4 and 5 taxa (unrooted trees).
    enum{TwoTaxaDim = 2, ThreeTaxaDim= 3, FourTaxaDim = 5, FiveTaxaDim = 7};

    //! Size of the character space for this model.
    enum{CharacterSpace = 4};

    /* members inherited from Fobj
    // a flag for working on the log(target shape) scale
    bool UsingLogDensity;  	
    // The initial collection of labelled boxes -- prior support
    vector<LabBox> LabDomainList;
    vector<real> LabDomainPriorIntegralList;
    // To specify uniform, exponential, user_defined.
    int PriorType; 					
    */
    
    //! The model space or number of taxa
    int tree_space;

    //! The number of edges in an unrooted tree
    int n_dimensions;

    //! The number of unrooted topology trees given the model space
    int topologies;

    //! The number of nodes in an unrooted tree given the number of taxa
    int tree_nodes;

    //! A struct for the sequence data read in
    struct DATA SeqData;

    //! A container of pointers to PhyloTrees
    mutable PhyloTreePtrs treeRoots;
    
    // the trees are in dynamic memory and so have to be
    // deleted when this object is destroyed

    //! Track number of interval function calls
    mutable int n_interval_calls;

    //! Track number of real function calls
    mutable int n_real_calls;

    // ------------------------- private member functions ----------------------

    //! Convert nucleotide characters to integers, for coding
    static int Char2Code(char ch);

    //! Convert coded nucleotide integers back to characters
    static char Code2Char(int c);

    //! A function to check and read a line from a file
    static void CheckReadLine(ifstream& ifs, string& line);

    //! A function to deal with an error found in a sequence input file
    static void FileInputError(ifstream& ifs, const string& msg);

    //! Reformats sequences read in to make them easy to analyse for patterns
    static void ReformatSequence(string& line);

    //! Reads sequences from a txt file
    // not static since uses data members, not const since alters this
    int ReadSequence(const string& s);

    /*! \brief Finds patterns in sequences read in and fills in datamembers 
      in SeqData struct
    */
    // not static since uses data members, not const since alters this
    int FindPattern(void);

    //! Print the summmary of results of ReadSequence() and FindPattern()
    void PrintSequence() const;

    /*! \brief Makes trees from a tree file.
      
      Reads in lines from the file, tries to construct a tree with each on, 
      puts trees into treeRoots
    \param s the name of the tree file
    \post treeRoots should be contain one pointer to a tree for every line in 
      the tree file 
    */
    // not static since uses data members, not const since alters this
    // should abort if any problems
    void ReadTrees(const string& s);

    //! Prints all the trees to console output
    // print all the trees
    void PrintTopologyTrees();

    //! Destroys the trees (dynamic memory), used by destructor
    void destroyRoots() const;

  public:
    //! Constructor
    FPhyloPOT(int ts, interval Domain, bool LogPi, int Prior);

    //! Destructor
    ~FPhyloPOT();

    //! interval function object
    interval operator()(const LabBox& lb) const;

    //! real function object
    real operator()(const LabPnt& lp) const;

    //! HessType function object
    HessType operator()(const HTvector& x, const int label = 0) const;

    //! Volume of rooted tree boxes inherited
    virtual real LabBoxVolume(const LabBox& LB)
    {
      return Fobj::LabBoxVolume(LB);
    }

    //! Get number of interval function calls
    int get_interval_calls()
    {
      return n_interval_calls;
    }

    //! Get number of real function calls
    int get_real_calls()
    {
      return n_real_calls;
    }

    //! Get the number of trees
    int getNoTrees()
    {
      return static_cast<int>(treeRoots.size());
    }

};                  // end of FPhyloPOT class declarations
#endif

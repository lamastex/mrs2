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

/*! \file      PhyloTree.hpp
\brief Declarations for a class of objects to represent Phylogenetic Trees
*/

#ifndef __PHYLOTREE_HPP__
#define __PHYLOTREE_HPP__

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
#include <numeric>  // to be able to use accumulate algorithm
#include <exception>
#include <algorithm>

#include "intvector.hpp"

#include <stdio.h>
#include <stdlib.h>

#include <gop.hpp>  // cxsc global optimisation, for HessTypes etc

using namespace std;
using namespace cxsc;

/** @name Forward class declaration
 */
//@{
class PhyloNode;    //!< a class for nodes of a phylogenetic tree
class PhyloTree;    //!< a class for phylogenetic trees
//@}

/** @name Typedefs for containers of nodes
 */
//@{
//! a container of pointers to PhyloNodes
typedef vector<PhyloNode*> PhyloPtrs;
//! iterator over a container of pointers to PhyloNodes
typedef PhyloPtrs::iterator PhyloPtrsItr;
//@}

/** @name Typedefs fro containers of trees
 */
//@{
//! a container of pointers to PhyloTrees
typedef vector<PhyloTree*> PhyloTreePtrs;
//! an iterator over a container of PhyloTrees
typedef PhyloTreePtrs::iterator PhyloTreePtrsItr;
//@}

/** @name Typedefs for nucleotide probabilities as real, interval or HessType
 */
//@{
//! a container of intervals for interval probabilities
typedef vector<interval> IntervalProbs;
//! iterator over container of interval probabilities
typedef IntervalProbs::iterator IntervalProbsIt;
//! a container of intervals for real probabilities
typedef vector<real> RealProbs;
//! iterator over container of real probabilities
typedef RealProbs::iterator RealProbsIt;
//! a container of HessTypes for HessType probabilities
typedef vector<HessType> HessProbs;
//! iterator over container of HessType probabilities
typedef HessProbs::iterator HessProbsIt;
//@}

/** @name Typedefs for function pointers for functions giving transition probabilities
 */
//@{
//! for probabilities as reals
typedef real (*RealTranProb_FctPtr)(const real&, const int, const int);
//! for probabilities as intervals
typedef interval (*IntervalTranProb_FctPtr)(const interval&, const int, const int);
//! for probabilities as HessTypes
typedef HessType (*HessTranProb_FctPtr)(const HessType&, const int, const int);
//@}

//! A struct to make a templatised function to delete objects pointed to by pointers
struct DeleteObject
{
  template<typename T>
    void operator()(const T* ptr) const
  {
    delete ptr;
  }
};

/** @name  Declarations for non-member functions used to calculate transition probabilities
 */
//@{
//! Jukes Cantor formula for transition from nucleotype i to j for an unrooted tree
/*!\param t is branchlength as a real
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as a real */
real PijofT_JC69R (const real& t, const int i, const int j );

//! Jukes Cantor formula for transition from nucleotype i to j for an unrooted tree
/*!\param t is branchlength as an interval
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as an interval */
interval PijofT_JC69I (const interval& t, const int i, const int j );

//! Jukes Cantor formula for transition from nucleotype i to j for an unrooted tree
/*!\param t is branchlength as a HessType
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as a HessType */
HessType PijofT_JC69H (const HessType& t, const int i, const int j );

//! CFN formula for transition from nucleotype i to j with character space 2
/*!\param t is branchlength as a real
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as a real */
real PijofT_CFNR(const real& t, const int i, const int j);

//! CFN formula for transition from nucleotype i to j with character space 2
/*!\param t is branchlength as a interval
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as a interval */
interval PijofT_CFNR(const interval& t, const int i, const int j);

//! CFN formula for transition from nucleotype i to j with character space 2
/*!\param t is branchlength as a HessType
\param i represents a nucleoide
\param j represents a nucleotide
\return probability of transition from nucleotide i to j as a HessType */
HessType PijofT_CFNR(const HessType& t, const int i, const int j);

//@}

//! A class for nodes of a phylogenetic tree
class PhyloNode
{

  private:
    int ibranch;    //!< signifies the edge number, ie the correspondence to the dimensions of the problem
    int seq_no;     //!< every leaf node should have a sequence number
    real time;      //!< some nodes may be given a time in the tree file
    int label;      //!< some nodes may be given a label in the input file
    string seq_name;//!< sequence name for leaf notes, derived from the sequence number and sequence file

                    //!< pointer to a node as the parent (NULL if no parent ie root)
    PhyloNode* parent;
                    //!< container of pointers to child nodes
    PhyloPtrs childNodes;

    //!  Copy constructor is private and not implemented
    /*! Prevents copying of the object. */
    explicit PhyloNode(const PhyloNode& other);

    //! Copy assignment operator is private and not implemented
    /*! Prevents copying of the object. */
    PhyloNode& operator=(const PhyloNode& rhs);

  public:
    //! Default constructor
    explicit PhyloNode() : ibranch(-1), seq_no(-1), time(0.0), label(0), parent(NULL)
    // vector childNodes is not initialised
      {}

    //! Destructor
    ~PhyloNode();

    //! Get the parent of this node
    /*! Returns a copy of the parent node pointer. */
    // this will mean that there are multiple pointers to the parent around - not very safe
    PhyloNode* getParent() const
    {
      return parent;
    }

    //! Get the children of this node
    /*! Returns a copy of the container of children. */
    // this will mean that there are multiple pointers to the children around - not very safe
    PhyloPtrs getChildren() const
    {
      return childNodes;
    }

    //! Get the ibranch member
    int getBranch() const
    {
      return ibranch;
    }

    //! Get the time member
    real getTime() const
    {
      return time;
    }

    //! Get the label member
    int getLabel() const
    {
      return label;
    }

    //! Get the seq_no member
    int getSeqNo() const
    {
      return seq_no;
    }

    //! Get the seq_name member
    string getSeqName() const
    {
      return seq_name;
    }

    //! Set the ibranch member
    void setBranch(const int branch)
    {
      ibranch = branch;
    }

    //! Set the time member
    void setTime(const real t)
    {
      time = t;
    }

    //! Set the label member
    void setLabel(const int n)
    {
      label = n;
    }

    //! Set the seq_no member
    void setSeqNo(const int n)
    {
      seq_no = n;
    }

    //! Set the seq_name member
    void setSeqName(const string& s)
    {
      seq_name = s;
    }

    //! Set the parent of this node
    void setParent(PhyloNode * newparent)
    {
      parent = newparent;
    }

    //! Add a child to this node
    void addChild(PhyloNode * const child);

    size_t noChildren() const;

    size_t noDescendents() const;

    //! Print node details
    /*! tab delimited */
    void printNode() const;

};                  // end of PhyloNode class

class PhyloTree
{

  private:

    // data members
    //! A pointer to a PhyloNode root of the tree
    PhyloNode* root;

    //! The size of the character space for the tree
    int CharacterSpace;
    //private functions

    //! Decodes and makes a tree from a line from the tree file
    /*! Expects sequence numbers to be delimited by commas */
    // returns a pointer to the root node of the new tree
    PhyloNode* DecodeTree(const string& line, const vector<string>& seqnames) const;

    //! A function to deal with an error found in a tree input file
    void TreeFileErrorExit(const string& msg, const PhyloNode * node) const;

    //! Extract a time given in a tree file if the identifier for a time has been found
    real findTime(const string& line, const size_t lineLength, size_t& pos, const PhyloNode* root) const;

    //! Extract a label given in a tree file if the identifier for a label has been found
    int findLabel(const string& line, const size_t lineLength, size_t& pos, const PhyloNode* root) const;

    //! Extract a sequence number from a tree file
    /*! If a sequence name rather than number is given, finds the corresponding sequence number using the names from sequences read in from sequence file. */
    int findSeqNo(const string& line, const size_t lineLength, size_t& pos, const PhyloNode* root, const vector<string>& seqnames) const;

    //! Prints a branch of the tree
    static void PrintBranch(const PhyloNode * const node, int level);

    //!  Copy constructor is private and not implemented
    /*! Prevents copying of the object. */
    explicit PhyloTree(const PhyloTree& other);

    //! Copy assignment operator is private and not implemented
    /*! Prevents copying of the object. */
    PhyloTree& operator=(const PhyloTree& rhs);

    //! Calculates the probabilities as reals of each character for a node given its children's state
    /*! Uses post order traversal of the tree to build probabilities from probabilities in child nodes.
    \param node the node for which to calculate probabilities
    \param x the rvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param rNucProb a vector of reals, passed by reference, to be filled in by the function
    \param tpf a function pointer for transition probabilities as reals
    \return a vector of reals, one for each character in the character space, returned by reference*/
    RealProbs& nodePOTreal(RealTranProb_FctPtr tpf, const PhyloNode * const node, const rvector& x, const vector<int>& pattern, RealProbs& rNucProb) const;

    //! Calculates the probabilities as an interval of each character for a node given its children's state
    /*! Uses post order traversal of the tree to build probabilities from probabilities in child nodes.
    \param node the node for which to calculate probabilities
    \param x the rvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param iNucProb a vector of intervals, passed by reference, to be filled in by the function
    \param tpf a function pointer for transition probabilities as intervals
    \return a vector of intervals, one for each character in the character space, returned by reference*/
    IntervalProbs& nodePOTinterval(IntervalTranProb_FctPtr tpf, const PhyloNode * const node, const ivector& x, const vector<int>& pattern, IntervalProbs& iNucProb) const;

    //! Calculates the probabilities as HessType of each character in a character space for a node given its children's state
    /*! Uses post order traversal of the tree to build probabilities from probabilities in child nodes.
    \param node the node for which to calculate probabilities
    \param x the rvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param hNucProb a vector of HessTypes, passed by reference, to be filled in by the function
    \param tpf a function pointer for probabilities as HessTypes
    \param tpf a function pointer for transition probabilities as HessTypes
    \return a vector of HessTypes, one for each character in the character space, returned by reference
    */
    HessProbs& nodePOTHess(HessTranProb_FctPtr tpf, const PhyloNode * const node, const HTvector& x, const vector<int>& pattern, HessProbs& hNucProb) const;

  public:
    //!Default constructor
    explicit PhyloTree() : root(NULL)
      {}

    //! Constructor using a string description, eg ((1 2) 3 4) and a set of sequence names
    explicit PhyloTree(int cs, string& line, const vector<string>& seqnames);

    //! Constructor using a string description, eg ((1 2) 3 4) but no sequence names
    explicit PhyloTree(int cs, string& line);

    //! Destructor
    ~PhyloTree()
    {               // deletes the node which is the root and all its children
      delete root;
    }

    //! print a tree starting with the root
    void PrintTree() const;

    //! Get the root of the tree
    PhyloNode* getRoot() const
    {
      return root;
    }

    //! Get the number of nodes in the tree, including the root
    size_t getNumberNodes() const;

    //! Get the number of branches in a tree
    size_t getNumberBranches() const
    {
      return getNumberNodes() - 1;
    }

    //! Return the probability of a given pattern for this tree
    /*! version dealing with rvectors, ie points
    // and applying JC69 transition probabilities
    \param x the rvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param rNucProb a vector of reals, passed by reference, to be filled in by the function
    \return a vector of reals, one for each character in the character space, returned by reference */
    RealProbs& fillProbRealJC69(const rvector& x, const vector<int>& pattern, RealProbs& rNucProb) const;

    //! Return the probability of a given pattern for this tree
    /*! version dealing with ivectors, ie boxes
    // and applying JC69 transition probabilities
    \param x the ivector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param iNucProb a vector of intervals, passed by reference, to be filled in by the function
    \return a vector of intervals, one for each character in the character space, returned by reference */
    IntervalProbs& fillProbIntervalJC69(const ivector& x, const vector<int>& pattern, IntervalProbs& iNucProb) const;

    //! Return the probability of a given pattern for this tree
    /*! version dealing with HTvectors
    // and applying JC69 transition probabilities
    \param x the HTvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param hNucProb a vector of HessTypes, passed by reference, to be filled in by the function
    \return a vector of HessTypes, one for each character in the character space, returned by reference */
    HessProbs& fillProbHessJC69(const HTvector& x, const vector<int>& pattern, HessProbs& hNucProb) const;

    //! Return the probability of a given pattern for this tree
    /*! version dealing with rvectors, ie points
    // and applying CFN transition probabilities
    \param x the rvector to calculate the probability of
    \param pattern the nucleotide pattern under which to calculate the probabilities (gives the probabilities at the leaf nodes)
    \param rNucProb a vector of reals, passed by reference, to be filled in by the function
    \return a vector of reals, one for each character in the character space, returned by reference */
    RealProbs& fillProbRealCFN(const rvector& x, const vector<int>& pattern, RealProbs& rNucProb) const;

};                  // end of PhyloTree class declarations

/*! \brief operator const LabPnt& X
Utility function to sum reals */
real realSum(real sumSoFar, const real r);

/*! \brief operator const LabBox& X
Utility function to sum intervals */
interval intervalSum(interval sumSoFar, const interval x);

/*! \brief operator const HTvector& X
Utility function to sum HessTypes */
HessType hessSum(HessType sumSoFar, const HessType x);
#endif

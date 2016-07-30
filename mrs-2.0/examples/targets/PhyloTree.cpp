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

/*! \file      PhyloTree.cpp
\brief Implementation for for a class of objects to represent Phylogenetic Trees
*/

#include "PhyloTree.hpp"

// functions which conform the the typedefs for transition probability functions

// conforms to typedef for *RealTranProb_FctPtr
// Jukes Cantor / CFN formula for transition from nucleotype 
// i to j for an unrooted tree
// t is branchlength
// returns a real
real PijofT_R (const real& t, const int i, const int j, const int k )
{
  if(k==4)
  {
    return (i==j) ? ( (0.25)+(0.75 * exp(-(4.0/3.0)*t)) ) : 
        	( (0.25)-(0.25 * exp(-(4.0/3.0)*t)) );
  }
  else if(k==2)
  {
    return (i==j) ? ( (0.5)+(0.5 * exp(-2.0*t)) ) : 
      		( (1.0 - exp(-2.0*t))/2.0 );
  }
  else 
  {
    cout << "unknown character space: neither JC69 nor CFN are applicable" << endl;
    exit(0);
  }
}

// conforms to typedef for *IntervalTranProb_FctPtr
// Jukes Cantor  / CFN formula for transition from nucleotype 
// i to j for an unrooted tree
// t is branchlength as an interval
// returns an interval
interval PijofT_I (const interval& t, const int i, const int j, const int k )
{
  if(k==4)
  {
    return (i==j) ? ( (0.25)+(0.75 * exp(-(4.0/3.0)*t)) ) : 
    		( (0.25)-(0.25 * exp(-(4.0/3.0)*t)) );
  }
  else if(k==2)
  {
    return (i==j) ? ( (0.5)+(0.5 * exp(-2.0*t)) ) : 
      		( (1.0 - exp(-2.0*t))/2.0 );
  }
  else
  {
    cout << "unknown character space: neither JC69 nor CFN are applicable" << endl;
    exit(0);
  }
}

// conforms to typedef for *GradTranProb_FctPtr
// Jukes Cantor / CFN formula for transition from nucleotype 
// i to j for an unrooted tree
// t is branchlength
// returns a GradType
GradType PijofT_G (const GradType& t, const int i, const int j, const int k )
{
  if(k==4)
  {
    return (i==j) ? ( (0.25)+(0.75 * exp(-(4.0/3.0)*t)) ) : 
    		( (0.25)-(0.25 * exp(-(4.0/3.0)*t)) );
  }
  else if(k==2)
  {
    return (i==j) ? ( (0.5)+(0.5 * exp(-2.0*t)) ) : 
      		( (1.0 - exp(-2.0*t))/2.0 );
  }
  else
  {
    cout << "unknown character space: neither JC69 nor CFN are applicable" << endl;
    exit(0);
  }
}

// conforms to typedef for *HessTranProb_FctPtr
// Jukes Cantor / CFN formula for transition from nucleotype 
// i to j for an unrooted tree
// t is branchlength
// returns a HessType
HessType PijofT_H (const HessType& t, const int i, const int j, const int k )
{
  //HessType p;
  //p = (i==j) ? ( (0.25)+(0.75 * exp(-(4.0/3.0)*t)) ) : 
  //  ( (0.25)-(0.25 * exp(-(4.0/3.0)*t)) );
  //return p;
  if(k==4)
  {
    return (i==j) ? ( (0.25)+(0.75 * exp(-(4.0/3.0)*t)) ) : 
    		( (0.25)-(0.25 * exp(-(4.0/3.0)*t)) );
  }
  else if(k==2)
  {
    return (i==j) ? ( (0.5)+(0.5 * exp(-2.0*t)) ) : 
      		( (1.0 - exp(-2.0*t))/2.0 );
  }
  else
  {
    cout << "unknown character space: neither JC69 nor CFN are applicable" << endl;
    exit(0);
  }
}


//------------------------------ PhyloNode class implementation

// Destructor
PhyloNode::~PhyloNode()
{
  // call DeleteObject to delete all the child nodes
  for_each(childNodes.begin(), childNodes.end(), DeleteObject());

}                   // end destructor

// Add a child to this node
void PhyloNode::addChild(PhyloNode * const child)
{
  child->setParent(this);
  childNodes.push_back(child);
}

size_t PhyloNode::noChildren() const
{
  return childNodes.size();
}

size_t PhyloNode::noDescendents() const
{
  size_t noDes = childNodes.size();
  //iterate over each child and repeat
  PhyloPtrs::const_iterator it;
  //PhyloPtrsItr it;
  for (it=childNodes.begin(); it < childNodes.end(); it++)
  {
    noDes += (*it)->noDescendents();
  }

  return noDes;
}

void PhyloNode::printNode() const
{
  // branch, sequence number, sequence name, time, label, newline
  std::cout << "ibranch is \t" << ibranch << "\t" << "seq_no \t" << seq_no
    << "\t" << "time \t" << time << "\t" << "label \t" << label  << "\t" 
    << "name \t" << seq_name << endl;
}

// --------------------------- end of PhyloNode implementation -----------------

// --------------------------- start PhyloTree implementation ------------------

// Constructor
PhyloTree::PhyloTree(int cs, string& line, const vector<string>& seqnames) 
: 
CharacterSpace(cs)
{
  // use decode tree to decode the line, returns a pointer to a PhyloNode
  root = DecodeTree(line, seqnames);

}                   // end constructor

// Constructor
PhyloTree::PhyloTree(int cs, string& line) : CharacterSpace(cs)
{

  // make an empty container
  vector<string> seqnames;
  
  // use decode tree to decode the line, returns a pointer to a PhyloNode
  root = DecodeTree(line, seqnames);

}                   // end constructor

// Get the number of nodes in the tree, including the root
size_t PhyloTree::getNumberNodes() const
{
  size_t noNodes = 0;

  if (root!=NULL)
  {

    noNodes++;      // count the root
    // and add on descendents of root
    noNodes += root->noDescendents();
  }

  return noNodes;
}

// print a tree starting with the root
void PhyloTree::PrintTree() const
{

  PrintBranch(root, 0);
  std:: cout << std::endl;
  std:: cout << std::endl;
}

// return the probability of a given pattern for this tree
// version dealing with rvectors, ie points
// and applying JC69/CFN transition probabilities
RealProbs& PhyloTree::fillProbReal(const rvector& x, 
                                       const vector<int>& pattern, 
                                       RealProbs& rNucProb) const
{

  // check we have some children
  if (root->noChildren() == 0)
  {
    std::cerr << "Error trying to calculate probabilities "
              << "over tree with no children" << std::endl;
    exit(1);
  }

  // check the dimensions of the point matches dimensions for this tree
  int d = Ub(x) - Lb(x) + 1;
  if(getNumberBranches() != static_cast<size_t>(d))
  {
    std::cerr << "tree dimensions (branches) " << getNumberBranches()
      << " does not match dimensions of point " << d << std::endl;
    exit(1);
  }

  // use the root as the node to start post order traversal
  // and pass in pointer to JC69 transition probabilites function
  rNucProb = nodePOTreal(&PijofT_R, root, x, pattern, rNucProb);

  // return the same vector but note that return is by reference
  return rNucProb;

}

// return the probability of a given pattern for this tree
// version dealing with ivectors, ie boxes
// and applying JC69/CFN transition probabilities
IntervalProbs& PhyloTree::fillProbInterval(const ivector& x, 
                                               const vector<int>& pattern, 
                                               IntervalProbs& iNucProb) const
{

  // check we have some children
  if (root->noChildren() == 0)
  {
    std::cerr << "Error trying to calculate probabilities "
              << "over tree with no children" << std::endl;
    exit(1);
  }

  // check the dimensions of the box matches dimensions for this tree
  int d = Ub(x) - Lb(x) + 1;
  if(getNumberBranches() != static_cast<size_t>(d))
  {
    std::cerr << "tree dimensions (branches) " << getNumberBranches()
      << " does not match dimensions of point " << d << std::endl;
    exit(1);
  }

  // use the root as the node to start post order traversal
  // and pass in pointer to JC69 transition probabilites function
  iNucProb = nodePOTinterval(&PijofT_I, root, x, pattern, iNucProb);

  // return the same vector but note that return is by reference
  return iNucProb;

}

// return the probability of a given pattern for this tree
// version dealing with GTvectors
// and applying JC69/CFN transition probabilities
GradProbs& PhyloTree::fillProbGrad(const GTvector& x, 
                                       const vector<int>& pattern, 
                                       GradProbs& gNucProb) const
{

  // check we have some children
  if (root->noChildren() == 0)
  {
    std::cerr << "Error trying to calculate probabilities "
              << "over tree with no children" << std::endl;
    exit(1);
  }

  // check the dimensions of the HTVector matches dimensions for this tree
  int d = x.Dim();
  if(getNumberBranches() != static_cast<size_t>(d))
  {
    std::cerr << "tree dimensions (branches) " << getNumberBranches()
      << " does not match dimensions of point " << d << std::endl;
    exit(1);
  }

  // use the root as the node to start post order traversal
  // and pass in pointer to JC69/CFN transition probabilites function
  gNucProb = nodePOTGrad(&PijofT_G, root, x, pattern, gNucProb);

  // return the same vector but note that return is by reference
  return gNucProb;

}

// return the probability of a given pattern for this tree
// version dealing with HTvectors
// and applying JC69/CFN transition probabilities
HessProbs& PhyloTree::fillProbHess(const HTvector& x, 
                                       const vector<int>& pattern, 
                                       HessProbs& hNucProb) const
{

  // check we have some children
  if (root->noChildren() == 0)
  {
    std::cerr << "Error trying to calculate probabilities "
              << "over tree with no children" << std::endl;
    exit(1);
  }

  // check the dimensions of the HTVector matches dimensions for this tree
  int d = x.Dim();
  if(getNumberBranches() != static_cast<size_t>(d))
  {
    std::cerr << "tree dimensions (branches) " << getNumberBranches()
      << " does not match dimensions of point " << d << std::endl;
    exit(1);
  }

  // use the root as the node to start post order traversal
  // and pass in pointer to JC69/CFN transition probabilites function
  hNucProb = nodePOTHess(&PijofT_H, root, x, pattern, hNucProb);

  // return the same vector but note that return is by reference
  return hNucProb;

}

// decodes a tree from a line
// expects all sequence numbers to be followed by ',' or )' 
// so that it can distinguish '12' from '1,2'
// returns a pointer to the root node of the new tree
PhyloNode* PhyloTree::DecodeTree(const string& line, 
                                 const vector<string>& seqnames) const
{

  // set up a pointer to node for a new tree
  PhyloNode* newtree = NULL;
  PhyloNode* currentParent = NULL;
  PhyloNode* currentNode = NULL;

  int nodesRead = 0;// to track number of nodes read in
  int level = 0;    // to track of tree levels
  int branch = -1;  // to track branch numbers

  string::iterator it;
  size_t pos = 0;   // a counter for characters in the line dealt with

                    // the length of the line
  size_t lineLength = line.size();

  while (pos < lineLength)
  {

    real time = 0.0;
    int label = 0;
    int sequenceNumber = 0;

    char ch = line[pos];
    // pos should show the position we are at in the string

    switch (ch)
    {
      case ' ':     // ignore spaces
        pos++;
        break;

      case ',':     // ignore commas
        // commas are used as delimiters between sequence numbers, 
        // but once found
        // and used to identify the end of the number, they can be ignored
        pos++;
        break;

      case ';':     // ignore semicolons
        pos++;
        break;

      case '(':     // go down a level in the tree, making new nodes 
                    // for each level
        // set up a root node for a new tree

        try
        {
                    // ibranch is -1, parent is NULL
          currentNode = new PhyloNode();
        }
        catch (bad_alloc&)
        {
          cout << "Error allocating memory in DecodeTree" << std::endl;
          delete newtree;
          exit(1);
        }

        if (newtree == NULL)
        {
          newtree = currentNode;
          currentParent = newtree;
          // currentParent is NULL
        }
        else        // there is already a root node
        {
          // make a new node
                    // set branch of currentNode
          currentNode->setBranch(branch);
                    // make the current node a child of the current parent
          currentParent->addChild(currentNode);
          currentParent = currentNode;
        }

        level++;
        branch++;
        nodesRead++;
        pos++;

        break;      // end '('

      case ')':     // go up a level in the tree

        level--;
                    // check we can go up
        if (level < 0)
        {
          string msg = "Error in tree file: reading too many ')'";
          // delete any tree that has been set up
          TreeFileErrorExit(msg, newtree);
        }
        else        // level is okay
        {
          // make the currentParent into the parent 
          // of the present current parent
          currentNode = currentParent;
          currentParent = currentParent->getParent();
          pos++;
        }

        // no change in nodesRead
        // no change in branch

        break;      // end ')'

      case ':':     // signifies start of a time given as a length, 
                    // eg ((1:0.5,2;0.5):0.25,3:0.75, 4:0.75);
        time = findTime(line, lineLength, pos, newtree);
        // this moves pos to be the first character following 
        // the time characters

        // make this the time for the current node
        currentNode->setTime(time);

        // no change in level
        // no change in branch

        break;      // end of ':'

      case '#':     // signifies the start of a label, which is expected to a 
                    // single digit integer
        label = findLabel(line, lineLength, pos, newtree);
        // findLabel() updates pos to be the character after 
        // the end of the label

                    // make the currentNode label the integer at this pos
        currentNode->setLabel(label);

        // no change in level
        // no change in branch

        break;      // end of '#'

      case '$' :    // signifies the start of a label, which is expected 
                    // to a single digit integer
        label = findLabel(line, lineLength, pos, newtree);
        // findLabel() updates pos to be the character after 
        // the end of the label

                    // make the currentNode label the integer at this pos
        currentNode->setLabel(label);

        // no change in level
        // no change in branch

        break;      // end of '$'

      default:      // not one of the recognised characters signifying 
                    // tree structure or a character to be ignored
        // so expected to be either a sequence number or a sequence name
        // we need to identify the sequence number and make a new node and add 
        // it to the tree

        sequenceNumber = findSeqNo(line, lineLength, pos, newtree, seqnames);
        // pos will be moved to first character following the sequence 
        // identifier (number or name)

        // make a new node
        try
        {
                    // ibranch is -1, parent is NULL
          currentNode = new PhyloNode();
        }
        catch (bad_alloc&)
        {
          cout << "Error allocating memory in DecodeTree" << std::endl;
          delete newtree;
          exit(1);
        }

                    // set branch of currentNode
        currentNode->setBranch(branch);
                    // set the sequence number of the current node
        currentNode->setSeqNo(sequenceNumber);

        if (!seqnames.empty())
        {
                    // and set the corresponding name
          currentNode->setSeqName(seqnames[sequenceNumber-1]);
        }
        // make it a child of the current parent
        currentParent->addChild(currentNode);

        // no change in level;
        branch++;
        nodesRead++;

        break;      // end default
    }               // end switch
  }                 // end while

  // we expect to find that level = 0 at the end
  // and something has been read in for the tree
  // and number of branches+1 is number of nodes
  if (level != 0 || newtree == NULL || (branch+1 != nodesRead))
  {
    // error in tree file

    string msg = 
    "Error in tree file, please check the number of branches and nodes:  Expecting unrooted tree";
    TreeFileErrorExit(msg, newtree);
  }
  return newtree;
}

void PhyloTree::TreeFileErrorExit(const string& msg, 
                                  const PhyloNode * node) const
{
  std::cerr << msg << std::endl << std::endl;
  delete node;
  exit(1);

}

// finds a double time number
// called when a ':' character has been found in a tree file
// changes the value of pos so that after routine has been called, pos points 
// to character after time
// aborts if expected time not found

real PhyloTree::findTime(const string& line, const size_t lineLength, 
                         size_t& pos, const PhyloNode* node) const
{
  pos++;
                    // characters accepted as part of a double number
  string doubleChars = ".0123456789";
                    // space characters
  string space = " ";

  // find the first non space and make sure it is acceptable 
  // as the start of a double
  pos = line.find_first_not_of(space, pos);

  if (pos==string::npos)
  {
    string msg = 
    "Error in tree file: no double for time following ':' character ";
    TreeFileErrorExit(msg, node);
  }

  size_t firstPos = line.find_first_of(doubleChars, pos);

  // check that first integer is the first non space character following ':'
  if ((firstPos==string::npos) || (firstPos != pos))
  {
    string msg = "Error in tree file: no integer following '#' character";
    TreeFileErrorExit(msg, node);
  }
  // pos is the same as firstPos

  // find last of the following doubleChars
  size_t posUpTo = line.find_first_not_of(doubleChars, pos);

                    // there could be no non-numeric character before the end
  if (posUpTo==string::npos)
  {
    posUpTo = lineLength;
  }

  // using the part of the string comprising the numbers
                    // convert to a stream
  istringstream sin(line.substr(pos, posUpTo-pos));
  double dblValue = 0.0;
  sin >> dblValue;  // extract number from the stream

  real retValue = _real(dblValue);

  // make this the time for the current node
  pos = posUpTo;    // move pos up the first non-numeric character

  return retValue;
}

// finds an integer label
// called when a '$' or '#' character has been found in a tree file
// changes the value of pos so that after routine has been called, 
// pos points to character after label
// aborts if expected label not found or more than 1 digit long
int PhyloTree::findLabel(const string& line, const size_t lineLength, 
                         size_t& pos, const PhyloNode* node) const
{
  pos++;            // move to the next character after the '#' or '$' character

                    // characters accepted as part of an integer
  string intChars = "0123456789";
                    // space characters
  string space = " ";

  // find the first non space and make sure it is an integer
  pos = line.find_first_not_of(space, pos);

  if (pos==string::npos)
  {
    string msg = 
    "Error in tree file: no integer label following '#' or '$' character";
    TreeFileErrorExit(msg, node);
  }

  size_t numerpos = line.find_first_of(intChars, pos);

  // check that first integer is the first non space character following ':'
  if ((numerpos==string::npos) || (numerpos != pos))
  {
    string msg = 
    "Error in tree file: no integer following '#' or '$' character ";
    TreeFileErrorExit(msg, node);
  }

  int retValue = 1000;

  // string to int
  string singleCharStr = line.substr(pos, 1);
  istringstream sin(singleCharStr);
  sin >> retValue;  // store integer at pos as sequenceNumber

  pos++;

  // now we expect a non integer character or the end of the line
  if (pos < lineLength)
  {
    numerpos = line.find_first_of(intChars, pos);
    if (numerpos == pos)
    {
      string msg = 
      "Sorry:  DecodeTree() can only cope with single-digit labels at present";
      TreeFileErrorExit(msg, node);
    }
  }                 // end check on next digit

  return retValue;
}

// expects sequence numbers to be delimited by commas
int PhyloTree::findSeqNo(const string& line, const size_t lineLength, 
                         size_t& pos, const PhyloNode* node, 
                         const vector<string>& seqnames) const
{

                    // characters accepted as part of an integer
  string intChars = "0123456789";
                    //characters accepted as part of an alphabetical label
  string alphaChars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

  // we know that there will be no spaces at the start since they are dealt 
  // with in DecodeTree()
  // we expect each distinct sequence number to be followed by a comma or a ), 
  // signifying the end of the number

  int sequenceNumber = -1;

  int noseq = static_cast<int>(seqnames.size());

  // check if character at this pos numeric?
  size_t numerpos = line.find_first_of(intChars, pos);

                    // the character at pos is numeric
  if (pos == numerpos)
  {
    // should be the start of our sequence number

    // find the non numeric character
                    // first occurence of something non-int
    size_t notnumerpos = line.find_first_not_of(intChars, pos);

    // the number we want starts at pos and has it's final character 
    // at notnumerpos-1
    // and so has length notnumberpos-1 - pos +1 = notnumererpos-pos
    string numberStr = line.substr(pos, (notnumerpos-pos));
    // string to int
    istringstream sin(numberStr);
                    // store the integer we are after as sequenceNumber
    sin >> sequenceNumber;

    //check sequenceNumber against number of sequence names if given
    if ((noseq > 0) && (noseq < sequenceNumber))
    {
      string msg = 
      "Error in tree file: A sequence number is \n greater than the number of sequences";
      TreeFileErrorExit(msg, node);
    }

                    // move pos to next character after the number
    pos = notnumerpos;
  }

  else              // the character at pos is not numeric
  {
    // should be the start of a recognised sequence name

    if (noseq == 0)
    {
      string msg = 
      "Error constructing tree: can deal only with sequences identified by number";
      TreeFileErrorExit(msg, node);
    }

    size_t posUpTo = line.find_first_not_of(alphaChars, pos);

    // but really there should be at least one non-alpha character 
    // before the end?
    if (posUpTo==string::npos)
    {
      posUpTo = lineLength;
    }

    // using the part of the string comprising the alpha characters
    string name = line.substr(pos, posUpTo-pos);

    // compare to the seqnames

    int i = 0;
    while ((sequenceNumber == -1) && (i < noseq))
    {

                    // find the ith sequence name
      string seqName = seqnames[i];

      // make both names all lower case
      std::transform(name.begin(), name.end(), name.begin(), ::tolower);
      std::transform(seqName.begin(), 
                     seqName.end(), seqName.begin(), ::tolower);

      if (seqName == name)
      {
        sequenceNumber = i+1;

      }
      i++;
    }               // end while

                    // did not find a match for a sequence name
    if (sequenceNumber == -1)
    {
      string msg = 
      "Error in tree file: a name in tree file did not match \n any in the sequence file (case insensitive matching)";
      TreeFileErrorExit(msg, node);
    }

    pos = posUpTo;  // move pos up the first non-alpha character
  }                 // end of alpha name for node identifier

  return sequenceNumber;
}

// print a level of the tree starting with the root
void PhyloTree::PrintBranch(const PhyloNode * const node, int level)
{
  PhyloPtrs theChildren = node->getChildren();
  PhyloPtrsItr it;

  for (int i = 0; i < level; i++)
  {
    std:: cout << "\t";
  }

  node->printNode();

  if(!theChildren.empty())
  {
    for (int i = 0; i < level; i++)
    {
      std::cout << "\t";
    }
    std::cout << "Children are:" << std::endl;
  }

  level++;
  for (it = theChildren.begin(); it < theChildren.end(); it++)
  {
    PrintBranch(*it, level);
  }
}

// return the probability of a given pattern on a tree with given root
// version dealing with rvectors, ie points
RealProbs& PhyloTree::nodePOTreal(RealTranProb_FctPtr tpf, 
                                  const PhyloNode * const node, 
                                  const rvector& x, 
                                  const vector<int>& pattern, 
                                  RealProbs& rNucProb) const
{
  rNucProb.clear(); // make sure the container is empty
  rNucProb.reserve(CharacterSpace);

                    // get the children of the node of this tree
  PhyloPtrs theChildren = node->getChildren();

  PhyloPtrsItr it;

  // CharacterSpace is the number of possible characters for nucleotides
  for (int thisNuc=0; thisNuc < CharacterSpace; thisNuc++)
  {
    real prob(1.0);

    // for each child
    for (it = theChildren.begin(); it < theChildren.end(); it++)
    {

      PhyloNode* child = *it;
                    // add one because rvectors are indexed 1 to n not 0 to n-1
      int branchIndex = child->getBranch() + 1;
                    // the element in x corresponding to the child's branch
      real ptElem = x[branchIndex];
      int sequenceNumber = child->getSeqNo();

                    // this child of this node has no children
      if (child->noChildren() < 1)
      {
        // get the coded character for this child's sequence 
        // number in the supplied pattern
        int childNuc = pattern[sequenceNumber-1];

        // get prob of transition from childNuc to thisNuc
        prob = prob * tpf(ptElem, thisNuc, childNuc, CharacterSpace);

      }             // end if child has no children
      else          // this child has children
      {
                    // fill this in by recursing on fillProbReal with child
        RealProbs childProbs;
        childProbs = nodePOTreal(tpf, child, x, pattern, childProbs);
        real tt(0.0);
        // run over the childProbs, accumulating
        for (int j = 0; j < CharacterSpace; j++)
        {
          tt = tt+ (tpf(ptElem, thisNuc, j, CharacterSpace) * childProbs[j]);
        }
        prob = prob * tt;
      }             // end case child has children
      // the child, whether it has children or not, has now 
      // multiplied prob by something
    }               // end of loop through children
    // we now have a prob for thisNuc to which each child has contributed
                    // copy of prob goes into rNucProb
    rNucProb.push_back(prob);
  }                 // end of loop through CharacterSpace
  // rNucProb should have had CharacterSpace intervals pushed into it
  // return the same vector but note that return is by reference
  return rNucProb;
}

// return the probability of a given pattern for a given node
// version dealing with ivectors, ie boxes
IntervalProbs& PhyloTree::nodePOTinterval(IntervalTranProb_FctPtr tpf, 
                                          const PhyloNode * const node, 
                                          const ivector& x, 
                                          const vector<int>& pattern, 
                                          IntervalProbs& iNucProb) const
{
  iNucProb.clear(); // make sure the container is empty
  iNucProb.reserve(CharacterSpace);

                    // get the children of the node
  PhyloPtrs theChildren = node->getChildren();

  PhyloPtrsItr it;

  // CharacterSpace is the number of possible characters for nucleotides
  for (int thisNuc=0; thisNuc < CharacterSpace; thisNuc++)
  {
    interval prob(1.0, 1.0);

    // for each child
    for (it = theChildren.begin(); it < theChildren.end(); it++)
    {

      PhyloNode* child = *it;
                    // add one because ivectors are indexed 1 to n not 0 to n-1
      int branchIndex = child->getBranch() + 1;
                    // the interval in x corresponding to the child's branch
      interval boxSide = x[branchIndex];
      int sequenceNumber = child->getSeqNo();

                    // this child of this node has no children
      if (child->noChildren() < 1)
      {

        // get the coded character for this child's sequence number in the 
        // pattern_i-th base sequence
        int childNuc = pattern[sequenceNumber-1];

        // get prob of transition from childNuc to thisNuc
        prob = prob * tpf(boxSide, thisNuc, childNuc, CharacterSpace);

      }             // end if child has no children

      else          // this child has children
      {

                    // fill this in by recursing on fillProbInterval with child
        IntervalProbs childProbs;
        childProbs = nodePOTinterval(tpf, child, x, pattern, childProbs);

        interval tt(0.0, 0.0);

        // run over the childProbs, accumulating
        for (int j = 0; j < CharacterSpace; j++)
        {
          tt = tt+ (tpf(boxSide, thisNuc, j, CharacterSpace) * childProbs[j]);
        }
        prob = prob * tt;
      }             // end case child has children

      // the child, whether it has children or not, has now 
      // multiplied prob by something
    }               // end of loop through children

    // we now have a prob for thisNuc to which each child has contributed
                    // copy of prob goes into iNucProb
    iNucProb.push_back(prob);
  }                 // end of loop through CharacterSpace

  // iNucProb should have had CharacterSpace intervals pushed into it

  // return the same vector but note that return is by reference
  return iNucProb;
}

// return the probability of a given pattern for given node
// version dealing with GradTypes, for centered-form range enclosures 
// but in the style of c-xsc's global optimisation eith HessTypes
GradProbs& PhyloTree::nodePOTGrad(GradTranProb_FctPtr tpf, 
                                  const PhyloNode * const node, 
                                  const GTvector& x, 
                                  const vector<int>& pattern, 
                                  GradProbs& gNucProb) const
{
  gNucProb.clear(); // make sure the container is empty
  gNucProb.reserve(CharacterSpace);

                    // get the children of the node
  PhyloPtrs theChildren = node->getChildren();

  PhyloPtrsItr it;

  int d = x.Dim();  // dimensions of x

  // CharacterSpace is the number of possible characters for nucleotides
  for (int thisNuc=0; thisNuc < CharacterSpace; thisNuc++)
  {
    GradType prob(d);
    prob = 1.0;

    // for each child
    for (it = theChildren.begin(); it < theChildren.end(); it++)
    {

      PhyloNode* child = *it;
                    // add one because ivectors are indexed 1 to n not 0 to n-1
      int branchIndex = child->getBranch() + 1;
                    // the element of x corresponding to the child's branch
      GradType ptElem(x[branchIndex]);
      int sequenceNumber = child->getSeqNo();

                    // this child of this node has no children
      if (child->noChildren() < 1)
      {

        // get the coded character for this child's sequence number in 
        // the pattern_i-th base sequence
        int childNuc = pattern[sequenceNumber-1];

        // get prob of transition from childNuc to thisNuc
        prob = prob * tpf(ptElem, thisNuc, childNuc, CharacterSpace);

      }             // end if child has no children

      else          // this child has children
      {

                    // fill this in by recursing on fillProbHess with child
        GradProbs childProbs;
        childProbs = nodePOTGrad(tpf, child, x, pattern, childProbs);

        GradType tt(d);
        tt = 0.0;

        // run over the childProbs, accumulating
        for (int j = 0; j < CharacterSpace; j++)
        {
          tt = tt + (tpf(ptElem, thisNuc, j, CharacterSpace) * childProbs[j]);
        }

        prob = prob * tt;
      }             // end case child has children

      // the child, whether it has children or not, 
      // has now multiplied prob by something
    }               // end of loop through children

    // we now have a prob for thisNuc to which each child has contributed
                    // copy of prob goes into hNucProb
    gNucProb.push_back(prob);
  }                 // end of loop through CharacterSpace

  // hNucProb should have had CharacterSpace intervals pushed into it

  return gNucProb;  // return the same vector but note that return is by reference
}

// return the probability of a given pattern for given node
// version dealing with HessTypes, for global optimisation
HessProbs& PhyloTree::nodePOTHess(HessTranProb_FctPtr tpf, 
                                  const PhyloNode * const node, 
                                  const HTvector& x, 
                                  const vector<int>& pattern, 
                                  HessProbs& hNucProb) const
{
  hNucProb.clear(); // make sure the container is empty
  hNucProb.reserve(CharacterSpace);

                    // get the children of the node
  PhyloPtrs theChildren = node->getChildren();

  PhyloPtrsItr it;

  int d = x.Dim();  // dimensions of x

  // CharacterSpace is the number of possible characters for nucleotides
  for (int thisNuc=0; thisNuc < CharacterSpace; thisNuc++)
  {
                    // HessType of same dimensions as x
    HessType prob(d);
    prob = 1.0;

    // for each child
    for (it = theChildren.begin(); it < theChildren.end(); it++)
    {

      PhyloNode* child = *it;
                    // add one because ivectors are indexed 1 to n not 0 to n-1
      int branchIndex = child->getBranch() + 1;
                    // the element of x corresponding to the child's branch
      HessType ptElem(x[branchIndex]);
      int sequenceNumber = child->getSeqNo();

                    // this child of this node has no children
      if (child->noChildren() < 1)
      {

        // get the coded character for this child's sequence number in 
        // the pattern_i-th base sequence
        int childNuc = pattern[sequenceNumber-1];

        // get prob of transition from childNuc to thisNuc
        prob = prob * tpf(ptElem, thisNuc, childNuc, CharacterSpace);

      }             // end if child has no children

      else          // this child has children
      {

                    // fill this in by recursing on fillProbHess with child
        HessProbs childProbs;
        childProbs = nodePOTHess(tpf, child, x, pattern, childProbs);

        HessType tt(d);
        tt = 0.0;

        // run over the childProbs, accumulating
        for (int j = 0; j < CharacterSpace; j++)
        {
          tt = tt+ (tpf(ptElem, thisNuc, j, CharacterSpace) * childProbs[j]);
        }

        prob = prob * tt;
      }             // end case child has children

      // the child, whether it has children or not, 
      // has now multiplied prob by something
    }               // end of loop through children

    // we now have a prob for thisNuc to which each child has contributed
                    // copy of prob goes into hNucProb
    hNucProb.push_back(prob);
  }                 // end of loop through CharacterSpace

  // hNucProb should have had CharacterSpace intervals pushed into it

  return hNucProb;  // return the same vector but note that return is by reference
}

// -------------------------------------- end PhyloTree implementation ---------

// utility function to sum reals
real realSum(real sumSoFar, const real r)
{
  return sumSoFar+r;
}

// utility function to sum intervals
interval intervalSum(interval sumSoFar, const interval x)
{
  return sumSoFar+x;
}

// utility function to sum GradTypes
GradType gradSum(GradType sumSoFar, const GradType x)
{
  return sumSoFar+x;
}

// utility function to sum HessTypes
HessType hessSum(HessType sumSoFar, const HessType x)
{
  return sumSoFar+x;
}

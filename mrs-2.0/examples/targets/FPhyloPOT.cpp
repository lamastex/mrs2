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

/*! \file      FPhyloPOT.cpp
  \brief Implementation for example function class FPhyloPOT 
  (Phylogenetic tree by post order traversal).
*/

#include "FPhyloPOT.hpp"

// Constructor
FPhyloPOT::FPhyloPOT(int ts, int cs, interval Domain, bool CentFrm, bool LogPi, int Prior):
tree_space(ts), CharacterSpace(cs), CenteredForm(CentFrm), n_interval_calls (0), 
n_real_calls(0), n_gradtype_calls(0)
{
  // an Fobj function, sets the Fobj member UsingLogDensity
  setUsingLogDensity (LogPi);
  PriorType = Prior;// PriorType is an inherited member from Fobj

  // file names for the tree and the sequences files
  string treeFileName;
  string seqFileName;

  //hard-coded the type of range enclosure here
  //CenteredForm=true;
  //CenteredForm=false;

  switch (tree_space)
  {

    // use just one treefile with multiple lines, one line for each topology

    /* need to refer to tree sequence file with full path
     * from mrs/trunk/examples so that it works with Gopt examples
     */

    case 2:         // 2 taxa tree
                    // Human-Neandertal data in Sainudiin and York, AMB, 2009
                    // Illustrates identifibality of the sum of two branches
                    // Similar to the 2 taxa example in Sainudiin, Ph.D., 2005,
                    // and that in Sainudiin and Yoshida, 2005
      topologies = TwoTaxaTop;
      n_dimensions = TwoTaxaDim;
      treeFileName = "../../targets/TreesSeqns/testtree2";
      seqFileName = "../../targets/TreesSeqns/testseq2";
      break;
      /*
      // another GOpt example with minimal sufficient sequence patterns
      // Human-Neandertal-Chimp Example Sainudiin and York, AMB, 2009
    case 3:         // 3 taxa tree (This is truly unrooted Tree as 
                    // opposed to the constrained unrooted tree in
                    // Sainudiin and York, AMB, 2009.)
      topologies = ThreeTaxaTop;
      n_dimensions = ThreeTaxaDim;
      treeFileName = "../../targets/TreesSeqns/testtree3";
      seqFileName = "../../targets/TreesSeqns/testseq3";
      break;
      */
    case 3: // 3 taxa tree Brown et al, 1982 example
            // Used in Yang 2000, Sainudiin, Ph.D., 2005, 
            // Sainudiin and Yoshida, 2005,
            // Sainudiin and York, AMB, 2009
      topologies = ThreeTaxaTop;
      n_dimensions = ThreeTaxaDim;
      treeFileName = "../../targets/TreesSeqns/testtree3";
      seqFileName = "../../targets/TreesSeqns/testseq3Brown"; 
      break;
      
    case 4:         // 4 taxa tree
      topologies = FourTaxaTop;
      n_dimensions = FourTaxaDim;
      treeFileName = "../../targets/TreesSeqns/testtree4";
      //seqFileName = "../../targets/TreesSeqns/testseq4";
      //seqFileName = "../../targets/TreesSeqns/testseq4bc2";
      seqFileName = "../../targets/TreesSeqns/testseq4ChorB100";
      break;

    case 5:         // 5 taxa tree -- way too long in this implementation
      topologies = FiveTaxaTop;
      n_dimensions = FiveTaxaDim;
      treeFileName = "../../targets/TreesSeqns/testtree5";
      seqFileName = "../../targets/TreesSeqns/testseq5";
      break;

    default:

      std::cerr<< "Does NOT contain 2, 3, 4 or 5 taxa";
      exit(1);

  }                 // end switch

  // for each topology
  // set up a labelled box of n_dimensions and put into the LabDomainList
  // LabDomainList is an inherited member from Fobj
  for (int j=0; j < topologies; j++)
  {

    LabBox Ldomain;
    ivector domain (1, n_dimensions);
    for (int i = 1; i <= n_dimensions; i++)
    {
      domain[i] = Domain;
    }

    Ldomain.L = j;
    Ldomain.Box = domain;

    LabDomainList.push_back(Ldomain);
  }

  // number of nodes is number of dimensions + 1;
  tree_nodes = n_dimensions + 1;

  //reads in the sequence file and initializes SeqData
  ReadSequence(seqFileName);

  //finds patterns and initializes them in SeqData.Baseseq
  FindPattern();

  PrintSequence();

  // read the tree file, which should have one line for each topology
  // and tries to make a tree from each line
  // trees made are pushed back into treeRoots
  ReadTrees(treeFileName);

  // print out the trees, one for each topology
  PrintTopologyTrees();

}                   // end constructor

// Destructor
FPhyloPOT::~FPhyloPOT()
{
  destroyRoots();   // makes sure the newed tree roots are destroyed

}                   // end destructor

void FPhyloPOT::destroyRoots() const
{
  // call DeleteObject to delete all the newed tree roots
  for_each(treeRoots.begin(), treeRoots.end(), DeleteObject());
  // and clear the container of pointers as well
  treeRoots.clear();
}

// -------------------------------------- functions ---------------------


//this is the original interval function
interval FPhyloPOT::operator() (const LabBox& lb) const
{
                    // no patterns or no treeRoots
  if(SeqData.BNo_pattern == 0 || treeRoots.empty())
  {
    cerr << "Function object contains either no patterns OR no trees" 
         << std::endl;
    exit(1);
  }

  ivector x(lb.Box);// the box we should be dealing with

  // identify which tree in the topology space from the label of the box
                    // labels run 0, 1, 2, etc
  int boxLabel = lb.L;
  if(boxLabel >= static_cast<int>(treeRoots.size()))
  {
    std::cerr << "Label on point " << boxLabel <<
      " greater outside bounds of available trees ("
      << treeRoots.size() << ")" << std::endl;
    exit(1);
  }
  PhyloTree* thisTree = treeRoots[boxLabel];

                    // for return value, to accumulate the likelihood we seek
  interval lkl(0.0, 0.0);

 //pure C0 enclosure only 
 if(!CenteredForm)
 {                 
  // for each pattern in the pattern counts
  for (int pattern_i = 0; pattern_i < SeqData.BNo_pattern; pattern_i++)
  {
			// set up a container for nucleotype probabilities
    IntervalProbs iNucProb;
    // use fillProbIntervalJC69 to fill this in for the pattern_ith member of 
    // the baseSequences so that it will have one member for each character in 
    // the CharacterSpace, each being
    // the probability of that character at the root node
    iNucProb = 
      thisTree->fillProbInterval(x, ((SeqData.baseSequences)[pattern_i]), 
                                     iNucProb);

    // check we have got the right number of probabilities back
    if(iNucProb.size() != static_cast<size_t>(CharacterSpace))
    {
      std::cerr << "Error, number of nucleotide probabilities "
                << "less than character space" << std::endl;
      exit(1);
    }

    interval prob(0.0, 0.0);
    interval initialvalue(0.0, 0.0);
    // combine the nucleotide probabilities
    // sum with the intervalSum function
    prob = std::accumulate(iNucProb.begin(), iNucProb.end(), 
                           initialvalue, intervalSum);

    // eg ln of (prob mult by 0.25 for JC69 4-character set/ by 0.5 for CFN 2-state set)
    prob = ln(prob/CharacterSpace);

    // add to the likelihood
    lkl = lkl + prob * (SeqData.BPatternCounts)[pattern_i];
  }

  // increment the number of interval() calls in the function object
  n_interval_calls++;
 }//end of pure C0 enclosure only 
//-------------------------new stuff to get centred form
 else
 {
 //interval fX, fC;
 rvector c = mid(x);
  if(0)// getting centered form range enclosure using hess_ari -- MUCH SLOWER!!
  {//This is mainly here for syntax and if one is doing hess_ari anyways...
   HessType hh;
   hh = (*this)(HessVar(c), boxLabel);
   interval lklc = fValue(hh);
   hh = (*this)(HessVar(x), boxLabel);
   lkl = fValue(hh);
   ivector Gradlkl = gradValue(hh);
   lkl = (lklc + Gradlkl*(x - c)) & lkl;
  }
  else// getting centered form range enclosure using grad_ari -- MUCH FASTER!!
  {
   GradType gg(Ub(x)-Lb(x)+1);
   gg = (*this)(GradVar(c), boxLabel);
   interval lklc = fValue(gg);
   //the above should be made faster by a interval evaluation instead
   gg = (*this)(GradVar(x), boxLabel);
   lkl = fValue(gg);
   ivector Gradlkl = gradValue(gg);
   lkl = (lklc + Gradlkl*(x - c)) & lkl;
  }
 }
//-------------------------new stuff to get centred form ends

 return lkl;
}

// the function's () operator
// version taking labeled point as parameter
real FPhyloPOT::operator() (const LabPnt& lp) const
{
                    // no patterns or no treeRoots
  if(SeqData.BNo_pattern == 0 || treeRoots.empty())
  {
    cerr << "Function object contains either no patterns OR no trees" 
         << std::endl;
    exit(1);
  }

  rvector x(lp.Pnt);// the point we should be dealing with

  // identify which tree in the topology space from the label of the point
                    // labels run 0, 1, 2, etc
  int ptLabel = lp.L;
  if(ptLabel >= static_cast<int>(treeRoots.size()))
  {
    std::cerr << "Label on point " << ptLabel <<
      " greater outside bounds of available trees ("
      << treeRoots.size() << ")" << std::endl;
    exit(1);
  }
  PhyloTree* thisTree = treeRoots[ptLabel];

  real lkl(0.0);    // for return value, to accumulate the likelihood we seek

  // for each pattern in the pattern counts
  for (int pattern_i = 0; pattern_i < SeqData.BNo_pattern; pattern_i++)
  {

                    // set up a container for nucleotype probabilities
    RealProbs rNucProb;
    // use fillProbRealJC69 to fill this in for the pattern_ith member of the 
    // baseSequences so that it will have one member for each character in the 
    // CharacterSpace, each being the probability of that character at the root 
    // node
    rNucProb = 
      thisTree->fillProbReal(x, ((SeqData.baseSequences)[pattern_i]), 
                                 rNucProb);

    // check we have got the right number of probabilities back
    if(rNucProb.size() != static_cast<size_t>(CharacterSpace))
    {
      std::cerr << "Error, number of nucleotide probabilities "
                << "less than character space" << std::endl;
      exit(1);
    }

    real prob(0.0);

    // combine the nucleotide probabilities
                    // sum with the realSum function
    prob = std::accumulate(rNucProb.begin(), rNucProb.end(), 
                           static_cast<real>(0.0), realSum);

    // eg ln of (prob mult by 0.25 for JC69 4-character set)
    prob = ln(prob/CharacterSpace);

    // add to the likelihood
    lkl = lkl + prob * (SeqData.BPatternCounts)[pattern_i];
  }

  // increment the number of real() calls in the function object
  n_real_calls++;   

  return lkl;
}

GradType FPhyloPOT::operator() (const GTvector& x, const int label) const
{

                    // no patterns or no treeRoots
  if(SeqData.BNo_pattern == 0 || treeRoots.empty())
  {
    cerr << "Function object contains either no patterns OR no trees" 
         << std::endl;
    exit(1);
  }

  // check the label is in our label set
  if (label >= static_cast<int>(treeRoots.size()))
  {
    std::cerr << "Tree label " << label
      << " does not identify a tree for this function object" << std::endl;
    exit(1);
  }

  // We do global optimisation for the tree labelled with label
  PhyloTree* thisTree = treeRoots[label];

  int d = x.Dim();
  GradType lkl(d);
  lkl = 0.0;        // for return value, to accumulate the likelihood we seek

  // for each pattern in the pattern counts
  for (int pattern_i = 0; pattern_i < SeqData.BNo_pattern; pattern_i++)
  {

    // set up a container for nucleotype probabilities
    GradProbs gNucProb;
    // use fillProbGradJC69 to fill this in for the pattern_ith member of the 
    // baseSequences so that it will have one member for each character in the 
    // CharacterSpace, each being
    // the probability of that character at the root node
    gNucProb = 
      thisTree->fillProbGrad(x, ((SeqData.baseSequences)[pattern_i]), 
                                 gNucProb);

    // check we have got the right number of probabilities back
    if(gNucProb.size() != static_cast<size_t>(CharacterSpace))
    {
      std::cerr << "Error, number of nucleotide probabilities less than " 
                << "character space" << std::endl;
      exit(1);
    }

    GradType prob(d);
    prob = 0.0;
    GradType initialvalue(d);
    initialvalue = 0.0;

    // combine the nucleotide probabilities
                    // sum with the gradSum function
    prob = std::accumulate(gNucProb.begin(), gNucProb.end(), 
                           initialvalue, gradSum);

    // eg ln of (prob mult by 0.25 for JC69 4-character set)
    prob = ln(prob/CharacterSpace);

    // add to the likelihood
    lkl = lkl + prob * (SeqData.BPatternCounts)[pattern_i];
  }

  // increment the number of gradtype() calls in the function object
  n_gradtype_calls++;   

  return lkl;
}

// the function's () operator
// version taking a HTvector as parameter
// only uses the first tree if the is more than one tree to choose from
// label defaults to 0 and identified which tree to do global optimisation on
HessType FPhyloPOT::operator() (const HTvector& x, const int label) const
{

                    // no patterns or no treeRoots
  if(SeqData.BNo_pattern == 0 || treeRoots.empty())
  {
    cerr << "Function object contains either no patterns OR no trees" 
         << std::endl;
    exit(1);
  }

  // check the label is in our label set
  if (label >= static_cast<int>(treeRoots.size()))
  {
    std::cerr << "Tree label " << label
      << " does not identify a tree for this function object" << std::endl;
    exit(1);
  }

  // We do global optimisation for the tree labelled with label
  PhyloTree* thisTree = treeRoots[label];

  int d = x.Dim();
  HessType lkl(d);
  lkl = 0.0;        // for return value, to accumulate the likelihood we seek

  // for each pattern in the pattern counts
  for (int pattern_i = 0; pattern_i < SeqData.BNo_pattern; pattern_i++)
  {

    // set up a container for nucleotype probabilities
    HessProbs hNucProb;
    // use fillProbHessJC69 to fill this in for the pattern_ith member of the 
    // baseSequences so that it will have one member for each character in the 
    // CharacterSpace, each being
    // the probability of that character at the root node
    hNucProb = 
      thisTree->fillProbHess(x, ((SeqData.baseSequences)[pattern_i]), 
                                 hNucProb);

    // check we have got the right number of probabilities back
    if(hNucProb.size() != static_cast<size_t>(CharacterSpace))
    {
      std::cerr << "Error, number of nucleotide probabilities less than " 
                << "character space" << std::endl;
      exit(1);
    }

    HessType prob(d);
    prob = 0.0;
    HessType initialvalue(d);
    initialvalue = 0.0;

    // combine the nucleotide probabilities
                    // sum with the hessSum function
    prob = std::accumulate(hNucProb.begin(), hNucProb.end(), 
                           initialvalue, hessSum);

    // eg ln of (prob mult by 0.25 for JC69 4-character set)
    prob = ln(prob/CharacterSpace);

    // add to the likelihood
    lkl = lkl + prob * (SeqData.BPatternCounts)[pattern_i];
  }

  return lkl;
}

// Converts the bases into numerical codes
// Gives out error message if it's not one of the T,C,A,G
int FPhyloPOT::Char2Code4(char ch)
{
                    //change to lower case
  char c = tolower(ch);
  int retValue = 0;
  switch (c)
  {
    case 't':
      retValue = 0;
      break;
    case 'c':
      retValue = 1;
      break;
    case 'a':
      retValue = 2;
      break;
    case 'g':
      retValue = 3;
      break;
    default:
      std::cerr<< "Contains non-nucleotide symbol " << c << "!" << std::endl;
      exit(1);

  }                 // end swtich

  return retValue;
}
// Converts the bases into numerical codes
// Gives out error message if it's not one of the T,C,A,G
// only purines and pyramidines {a,g} -> 0 {t,c}->1
int FPhyloPOT::Char2Code2(char ch)
{
                    //change to lower case
  char c = tolower(ch);
  int retValue = 0;
  switch (c)
  {
    case 't':
      retValue = 1;
      break;
    case 'c':
      retValue = 1;
      break;
    case 'a':
      retValue = 0;
      break;
    case 'g':
      retValue = 0;
      break;
    default:
      std::cerr<< "Contains non-nucleotide symbol " << c << "!" << std::endl;
      exit(1);
  }                 // end swtich

  return retValue;
}

char FPhyloPOT::Code2Char(int c)
{
  char retValue = 'x';

  switch (c)
  {
    case 0:
      retValue = 't';
      break;
    case 1:
      retValue = 'c';
      break;
    case 2:
      retValue = 'a';
      break;
    case 3:
      retValue = 'g';
      break;
    default:
      std::cerr<< "Contains non-nucleotide symbol " << c << "!" << std::endl;
      exit(1);
  }                 // end switch

  return retValue;
}

void FPhyloPOT::CheckReadLine(ifstream& ifs, string& line)
{
  if(!ifs.fail() && !ifs.bad())
  {
    getline(ifs, line);
  }
  else
  {
    FileInputError(ifs, " Error reading from file ");
  }
}

void FPhyloPOT::FileInputError(ifstream& ifs, const string& msg)
{
  std::cerr << msg << endl;
  ifs.close();
  exit(-1);
}

void FPhyloPOT::ReplaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) 
{
  size_t pos = 0;
  while((pos = subject.find(search, pos)) != std::string::npos) {
     subject.replace(pos, search.length(), replace);
     pos += replace.length();
  }
}

// reformat a line
// take out white space and \0, replace u with t, and make it all lower case
void FPhyloPOT::ReformatSequence(string& line)
{

  size_t pos;

  string toErase(" \r\t\0");
  string toReplace("u");

  // erase unwanted characters
  pos = line.find_first_of(toErase);
  while (pos!=string::npos)
  {
    line.erase(pos);
    pos = line.find_first_of(toErase, pos);
  }

  //turn it all to lower case
  std::transform(line.begin(), line.end(), line.begin(), ::tolower);

  // if a character is 'u', replace by 't'
  pos = line.find_first_of(toReplace);
  while (pos!= string::npos)
  {
    line.replace(pos, 1, "t");
    pos = line.find_first_of(toReplace, pos+1);
  }
}

// read in sequences from a given file
// fills in part of the data member SeqData which is a struct of type data
// returns an int 0 when completed
int FPhyloPOT::ReadSequence(const string& s)
{
  // reset
  SeqData.No_seq = 0;
  SeqData.Seq_length = 0;
  // make sure the sequence name and sequence containers are empty
  (SeqData.seqNames).clear();
  (SeqData.rawSequences).clear();

                    // make and open an ifstream object
  ifstream seqFile(s.c_str());

  int i;            // for counters

  string line;      // for a line read in

  // read the first line
  CheckReadLine(seqFile, line);

                    // convert to a stream
  istringstream sin(line);
  // the first line which should be the number of sequences 
  // and the sequence length

                    // extraction operator on a stream
  sin >> SeqData.No_seq >> SeqData.Seq_length;

  //  check we read two positive two integers from the first line
  if ((SeqData.No_seq <= 0) || SeqData.Seq_length <= 0)
  {
    FileInputError(seqFile, "Check your sequence data file first line");
  }

  // for each of the expected number of sequences
  // if there are more sequences they will not be read
  // and if there are less the program will give error message and exit
  // it would be easy to have a more flexible format without having to 
  // pre-specify number of sequences and length
  for (i = 0; i < SeqData.No_seq; i++)
  {
    if (!seqFile.eof())
    {
      //get the next line in the file, which should be a sequence name
      CheckReadLine(seqFile, line);

      // check line has some characters and less than 
      // maxCharInSeqName characters
      if (line.empty() || line.size() > maxCharInSeqName)
      {
        FileInputError(seqFile, 
                       "Check sequence file: expected a sequence name");
      }
      else          // line is okay
      {
        (SeqData.seqNames).push_back(line);
      }
    }
    else
    {
      FileInputError(seqFile, 
      "Check sequence file: found less than expected number of sequence names");
    }

    if (!seqFile.eof())
    {
      //get the next line in the file, which should be a sequence
      CheckReadLine(seqFile, line);

      // check line has some characters
      if (line.empty())
      {
        FileInputError(seqFile, "Check sequence file: expected a sequence");
      }
      else          // line is okay
      {

        //check the line

                    // reformats the sequence line
        ReformatSequence(line);
        // and put it into the container of sequences

        if(line.size() == static_cast<unsigned int>(SeqData.Seq_length))
        {
          // copy the string into vector of chars and copy this 
          // into the rawSequences
          vector<char> tempCharVec(line.begin(), line.end());
          (SeqData.rawSequences).push_back(tempCharVec);
        }

        else
        {
          FileInputError(seqFile, 
                    "Check sequence file: a sequence is shorter than expected");
        }
      }
    }
    else
    {
      FileInputError(seqFile, 
           "Check sequence file: found less than expected number of sequences");
    }

  }                 // end of loop reading sequences and sequence names

  seqFile.close();

  return (0);
}

/* Finds the pattern in the sequences to save computational time */
int FPhyloPOT::FindPattern(void)
{
  int i, j;

  // make sure the vectors for sequences anlysis are clear
  (SeqData.baseSequences).clear();
  //  this destructs all the component vectors if any
  (SeqData.BPatternCounts).clear();

  (SeqData.BPatterns).clear();

  SeqData.BNo_pattern = 0;

  // multiple sequence alignment sets
  vector<string> msaSets;

  // form the msaSets
  for (i = 0; i < SeqData.Seq_length; i++)
  {
    // make a string from the ith element of the first sequences's data
    string s(1, ((SeqData.rawSequences)[0])[i]);
    for (j = 1; j<SeqData.No_seq; j++)
    {
      char add = ((SeqData.rawSequences)[j])[i];
      s+=add;       // add on the equivalent character from the other sequences
    }
    if (CharacterSpace==2)
    {
      ReplaceStringInPlace(s, "g", "a");
      ReplaceStringInPlace(s, "t", "c");
    }
    msaSets.push_back(s);
  }

  // should now have SeqData.Seq_length strings in msaSets

                    // sort the msaSets
  sort (msaSets.begin(), msaSets.end());

  vector<string>::iterator sit;

  int patternsFound = 0;

  sit = msaSets.begin();
  // vector<int> BPatternCounts;

  int seqPos = 0;   // to keep track of where we are in the sequence

  // while we have not found all the patterns
  while (seqPos < SeqData.Seq_length)
  {
    string found = *sit;

    SeqData.BNo_pattern++;
    (SeqData.BPatterns).push_back(found);

    int thisCount = 0;

    // move to the next non matching set  this works 
    // because msaSets is is sorted
    while (*sit==found && (seqPos < SeqData.Seq_length-1))
    {
      thisCount++;
      seqPos++;
      sit++;
    }

    // move to the next non matching set  this works 
    //bbecause msaSets is is sorted
    //if (*sit==found && (seqPos = SeqData.Seq_length -1))//Jenny's
    if (*sit==found && (seqPos == SeqData.Seq_length -1))
    {
      thisCount++;
      seqPos++;
    }

    (SeqData.BPatternCounts).push_back(thisCount);
    patternsFound += thisCount;

    string::iterator it;
    vector<int> tempvec;

    // put the found sequence in coded form into baseSequences
    for (it = found.begin(); it< found.end(); it++)
    {
      //char c = *it;
      if (CharacterSpace==2)
        tempvec.push_back(Char2Code2(*it));
      if (CharacterSpace==4)
        tempvec.push_back(Char2Code4(*it));
    }

    (SeqData.baseSequences).push_back(tempvec);

  }                 // return to the top of the loop

  return (0);
}

void FPhyloPOT::PrintSequence() const
{
  vector<int>::iterator it;

  //printing sequence patterns
  std::cout << std::endl;
  std::cout << "For the " << SeqData.No_seq << " sequences, there are " <<
    SeqData.BNo_pattern << " patterns " << std::endl;

  std::cout << std::endl;
  std::cout << "The pattern counts are:" << std::endl;
  for (int i = 0; i < SeqData.BNo_pattern; i++)
  {
    std::cout << (SeqData.BPatternCounts)[i] << std::endl;
  }

  std::cout << std::endl;
  //printing the nucleotide code patterns
  std::cout << "The patterns in code are:" << std::endl;
  for (int i = 0; i < SeqData.BNo_pattern; i++)
  {
    vector<int> tempvec = (SeqData.baseSequences)[i];
    for (it = tempvec.begin(); it < tempvec.end(); it++)
    {
      std::cout << *it;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
  //printing the nucleotide patterns
  std::cout << "The patterns in nucleotides are:" << std::endl;
  for (int i = 0; i < SeqData.BNo_pattern; i++)
  {
    std::cout << (SeqData.BPatterns)[i] << std::endl;
  }
}

// read trees from file
// and tries to make a tree out of each line
// trees are pushed back into tree roots
void FPhyloPOT::ReadTrees(const string& s)
{
                    // a container for treeLines read in
  vector<string> treeLines;

  // we have treefile with multiple lines, one line for each topology
  // read the file and store each line
                    // make and open an ifstream object
  ifstream treeFile(s.c_str());

  // get each line and push back into treeLines container
  while(treeFile.good())
  {
    string line;
    getline(treeFile, line);

    // erase trailing white space
    string whitespaces (" \t\f\v\n\r");
    size_t found;

    found=line.find_last_not_of(whitespaces);
    if (found!=string::npos)
    {
                    // erase everything from found onwards
      line.erase(found+1);
      treeLines.push_back(line);
    }
    else
    {
      line.clear(); // str is all whitespace
      std::cerr << "Line in treefile is all whitespace" << std::endl;
    }
  }

  treeFile.close(); // finished with the treeFile

                    // count the lines read in
  size_t countLines = treeLines.size();

  // checks on the treefile lines read in
  if (countLines == 0)
  {
    std::cerr << "Error in tree file: no trees read in" << std::endl << std::endl;
    exit(1);
  }
  if (countLines!=static_cast<size_t>(topologies))
  {
    std::cerr << "Error in tree file: number of lines for trees is not the "
              << "expected number of topologies " << topologies 
              << std::endl << std::endl;
    exit(1);
  }

  // should now have topologies lines read in 
  // (eg topologies = 1 for 3 taxa case, 3 for 4 taxa, etc)

  // for each line, ie each topology, set up a new tree for this topology
  for (size_t i = 0; i < static_cast<size_t>(topologies); i++)
  {

    PhyloTree* newtree;

    try
    {
      // use PhyloTree constructor with the 
      // line, sequence names, and character space
      newtree = new PhyloTree(CharacterSpace, treeLines[i], SeqData.seqNames);
    }

    catch (bad_alloc&)
    {
      std::cout << "Error allocating memory in ReadTrees" << std::endl;
    }

    // we expect to find the tree_nodes number of nodes
    if (static_cast<int>((newtree->getNumberNodes())) != tree_nodes)
    {
      // error in tree file

      string msg = "Error in tree file, please check the number of nodes";
      delete newtree;
      destroyRoots();
      exit(1);
    }

                    // pushes a copy of the newtree pointer into treeRoots;
    treeRoots.push_back(newtree);
  }
  // we now have treeRoots full of pointers to PhyloTrees
}

void FPhyloPOT::PrintTopologyTrees()
{

  //describe the generic tree
  std::cout << endl;
  std::cout << "The trees will have " << tree_nodes
    << " nodes and " << n_dimensions << " edges" << std::endl << std::endl;

  PhyloTreePtrsItr it;
  int i = 0;

  //for each tree in treeRoots
  for (it = treeRoots.begin(); it < treeRoots.end(); it++)
  {

    std::cout << "Tree number " << i << " is:" <<std::endl;
    (*it)->PrintTree();
    std::cout << std::endl << std::endl;
    i++;
  }
}

/*
 * Copyright (C) 2009 Raazesh Sainudiin and Gloria Teng
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

/*! \file CatalanCoeff.cpp
  \brief Definitions for Computing Catalan Coefficients and their Frequencies.
*/

// Purpose: To compute the Catalan Coefficients of a spatial binary branching 
// process based on the SAGE worksheet 
// (http://sage.math.canterbury.ac.nz/home/pub/15/) authored by 
// Raazesh Sainudiin, 2009.

#include "CatalanCoeff.hpp"

/*! Function to convert int state to string state
*/
inline std::string getStringState(vector<int> intState)
{
  string stateString;
  vector<int>::iterator iterIntState;

  iterIntState = intState.begin();
  while (iterIntState != intState.end())
  {
    int c =  *iterIntState++;
    stateString += to_string(c);
  }
  return stateString;
}

/*! Function to get next state after split
*/
inline vector<int> getNextState(vector<int> currSState, int nodeI)
{
  vector<int> nextSState;
  int i = nodeI;
  int j;
  int currSStateSize = currSState.size();
  int currNode, nextNode;
  currNode = currSState[i];
  nextNode = currNode + 1;
  //concatenate to get nextSState
  for (j=0; j <i; j++)
  {
    nextSState.push_back(currSState[j]); 
  //cout << nextSState[j] << endl;
  }

  nextSState.push_back(nextNode);
  nextSState.push_back(nextNode);

  for (j=(i+1); j < currSStateSize; j++)
  {
    nextSState.push_back(currSState[j]); 
  //cout << nextSState[j] << endl;
  }
  //while(!nextSState.empty())
  //nextSState.pop_back();
  return nextSState;
}


/*! Function to compute factorials
should use GSL factorial
*/
inline unsigned long int factorial(int num)
{
  if (num==1 || num==0)
    return 1;
  else
    return factorial(num-1)*num;
}

/*! Function to compute the Catalan number
*/
//factorial(2*k)/factorial(k+1)/factorial(k);
inline long double Catalan(int num)
{
  int endNum = num-1;
  int i;
  long double CatNum = 1;
  if(num==0 || num==1)
  return 1;
  for (i=0; i < endNum; i++)
  {
    CatNum *= (2.0*num - i)/(1.0*num - i);
  }
  return CatNum;
}

/*! Function to empty vector< vector<int> > States
*/
inline vector< vector<int> > emptyStates(vector< vector<int> > SStates)
{
  while ( !SStates.empty())
  SStates.pop_back();
  /*
  if (SStates.empty())
  cout << "state space emptied" << endl;
  */
  return SStates;
}

int main()
{
  int i; //for 'for' loops

  //level to stop splitting
  int L, L1; //this means each state has L leaves resulting from k splits
  cout << "Number of splits? ";
  cin >> L1;
  L = L1 + 1;

  //number of splits
  int k;

  //number of distinct states or Catalan number
  long double Csplit[L];

  //total number of possible paths to reach any state at level L
  unsigned long int totalPaths[L];

  //container for states
  States currStates, nextStates;
  States::iterator iterCurrStates;
  State currState, nextState;
  State::iterator iterCurrState;

  // a map for the Catalan Coefficients
  typedef map<State, int, LexicoSorting<State> > StatesMap;
  StatesMap StateCount;
  StatesMap::iterator iterStateCount;
  std::pair<StatesMap::iterator, bool> boolStateCount;

  //state space at 0 split
  State state0;
  state0.push_back(0);

  StateCount.insert(make_pair(state0,1));

  Csplit[0]=1;
  totalPaths[0]=1;

  //ofstream to output CatalanCoefficient and CatalanFrequencies to .txt.
  ofstream oss;

  clock_t start, end;
  double timeTaken;

  //Loop to get Catalan coefficients and frequencies using L-R leaf-depth encoding
  for (k=1; k < L; k++)  //for a given number of splits
  {
    //Set up string for filename to output the Catalan coefficients and frequencies
    string fileNameCC, fileNameCF;
    fileNameCC = "CatalanCoefficientSplit";
    fileNameCF = "CatalanFrequencySplit";
    std::ostringstream stm1;  //convert k to std::ostringstream

    start=clock();  //record the time taken to get the Catalan coefficient and frequencies 

    //Get number of nodes, Catalan number and total possible paths for each split
    Csplit[k] = Catalan(k); 
    totalPaths[k] = factorial(k);
    cout << "========================================" << endl;
    cout << "At " << k << " splits:" << endl;
    cout << Csplit[k] << " distinct states " << "\t" << totalPaths[k] << " possible paths " << endl;

    //Set up an empty map to store new states
    StatesMap newStateCount;
    StatesMap::iterator iterNewStateCount;
    std::pair<StatesMap::iterator, bool> boolNewStateCount;

    //for each state in current states map, StateCount
    for (iterStateCount=StateCount.begin(); iterStateCount != StateCount.end(); iterStateCount++)
    {
      currState = iterStateCount->first;  //current state in current states map
      //cout << "Current state:" << getStringState(currState) << endl;

      //for each node in current state
      for(int node = 0; node < currState.size(); node++)
      {
	//get state for each node in current state
	nextState = getNextState(currState, node);
	//cout << "Next state: " << getStringState(nextState) << endl;

	//insert the new state into newStateMap
	boolNewStateCount = newStateCount.insert(make_pair(nextState, 1));
	//Check if insertion is successful - a unique set will render a successful insertion.
	if(!(boolNewStateCount.second)) //if state is not unique
	{
	  //increment by the Catalan coefficient of the current state producing this next state
          newStateCount[nextState] += StateCount[currState];	
	  //cout << "this state already exist"  << endl;
	}  //end of "non-unique state"
	else
	{
          //set the count of this next state as the Catalan coefficient of current state
          newStateCount[nextState] = StateCount[currState];
          //cout << "this is a unique state" << endl;
	}  //end of "unique state"
      } //end of looping over node in current state
    } //end of looping through all Ck states in current states map

    cout << "There are " << newStateCount.size() << endl;
    //cout << " distinct states with the following breakdown:" << endl; 

    //Catalan Coefficients
    //Output states with corresponding counts in .txt.
    stm1 << k; //name file according to the number of splits k
    fileNameCC += stm1.str();
    fileNameCC += ".txt";
    oss.open(fileNameCC.c_str());

    //set up a vector for the number of paths leading to each state in newStateCount
    vector<int> StatePathVec;
    vector<int>::iterator iterStatePathVec;

    //set up a map for the frequencies corresponding to state paths
    map<int, int, less<int> > StateFreqMap;
    map<int, int, less<int> >::iterator iterStateFreqMap;
    std::pair<map<int, int, less<int> >::iterator, bool> boolStateFreqMap;
    int Path;

    for(iterNewStateCount=newStateCount.begin(); 
        iterNewStateCount != newStateCount.end(); iterNewStateCount++)
    {
      State thisState = iterNewStateCount->first;
      //cout << getStringState(thisState) << "\t\t" << iterNewStateCount->second << endl;
      oss << getStringState(thisState) << "\t" << (iterNewStateCount->second)<< "\n";

      //Find frequencies for state paths
      Path = iterNewStateCount->second;
      boolStateFreqMap = StateFreqMap.insert(make_pair(Path, 1));

      if (!(boolStateFreqMap.second)) //if there is more than one state with Path
      {
        StateFreqMap[Path] += 1;
      }
    }

    //close ostream for CatalanCoefficientSplitk.txt
    oss << flush;
    oss.close();
    cout << "State counts output to " << fileNameCC << endl;


    //output Catalan Frequencies to.txt
    fileNameCF += stm1.str();
    fileNameCF += ".txt";
    oss.open(fileNameCF.c_str());
    oss << "Paths" << "\t\t" << "States" << "\n";
    //Output to .txt
    for (iterStateFreqMap = StateFreqMap.begin(); iterStateFreqMap != StateFreqMap.end(); iterStateFreqMap++)
    {
      //cout << "There are " << iterStateFreqMap->second << " states with " 
      //     << iterStateFreqMap->first << " paths." << endl;
      oss << iterStateFreqMap->first << "\t\t" << iterStateFreqMap->second << "\n";
    }

    oss << flush;
    oss.close();
    cout << "State frequencies output to " << fileNameCF << endl;


    //Empty current states map
    StateCount.clear();

    if ( StateCount.empty())
    {
      //set nextStateCount map to be StateCount to continue loop
      //cout << "space emptied" << endl;
      StateCount = newStateCount;
    }
    else
    {
      cout << "space not emptied. Should not proceed." << endl;
      break;
    }

    end=clock();
    timeTaken = static_cast<double>(end-start)/CLOCKS_PER_SEC;
    cout << "Computing time : " <<timeTaken<< " s." << endl;
  }

  return 0;
}

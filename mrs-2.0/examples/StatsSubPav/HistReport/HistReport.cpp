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


/*! \file HistReport.cpp
\brief Testing histograms with histrogram reports distribution examples
*/
// to use std::vector
#include <vector>
// to use iterators
#include <iterator>
// to use valarray
#include <valarray>
#include<map>

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams
#include <unistd.h> // needed for getpid() call

#include "toolz.hpp"    // toolz headers
#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data

using namespace cxsc;
using namespace std;
using namespace subpavings;

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

int main()
{
    // ------- prepare to generate some data for the tests -----------

    // set up a random number generator for uniform rvs
    const gsl_rng_type * T;
    gsl_rng * r;
    //create a generator chosen by the environment variable GSL_RNG_TYPE
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    long s = time (NULL) * getpid();
    gsl_rng_set(r, s);


    string fileName; // create a name for the file to use

    // ----------------   example to create multiple histograms -------------------
    int numHist = 100000; // the number of histograms to make
    int minLeaves = 1; // for number of leaves pq splitting stopping criteria
    int MINLeaves = 4; // for number of leaves pq splitting stopping criteria

    // a map for counting visits to states of Histograms encoded
    // as OrderedLeafDepths
    typedef map<vector<int>,int,LexicoSorting<vector<int> > > OrdLeafDepthsMap;
    vector< OrdLeafDepthsMap > AnOrdLfDpMapVec(MINLeaves-minLeaves+1);
    //OrdLeafDepthsMap AnOrdLfDpMap;
    OrdLeafDepthsMap::iterator OrdLeafDepthsMap_Iter;
    std::pair<OrdLeafDepthsMap::iterator, bool> OrdLeafDepthsMap_bool;

    string outputFileName;// for output file
    ofstream oss;         // ofstream object
    oss << scientific;  // set formatting for input to oss
    oss.precision(5);
    bool successfulInsertion = false;
    bool successfulPQSplit = false;
    bool successfulPQMerge = false;
    bool Mixture_bool = false;//using mixture or not
    //----data generating partition----------------------------------------
        int d = 1; // dimension of the uniform hypercube to sample data from
        ivector pavingBox(d);

        // domain is a hypercube
        interval pavingInterval(0,1);
        for(int i=1; i <= d; i++) pavingBox[i] = pavingInterval;

        // data sampled as uniform mixture over leaves of sub-paving myPart
        AdaptiveHistogram myPart(pavingBox);
        // a container for the boxes
        vector<ivector> Pboxes;
        size_t PartSize;
if (Mixture_bool) {
        //not working bool successfulInstruction = myPart.splitToShape("0");// uniform
        bool successfulInstruction = myPart.splitToShape("1,1");// uniform mixture
        //bool successfulInstruction = myPart.splitToShape("2,2,1");// uniform mixture
        //bool successfulInstruction = myPart.splitToShape("3,3, 2, 1");
        cout << "result is " << successfulInstruction << endl;
        cout << "Level string for new partition is "
             << myPart.getLeafLevelsString() << endl;
        // we can keep splitting if the state is reachable
        /*
        if (successfulInstruction) {
          successfulInstruction = myPart.splitToShape("3, 4,4, 2,2, 3,3");
          cout << "result is " << successfulInstruction << endl;
          cout << "Level string for new partition is "
               << myPart.getLeafLevelsString() << endl;
        // We could also have put the data from the container into
        // the histogram with the splitToShape partition, no splitting here
        // successfulInsertion = myPart.insertFromRVec(theData);
        }
        */
        SPSnodePtrs Pleaves; // set up empty container for leaf node pointers
        SPSnodePtrsItr it; // and an iterator over the container
        myPart.getSubPaving()->getLeaves(Pleaves); // fill the container
        // container is filled by reading leaves off tree from left to right
        for(it = Pleaves.begin(); it < Pleaves.end(); it++) {
            // remember that it points to a pointer, so *it is still a ptr
            //get the counts in all the Pleaves
            //get the boxes from all the Pleaves
            Pboxes.push_back((*it)->getBox());
        }
        PartSize = Pboxes.size();
    //----end of data generating partition----------------------------------------
}
    for (int j = 1; j <= numHist; j++)
    { // loop to make histograms

        //cout << "Doing hist number " << j << endl;
        RVecData theData;   // a container for all the points generated
        // make a simulated data set allData to sample from
        int n = 4; // total points from random number generator
        // data sampled as uniform equi-mixture over leaves of sub-paving myPart
        for (int i = 0; i < n; i++) {
            rvector thisrv(d);
            if(Mixture_bool) {
              size_t RndBoxNum = floor(PartSize*gsl_rng_uniform(r));
              //cout << Pboxes[RndBoxNum] << endl;
              thisrv = DrawUnifBox(r,Pboxes[RndBoxNum]);
            }
            else {
              for(int i=1; i <= d; i++) {
                  thisrv[i]  = gsl_rng_uniform(r);
              }
            }

            // put points generated into container
            //cout << thisrv;
            theData.push_back(thisrv);
        }  // data  should be in theData

        // make an Adaptive Histogram object with a specified box
        AdaptiveHistogram myHist(pavingBox);

        // put the data from the container into the histogram, no splitting here
        successfulInsertion = myHist.insertFromRVec(theData);

        for(int L=minLeaves; L<=MINLeaves; L++)
        {

          CompCount nodeCompCount;
          CritLeaves_GTE critLeavesGTE(L);
          if (successfulInsertion) {
              // now split with priority queue
              // split node wth most pointsin first (compCount)
              // until leaves >= minLeaves (critLeaves_GTE)
              // no minPoints or minVolB limitations on splittable nodes
              successfulPQSplit = myHist.prioritySplit(nodeCompCount,
                                critLeavesGTE, NOLOG, r); // no logs
          }

          if (successfulPQSplit) {
            //
            // optional - if you want to get a txt output of each histogram
            // create a name for the file to output
            // fileName = "Hist";
            //convert j to a string
            // std::ostringstream stm2;
            // stm2 << j;
            // add the stringed j to the filename
            // fileName += stm2.str();
            // fileName += ".txt"; // and finish the filename
            // To realize a file output
            // myHist.outputToTxtTabs(fileName);
            //myHist.outputToTxtTabsWithEMPs(fileName);

            SPSnodePtrs leaves; // set up empty container for leaf node pointers
            SPSnodePtrsItr it; // and an iterator over the container
            myHist.getSubPaving()->getLeaves(leaves); // fill the container
            // container is filled by reading leaves off tree from left to right

            // a container for the counts
            IntVec counts;  // IntVec is a typedef for vector<int>
                            // the iterator is typedefed as IntVecItr
            // a container for the boxes
            vector<ivector> boxes;
            // a container for the volumes
            vector<double> volumes;
            // a valarray container for the node levels, sized to fit
            valarray<int> levels(spLeaves(myHist.getSubPaving()));

            int v=0;
//            for(it = leaves.begin(); it < leaves.end(); it++) {
//                // remember that it points to a pointer, so *it is still a ptr
//                //get the counts in all the leaves
//                counts.push_back((*it)->getCounter());
//                //get the boxes from all the leaves
//                boxes.push_back((*it)->getBox());
//                //get the volumes of all the leaves
//                volumes.push_back((*it)->nodeVolume());
//                //get the levels of the leaves
//                levels[v] = (*it)->getNodeDepth();
//                v++;
//            }
//
            IntVec altLevels = myHist.getLeafLevels();
//
//            copy (altLevels.begin(), altLevels.end(),
//                    ostream_iterator<int>(cout, "\t"));
//            cout << '\n';
//
            OrdLeafDepthsMap_bool = AnOrdLfDpMapVec[L-minLeaves].insert(make_pair(altLevels,1));
            if(!(OrdLeafDepthsMap_bool.second)) {
                (OrdLeafDepthsMap_bool.first)->second +=1;
            }
            //  This is where you'd have to use/manipulate/store/whatever
            //    the stuff you have from the containers for each hist
            //
          }//end op PQSplit
        }//end of leaf level pecific splitting loop
    } // end of loop for histograms
    // free the random number generator
    gsl_rng_free (r);


    cout << "empirical histogram frequencies\n";
    for(int L=minLeaves; L<=MINLeaves; L++)
    {
      cout << "L : " << '\t' << L << '\n';
      for(OrdLeafDepthsMap_Iter=AnOrdLfDpMapVec[L-minLeaves].begin();
          OrdLeafDepthsMap_Iter != AnOrdLfDpMapVec[L-minLeaves].end();
          ++OrdLeafDepthsMap_Iter){
              IntVec altLevels = OrdLeafDepthsMap_Iter->first;
              cout << OrdLeafDepthsMap_Iter->second << " : " << '\t';
              copy (altLevels.begin(), altLevels.end(),
                      ostream_iterator<int>(cout, ";"));
              cout << '\n';
      }
    }
    //  All the histograms done, you've somehow stored some summary from each,
    //  you can now do something with the summary

    return 0;

} // end of histrogram report test program

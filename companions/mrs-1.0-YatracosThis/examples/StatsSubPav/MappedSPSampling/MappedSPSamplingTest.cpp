/*! \file
\brief MappedSPnode example for Gaussian objects main.
*/

#include "RosenFobj2D.hpp"
#include "RosenFobj10D.hpp"
#include "GaussianFobj1D.hpp"
#include "GaussianFobj2D.hpp"
#include "GaussianFobj9D.hpp"
#include "GaussianFobj10D.hpp"
#include "GaussianFobj100D.hpp"
//#include "GaussianFobj1000D.hpp"

#include "mappedspnode.hpp"
#include "realmappedspnode.hpp"
#include "mappedspnodevisitor_expand.hpp"

#include "histall.hpp"  // headers for the histograms
#include "dataprep.hpp" // headers for getting data
#include "MCMCGRtools.hpp" // tools to help 

#include <time.h>   // clock and time classes
#include <fstream>  // input and output streams
#include <sstream>  // to be able to manipulate strings as streams
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions
#include <functional> // mutliplies<>
#include <algorithm> // transform

#include <gsl/gsl_randist.h> // to use gsl_ran_discrete_preproc
#include <valarray> 
#include "toolz.hpp" //draw unif box

#include "auto_tools.hpp"

//#define NDEBUG // uncomment this to turn off assertion checking and all extra debugging

#ifndef NDEBUG // ie only define these if we have not defined NDEBUG for no debugging
//#define MYDEBUG_OUTPUT // extra console output etc for debugging - only use for small examples!
//#define MYDEBUG_CALCS // extra console output for calculations
#define MYDEBUG // extra files for collations, averages and diffs to av as chains develop

//#define FORCEFAILINSERTION // debugging flag to force a failure during insertion of data

//#define FORCEFAILMCMCLOOP // debugging flag to force a failure during an MCMC loop

#endif

using namespace cxsc;
using namespace subpavings;
using namespace std;

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

//==========Functions for MappedSPnode===================================//
//to output MappedSPnode to .txt file
void output(string& filename,  const SPnode& node)
{
   // To generate a file output
   ofstream os(filename.c_str());         // Filename, c-string version
   if (os.is_open()) {
      node.leavesOutputTabs(os); // the output
      std::cout << "The output of the estimated function"
               << " has been written to " << filename << std::endl << std::endl;
         os.close();
      }
   else {
      std::cerr << "Error: could not open file named "
         << filename << std::endl << std::endl;
   }
}

//function to iterate through the leaves and get weights and boxes
void getAllWeights(RealMappedSPnode* thisNodePtr, vector<double>& WeightsVector,
							vector<interval>& WeightsInt)
{
	if (!(thisNodePtr->isEmpty()) && thisNodePtr->isLeaf()) { // this is a non-empty leaf
		 //get the weights
		 RangeCollectionClass<real> myContainer;
		 myContainer = thisNodePtr->getRangeCollection();
		 myContainer.getWeights(WeightsVector, WeightsInt, thisNodePtr->nodeVolume());
	}

	//recurse on the children
	if (thisNodePtr->hasLCwithBox()) {
		getAllWeights(thisNodePtr->getLeftChild(), WeightsVector, WeightsInt);
	}
   if (thisNodePtr->hasRCwithBox()) {
		getAllWeights(thisNodePtr->getRightChild(), WeightsVector, WeightsInt);
   }
}

//function to iterate through the leaves and get heights and boxes
void getHeightAndBox(RealMappedSPnode* thisNodePtr, vector<ivector>& BoxVector,
					 vector<real>& HeightsVector)
{
	if (!(thisNodePtr->isEmpty()) && thisNodePtr->isLeaf()) { // this is a non-empty leaf
		 //push back this box into the BoxVector
		 BoxVector.push_back(thisNodePtr->getBox());
		 
		 //get the heights
		 RangeCollectionClass<real> myContainer;
		 myContainer = thisNodePtr->getRangeCollection();
		 myContainer.getHeight(HeightsVector);
	}
  //recurse on the children
  if (thisNodePtr->hasLCwithBox()) {
		getHeightAndBox(thisNodePtr->getLeftChild(), BoxVector, HeightsVector);
	}
   if (thisNodePtr->hasRCwithBox()) {
		getHeightAndBox(thisNodePtr->getRightChild(), BoxVector, HeightsVector);
   }
}

// normalize the heights
void normHeights(RealMappedSPnode* thisNodePtr, double totalArea, 
					vector< RangeCollectionClass<real> >& heightNorm)
{
	if (!(thisNodePtr->isEmpty()) ) { // this is non-empty
		 RangeCollectionClass<real> myContainer;
		 myContainer = thisNodePtr->getRangeCollection();
		 real newHeight = myContainer.normNodeHeight(totalArea);
		 
		 RangeCollectionClass<real> height(newHeight);
		heightNorm.push_back(height);
	}
  //recurse on the children
  if (thisNodePtr->hasLCwithBox()) {
		normHeights(thisNodePtr->getLeftChild(), totalArea, heightNorm);
	}
   if (thisNodePtr->hasRCwithBox()) {
		normHeights(thisNodePtr->getRightChild(), totalArea, heightNorm);
   }
}
//=======================end of functions====================================//

int main(int argc, char* argv[])
{
	//========user-defined parameters====================//
	size_t n=atoi(argv[1]);  // number of datapoints to generate for each histogram
	int d = atoi(argv[2]); // dimensions
	size_t numHist = atoi(argv[3]); // number of repetitions for simulation purposes
	
	//	for generating samples from MappedSPnode 
	// ensure max leaves is < 1E6 or something reasonable
	size_t maxLeaves = atoi(argv[4]);
	
	// for the MCMC run
	int maxLoops = atoi(argv[5]); // maximum changes of state from initial state to try
	int samplesNeeded = atoi(argv[6]); // how many samples do we want (ie once chains have burned in)
	int thinout = atoi(argv[7]); // sample every thinout state, ie thinout-1 states between samples
	
	real tolerance = atof(argv[8]);
	cxsc::real tol(tolerance); //tolerance for automated burn in criteria
	
	size_t minPoints = atoi(argv[9]); 

	int dataSeed = atoi(argv[10]);
	
	double maxLeaf = atof(argv[11]);
	
	// should really do more checks on parameters, but just check thinout here
	if (thinout < 1 ) {
		throw std::invalid_argument("Invalid thinout argument");
	}

	// use the cxsc manipulators for changing printing of cxsc::reals to console
	int prec = 15;
	cout << cxsc::SaveOpt;
	cout << cxsc::Variable;
	cout << cxsc::SetPrecision(prec+2, prec);

	//string formatting
	ofstream oss;         // ofstream object
   oss << scientific;  // set formatting for input to oss
   oss.precision(10);

	//=========set up to estimate the function==============================// 
	// Function object
//	GaussianFobj1D realF;
//	GaussianFobj2D realF;
	 GaussianFobj9D realF;
//	 GaussianFobj10D realF;
	//RosenFobj2D realF;
	//RosenFobj10D realF;

	//make a root box
	ivector pavingBox(d);
//	interval pavingInterval(-3,3);
//	interval pavingInterval(-10,10);
	interval pavingInterval(-6.5,7.5);
	for(int k=1; k <= d; k++) pavingBox[k] = pavingInterval;
	
	 RealMappedSPnode nodeEst(pavingBox); // make a MappedSPnode object
    // estimate the function
	 MappedSPnodeVisitorExpand expander(realF, 0);
	 
	 vector<real> epsVec;
	 
	 //RealMappedSPnode nodeEst1(pavingBox);
	 //nodeEst1.accept(expander);
	 
	 nodeEst.priorityAccept(expander, maxLeaves, epsVec);

 	string avgL1FileName = "Eps";
	avgL1FileName += ".txt";
	oss.open(avgL1FileName.c_str());
		for (size_t i = 0; i < epsVec.size(); i++) { 
			//cout << epsVec[i] << endl;
			oss << epsVec[i] << "\n";
		}
		oss << flush;
		oss.close();

 	 cout << "Estimate function has " << nodeEst.getNumLeaves() << " leaf nodes." << endl;
	 
	 //RealMappedSPnode nodeEst1;
	 //nodeEst1 = nodeEst;
	 
	 //RealMappedSPnode diff = nodeEst1 - nodeEst;
	 
	 //output to .txt  
 	//tring thefilename = "Est.txt";
	//utput(thefilename, nodeEst);

	//=======================================================================//
	
	//==================Get the weights of the boxes=========================//
	 cout << "Getting boxes and weights:" << endl;
    vector<ivector> BoxVector;
	 vector<real> HeightsVector;
	 RealMappedSPnode* nodePtr;
	 nodePtr = &nodeEst;
	 vector<double>* WeightsVectorPtr;
	 WeightsVectorPtr = new vector<double>;
	 vector<interval>* WeightsIntPtr;
	 WeightsIntPtr = new vector<interval>;
	 
	 // iterate through the leaf nodes to get boxes and heights and weights
	 getHeightAndBox(nodePtr, BoxVector, HeightsVector);
	 getAllWeights(nodePtr, *WeightsVectorPtr, *WeightsIntPtr);
	 
	 //now put elements of WeightsVector into an array of doubles
	 size_t sizeWeight =(*WeightsVectorPtr).size();
	 //check that number of boxes < 10^6
	 if (sizeWeight > pow(10,7)) { 
			cerr << "Too many boxes (" << sizeWeight << ")." << endl;
			exit(1);
	 }

	// normalize using heights
	interval areaInt = interval(0);
	 //normalize the heights so that the function integrates to 1
	 for (size_t i = 0; i < sizeWeight; i++) {
		areaInt = areaInt + (*WeightsIntPtr)[i];
		
		//cout << (*WeightsVectorPtr)[i] << "\t" << (*WeightsIntPtr)[i] << endl;
	 }
	cout << "Total area: " << mid(areaInt) << endl; 
	 
	 double totalArea = _double(mid(areaInt));
	 
	 // very important - normalize the heights in nodeEst
	 vector< RangeCollectionClass<real> >* heightNorm = new vector< RangeCollectionClass<real> >;
	 normHeights(nodePtr, totalArea, *heightNorm);
	 nodeEst.allocateRanges(*heightNorm, 0);
	 string filename = "EstFunctionAfterNormalized.txt";
	 output(filename, nodeEst);
	 
	 //need to check that the weights equal to 1
	 double densityCheck = 0.0;
    //convert vector to array
	 double WeightsArray[sizeWeight];
	 for (size_t i = 0; i < sizeWeight; i++) {
			WeightsArray[i] = (*WeightsVectorPtr)[i];
			densityCheck += (*WeightsVectorPtr)[i]/totalArea;
	 }
	 
	cout << "Total area after normalizing: " << densityCheck << endl;
	/*
	if ( (densityCheck != 1.0) ) {
		cout << densityCheck << endl; 
		cerr << "Function does not integrate to 1. Need to normalize." << endl; 
		exit(0); 
	}
	*/
	
	 //return to the system the memory that was pointed to by WeightsVectorPtr
	 // and heightNorm
	 delete WeightsVectorPtr;
	 delete heightNorm;
	 
    //now get the lookup table returned from gsl_ran_discrete_preproc
	 //i.e. the box indices with their weights
	 gsl_ran_discrete_t* gslpdfstruct;
	 gslpdfstruct = gsl_ran_discrete_preproc(sizeWeight, WeightsArray);
	 //===================end of getting box weights=======================//

	//===========preliminaries before simulations========================//
	// set up a random number generator to draw from weighted boxes
	const gsl_rng_type * T;
	gsl_rng * r;

	//create a generator chosen by the environment variable GSL_RNG_TYPE
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	// set a seed for the data
	
	//===========end of setting up preliminaries=======================//

	for (int k = 1; k <= numHist; k++) {
		cout << "Data set " << k << endl;	
		dataSeed = k;	
	
		gsl_rng_set(r, dataSeed);

	//-------------generate data--------------------------------------//
	//now sample n data points from boxes given by the proposed indices
	cout << "Sample data points using weighted boxes:" << endl;
	RVecData theData;   // a container for all the points generated
  // make a simulated data set
	// data sampled as weighted-uniform-mixtures
	for (size_t i = 0; i < n; i++) {
		rvector thisrv(d);
		size_t proposedIndex = gsl_ran_discrete(r, gslpdfstruct);
		//int proposed_index = static_cast<int>(gsl_ran_discrete(r, gslpdfstruct));
		thisrv = DrawUnifBox(r, BoxVector[proposedIndex]);
		// put points generated into container
		theData.push_back(thisrv);
	}  // data  should be in theData
	
	cout << (theData).size() << " points generated" << endl;
	
	string dataFileName = "MappedData";
	ostringstream stm;
	stm << dataSeed;
	dataFileName += stm.str(); 
	dataFileName += ".txt";
	
	oss.open(dataFileName.c_str());
	for (size_t i = 0; i < n; i++) { 
		for (size_t j = 1; j <= d; j++) {
			//cout << (theData)[i][j] << "\t"; 
			oss << (theData)[i][j] << "\t";
		}
		oss << "\n";
		//cout << "\n";
	}
	oss << flush;
	oss.close();
	
	cout << "Mapped data written to  " << dataFileName << endl;

	//=================generate Gaussian data========================

	cout << "Generating Gaussian data: " << endl;
	
	const gsl_rng_type * T1;
	gsl_rng * r1;
	gsl_rng_env_setup();
	T1 = gsl_rng_default;
	r1 = gsl_rng_alloc (T1);
	gsl_rng_set(r1, dataSeed);

	RVecData actualData;
	
	for (size_t i = 0; i < n; i++) {
		rvector thisrv(d);
		for (size_t j = 1; j <= d; j++) {
			double z = gsl_ran_gaussian(r1, 1.0); // generate a normal r.v.
			thisrv[j] = _real(z);
		}
		//cout << thisrv << endl;
		actualData.push_back(thisrv);
	}

	cout << (actualData).size() << " points generated" << endl;
	
	dataFileName = "ActualData";
	dataFileName += stm.str(); 
	dataFileName += ".txt"; 
	oss.open(dataFileName.c_str());
	for (size_t i = 0; i < n; i++) { 
		for (size_t j = 1; j <= d; j++) {
				//cout << (actualData)[i][j] << "\t";
				oss << (actualData)[i][j] << "\t";
		}
		oss << "\n";
		//cout << "\n";
	}
	oss << flush;
	oss.close();

	cout << "Actual data written to  " << dataFileName << endl;
	
	}
	
	gsl_rng_free(r);
	//gsl_rng_free(r1);
	gsl_ran_discrete_free (gslpdfstruct);
	
	return(0);
	
} // end of MCMC test program

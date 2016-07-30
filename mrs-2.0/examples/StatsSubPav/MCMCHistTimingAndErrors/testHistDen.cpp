/*
* Copyright (C) 2014 Jennifer Harlow
*
*/


/*! \file
\brief Example for an mvn mixture distribution
making hist with MCMCGR multiple chains, no kde, 

showing how to make log files for timing, leaves, errors, etc
 */

#include "testDenCommon.hpp"
#include "testDenTools.hpp"

#include "mixture_mvn.hpp"



#include <vector>
#include <string>
#include <iostream>  // input and output streams
#include <fstream>  // file streams
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept> // throwing exceptions
#include <unistd.h>


using namespace std;
using namespace subpavings::kde;


void test1();



int main(int argc, char* argv[])
{
	
		test1();
		
} 


void test1()
{
    /* dimensions of data */
	size_t dim = 2;
	
	long unsigned int seed = 1234;
	
	std::vector < size_t > n_vec; 
	{
        /* number of data points in sample to make hist from */
		//size_t tmp[] = {5000, 10000, 50000, 100000,500000,1000000}; //, ,500000,1000000
		size_t tmp[] = {100}; 
		n_vec.insert (n_vec.begin(), tmp, tmp+sizeof(tmp) / sizeof(size_t));
	}
	
    /* replications (for looking at timings, leaves, etc, over */
	int reps = 3;
	
	size_t minPoints = 1;
	
	/* points to use to try to assess approximation accuracy - the errors*/
	size_t intNbase = 1000000;
	size_t intN = (dim > 2 ? dim : 1) * intNbase;
	
	cout << "pid = " << getpid() << endl;
	

	std::string baseOutputDir("output/");
	
	std::string thisDir;
	{
		ostringstream oss;
		oss << "HistDen" << dim << "DAlt";
		thisDir = oss.str();
	}
	
	// make the output dir and get back path
	std::string path = makeDir(baseOutputDir, thisDir);
	
	/* empty here (not used) */
	std::vector < std::vector < int > > vecSliceDims;
	std::vector < std::vector < double > > vecSlicePts;
	
	/* setting up filenames to use for output of histograms and the log files
     * the program will add dimensions, number of samples, etc, to these base strings */
	string histFilenameBase = path +"Hist"; // base file name for histogams
	string logFilenameBase = path + "LogMixTwo"; // base file name for log file
	
	
	for (size_t ni = 0; ni < n_vec.size(); ++ni) {
		
		size_t n = n_vec[ni];
		
		cout << "\n" << dim << "-d mixture example to get histogram errors, n = " << n << endl;
		
		MixtureMVN* mixMVNptr = makeMixture(dim , seed); // seed gets changed in doDenEst

		doDenEst(
			dim,
			seed,
			reps,
			minPoints,
			n,
			intN,
			histFilenameBase,
			logFilenameBase,	
			mixMVNptr,
			vecSliceDims,
			vecSlicePts);
		

		/*delete the data generator */
		delete mixMVNptr;
		
	}
}	




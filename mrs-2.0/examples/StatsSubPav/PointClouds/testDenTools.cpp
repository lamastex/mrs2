/*
* Copyright (C) 2012 Jennifer Harlow
*
*/


/*! \file
\brief tools for comparing mcmc and kde

 */

#include "testDenTools.hpp"

#include "sptools.hpp"


#include <gsl/gsl_qrng.h>

#include <ctime>   // clock and time classes
#include <iostream>  // input and output streams
#include <fstream>  // file streams
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator>  // output iterator
#include <cassert> // for assertions
#include <stdexcept> // throwing exceptions

using namespace cxsc;
using namespace subpavings;
using namespace subpavings::kde;
using namespace std;

std::vector < std::vector < double > >& getQuasiRandomPointsUniformBox(
			std::vector < std::vector < double > >& intPts,
			size_t dim, size_t N);


void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < real > >& repResults);

void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < double > >& repResults);

void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < size_t > >& repResults);
					
void outputRepResultsMeans(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < double > >& repResults);


/* standard multivariate normal */
MixtureMVN* makeStandard(size_t dim, long unsigned int dataseed)
{
	vector < vector < double > > means;
	vector < vector < vector < double > > > scales;
	vector < double > mixes;
	
	{
		vector < double > mean(dim, 0.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double var = 1.0;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 0.0);
			scalerow[i] = var;
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(1.0);
	}
		
	MixtureMVN* mixMVNptr = new MixtureMVN(means, scales, mixes, dataseed);
	
	return mixMVNptr;
	
}



/* mixture of multivariate normals (if dim=2, as in Zhang, 2006, density A) */
MixtureMVN* makeMixture(size_t dim, long unsigned int dataseed)
{
	vector < vector < double > > means;
	vector < vector < vector < double > > > scales;
	vector < double > mixes;
	
	{
		vector < double > mean(dim, 2.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = -0.9;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.0);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.5);
	}
	{
		vector < double > mean(dim, -1.5);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = 0.3;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.0);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.5);
	}
		
	MixtureMVN* mixMVNptr = new MixtureMVN(means, scales, mixes, dataseed);
	
	return mixMVNptr;
	
}


/* mixture of 4 multivariate normals */
MixtureMVN* makeMixture4(size_t dim, long unsigned int dataseed)
{
	vector < vector < double > > means;
	vector < vector < vector < double > > > scales;
	vector < double > mixes;
	
	{
		vector < double > mean(dim, 4.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = -0.9;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.0);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.25);
	}
	{
		vector < double > mean(dim, 5.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = -0.7;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.3);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.25);
	}
	{
		vector < double > mean(dim, 7.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = -0.3;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.7);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.25);
	}
	{
		vector < double > mean(dim, 9.0);
		means.push_back(mean);
		
		vector < vector < double > > scale;
		scale.reserve(dim);
		double covar = 0.5;
		for (int i = 0; i < dim; ++i) {
			vector < double > scalerow(dim, 1.0);
			for (int j = 0; j < dim; ++j) {
				if (j != i) 
					scalerow[j] = std::pow(covar, (j > i ? j-i : i-j ));
				}
			
			scale.push_back(scalerow);
		}
		
		scales.push_back(scale);
		
		mixes.push_back(0.25);
	}
		
	MixtureMVN* mixMVNptr = new MixtureMVN(means, scales, mixes, dataseed);
	
	return mixMVNptr;
	
}



std::vector < real >& getTrueDensities(const MixtureMVN& mixMVN,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& trueDensities)
{
	std::vector < double > timing;
	return getTrueDensities(mixMVN,
						intPts,
						trueDensities,
						timing);
	
}

std::vector < real >& getTrueDensities(const MixtureMVN& mixMVN,
						const std::vector < std::vector < real > >& intPts,
						std::vector < real >& trueDensities,
						std::vector < double >& timing)
{
	
	size_t N = intPts.size();
	
	std::vector < real > tmp(N);
		
	clock_t start = clock();
	
	for (size_t i = 0; i < N; ++i) { 
	
		tmp[i] = mixMVN.f(intPts[i]);
		
	}
	
	clock_t end = clock();
	timing.push_back(static_cast<double>(end-start)/CLOCKS_PER_SEC);
		
	trueDensities.swap(tmp);
	
	return trueDensities;
	
}

real avLogDen(const std::vector < real >& densities)
{
	size_t N = densities.size();
	
	real sumLogDen(0.0);
		
	for (size_t i = 0; i < N; ++i) { 
	
		sumLogDen += cxsc::ln(densities[i]);
		
	}
	
	return ( N > 0 ? sumLogDen/(1.0*N) : real(0.0));
	
}

real avAbsDiffDen(const std::vector < real >& densities1,
				const std::vector < real >& densities2)
{
	size_t N = densities1.size();
	
	real sumAbsDenDiff(0.0);
		
	for (size_t i = 0; i < N; ++i) { 
	
		real denDiff = densities1[i] - densities2.at(i);
		if (denDiff < 0) sumAbsDenDiff -= denDiff;
		else sumAbsDenDiff += denDiff;
		
	}
	
	return ( N > 0 ? sumAbsDenDiff/(1.0*N) : real(0.0));
	
}




std::vector < std::vector < real > >& getQuasiRandomPoints(
			const cxsc::ivector& box,
			std::vector < std::vector < real > >& intPts,
			size_t N)
{
	size_t dim = VecLen(box);
	
	std::vector < std::vector < double > > tmpD;
	
	getQuasiRandomPointsUniformBox(tmpD, dim, N);
	
	std::vector < std::vector < real > > tmpR(N, std::vector <real >(dim));
	
	for (size_t j = 0; j < dim; ++j) {
		
		real len = diam(box[j+1]);
		real lb = Inf(box[j+1]);
		
		for (size_t i = 0; i < N; ++i) tmpR[i][j] = tmpD[i][j]*len + lb;

	}
 
	intPts.swap(tmpR);
	
	return intPts;
}

std::vector < std::vector < double > >& getQuasiRandomPointsUniformBox(
			std::vector < std::vector < double > >& intPts,
			size_t dim, size_t N)
{
	gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, dim);

	std::vector < std::vector < double > > tmp(N, std::vector < double> (dim, 0.0));
  
	for (size_t i = 0; i < N; ++i) 	{
		gsl_qrng_get (q, &(tmp[i][0]));
	}
 
	intPts.swap(tmp);
	gsl_qrng_free (q);
	
	return intPts;
}



std::string makeDir(const std::string& baseOutputDir, 
				const std::string& thisDir)
	{
		std::string dStr; // string identifying this directory
		{	
			std::ostringstream oss;
			oss << baseOutputDir << thisDir;
			dStr = oss.str(); 
		}
		
		/*command make the output directory (-p switch makes intermediates
		 *  as well, okay if already exists */
		std::string  mkdirCmd;
		{
			std::ostringstream oss;
			oss << "mkdir -p " 
				<< dStr;
			 mkdirCmd = oss.str(); 
		}	
		// and system call to make it
		system( mkdirCmd.c_str() );
		
		return dStr+"/";
	}

void outputResults(const std::string& logFilename,
			const std::vector < std::vector < double > >& timingMake,
			const std::vector < std::vector < double > >& timingIntDensities,
			const std::vector < std::vector < real > >& avLogDens,
			const std::vector < std::vector < real > >& avLogDenRatios,
			const std::vector < std::vector < real > >& estL1ErrorsQR,
			const std::vector < std::vector < size_t > >& leaves,
			const std::vector < std::vector < double > >& topt)
{
	outputFileStart(logFilename);
	
	{ // optimal temperature
		const std::string intro("Optimal Temperature");
		outputRepResults(logFilename, intro, topt);
	}
	{ // timings for make
		const std::string intro("Make\nrep\tMCMC\tKDEmcmc\tKDEnorm\tFunEst");
		outputRepResults(logFilename, intro, timingMake);
	}
	{ // leaves)
		const std::string intro("Leaves");
		outputRepResults(logFilename, intro, leaves);
	}
	{ // timings for qr densities
		const std::string intro(
			"Timings for densities of random points\nrep\tTrue\tMCMC\tKDEmcmc\tKDEnorm\tFunEst");
		outputRepResults(logFilename, intro, timingIntDensities);
	}
	{ // average log densities
		const std::string intro("Average log densities");
		outputRepResults(logFilename, intro, avLogDens);
	}
	{ // average log density ratios (KL information))
		const std::string intro("Average log density ratios (KL information)");
		outputRepResults(logFilename, intro, avLogDenRatios);
	}
	{ // estimated L1 error using quasi-random leaves)
		const std::string intro("Approx L1 error by quasi-random sampling");
		outputRepResults(logFilename, intro, estL1ErrorsQR);
	}
	
	
}

void outputResults(const std::string& logFilename,
			const std::vector < std::vector < double > >& timings,
			const std::vector < std::vector < real > >& avLogDens,
			const std::vector < std::vector < real > >& avLogDenRatios,
			const std::vector < std::vector < real > >& estL1ErrorsQR,
			const std::vector < std::vector < size_t > >& leaves)
{
	outputFileStart(logFilename);
	

	{ // timings for 
		const std::string intro("Timings");
		outputRepResults(logFilename, intro, timings);
	}
	{ // leaves)
		const std::string intro("Leaves");
		outputRepResults(logFilename, intro, leaves);
	}
	{ // average log densities
		const std::string intro("Average log densities");
		outputRepResults(logFilename, intro, avLogDens);
	}
	{ // average log density ratios (KL information))
		const std::string intro("Average log density ratios (KL information)");
		outputRepResults(logFilename, intro, avLogDenRatios);
	}
	{ // estimated L1 error using quasi-random leaves)
		const std::string intro("Approx L1 error by quasi-random sampling");
		outputRepResults(logFilename, intro, estL1ErrorsQR);
	}
	
	
}


/* output real results for reps*/
void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < real > >& repResults)
					
{
	std::vector < std::vector < double > > tmp(repResults.size());
	size_t m = repResults.size(); // reps
	for (size_t i = 0; i < m; ++i) {
		size_t d = repResults[i].size();
		tmp[i] = std::vector < double >(d);
		for (size_t j = 0; j < d; ++j) tmp[i][j]=_double(repResults[i][j]);
	}
	
	outputRepResults(filename, intro, tmp);
}

/* output double results for reps*/
void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < double > >& repResults)
					
{
	ofstream os(filename.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << "\n\n" << intro << endl;
		size_t m = repResults.size(); // reps
		for (size_t i = 0; i < m; ++i) {
			os << (i+1) << "\t";
			ostream_iterator<double> out_it (os,"\t");
			copy ( repResults[i].begin(), repResults[i].end(), out_it );
			
			os << endl;
		}
				
		os.close();
		
		if(m > 1) outputRepResultsMeans(filename, intro, repResults);
	}
	else {
		std::cerr << "Error: could not open file named "
			<< filename << std::endl << std::endl;
	}
}

/* output size_t results for reps*/
void outputRepResults(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < size_t > >& repResults)
					
{
	std::vector < std::vector < double > > tmp(repResults.size());
	
	ofstream os(filename.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << "\n\n" << intro << endl;
		size_t m = repResults.size(); // reps
		for (size_t i = 0; i < m; ++i) {
			
			size_t d = repResults[i].size();
			tmp[i] = std::vector < double >(d);
			os << (i+1) << "\t";
			for (size_t j = 0; j < d; ++j) {
				
				tmp[i][j]=1.0*repResults[i][j];
				os << repResults[i][j] << "\t";
			}
			os << endl;
		}
				
		os.close();
		
		if(m > 1) outputRepResultsMeans(filename, intro, tmp);
	}
	else {
		std::cerr << "Error: could not open file named "
			<< filename << std::endl << std::endl;
	}
}

/* output column means for reps*/
void outputRepResultsMeans(const std::string& filename,
					const std::string& intro,
					const std::vector < std::vector < double > >& repResults)
					
{
	size_t m = repResults.size(); // reps
	
	if (m) {
		std::vector < double > means = repResults[0];
		size_t dim = means.size();
		for (size_t i = 1; i < m; ++i) {
			for (size_t j = 0; j < dim; ++j) means[j]+=repResults[i][j];
		}
		for (size_t j = 0; j < dim; ++j) means[j]/=(1.0*m);
	
		ofstream os(filename.c_str(), ios::app);         // append
		if (os.is_open()) {
			os << "Mean";
			for (size_t j = 0; j < dim; ++j) os << "\t" << means[j];
			os << endl;		
			os.close();
		}
		else {
			std::cerr << "Error: could not open file named "
				<< filename << std::endl << std::endl;
		}
	}
}

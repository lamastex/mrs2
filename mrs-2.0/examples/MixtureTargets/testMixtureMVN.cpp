/*
* Copyright (C) 2012 Jennifer Harlow
*
*/


/*! \file
\brief test mixture mvn

 */


#include "mixture_mvn.hpp"

#include "cxsc.hpp"

#include <iostream>
//#include <string>
#include <vector>
#include <cstddef> // size_t
//#include <sstream>
//#include <fstream>
#include <stdexcept>
#include <cassert>

#include <gsl/gsl_randist.h>



#if(0)
std::vector< std::vector < double > >& fileToVector(
	std::vector< std::vector < double > >& vec,
	const std::string& inputFilename);

std::vector< std::vector < cxsc::real > >& fileToVectorReal(
	std::vector< std::vector < cxsc::real > >& vec,
	const std::string& inputFilename);
#endif

using namespace std;

void testUnivariateNormal1();
void testUnivariateNormal2();

void testBivariateNormal1();
void testBivariateNormal2();

void test(const subpavings::kde::MixtureMVN& mixMVN,
	const vector < double >& mean,
		const vector < vector < double > >& scale);

int main()
{
	
	//testUnivariateNormal1();
	
	testUnivariateNormal2();
	
	//testBivariateNormal1();
	
	testBivariateNormal2();
	
	
	return 0; 
}  

void testUnivariateNormal1()
{
	cout << "\n\nTestUnivariate Standard Normal " << endl;
	size_t dim = 1;
	
	vector < vector < double > > means;
	
	vector < double > mean1(dim, 0.0);
	means.push_back(mean1);
	
	vector < vector < vector < double > > > scales;
	
	vector < vector < double > > scale1;
	vector < double > scale1row1(1, 1.0);
	scale1.push_back(scale1row1);
	scales.push_back(scale1);
	
	vector < double > mixes(1, 1.0);
	unsigned long int seed = 1234;
	
	subpavings::kde::MixtureMVN mixMVN(means, scales, mixes, seed);
	
	test(mixMVN, mean1, scale1);
	
}  

void testBivariateNormal1()
{
	cout << "\n\nTestBiivariate Standard Normal " << endl;
	size_t dim = 2;
	
	vector < vector < double > > means;
	
	vector < double > mean1(1, 0.0);
	mean1.push_back(0.0);
	means.push_back(mean1);
	
	vector < vector < vector < double > > > scales;
	
	vector < vector < double > > scale1;
	vector < double > scale1row1(1, 1.0);
	double covar = 0.0;
	scale1row1.push_back(covar);
	vector < double > scale1row2(1, covar);
	scale1row2.push_back(1.0);
	scale1.push_back(scale1row1);
	scale1.push_back(scale1row2);
	scales.push_back(scale1);
	
	vector < double > mixes(1, 1.0);
	unsigned long int seed = 1234;
	
	subpavings::kde::MixtureMVN mixMVN(means, scales, mixes, seed);
	
	test(mixMVN, mean1, scale1);
	
}  

void testUnivariateNormal2()
{
	cout << "\n\nTestUnivariate Non-Standard Normal " << endl;
	size_t dim = 1;
	
	vector < vector < double > > means;
	
	vector < double > mean1(dim, 3.5);
	means.push_back(mean1);
	
	vector < vector < vector < double > > > scales;
	
	vector < vector < double > > scale1;
	vector < double > scale1row1(1, 2.0);
	scale1.push_back(scale1row1);
	scales.push_back(scale1);
	
	vector < double > mixes(1, 1.0);
	unsigned long int seed = 1234;
	
	subpavings::kde::MixtureMVN mixMVN(means, scales, mixes, seed);
	
	test(mixMVN, mean1, scale1);
	
}  

void testBivariateNormal2()
{
	cout << "\n\nTestBiivariate Non-Standard Normal " << endl;
	size_t dim = 2;
	
	vector < vector < double > > means;
	
	vector < double > mean1(1, -3.5);
	mean1.push_back(2.5);
	means.push_back(mean1);
	
	vector < vector < vector < double > > > scales;
	
	vector < vector < double > > scale1;
	vector < double > scale1row1(1, 2.0);
	double covar = 0.6;
	scale1row1.push_back(covar);
	vector < double > scale1row2(1, covar);
	scale1row2.push_back(0.5);
	scale1.push_back(scale1row1);
	scale1.push_back(scale1row2);
	scales.push_back(scale1);
	
	vector < double > mixes(1, 1.0);
	unsigned long int seed = 1234;
	
	subpavings::kde::MixtureMVN mixMVN(means, scales, mixes, seed);
	
	test(mixMVN, mean1, scale1);
	
}  


void test(const subpavings::kde::MixtureMVN& mixMVN,
	const vector < double >& mean,
		const vector < vector < double > >& scale)
{
	size_t dim = scale.size();
	{
		
		vector < real > test(mean.begin(), mean.end());
		real den = mixMVN.f(test);
		cout << "\nDensity at mean (";
		for (size_t j = 0; j < dim; ++j) cout << " " << _double(test[j]);
		cout << ") = " << _double(den) << endl;
		if (dim == 1) {
			double check = gsl_ran_gaussian_pdf(_double(test[0] - mean[0]),
													std::sqrt(scale[0][0]));
			cout << "GSL N(mu, sigma) density at mean " << check << endl;
		}
		if (dim == 2) {
			double check = gsl_ran_bivariate_gaussian_pdf(_double(test[0]  - mean[0]), 
						_double(test[1] - mean[1]),
						std::sqrt(scale[0][0]), 
						std::sqrt(scale[1][1]), 
						(scale[0][1]/(scale[0][0] * scale[1][1])));
			cout << "GSL BiN(0, sigma) density at this point = " << check << endl;
		}
	}
	
	{
		size_t testN = 5;
		cout << "\nGenerate " << testN << " random values:" << endl;
		std::vector < std::vector < real > > rvs;
		mixMVN.prn(rvs, testN);
	
		for (size_t i = 0; i < testN; ++i) {
			cout << "Value " << (i+1) << " is:" << endl;
			for (size_t j = 0; j < dim; ++j) {
				cout << "\t" <<  _double(rvs[i][j]);
				
			}
			cout << endl;
		
			real den = mixMVN.f(rvs[i]);
			cout << "Density at this point = " << _double(den) << endl;
			if (dim == 1) {
				double check = gsl_ran_gaussian_pdf(_double(rvs[i][0] - mean[0] ), 
						std::sqrt(scale[0][0]));
				cout << "GSL N(0, sigma) density at this point  = " << check << endl;
			}
			if (dim == 2) {
				double check = gsl_ran_bivariate_gaussian_pdf(_double(rvs[i][0]  - mean[0]), 
							_double(rvs[i][1] - mean[1]),
							std::sqrt(scale[0][0]), 
							std::sqrt(scale[1][1]), 
							(scale[0][1]/(scale[0][0] * scale[1][1])));
				cout << "GSL BiN(0, sigma) density at this point = " << check << endl;
			}
		}
			
		
	}
	
}  


#if(0)
std::vector< std::vector < double > >& fileToVector(
	std::vector< std::vector < double > >& vec,
	const std::string& inputFilename)
{
    std::ifstream input (inputFilename.c_str());
    std::string lineData;
	
	std::vector< std::vector < double > > result;
	result.reserve(1000000);

    if (input.is_open()) {
		while (getline(input, lineData)) {
			double d;
			std::vector< double > row;
			row.reserve(10);
			std::stringstream lineStream(lineData);

			while (lineStream >> d)
				row.push_back(d);

			result.push_back(row);
		}
	}
	else { // dataFile not open
		std::cerr << "Unable to open file " << inputFilename << std::endl;
		throw std::invalid_argument("File name invalid");
	}

    vec.swap(result);
	return vec;
}

std::vector< std::vector < cxsc::real > >& fileToVectorReal(
	std::vector< std::vector < cxsc::real > >& vec,
	const std::string& inputFilename)
{
    std::ifstream input (inputFilename.c_str());
    std::string lineData;
	
	std::vector< std::vector < cxsc::real > > result;
	result.reserve(1000000);

    if (input.is_open()) {
		while (getline(input, lineData)) {
			double d;
			std::vector< cxsc::real > row;
			row.reserve(10);
			std::stringstream lineStream(lineData);

			while (lineStream >> d)
				row.push_back(d);

			result.push_back(row);
		}
	}
	else { // dataFile not open
		std::cerr << "Unable to open file " << inputFilename << std::endl;
		throw std::invalid_argument("File name invalid");
	}

    vec.swap(result);
	return vec;
}
#endif

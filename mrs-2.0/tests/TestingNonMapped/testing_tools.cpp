/*
* Copyright (C) 2011 Jennifer Harlow
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


/*! \file
\brief Testing tools
 */

#include "testing_tools.hpp"
#include "toolz.hpp"

#include <gsl/gsl_statistics.h>

#include <iterator>
#include <fstream>  // input and output streams
#include <cfloat> // for DBL_EPSILON


using namespace std;
using namespace subpavings;



bool checkFileLines(const std::string& s, std::size_t expectedLines)
{
	bool retValue = false;
	
	std::ifstream dataFile(s.c_str());
	
	if (dataFile.is_open())
	{
		std::size_t readNonBlankLines = 0;
		std::string line;
	
		while (dataFile.good() )
		{
			getline (dataFile,line);
			if (!line.empty()) readNonBlankLines++;
		}
		retValue = (readNonBlankLines == expectedLines);	
		
	}
	else { // dataFile not open
		std::cerr << "Unable to open file " << s << std::endl;
	}

	return retValue;
}

RVecData combineData(const RVecData& v1, const RVecData& v2)
{
	RVecData retValue = v1;
	retValue.insert(retValue.end(), v2.begin(), v2.end() );
	
	return retValue;
}

RVecData combineData(const RVecData& v1, const RVecData& v2, const RVecData& v3)
{
	RVecData retValue = v1;
	retValue.insert(retValue.end(), v2.begin(), v2.end() );
	retValue.insert(retValue.end(), v3.begin(), v3.end() );
	
	return retValue;
}


void outputADH(const std::string& s, const AdaptiveHistogram& adh)
{
	adh.outputToTxtTabs(s, 5, true);
}

void outputSPS(const std::string& s, const SPSnode& spn, int prec)
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		spn.leavesOutputTabs(os, prec);
		os.close();
		cout << s << " output to file" << endl;
	}
	else cout << "Error opening file " << s << endl;
}

void swapCheckOutput(const std::string& s, const AdaptiveHistogram& adh) 
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << adh.stringSummary();
	
		os.close();
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


void swapCheckOutput(const std::string& s, const SPSnode& spn, bool append)
{
	ofstream os;
	if (!append) os.open(s.c_str()); // don't append
	else os.open(s.c_str(), ios_base::app); //append
	if (os.is_open()) {

		os << spn.nodeStringSummary() << endl;
		
		os << endl;
		os.close();
		
		// recurse
		if ( spn.getLeftChild() ) swapCheckOutput(s, *spn.getLeftChild(), true);
		if ( spn.getRightChild() ) swapCheckOutput(s, *spn.getRightChild(), true);
		
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

void checkOutput(const std::string& s, const AdaptiveHistogram& adh) 
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << adh.stringSummary() << endl;
		
		os << "SPS is " << endl;
			
		os.close();
		swapCheckOutput(s, *(adh.getSubPaving()), true); //append
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


void doubleCheckOutput(const std::string& s, const SPSnode& spn, bool append)
{
	ofstream os;
	if (!append) os.open(s.c_str()); // don't append
	else os.open(s.c_str(), ios_base::app); //append
	if (os.is_open()) {

		os << spn.doubleCheckStringSummary() << endl;
		
		os.close();
		
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

void doubleCheckOutput(const std::string& s, const AdaptiveHistogram& adh) 
{
	ofstream os;
	os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << adh.doubleCheckStringSummary() << endl;
		
		os.close();
	}
	else {
		std::cout << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

void cov_calculate_gsl(gsl_matrix *r, gsl_matrix *m)
{
	gsl_vector_view a, b;
	size_t i, j;
	
	for (i = 0; i < m->size2; i++) {
		a = gsl_matrix_column (m, i);
		for (j = 0; j < m->size2; j++) {
			
			double v;
			b = gsl_matrix_column (m, j);
			//gsl_stats_cov does the sample var(covar), ie the unbiased estimator (div by n-1), not the mle
			v = gsl_stats_covariance (a.vector.data, a.vector.stride,
				b.vector.data, b.vector.stride, a.vector.size);
			gsl_matrix_set (r, i, j, v);
		}
	}
}


RealVec checkVarCov(const RVecData& data) 
{
	// each rvector goes into a row
	// so that the number of columns is the data dimensions
	// number of rows is the number of rvectors
	
	// dress the data up as a gsl matrix
	int rows = 0;
	int cols = 0;
	if (!data.empty()) {
		rows = data.size();
		cols = VecLen(data.at(0));
	}
	
	RealVec retValue;
	
	if (rows > 1 && cols) {
		
		gsl_matrix * m = NULL;
		gsl_matrix * r = NULL;
		m = gsl_matrix_alloc (rows, cols);
		r = gsl_matrix_alloc (cols, cols);
		
		for (int i = 0; i < rows; ++i) {
			cxsc::rvector thisrv = data.at(i);
			int lb = Lb(thisrv);
			for (int j = 0; j < cols; ++j) {
				gsl_matrix_set ( m, i, j, _double(thisrv[lb+j]) );
			} 
		}
			
		if (rows && cols) cov_calculate_gsl(r, m);
		
		// now get the covariances out again
		for (int i = 0; i < cols; ++i) {
			for (int j = 0; j < cols; ++j) {
				retValue.push_back( cxsc::real(gsl_matrix_get(r, i,j)) );
			} 
		}

		gsl_matrix_free (m);
		gsl_matrix_free (r);
	}
	if (rows < 2 && cols) {
		
		// now get the covariances out again
		for (int i = 0; i < cols*cols; ++i) {
			retValue.push_back( SignalingNaN );
		}
	}

	return retValue;
	
}

cxsc::rvector checkMean(const RVecData& data) 
{
	// each rvector goes into a row
	// so that the number of columns is the data dimensions
	// number of rows is the number of rvectors
	
	// dress the data up as a gsl matrix
	int rows = 0;
	int cols = 0;
	if (!data.empty()) {
		rows = data.size();
		cols = VecLen(data.at(0));
	}
		
	cxsc::rvector retValue(cols);
	
	if (rows) {
	
		for (int j = 0; j < cols; ++j) {
			retValue[1+j] = 0.0;
		} 
		
		for (int i = 0; i < rows; ++i) {
			retValue += data.at(i);
		}
		
		retValue /= (rows * 1.0);
	}
	else {
		for (int j = 0; j < cols; ++j) {
			retValue[1+j] = SignalingNaN;
		}
	}
	
	return retValue;
}

bool checkAllNaN(const RealVec& vcov) 
{
	int n = vcov.size();
	
	bool retValue = true;
	
	int i = 0;
	
	while (retValue && i < n) {
		retValue = IsSignalingNaN( vcov.at(i) );
		i++;
	}
	return retValue;
}

bool checkAllNaN(const cxsc::rvector& rv) 
{
	int d = VecLen(rv);
	int lb = Lb(rv);
	bool retValue = true;
	
	int i = 0;
	
	while (retValue && i < d) {
		retValue = IsSignalingNaN( rv[lb+i] );
		i++;
	}
	return retValue;
}

bool checkSame(const cxsc::rvector& rv1, const cxsc::rvector& rv2, int n) 
{
	int d = VecLen(rv1);
	if (VecLen(rv2) != d) {
		return false;
	}
	
	int lb1 = Lb(rv1);
	int lb2 = Lb(rv2);
	
	bool retValue = true;
	
	int i = 0;
	
	//cout << cxsc::SaveOpt << cxsc::SetPrecision(23,16);
	while (retValue && i < d) {
		
		//cout << "pair is " << rv1[lb1+i] << "\t" << rv2[lb2+i] << endl;
		//cout <<  "cxsc::abs(rv1[lb1+i] - rv2[lb2+i]) = " << cxsc::abs(rv1[lb1+i] - rv2[lb2+i]) << endl;
		//cout <<  "n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(rv1[lb1+i]), cxsc::abs(rv2[lb2+i]))) = " << n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(rv1[lb1+i]), cxsc::abs(rv2[lb2+i]))) << endl;
			
		bool isNotEqual = cxsc::abs(rv1[lb1+i] - rv2[lb2+i]) 
			> n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(rv1[lb1+i]), cxsc::abs(rv2[lb2+i])));
		//cout << "isNotEqual is " << isNotEqual;	
		retValue = !isNotEqual;
		//cout << "retValue is " << retValue << endl;
		i++;
	}
	//cout << cxsc::RestoreOpt;
	return retValue;
}

bool checkSame(cxsc::real r1, const cxsc::real r2, int n) 
{
	
	return !(cxsc::abs(r1 - r2) 
			> n* DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(r1), cxsc::abs(r2))));
	
}

bool checkSame(const subpavings::RealVec& vec1, const subpavings::RealVec& vec2, int n) 
{
	int d = vec1.size();
	if (vec2.size() != d) {
		return false;
	}
	
	bool retValue = true;
	
	int i = 0;
	
	cout << cxsc::SaveOpt << cxsc::SetPrecision(23,16);
	while (retValue && i < d) {
		//cout << "pair is " << vec1.at(i) << "\t" << vec2.at(i) << endl;
		//cout <<  "cxsc::abs(vec1.at(i) - vec2.at(i)) = " << cxsc::abs(vec1.at(i) - vec2.at(i)) << endl;
		//cout <<  "n * n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(vec1.at(i)), cxsc::abs(vec2.at(i)))) = " 
				//<< n * n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(vec1.at(i)), cxsc::abs(vec2.at(i)))) << endl;
			
		bool isNotEqual = cxsc::abs(vec1.at(i) - vec2.at(i)) 
			> n * n * DBL_EPSILON * cxsc::max(1.0, cxsc::max(cxsc::abs(vec1.at(i)), cxsc::abs(vec2.at(i))));
		//cout << "isNotEqual is " << isNotEqual << endl;	
		
		retValue = !isNotEqual;
		i++;
	}
	cout << cxsc::RestoreOpt;	
	return retValue;
}

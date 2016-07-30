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

/*!/ \file
\brief Testing input functions
*/


#include "sptools.hpp"
#include "toolz.hpp"

// include fstream so as to be able to output a file
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>
#include <iterator>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace subpavings;


void theProblem();
void types();
void testCountNumbers();
void testCheckString();
void testReadRvectorsFromTxt();
void testGetRvectorsFromVecDbl();
void testReadRvectorsFromRVec();
void readFileToCout(const std::string& s);
void testReadVecRealVecFromTxt();
void testReadRvectorsFromTxtReqDims();
void testReqDimsOrd(RVecData data, const string& s, size_t expectedDataLines, 
				vector <int>& reqDims, const string& desc, const string& should, int headers = -1);
void testReqDimsParanoid(RVecData data, const string& s, size_t expectedDataLines, 
				vector <int>& reqDims, const string& desc, const string& should, int headers);
void testReqDims(int checkLevel, RVecData data, const string& s, size_t expectedDataLines, 
		vector <int>& reqDims, const string& desc, const string& should, int headers);

int main()
{
	//theProblem();
	//types();
	//testCountNumbers();
	//testCheckString();
	testReadRvectorsFromTxt();
	//testGetRvectorsFromVecDbl();
	//testReadRvectorsFromRVec();
	//testReadVecRealVecFromTxt();
	//testReadRvectorsFromTxtReqDims();
	
	return 0;
}

void testReqDimsOrd(RVecData data, const string& s, size_t expectedDataLines, 
vector <int>& reqDims, const string& desc, const string& should, int headers)
{
	testReqDims(1, data, s, expectedDataLines, reqDims, desc, should, headers);

}

void testReqDimsParanoid(RVecData data, const string& s, size_t expectedDataLines, 
vector <int>& reqDims, const string& desc, const string& should, int headers)
{
	testReqDims(2, data, s, expectedDataLines, reqDims, desc, should, headers);

}

void testReqDims(int checkLevel, RVecData data, const string& s, size_t expectedDataLines, 
vector <int>& reqDims, const string& desc, const string& should, int headers)
{
	if (checkLevel >= 2 && headers < 0) throw invalid_argument("Paranoid checking needs headers specified");

	size_t existingSize = data.size();

	if (checkLevel == 1) cout << "Ordinary level checking" << endl;
	if (checkLevel == 2) cout << "Paranoid level checking" << endl;
	cout << desc << ", " << " filename " << s <<endl;
	if (headers >= 0) cout << headers << " headers specified, ";
	else cout << "no headers specified, ";
	cout << "reqDims is " << reqDims << endl;
	cout << "expect " << expectedDataLines << " lines" << endl;
	cout << "\n" << should << "\n" << endl;
	
	if (data.empty()) {
		cout << "Before operation data vector is empty" << endl;
	}
	else {
		cout << "Before operation, data vector is:" << endl;
		for (RVecDataItr it = data.begin(); it < data.end(); ++it) {
			prettyPrint(cout, (*it));
			cout << endl;
		}
	}
	cout << endl;
	
	cout << "File to be read in is" << endl;
	readFileToCout(s);
	
	bool result = false;		
	if (checkLevel == 1 && headers >= 0) result = readRvectorsFromTxtOrd(data, s, reqDims, headers); 
	if (checkLevel == 1 && headers < 0) result = readRvectorsFromTxtOrd(data, s, reqDims); 
	if (checkLevel == 2) result = readRvectorsFromTxtParanoid(data, s, reqDims, headers); 
	
	
	if (result) {
		cout << "\nOperation returned true: after operation, data vector is:" << endl;
		for (RVecDataItr it = data.begin(); it < data.end(); ++it) {
			prettyPrint(cout, (*it));
			cout << endl;
		}
		cout << endl;

		if (!data.empty()) assert(VecLen(data.back()) == reqDims.size());
	}
	else {
		cout << "Operation returned false" << endl;
		
	}
	
	assert(data.size() == existingSize + expectedDataLines);
	cout << "(end test )\n" << endl;
	//cout << "press enter\n" << endl;
	//if (checkLevel ==2) getchar();
}

void testReadRvectorsFromTxtReqDims()
// Ord reading 
{
	
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		vector <int> reqDims;
		
		string desc("Test basic file of 1-d data, empty reqDims");
		string should("should throw illegal argument exception");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {0, 4, 10}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 1-d data, reqDims with negative numbers");
		string should("should throw illegal argument exception");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 1-d data, reqDims with all dims");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 1-d data, reqDims out of range");
		string should("should abort");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 1-d data, ");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 3;
		
		int int_array[] = {5, 2, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 5-d data with DOS line endings, ");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 0;
		
		int int_array[] = {5, 0, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 5-d data with DOS line endings, ");
		string should("should throw an invalid argument exception");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 0;
		
		int int_array[] = {5, -2, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 5-d data with DOS line endings, ");
		string should("should throw an invalid argument exception");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 0;
		
		int int_array[] = {6, 2, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 5-d data with DOS line endings, ");
		string should("should abort because specified dim out of range");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {1, 2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 1-d data, ");
		string should("should abort because number of dims required > data dims");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test basic file of 2-d data, but reqDims only specifies one of these");
		string should("should read in second column of given data");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_ints.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data as ints with extra blank line, reqDims reverses order");
		string should("should read okay with column order reversed");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 3;
		
		int int_array[] = {1, 2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with no lines following headers");
		string should("should abort because data dims cannot be found");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 1;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);

	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, 3 lines of headers, first line just numbers");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers, first line just numbers");
		string should("should abort because first line is regarded as data and dims don't match");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 11;
		size_t headers = 0;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers");
		string should("will read in 1-d data, including part of the date line");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers");
		string should("should abort because it uses header line and thinks data is only 1-d");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {1,2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, actually 3 lines of headers");
		string should("should abort because first data line found does not match expected dims from reqDims");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;

		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should read in first elements on each line as 1-d data");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should abort because data dims from first line taken as 1-d");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should abort because data dims from first line taken as 1-d");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should read in first elements on each line");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should read in second column, replacing missing elements with 0.0");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements in some lines");
		string should("should read in, replacing missing elements in second column of txt file with 0.0");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 7;
		size_t headers = 0;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements and illegals in some lines");
		string should("should spit out 3 data lines with illegal formatting, and read in first element on each line");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);

	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 7;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements and illegals in some lines");
		string should("should spit out 3 data lines with illegal formatting but replace missing values with 0.0");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);

	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 7;
		size_t headers = 0;
		
		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data, missing elements and illegals in some lines");
		string should("should spit out 3 data lines with illegal formatting but replace missing values in second column of txt file with 0.0");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);

	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile3d1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;

		int int_array[] = {3, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with comma delimiters");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		rvector r(2); // 2-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		data.push_back(r); // add 2-d vector to the data
		size_t existingSize = data.size();
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;

		int int_array[] = {3, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with comma delimiters, reqDims only wanting 2 of these dimensions");
		string should("should be okay");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		rvector r(2); // 2-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		data.push_back(r); // add 2-d vector to the data
		size_t existingSize = data.size();
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with comma delimiters, reqDims only wanting one of these");
		string should("should abort because number of required dimensions does not match existing data");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		rvector r(2); // 2-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		data.push_back(r); // add 2-d vector to the data
		size_t existingSize = data.size();
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {3, 2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with comma delimiters, reqDims wanting all of these in a different order");
		string should("should abort because required number of dimensions does not match existing data");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	
	
	try {
		RVecData data;
		rvector r(4); // 4-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		r[3] = 3.3;
		r[4] = 4.4;
		data.push_back(r); // add 4-d vector to the data
		size_t existingSize = data.size();
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {3, 2, 1, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with comma delimiters, reqDims size 4 by repeating dim 1");
		string should("should abort because fewer columns in data than seems to be required (ie ignores repeats)");
		testReqDimsOrd(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	// Paranoid reading
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with first header line just numbers");
		string should("should be okay");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);

	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with first header line just numbers");
		string should("should abort because first 'data' line does not match expected dimensions");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 11;
		size_t headers = 0;
		
		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with first header line just numbers");
		string should("will read in first column of data including a line for the assumed data in the first header ");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with first header line just numbers");
		string should("should abort because required dim is larger than dim of data assumed from first 'data' line read");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		int int_array[] = {0}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should throw illegal argument exception");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
		throw logic_error("Should not be able to do this");
		
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;

		int int_array[] = {1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should read in first value found in each line of txt file");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with first header line just numbers");
		string should("abort since data does not seem to match expected dimensions");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 6;
		size_t headers = 0;

		int int_array[] = {2, 1}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should reject 4 lines");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile3d_missingelements1.txt");
		size_t expectedDataLines = 2;
		size_t headers = 0;
		int dim = 3;

		int int_array[] = {2, 1, 3}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 3-d data with missing elements");
		string should("should reject 8 lines with missing elements");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 6;
		size_t headers = 3;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should reject 4 lines");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 6;
		size_t headers = 0;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should spit out headers and 4 data lines");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 4;
		size_t headers = 0;

		int int_array[] = {2}; 
		std::vector <int> reqDims (int_array, int_array + sizeof(int_array) / sizeof(int) );
		
		string desc("Test file of 2-d data with missing elements");
		string should("should reject headers and 6 data lines");
		testReqDimsParanoid(data, s, expectedDataLines, reqDims, desc, should, headers);
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
}

void testReadVecRealVecFromTxt()

{
	
	try {
		vector < RealVec > data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = -1;
		
		cout << "Test basic file of 1-d data, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould throw illegal argument exception" << endl;
		bool result = readVecRealVecFromTxt(data, s, headers, dim); 
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test basic file of 1-d data, " << headers << " headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		cout << "File is" << endl;
		
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		int dim = 1;
		
		cout << "Test basic file of 1-d data, no headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, 0); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 3;
		size_t headers = 0;
		int dim = 5;
		
		cout << "Test basic file of 5-d data with DOS line endings, no headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test basic file of 1-d data, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dim); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		vector < RealVec > data;
		string s("inputFile2d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test basic file of 2-d data with, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nreads in 1-d data" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dim); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test basic file of 2-d data, " << headers << " headers specified, no dim specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_ints.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim =2;
		
		cout << "Test basic file of 2-d data, " << headers << " headers specified, extra blank line at end, data as ints, no dims specified ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers ); //const std::size_t headerlines)
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 3;
		
		cout << "Test file with 2 column headings containing , " << headers << " headers, no lines following headers, no dims specified";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be unable to find any data" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dim = 2;
		
		cout << "Test basic file of 2-d data, " << headers << " headers specified, with data lines following headers, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 1;
		int dim = 2;
		
		cout << "Test basic file of 2-d data, " << headers << " headers specified, but actual headers = 3, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay because  first non-blank line read as data is the header which is ignored due to illegal chars" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		int dim = 2;
		
		cout << "Test basic file of 2-d data, no headers specified, but actual headers = 3, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, 0); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, and actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dims); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dims);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould abort because dims found from first line (1) < specified dims" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dims); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 30;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, but actual headers = 3, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << "\nwill import messed up data because of misleading first line read" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 15;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test malformed file of 2-d data, " << headers << " headers, missing elements in some lines, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill read in all data as 1-d" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		vector < RealVec > data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test malformed file of 2-d data, " << headers << " headers, missing elements in some lines, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill abort read because dims found from data (1) < dims specified" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dims); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 8;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data, " << headers << " headers, missing elements in some lines, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill read data incorrectly, muddling elements on different lines" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
			cout << (*it) << endl;
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		vector < RealVec > data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data, " << headers << " headers, missing elements & illegals in some lines, dims not specified";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill read in data incorrectly because of illegal characters" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		string s("inputFile3d1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 3;
		
		cout << "Test file of 3-d data with comma delimiters, " << headers << " headers, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould be okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		cout << "\n data vector is" << endl;
		
		for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
		{
			cout << (*it) << endl;
		}
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		data.push_back(RealVec());
		data.back().push_back(1.1);
		data.back().push_back(2.2);
		
		size_t existingSize = data.size();
		
		{
			cout << "\n data vector at start is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test file of 3-d data with comma delimiters, " << headers << " headers, no dims specified, existing data does not match new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill only read 2 values from each line and discard rest" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		{
			cout << "\n data vector at end is" << endl;
		
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it)
			{
				cout << (*it) << endl;
			}
		
		cout << "\n" << endl;
		
		}
		
		assert(data.size() == expectedDataLines + existingSize);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		data.push_back(RealVec());
		data.back().push_back(1.1);
		data.back().push_back(2.2);
		data.back().push_back(3.3);
		data.back().push_back(4.4);
		size_t existingSize = data.size();
		
		{
			
			cout << "\n data vector at start is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		cout << "Test file of 3-d data with comma delimiters, " << headers << " headers, no dims specified, will expect more 4-d data";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill abort read because existing data dims > dims found from file " << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		{
			cout << "\n data vector at end is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		data.push_back(RealVec());
		data.back().push_back(1.1);
		data.back().push_back(2.2);
		data.back().push_back(3.3);
		size_t existingSize = data.size();
		
		{
			cout << "\n data vector at start is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 5;
		int dim = 3;
		
		cout << "Test file of 3-d data with comma delimiters, " << headers << " headers, no dims specified, existing data matches new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << "\n\nshould b okay" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers); 
		
		{
			cout << "\n data vector at end is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines + existingSize);
		assert(data.back().size() == dim);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		vector < RealVec > data;
		data.push_back(RealVec());
		data.back().push_back(1.1);
		data.back().push_back(2.2);
		data.back().push_back(3.3);
		size_t existingSize = data.size();
		
		{
			cout << "\n data vector at start is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 5;
		int dim = 2;
		
		cout << "Test file of 3-d data with comma delimiters, " << headers << " headers, dims specified as " << dim << ", existing data matches new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould abort because dims specified does not match data" << endl;
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		bool result = readVecRealVecFromTxt(data, s, headers, dim); 
		
		{
			cout << "\n data vector is" << endl;
			for (vector < RealVec >:: iterator it = data.begin(); it < data.end(); ++it) {
				cout << (*it) << endl;
			}
			cout << "\n" << endl;
		
		}
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	cout << "\nEnd all testing\n" << endl;
}


void testReadRvectorsFromRVec()
{
	try {
			
		RVecData data;
		
		RVecData inputData;
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test empty input vector, ";
		cout << expectedDataLines << " lines " << " \n\nshould (silently) fail to read anything\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		RVecData inputData;
		inputData.push_back(cxsc::rvector());
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		
		cout << "Test input vector containing just one empty vector, ";
		cout << expectedDataLines << " lines " << " \n\nshould throw illegal_argument exception\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		throw logic_error("should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		
		RVecData inputData;
		inputData.push_back(cxsc::rvector());
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		
		cout << "Test input vector containing one empty vector and then two 2-d vectors";
		cout << expectedDataLines << " lines " << " \n\nshould throw illegal_argument exception\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		throw logic_error("should not be able to do this");
		
		
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		
		RVecData inputData;
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test input vector of all 2-d data, ";
		cout << expectedDataLines << " lines " << " \n\nshould be okay\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		RVecData inputData;
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test input vector with first 'line' 2-d, rest 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		RVecData inputData;
		
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();;
		size_t existingSize = data.size();
		
		cout << "Test input vector with first line 3-d, second and third 2-d, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	try {
			
		RVecData data;
		
		RVecData inputData;
		
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		size_t expectedDataLines = inputData.size();;
		size_t existingSize = data.size();
		
		cout << "Test input vector with first line 2-d, second 3-d and third 2-d, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		
		RVecData inputData;
		
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 1;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test input vector with first line 3-d, second and third 2-d, fourth 1-d";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(3);
			rv[1] = 1.1;
			rv[2] = 4.4;
			rv[3] = 6.6;
			data.push_back(rv);
		}
		
		RVecData inputData;
		
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		
		cout << "Test 4 lines of existing 3-d data, input vector with first and second line 2-d, third 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nshould abort: existing data dimension not same as first in input \n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(1);
			rv[1] = 1.1;
			data.push_back(rv);
		}
		
		RVecData inputData;
		
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 1;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			cxsc::rvector rv;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test 4 lines of existing 1-d data, input vector with first line 1-d, second line empty, third line 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
	for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(2);
			rv[1] = 1.1;
			rv[2] = 4.4;
			data.push_back(rv);
		}		
		RVecData inputData;
		
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();;
		size_t existingSize = data.size();
		
		cout << "Test 4 lines of existing 2-d data, input vector with first and second line 2-d, third 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

// with checks on dims
try {
			
		RVecData data;
		
		RVecData inputData;
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test input vector of all 2-d data, with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nshould be okay\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		RVecData inputData;
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 2;
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test input vector with first 'line' 2-d, then 3-d, 2-d, 3-d with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in only 2-d data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		RVecData inputData;
		
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			cxsc::rvector rv;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		
		size_t expectedDataLines = 3;
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test input vector with first and second line 3-d, third 0-d,fourth 3-d with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in only 3-d data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	
	
	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(3);
			rv[1] = 1.1;
			rv[2] = 4.4;
			rv[3] = 6.6;
			data.push_back(rv);
		}
		
		RVecData inputData;
		
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test 4 lines of existing 3-d data, input vector with first and second line 2-d, third 3-d, with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nshould abort: existing data dimension not same as first in input \n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims);  
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(1);
			rv[1] = 1.1;
			data.push_back(rv);
		}
		
		RVecData inputData;
		
		int n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 1;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			cxsc::rvector rv;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 1;
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test 4 lines of existing 1-d data, input vector with first line 1-d, second line empty, third line 3-d, with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in just first line data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
	for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(2);
			rv[1] = 1.1;
			rv[2] = 4.4;
			data.push_back(rv);
		}		
		RVecData inputData;
		
		int n = 2;
		for (int i = 0; i < n; ++i) {
			int d = 2;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			int d = 3;
			cxsc::rvector rv(d);
			for (int k = 1; k <= d; ++k) rv[k] = 1.0;
			inputData.push_back(rv);
			
		}
		n = 1;
		for (int i = 0; i < n; ++i) {
			cxsc::rvector rv;
			inputData.push_back(rv);
			
		}
		
		size_t expectedDataLines = 2;
		size_t existingSize = data.size();
		bool checkDims = true;
		
		cout << "Test 4 lines of existing 2-d data, input vector with first and second line 2-d, third 3-d, fourth empty, with checkDims = true, ";
		cout << expectedDataLines << " lines " << " \n\nwill read in all data\n" << endl;
		cout << "input data is" << endl;
		
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( inputData.begin(), inputData.end(), out_it1 );
		cout << "\n" << endl;
		
		size_t result = getRvectorsFromRVec(data, inputData, checkDims); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}


}


void testGetRvectorsFromVecDbl()
{
	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test empty input vector, ";
		cout << expectedDataLines << " lines " << " \n\nshould (silently) fail to read anything\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(result == expectedDataLines);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		
		cout << "Test input vector containing just one empty vector, ";
		cout << expectedDataLines << " lines " << " \n\nshould throw illegal_argument exception\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		throw std::logic_error("should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		cxsc::rvector rv;
		data.push_back(rv);
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		
		cout << "Test existing data with one empty rv in it, ";
		cout << expectedDataLines << " lines " << " \n\nshould throw an invalid argument exception\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		throw std::logic_error("should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		int dim = 0;
		
		cout << "Test dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines " << " \n\nshould abort with an exception\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData, dim); 
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}

	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test input vector of all 2-d data, ";
		cout << expectedDataLines << " lines " << " \n\nshould be okay\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.back().push_back(13);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		inputData.back().push_back(23);
		
		size_t expectedDataLines = inputData.size();
		size_t existingSize = data.size();
		
		cout << "Test input vector with first 'line' 2-d, rest 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nshould read in all data but only 2-d\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		
		size_t expectedDataLines = 1;
		size_t existingSize = data.size();
		
		cout << "Test input vector with first line 3-d, second and third 2-d, ";
		cout << expectedDataLines << " lines " << " \n\nshould read in only 3-d data\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		
		size_t expectedDataLines = 1;
		size_t existingSize = data.size();
		size_t dim = 3;
		
		cout << "Test input vector with first line 2-d, second 3-d and third 2-d, ";
		cout << "dim specified as " << dim << endl;
		cout << expectedDataLines << " lines " << " \n\nshould read in only data matching dim\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData, dim); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		inputData.push_back(VecDbl());
		inputData.back().push_back(31);
		
		size_t expectedDataLines = 3;
		size_t existingSize = data.size();
		size_t dim = 2;
		
		cout << "Test input vector with first line 3-d, second and third 2-d, fourth 1-d";
		cout << "dim specified as " << dim << endl;
		cout << expectedDataLines << " lines " << " \n\nshould read in all as 2-d data and fail to read 1-d\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData, dim); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	

	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(3);
			rv[1] = 1.1;
			rv[2] = 4.4;
			rv[3] = 6.6;
			data.push_back(rv);
		}
		
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		
		size_t expectedDataLines = 1;
		size_t existingSize = data.size();
		
		cout << "Test 4 lines of existing 3-d data, input vector with first and second line 2-d, third 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nshould read in only data 3-d data\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it1 (cout,"\n");
		copy ( data.begin(), data.end(), out_it1 );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
			
		RVecData data;
		for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(1);
			rv[1] = 1.1;
			data.push_back(rv);
		}
		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		
		size_t expectedDataLines = 2;
		size_t existingSize = data.size();
		
		cout << "Test 4 lines of existing 1-d data, input vector with first line 2-d, second line empty, third line 3-d, ";
		cout << expectedDataLines << " lines " << " \n\nshould read in only at least 1-d data\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
			
		RVecData data;
	for (int i = 0; i < 4; ++i) {
			cxsc::rvector rv(3);
			rv[1] = 1.1;
			rv[2] = 4.4;
			rv[3] = 6.6;
			data.push_back(rv);
		}		
		std::vector < VecDbl > inputData;
		inputData.push_back(VecDbl());
		inputData.back().push_back(11);
		inputData.back().push_back(12);
		inputData.push_back(VecDbl());
		inputData.back().push_back(21);
		inputData.back().push_back(22);
		inputData.push_back(VecDbl());
		inputData.back().push_back(1);
		inputData.back().push_back(2);
		inputData.back().push_back(3);
		
		size_t expectedDataLines = 0;
		size_t existingSize = data.size();
		int dim = 2;
		
		cout << "Test 4 lines of existing 3-d data, input vector with first and second line 2-d, third 3-d, ";
		cout << "dim specified as " << dim << endl;
		cout << expectedDataLines << " lines " << " \n\nshould abort because expected dim does not match existing data\n" << endl;
		cout << "input data is" << endl;
		
		{
			for (std::vector < VecDbl >::iterator it = inputData.begin();
				it < inputData.end();
				++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
		}
		cout << "\n" << endl;
		
		bool result = getRvectorsFromVectorOfVecDbl(data, inputData, dim); 
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

}


void testReadRvectorsFromTxt()

{
	
// Ord reading 

	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = -1;
		
		cout << "Test basic file of 1-d data with ordinary input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should throw illegal argument exception***" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dim); 
		
		throw logic_error("Should not be able to do this");
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test basic file of 1-d data with ordinary input checking, " << headers << " headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should be okay***" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		
		cout << "Test basic file of 1-d data with ordinary input checking, no headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should be okay***" << endl;
		readRvectorsFromTxtOrd(data, s); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFileDosCR1.txt");
		size_t expectedDataLines = 3;
		
		cout << "Test basic file of 5-d data with DOS line endings ordinary input checking, no headers specified, dim not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should be okay***" << endl;
		readRvectorsFromTxtOrd(data, s); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test basic file of 1-d data with ordinary input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***give an error that data dims < expected dimensions***" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dim); 
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	
	try {
		RVecData data;
		string s("inputFile1d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test basic file of 1-d data with ordinary input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should be okay***" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 1;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***does not give an error about dimensions: reads in 1-d data***" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers specified, no dim specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should be okay***" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_ints.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers specified, extra blank line at end, data as ints, no dims specified ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtOrd(data, s, headers ); //const std::size_t headerlines)
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 3;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers, no lines following headers, no dims specified";
		cout << expectedDataLines << " lines, filename " << s << " (should be unable to get dimensions)" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers specified, with data lines following headers, ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 1;
		
		cout << "Test basic file of 2-d data with ordinary input checking, " << headers << " headers specified, but actual headers = 3, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers3.txt");
		size_t expectedDataLines = 10;
		
		cout << "Test basic file of 2-d data with ordinary input checking, no headers specified, but actual headers = 3, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtOrd(data, s); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, ordinary input checking, " << headers << " headers specified, and actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (should abort because first line is regarded as data and dims don't match)" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 11;
		size_t headers = 0;
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, but actual headers = 3, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " (will read in 1-d data, including part of the date line)" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 11;
		size_t headers = 0;
		int dims = 1;
		
		
		cout << "Test file of 2-d data with first header line just numbers, " << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (will read in 1-d data including headers)" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test malformed file of 2-d data with ordinary input checking, " << headers << " headers, missing elements in some lines, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will read in all data as 1-d***" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test malformed file of 2-d data with ordinary input checking, " << headers << " headers, missing elements in some lines, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will abort as first data line < expected dims***" << endl;
		readRvectorsFromTxtOrd(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test malformed file of 2-d data with ordinary input checking, " << headers << " headers, missing elements in some lines, dims not specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will read in 2-d data incorrectly***" << endl;
		readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 7;
		size_t headers = 0;
		
		cout << "Test malformed file of 2-d data with ordinary input checking, " << headers << " headers, missing elements & illegals in some lines, dims not specified";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out 3 data lines with illegal formatting - others will be wrong***" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile3d1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test file of 3-d data with comma delimiters, ordinary input checking, " << headers << " headers, no dims specified, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould be okay" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		rvector r(2); // 2-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		data.push_back(r); // add 2-d vector to the data
		size_t existingSize = data.size();
		
		{
			
			cout << "\n data vector before new data (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
			cout << endl;
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		
		cout << "Test file of 3-d data with comma delimiters, ordinary input checking, " << headers << " headers, no dims specified, existing data does not match new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill expect more 2-d data and read all new data as 2-d" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		{
			cout << "\n data vector is (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		rvector r(4); // 4-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		r[3] = 3.3;
		r[4] = 4.4;
		data.push_back(r); // add 4-d vector to the data
		size_t existingSize = data.size();
		
		{
			
			cout << "\n data vector before new data (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
			cout << endl;
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		
		cout << "Test file of 3-d data with comma delimiters, ordinary input checking, " << headers << " headers, no dims specified, will expect more 4-d data";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould abort" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		{
			cout << "\n data vector is (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		rvector r(3); // 3-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		r[3] = 3.3;
		data.push_back(r); // add 3-d vector to the data
		size_t existingSize = data.size();
		
		{
			
			cout << "\n data vector before new data (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
			cout << endl;
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 5;
		
		cout << "Test file of 3-d data with comma delimiters, ordinary input checking, " << headers << " headers, no dims specified, existing data matches new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould b okay" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		{
			cout << "\n data vector is (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		rvector r(3); // 3-d rvector
		r[1] = 1.1;
		r[2] = 2.2;
		r[3] = 3.3;
		data.push_back(r); // add 3-d vector to the data
		size_t existingSize = data.size();
		
		{
			
			cout << "\n data vector before new data (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
			cout << endl;
		}
		
		string s("inputFile3d2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 5;
		int dim = 2;
		
		cout << "Test file of 3-d data with comma delimiters, ordinary input checking, " << headers << " headers, dims specified as " << dim << ", existing data matches new data dimensions";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould abort because dims specified does not match data" << endl;
		bool result = readRvectorsFromTxtOrd(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		{
			cout << "\n data vector is (using cxsc formatting)" << endl;
			ostream_iterator<rvector> out_it (cout,"\n");
			copy ( data.begin(), data.end(), out_it );
		}
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines + existingSize);
		
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
// Paranoid reading
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, paranoid input checking, " << headers << " headers specified, and actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, paranoid input checking" << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (should abort because first 'data' line does not match expected dimensions)" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 1;
		size_t headers = 0;
		int dims = 1;
		
		
		cout << "Test file of 2-d data with first header line just numbers, paranoid input checking, " << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (will only read in data from first 'data' line, ie the numerical header)" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = -1;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should throw illegal argument exception***" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		throw logic_error("Should not be able to do this");
		
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should abort since data does not seem to match expected dimensions***" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 6;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers, missing elements in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out 4 lines***" << endl;
		readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = 3;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers, missing elements in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out all lines***" << endl;
		bool result = readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 6;
		size_t headers = 3;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers, missing elements in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out 4 lines***" << endl;
		bool result = readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 6;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers, missing elements in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out incorrect headers and 4 data lines***" << endl;
		bool result = readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 4;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with paranoid input checking, " << headers << " headers, missing elements & illegals in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should spit out incorrect headers and 6 data lines***" << endl;
		bool result = readRvectorsFromTxtParanoid(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
// fast reading

	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, fast input checking, " << headers << " headers specified, and actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (should be okay)" << endl;
		readRvectorsFromTxtFast(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 12;
		size_t headers = 0;
		int dims = 2;
		
		cout << "Test file of 2-d data with first header line just numbers, fast input checking" << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (will read in any non blank row and make the best of the input that it can)" << endl;
		readRvectorsFromTxtFast(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	try {
		RVecData data;
		string s("inputFile2d_headers4.txt");
		size_t expectedDataLines = 12;
		size_t headers = 0;
		int dims = 1;
		
		
		cout << "Test file of 2-d data with first header line just numbers, fast input checking, " << headers << " headers specified, but actual headers = 3, dims specified as " << dims << ", ";
		cout << expectedDataLines << " lines, filename " << s << " (will read in any non blank row and make the best of the input that it can)" << endl;
		readRvectorsFromTxtFast(data, s, headers, dims); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}




	try {
		RVecData data;
		string s("inputFile2d_blanklines2.txt");
		size_t expectedDataLines = 11;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test file of 2-d data with fast input checking, " << headers << " headers, extra blank lines in file, and line with just delimiter characters ";
		cout << expectedDataLines << " lines, filename " << s << " (will iclude line with just delimiters incorrectly)" << endl;
		readRvectorsFromTxtFast(data, s, headers, dim);
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 0;
		size_t headers = 0;
		int dim = -1;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***should throw illegal argument exception***" << endl;
		readRvectorsFromTxtFast(data, s, headers, dim); 
		
		throw std::logic_error("Should not be able to do this");
		
	}
	catch (std::invalid_argument& ee) {
		cout << "Exception: " << ee.what() << "\n"<< endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will not spit out any lines 5 read in incorrectly***" << endl;
		readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nwill not spit out any lines 4 read in incorrectly" << endl;
		readRvectorsFromTxtFast(data, s, headers, dim); 
		
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}	
	try {
		RVecData data;
		string s("inputFile2d_missingelements2.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 3;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***all lines read in incorrectly (dims padded)***" << endl;
		bool result = readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***read in all data, 4 lines incorrectly***" << endl;
		bool result = readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile2d_missingelements3.txt");
		size_t expectedDataLines = 12;
		size_t headers = 0;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers specified, dim specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will read in whole file including 2 extra lines from non-blank headers, and 4 other incorrect lines***" << endl;
		bool result = readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}

	try {
		RVecData data;
		string s("inputFile2d_missingelements_and_illegalformatting.txt");
		size_t expectedDataLines = 10;
		size_t headers = 3;
		int dim = 2;
		
		cout << "Test malformed file of 2-d data with fast input checking, " << headers << " headers, missing elements & illegals in some lines, ";
		cout << expectedDataLines << " lines, filename " << s << " \n\n***will read in all data including lines with illegal formatting***" << endl;
		bool result = readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	try {
		RVecData data;
		string s("inputFile3d1.txt");
		size_t expectedDataLines = 10;
		size_t headers = 0;
		int dim = 3;
		
		cout << "Test file of 3-d data with comma delimiters, fast input checking, " << headers << " headers, dims specified as " << dim << ", ";
		cout << expectedDataLines << " lines, filename " << s << " \n\nshould be okay" << endl;
		bool result = readRvectorsFromTxtFast(data, s, headers, dim); 
		
		cout << "File is" << endl;
		readFileToCout(s);
		
		cout << "\n data vector is (using cxsc formatting)" << endl;
		ostream_iterator<rvector> out_it (cout,"\n");
		copy ( data.begin(), data.end(), out_it );
		
		cout << "\n" << endl;
		
		assert(data.size() == expectedDataLines);
		cout << "(end test)\n" << endl;
		
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
		throw;
	}
	
	cout << "\nEnd of all tests\n" << endl;
}

void readFileToCout(const std::string& s)
{
	ifstream is;
	is.open (s.c_str() );

	while( is.good() ) {
		std::string line;
		getline(is, line);
		cout << line << endl;
	}
}



void testCheckString()

{
	//string legal("eE+-.0123456789 \t,");
	
	try {
		string stringvalues("");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" ");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("\t");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("eE+-.0123456789 \t");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	
	try {
		string stringvalues("\n");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("123.56e-10");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("-123.56E+10");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("a-123.56e+10");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("-123.56E+10a");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(",");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues("12.0, 3.4");
		
		bool ck = checkString(stringvalues);
		
		cout << "string is " << stringvalues << " checkString returns " << ck << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
}


void testCountNumbers()
{
	try {
		string stringvalues("");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" \t, \n");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" 123");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(".");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" 1.23\t");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" -4.0 +1\t");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" e");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
		try {
		string stringvalues("1e");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
		try {
		string stringvalues("e1");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}

	try {
		string stringvalues("e-04");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" 1.23 1e-04\t");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
	try {
		string stringvalues(" 1a23 1e-04\t");
		
		int cn = countNumbers(stringvalues);
		
		cout << "string is " << stringvalues << " cn is " << cn << endl;
	}
	catch (exception& e) {
		cout << "Exception: " << e.what() << "\n"<< endl;
	}
}

void types()
{
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues(" ");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("a");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("a1");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1a");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("e");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("1e");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("e1");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues(".");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("..");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0 a");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 e");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("a 2.0");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("e 2.0");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1e-3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1E3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0 1e3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0,1e3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0\t1e3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0\n1e3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("\t1e3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues(" 1 2 3,");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("s1s2s3,");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
}


void theProblem()
{
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
		
		string stringvalues("2.2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues(".");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("..");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
				
		string stringvalues("..");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
				
		string stringvalues(". .");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 1;
		cxsc::rvector rv1(dim);
				
		string stringvalues("1.0 2.2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0 2.2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2.2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0 2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 2;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1.0");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2 3");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2 3,");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2 3 ");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2 3.0");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
	{
		int dim = 3;
		cxsc::rvector rv1(dim);
		
		string stringvalues("1 2 3.0,");
		istringstream iss (stringvalues,istringstream::in);
		
		iss >> rv1;
		
		cout << "dim is " << dim << ", string is " << stringvalues << ": rv1 is " << rv1 << endl;
		cout << "good is " << iss.good() << ", eof is " << iss.eof() << ", fail is " << iss.fail() << ", bad is " << iss.bad() << "\n(end test)\n" <<endl;
	}
}




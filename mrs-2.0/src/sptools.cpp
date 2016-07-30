/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
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
\brief Implementation of sptools functions
*/


#include "sptools.hpp"

#include "toolz.hpp"

#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <unistd.h>    

using namespace std;
using namespace subpavings;

// make a unique file name that has a time stamp number appended to it
string subpavings::getUniqueFilename(const std::string& baseFileName, 
				const std::string& suffix)
{
	string s;
	bool newfilename = false;

	while (!newfilename) {
		stringstream out;
		out << time(NULL);
		s = baseFileName + out.str() + suffix;
		//try opening a file with this name for input
		ifstream inp;
		inp.open(s.c_str(), ifstream::in);
		inp.close();
		if(inp.fail()) // not already a file of this name
		{
			inp.clear(ios::failbit);
			newfilename = true;
		}
	} // we should now have a unique filename
	return s;

}

// make a file name that has this process ID (PID) appended to it
string subpavings::getPidFilename(const std::string& baseFileName, 
				const std::string& suffix)
{
	string s;
	
	ostringstream out;
	out << getpid();
	return (baseFileName + out.str() + suffix);
	
}

// Method to add a line to a file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, const std::string& line, bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << endl;
		os << line << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

//UPDATED JUNE 2012 for logposteriors
// Method to add a string to a file
// Output goes to file named according to argument s
// does not add any new lines
void subpavings::outputFileString(const std::string& s, const std::string& line, bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << line << flush;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to put opening line into a log file
void subpavings::outputFileStart(const std::string& s)
{
	// Make a string with filename and timestamp to start log file
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	ofstream os(s.c_str());         // replace data
	if (os.is_open()) {
		os << "File " << s << " created " <<  asctime (timeinfo) << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

//method to strip path from filename
std::string subpavings::stripPath(const std::string& filename)
{
	/* Search for the last '/' in the filename, break
	 * it there*/

	std::string file = filename;

	size_t found;
	found=filename.find_last_of("/");
	if (found!=string::npos) {
		file = filename.substr(found+1);
	}
	return file;
}

// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, RealVec& vals, int i)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << "Pass " << i << std::endl; // numbering
		ostream_iterator<real> out_it (os,"\t");
		copy ( vals.begin(), vals.end(), out_it );
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

//NEWAPRIL2012
// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, RealVec& vals, int i,
		int nodeType, int accepted, size_t nCount)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << i << "\t" << nodeType << "\t" << accepted << "\t" << nCount << "\t" << std::flush;
		ostream_iterator<real> out_it (os,"\t");
		copy ( vals.begin(), vals.end(), out_it );
		os << endl;
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

//NEWAPRIL2012
// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, int i, int val)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << i << "\t" << val << std::endl;
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

void subpavings::outputFile(const std::string& s, int i, size_t val)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << i << "\t" << val << std::endl;
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, const std::string& intro,
				RealVec& vals)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << intro << "\t";
		ostream_iterator<real> out_it (os,"\t");
		copy ( vals.begin(), vals.end(), out_it );
		vals.clear();

		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, const std::string& intro,
				subpavings::IntVec& vals)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << intro << "\t";
		ostream_iterator<int> out_it (os,"\t");
		copy ( vals.begin(), vals.end(), out_it );
		os << "\n";

		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to append values to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, const std::string& intro,
RVecData& vals)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << intro << "\t";
		ostream_iterator<rvector> out_it (os,"\t");
		copy ( vals.begin(), vals.end(), out_it );
		os << "\n";

		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to append strings to output log file
// Output goes to file named according to argument s
void subpavings::outputFile(const std::string& s, vector<std::string>& strings, int i)
{
	ofstream os(s.c_str(), ios::app);         // append
	if (os.is_open()) {
		os << "Pass " << i << "\t" << std::endl; // numbering
		ostream_iterator<string> out_it (os,"\t");
		copy ( strings.begin(), strings.end(), out_it );
		strings.clear();
		os << endl;

		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// parse a tree string to make a .dot file
bool subpavings::parseForGraphDot(const std::string& s, std::string toParse)
{
	bool success = false;
	if (toParse.length() > 1) { // min of "XL" or "XR"

		string sep = " \t,";
		size_t startpos = toParse.find_first_not_of(sep); // start of segment
		size_t endpos = string::npos;
		string segment = "";
		if (startpos != string::npos) { // find start next sep
			endpos = toParse.find_first_of(sep, startpos+1);
			if (endpos != string::npos) {  // get the segment
				segment = toParse.substr(startpos, endpos-startpos);
				toParse = toParse.substr(endpos); // what is left

			}
			else { // must be last segment
				segment = toParse.substr(startpos);
				toParse = "";
			}
		}
		if (segment.length() > 1) {
			string parent = segment.substr(0,segment.length()-1);
			string line = "\t " + parent + " -> " + segment + ";";
			outputFile(s, line);
			success = true;
		}
		if (success && toParse.length() > 0) {
			//recurse
			success = parseForGraphDot(s, toParse);
		}
	} // end if length at start > 1

	return success;
}


// make a Dot graph image given a dot file
void subpavings::makeDotImage(const std::string& s)
{
	string baseFileName = "graph";
	string suffix = ".png";
	string graphName = getUniqueFilename(baseFileName, suffix);

	// make the graph
	string commandLine = "dot -Tpng " + s + " -o " + graphName;
	system(commandLine.c_str());
	std::cout << endl;
	std::cout << "Graph output to " << graphName << std::endl << std::endl;

}

// method to count lines in a txt file
int subpavings::countLinesInTxt(const std::string& s, const size_t headerlines)
{
	//set up the file and read input line by line
	// we need to convert the string argument to a c-string for ifstream
	ifstream dataFile(s.c_str());
	string line;    // a string object to use in counting lines
	int howMany = -headerlines; // how many lines in the file

	if (dataFile.is_open())
	{
		// count the lines in the file
		while (dataFile.good() )
		{
			getline (dataFile,line);
			howMany++;  // count the number of lines in file
		}
		dataFile.close();
	}

	else {
		std::cerr << "Error in countLinesInTxt. "
			<< "Unable to open file" << std::endl;
		
	}

	return howMany;
}

// make a box to fit some datat
ivector subpavings::makeBox(const RVecData& theData, real padding)
{
	if (theData.empty()) throw std::invalid_argument(
		"subpavings::makeBox(const RVecData&, real): theDate.empty()");
	
	int dim = VecLen(theData.back());
	
	return makeBox(theData, dim, padding);
}

// make a box to fit some datat
ivector subpavings::makeBox(const RVecData& theData, int dim, 
							real padding)
{
	if (theData.empty()) throw std::invalid_argument(
		"subpavings::makeBox(const RVecData&, int, real): theData.empty()");
	
	// set up a vector of maxes
    vector<real> maxs;

    // give maxs starting values from the first element in the rvectors
    rvector first = *theData.begin();

    for (int i = 1; i <=dim; ++i) {
        maxs.push_back(first[i]);
    }

    // make mins the same as maxes to start with
    vector<real> mins = maxs;

    RVecDataCItr cit;

    // go over the rest of the container
    for(cit = theData.begin()+1; cit < theData.end(); ++cit) {
        for (int i = 1; i <= dim; ++i) {
            real r = (*cit)[i];
            // vectors indexed 0 - n-1, rvectors ndexed 1 - n
            if(r < mins[i-1]) {
                mins[i-1] = r;
            }
            if(r > maxs[i-1]) {
                maxs[i-1] = r;
            }
        } // end going through rvector elements
    } // end going through rvectors

    ivector retVal(dim);    // set up an ivector to become the return value

    // and make each interval the (min, max) of the corresponding elements
    // of the rvectors -/+ some padding

    // make intervals and make them elements of the ivector
    for (int i = 1; i <=dim; ++i) {
        interval myInterval(mins[i-1]-padding, maxs[i-1]+padding);
        retVal[i]=myInterval;
    }
   
    return retVal;

}

// make a box to fit some datat
ivector subpavings::makeBox(
					const std::vector < std::vector <double> >& theData,
					real padding)
{
	if (theData.empty()) throw std::invalid_argument(
		"subpavings::makeBox(const std::vector < std::vector <double> >&, real): theData.empty()");
	
	size_t dim = theData.back().size();
	// set up a vector of maxes
	// give maxs starting values from the first element in the theData
    vector < real > maxs(theData.begin()->begin(), theData.begin()->end()); 

	// make mins the same as maxes to start with
    vector<real> mins = maxs;

    std::vector < std::vector <double> >::const_iterator cit;

    // go over the rest of the container
    for(cit = theData.begin()+1; cit < theData.end(); ++cit) {
		
        for (size_t i = 0; i < dim; ++i) {
            real r = (*cit)[i];
            if(r < mins[i]) {
                mins[i] = r;
            }
            if(r > maxs[i]) {
                maxs[i] = r;
            }
        } // end going through inner vector elements
    } // end going through inner vectors

    ivector retVal(dim);    // set up an ivector to become the return value

    // and make each interval the (min, max) of the corresponding elements
    // of the rvectors -/+ some padding

    // make intervals and make them elements of the ivector
    for (size_t i = 1; i <=dim; ++i) {
        interval myInterval(mins[i-1]-padding, maxs[i-1]+padding);
        retVal[i]=myInterval;
    }
   
    return retVal;

}


bool subpavings::checkBox(const cxsc::ivector& box)
{
	int low_index = Lb(box);
	int upp_index = Ub(box);

	bool retValue = (upp_index >= low_index);

	if (retValue) {
		int d = upp_index - low_index +1;
		int i = 1;
		while (i <= d && retValue) {
			retValue = !(succ(Inf(box[i])) > Sup(box[i]));
			++i;
		}
		if (retValue) retValue = (realVolume(box) >= cxsc::MinReal);
	}	
	return retValue;
}

// the first set of these functions use only base class attributes

// return TRUE if volume of a < volume of b
bool subpavings::volCompare(const ivector &a, const ivector &b)
{
	bool returnValue = 0;

	// Make sure vectors have same number of elements
	// and at least one element each
	if( (Ub(a)-Lb(a)) == (Ub(b)-Lb(b)) && (Ub(a)-Lb(a))>=0 ) {
		// compare the two volumes
		returnValue = ((Volume(a)<Volume(b)));

	}
	else {
		std::cerr << "Error in volCompare : comparing "
			<< "ivectors of different dimensions"
			<< std::endl;
	}

	return returnValue;
}


real subpavings::tryExp(real r)
{
	real result = 0.0;
	try
	{
		result = _real(exp(_double(r)));

	}
	catch (...)
	{
		// exponentiation error

	}
	return result;

}


// turn a double into an rvector
rvector subpavings::makeDoubleIntoRvector (const double d)
{
	rvector newdata = _rvector(d);
	// 1-d rvector, only element is the implicit cast of d to real
	// note that cxsc allows an implicit cast from double to real
	return newdata;

}


// remove carriage returns
std::string& subpavings::reformatLine(std::string& line)
{
	string carriageRet("\r"); 
	
	size_t cr_pos = line.find_first_of(carriageRet);
	while (cr_pos != std::string::npos) {
		size_t next = line.find_first_not_of(carriageRet, cr_pos);
		line.erase(cr_pos, 1);
		if (next != std::string::npos) {
			cr_pos = line.find_first_of(carriageRet, next-1);
		}
		else cr_pos = next;
	}
	return line;
}

// quick check on a string expecting only numbers or spaces or tabs or 
// decimal points or commas or E or e characters
// used for checking text file input
bool subpavings::checkString(const std::string& s)
{
	// check for illegals, as anything but numbers, ".", e, E, ",", space or tab or carriage return
	string legal("eE+-.,0123456789 \t\r");
	string numeric("0123456789");
	string carriageRet("\r"); // only at the end
	bool retValue = ( (s.find_first_not_of(legal) == string::npos)
		&& (s.find_first_of(numeric) != string::npos));
	
	return retValue;
	// return true there are no illegal characters
}

// find number of blocks of numbers in string of numbers
// returns the number of blocks that have at least one 'number' in them
// blocks are separated by delimiters  space ' ', comma ',' or tab '\t'
int subpavings::countNumbers(const std::string& s)
{
	
	int countFound = 0;     // to count the number of finds
	
	// specify what to look for as numbers or decimal point or + or -
	string toFind(".0123456789");
	// and delimiters
	string delimFind(" \t,");
	// and illegals, as anything but numbers, ".", e, =, - space or tab
	//string illegalFind("e+-.0123456789 \t");

	size_t numberFound = s.find_first_of(toFind);    // first toFind
	size_t delimFound = 0;      //  first delimiter after number
	//size_t illegalFound = s.find_first_not_of(illegalFind);// first illegal

	while (numberFound != string::npos && delimFound != string::npos) {
		
		countFound++;
		
		delimFound = s.find_first_of(delimFind, numberFound+1);
		numberFound = s.find_first_of(toFind, delimFound); 
		
	}

	return countFound;
}

cxsc::rvector& subpavings::readRV(cxsc::rvector& r, const std::string& line)
{
	// convert line to an istream type
	istringstream sin(line);
	// c-xsc can convert this to an rvector
	sin >> r;
	// put r into the container
	return r;
}

cxsc::rvector& subpavings::readRV(cxsc::rvector& r, const std::string& line,
									const std::vector < int >& reqDims,
									int lineDim)
{
	// c-xsc can convert this to an rvector
	int n = reqDims.size();
	
	cxsc::rvector tmp(lineDim);
	
	tmp = readRV(tmp, line);
	
	for (int i = 1; i <= n; ++i) {
		r[i] = tmp[reqDims[i-1]];
	}

	return r;
}

void subpavings::lineErrorLogger(const std::string& line, int lineNumber)
{
	std::cerr << "Error in data input file, "
					<< "ignored line " << lineNumber
					<< ".  Data ignored is:  "
					<< line << std::endl;

}



bool subpavings::readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim)
{
	int checkLevel = 1;
	if (dim < 1) throw std::invalid_argument(
		"readRvectorsFromTxtOrd(RVecData&, const string&, std::size_t, int) : dim < 1");
	std::vector <int> reqDims;
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dim);
}

bool subpavings::readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines)
{
	int checkLevel = 1;
	int dataDim = -1;
	std::vector <int> reqDims;
	
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dataDim);
}

bool subpavings::readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							const std::vector < int >& reqDims,
							std::size_t headerlines)
{
	if (reqDims.empty()) {
		throw std::invalid_argument(
		"readRvectorsFromTxtOrd(RVecData&, const string&, const std::vector < int >&, std::size_t) : reqDims.empty()");
	}
	if ( (*std::min_element(reqDims.begin(), reqDims.end())) < 1) {
		throw std::invalid_argument(
		"readRvectorsFromTxtParanoid(RVecData&, const string&, const std::vector < int >&, std::size_t) : min element < 1");
	}
	
	int checkLevel = 1;
	int dataDim = reqDims.size();
	
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dataDim);
}

bool subpavings::readRvectorsFromTxtParanoid(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim)
{
	int checkLevel = 2;
	if (dim < 1) throw std::invalid_argument(
		"readRvectorsFromTxtParanoid(RVecData&, const string&, std::size_t, int) : dim < 1");
	std::vector <int> reqDims;
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dim);
}

bool subpavings::readRvectorsFromTxtParanoid(subpavings::RVecData& theData, 
							const std::string& s,
							const vector < int > & reqDims,
							std::size_t headerlines)
{
	if (reqDims.empty()) {
		throw std::invalid_argument(
		"readRvectorsFromTxtParanoid(RVecData&, const string&, const std::vector < int >&, std::size_t) : reqDims.empty()");
	}
	
	if ( (*std::min_element(reqDims.begin(), reqDims.end())) < 1) {
		throw std::invalid_argument(
		"readRvectorsFromTxtParanoid(RVecData&, const string&, const std::vector < int >&, std::size_t) : min element < 1");
	}
	
	int checkLevel = 2;
	int dataDim = reqDims.size();
	
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dataDim);
}

bool subpavings::readRvectorsFromTxtFast(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim)
{
	int checkLevel = 0;
	if (dim < 1) throw std::invalid_argument(
		"readRvectorsFromTxtFast(RVecData&, const string&, std::size_t, int) : dim < 1");
	std::vector <int> reqDims;
	return _readRvectorsFromTxt(theData, 
							s,
							headerlines,
							checkLevel,
							reqDims,
							dim);
}

// method to read rvectors from a txt file
bool subpavings::_readRvectorsFromTxt(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int checkLevel,
							const std::vector <int>& reqDims,
							int dim)
{
	bool retValue = false;
	bool cancontinue = true;
	
	int maxDim = 0;
	
	if (!reqDims.empty()) {
		dim = reqDims.size();
		maxDim = *max_element(reqDims.begin(), reqDims.end());
	}

	if (!theData.empty()) { // data already in theData, 
		std::cout
		//print a warning
			<< "Warning: adding to existing data "
			<< "- mixing datasets"
			<< std::endl;
		
		// and if dim is -1, set dim to be size of existing data
		if (dim < 0) dim =VecLen( (*theData.begin()) );
		if (dim < 1) throw std::invalid_argument(
		"_readRvectorsFromTxt(RVecData&, const std::string&, std::size_t, int, int) : VecLen(*theData.begin()) < 1");	
	}

	// we can read rvectors as strings in c-xsc format
	// and then convert to rvectors with c-xsc::>>

	//set up the file and read input line by line
	// we need to convert the string argument to a c-string for ifstream
	ifstream dataFile(s.c_str());
	int howMany = 0; // how many datapoints, to be read from file

	// if dim = -1, file will attempt to find it

	string line;
	
	std::streampos filePos = 0; // record where we are in file
	int lineDim = 0;
		
	if (dataFile.is_open())
	{
		std::size_t headers = headerlines;
		// skip headerlines
		while ((headers > 0) && (dataFile.good()) )
		{
			getline (dataFile,line);
			headers--;
		}
		filePos = dataFile.tellg( );
		
		if ( (dim == -1) || checkLevel || !reqDims.empty()) { // try to find dim from file
			// get the first data line from the file
			if (dataFile.good()) {
				
				getline (dataFile,line); 

				// find the number of blocks of numbers in first 
				// line that does not have illegal characters in
				while ( ( line.empty() || !checkString(line) ) && dataFile.good()) { // problem try next line
					filePos = dataFile.tellg( );
					getline(dataFile, line);
				}
				if ( !line.empty() && checkString(line) ) {
					lineDim = countNumbers(line);
					howMany++; // just read a valid line
				}
			}
		}

		if ((dim == -1) && (lineDim == 0)) {     // failed to find a valid line
			std::cerr << "Error in "
				<< "_readRvectorsFromTxt: "
				<< "could not establish data dimension from file "
				<< "- aborting read of file " << s << std::endl;
			dataFile.close();
			cancontinue = false;
		}
		else if ((checkLevel) && (lineDim < dim) ) { 
			// check first valid line matches given dimensions
			std::cerr << "Error in "
				<< "_readRvectorsFromTxt: "
				<< "data dimensions seems to be < given dimensions "
				<< "- aborting read of file " << s << std::endl;
			dataFile.close();
			cancontinue = false;
		}
			
		if (dim == -1) dim = lineDim;
		if (!maxDim) maxDim = dim;
		
		// dim now becomes the number of dimensions we expect
		
		// filePos is at the start of the first valid data line
	}

	else { // dataFile not open
	// todo should have this as our own subpavings io exception 
		std::cerr << "Error in "
			<< "_readRvectorsFromTxt."
			<< "Unable to open file " << s << std::endl;
		cancontinue = false;
	}

	// need to check we have enough values in the data
	if (cancontinue && !reqDims.empty()) {
		if (lineDim < maxDim) {
			std::cerr << "Error in "
				<< "_readRvectorsFromTxt: "
				<< "data dimensions seem to be < largest required dimension "
				<< "- aborting read of file " << s << std::endl;
			dataFile.close();
			cancontinue = false;
		}
	}

	// if there is already data, check dimensions match
	if (cancontinue && !theData.empty()) {

		//find the data dimensions from the first datapoint
		int dataDim = Ub(*theData.begin()) -
			Lb(*theData.begin()) + 1;
		if (dim != dataDim) {
			std::cerr << "Error in "
				<< "_readRvectorsFromTxt: "
				<< "Existing data different "
				<< "dimension to data to be read in "
				<< "- aborting read of file " << s << std::endl;
			dataFile.close();
			cancontinue = false;
		}
	}

	if (cancontinue) {
		// file is still open and we have just read the first valid line, for dims
		// or are positioned at the start, if dims are given
		
		// count remaining lines from here
		
		// count the lines still to go, without checking if they are good
		while (dataFile.good() )
		{
			getline (dataFile,line);
			howMany++;  // count number of lines in the file
		}	
			
		// reserve space in vector, assuming all remaining files will be good	
		theData.reserve(howMany + theData.size());

		dataFile.clear(); // reset the flags on the file
		dataFile.seekg(filePos); // and put file pointer to read first [valid] line

		int countIn = 0;
		int countLines = 0;

		while (dataFile.good() )
		{
			// get from the file line by line
			getline (dataFile, line);

			line = reformatLine(line); 

			countLines++;

			// validity checks only if checkLevel != 0
			if (!line.empty() && (!checkLevel || checkString(line)) ) { // quick validity check

				rvector r(dim); // empty rvector of correct dimension

				// the cxsc insertion operator >> is too forgiving
				// and idiosyncratic.  It will not just fail if the conversion
				// was not perfect, and it produces completely the wrong result if 
				// the final number is an integer not followed by anything
				
				bool proceed = (checkLevel < 2);
				int giveDim = lineDim;
				
				// only check dims using countNumbers if checkLevel >= 2
				if (!proceed) {
					
					giveDim = countNumbers(line);
					proceed = (reqDims.empty() ? (dim == giveDim) : maxDim <= (giveDim));
					
				}
				
				if ( proceed ) {
				
					// work around for the problem of no final decimal place
					// by just adding a space at the end - cxsc is okay with this
					line += " "; 

					if (reqDims.empty()) {
						r = readRV(r, line);
					}
					else {
						r = readRV(r, line, reqDims, giveDim);
					}
					theData.push_back(r);
					countIn++;
				}
				else {  // invalid line, spit out
					lineErrorLogger(line, countLines);
				}
			}
			else if (!line.empty() ) {  // invalid line, spit out
				lineErrorLogger(line, countLines);
			}
		}
		dataFile.close();

		// confirm the amount of data read in
		std::cout << "End of reading data input file: "
			<< countIn << " valid data points read in"
			<< std::endl;

		if (countIn > 0) retValue = true; // some data successfully read in
	}

	return retValue;
}




// method to read rvectors from a vector < vector < double > >
bool subpavings::getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
							const std::vector < subpavings::VecDbl > & inputData)
{
	int dim = -1;
	return _getRvectorsFromVectorOfVecDbl(theData, inputData, dim);
}

bool subpavings::getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
							const std::vector < subpavings::VecDbl > & inputData,
							int dim)
{
	if (dim < 1) throw std::invalid_argument(
		"readRvectorsFromVectorOfVecDbl(RVecData&, const std::vector < VecDbl >&, int) : dim < 1");
	return _getRvectorsFromVectorOfVecDbl(theData, inputData, dim);
}

bool subpavings::_getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
							const std::vector < subpavings::VecDbl > & inputData,
							int dim)
{
	bool retValue = false;
	bool cancontinue = true;
	if (!theData.empty()) { // data already in theData, print a warning
		std::cerr
			<< "Warning: adding to existing data "
			<< "- mixing datasets"
			<< std::endl;
		
		if (dim < 0) dim = VecLen(*theData.begin());
		
		if (dim < 1) throw std::invalid_argument(
		"_readRvectorsFromVectorOfVecDbl(RVecData&, const std::vector < VecDbl >&, int) : VecLen(*theData.begin()) < 1");	
	}

	/* we could convert to a stream and then use cxsc input operator
	 * but would then be subject to precision in conversion to stream,
	 * I think, so leave it as the explicit way for the moment . */

	if (inputData.empty()) {
		cancontinue = false;
	}
	
	// if we need to get data dim ourselves	
	if (cancontinue && dim < 0) {
		dim = inputData[0].size(); // size of the first one
		if (dim < 1) throw std::invalid_argument(
		"_readRvectorsFromVectorOfVecDbl(RVecData&, const std::vector < VecDbl >&, int) : inputData[0].empty()");	
		// dim now becomes the number of dimensions we expect
		
	}

	// if there is already data, check dimensions match
	if (cancontinue && !theData.empty()) {

		//find the existing data dimensions from the first datapoint
		int dataDim = VecLen(*theData.begin());
		if (dim != dataDim) {
			std::cerr
				<< "Existing data dimensions != "
				<< "dimensions of data to be read in" << std::endl;
			cancontinue = false;
		}
	}

	if (cancontinue) {
		
		size_t countIn = 0;
		size_t countLines = 0;
		
		theData.reserve(inputData.size() + theData.size());

		for (vector < VecDbl>::const_iterator it = inputData.begin();
			it < inputData.end();
			++it) {
				
			countLines++; // start at 1	
			
			if (it->size() >= static_cast<size_t>(dim)) {	
			
				rvector r(dim);
				theData.push_back(r);
				int index = 1;
				
				for (VecDbl::const_iterator dit = it->begin();
					dit < it->end(), index <= dim;
					++ dit) {
						
					theData.back()[index] = *dit;
					index++;
				}
					
				countIn++;
			}
			else {  // invalid input, spit out
				std::string line = toString(*it);
				lineErrorLogger(line, countLines);
			}
			
		}
		
		// confirm the amount of data read in
		std::cout << "End of getting data from vector: "
			<< countIn << " valid data points read in"
			<< std::endl;

		if (countIn > 0) retValue = true; // some data successfully read in
	}

	return retValue;
}


// method to get all rvectors from a container of rvectors
size_t subpavings::getRvectorsFromRVec(subpavings::RVecData& theData,
							const subpavings::RVecData& inputData,
							bool checkDims)
{

	size_t retValue = 0;
	bool cancontinue = true;
	int dim = 0;
	size_t existingSize = 0;
	
	if (!theData.empty()) { // data already in theData, print a warning
		existingSize = theData.size();
		std::cout
			<< "Warning: adding to existing data "
			<< "- mixing datasets"
			<< std::endl;

	}

	if (inputData.empty()) {
		cancontinue = false;
	}
	
	// check first point in file
	if (cancontinue) {
		dim = VecLen(*(inputData.begin())); // size of the first one
		if (dim < 1) throw std::invalid_argument(
		"_getRvectorsFromRVec(RVecData&, RVecData&, int) : dimensions of first input element < 1");	
		
	}

	if (cancontinue && !theData.empty()) { // check dimensions with first existing data point
		int dataDim = VecLen(*theData.begin());
		if (dim != dataDim) {
			std::cerr
				<< "Aborting:: Existing data dimensions !=  "
				<< "dimensions of first data point data to be read in" << std::endl;
			cancontinue = false;
		}
	}

	// dim is now irrelevant if checkDims = false
	
	// if cancontinue is still true we know there is data we can take
	if (cancontinue) {
		
		theData.reserve(inputData.size() + theData.size());

		if (!checkDims) {
			
			// just transfer the whole lot
			theData.insert(theData.end(), inputData.begin(),inputData.end());
		}
		else {
			int countLines = 0;
			for (RVecDataCItr it = inputData.begin();
					it < inputData.end();
					++it) {
				countLines++;
				if ( VecLen(*it) == dim ) theData.push_back(*it);
				else {
					std::string line = ::toString(*it);
					lineErrorLogger(line, countLines);
				}
			}
		}
		assert(theData.size() >= existingSize);
		
		retValue = theData.size() - existingSize;    // data read in
	}

	return retValue;
}


// method to get a sample of rvectors from a container
// the container data is assumed to be empty but
// should work even if it is not since we are only pushing back
// provided data dimensions match existing data if any
// takes a set-up rgsl
size_t subpavings::getSampleRvectorsFromRVec(subpavings::RVecData& data,
			gsl_rng * rgsl, size_t samplesize, const subpavings::RVecData& rvec)
{
	getSampleFromContainer(samplesize, rgsl, rvec, data);

	return data.size();
}


// method to get a sample of rvectors from an RSSample object
// the container data is assumed to be empty but
// should work even if it is not since we are only pushing back
// provided data dimensions match existing data if any
// takes a set-up rgsl
size_t subpavings::getSampleRvectorsFromRSSample(subpavings::RVecData& data,
			gsl_rng * rgsl, size_t samplesize, const RSSample& rss, int label)
{

	RVecData allData; // a container for all the rvectors from rss.Samples

	//use getRvectorsFromRSSampleForSampling to put rvectors from labeled
	//points in rss.Samples into allData where the labeled point label
	//matches label
	size_t numberAll = getRvectorsFromRSSampleForSampling(allData,
					data, samplesize, rss, label);

	// get the sample
	if (numberAll > 0) {
		getSampleFromContainer(samplesize, rgsl, allData, data);
	}

	return data.size();
}


// method to get all rvectors from an RSSample object
// the container data is assumed to be empty but
// should work even if it is not since we are only pushing back
// provided data dimensions match existing data if any
size_t subpavings::getRvectorsFromRSSample(subpavings::RVecData& data,
						const RSSample& rss, int label)
{

	size_t retValue = 0;
	bool cancontinue = true;

	if (!data.empty()) { // data already in data, print a warning
		std::cout
			<< "Warning: adding to existing data "
			<< "- mixing datasets"
			<< std::endl;

		// assume we will be adding all the data from rss.Samples
		data.reserve((rss.Samples).size() +
							data.size());
	}

	else {
		data.reserve((rss.Samples).size());
	}

	//get the data with the appropriate label out of rss

	// transform would be better but we would have to assume that
	// we wanted all the data in rss.Samples, but we know that
	// an RSSample object can contain samples with with different labels
	// and with different dimensions

	bool foundfirst = false;
	rvector first;

	if (cancontinue) {
		//find the first point with the right label and
		// check dimensions
		vector<LabPnt>::const_iterator cit = rss.Samples.begin();

		while (!foundfirst && (cit<(rss.Samples).end())) {
			if ((*cit).L == label) {
				foundfirst = true;
				first = cit->Pnt;
			}
			cit++;
		}// end while
	}

	if (!foundfirst) { // could not find points with the right label
		std::cout << "Could not find any points in RSSample object "
			<< " with the right label - aborting "
			<< std::endl;
		cancontinue = false;
	}

	if (cancontinue && !data.empty()) {
		//find the data dimensions from the first datapoint fount
		int dim = Ub(first) - Lb(first) + 1;
		//find the data dimensions from the existing data
		int dataDim = Ub(*data.begin()) - Lb(*data.begin()) + 1;
		if (dim != dataDim) {
			std::cout
				<< "Existing data different dimension "
				<< " to data in RSSample object - "
				<< " aborting"
				<< std::endl;
			cancontinue = false;
		}
	}

	// if cancontinue is still true we know there is data we can take

	if (cancontinue) {

		size_t countIn = 0; // track points inserted

		//iterate through the rss.Samples and take points
		//that match the supplied label
		vector<LabPnt>::const_iterator cit;

		for (cit=(rss.Samples).begin();cit<(rss.Samples).end();cit++){
			if (cit->L == label) {
				data.push_back(cit->Pnt);
				countIn++;
			}
		}

		retValue = countIn;    // some data successfully read in
	}

	return retValue;
}



// method to get data from an RSSample object to take samples from
size_t subpavings::getRvectorsFromRSSampleForSampling(subpavings::RVecData& allData,
					subpavings::RVecData& sampleData, size_t samplesize,
					const RSSample& rss, int label)
{

	size_t retValue = 0;

	//use getRvectorsFromRSSample to put rvectors from labeled points in
	// rss.Samples into allData where the labeled point label matches label
	size_t numberFound = getRvectorsFromRSSample(allData, rss, label);

	bool cancontinue = (numberFound > 0);
	// cancontinue will be false if there was a problem getting data points
	// if cancontinue is true data should contain at least some data points

	if (cancontinue && (allData.size() < samplesize)) {
		std::cout << "Warning: Sample size required is greater than "
			<< allData.size() << " the number of available data points "
			<< "matching label " << label << "in the RSSample object's Samples"
			<< std::endl;
	}

	// data already in data, print a warning
	if (cancontinue && !sampleData.empty()) {
		std::cout
			<< "Warning: adding to existing data "
			<< "- mixing datasets"
			<< std::endl;

		 //find the data dimensions from the first datapoint in allData
		int dim = Ub(*allData.begin()) - Lb(*allData.begin()) + 1;
		//find the data dimensions from the existing data
		int dataDim = Ub(*sampleData.begin()) - Lb(*sampleData.begin()) + 1;
		if (dim != dataDim) {
			std::cout
				<< "Existing data different dimension "
				<< " to data in RSSample object - aborting" << std::endl;
			cancontinue = false;
		}
	}

	if (cancontinue) {

		retValue = numberFound;    // some data successfully read in
	}

	return retValue;
}


void subpavings::getSampleFromContainer(size_t samplesize,
			gsl_rng * rgsl, const subpavings::RVecData& allData,
			subpavings::RVecData& sampleData)
{
	// make space to add the sample in addition to existing data
	if (sampleData.empty()) {
		sampleData.reserve(samplesize);
	}
	else {
		sampleData.reserve(sampleData.size() + samplesize);
	}

	//get a sample of the data out of allData
	for (size_t i = 0; i < samplesize; i++) {
		//draw a random number in [0,1)
		double rand = gsl_rng_uniform(rgsl);
		//turn this into an index in [0, samplesize-1]
		int index = static_cast<int>(ceil(rand*samplesize - 1));
		//put element in allData indexed into data
		sampleData.push_back(allData[index]);
	}
}

//new
// method to read rvectors from a txt file
bool subpavings::readVecRealVecFromTxt(std::vector < subpavings::RealVec >& theData, 
							const std::string& s, std::size_t headerlines, int dims)
{
	if (dims < 0 ) throw invalid_argument("readVecRealVecFromTxt(...) : dims < 0");
	
	size_t currentSize = 0;
	int existingDims = -1;
	bool canContinue = true;
	
	if (!theData.empty()) { // data already in theData, 
		std::cout << "\nWarning: adding to existing data " << std::endl;
		currentSize = theData.size();
		existingDims = theData.back().size();
	
		if (dims && dims != existingDims) {
			std::cerr
				<< "\nExisting data different "
				<< "dimension to dims given "
				<< "- aborting read of file " << s << std::endl;
			canContinue = false;
		}
		if (!dims) dims = existingDims;
	}

	if (canContinue) {
		ifstream dataFile(s.c_str());
		
		if (dataFile.is_open())
		{
			std::size_t headers = headerlines;
			// skip headerlines
			string line;
			
			while ((headers > 0) && (dataFile.good()) )
			{
				getline (dataFile,line);
				headers--;
			}
			std::streampos filePos = dataFile.tellg( );
		
			int dataDims = 0;
		
			while (dataFile.good() && !dataDims) {
				
				filePos = dataFile.tellg( );
				// get from the file line by line to find some dims
				getline (dataFile, line);

				if (!line.empty() && checkString(line)) dataDims = countNumbers(line);
			}
			
			canContinue = (dataDims > 0);
			
			if (!canContinue) { 
				std::cerr << "\nCould not find any data in the file to read in" << std::endl;
			}
			
			else if (dataDims < dims) {
				std::cerr
					<< "\nDimensions of data to be read in < expected data dimensions "
					<< "- aborting read of file " << s << std::endl;
				canContinue = false;
			}
			
			if (canContinue) {
				
				if (!dims) {
					dims = dataDims;
				}
				// otherwise dims will be as given or found from existing dims
				
				int howMany = 0;

				while (dataFile.good() )
				{
					getline (dataFile,line);
					howMany++;  // count number of lines in the file
				}
				
				// reserve space in vector, assuming all remaining files will be good	
				theData.reserve(howMany + currentSize);

				dataFile.clear(); // reset the flags on the file
				dataFile.seekg(filePos); // and put file pointer to read first non-empty line

				RealVec tmp (dims);

				while (dataFile.good() ) {
					dataFile >> ws;
				
					theData.push_back(tmp);
				
					for (int i = 0; i < dims; ++i) {
						dataFile >> (theData.back())[i];
					}
					for (int i = dims; i < dataDims; ++i) {
						cxsc::real tmp;
						dataFile >> tmp;
					}
				
					dataFile >> ws;
				}
			}
			
			dataFile.close();
		}
		else {
			std::cerr << "\nCould not open the file " << s << std::endl;
		}
	}
	
	return (theData.size() - currentSize > 0);
}


/* output a vector of vectors of reals to a file, with each inner vector
* as a column, headed by the column names in colNames,
* using precision given by prec.
* Checks that size of colNames matches size of dataPtrs,
* but does not check that vectors within dataPtrs are all the same length*/
void subpavings:: outputToFileVertical( std::vector < const subpavings::RealVec* >& dataPtrs, 
							const std::vector < std::string >& colNames,
							const std::string& filename,
							int prec)
{
	
	size_t n = colNames.size();
	if (n != dataPtrs.size()) {
		throw std::length_error("Size of colNames does not match size of dataPtrs");
	}
	
	outputFileStart(filename);
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	

	if (os.is_open()) {
		
		// use cxsc manipulators to set the precision for output
		os << cxsc::Variable;
		os << cxsc::SetPrecision(prec+2, prec);
		
		// colNames
		os << "line\t";
		ostream_iterator< std::string > out_s_it (os,"\t");
		copy ( colNames.begin(), colNames.end(), out_s_it );
		os << std::endl;
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < RealVec::const_iterator > dataIterators;
		
		for (std::vector < const RealVec* >::iterator it = dataPtrs.begin();
				it < dataPtrs.end(); ++it) {
					
			dataIterators.push_back( (*it)->begin() );
			if ((*it)->size() > longest) longest = (*it)->size();
		}
		// use this to output the data
		for (size_t j = 0; j < longest; ++j) {
			os << j+1;
			for (size_t i = 0; i < n; ++i) {
				if ( dataIterators[i] < dataPtrs[i]->end() ) {
					os << "\t" << (*dataIterators[i]);
					dataIterators[i]++;
				}
				else {
					os << "\t"; // not sure about this - may not be a good idea?
				}
				
			}
			// if it is not the last row, newline
			if ( j < longest - 1) os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}

/* output a vector of vectors of reals to a file, with each inner vector
* as a column, headed by the column names in colNames,
* using lhsLines as the line numbers on the left hand side,
* using precision given by prec.
* Checks that size of colNames matches size of dataPtrs,
* but does not check that vectors within dataPtrs are all the same length*/
void subpavings:: outputToFileVertical( const std::vector < size_t >& lhsLines,
							std::vector < const subpavings::RealVec* >& dataPtrs, 
							const std::vector < std::string >& colNames,
							const std::string& filename,
							int prec)
{
	
	
	
	size_t n = colNames.size();
	if (n != dataPtrs.size()) {
		throw std::length_error("Size of colNames does not match size of dataPtrs");
	}
	
	outputFileStart(filename);
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	

	if (os.is_open()) {
		
		// use cxsc manipulators to set the precision for output
		os << cxsc::Variable;
		os << cxsc::SetPrecision(prec+2, prec);
		
		// colNames
		os << "line\t";
		ostream_iterator< std::string > out_s_it (os,"\t");
		copy ( colNames.begin(), colNames.end(), out_s_it );
		os << std::endl;
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < RealVec::const_iterator > dataIterators;
		
		for (std::vector < const RealVec* >::iterator it = dataPtrs.begin();
				it < dataPtrs.end(); ++it) {
					
			dataIterators.push_back( (*it)->begin() );
			if ((*it)->size() > longest) longest = (*it)->size();
		}
		if (lhsLines.size() < longest)
			throw std::length_error("lhsLines is too short");
		// use this to output the data
		for (size_t j = 0; j < longest; ++j) {
			
			os << lhsLines[j];
			
			for (size_t i = 0; i < n; ++i) {
				if ( dataIterators[i] < dataPtrs[i]->end() ) {
					os << "\t" << (*dataIterators[i]);
					dataIterators[i]++;
				}
				else {
					os << "\t"; // not sure about this - may not be a good idea?
				}
				
			}
			// if it is not the last row, newline
			if ( j < longest - 1) os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}


/* output a vector of vectors of reals to a file, with each inner vector
* as a column, headed by the column names in colNames,
* using precision given by prec.
* Checks that size of colNames matches size of dataPtrs,
* but does not check that vectors within dataPtrs are all the same length*/
void subpavings:: outputToFileVertical( 
						std::vector < const subpavings::Size_tVec* >& dataPtrs, 
						const std::vector < std::string >& colNames,
						const std::string& filename)
{
	
	size_t n = colNames.size();
	if (n != dataPtrs.size()) {
		throw std::length_error("Size of colNames does not match size of dataPtrs");
	}
	
	outputFileStart(filename);
	
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	

	if (os.is_open()) {
		
		
		// colNames
		os << "line\t";
		ostream_iterator< std::string > out_s_it (os,"\t");
		copy ( colNames.begin(), colNames.end(), out_s_it );
		os << std::endl;
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < Size_tVec::const_iterator > dataIterators;
		for (std::vector < const Size_tVec* >::iterator it = dataPtrs.begin();
				it < dataPtrs.end(); ++it) {
			dataIterators.push_back( (*it)->begin() );
			if ((*it)->size() > longest) longest = (*it)->size();
		}
		// use this to output the data
		for (size_t j = 0; j < longest; ++j) {
			os << j+1;
			for (size_t i = 0; i < n; ++i) {
				if ( dataIterators[i] < dataPtrs[i]->end() ) {
					os << "\t" << (*dataIterators[i]);
					dataIterators[i]++;
				}
				else {
					os << "\t"; // not sure about this - may not be a good idea?
				}
				
			}
			// if it is not the last row, newline
			if ( j < longest - 1) os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}

void subpavings:: outputToFileVertical( 
							std::vector < const subpavings::Size_tVec* >& sizePtrs,
							std::vector < const subpavings::RealVec* >& realPtrs, 
							const std::vector < std::string >& colNamesSize,
							const std::vector < std::string >& colNamesReal,
							const std::string& filename,
							int prec)
{
	size_t n1 = colNamesSize.size();
	if (n1 != sizePtrs.size()) {
		throw std::length_error("Size of colNamesSize does not match size of sizePtrs");
	}
	size_t n2 = colNamesReal.size();
	if (n1 != n2) {
		throw std::length_error("Size of colNamesReal does not match size of colNamesSize");
	}
	if (n2 != realPtrs.size()) {
		throw std::length_error("Size of colNamesReal does not match size of realPtrs");
	}
	size_t n = n1;
	
	outputFileStart(filename);
	
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	

	if (os.is_open()) {
		
		// use cxsc manipulators to set the precision for output
		os << cxsc::Variable;
		os << cxsc::SetPrecision(prec+2, prec);
		
		// colNames
		os << "line\t";
		for (size_t i = 0; i < n-1; ++i) { // all except last
			os << colNamesSize[i] << "\t" << colNamesReal[i] << "\t";
		}
		// last one
		os << colNamesSize[n-1] << "\t" << colNamesReal[n-1] << endl;
			
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < RealVec::const_iterator > realIterators;
		for (std::vector < const RealVec* >::iterator it = realPtrs.begin();
				it < realPtrs.end(); ++it) {
			realIterators.push_back( (*it)->begin() );
			if ((*it)->size() > longest) longest = (*it)->size();
		}
		
		std::vector < Size_tVec::const_iterator > sizeIterators;
		for (std::vector < const Size_tVec* >::iterator it = sizePtrs.begin();
				it < sizePtrs.end(); ++it) {
			sizeIterators.push_back( (*it)->begin() );
			if ((*it)->size() > longest) longest = (*it)->size();
		}
		
		// use this to output the data
		for (size_t j = 0; j < longest; ++j) {
			os << j+1; // line number
			for (size_t i = 0; i < n; ++i) {
				if ( sizeIterators[i] < sizePtrs[i]->end() ) {
					os << "\t" << (*sizeIterators[i]);
					sizeIterators[i]++;
				}
				else {
					os << "\t"; 
				}
				if ( realIterators[i] < realPtrs[i]->end() ) {
					os << "\t" << (*realIterators[i]);
					realIterators[i]++;
				}
				else {
					os << "\t";
				}
				
			}
			// newline
			os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
	
}

/* output to file vertical assuming we are adding to an existing file*/
void subpavings:: outputToFileVertical( 
							std::vector < const subpavings::RealVec* >& realPtrs, 
							const std::string& filename,
							size_t startPos, 
							int prec)
{
	size_t n = realPtrs.size();
	
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	
	if (os.is_open()) {
		
		// use cxsc manipulators to set the precision for output
		os << cxsc::Variable;
		os << cxsc::SetPrecision(prec+2, prec);
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < RealVec::const_iterator > realIterators;
		for (std::vector < const RealVec* >::iterator it = realPtrs.begin();
				it < realPtrs.end(); ++it) {
			RealVec::const_iterator rit = (*it)->begin();
			size_t thisSize = (*it)->size();
			if (thisSize > startPos) advance(rit, startPos);		
			else rit = (*it)->end();	
					
			realIterators.push_back( rit );
			if (thisSize > longest) longest = (*it)->size();
		}
		
		if (longest < startPos)
			throw std::invalid_argument(
						"outputToFileVertical(...) : startPos");
		
		// use this to output the data
		for (size_t j = 0; j < longest-startPos; ++j) {
			os << j+1+startPos; // line number
			for (size_t i = 0; i < n; ++i) {
				if ( realIterators[i] < realPtrs[i]->end() ) {
					os << "\t" << (*realIterators[i]);
					realIterators[i]++;
				}
				else {
					os << "\t";
				}
				
			}
			// newline
			os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
	
}


/* output to file vertical assuming we are adding to an existing file*/
void subpavings:: outputToFileVertical( 
							std::vector < const subpavings::Size_tVec* >& sizePtrs,
							std::vector < const subpavings::RealVec* >& realPtrs, 
							const std::string& filename,
							size_t startPos, 
							int prec)
{
	size_t n = realPtrs.size();
	
	std::ofstream os;
	
	os.open(filename.c_str(), ios::app);         // append
	
	if (os.is_open()) {
		
		// use cxsc manipulators to set the precision for output
		os << cxsc::Variable;
		os << cxsc::SetPrecision(prec+2, prec);
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < RealVec::const_iterator > realIterators;
		for (std::vector < const RealVec* >::iterator it = realPtrs.begin();
				it < realPtrs.end(); ++it) {
			RealVec::const_iterator rit = (*it)->begin();
			size_t thisSize = (*it)->size();
			if (thisSize > startPos) advance(rit, startPos);		
			else rit = (*it)->end();	
					
			realIterators.push_back( rit );
			if (thisSize > longest) longest = (*it)->size();
		}
		
		std::vector < Size_tVec::const_iterator > sizeIterators;
		for (std::vector < const Size_tVec* >::iterator it = sizePtrs.begin();
				it < sizePtrs.end(); ++it) {
					
			Size_tVec::const_iterator sit = (*it)->begin();
			size_t thisSize = (*it)->size();
			if (thisSize > startPos) advance(sit, startPos);		
			else sit = (*it)->end();	
					
			sizeIterators.push_back( sit );
			if (thisSize > longest) longest = (*it)->size();
		}
		
		if (longest < startPos)
			throw std::invalid_argument(
						"outputToFileVertical(...) : startPos");
		
		// use this to output the data
		for (size_t j = 0; j < longest-startPos; ++j) {
			os << j+1+startPos; // line number
			for (size_t i = 0; i < n; ++i) {
				if ( sizeIterators[i] < sizePtrs[i]->end() ) {
					os << "\t" << (*sizeIterators[i]);
					sizeIterators[i]++;
				}
				else {
					os << "\t"; 
				}
				if ( realIterators[i] < realPtrs[i]->end() ) {
					os << "\t" << (*realIterators[i]);
					realIterators[i]++;
				}
				else {
					os << "\t";
				}
				
			}
			// newline
			os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
	
}


std::vector < const RealVec* >& subpavings::addDataPtrs(
								std::vector < const subpavings::RealVec* >& container,
								const std::vector < subpavings::RealVec >& toAdd)
{
	for (std::vector < RealVec >::const_iterator it = toAdd.begin();
					it < toAdd.end();
					++it) {
		
		container.push_back( &(*it) );
						
	}
	
	return container;

}

std::vector < const Size_tVec* >& subpavings::addDataPtrs(
								std::vector < const subpavings::Size_tVec* >& container,
								const std::vector < subpavings::Size_tVec >& toAdd)
{
	for (std::vector < Size_tVec >::const_iterator it = toAdd.begin();
					it < toAdd.end();
					++it) {
		
		container.push_back( &(*it) );
						
	}
	
	return container;

}

std::ostream & subpavings::operator<<(std::ostream &os, const subpavings::IntVec& vec)
{
	if (!vec.empty()) {
		std::ostream_iterator< int > out_it (os, "\t");
		copy ( vec.begin(), vec.end(), out_it );
	}
	return os;
}

std::ostream & subpavings::operator<<(std::ostream &os, const subpavings::RealVec& vec)
{
	if (!vec.empty()) {
		std::ostream_iterator< cxsc::real > out_it (os, "\t");
		copy ( vec.begin(), vec.end(), out_it );
	}
	return os;
}

std::ostream & subpavings::operator<<(std::ostream &os, const subpavings::VecDbl& vec)
{
	if (!vec.empty()) {
		std::ostream_iterator< double > out_it (os, "\t");
		copy ( vec.begin(), vec.end(), out_it );
	}
	return os;
}

std::string subpavings::toString(const subpavings::IntVec vec, 
														bool compact)
{
	std::ostringstream oss;
	
	std::string delim(", ");
	if (compact) delim = ",";
	
	std::ostream_iterator<int> out_it (oss, delim.c_str());
	copy ( vec.begin(), vec.end(), out_it );
	
	std::string result = oss.str();

	result = result.substr(0, result.size()-delim.size());
	
	return result;
}

std::string subpavings::toString(const subpavings::RealVec vec,
														bool compact)
{
	std::ostringstream oss;
	
	std::string delim(", ");
	if (compact) delim = ",";
	
	std::ostream_iterator<cxsc::real> out_it (oss, delim.c_str());
	copy ( vec.begin(), vec.end(), out_it );
	
	std::string result = oss.str();

	result = result.substr(0, result.size()-delim.size());
	
	return result;
}

std::string subpavings::toString(const subpavings::VecDbl vec,
														bool compact)
{
	std::ostringstream oss;
	
	std::string delim(", ");
	if (compact) delim = ",";
	
	std::ostream_iterator<double> out_it (oss, delim.c_str());
	copy ( vec.begin(), vec.end(), out_it );
	
	std::string result = oss.str();

	result = result.substr(0, result.size()-delim.size());
	
	return result;
}

//new AHABC
// a sorting order for rvectors
// sorts by element by element comparison
bool RvecComp::operator() (const cxsc::rvector& lhs, const cxsc::rvector& rhs) const
{
	int len = VecLen(lhs);
	if ( len != VecLen(rhs) ) {
		return ( len < VecLen(rhs) );
	}
	// they are the same length
	
	if (len <= 0) {
		return false;
	}
	
	// same length and at least one element
	
	int index = 1;
	while (index <= len) {
		if (lhs[index] < rhs[index]) return true;
		if (rhs[index] < lhs[index]) return false;
		index ++;
	}
	// only get to here if everything equal
	return false;

}

//new AHABC
// Point mass filtering from data in a RVecData container.
std::map<rvector, size_t, subpavings::RvecComp >& subpavings::pointMassFilter(
		const subpavings::RVecData& theData, 
		std::map<rvector, size_t, subpavings::RvecComp > & countsMap)
{
	RVecDataCItr dataIt;
	std::pair< map< rvector, size_t, RvecComp >::iterator, bool> boolCount;
	
	//go through the dataset
	for (dataIt = theData.begin(); dataIt < theData.end(); ++dataIt) {
		
		//insert into map
		boolCount = countsMap.insert(make_pair(*dataIt,1));		
		//Check if insertion is successful - a unique set will render a successful insertion.
		if(!(boolCount.second)) { //if there are repeated points
		
			boolCount.first->second++;
		} 
	} // end of going through each data in this dataset
	return countsMap;
}


//new AHABC
// Create an RSSample object from a counts Map
RSSample& subpavings::labelDataFromFilter(RSSample&  labData,
		const std::map<rvector, size_t, subpavings::RvecComp >& countsMap,
		int uniqueLabel, int ptMassLabel)
{
	
	std::map<rvector, size_t, RvecComp >::const_iterator it;
	
	std::vector<LabPnt> tmp;
	
	tmp.reserve(countsMap.size() * 3);  // this is just an approximate guess
	
	/* for every point in the map, add to the sample, count times */
	for (it = countsMap.begin(); it != countsMap.end(); ++it) {
		LabPnt data0;
		data0.Pnt = it->first; // the rvector
		data0.L = uniqueLabel;
		size_t count = it->second;
		if (count > 1) data0.L = ptMassLabel;
		
		for (size_t i = 0; i < count; ++i) {
			tmp.push_back(data0);
		}
	}
	
	labData.Samples.swap(tmp);
	
	return labData;
}
	


//new AHABC
// Fills in an EMFMap from a countsMap
// return the total point mass weight
double subpavings::makeEMFMap(
		std::map < rvector, double, subpavings::RvecComp >& EMFMap,
		const std::map<rvector, size_t, subpavings::RvecComp >& countsMap)
{
	double pmWeight = 0.0;
	
	std::map<rvector, size_t, RvecComp >::const_iterator it;
	
	std::map < rvector, double, RvecComp > tmp;
	
	size_t total = 0;
	for (it = countsMap.begin(); it != countsMap.end(); ++it) {
		total += it->second;
	}
		
	/* for every non-unique point in the map, add to the EMF map */
	for (it = countsMap.begin(); it != countsMap.end(); ++it) {
		size_t count = it->second;
		if (count > 1 ) {
			double weight = ( static_cast<double>(count) )/total;
			pmWeight += weight;
			tmp.insert( std::make_pair(it->first, weight) );
			
		}
	}
	tmp.swap(EMFMap);
	
	return pmWeight;
}

//mean and sample standard deviation	
std::pair < cxsc::real, cxsc::real > subpavings::calcMeanAndSD(
							const RealVec& vec) 
{
	size_t n = vec.size();
	
	cxsc::dotprecision sum_sqs(0.0);
	cxsc::real sum = 0.0;
	for (std::vector < cxsc::real >::const_iterator it = vec.begin();
			it < vec.end();
			++it)
	{
		accumulate(sum_sqs, (*it), (*it));
		sum += (*it);
	}
	
	cxsc::real mean = sum/cxsc::real(1.0*n);
	cxsc::real var = 0.0;
	
	if (n > 1) {
		accumulate(sum_sqs, -sum, mean);
		var = rnd(sum_sqs)/cxsc::real(n-1.0);
		
	}
	return make_pair( mean, cxsc::sqrt(var) );
	
}


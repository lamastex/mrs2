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

/*!/ \file:     sptools.cpp
\brief Implementation of sptools functions
*/

#include "sptools.hpp"

// include fstream so as to be able to output a file
#include <fstream>

// to be able to manipulate strings as streams
#include <sstream>

#include <stdexcept>
#include <cassert>

// general interval tools
#include "toolz.hpp"

// to use LabBox and RSSample objects
#include "SmallClasses.hpp"

using namespace std;
using namespace subpavings;

    // make a unique file name that has a time stamp number appended to it
    string subpavings::getUniqueFilename(string baseFileName,string suffix)
    {
        string s;
        bool newfilename = false;

        while (!newfilename) {
            stringstream out;
            out << time(NULL);
            s = baseFileName + out.str() + suffix;
            //try opening a file with this name
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

    // Method to add a line to a file
    // Output goes to file named according to argument s
    void subpavings::outputFile(const std::string& s, const string line, bool append)
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
            std::cout << "Error: could not open file named "
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
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
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
            vals.clear();

            os.close();
        }
        else {
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }

    // Method to append values to output log file
    // Output goes to file named according to argument s
    void subpavings::outputFile(const std::string& s, const std::string intro,
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
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }


    // Method to append values to output log file
    // Output goes to file named according to argument s
    void subpavings::outputFile(const std::string& s, const std::string intro,
                    IntVec& vals)
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
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }

    // Method to append values to output log file
    // Output goes to file named according to argument s
    void subpavings::outputFile(const std::string& s, const std::string intro,
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
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }

    // Method to append strings to output log file
    // Output goes to file named according to argument s
    void subpavings::outputFile(const std::string& s, vector<string>& strings, int i)
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
            std::cout << "Error: could not open file named "
                << s << std::endl << std::endl;
        }
    }

    // parse a tree string to make a .dot file
    bool subpavings::parseForGraphDot(string s, string toParse)
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
    void subpavings::makeDotImage(string s)
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
    int subpavings::countLinesInTxt(const string& s)
    {
        //set up the file and read input line by line
        // we need to convert the string argument to a c-string for ifstream
        ifstream dataFile(s.c_str());
        string line;    // a string object to use in counting lines
        int howMany = 0; // how many lines in the file

        if (dataFile.is_open())
        {
            // count the lines in the file
            while (dataFile.good() )
            {
                getline (dataFile,line);
                howMany++;  // count the number of lines in file
            }
        }

        else {
            std::cout << "Error in AdaptiveHistogram::countLinesInTxt. "
                << "Unable to open file" << std::endl;
            exit(1);
        }

        dataFile.close();

        return howMany;
    }


    // parse a vtk file dimension spacings line for spacings
    IntVec subpavings::parseSpacings(string line, IntVec& spacings)
    {
        // specify what to look for
        string nums("0123456789");
        string space(" \t");

        int found = 0;
        istringstream sin;

        size_t firstFound = line.find_first_of(nums);// first number

        // go through the line extracting integer numbers
        while (firstFound!=string::npos) {
            size_t endpos = line.find_first_of(space, firstFound);
            if (endpos!=string::npos) {
                istringstream sin(line.substr(firstFound, endpos-firstFound));
                sin >> found;  // extract number from the stream
                firstFound = line.find_first_of(nums, endpos);
            }
            else {
                istringstream sin(line.substr(firstFound));
                sin >> found;  // extract number from the stream
                firstFound = string::npos;
            }

            spacings.push_back(found);
            found = 0;
        }

        return spacings;
    }

    bool findBlack(string line, string seek)
    {
        bool found = false;
        size_t start = line.find(seek.substr(0,1));
        if (start != string::npos) {
            size_t startseek = 1;
            string seekrest = seek.substr(startseek,string::npos);
            found = true;
            while (seekrest.length() > 1 && found) {
                if (line.find(seekrest.substr(1,string::npos)) != start+1)
                    found = false;
                startseek ++;
                seekrest = seek.substr(startseek,string::npos);
            }
        }
        return found;
    }

    // method to read coordinates from a .vtk file
    // there is some checking on the data input
    // expects structured point format data for 3 dimensions
    IntVec subpavings::getCoordinatesFromVtk(IntVec& Xs, IntVec& Ys, IntVec& Zs,
                                                const string& s)
    {
        size_t resCap = 1000; // guess about how much capacity to reserve
        int headerLines = 10; // number of header lines expected
        int spacingHeader = 4; // line number on which spacings expected
        IntVec spacings; // to hold dimensions
        size_t x_dim = 0;  // x-dimensions
        size_t y_dim = 0;  // y-dimensions
        size_t z_dim = 0;  // z-dimensions

        size_t expectedDims = 3;     // expected dimensions
        string seek = "255"; // the string indicating a voxel in the image

        bool retValue = false;
        bool cancontinue = true;

        //set up the file and read input line by line
        // we need to convert the string argument to a c-string for ifstream
        ifstream dataFile(s.c_str());
        size_t howMany = 0; // how many datapoints, to be read from file

        string line;

        if (!dataFile.is_open())
        {
            std::cout << "Error in "
                << "SPnode::getCoordinatesFromVtk: "
                << "Unable to open file " << s << std::endl;
            cancontinue = false;
        }

        int lineNumber = 0;
        if (cancontinue) {

            // find the spacings
            while (lineNumber < spacingHeader && dataFile.good()) {
                getline(dataFile, line);
                lineNumber++;
            }
            if (lineNumber == spacingHeader)
            {
                getline(dataFile, line);
                lineNumber++;
                // line should contain spacing data
                // parse out the dimensions
                spacings = parseSpacings(line, spacings);
            }
            if (spacings.size() != expectedDims) cancontinue = false;
        }

        if (cancontinue) {

            // spacings should now have dimensions (x_dim, y_dim, z_dim)
            x_dim = spacings[0]; //x-dimensions
            y_dim = spacings[1]; //y-dimensions
            z_dim = spacings[2]; //z-dimensions

            // get to the top of the data line
            while (lineNumber < headerLines && dataFile.good()) {
                getline(dataFile, line);
                lineNumber++;
            }
        }
        if (cancontinue && (lineNumber == headerLines) && dataFile.good()) {
            // reserve capacity in vectors (this is just an over-estimate guess)
            Xs.reserve(resCap);
            Ys.reserve(resCap);
            Zs.reserve(resCap);

            int coordNumber = 0;

            // while the datafile is good, extract coordinates
            while (dataFile.good() )
            {
                getline (dataFile,line);
                // check if line contains the string we are seeking
                if (findBlack(line, seek)) {
                    // using integer division here
                    int x_coord = coordNumber/(y_dim*z_dim);
                    int y_coord = (coordNumber - x_coord*y_dim*z_dim)/z_dim;
                    int z_coord = coordNumber - x_coord*y_dim*z_dim
                                - y_coord*z_dim;
                    Xs.push_back(x_coord);
                    Ys.push_back(y_coord);
                    Zs.push_back(z_coord);

                }
                coordNumber++;
            }

            // check the amount of data read in
            if (coordNumber - 1 == x_dim*y_dim*z_dim) retValue = true;
        }

        dataFile.close();

        return spacings;
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
            std::cout << "Error in volCompare : comparing "
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


    // quick check on a string expecting only numbers white space or decimal points
    // checks for illegal characters and check number of decimal points found is n
    // used for checking text file input
    bool subpavings::checkString(const string& s, const int n)
    {
        int countFound = 0;     // to count the number of finds

        string toFind("."); // look for occurrences of a decimal point

        size_t posFound = s.find(toFind); // find first occurrence of toFind

        while(posFound!=string::npos) {  // count the occurrences of toFind
            countFound++;
            posFound = s.find(toFind, (posFound + 1));
        }

        // check for illegals, as anything but numbers, ".", e, space or tab
        string illegalFind("e+-.0123456789 \t");
        if ((s.find_first_not_of(illegalFind))!=string::npos) {
            countFound = 0;
        }

        return (countFound == n);
        // return true if the number of decimal points found is n and
        // there are no illegal characters
    }

    // find number of blocks of numbers in a properly formatted string of numbers
    // properly formatted means each number includes a decimal point
    // returns 0 if not all numbers have a decimal point
    // or if the string contains anything not in ".123456789 \t"
    // used for checking text file input
    int subpavings::countNumbers(const string& s)
    {
		
		int countFound = 0;     // to count the number of finds
        int countDec = 0;   // to count the number of decimal places

        // specify what to look for as numbers or decimal point or + or -
        string toFind("e+-.0123456789");
        // and decimals only
        string decFind(".");
        // and illegals, as anything but numbers, ".", e, =, - space or tab
        string illegalFind("e+-.0123456789 \t");

        size_t firstFound = s.find_first_of(toFind);    // first toFind
        size_t decFound = s.find(decFind);      //  first decFind
        size_t illegalFound = s.find_first_not_of(illegalFind);// first illegal

        while (firstFound!=string::npos && illegalFound==string::npos) {

            if (decFound!=string::npos) {
                countDec++;
            }

            countFound++;

            // find first something else
            size_t otherFound = s.find_first_not_of(toFind,
                                                    (firstFound + 1));

            // if not already at end search again for toFind and  decFind
            if (otherFound!=string::npos) {
                firstFound = s.find_first_of(toFind, (otherFound + 1));
                decFound = s.find(decFind, (otherFound + 1));
            }
            else {
                // terminate while loop if otherFound==string::npos
                firstFound = string::npos;
            }
        }

        // if illegal character in s, countFound and decFound are still be 0

        // check if # of decimal points matches the number of number blocks
        if (countDec!=countFound) {
            countFound = 0; // if not, countFound will be 0
        }
		
		  

        return countFound;
        // return the number of blocks of numbers in a legal string
        // returns 0 if any numbers don't have decimal points
        // or if there are decimal points without numbers
        // of if there are any illegal characters in the string
    }


    // method to read 1-d data from a txt file and convert to 1-d rvectors
    // there is some checking on the data input
    // Only suitable for for 1-d data
    // expects values for a single integer or double on one line
    // Will ignore anything else on the line
    // the container theData is assumed to be empty but
    // should work even if it is not since we are only pushing back
    bool subpavings::readOneDimDataFromTxt(RVecData& theData, 
									const string& s,
									const size_t headerlines)
    {
        bool retValue = false;
        bool cancontinue = true;
        int expectedDims = 1;

        if (!theData.empty()) { // data already in theData, print a warning
            std::cout
                << "Warning: adding to existing data "
                << "- mixing datasets"
                << std::endl;

            //find the data dimensions from the first datapoint
            int dataDim = Ub(*theData.begin()) -
                Lb(*theData.begin()) + 1;
            if (dataDim != expectedDims) {
                std::cout
                    << "Existing data not one-dimensional"
                    << std::endl;
                cancontinue = false;
            }
        }

        if (cancontinue) {

            //set up the file and read input line by line
            // we need to convert the string argument to a c-string for ifstream
            ifstream dataFile(s.c_str());
            string line;

            if (!dataFile.is_open())
            {
                std::cout << "Error in "
                    << "AdaptiveHistogram::readOneDimDataFromTxt."
                    << "Unable to open file" << std::endl;
                cancontinue = false;
            }

            else {
				
				std::size_t headers = headerlines;
				// skip headerlines
                while ((headers > 0) && (dataFile.good()) )
                {
                    getline (dataFile,line);
					headers--;
                }

                size_t howMany = 0;

                // count the lines in the file
                while (dataFile.good() )
                {
                    getline (dataFile,line);
                    howMany++;  // count number of lines in the file
                }

                theData.reserve(howMany);

                dataFile.clear(); // reset the flags on the file
                dataFile.seekg(0, ios::beg); // and put file pointer to start
				
				headers = headerlines;
				// skip headers again
                while ((headers > 0) && (dataFile.good()) )
                {
                    getline (dataFile,line);
					headers--;
                }

                int countIn = 0; // to keep track of lines successfully read in
                int countLines = 0; // to keep track of line numbers

                while (dataFile.good() )
                {
                    // get from the file line by line
                    getline (dataFile, line);
                    countLines++;

                    // convert to an istream type
                    istringstream sin(line);

                    int intIn = 0; // to hold int expected
                    bool goodFind = false; // boolean for int or double found
                    rvector r(expectedDims); // to hold conversion to rvector

                    sin >> intIn; // first number on line only read
                    if (sin) // no error extracting line content to an int
                    {
                        // convert this to an rvector via a real
                        real realIn(intIn);
                        r[1] = realIn;
                        goodFind = true;

                    }

                    else { // error so try a double instead
                        double doubleIn = 0.0;

                        sin.clear(std::ios::goodbit);
                        sin >> doubleIn;

                        if (sin) // no error this time
                        {
                            r(doubleIn); // convert to rvector
                            goodFind = true;
                        }
                    }

                    if (goodFind) {
                        // put r into the container
                        theData.push_back(r);
                        countIn++;

                    }

                    else if (line != "") {  // invalid line, spit out
                        std::cout << "Error in data input file, "
                            << "ignored line " << countLines
                            << ".  Data ignored is:  "
                            << line << std::endl;
                        // better make this output to an input error
                        //log txt file?
                    }
                }
                dataFile.close();

                // confirm the amount of data read in
                std::cout << "End of reading data input file: "
                    << countIn << " valid data points read in"
                    << std::endl;

                if (countIn > 0) retValue = true;    // some data successfully read in
            }
        }

        return retValue;
    }



    // method to read rvectors from a txt file
    // there is some checking on the data input
    // This works for 1-d data
    // expects values for a single rvector on one line, separated by white space
    // Will spit out data which does not match this, but keep reading the rest
    // the container theData is assumed to be empty but
    // should work even if it is not since we are only pushing back
    bool subpavings::readRvectorsFromTxt(RVecData& theData, 
								const string& s,
								const std::size_t headerlines)
    {
        bool retValue = false;
        bool cancontinue = true;

        if (!theData.empty()) { // data already in theData, print a warning
            std::cout
                << "Warning: adding to existing data "
                << "- mixing datasets"
                << std::endl;
        }

        // we can read rvectors as strings in c-xsc format
        // and then convert to rvectors with c-xsc::>>

        //set up the file and read input line by line
        // we need to convert the string argument to a c-string for ifstream
        ifstream dataFile(s.c_str());
        int dim = 0; // dimensions, to be assessed from the file
        int howMany = 0; // how many datapoints, to be read from file

        string line;

        if (dataFile.is_open())
        {
			std::size_t headers = headerlines;
			// skip headerlines
			while ((headers > 0) && (dataFile.good()) )
			{
				getline (dataFile,line);
				headers--;
			}
			
			// get the first data line from the file
            if (dataFile.good()) {
				getline (dataFile,line); 

				if(line.empty()) { // if no characters extracted
					std::cout
						<< "Error in "
						<< "AdaptiveHistogram::readRvectorsFromTxt: "
						<< "no data in input file " << s << std::endl;
					dataFile.close();
					cancontinue = false;
				}
            }

        }

        else { // dataFile not open
            std::cout << "Error in "
                << "AdaptiveHistogram::readRvectorsFromTxt."
                << "Unable to open file " << s << std::endl;
            cancontinue = false;
        }

        if (cancontinue) {

            //  find the number of blocks of numbers with decimal points
            // countNumbers will return 0 if the line contains illegal
            // characters or any number is missing its decimal point

            dim = countNumbers(line);
            while (dim == 0 && dataFile.good()) { // problem try next line
                getline(dataFile, line);
                dim = countNumbers(line);
            }

            if (dim == 0) {     // failed to find a valid line
                std::cout << "Error in "
                    << "AdaptiveHistogram::readRvectorsFromTxt: "
                    << "all lines of input file " 
					<< s << " contain illegal formatting" << std::endl;
                dataFile.close();
                cancontinue = false;
            }

            // dim now becomes the number of dimensions we expect
        }

        // if there is already data, check dimensions match
        if (cancontinue && !theData.empty()) {

            //find the data dimensions from the first datapoint
            int dataDim = Ub(*theData.begin()) -
                Lb(*theData.begin()) + 1;
            if (dim != dataDim) {
                std::cout
                    << "Existing data different "
                    << "dimension to data to be read in "
                    << "- aborting read of file " << s << std::endl;
                dataFile.close();
                cancontinue = false;
            }
        }

        if (cancontinue) {

            dataFile.clear(); // reset the flags on the file
            dataFile.seekg(0, ios::beg); // and put file pointer to start
			
			std::size_t headers = headerlines;
			// skip headers again
			while ((headers > 0) && (dataFile.good()) )
			{
				getline (dataFile,line);
				headers--;
			}

            rvector r(dim);

            // count the lines in the file
            while (dataFile.good() )
            {
                getline (dataFile,line);
                howMany++;  // count number of lines in the file
            }

            theData.reserve(howMany);

            dataFile.clear(); // reset the flags on the file
            dataFile.seekg(0, ios::beg); // and put file pointer to start
			
			headers = headerlines;

			// skip headers again
			while ((headers > 0) && (dataFile.good()) )
			{
				getline (dataFile,line);
				headers--;
			}

            int countIn = 0;
            int countLines = 0;

            while (dataFile.good() )
            {
                // get from the file line by line
                getline (dataFile, line);
                countLines++;

                if (checkString(line, dim)) { // quick validity check

                    // could replace all this by checking if sin is true (good)
                    // since it will not be if failed to convert input to r
                    // but should check this before implementing

                    // convert to an istream type
                    istringstream sin(line);
                    // c-xsc can convert this to an rvector
                    sin >> r;
                    // put r into the container
                    theData.push_back(r);
                    countIn++;
                }
                else if (line != "") {  // invalid line, spit out
                    std::cerr << "Error in data input file, "
                        << "ignored line " << countLines
                        << ".  Data ignored is:  "
                        << line << std::endl;
                    // better make this output to an input error
                    //log txt file?
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

	
    // method to get all rvectors from a container of rvectors
    // the container data is assumed to be empty but
    // should work even if it is not since we are only pushing back
    // provided data dimensions match existing data if any
    size_t subpavings::getRvectorsFromRVec(RVecData& data,
                                const RVecData& rvec)
    {

        size_t retValue = 0;
        bool cancontinue = true;

        if (!data.empty()) { // data already in data, print a warning
            std::cout
                << "Warning: adding to existing data "
                << "- mixing datasets"
                << std::endl;

            // assume we will be adding all the data from rss.Samples
            data.reserve(rvec.size() +
                                data.size());
        }

        else {
            data.reserve(rvec.size());
        }

        if (cancontinue && !data.empty()) { // check dimensions
            //find the data dimensions from the first datapoint
            int dim = Ub(*rvec.begin()) - Lb(*rvec.begin()) + 1;
            //find the data dimensions from the existing data
            int dataDim = Ub(*data.begin()) - Lb(*data.begin()) + 1;
            if (dim != dataDim) {
                std::cout
                    << "Existing data different dimension "
                    << " to data in container: aborting insertion"
                    << std::endl;
                cancontinue = false;
            }
        }

        // if cancontinue is still true we know there is data we can take
        if (cancontinue) {

            data.insert(data.end(), rvec.begin(),rvec.end());

            retValue = rvec.size();    // data read in
        }

        return retValue;
    }

    // method to get a sample of rvectors from a container
    // the container data is assumed to be empty but
    // should work even if it is not since we are only pushing back
    // provided data dimensions match existing data if any
    // takes a set-up rgsl
    size_t subpavings::getSampleRvectorsFromRVec(RVecData& data,
                gsl_rng * rgsl, size_t samplesize, const RVecData& rvec)
    {
        getSampleFromContainer(samplesize, rgsl, rvec, data);

        return data.size();
    }


    // method to get a sample of rvectors from an RSSample object
    // the container data is assumed to be empty but
    // should work even if it is not since we are only pushing back
    // provided data dimensions match existing data if any
    // takes a set-up rgsl
    size_t subpavings::getSampleRvectorsFromRSSample(RVecData& data,
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
size_t subpavings::getRvectorsFromRSSample(RVecData& data,
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
            << " with label " << label << "- aborting "
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
    size_t subpavings::getRvectorsFromRSSampleForSampling(RVecData& allData,
                        RVecData& sampleData, size_t samplesize,
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
                gsl_rng * rgsl, const RVecData& allData,
                RVecData& sampleData)
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

	//gloria's additions
	// Point mass filtering from data in a RVecData container.
	void subpavings::pointMassFilter(RVecData& theData, std::map<rvector, size_t, std::less<rvector> > & CountsMap)
	{
		RVecDataItr dataIt;
		std::pair<map<rvector, size_t, less<rvector> >::iterator, bool> boolCount;
	   //go through the dataset
		for (dataIt = theData.begin(); dataIt < theData.end(); dataIt++) {
			//insert into map
			boolCount = CountsMap.insert(make_pair(*dataIt,1));		
			//Check if insertion is successful - a unique set will render a successful insertion.
			if(!(boolCount.second)) { //if there are repeated points
		      //increment count for this key
		      CountsMap[(*dataIt)] += 1;	
			}  
		} // end of going through each data in this dataset
	}
	
	// Labels an RVecData object and store as an RSSample object.
	void subpavings::labelDataFromFilter(RVecData& theData, RSSample&  labData,
				  std::map<rvector, size_t, std::less<rvector> > & CountsMap,
				  std::map<rvector, double, std::less<rvector> > & EMFMap)
	{
		RVecDataItr dataIt;
		for (dataIt = theData.begin(); dataIt < theData.end(); dataIt++){
			if (CountsMap[(*dataIt)] > 1 ) { //for non-unique points
			//if (CountsMap[(*dataIt)] > int(theData.size()*0.01) ) { //for loosely non-unique points
				// insert into RSSample object with label 0
				LabPnt data0;
				data0.Pnt = *dataIt;
				data0.L = 0;
           // cout << "non-unique point: "; 
				//data0.Print(cout);
				labData.Samples.push_back(data0);
				//Also store the EMF of (*dataIt)
                                //TODO: this is crazy GT?? Every time you insert multiple pointmass events at the same point you update the same thing???this EMFMap should be removed and weights obtained from CountsMap directly (perhaps replace size_t with double for value in CountsMap map)
				EMFMap[(*dataIt)] = (double(CountsMap[(*dataIt)])*1.0)/(double(theData.size())*1.0);
			}
			else {
				// insert into RSSample object with label 1
				LabPnt data1;
				data1.Pnt = *dataIt;
				data1.L = 1;
				//cout << "unique point: ";
				//data1.Print(cout);
				labData.Samples.push_back(data1);				
			}			
		} // end of iterating through data
	}

//--src_trunk_0701

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
//--src_trunk_0701


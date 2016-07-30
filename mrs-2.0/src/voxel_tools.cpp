/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief Implementation of voxel tools functions
*/


#include "voxel_tools.hpp"
//#include "sptools.hpp"
//#include "toolz.hpp"

#include <fstream>
#include <sstream>
#include <iterator>
#include <stdexcept>
#include <cassert>

using namespace std;
using namespace subpavings;

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

bool subpavings::findBlack(string line, string seek)
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
	int x_dim = 0;  // x-dimensions
	int y_dim = 0;  // y-dimensions
	int z_dim = 0;  // z-dimensions

	size_t expectedDims = 3;     // expected dimensions
	string seek = "255"; // the string indicating a voxel in the image

	bool retValue = false;
	bool cancontinue = true;

	//set up the file and read input line by line
	// we need to convert the string argument to a c-string for ifstream
	ifstream dataFile(s.c_str());
	
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



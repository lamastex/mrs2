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

#ifndef ___SPTOOLS_HPP__
#define ___SPTOOLS_HPP__

#include "sptypes.hpp"
#include <gsl/gsl_rng.h>    // to use the gsl random number generator

/*! \file sptools.hpp
\brief General tools for subpavings

*/

//! Forward class declarations
    class RSSample;

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    /*! \brief Construct a unique filename from base and timestamp number.

    Waits until timestamp number combined with baseFileName gives
    a name for a file that does not yet exist.

    \param baseFileName the base to use as the filename
    \param suffix to put on filename, defaults to ".txt"
    \return a unique filename with baseFileName suffix, timestamp number prefix
    */
    string getUniqueFilename(string baseFileName,
                string suffix=".txt");

    /*! \brief Append a line to a txt file.

    \param s the name of the txt file to send output to.
    \param line the line to append to the file.
    \param append controls whether line is appended to file, defaults to true.
    */
    void outputFile(const string& s, const string line, bool append = true);


    /*! \brief Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.
    */
    void outputFileStart(const string& s);

    /*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param vals is a vector of values to add to the file
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process
    */
    void outputFile(const std::string& s, RealVec& vals, int i);

    
    /*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string intro,
                    RealVec& vals);

	
	/*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string intro,
					IntVec& vals);

    /*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string intro,
					RVecData& vals);

    /*! \brief Append a tab delimited list seq of strings as a line in log file.
    \param s the name of the txt file to send output to.
    \param strings is a vector of strings to add to the file
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process
    */
    void outputFile(const std::string& s, std::vector<string>& strings, int i);



    /*! \brief Parse a string to make lines for a dot graph.

    \param s the name of the file to output to
    \param toParse the string to parse
    \return true if parsing worked, false if problems encountered
    */
    bool parseForGraphDot(string s, string toParse);

    /*! \brief make a Dot graph png image given a dot file.

    Prints out a line to console with name of image file

    \param s the name of the .dot file to make the graph out of.
    */
    void makeDotImage(string s);


    /*! \brief Method to count lines in a txt file.

    Counts all lines with no data checking.  Used for counting lines before
    setting up an adaptive histogram.

    \param s is the name of the txt file in which to count the lines.
    \return the number of lines in the file.
    */
    int countLinesInTxt(const string& s);



    /*! \brief parse a .vtk header line for spacings data.

    Expects a line of the form "DIMENSIONS 32 32 32"
    \param line is the line to be parsed.
    \param spacings is an empty container for the spacings.
    \return the container of spacings filled up from the data in the line.
    */
    IntVec parseSpacings(string line, IntVec& spacings);

    /*! \brief Get coordinates from a .vtk file.

    Expects structured point format data.
    @param Xs a container of integers to put x-coordinate data into.
    @param Ys a container of integers to put y-coordinate data into.
    @param Zs a container of integers to put z-coordinate data into.
    @param s the filename to get the data from.
    @pre empty containers for Xs, Ys, Zs and a properly formated file s.
    @post containers Xs, Ys, Zs filled with coordinate data.
    @return a vector of the coordinate spacing in the x, y, z direction.
    */
    IntVec getCoordinatesFromVtk(IntVec& Xs, IntVec& Ys, IntVec& Zs,
                                            const string& s);


    /*! \brief A function for comparing ivectors based on volume.

    Used in sorting a list of ivectors ordered by volume.

    Uses Volume() as defined in toolz.
    \return TRUE if the Volume of a is strictly less than Volume of b.
    */
    bool volCompare(const ivector &a, const ivector &b);


    /*! \brief Work arround for c-xsc math library for exponentiation of reals.

    The library handling exponentiation of reals cannot deal with large reals.
    This function tries the exponentiation after conversion to a double and
    converting the results back to a real if success.

    */
    real tryExp(real r);

    /*! \brief Make a double into an rvector.

    Used for processing input to an AdaptiveHistogram object.
    */
    rvector makeDoubleIntoRvector (const double d);

    /*! \brief A quick check on a data string: expecting only numbers
    white space or decimal points.

    Checks for illegal characters and checks number of decimal points.
    Used for checking txt file input.  For example, a string
    "12.04 1.00005e-10 -30.0006" contains no illegal characters and 3
    decimal points and would pass the test if the value of n is 3.
    For example, a string "12.04ab 1.00005e-10 -30" only has 2 decimal
    points ('-30' is not a valid format since it does not have its decimal
    point) and the string also has illegal characters 'a' and 'b'.
    Allowing 'e' could lead to problems but is necessary to be be able to
    accept numbers in floating point format.

    \param s a reference to the string to check.
    \param n the number of decimal points to expect, which is taken to
    indicate how many blocks of numbers the string should contain since
    each number block is expected to contain a decimal point.
    \return true if the string passed the test, false otherwise.
    */
    bool checkString(const string& s, const int n);

    /*! \brief Find the number of 'blocks' of numbers in a properly
    formatted string of numbers.

    Properly formatted means each number includes a decimal point.
    Blocks of numbers are assessed by looking for numerical characters
    including 'e', '+', '-' and also looking for decimal points '.'.
    Used for checking txt file input.  For example, a
    string "12.04 1.00005e-10 -30.0006" contains 3 blocks of numbers.
    Allowing 'e' could lead to problems but is necessary to be be able to
    accept numbers in floating point format.
    Returns 0 if not all numbers have a decimal point or if the string
    contains anything not in ".123456789 \t"

    \param s a reference to the string to check.
    \return the number of 'blocks' of numbers found.
    */
    int countNumbers(const string& s);


    /*! \brief Read in one-dimensional data from a txt file.

    Reads in lines of data representing integers or doubles from a txt file.
    Expects one line per integer or double and will ignore any further numbers
    on a line.
    Will ignore lines where first 'word' (block of characters terminated by
    white space or endofline character) found cannot be converted to an
    integer or double, and print a message to console about ignored line.

    For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    12.04.  A string "30 abc 12.04.0006" will be read as  30, and a string
    "abc 30 12.04.0006" will be rejected.

    The container theData is assumed to be empty but this will work even if
    it is not empty provided that the data already there is also 1-dimensional.

    \param theData is a reference to a container of data into which to
    place the data once read.
    \param s is the name of the txt file from which to read the data.
    \param headerlines is number of headerlines to skip before reading 
	data.  Defaults to 0.
    \pre A container theData, assumed to be empty, in which to store data.
    theData should support push_back() method.
    \pre A file with filename s in the same directory as main() or with the
    filename incorporating directory location.
    \post The container theData contains the valid rvectors read in.
    \return true if at least some data was read in, false if the file could
    not be opened, if the file contained no data, or if no valid data was
    found in the file.
    */
    bool readOneDimDataFromTxt(RVecData& theData, 
								const string& s,
								const std::size_t headerlines = 0);


    /*! \brief Read in rvectors from a txt file.

    Reads in lines of data representing rvectors from a txt file.
    The dimensions of the rvector are deduced from the input format and all
    the data then is expected to have the same dimension.  Any data not
    matching the expected dimensions, based on assessing the first valid
    line found, will be rejected.  Expects one line per rvector with the
    elements separated by white space (space or tabs), with no non-numeric
    characters.  Carries out basic data checking through checkString().
    Input lines which do not pass are printed to standard output with an
    error message but the entire file will continue to be processed and
    valid lines converted to rvectors which are stored in theData.
    Conversion to rvectors is via the cxsc::operator<< which allows an
    rvector to be constructed from a stream. Can read 1-d rvector data from
    doubles or floats but insists on a decimal point in each number
    (otherwise the number is rejected).

    For example, a string "12.04 1.00005e-10 -30.0006" will be read as a
    3-dimensional rvector, a string "-30.0006" will be read as a
    1-dimensional rvector, and a string "30" will be rejected.

    The container theData is assumed to be empty but this will work even if
    it is not empty provided that the data added has the same dimension as
    the existing data, since we are only pushing back: new data will just
    get put in after the stuff already there.

    \param theData is a reference to a container of data into which to
    place the data once read.
    \param s is the name of the txt file from which to read the data.
    \param headerlines is number of headerlines to skip before reading 
	data.  Defaults to 0.
    \pre A container theData, assumed to be empty, in which to store data.
    theData should support push_back() method.
    \pre A file with filename s in the same directory as main() or with the
    filename incorporating directory location.
    \post The container theData contains the valid rvectors read in.
    \return true if at least some data was read in, false if the file could
    not be opened, if the file contained no data, or if no valid data was
    found in the file.
    */
    bool readRvectorsFromTxt(RVecData& theData, 
								const string& s,
								const std::size_t headerlines = 0);
								
   
    /*! \brief Get all rvectors from a container of rvectors

    The container theData is assumed to be empty but this will work even if
    it is not empty provided that the data added has the same dimension as
    the existing data, since we are only pushing back: new data will just
    get put in after the stuff already there.

    \param data is a reference to a container of data into which to
    place the data once read.
    \param rvec is the container of rvectors with the data we want.
    \return the number of data points read in from the container rvec,
    which should be the number of data points in rvec provided that
    they match the dimension of data points already in data if any.
    Returns 0 if dimensions did not match existing data points in data.
    */
    size_t getRvectorsFromRVec(RVecData& data, const RVecData& rvec);

    /*! \brief Get a sample of rvectors from an a container.

    Samples with replacement to find required sample size.

    \param data is a reference to a container of data into which to
    place the data once read.
    \param rgsl is a pointer to a random number generator
    \param samplesize is the size of the required sample
    \param rvec is the container with the data we want to sample from.
    \return the number of data points taken as a sample from the container.
    Returns 0 if dimensions did not match existing data points in data.
    */
    size_t getSampleRvectorsFromRVec(RVecData& data, gsl_rng * rgsl,
                size_t samplesize, const RVecData& rvec);


    /*! \brief Get a sample of rvectors from an RSSample object.

    A function to get a sample of rvectors from the
    labeled point sample data held by an RSSample object

    Uses getRvectorsFromSample to get something to sample from.  Then
    samples with replacement to find required sample size

    \param data is a reference to a container of data into which to
    place the data once read.
    \param rgsl is a pointer to a random number generator
    \param samplesize is the size of the required sample
    \param rss is the RSSample object with the data we want.
    \param label is the label for the points that we want from rss.
    \pre A container theData, assumed to be empty, in which to store data
    sampled from rss. theData should support push_back() method.
    \post The container theData contains a sample of rvectors from rss.
    \return the number of data points taken as a sample from rss.Samples.
    Returns 0 if no data points with the correct label found in rss.Samples
    or if dimensions did not match existing data points in data.
    */
    size_t getSampleRvectorsFromRSSample(RVecData& data, gsl_rng * rgsl,
                    size_t samplesize, const RSSample& rss, int label);


    /*! \brief Get all rvectors from an RSSample object.

    Takes the rvector part of labeled points in the Samples data member
    of an RSS object, where the label on the point matches the label
    supplied as a function argument.

    RSSample objects can have Samples containing points with different
    labels and points with different labels can have different dimensions.

    The container theData is assumed to be empty but this will work even if
    it is not empty provided that the data added has the same dimension as
    the existing data, since we are only pushing back: new data will just
    get put in after the stuff already there.

    \param data is a reference to a container of data into which to
    place the data once read.
    \param rss is the RSSample object with the data we want.
    \param label is the label for the points that we want from rss.
    \pre A container theData, assumed to be empty, in which to store data.
    theData should support push_back() method.
    \post The container theData contains the valid rvectors read in.
    \return the number of data points matching label, and matching dimensions
    data points already in data if any, found in rss.Samples.  Returns 0 if
    no data points with the correct label found or if dimensions did not
    match existing data points in data.
    */
    size_t getRvectorsFromRSSample(RVecData& data, const RSSample& rss,
                                int label);



    /*! \brief Get data from an RSSample object to take samples from

    \param allData is the container to put all the data found into.
    \param sampleData is the container which the sample wll eventually go into.
    \param samplesize is the size of the sampel to be drawn.
    \param rss is the RSSample object to sample data from.
    \param label is the label for the labeled points we want to find in rss.
    \return the number of data points taken as a sample from rss.Samples.
    Returns 0 if no data points with the correct label found in rss.Samples
    or if dimensions did not match existing data points in data.
    */
    size_t getRvectorsFromRSSampleForSampling(RVecData& allData,
                    RVecData& sampleData, size_t samplesize, const RSSample& rss,
                    int label);


    /*! \brief Get a sample from data in a container.

    A function to get a random sample of of a specified size from a given
    container of data.  Sampling is with replacement.

    \param samplesize is the size of the sample to get.
    \param rgsl is a pointer to a random number generator.
    \param allData is the container of data points from which to draw
    the sample.
    \param sampleData is a container into which to put the sample drawn

    \post sampleData has samplesize data points from allData, drawn with
    replacement.
    */
    void getSampleFromContainer(size_t samplesize, gsl_rng * rgsl,
                const RVecData& allData, RVecData& sampleData);
	 

	//gloria's additions
	/*! \brief Point mass filtering from data in a RVecData container.
	Takes an RVecData container and sieves through the data to obtain "point mass"
	data, i.e. data that falls on the same location with precision up to
	instrument settings or user's definition.
	
	\param theData is the container to put all the data found into.
	\param CountsMap is passed in by reference and stores the empirical mass 
	function of the data.
	*/
	void pointMassFilter(RVecData& theData, 
						std::map<rvector, size_t, std::less<rvector> > & CountsMap);
							
	/*! \brief Labels an RVecData object and store as an RSSample object.
	Takes an RVecData container and its corresponding map to label the points.
	Point mass data gets label 0 and the rest gets label 1. The labelled data is
	then stored in an RSSample object.
	
	The EMF of the point mass data is also returned and stored in a map.
	
	\param theData is the container to put all the data found into.
	\param CountsMap is passed in by reference and stores the empirical mass 
	function of the data.
	*/
	void labelDataFromFilter(RVecData& theData, RSSample& labData, 
						std::map<rvector, size_t, std::less<rvector> > & CountsMap,
						std::map<rvector, double, std::less<rvector> > & EMFMap);
	
	
	//src_trunk_0701
	
	/*! @name String representations of collections.
	 * 
	 Formatted to have ", " or "," (\a compact = true) between the 
	 elements of the collection.
	 
	 * \param vec the collection to describe in the output.
	 * \param compact indicator for whether format should be compact.
	 * \return string representation of \a vecints.*/
	//@{
	/*! \brief String representation of an IntVec. */
	std::string toString(const subpavings::IntVec vec,
												bool compact = false);


	/*! \brief String representation of a RealVec.*/
	std::string toString(const subpavings::RealVec vec,
												bool compact = false);
	
	/*! \brief String representation of a VecDbl.*/
	std::string toString(const subpavings::VecDbl vec,
												bool compact = false);
	//@}

	
} // end namespace subpavings

#endif

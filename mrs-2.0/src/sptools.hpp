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



/*! \file
\brief General tools for subpavings

*/

#ifndef ___SPTOOLS_HPP__
#define ___SPTOOLS_HPP__

#include "sptypes.hpp"
#include "SmallClasses.hpp"

#include <gsl/gsl_rng.h>    // to use the gsl random number generator

#include <map>
#include <utility>

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {


    /*! \brief Construct a unique filename from base and timestamp number.

    Waits until timestamp number combined with baseFileName gives
    a name for a file that is unique.
	
	\param baseFileName the base to use as the filename
    \param suffix to put on filename, defaults to ".txt"
    \return a unique filename
	constructed as  baseFileName + timestamp number + suffix. */
    std::string getUniqueFilename(const std::string& baseFileName, 
				const std::string& suffix = ".txt");

	/*! \brief Construct a filename from base and process ID (PID).

    This does not necessarily give a unique filename, but may be more
	helpful than getUniqueFilename in matching output files to process 
	provided that only one output file with base filename 
	\a baseFileName will be generated in one process. 
	
    \param baseFileName the base to use as the filename
    \param suffix to put on filename, defaults to ".txt"
     \return afilename
	constructed as  baseFileName + PID + suffix. */
    std::string getPidFilename(const std::string& baseFileName, 
				const std::string& suffix = ".txt");

    /*! \brief Append a line to a txt file.
	
	Puts a new line before and after the \a line. 

    \param s the name of the txt file to send output to.
    \param line the line to append to the file.
    \param append controls whether line is appended to file, defaults to true.
    */
    void outputFile(const std::string& s, const std::string& line, bool append = true);

	/*! \brief Append a string to a txt file.
	
	Does not put in any new line characters.

    \param s the name of the txt file to send output to.
    \param line the string to append to the file.
    \param append controls whether line is appended to file, defaults to true.
    */
    void outputFileString(const std::string& s, const std::string& line, bool append = true);

    /*! \brief Opening line of a txt log file.

    Starts the log file with file name and date and time.
    \param s the name of the txt file to send output to.    */
    void outputFileStart(const std::string& s);
    
    /*! \brief Strip the path from a filename.

    If the filename is in the form "path/file" then just the
    string file is returned.  If the filename is in the form file
    (ie there is no path) then the string returned is file.
    
    \note The path and file are assumed to be delimited by "/".
    
    \param filename the name of the file to strip path from.
    \return a string comprising \a filename without path.*/
    std::string stripPath(const std::string& filename);

    /*! \brief Append a tab delimited seq of values as a line in log file.
	 
	This routine is intended to help with logging a process like MCMC.
	
    \param s the name of the txt file to send output to.
    \param vals is a vector of values to add to the file.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.    */
    void outputFile(const std::string& s, subpavings::RealVec& vals, int i);

	/*! \brief Append a tab delimited seq of values as a line in log file.
     
	This routine is intended to help with logging a process like MCMC.
	
    \param s the name of the txt file to send output to.
    \param vals is a vector of values to add to the file.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
	\param nodeType an indictator for the type of node (eg leaf, cherry).
	\param accepted an indicator for whether the move was accepted.
	\param nCount number of points in the node concerned.    */
    void outputFile(const std::string& s, RealVec& vals, int i,
			int nodeType, int accepted, size_t nCount);

	/*! \brief Append a value to a log file.
     
	\param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
    \param val a value to log.*/
    void outputFile(const std::string& s, int i, int val);
   
	/*! \brief Append a value to a log file.
     
	\param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
    \param val a value to log.*/
    void outputFile(const std::string& s, int i, size_t val);
   
    
    /*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string& intro,
                    subpavings::RealVec& vals);

	
	/*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string& intro,
					subpavings::IntVec& vals);

    /*! \brief Append a tab delimited list seq of values as a line in log file.
    \param s the name of the txt file to send output to.
    \param intro is the introduction before the rest is put in, all on one line.
    \param vals is a vector of values to add to the file
    */
    void outputFile(const std::string& s, const std::string& intro,
					subpavings::RVecData& vals);

    /*! \brief Append a tab delimited list seq of strings as a line in log file.
    \param s the name of the txt file to send output to.
    \param strings is a vector of strings to add to the file
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process
    */
    void outputFile(const std::string& s, 
					std::vector< std::string >& strings, int i);



    /*! \brief Parse a string to make lines for a dot graph.

    \param s the name of the file to output to
    \param toParse the string to parse
    \return true if parsing worked, false if problems encountered
    */
    bool parseForGraphDot(const std::string& s, std::string toParse);

    /*! \brief make a Dot graph png image given a dot file.

    Prints out a line to console with name of image file

    \param s the name of the .dot file to make the graph out of.
    */
    void makeDotImage(const std::string& s);


    /*! \brief Method to count lines of data in a txt file.

    Counts all lines with no data checking.  Used for counting lines before
    setting up an adaptive histogram.

    \param s is the name of the txt file in which to count the lines.
	\param headerlines is the number of lines to ignore.
    \return the number of lines of data, excluding \a headerlines, 
	in the file.
    */
    int countLinesInTxt(const std::string& s, const size_t headerlines = 0 );
	
	/*! @name Make a box to contain all the data in a collection.
	
	The box is made to fit all of the data in \a theData, with padding
	\a padding around the minimums and maximums of the data on each 
	dimension.

    \param theData a reference to a container of rvector data.
    \param padding the padding to put around the minimums and maximums
	of the data in \a theData.
    \return an ivector same dimensions as the data and to fit all the data
    including an allowance for padding.    
	\pre \a theData is not empty.*/
	//@{
    cxsc::ivector makeBox(const RVecData& theData, cxsc::real padding);
	
	cxsc::ivector makeBox(
					const std::vector < std::vector <double> >& theData,
					cxsc::real padding);
	//@}
	/*! \brief Make a box to contain all the data in a collection.
	
	The box is made to fit all of the data in \a theData, with padding
	\a padding around the minimums and maximums of the data on each 
	dimension.

    \param theData a reference to a container of rvector data.
    \param dim the dimensions of the data.
	\param padding the padding to put around the minimums and maximums
	of the data in \a theData.
    \return an ivector same dimensions as the data and to fit all the data
    including an allowance for padding.
    \pre \a theData is not empty.
	\pre \a dim is the dimension of all the data in \a theData.*/
    cxsc::ivector makeBox(const RVecData& theData, int dim, 
												cxsc::real padding);

	bool checkBox(const cxsc::ivector& box);

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

	/*! \brief Reformat a line to be suitable for reading into an rvector.
	
	Strips out dos carriage returns.
	
	\param line a reference to a string to reformat.
	\return the same \a line reformatted.	*/
	std::string& reformatLine(std::string& line);

    /*! \brief A quick check on a data string: expecting only numbers
    white space or decimal points.

    Checks for illegal characters;  Used for checking txt file input. 
	
	Legal characters are any characters in the string 
	"eE+-.,0123456789 \\t\\r" 
	
	Note that legal characters includes dos carriage return '\\r'.
	
	For example, a string
    "12.04 1.00005e-10 -30.0006" contains no illegal characters.
	
    For example, a string "12.04ab 1.00005e-10 -30" has illegal 
	characters 'a' and 'b'.
	
    \param s a reference to the string to check.
    \return true if the string passed the test, false otherwise.
    */
    bool checkString(const std::string& s);

    /*! \brief Find the number of 'blocks' of numbers in a 
    string of numbers.

    A number is identified as a block of characters containing at least 
    one of the characters 0123456789
    * 
    Blocks of numbers delimited by space ' ', comma ',' or tab '\\t'
    * 
    This function does not check for non-numeric characters, for example, 
    the string "a12.3 2E+03" will be assessed as containing 2 blocks of
    numbers, "e" will be assessed as containing zero blocks of numbers,
    "1e" will be assessed as containing one block of numbers,
    "." will be assessed as containing zero block of numbers,
    ".0" will be assessed as containing one block of numbers,
    "a.0" will be assessed as containing one block of numbers.
    
    \param s a reference to the string to check.
    \return the number of 'blocks' of numbers found.
    */
    int countNumbers(const std::string& s);

	/*! \brief Log an error from reading in data.
	 * 
	 * At present this just spits the error message out to cerr.
	 * 
	 * \param line is the line that caused the problem.
	 * \param lineNumber is the number of the line in the input file.
	 * */							
	void lineErrorLogger(const std::string& line, int lineNumber);

	//new
	/*! \brief Read an rvector in from a string.
	 * 
	 * \param r is a reference to the rvector to read into.
	 * \param line is the line to read into \a r.
	 * \return the reference to \a r, with \a line read in. 
	 * \pre The length of \a r should be the same as the
	 * number of values in \a line.*/
	cxsc::rvector& readRV(cxsc::rvector& r, const std::string& line);

	//new
	/*! \brief Read specified elements from a string into an rvector.
	 * 
	 * The elements are specified in \a reqDims.
	 * 
	 * Eg reqDims = 5, 3, 7 would create a 3-d rvector r where 
	 * r[1] is the 5th value in the line, r[2] is the 3rd value in the 
	 * line, r[3] is the 7th value in the line.   
	 * 
	 * \param r is a reference to the rvector to read into.
	 * \param line is the line to read into \a r.
	 * \param reqDims is a container of dimensions to read from \a line.
	 * \param lineDim is the number of values the line is meant to contain.
	 * \return the reference to \a r, with \a line read in.
	 * \pre \a line must have at least as many separate values in 
	 * as the maximum value in \a reqDims.  The length of \a r
	 * should be the same as the size of reqDims. 
	 * \a line contains \a lineDim values.*/
	cxsc::rvector& readRV(cxsc::rvector& r, const std::string& line,
									const std::vector < int >& reqDims,
									int lineDim);


	/*! @name Routines to read in rvectors from a txt file and 
	 * store in a given data container \a theData.

    Reads in lines of data representing rvectors from a txt file.
    The dimensions of the data may given or may be deduced 
	from existing data in \a theData or from the data in the file.
	
	Expects one line per rvector with the elements separated by 
	white space (space or tabs), or by commas.
	
	Input files can be checked to try to ensure data conforms with 
	some given dimensions and/or has no illegal characters.  
	* 
	Illegal characters are any characters other than those in the string
	"eE+-.,0123456789 \\t".
	* 
	Note that this allows numbers in scientific format, and with space,
	tab, or comma delimiters.  However, allowing the characters 'e'
	and 'E' as part of scientific format, and the delimiters,
	also allows the possibility that these are erroneously present.  
	The level of checking for illegal characters is not rigorous and will
	not pick up the inadvertent or erroneous presence of these non-numeric
	characters.  This includes the presence of delimiters before the first 
	numeric character, which may result in 0.0 being inserted into the
	first dimension(s) of the rvector read from the relevant data line.
	* 
	<strong>Ordinary checking</strong> means that:<ul>
	<li>The first valid line is checked against the 
	given data dimensions if any.  Oherwise the first valid data line
	is used to find data dimensions.</li>
	<li> If data dimensions are given and the dimensions found in the 
	first valid data line are less than those given dimensions, an error
	message is given and reading of all data is aborted. Any data 
	in excess of the given dimensions will however just be ignored.    
	<li> The first valid data line, as described above, is the 
	first line after the given headers not containing illegal characters.
	All lines after that are also checked to ensure that they do not 
	contain illegal characters.    
	</ul>
	Note that under ordinary checking,   data lines after the first 
	valid line are not checked to ensure that they conform to 
	expected dimensions: data in excess of the given dimensions will 
	be ignored and 'missing' dimensions will be padded with 0.0's.    
	 
	<strong>Paranoid checking</strong> means that expected 
	data dimensions must be given and checking is carried out as for 
	Ordinary checking plus further checking:<ul>
	<li> All lines after the first valid data line are also checked to
	to check that the number of blocks of numbers is equal to the
	expected dimensions.  </li> 
	</ul>
	* 
	<strong>Fast checking</strong> effectively means no checking (fast
	reading of data).  Data is assumed to conform to the given dimensions:
	data in excess of the given dimensions will be ignored and
	all data lines which are less than the given dimensions will result in
	an rvector "padded" with 0.0 in the additional dimensions.  
	Illegal characters in the input, including initial delimiters before
	the first numeric values, may be interpreted as part of a
	real vector, or as delimiters, or ignored.  No guarantees are given 
	about the resulting data read in if the data does not conform to the
	values given for the parameters for headerlines and data dimensions. 
	* 
	While reading of the file continues (has not been aborted, for 
	example because the dimension of the data does not seem to be as 
	expected) , input lines which do not pass checks are rejected and logged
	but the rest of the file will continue to be processed and
    valid lines converted to rvectors to be stored in \a theData.
	
    Blank lines are ignored.

    These routines will work if \a theData it is not empty, provided that 
    if \a dim is given, this matches the dimensions of the existing data.

	If no data dimension is specified by the user (\a dim = -1), 
	then the routine will attempt to find it.  If there is existing data
	in \a theData, the expected dimension \a dim will be taken as the 
	dimension of the first element in \a theData.  If \a theData is empty,
	the expected dimension \a dim will be taken as the dimension of the
	first valid data line read in from the file.

	Throws an invalid_argument exception if \a dim is supplied by the
	user and \a dim < 1.
	
	Throws an invalid_argument exception if there is existing data in
	\a theData and the dimensions of this data are < 1.   

    \param theData is a reference to a container of data into which to
    place the data once read.
    \param s is the name of the txt file from which to read the data.
    \param headerlines is number of headerlines to skip before reading 
	data.
	\param dim expected data dimension.  
	
   \return true if at least some data was read in, else false (if
    the file could not be opened, if the file contained no data, 
    or if no valid data was found in the file).
    \pre A file with filename s (s can include path and name), if 
	\a dim is supplied then \a dim >= 1, and if there is existing
	data in \a theData, the dimensions of that data are >= 1.
    \post The container theData contains all data that it contained 
	before the operation and the valid rvectors read in.
    \todo Even paranoid checking does not check that the 
	text can be properly represented in floating point format.
	See "C++ Toolbox for Verified Computing" book pp. 46-50.*/
	
	//@{

	/*! \brief Read in rvectors from a txt file with Ordinary level
	 checking, with data dimensions given. */
	bool readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim);
	
	/*! \brief Read in rvectors from a txt file with Ordinary level
	 checking.
	 *  
	 * Data dimensions are deduced from the first valid input line 
	 * read after headerlines.  
	 */						
	bool readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines = 0);

	/*! \brief Read in rvectors from a txt file with Ordinary level
	 checking.
	 *  
	 * Data dimensions are deduced from the first valid input line 
	 * read after headerlines.  
	 */						
	bool readRvectorsFromTxtOrd(subpavings::RVecData& theData, 
							const std::string& s,
							const std::vector < int > & reqDims,
							std::size_t headerlines = 0);

	/*! \brief Read in rvectors from a txt file with Paranoid level
	 checking.*/ 
	bool readRvectorsFromTxtParanoid(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim);

	/*! \brief Read in rvectors from a txt file with Paranoid level
	 checking.*/ 
	bool readRvectorsFromTxtParanoid(subpavings::RVecData& theData, 
							const std::string& s,
							const std::vector < int > & reqDims,
							std::size_t headerlines);

	/*! \brief Read in rvectors from a txt file with Fast level
	 checking. */
	bool readRvectorsFromTxtFast(subpavings::RVecData& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dim);
	
	//@}
	
	/*! \brief Get rvectors from text.
	
	\note Use getRvectorsFromTxtOrd, or
	getRvectorsFromTxtParanoid, or getRvectorsFromTxtFast, rather than this.
	
	\param theData is a reference to a container of data into which to
    place the data once read.
    \param s is the name of the txt file from which to read the data.
    \param headerlines is number of headerlines to skip before reading 
	data.
	\param checkLevel the level of checking to be carried out.
	\param reqDims is a container specifying the positions of the 
	required values in the line in the file.  If \a reqDims is empty 
	and \a dim is not supplied, it will be assumed that all values from
	a line are required to make one rvector.  If \a reqDims is empty 
	and \a dim is supplied, it will be assumed that the first \a dim
	values from a line of the file are required to make one rvector.
	\param dim expected data dimension.  
	\return true if at least some data was read in.
	*/							
	bool _readRvectorsFromTxt(subpavings::RVecData& theData, 
								const std::string& s,
								std::size_t headerlines,
								int checkLevel,
								const std::vector <int>& reqDims,
								int dim = -1);
	
								

	/*! @name Read in rvectors from a vector of vectors of doubles.

    Reads in vectors of doubles, representing rvectors, from a std::vector.
    The dimensions of the rvector can be given or (if \a dim is not given) 
	deduced from existing data in \a theData or from the 
	first vector in the input vector \a inputData.
   
    The rest of the data then is expected to have the same dimension as
	that given or found from the frst vector in the input vector.  
	Any data with dimensions less than the expected dimensions will be 
	rejected.  If data has more than the expected dimensions, only the
	values on the first expected dimensions will be read in (the rest 
	will be ignored).
	   
    Expects one inner vector per rvector.

    If the container theData is not empty this will work 
	provided that the expected data dimension is the same 
	as the existing data has the same dimension as
    the existing data (if \a dim is not specified by the user and there
	is no existing data, then the expected data dimension will have
	been found from the first element in the input vector \a inputData) 

	Throws an invalid_argument exception if \a dim is supplied by the
	user and \a dim < 1.
	
	Throws an invalid_argument exception if there is existing data in
	\a theData and the dimensions of this data are < 1.   

    \param theData is a reference to a container of data into which to
    place the data once read.
    \param inputData is the vector of vectors of doubles to be input.
    \param dim is the dimensions of the data to expect.  If \a dim = -1 then
	the method will try to establish the dimensions from the first 
	data point in the \a inputData.
    \return true if at least some data was put into theData,
    false otherwise.
	\pre If \a dim is supplied then \a dim >= 1, and if there is existing
	data in \a theData, the dimensions of that data are >= 1.
    \post The container theData contains all data that it contained 
	before the operation and the valid rvectors read in.
    */
    //@{
	/*! \brief Version with no expected dimensions specified.*/	
	bool getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
				const std::vector < subpavings::VecDbl > & inputData);
    
	/*! \brief Version with expected dimensions specified by user.*/	
	bool getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
				const std::vector < subpavings::VecDbl > & inputData,
				int dim);
		
	/*! \brief Internal version:  use getRvectorsFromVectorOfVecDbl rather than this.*/							
	bool _getRvectorsFromVectorOfVecDbl(subpavings::RVecData& theData, 
								const std::vector < subpavings::VecDbl > & inputData,
								int dim = -1);
    //@}
   
    /*! \brief Get rvectors from a container of rvectors.
	
	If \a checkDims == false, appends the entire contents of
	\a rvec to the end of \a data unless the first element in \a rvec
	is an empty (0-d) rvector.
	
	If \a checkDims == true, checks the dimensions of each 
	element in \a inputData and
	only add to \a theData those elements where with dimension equal to
	that of the first element in \a inputData .
	
	Aborts with a message to standard error output if there 
	is existing data in \a theData and the dimensions of the 
	first element of \a theData do not match the dimensions of the first
	element of \a inputData.
		
	Throws an illegal_argument exception if the first element
	in \a inputData is an empty rvector.
	
	\param theData is a reference to a container of data into which to
    place the data once read.
    \param inputData is the container of rvectors we want to transfer
	to \a theData.
    \param checkDims is an indicator for whether the dimensions
	each of the elements of \a inputData are compared to the dimensions
	of the first (true) or not (false). 
	\return the number of data points read in from \a inputData.
    Returns 0 if no data points read in from \a inputData.
	\pre  If there is existing
	data in \a theData, the dimensions of that data are >= 1.
    */
    size_t getRvectorsFromRVec(subpavings::RVecData& theData, 
		const subpavings::RVecData& inputData,
		bool checkDims = false);


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
    size_t getSampleRvectorsFromRVec(subpavings::RVecData& data, 
									gsl_rng * rgsl,
									size_t samplesize, 
									const subpavings::RVecData& rvec);


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
    size_t getSampleRvectorsFromRSSample(subpavings::RVecData& data, 
									gsl_rng * rgsl,
									size_t samplesize, 
									const RSSample& rss, int label);


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
    size_t getRvectorsFromRSSample(subpavings::RVecData& data, 
								const RSSample& rss,
                                int label);



    /*! \brief Get data from an RSSample object to take samples from

    \param allData is the container to put the all the data found into.
    \param sampleData is the container which the sample wll eventually go into.
    \param samplesize is the size of the sampel to be drawn.
    \param rss is the RSSample object to sample data from.
    \param label is the label for the labeled points we want to find in rss.
    \return the number of data points taken as a sample from rss.Samples.
    Returns 0 if no data points with the correct label found in rss.Samples
    or if dimensions did not match existing data points in data.
    */
    size_t getRvectorsFromRSSampleForSampling(subpavings::RVecData& allData,
										subpavings::RVecData& sampleData, 
										size_t samplesize, 
										const RSSample& rss,
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
							const subpavings::RVecData& allData, 
							subpavings::RVecData& sampleData);

	//new
	/*! \brief Read from a txt file into a container of containers
	of reals.
	
	Will print error messages to standard error output if the file
	cannot be opened or if no numeric data was found after the
	headerlines.  
	* 
	If \a dims is not given, the expected dimensions will be 
	taken as the number of 
	separate numeric values on the first non-blank line 
	that passes the check checkString() read after
	the headerlines.
	
	Numeric data is expected to be presented as one line per RealVec
	to be read in, with (at least) as many values as the expected
	dimensions (the excess will be silently ignored).
	 
	No checks for illegal characters or bad formatting are carried out: 
	the results of this are unpredictable but will almost certainly
	result in the data read in being incorrect.
	
	Throws an invalid_argument exception if \a dims < 0.
	
	\param theData a reference to a container to put the data in.
	\param s the name of the file to extract from.
	\param headerlines how many lines to skip at the start of the file.
	\param dims the dimensions of the data to expect (defaults to 0).
	\return true if at least one RealVec was added to \a theData,
	false otherwise. */
	bool readVecRealVecFromTxt(std::vector < subpavings::RealVec >& theData, 
							const std::string& s,
							std::size_t headerlines,
							int dims = 0);

	/*! \brief Output multiple containers of data in vertical format.
	 
	Opening line of file gives creation details.
	 
	Puts data from each container in adjacent columns, right to left.
	
	Checks that size of colNames matches size of dataPtrs,
	but does not check that vectors within dataPtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param dataPtrs a container of pointers to containers for the data.
	\param colNames headers for the columns.
	\param filename filename to use for the output.*/
	void outputToFileVertical( std::vector < const subpavings::RealVec* >& dataPtrs, 
							const std::vector < std::string >& colNames,
							const std::string& filename,
							int prec = 5);
	
	/*! \brief Output multiple containers of data in vertical format.
	 
	Opening line of file gives creation details.
	 
	\a lhsLines gives line or state numbers for left hand side column.
	 
	Puts data from each container in adjacent columns, right to left.
	
	\alhsLines must be at least as long as the longest element in 
	\a dataPtrs.  Checks that size of colNames matches size of dataPtrs,
	but does not check that vectors within dataPtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param lhsLInes a container of line or state numbers for left hand
	side column.
	\param dataPtrs a container of pointers to containers for the data.
	\param colNames headers for the columns.
	\param filename filename to use for the output.*/
	void outputToFileVertical( const std::vector < size_t >& lhsLines,
							std::vector < const subpavings::RealVec* >& dataPtrs, 
							const std::vector < std::string >& colNames,
							const std::string& filename,
							int prec);
	
	/*! \brief Output multiple containers of data in vertical format.
	 
	Opening line of file gives creation details.
	 
	Puts data from each container in adjacent columns, right to left.
	
	Checks that size of colNames matches size of dataPtrs,
	but does not check that vectors within dataPtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param dataPtrs a container of pointers to containers for the data.
	\param colNames headers for the columns.
	\param filename filename to use for the output.
	\param prec precision to use for outputting data (cxsc formatting,
	Variable type output).	*/
	void outputToFileVertical( std::vector < const subpavings::Size_tVec* >& dataPtrs, 
							const std::vector < std::string >& colNames,
							const std::string& filename);
	
	/*! \brief Output multiple containers of data in vertical format.
	 
	Opening line of file gives creation details.
	 
	Puts data from each container in adjacent columns, alternating
	between successive columns from \a sizePts and \a dataPtrs,
	right to left.
	
	Checks that size of colNames matches size of collection of ptrs 
	for both size and real, and that size of size data is same as
	size of real data,
	but does not check that vectors within dataPtrs and sizePtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param sizePtrs a container of pointers to containers for the data
	as size_t types.
	\param dataPtrs a container of pointers to containers for the data
	as reals.
	\param colNamesSize headers for the size data columns.
	\param colNamesData headers for the real data columns.
	\param filename filename to use for the output.
	\param prec precision to use for outputting data (cxsc formatting,
	Variable type output).	*/
	void outputToFileVertical( std::vector < const subpavings::Size_tVec* >& sizePtrs,
							std::vector < const subpavings::RealVec* >& realPtrs, 
							const std::vector < std::string >& colNamesSize,
							const std::vector < std::string >& colNamesReal,
							const std::string& filename,
							int prec = 5);

	/*! \brief Add output multiple containers of data in vertical format
	 to an existing file.
	 
	Puts last (see \a startPos) data from each container in adjacent 
	columns, alternating
	between successive columns from \a sizePts and \a dataPtrs,
	right to left.
	
	Does not check that vectors within realPtrs and sizePtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param sizePtrs a container of pointers to containers for the data
	as size_t types.
	\param realPtrs a container of pointers to containers for the data
	as reals.
	\param filename filename to use for the output.
	\param startPos position in collections inside \a sizePtrs, 
	\a realPtrs, at which to start finding data to output.
	\param prec precision to use for outputting data (cxsc formatting,
	Variable type output).	*/
	void outputToFileVertical( 
				std::vector < const subpavings::Size_tVec* >& sizePtrs,
				std::vector < const subpavings::RealVec* >& realPtrs, 
				const std::string& filename,
				size_t startPos, 
				int prec);

	/*! \brief Add output from multiple containers of data in vertical format
	 to an existing file.
	 
	Puts last (see \a startPos) data from each container in adjacent 
	columns, right to left.
	
	Does not check that vectors within realPtrs are all the
	same length - if this is so, the resulting columns will simply
	be of unequal lengths.

	\param realPtrs a container of pointers to containers for the data
	as reals.
	\param filename filename to use for the output.
	\param startPos position in collections inside \a realPtrs,
	at which to start finding data to output.
	\param prec precision to use for outputting data (cxsc formatting,
	Variable type output).	*/
	void outputToFileVertical( 
				std::vector < const subpavings::RealVec* >& realPtrs, 
				const std::string& filename,
				size_t startPos, 
				int prec);

	/*! \brief Add a pointer to a container to a container of pointers.
	
	\param container a reference to a container to add the pointer to.
	\param toAdd a reference to the container for which the pointer 
	is to be added to \a container.
	\return a reference to \a container.*/
	std::vector < const subpavings::RealVec* >& addDataPtrs(
						std::vector < const subpavings::RealVec* >& container,
						const std::vector < subpavings::RealVec >& toAdd);

	/*! \brief Add a pointer to a container to a container of pointers.
	
	\param container a reference to a container to add the pointer to.
	\param toAdd a reference to the container for which the pointer 
	is to be added to \a container.
	\return a reference to \a container.*/
	std::vector < const subpavings::Size_tVec* >& addDataPtrs(
						std::vector < const subpavings::Size_tVec* >& container,
						const std::vector < subpavings::Size_tVec >& toAdd);

	/*! \brief Output the contents of a collection.
	
	Outputs colections as values separated by tabs.
	
	* \param os the stream to send the output to.
	* \param vec the collection to output.
	* \return the stream with the new output.*/
	//@{
	/*! \brief Output the contents of an IntVec. */
	std::ostream & operator<<(std::ostream &os, 
										const subpavings::IntVec& vec);
	
	/*! \brief Output the contents of a RealVec. */
	std::ostream & operator<<(std::ostream &os, 
										const subpavings::RealVec& vec);

	/*! \brief Output the contents of a VecDbl.*/
	std::ostream & operator<<(std::ostream &os, 
										const subpavings::VecDbl& vec);
	//@}
	
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
	
	//new AHABC
	/*! \brief A struct for comparing rvectors.
	 */
	struct RvecComp
	{
		bool operator() (const cxsc::rvector& lhs, const cxsc::rvector& rhs) const;
	};
	
	//new AHABC
	/*! \brief Create a map of counts of points in an rvector.
	 */
	std::map<rvector, size_t, subpavings::RvecComp > & pointMassFilter(
				const subpavings::RVecData& theData,
				std::map<rvector, size_t, subpavings::RvecComp > & countsMap);
	
	/*! \brief Label data from a filter.
	 */
	RSSample& labelDataFromFilter(RSSample&  labData,
		const std::map<rvector, size_t, subpavings::RvecComp >& countsMap,
		int uniqueLabel, int ptMassLabel);
	
	/*! \brief Make an empirical mass filter map.
	 */
	double makeEMFMap(
		std::map < rvector, double, subpavings::RvecComp >& EMFMap,
		const std::map<rvector, size_t, subpavings::RvecComp >& countsMap);
	
	/*! \brief Calculate mean and sample standard deviation from \a vec.
	 */
	std::pair < cxsc::real, cxsc::real > calcMeanAndSD(
							const RealVec& vec);
} // end namespace subpavings

#endif

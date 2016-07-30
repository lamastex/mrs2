/*
* Copyright (C) 2012 Jennifer Harlow
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
\brief FunctionEstimatorReal declarations.
*/

#ifndef ___FUNCTIONESTIMATOR_REAL_HPP__
#define ___FUNCTIONESTIMATOR_REAL_HPP__

#include "piecewise_constant_function.hpp"
#include "realmappedspnode.hpp"
#include "spnode.hpp"
#include "mappedFobj.hpp"


namespace subpavings {


/*! \brief A wrapper or manager for an RealMappedSPnode
tree to be used for function estimation.

*/

class FunctionEstimatorReal {

    public:
	
	/*! \brief Initialised constructor.

    Initialised  with domain box.
    
	Throws a MalconstructedBox_Error if the box is not suitable as the
    basis of a subpaving (eg, box has no dimensions, or the box has
    a thin interval on at least one dimension).

    Ideal constructor when the support domain of the function is
	set a priori.
	\param v The box to use for the subpaving to be managed.
	\param f The function to be estimated by this.
	\param lab The label for this (defaults to 0).    
	\post The estimator constructed has subpaving that consists 
	of single leaf node (the root) with a box like \a v) 
	and the range on that single node is the 
	the real image of \a v under the function represented by \a f.*/
    FunctionEstimatorReal(const ivector& v,
										const MappedFobj& f,
										int lab = 0);

    /*! \brief Initialised constructor.

    Initialised  with a subpaving.
    
	\param spn A subpaving to copy as the subpaving to be managed.
	\param f The function to be estimated by this.
	\param lab The label for this (defaults to 0).    
	\pre \a spn has a box.
	\post The estimator constructed has label \a lab, 
	a subpaving that that is a copy 
	of spn and the range on each box in that subpaving is 
	the real mid-image of that box
	under the function represented by \a f.*/
    FunctionEstimatorReal(const SPnode& spn,
							const MappedFobj& f, 
							int lab = 0);

	/*! \brief  Copy constructor.
    */
    FunctionEstimatorReal(const FunctionEstimatorReal& other);

	
    //! Destructor
    ~FunctionEstimatorReal();

	/*! \brief Get the reference to the function object used by this.
	
	\return the function object reference used by this.*/
    const MappedFobj& getFobjReference() const;

	/*! \brief Get the label.
	
	\return the label for this.*/
    int getLabel() const;

	/*! \brief Set the label.*/
    void setLabel(int lab);

	/*! \brief Get whether this has a subpaving to manage.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return true if this has a subpaving to manage.
	false otherwise.*/
    bool hasSubPaving() const;
	
	/*! \brief Get whether the subpaving managed by this 
	 has negative range values.
	
	\return true if the subpaving managed by this has negative ranges,
	false otherwise.*/
    bool hasNegativeFunctionEstimates() const;
	
	/*! \brief Get the box of the subpaving managed by this.
	
	\note with the present constructors, it is impossible for
	this to have a subpaving but for the subpaving to have no box.

    \return copy of the box of the subpaving managed by this.
	\pre hasSubPaving() == true.*/
	cxsc::ivector getRootBox() const;
   
	/*! \brief get the dimensions of the subpaving this manages.

    \return 0 if this does not have a subpaving, else returns the
    dimensions of the subpaving.*/
    int getDimensions() const;
	
	/*! \brief get volume of the root box of the subpaving this manages.

    \return volume of the root box of the subpaving this manages, or
	0.0 if this has no subpaving.*/
    cxsc::real getDomainVolume() const;
    
    /*! \brief Gets number of leaf nodes in the root paving.
    
	Throws NullSubpavingPointer_Error is the subpaving that this 
    manages is a NULL pointer. 
	
	\return the total number of leaves in the subpaving managed.
    \pre hasSubPaving() == true.	     */
    size_t getRootLeaves() const;

    /*! Get a vector of the leaf node levels.

    Root is level 0, next level down is 1, etc.

    \return a vector of leaf levels, left to right order,
    */
    IntVec getLeafLevels() const;


    /*! Get a string of the leaf node levels.

    Root is level 0, next level down is 1, etc.
    Example return string "3,3,2,1"

    \return a comma separated string of leaf levels, left to right order
    */
    std::string getLeafLevelsString() const;
	
	/*! \brief Make a PiecewiseConstantFunction out of this estimator.
	  
	 \return A %PiecewiseConstantFunction with the same subpaving
	 as this with the same function values mapped onto the subpaving
	 and the same label as this. */
	subpavings::PiecewiseConstantFunction 
						makePiecewiseConstantFunction() const;

	/*! \brief Create an estimate of the function described by fobj using a 
	depth-first brute force approach.
	
	The over all aim to create an estimate of the described by fobj
	such that each of the leaf nodes of the subpaving managed by this
	has a box associated with it
	that meets the <b>interval image tolerance requirement</b> specified
	by \a tolerance.  

	The interval image tolerance requirement
	is that the diameter of the interval image of the box associated
	with any leaf node of of the subpaving managed by this
	should be less than or equal to the tolerance.
	 
	But in some cases this can result in boxes being split beyond the limits
	of the real number screen.  If a split of a leaf node that does not
	meet the interval image tolerance requirement would result in
	the child nodes' volumes being < the value of the cxsc::MinReal (the 
	smallest representable real number) then that leaf will not be split.

	This means that after the operation, the subpaving managed by this
	may have leaf nodes where the interval image tolerance requirement 
	is not met. 
	
	This operation works depth first, ie working from a leaf root
	it will successively split the root then the root's left child
	and that child's, left child, etc.  It will only start
	splitting the right child of any node already split when all the
	descendents of the left child meet the interval image tolerance
	requirement (or they cannot 
	be split further).  This may cause the machine to run out of memory
	before the process is completed.   
	 
	\note If this operation is applied to this when the subpaving managed
	by this already has leaf nodes which exceed the interval
	image tolerance requirement, those leaf nodes will be left unaltered
	by this operation.
	 
	\param tolerance describes the tolerance to be used in making the 
			estimate.
	\pre hasSubPaving() == true.
	\pre hasSubPaving() == true.
	\pre \a tol >= cxsc::MinReal.
	\post the leaf nodes of the subpaving managed by this have a
	range which is the real image of the mid-point of the box of the node under
	the function represented by this and either the interval image tolerance
	requirement is met for each leaf node or that leaf node is not
	splittable.*/
	void bruteForceEstimate(cxsc::real tolerance);
	
		
	/*! \brief Split an estimator to a specified shape.
    
   	Throws a NullSubpavings_Error if the subpaving that this manages
	is a NULL pointer.
		
	Throws a NoBox_Error if the subpaving box is empty.

	Prints a message to the standard error output if the instruction
	could not be carried out. 

    \param instruction specifies the required shape, eg "3, 3, 2, 1"
	\return true if the split was successful, false otherwise
    \pre hasSubPaving() == true.
	\post this has subpaving with shape specified by \a instruction
	and the value mapped to each node in the subpaving is the 
	real mid-image of that node's box under the function represented
	by this.*/
    bool splitToShape(std::string instruction);

			
		
	/*! \brief  Get the total integral of the function as estimated
	by this.
	
	The integral is calculated as the sum over
	all the leaves of the subpaving managed by this of the absolute
	value of the real range
	on the leaf multiplied by the volume of the 
	box represented by the leaf.
	
	\return total area between function and 0 for the function as
	estimated by this.
	\pre This must have a subpaving to manage.*/
	cxsc::real getTotalIntegralOfRealEstimate() const;
	
	/*! \brief Output the subpaving managed by this to a given stream.

    Format is a tab-delimited data giving details of leaf nodes.

    \param os is a reference to the stream to output the estimator to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \return a reference to the given stream.   */
    std::ostream & outputToStreamTabs(std::ostream & os, int prec = 5) const;

	/*! @name Output the subpaving managed by this to a txt file.

    Format is a tab-delimited file of numeric data starting with nodeName, then
    the node box volume, then the node counter, then the description of the
    node box as a tab-delimited list of interval upper and lower bounds.

    \param s the name of the txt file to send output to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
	\param confirm is a boolean controlling whether confirmation goes to
    console output.     */
	//@{
	void outputToTxtTabs(const std::string& s, int prec = 5) const;

	void outputToTxtTabs(const std::string& s, int prec, bool confirm) const;

	//@}

	/*! @name Output details of full sample (from root) to txt file.

    Format is a mixture of alpha and  numeric data.

    \param s the name of the txt file to send output to.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    \param confirm is a boolean controlling whether confirmation goes to
    console output.
    */
	//@{
    void outputRootToTxt(const std::string& s, int prec = 5) const;

    void outputRootToTxt(const std::string& s, 
									int prec, bool confirm) const;


	//@}
	
	/*! @name Output all nodes of the subpaving managed by this to
	a given stream.

    Format is a tab-delimited data giving details of all nodes.

    \param os is a reference to the stream to output to.
	\param prec the precision for output formatting. ie, number
	of decimal places, defaulting to 5.
    \return a reference to the given stream.   */
    std::ostream & outputRootToStreamTabs(std::ostream & os,
													int prec = 5) const;
	
	/*! \brief Append current state of estimator to a txt log file.

    Format is a tab-delimited file of numeric data.
    Output includes node contributions to unscaled EMP under COPERR and AIC
    and the changes in EMP that would result from splitting the node.

    \param s the name of the txt file to send output to.
    \param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
	\param prec the precision for output formatting. ie, number
	of decimal places.
    */
    void outputLog(const std::string& s, int i, int prec = 5) const;
	
	
    /*! \brief Get a string summary of this estimator's properties.
		
	A string description of this. Includes the address of the subpaving
	managed but not details of that subpaving.
	
	\return the string summary.		*/
	std::string stringSummary() const;
	
	
	private:
		
	
	/*No argument constructor is private and not implemented.*/
    explicit FunctionEstimatorReal();

	 /*Copy assignment operator is private an not implemented.    */
    FunctionEstimatorReal& operator=(const FunctionEstimatorReal& rhs);


	/*Return a pointer to the RealMappedSPnode 
	this manages. */
    RealMappedSPnode* getSubPaving() const;
	
	/*Opening line of a txt log file.

    Starts the log file with file name and date and time
    \param s the name of the txt file to send output to.    */
    void outputLogStart(const std::string& s) const;

	
							
	/*Check that the box is okay as the basis for a subpaving.
	 * 
	 * \return true if the box has at least one dimension and that
	 * no dimension of the box has a thin interval, false otherwise. */
	static bool checkBox(const cxsc::ivector& box);

	/* Handle exceptions thrown in splitting root to a specific shape. */
	void handleSplitToShapeError(RealMappedSPnode& spn);
	
	/*Handle exceptions thrown in constructors. */
	void constructor_error_handler();
	
	// data members
	/*! \brief Pointer to the root node of the subpaving tree.

    */
    RealMappedSPnode* rootPaving;

    /*! \brief A representation of the function to estimate.   */
    const MappedFobj& fobj;
	
	/*! The label.*/
	int label;

    

}; // end of FunctionEstimatorReal class declarations




	// ----------  declarations of non-member functions ----------------------


	/*! \brief Output the contents of an FunctionEstimatorReal object.

	Verbose output for an FunctionEstimatorReal object, including all boxes
	(not just leaves), data, and summary statistics.
	*/
	std::ostream & operator<<(std::ostream &os, const subpavings::FunctionEstimatorReal& fei);

} // end namespace subpavings



#endif


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
\brief RealMappedShapeNode (SubPaving) and associated non-member functions declarations

*/

#ifndef ___REALMAPPEDSHAPENODE_HPP__
#define ___REALMAPPEDSHAPENODE_HPP__

#include "realmappedspnode.hpp"
#include "cxsc.hpp"

namespace subpavings {
	
	/*! \brief RealMappedShapeNodes are nodes in the representation of a
	non-regular non-rectangular non-axis aligned polygon-paving 
	("shape") as a
	binary tree, with each node having a real value mapped to it.
	
    A %RealMappedShapeNode represents a non-regular non-rectangular 
	polygon shape with a real value mapped to it.
	RealMappedShapeNodes are linked together
    to form the tree.  The initial polygon of the non-regular 
	non-rectangular paving is the "root box"
    represented by the root node of the tree.  A polygon which has been
	split
    will be represented as node with one or two children.

    RealMappedShapeNodes have limited functionality compared to 
	SPnode types and
	are intended to be able to display the results of transforming 
	a RealMappedSPnode using transformations
	such as scaling, rotation and reflection, and also translation 
	of the origin of the axes.    */
    class RealMappedShapeNode {

    public:
	
		//typedefs
		typedef std::vector <RealMappedShapeNode* >
			Ptrs;
		typedef Ptrs::iterator
			PtrsItr;
		typedef Ptrs::const_iterator
			PtrsConstItr;
		typedef std::vector < const RealMappedShapeNode* >
			ConstPtrs;
		typedef ConstPtrs::const_iterator
			ConstPtrsItr;
		
        /*! \brief Constructor with a transformation matrix.
		
		This is constructed by applying a volume and angle-preserving
		back-transformation to the boxes in the paving 
		represented by \a rmsp, then scaling as specifed in \a scale,
		then shifting as specified in \a shift.
		
		The back-transformation
		is specified by the matrix represented by \a backtransform.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The tree rooted at this
		will have as many descendent nodes as \a rmsp and each 
		node in the tree will have a box that is the result of 
		applying the back-transformation matrix \a backtransform
		to the box of the equivalent node in \a rmsp, scaling
		it and shifting it, and then
		giving it a range value equal to the range value of
		the equivalent node in \a rmsp divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the box that results from scaling the shapes in 
		this).
		
		\param rmsp the RealMappedSPnode to use to provide the tree
		structure for this and the values to be adjusted 
		for volume scaling to give the values mapped to nodes in this.
		\param backtransform A representation of a volume preserving
		and angle preserving back-transformation matrix specified
		in row major order, ie 
		the inner vectors of \a backtransform are the rows of 
		the back-transformation matrix.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shiftg to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\pre \a backtransform should represent a square matrix 
		of dimension d where d is also the dimension of 
		the subpaving represented by \a rmsp. 
		\pre \a backtransform should be a volume-preserving 
		transformation matrix. 
		\pre \a scale and \a shift should both also be of dimension d. 
		\pre \a scale should contain only strictly positive values. */
        RealMappedShapeNode(const RealMappedSPnode& rmsp, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift);
		
		/*! \brief Constructor with no transformation matrix.
		
		This is constructed by scaling and shifting the subpaving
		represented in \a rmsp.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The tree rooted at this
		will have as many descendent nodes as \a rmsp and each 
		node in the tree will have a box that is the result of 
		scaling and shifting it the box of the equivalent node in 
		\a rmsp, and then
		giving it a range value equal to the range value of
		the equivalent node in \a rmsp divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the box that results from scaling the shapes in 
		this).
		
		\param rmsp the RealMappedSPnode to use to provide the tree
		structure for this and the values to be adjusted 
		for volume scaling to give the values mapped to nodes in this.
		 \param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shiftg to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\pre \a scale and \a shift should both also be of dimension d
		where d is the dimension of the subpaving represented by \a rmsp. 
		\pre \a scale should contain only strictly positive values. */
        RealMappedShapeNode(const RealMappedSPnode& rmsp, 
				const std::vector < double >& scale,
				const std::vector < double >& shift);
				
		/*! \brief Constructor with a transformation matrix, box and range.
		
		This is constructed by applying a volume and angle-preserving
		back-transformation to the box\a box,
		then scaling as specifed in \a scale,
		then shifting as specified in \a shift.
		
		The back-transformation
		is specified by the matrix represented by \a backtransform.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The tree rooted at this
		will consist of a single leaf node. That node will have a box
		that is the result of 
		applying the back-transformation matrix \a backtransform
		to the box \a box, scaling
		it and shifting it, and then
		giving it a range value equal to the value of \a r
		divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the box that results from scaling the shapes in 
		this).
		
		\param box the box to use to find the vertices of this 
		after scaling and shifting and back transforming.
		\param backtransform A representation of a volume preserving
		and angle preserving back-transformation matrix specified
		in row major order, ie 
		the inner vectors of \a backtransform are the rows of 
		the back-transformation matrix.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\param r The value to be used for the range of this
		after adjusting for changes in scaling.
		\pre \a backtransform should represent a square matrix 
		of dimension d where d is also the dimension of 
		the subpaving represented by \a rmsp. 
		\pre \a backtransform should be a volume-preserving 
		transformation matrix. 
		\pre \a scale and \a shift should both also be of dimension d. 
		\pre \a scale should contain only strictly positive values. */
        RealMappedShapeNode(const cxsc::ivector& box, 
				const std::vector < std::vector < double > >& backtransform,
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real r = 0.0);
		
		/*! \brief Constructor with a scale, shift, box and range.
		
		This is constructed scaling as specifed in \a scale,
		then shifting as specified in \a shift.
		
		The scaling is specifed by the elements of \a scale.  The
		ith element of scale is the scalar to be applied on the ith
		dimension of this. 
		 
		The shifting is specifed by the elements of \a shift.  The
		ith element of shidft is the shift to be applied on the ith
		dimension of this. 
		
		The tree rooted at this
		will consist of a single leaf node. That node will have a box
		that is the result of scaling \a box
		it and shifting it, and then
		giving it a range value equal to the value of \a r
		divided by the product of the
		elements in \a scale (ie adjusting for the change in the
		volume of the box that results from scaling the shapes in 
		this).
		
		\param box the box to use to find the vertices of this 
		after scaling and shifting and back transforming.
		\param scale The scaling to be applied on each dimension.  
		The ith element of \a scale is the scalar on the ith
		dimension of this.
		\param shift The shift to be applied on each dimension.  
		The ith element of \a shift is the translation on the ith
		dimension of this.
		\param r The value to be used for the range of this
		after adjusting for changes in scaling.
		\pre \a scale and \a shift should both also be of dimension d
		where d is the dimension of the subpaving represented by \a rmsp. 
		\pre \a scale should contain only strictly positive values. */
        RealMappedShapeNode(const cxsc::ivector& box, 
				const std::vector < double >& scale,
				const std::vector < double >& shift,
				cxsc::real r = 0.0);
		
		
		/*! \brief Copy constructor.       */
        RealMappedShapeNode(const RealMappedShapeNode& other);

        /*! \brief Destructor.    */
        virtual ~RealMappedShapeNode();

        /*! \brief Copy assignment operator.     */
        RealMappedShapeNode& operator=(RealMappedShapeNode rhs);

		
        /*! \brief Accessor for the dimension of theBox of a node.
        */
        size_t getDimension() const;

        /*! \brief Accessor for vertices of the box associated with the node.
		 */
        std::vector < std::vector < real > >  getVertices() const;
		
		/*! \brief Accessor for the value mapped to this.*/
        real getRange() const;

        /*! @name Accessors for links between the nodes.
        \note Pointers for parent, leftChild, and rightChild are not
        reference counted so there could potentially be problems with
        the use of returned pointers (for instance, being used to
        delete  nodes). These pointers might be better implemented
        with boost::shared_ptr.
        */
        //@{
        /*! \brief Accessor for the parent of a node.

        Returns a pointer to parent node.
        */
        RealMappedShapeNode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Returns a pointer to leftChild node.
        */
        RealMappedShapeNode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

		Returns a  pointer to rightChild node.
        */
        RealMappedShapeNode* getRightChild() const;
        //@}

        /*! \brief Get the node name.
        */
        std::string getNodeName() const;

       
        /*! \brief Check if this RealMappedShapeNode is a leaf.
		
		\return true if this has no children, false otherwise.
        */
        bool isLeaf() const;

        
		/*! @name Output for for <b>all leaves</b> of a binary tree.

        Output intended for a txt file, in numeric form only.

        \param os is the stream to send to
        \param prec is the precision used for printing.   */
        //@{
		std::ostream& leavesOutputTabs(std::ostream &os) const;
		std::ostream& leavesOutputTabs(std::ostream &os, 
												int prec) const;
		//@}
		
        /*! @name Output for <b>all nodes</b> in binary tree.

        Output for all nodes in binary tree representing a subpaving
        in tab-delimited form.

        Output intended for console output.*/
		//@{
        std::ostream& nodesAllOutput(std::ostream &os,
												int level) const;
		std::ostream& nodesAllOutput(std::ostream &os,
                                    int level, int prec) const;
		//@}
									
		


		/*! \brief Swap this and another node.
		
		Swaps all the data members of this with the other node. 
				
		\param spn a reference to the node to swap with
		\post this is identical,in terms of its data members, 
		to spn before the swap, and spn is
		identical to this after the swap.*/
		void swap (RealMappedShapeNode& rmsn); //throw()
		
		
	private:
	
		explicit RealMappedShapeNode();
	
        /*! \brief Add a left child to this.
		
		\param lChild a pointer to the new left child node to add.
        \post getLeftChild() == lChild. 
       */
        void nodeAddLeft(RealMappedShapeNode *lChild);

        /*! \brief Add a right child to this.
		
		\param rChild a pointer to the new right child node to add.
        \post getRightChild() == rChild. */
        void nodeAddRight(RealMappedShapeNode *rChild);
		
				
        /*! @name Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here shows the box
        and produces numeric tab-delimited data.

        The format when the box is an n-dimensional interval vector is:

        name [tab] volume [tab] range[tab]
        Inf(ivector[1])[tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]        */
		//@{
        std::ostream& leafOutputTabs(std::ostream &os) const;
		std::ostream& leafOutputTabs(std::ostream &os, 
												int prec) const;
		//@}
		
		void extractVertices(const ivector& box,
			const std::vector < std::vector < double > >& backtransform);
		
		void extractVertices(const ivector& box);
		
		void scaleAndShiftShape(const std::vector < double >& scale,
								const std::vector < double >& shift);
		
		static std::vector < std::vector < real > >& 
			buildVertices(const ivector& box, 
							int d,
							std::vector < std::vector < real > >& verts,
							std::vector < real > vertex);

		
		// data members
		/*! \brief The vertices of the box the node represents
        */
        std::vector < std::vector < real > >* verticesPtr;

		/*! \brief The real value mapped to the node.
        */
        real range;
        
		
        /*! \brief The name given to the node.
        */
        std::string nodeName;

		/*! \brief The node's parent.

        */
        RealMappedShapeNode* parent;

        /*! \brief The node's left child.
        */
        RealMappedShapeNode* leftChild;

        /*! \brief The node's right child.
        */
        RealMappedShapeNode* rightChild;
	
	
		/*! \brief Handle exceptions in the construction of a node.
        */
		void constructor_error_handler(); // throw()
        
    };
    // end of RealMappedShapeNode class
	


    /*! \brief Output the contents of a SubPaving.

    Output for all nodes in the subPaving.
	
	This is intended for console output or output
    to a mixed alpha and numeric file.    */
    std::ostream & operator<<(std::ostream &os, const RealMappedShapeNode& spn);



} // end namespace subpavings

/*! A specialisation of std::swap for RealMappedShapeNode types.*/
namespace std
{
	template <>
	void swap(subpavings::RealMappedShapeNode & s1, 
			subpavings::RealMappedShapeNode & s2); // throw ()
	
}



#endif

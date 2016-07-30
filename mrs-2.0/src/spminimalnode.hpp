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

#ifndef ___SPMINIMALNODE_HPP__
#define ___SPMINIMALNODE_HPP__

#include "spnode.hpp"



/*! \file 
\brief SPMinimalnode (minimal SubPaving) 
and associated non-member functions declarations.

*/


/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    /*! \brief SPMinimalnodes are nodes in the representation of a 
	subpaving as a minimal binary tree.
	
	A minimal tree means that a node will have at most one leaf child
	(ie, no node has two leaf or terminal children).

    A node represents a box (interval vector). SPMinimalnodes are linked together
    to form the tree.  The initial box of the subpaving is the box
    represented by the root node of the tree.  A box which has been split
    will be represented as node with one or two children.

    A subpaving of [<b>x</b>] (union of non-overlapping subboxes of
    [<b>x</b>]) is represented by the leaves (degenerate/ child-less)
    nodes in the tree.

    SPMinimalnodes are an adaptation of the AIASPnode class [from Jaulin,
    Kieffer, Didrit and Walter, Applied Interval Analysis, Springer, 2001,
    p. 336-348] including changes in coding, class members, and class
    structure.  SPMinimalnodes are implemented under C-XSC.

    This class replicates the set computation functionality of the
    subpaving nodes developed by [Jaulin, Kieffer, Didrit and Walter,
    Applied Interval Analysis, Springer, 2001] but also provides 
	some additional functionality.
    */
    class SPMinimalnode  : public SPnode {

    
        
    public:
        /*! \brief Default constructor.
        */
        SPMinimalnode();

        /*! \brief Initialised constructor with box.

        Initialised with one interval vector for the box.       */
        explicit SPMinimalnode(ivector& v);

        /*! \brief Initialised constructor.
        Initialised with a LabBox (labeled box).
        */
        explicit SPMinimalnode(LabBox& lb);

        /*! \brief Copy constructor.
        */
        explicit SPMinimalnode(const SPMinimalnode& other);

        /*! \brief Destructor.
        Declare as virtual because we may use SPMinimalnode as a base class.
        */
        virtual ~SPMinimalnode();

        /*! \brief Copy assignment operator.
        */
        SPMinimalnode& operator=(SPMinimalnode rhs);

        /*! @name Accessors for links between the nodes.
        These accessor methods shadow equivalent methods in the base
        class.  Thus the method used is determined  at compile time,
        not run time as would be the case if virtual methods were used.
        Because the pointers to parents and children are part of the
        base class definition, the methods have to cast the base class
        form to the derived class form in order for the pointer
        returned to be able to be used with derived class members.

        Note that pointers for parent, leftChild, and rightChild are
        not reference counted so there could potentially be problems
        with the use of returned pointers (for instance, being used to
        delete nodes). These pointers might be better implemented with
        boost::shared_ptr .
        */

        //@{
		/*! \brief Accessor for the parent of a node.
		
		Hides the base class version of this method.

        Returns a copy of the pointer to parent node.
        */
        SPMinimalnode* getParent() const;

        /*! \brief Accessor for the left child of a node.
		 
		Hides the base class version of this method.

        Returns a copy of the pointer to leftChild node.
        */
        SPMinimalnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.
		 
		 Hides the base class version of this method.
		 
		 Returns a copy of the pointer to rightChild node.
        */
        SPMinimalnode* getRightChild() const;
        //@}

		
        /*! \brief Tries to reunite nodes to form a minimal subpaving.

        Note that the nodes provided, lChild and rChild, are not
        the actual children of this, they are potential children
        which we are trying to either totally bring into
        this (if there are two of them) or to graft onto this
        if there is only one of them or if they are not both leaves,
		ie creating a minimal subpaving.
		
		This is typically a new, part-formed
        node whose formation can be completed by reuniting already two
        already-formed nodes into it or by adding on one child if only
        one is available.  nodeReunite is used in building a tree upwards
        (rather than in pruning leaves of formed tree from the bottom up).

        If two potential children are provided and they are both
        leaves, it combines the two leaf siblings into this.  If the
        potential children are not leaves or if only one potential
        child is provided, it grafts the potential child/children
        onto this as its child/children.
        */
        void nodeReunite(SPMinimalnode *lChild, SPMinimalnode* rChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a left child rather than attempting.
        to reunite two children into this.
        */
        void nodeAdoptLeft(SPMinimalnode *lChild);

        /*! \brief Builds a higher level of a tree from existing nodes.

        This adopts a right child rather than attempting.
        to reunite two children into this.
        */
        void nodeAdoptRight(SPMinimalnode *rChild);
		
		/*! \brief Check if ivector z is contained in this or children.

        Note that this checks not only the box represented by this node
        but the children as well.  Returns a BOOL_INTERVAL type, ie can
        be true, false, or indeterminate (indeterminate if the ivector
        covers the border of the subpaving).
        */
        virtual BOOL_INTERVAL spContains(const ivector& z) const;

        /*! \brief Check if rvector p is contained in this node or
        any of its children.

        Note that this checks not only the box represented by this node
        but the children as well.  Returns a BOOL_INTERVAL type but
        this can actually only be BI_TRUE or BI_FALSE
        (not BI_INDET = indeterminate).
        */
        virtual BOOL_INTERVAL spContains(const rvector& p) const;

        /*! Return container of boxes represented by commonality of 2 subpavings.

        Finds the boxes represented by the deepest common level of nodes between
        two SPMinimalnode trees  i.e. the 'outer jacket' that is the finest subpaving
        that fits around both of the subpavings represented by the given trees.

        \param boxes is a reference to a container of boxes to be filled.
        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return a container of boxes that is the outer jacket of spn1, spn2.
        */
        static BoxVec& vecLeafBoxOuterJacket(BoxVec& boxes,
                        const SPMinimalnode * const spn1, 
						const SPMinimalnode * const spn2);

        /*! Return an SPMinimalnode tree representing the outer jacket of 2 subpavings.

        The outer jacket tree is the nodes in common between two SPMinimalnode trees.
        i.e. the 'outer jacket' that fits both of the subpavings represented by
        the given trees.

        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return root of a tree whose leaves represent the boxes that cover both
        subpavings.
        */
        static SPMinimalnode* spLeafBoxOuterJacket(
										const SPMinimalnode * const spn1,
                                        const SPMinimalnode * const spn2);


        /*! Return the sum of the volume of outer jacket of two subpavings.

        The outerjacket is found using vecLeafBoxOuterJacket.

        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return the volume of the smallest subpaving containing both subpavings.
        */
        static double volOuterJacket(const SPMinimalnode * const spn1,
                                    const SPMinimalnode * const spn2);

        /*! Return container of boxes represented by intersection of 2 subpavings.

        The intersection is the boxes covering space in the subpavings
        represented by both trees.

        \param boxes is a reference to a container of boxes to be filled.
        \param spn1 root of tree representing a subpaving (can be NULL).
        \param spn2 root of tree representing another subpaving (can be NULL).
        \return a container of boxes that is the intersection of the boxes
        represented by the leaves of spn1 and the boxes represented by the
        leaves ofspn2.
        */
        static BoxVec& vecLeafBoxIntersection(BoxVec& boxes,
                        const SPMinimalnode * const spn1, 
						const SPMinimalnode * const spn2);

        /*! Return subpaving representing intersection of two subpavings.

        The intersection is the boxes covering space that that is represented by
        both trees.

        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return root of a tree whose leaves represent the boxes that cover
        space in common between both subpavings.
        */
        static SPMinimalnode* spLeafBoxIntersection(
										const SPMinimalnode * const spn1,
                                        const SPMinimalnode * const spn2);


        /*! Return the sum of the volume of intersection between two subpavings.

        The intersection is taken on leaf nodes using vecLeafBoxIntersection.

        \param spn1 first subpaving to intersect, can be NULL.
        \param spn2 second subpaving to intersect, can be NULL.
        \return the total volume of the space in common between the
        two subpavings.
        */
        static double volIntersection(const SPMinimalnode * const spn1,
                                    const SPMinimalnode * const spn2);


        /*! Return boxes of space that is the difference between two subpavings.

        The difference is the boxes covering space that is represented by spn1
        but is not in the subpaving represented by spn2.

        \param boxes is a reference to a container of boxes to be filled.
        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return boxes of space represented by spn1 that is not represented
        by spn2.
        */
        static BoxVec& vecLeafBoxDifference(BoxVec& boxes,
                        const SPMinimalnode * const spn1, 
						const SPMinimalnode * const spn2);

        /*! Return boxes of space that is difference between box and subpaving.

        The difference is the boxes covering space that is in box1 but
        is not in the subpaving represented by spn2.

        \param boxes is a reference to a container of boxes to be filled.
        \param box1 a box.
        \param spn2 root of tree representing a subpaving.
        \return boxes in box1 that is not in subpaving represented by spn2.
        */
        static BoxVec& vecBoxNodeDifference(BoxVec& boxes,
                        ivector box1, const SPMinimalnode * const spn2);

        /*! Return a tree representing difference between two subpavings.

         The difference is the boxes covering space that is represented by spn1
         but is not in the subpaving represented by spn2.

        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return a root of a tree representing a subpaving which is the space
        represented by spn1 that is not represented by spn2.
        */
        static SPMinimalnode* spLeafBoxDifference(
										const SPMinimalnode * const spn1,
                                        const SPMinimalnode * const spn2);


        /*! Return the sum of the volume of the difference between 2 subpavings.

        The difference is taken on leaf nodes using vecLeafBoxDifference.

        \param spn1 root of tree representing a subpaving.
        \param spn2 root of tree representing another subpaving.
        \return the total volume of the space represented by spn1 that is not
        represented by spn2.
        */
        static double volDifference(const SPMinimalnode * const spn1,
                                    const SPMinimalnode * const spn2);


        /*! \brief Forms a minimal SPMinimalnode subpaving from leaf boxes

        Make a minimal subpaving tree from a list of interval vectors which
        are represented by leaves of the tree.  The root of the subpaving tree
        will have Box = root, and the leaves in the list will be formable by a
        series of bisections of the given root.

        MakeTreeFromLeaves is applied recursively on bisected halves of root
        and new lists there are no boxes in the list.

        Uses Reunite() and recursive calls to MakeTreeFromLeaves() to work
        upwards to form a minimal subpaving tree.

        \param root the root of the tree.
        \param leafList a collection of non-overlapping interval vector leaves
        each of which can be formed by series of bisections of the given root.
        \return a regular minimal subpaving with rootbox root.
        */
        static SPMinimalnode* makeTreeFromLeaves(ivector& root, 
												ImageList& leafList);

        /*! \brief Forms a minimal SPMinimalnode subpaving from voxel boxes

        Make a minimal subpaving tree from a list of interval vectors which
        approximate the leaves of the tree.  The root of the subpaving tree will
        have Box = root, and the leaves in the list will have the same width in
        each dimension and be formable by a series of bisections of the given
        root.  i.e. each leaf will be the same size as each other leaf and will
        be a equal-sided hypercube.

        MakeTreeFromVoxels is applied recursively on bisected halves of root
        and new lists until either there are no boxes in the list or the
        largest diameter of the root the smallest we would expect given spacing
        or the root is equal to the largest box in the list.

        Uses Reunite() and recursive calls to MakeTreeFromVoxels() to work
        upwards to form a minimal subpaving tree

        \param root the root of the tree.
        \param leafList a collection of non-overlapping interval vector leaves
        each of which is approximately the same as the boxes formed by series
        of bisections of the given root.
        \param spacing is the spacing from the voxel data, ie the number of
        voxels in any dimension (assumed to be the same for all dimensions).
        \param dim is the number of dimensions we are using.
        \return a regular minimal subpaving with rootbox root.
        */
        static SPMinimalnode* makeTreeFromVoxels(ivector& root, 
						ImageList& leafList,
                        double spacing, size_t dim);

        /*! \brief Make a subpaving from vtk file data

        Expects a file of structured point vtk data, with 10 header lines,
        spacing on 4th line, and 255 signifying voxel in the image.
        Uses MakeTreeFromVoxels(...)

        \param filename is the name of the file to get the vtk data from.
        \return a pointer to the root node of a new tree (can be null).
        */
        static SPMinimalnode* vtkPaving(const std::string filename);
		
		void swapMin(SPMinimalnode& spn);

	protected:
        /* theBox, parent, leftChild,
        rightChild and nodeName are inherited from base class */



    };
    // end of SPMinimalnode class

	/*! \brief Check for containment of interval vector in the SubPaving.
    \param z the interval vector we are testing for containment
    in a subpaving
    \param spmn a pointer to the subpaving node we are testing to see
    if it contains z
    \return a BOOL_INTERVAL, ie BI_TRUE, or BI_FALSE,
    or BI_INDET if z overlaps the boundary of the subpaving spn.
    */
    BOOL_INTERVAL operator<=(const ivector& z, const SPMinimalnode * const spmn);


} // end namespace subpavings

// Full specializations of the templates in std namespace can be added in std namespace.
namespace std
{
	template <>
	void swap(subpavings::SPMinimalnode & s1, 
			subpavings::SPMinimalnode & s2); // throw ()
	
}

#endif

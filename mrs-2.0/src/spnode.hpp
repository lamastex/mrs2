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
\brief SPnode (SubPaving) and associated non-member functions declarations

*/

#ifndef ___SPNODE_HPP__
#define ___SPNODE_HPP__

// to use LabBox and RSSample objects
#include "MCMCPartitionGenerator.hpp"

#include "SmallClasses.hpp"

#include "spnodevisitor.hpp"

#include "sp_check_visitor.hpp"

#include "sptypes.hpp"




/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {
	
	/*! \brief SPnodes are nodes in the representation of a subpaving as a
    binary tree.

    A node represents a box (interval vector). SPnodes are linked together
    to form the tree.  The initial box of the subpaving is the box
    represented by the root node of the tree.  A box which has been split
    will be represented as node with one or two children.

    A subpaving of [<b>x</b>] (union of non-overlapping subboxes of
    [<b>x</b>]) is represented by the leaves (degenerate/ child-less)
    nodes in the tree.

    SPnodes are an adaptation of the AIASPnode class [from Jaulin,
    Kieffer, Didrit and Walter, Applied Interval Analysis, Springer, 2001,
    p. 336-348] including changes in coding, class members, and class
    structure.  SPnodes are implemented under C-XSC.

    This class replicates the set computation functionality of the
    subpaving nodes developed by [Jaulin, Kieffer, Didrit and Walter,
    Applied Interval Analysis, Springer, 2001] but also provides additional
    functionality to allow the class to become the basis for derived
    classes of subpaving nodes that can be used for
    <b>statistical set processing</b>.
    */
    class SPnode {

    public:
        /*! \brief Default constructor.
        */
        SPnode();

        /*! \brief Initialised constructor with box.
		 
		Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

        \param v interval vector for the box this represents.
		\pre v should be a proper ivector, with length >= 1.
        */
        explicit SPnode(const ivector& v);

        /*! \brief Initialised constructor.
        Initialised with a LabBox (labeled box).
		
		Throws a MalconstructedBox_Error if \a v has no practical
		(0 or negative) dimensions.

		\param lb a labeled box which supplies information about 
		box for this.
		\pre the box of the labeled box should be a proper box,
		with length >= 1.
        */
        explicit SPnode(const LabBox& lb);

        /*! \brief Copy constructor.       */
        SPnode(const SPnode& other);

        /*! \brief Destructor.    */
        virtual ~SPnode();

        /*! \brief Copy assignment operator.     */
        SPnode& operator=(SPnode rhs);

        /*! \brief Recursively rename children based on this node's nodeName.
        */
        void recursiveRename();

		/*! Return boolean to indicate if node is splittable.
		* 
		A node is splittable if there is at least one value in the 
		basic number screen between the inf and sup of the interval
		on the coordinate with maximum diameter (ie, the box 
		can be split into two child boxes which will both be considered
		smaller than the parent box) \b and
		the node volume is >= 2 * cxsc::MinReal (the
		smallest representable real number).*/
		virtual bool isSplittableNode() const;
		
        /*! \brief Accessor for the dimension of theBox of a node.
        */
        size_t getDimension() const;

        /*! \brief Accessor for theBox of a node.
		
		Throws a NoBox_Error if theBox is NULL. 

        Returns a copy of the object pointed to by theBox of a node.
        */
        ivector getBox() const;

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
        SPnode* getParent() const;

        /*! \brief Accessor for the left child of a node.

        Returns a pointer to leftChild node.
        */
        SPnode* getLeftChild() const;

        /*! \brief Accessor for the right child of a node.

		Returns a  pointer to rightChild node.
        */
        SPnode* getRightChild() const;
        //@}

        /*! \brief Get the node name.
        */
        std::string getNodeName() const;

        /*! \brief Set the node name.
        */
        void setNodeName(std::string newname);


        /*! \brief Return the volume of the box as a double.
		
		Throws a NoBox_Error if theBox is NULL.
		\return the volume of the box associated with the node.
		\pre !isEmpty(). 
		*/
        double nodeVolume() const;

        /*! \brief Return the volume of the box as a real.
        
		Throws a NoBox_Error if theBox is NULL. 
		\return the volume of the box associated with the node.
		\pre !isEmpty(). 
		*/
        real nodeRealVolume() const;

        /*! \brief Check if this SPnode is empty.

        Can only check if an actual subpaving object is empty,
		\return true if this has no box (theBox == NULL), false otherwise.
        */
        bool isEmpty() const;

        /*! \brief Check if this SPnode is a leaf.
		
		\return true if this has no children, false otherwise.
        */
        bool isLeaf() const;

        /*! \brief Check if this has a non-empty left child.
        
		\return true if this has a left child with a box,
		false otherwise.*/
        bool hasLCwithBox() const;

        /*! \brief Check if this has a non-empty right child.
        
		\return true if this has a right child with a box,
		false otherwise.*/
        bool hasRCwithBox() const;

        /*! \brief find if this node is a subleaf node.
        \return true if this is not a leaf and has any children
		of this are leaves, false otherwise.*/
        bool isSubLeaf() const;

		/*! \brief Get the split dimension for this.
		@todo We could just keep the splitDim for all spnode 
		types.  But then we'd have to adjust this when nodes
		were marginalised etc.  but this is expensive ...
		* See RealMappedSPnodeSingleRange etc as well*/
		virtual int getSplitDim() const;
		
		/*! \brief Get the split value.
        */
        virtual real getSplitValue() const;
		
		/*! \brief Get the total depth of the leaf nodes of the tree
		 * rooted at this.
		
		\return Total accumulated depths of all the leaf node
		descendants of this (depth of this if this is a leaf).*/
		unsigned long int getTotalLeafDepth() const;
		
		/*! \brief Get number of leaf descendents of this.
		
		\return Number of leaf descendents of this (1 if this is a leaf).*/
		virtual size_t getNumberLeaves() const;
		
		/*! \brief Get number of cherry descendents of this.
		
		\return Number of cherry descendents of this (1 if this is a cherry).*/
		virtual size_t getNumberCherries() const;


        /*! @name Get a container of all descendent leaf nodes as SPnodes.

        Will contain just this if this is a leaf.

        When this function is inherited in derived classes it will return the
        leaf nodes of a tree of the derived class nodes as plain SPnodes.

		\param leaves a reference to a container of node pointers to fill in.
        \return a reference to the container \a leaves 
			filled with pointers to leaf nodes.   */
		//@{
		/*! \brief Returns pointers to non-const nodes.
		 * \todo This is bad bad bad and should not happen.  Change when 
		need for it for testing is over...
		*/ 
		SPnodePtrs& getSPnodeLeaves(SPnodePtrs& leaves);

		/*! \brief Return pointers to const nodes.*/ 
		SPnodeConstPtrs& getConstSPnodeLeaves(SPnodeConstPtrs& leaves) const;
		//@}
		
		/*! @name Get a container of all sub-leaf nodes.

        Sub-leaf nodes (aka 'cherry nodes) have at least one child
		but any child must be a leaf,
        ie sub-leaves are the parents only of leaf nodes.

		Will just contain this if this is a sub-leaf.

        \param subleaves a reference to a container of node pointers to fill in.
        \return a reference to the container \a subleaves 
			filled with pointers to sub-leaf nodes.   		  */
		
		//@{
		
		/*! \brief Returns pointers to non-const nodes.
		\todo This is bad bad bad and should not happen.  Change when 
		need for it for testing is over... */ 
        SPnodePtrs& getSPnodeSubLeaves(SPnodePtrs& subleaves);

		/*! \brief Returns pointers to const nodes.*/ 
		SPnodeConstPtrs& getConstSPnodeSubLeaves(
							SPnodeConstPtrs& subleaves) const;
		
		//@}

		/*! \brief Fill in a vector of leaf node volumes.

        \param vols a container in which to store volumes.
        \return a reference to the container of volumes to which the volumes of
        all leaves descended from this have been added.
        */
        RealVec& getLeafNodeVolumes(RealVec& vols) const;


        /*! \brief Fill in a vector of leaf node levels.

        Levels are taken from length of node name.  nodeName for root is X.

        \param levels a container in which to store levels.
		\param level the level to allocate to this.
        \return a reference to the container of levels.
        */
        IntVec& getLeafNodeLevels(IntVec& levels, int level = 0) const;


        /*! \brief Get a string of leaf node levels.

        Levels start at 0 for root.

        \return a comma-separated string of levels.
        */
        std::string getLeafNodeLevelsString() const;

        /*! \brief Get the node depth from the root.

        Root has depth 0.
		
		\return the number of ancestors the node has between it and the
		root (a root returns 0);
        */
        int getNodeDepth() const;

		
        /*! \brief Get height of tree rooted at this node.

        Height of tree rooted at this is effectively the depth 
		of the deepest leaf child descendent of this,
		relative to the depth of this.  A leaf has tree height 0. 
		The ancestor of a pair of leaves has tree height 1, etc.
        */
        int getTreeHeight() const;


        /*! \brief Get the volume of the leaf with the smallest volume.

        \return the volume of the smallest (by vol) leaf node of this.
        */
        double getSmallestLeafVol() const;

        /*! \brief Get the volume of the leaf with the largest volume.

        \return the volume of the largest (by vol) leaf node of this.
        */
        double getLargestLeafVol() const;

		/*! \brief Check tree rooted at this is legal with respect
		to isSplittableNode().
		* 
		'Legal' means that all non-leaf nodes in the tree are 
		splittable, ie return isSplittableNode() = true; 
        
		\return true if all non-leaf nodes in the tree rooted at
		this are splittable, false if any non-leaf node in the tree 
		rooted at this is not splittable.*/
        virtual bool checkTreeStateLegal() const;

        /*! \brief Check if this has a leaf sibling.
        
		\return true if this has a leaf sibling, false otherwise.*/
        bool hasLeafSibling() const;

		/*! \brief Get the volume of the root node of the tree for this.
		
		\return The volume of the root node that is the ultimate
		ancestor of this node.        */
        cxsc::real getRootVolume() const;
									
		/*! \brief Accept a visitor of the type SPnodeVisitor.
		
		Recursively asks children to accept visitor too.
		 
		 \param visitor The visitor to accept.
        */
        virtual void acceptSPnodeVisitor(const SPnodeVisitor& visitor);
		
		/*! \brief  Accept an SPCheckVisitor.
		
		\param visitor a visitor capable of
		performing some check on this.
		\post this is unchanged.*/
		virtual void acceptSPCheckVisitor(const SPCheckVisitor& visitor) const;


		/*! @name Output for for <b>all leaves</b> of a binary tree.

        Output intended for a txt file, in numeric form only.

        \param os is the stream to send to
        \param prec is the precision used for printing.   */
        //@{
		virtual std::ostream& leavesOutputTabs(std::ostream &os) const;
		virtual std::ostream& leavesOutputTabs(std::ostream &os, 
												int prec) const;
		//@}
		
        /*! @name Output for <b>all nodes</b> in binary tree.

        Output for all nodes in binary tree representing a subpaving
        in tab-delimited form.

        Output intended for console output.*/
		//@{
        virtual std::ostream& nodesAllOutput(std::ostream &os,
												int level) const;
		virtual std::ostream& nodesAllOutput(std::ostream &os,
                                    int level, int prec) const;
		//@}
									
		/*! \brief Output details of a specific node in a tree.

        This is intended for console output or output
        to a mixed alpha and numeric file.        */
		//@{
        virtual std::ostream& nodePrint(std::ostream &os) const;
		virtual std::ostream& nodePrint(std::ostream &os,
									int prec) const;

		//@}

		/*! @name Append current state to a txt log file.

		Output is in format used by leavesOutputTabs.

		\param s the name of the txt file to send output to.
		\param i the number of pass (ie, 0, 1, 2, 3 etc) in process.
		\param prec the precision for output formatting.*/
		//@{
		void outputLog(const std::string& s, int i) const;
		
		void outputLog(const std::string& s, int i, int prec) const;
		//@}

        /*! \brief Make a .dot graph file from an SPnode tree structure.

        Makes a simple .dot graph from the SPnode tree using node names and the
        .png image for this graph.

        \pre a constructed SPnode tree
        \post a .dot file and a .png in the same directory as the program
        creating it was run in.
        */
        virtual bool outputGraphDot() const;

		/*! \brief Check if the box a node represents contains a datapoint p.

        \param p the value of the data point being tested for containment in
        the box represented by this node.
        \param childInd indicates whether this should be considered
        to be a left child or a right child (ie where we need to take
        splitting dimension and value into account) or as a parent
        node (default).
		
		If a node is being considered as a child node (ie childInd is 
		ON_RIGHT or ON_LEFT) and actaully has a parent node, 
		then is is assumed that the data 
		point would have been in the box associated with the parent
		node, and the question now is just whether it falls into 
		the left or right child. Thus no test for full containment of
		the data point in the box is carried out, and this method
		just checks whether the data point would have fallen 
		from the parent into the left or right child.
		 
		If childInd is ON_PARENT, or if the node has no parent node,
		then the method checks for full containment of the data point
		within the box associated with this node.

        Note that the interval on the parent's split dimension
        of the right child's box is closed at the split value and
        the interval of the left child's box is open.  A datapoint
        whose element in the split dimension is exactly the split
        value should be assessed to be in the right child's box but
        not the left child's box.
        */
        virtual bool nodeContains(const rvector& p,
                        OPERATIONS_ON childInd = ON_PARENT) const;

        
        /*! \brief Expand a leaf node to make two leaves as children.

        Equivalent to bisecting a box in a regular subpaving.  Makes
        two new sibling child nodes of this one and grafts them on.
        
		Throws a NoBox_Error if theBox is NULL. 
		If this is not a leaf, no expansion will take place.

		\param comp is the dimension on which to bisect theBox.
		\pre this has a non-empty box.
		\post this will not be a leaf.
        */
        virtual void nodeExpand(int comp);

        /*! \brief Expand a leaf node to have two children and pass
        data  down to the children with no further splitting.

        Finds the splitting dimension as the first longest dimension
		of this's box.
		
		Throws a NoBox_Error if theBox is NULL. 
		If this is not a leaf, no expansion will take place.

		\pre this has a non-empty box.
		\post this will not be a leaf.*/
        virtual void nodeExpand();

        /*! \brief Reabsorbs both children of the node.

        \post this will be a leaf node. */
        virtual void nodeReabsorbChildren();

        /*! \brief Add a left child to this.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if this has no box or if the
		prospective child has no box.</li>
		<li>Throws a logic_error if there is an existing
		child on this side.</li>
		<li>Throws a UnfulfillableRequest if the prospective
		child's box is the same as this's box.  Note: this can happen 
		if the dimension split on is very small already, so that 
		using the available precision, one 'half' is the whole.</li>
		<li>Throws an IncompatibleDimensions_Error if the prospective
		child's box is incompatible with this's box or, if there
		is another child, with that other child's box, if any.</li>
		</ul>

		\param lChild a pointer to the new left child node to add.
        \pre This is not empty, the prospective child is not empty, this
		has no existing left child, the propective child's box
		is compatible in size and dimension with this's box and with
		the box of any already existing right child.
		\post getLeftChild() == lChild. 
       */
        void nodeAddLeft(SPnode *lChild);

        /*! \brief Add a right child to this.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if this has no box or if the
		prospective child has no box.</li>
		<li>Throws a logic_error if there is an existing
		child on this side.</li>
		<li>Throws a UnfulfillableRequest if the prospective
		child's box is the same as this's box.  Note: this can happen 
		if the dimension split on is very small already, so that 
		using the available precision, one 'half' is the whole.</li>
		<li>Throws an IncompatibleDimensions_Error if the prospective
		child's box is incompatible with this's box or, if there
		is another child, with that other child's box, if any.</li>
		</ul>

        \param rChild a pointer to the new right child node to add.
        \pre This is not empty, the prospective child is not empty, this
		has no existing right child, the propective child's box
		is compatible in size and dimension with this's box and with
		the box of any already existing left child.
		\post getRightChild() == rChild. */
        void nodeAddRight(SPnode *rChild);

        /*! \brief Gets a string of child nodes names.

        Left to right order.

        \return a string of child nodes names in pairs, left to right.
        */
        std::string getChildNodeNames() const;
        
		/*! \brief Reshape so that the tree rooted at this has shape that
		is the union of this shape and the shape of another tree.
		
		Throws a NoBox_Error if this has no box or if \a other has no box. 
		
		Throws an IncompatibleDimensions_Error if boxes of this and \a other
		are not the same.
		
		Throws an std::runtime_error if \a other has an illegal state
		(see checkTreeStateLegal()).
		* 
		\param other is the tree to make the union against.
		\pre This has a box and that box is identical to the box of \a other. 
		\post the tree rooted at this has shape that is the
		union of the shape of this before the operation and the shape of 
		\a other.  \a other is unchanged.      */
		virtual void reshapeToUnion(const SPnode& other);
	
		/* docs */
		virtual bool randomSplitRootAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						bool saveInstructions = false);
	
		virtual bool randomNaturalSplitRootAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner,
						bool saveInstructions = false);
	
		
	
        /*! \brief Split a root paving to a specified shape.
		
		\note:  This method will create a non-minimal tree only, ie 
		a split will result in two child nodes. There is no equivalent
		method for creating a minimal tree.  

        The instruction specifies the required shape in terms of the depth of
        the leaf nodes, in left to right order.  The depth of a leaf node is
        equivalent to the number of bisections of the root box required to make
        the box represented by that leaf.  i.e., the root has depth 0 and if
        that were bisected, the root node would have two child nodes each of
        level 1.  Bisection of any one of the boxes represented by a
        child would give two more children, each of level 2 (2 bisections), etc.

        Leaf levels in \a instruction can be separated by commas, spaces, or 
		both.
		
		For example, an instruction string  "3, 3, 2, 1" would give an SPSnode
        tree with 4 leaves, 2 of which would be level 3 (i.e. representing
        boxes resulting from 3 bisections of the root box, each of which would
        have volume 1/8 the size of the root box).  Another leaf would represent
        a box resulting from 2 bisections of the root box (volume 1/4 that of
        the root box) and the 'right-most' leaf (in a drawing of the tree) would
        be the result of a single bisection of the root box and would have half
        the volume of the root box.  This is a valid instruction string because
        it is possible to get leaves of those levels by a series of successive
        bisections of the root box and the volume of the leaves adds up to the
        volume of the root box.
		
		Throws the following exceptions:
		<ul>
		<li>Throws a NoBox_Error if the box is NULL.</li>
		<li>Throws a NonRootNode_Error if this is not a root node
		(ie if this has a parent node).</li>
		<li>Throws an std::invalid_argument exception if the \a
		instruction constain invalid characters (anything other
		than digits, whitespace or commas).</li> 
		<li>Throws an std::logic_error if the instruction is valid but does
		not describe an achievable tree shape.</li>
		</ul>

        \param instruction specifies the required shape, eg "3, 3, 2, 1".
        \return true if the instruction could be successfully carried out,
        false if the instruction could not be carried out successfully.
		\pre The instruction must be valid and describe an achievable
		tree shape.  This must be a non-empty root node.  
		\post getLeafNodeLevelsString() == instruction */
        virtual bool splitRootToShape(std::string instruction);

		/*! \brief Try to get tree rooted at this into shape at least as 
		deep as described by instruction \a reqDepths.
		
		Tries to ensure that this has at least shape of \a reqDepths, but 
		allows the tree rooted at this to also be more split at some 
		or all nodes.
		
		This can only follow the instruction if the nodes it is required
		to split are splittable nodes (to be split a node 
		must return isSplittableNode() = true).
		
		\param reqDepths a collection of leaf levels to try to split
		to.
		\return true if this was able to follow the instruction, false
		otherwise (ie if instruction could not be followed because
		nodes were not splittable). 
		\pre this has a box.
		\pre this is a root node.
		\pre \a reqDepths describes a valid string (including is not empty).
		\post If returns true, the tree rooted at this has at least the
		shape described by \a reqDepths; if returns false, the tree rooted
		at this will have been split according to \a reqDepths, working
		from the right, until a node that cannot be split needs to 
		be split (all attempts to split will then have ceased).*/
		virtual bool splitRootAtLeastToShape(std::vector < size_t > reqDepths);
		
		/* assumes a binary tree type instruction, ie any non-root node has a sibling
		 * will expand or prune as required but if it expands, will use nodeExpand for this 
		 * type of node.
		 * No checking that a node is expandable here: instructions are assumed to be valid
		 * for a node like this
		 * Make bool return type and allow subtypes to override if necessary with 
		 * a version that allows for reshape to fail*/
		virtual bool reshapeToMirrorValidSplitPartitions(
						const std::vector < unsigned long int >& instructions,
						unsigned long int numLeaves);

		unsigned long int getPartitions(
				std::vector < unsigned long int >& partitions) const;

		/*! \brief Swap this and another node.
		
		Swaps all the data members of this with the other node. 
				
		\param spn a reference to the node to swap with
		\post this is identical,in terms of its data members, 
		to spn before the swap, and spn is
		identical to this after the swap.*/
		void swap (SPnode& spn); //throw()
		
		/*! \brief Get a string summary of this node's properties.
		
		No recursion, ie just this node.
		
		\return the string summary.		*/
		virtual std::string nodeStringSummary() const;
		
		
	protected:
	
		/*! \brief Splitting according to instruction string.

        Gets the first level in the string, checks if it is already at that
        level, if not splits and recurses to children in order left then right,
        if it is already at that level strips this off the string and returns
        the remaining string to the caller.
        \param instruction is a string specifiying shape, eg "3, 3, 2, 1"
        */
        std::string splitLeft(std::string instruction);
		
		/*! \brief Try to get tree rooted at this into shape at least as 
		deep as described by instruction \a reqDepths.
		
		Tries to ensures that this has at least shape of \a reqDepths, but 
		allows the tree rooted at this to also be more split at some 
		or all nodes.
		
		This can only follow the instruction if the nodes it is required
		to split are splittable nodes (to be split a node 
		must return isSplittableNode() = true).
		
		\param reqDepths a collection of leaf levels to try to split
		to.
		\param myDepth the depth this node is in the tree.
		\return true if this was able to follow the instruction, false
		otherwise (ie if instruction could not be followed because
		nodes were not splittable). 
		\pre \a reqDepths is not empty*/
		virtual bool _splitAtLeastToShape(std::vector < size_t >& reqDepths,
								size_t myDepth);

		/* assumes a binary tree type instruction, ie any non-root node has a sibling
		 * will expand or prune as required but if it expands, will use nodeExpand for this 
		 * type of node */
		virtual size_t _reshapeToMirrorValidSplitPartitions(
						const std::vector < unsigned long int >& partitions,
						unsigned long int numLeaves,
						size_t currentIndex);
		
		/*docs */
		unsigned long int _getMyPartitions(
			std::vector < unsigned long int >& partitions) const;
		
		/*! \brief Reshape the tree rooted at this so that it has the
		shape that is the non-minimal union of the tree orginally
		rooted at this and the tree rooted at the node pointed to
		by \a other.*/
		virtual void _reshapeToUnion(const SPnode * const other);
		
		/*! \brief Randomly partition the tree rooted at this so
		that it has at least numLeaves.  If this is already more
		split than the partitioner directs that is okay.  If this 
		or any of its children
		cannot split as directed by the partitioner (because of 
		isSplittableNode()) then return false.  Otherwise return true.*/
		virtual bool _randomSplitAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner);

		virtual bool _randomNaturalSplitAtLeast(
						unsigned long int numLeaves,
						const MCMCPartitionGenerator& partitioner);
		
		unsigned long int _getTotalLeafDepth(
				unsigned long int thisDepth) const;
		
        /*! @name Output for a node in a binary tree, tab-delimited

        Output intended for a txt file, in numeric form only.

        Replaces the format that that the cxsc::<< operator produces
        for interval vectors.   The format used here shows the box
        and produces numeric tab-delimited data.

        The format when the box is an n-dimensional interval vector is:

        name [tab] volume [tab] Inf(ivector[1])
        [tab] Sup(ivector[1].[tab] . .
        [tab] Inf(ivector[n]) [tab] Sup(ivector[n]        */
		//@{
        virtual std::ostream& leafOutputTabs(std::ostream &os) const;
		virtual std::ostream& leafOutputTabs(std::ostream &os, 
												int prec) const;
		//@}
		
		// data members
		/*! \brief Pointer to the interval vector the node represents
        */
        ivector* theBox;

        
        /*! \brief The name given to the node.
        */
        std::string nodeName;

		/*! \brief The node's parent.

        This data member is not present in the AIASPnode class.
        This allows us to go up from the child node to its parent.
        */
        SPnode* parent;

        /*! \brief The node's left child.
        */
        SPnode* leftChild;

        /*! \brief The node's right child.
        */
        SPnode* rightChild;
	
	private:
	
		/*! \brief Handle exceptions in the construction of a node.
        */
		void constructor_error_handler(); // throw()
        
    };
    // end of SPnode class
	


    /*! \brief Output the contents of a SubPaving.

    Output for all nodes in the subPaving.
	
	This is intended for console output or output
    to a mixed alpha and numeric file.    */
    std::ostream & operator<<(std::ostream &os, const SPnode& spn);

    /*! \brief Check if a SubPaving is empty.
    \return true if the given node is NULL or
    SPnode::isEmpty() is true */
    bool isEmpty(const SPnode * const spn);

    /*! \brief Check if a node is a leaf.
	
    \return true if the given node is not NULL
    and SPnode::isLeaf() is true.    */
    bool isLeaf(const SPnode * const spn);


    /*! \brief Get the volume of the subpaving represented by spn.    */
    double spVolume(const SPnode * const spn);

    /*! \brief Get the number of leaves of a tree (boxes in the 
	subpaving). */
    size_t spLeaves(const SPnode * const spn);

    

} // end namespace subpavings

/*! A specialisation of std::swap for SPnode types.*/
namespace std
{
	template <>
	void swap(subpavings::SPnode & s1, 
			subpavings::SPnode & s2); // throw ()
	
}



#endif

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

/*! \file AIAsubpaving.hpp
\brief AIASPnode (AISSubPaving) declarations

*/

#include "cxsc.hpp" //to use ivectors
#include "toolz.hpp"//to use methods declared in toolz.hpp
#include <iostream> // to use standard input and output
#include <fstream>  // include fstream to output a file from spImage
#include <list>     // to use std::list

/*! \brief Define type "Interval Booleans"
This is used to extend the usual booleans TRUE and FALSE for use with 
intervals. 	With intervals the result of a test may be indeterminate at a 
particular level of accuracy rather than clearly true or false.
 */
typedef enum {BI_TRUE,BI_FALSE,BI_INDET}
AIA_BOOL_INTERVAL;

/** @name Typedefs for function pointers
*/
//@{
/*! \brief Define type "Pointer to an interval boolean test"

The test is for containment of the inclusion function image of the given 
interval vector in some subpaving Y.  The inclusion function is specified in 
the body of the interval boolean test.  The test returns BI_INDET if the 
inclusion function image overlaps the borders of the subpaving.  As implemented 
here, the test relies on the use of globals for the subpaving Y so this is not 
passed as parameter to the function but is identified within the body of the 
function.
*/
typedef AIA_BOOL_INTERVAL (*AIA_PIBT)(const ivector&);

/*! \brief Define type "Pointer to an interval vector function"

PIVF is an interval vector inclusion function, ie a function which returns the 
interval vector which encloses f(x) for f as specified in the function and 
interval vector given in the function parameter. 
*/
typedef ivector (*AIA_PIVF)(const ivector&);

//! forward class declarations
class AIASPnode;

//! Makes AIASubPaving alias of pointer to a AIASPnode
typedef AIASPnode* AIASubPaving;
//@}

/*! \brief AIASubPaving node class

A class based on the Sub-paving node class from 
[AIA2001, p. 336-348] implemented over C-XSC.  AIASPnode is a node in the tree 
representation of a regular subpaving.  A node represents a box 
(interval vector). SPnodes are linked to children to form a binary tree. A 
subpaving of [<b>x</b>] (union of non-overlapping subboxes of [<b>x</b>]) is 
represented by the leaves (degenerate/ child-less) nodes in the tree.  This 
class replicates the set computation functionality of the subpaving nodes 
developed in AIA2001.
*/
class AIASPnode
{
  private:
    //! The interval vector the node represents
    ivector* theBox;
    //! The node's left child
    AIASubPaving leftChild;
    //! The node's right child
    AIASubPaving rightChild;

  public:
    //! Default constructor
    AIASPnode()
      {theBox=NULL; leftChild=NULL; rightChild=NULL;}

    //! Initialized constructor
    AIASPnode(ivector& v)
      {theBox=new ivector(v); leftChild=NULL; rightChild=NULL;}

    //! Copy constructor
    AIASPnode(const AIASPnode& n);

    //! Destructor
    ~AIASPnode()
      {delete theBox; delete leftChild; delete rightChild;}

    /*! @name Friend functions
    We have followed the implementation in [AIA2001, p. 336-348], which uses 
    friend functions in order to give non-member functions access to the 
    private data members of the class. */
    //@{
    //! Returns a copy of theBox of an AIASubPaving. */
    friend ivector Box(AIASubPaving a)
      { return (*a->theBox); }

    //! Check if the AIASubPaving is empty
    friend bool IsEmpty(AIASubPaving a)
      {return ( (a==NULL) || (a->theBox==NULL) ); }

    //! Check if the AIASubPaving is a leaf
    friend bool IsLeaf(AIASubPaving a)
    {
      return ( !IsEmpty(a) && IsEmpty(a->leftChild)
        && IsEmpty(a->rightChild) );
    }

    /*! \brief Volume of an AIASubPaving
    
    Sums the volumes of the leaves of the AIASubPaving. 
    */
    friend double Volume(AIASubPaving a);

    //! Number of leaves of an AIASubPaving
    friend int NbLeaves(AIASubPaving a);

    //! Check for containment of interval vector or box in the AIASubPaving
    friend AIA_BOOL_INTERVAL operator<=(const ivector&, AIASubPaving);

    /*! \brief Output the contents of an AIASubPaving
    
    This output format is not the same as used in the example code 
    provided for [AIA2001], but is adapted for easy reading and rendering 
    in Matlab.
    */
    friend std::ostream & operator<< (std::ostream &, AIASubPaving);

    /*! \brief Set Inversion Via Interval Analysis method taken from AIA2001.
    
    \return a minimal regular subpaving covering the reciprocal image of a 
      subpaving under some function f
    \param BoolTest a (*AIA_PIBT)() which specifies the inclusion function of 
      f and tests the inclusion function image [f][x] of an interval vector 
      [x] for containment in an SPnode.  In this implementation, the body of 
      BoolTest specifies both the subpaving we want the image of and the 
      function f forming the image
    \param A a subpaving whose box forms the interval vector passed to BoolTest
    \param eps the precision with which the returned subpaving should be formed
    SIVIA progressively subdivides the boxes of the initial search subpaving 
    and calls itself recursively to select or reject or retest the resulting 
    subpavings until the desired level of precision, specified by eps, in 
    forming the subpaving to be returned has been achieved.
    */
    friend AIASubPaving Sivia(AIA_PIBT BoolTest, AIASubPaving A, double);

    /*! \brief Expand a leaf node to have two child nodes.
    
    \param A a pointer to the node to be expanded.
    \param comp is the dimension to split along.
    */
    friend void Expand(AIASubPaving A, int comp);

    /*! \brief Return a minimal subpaving from two sibling subpavings. 
    
    If two potential children are provided and they are both leaves, then 
    Reunite combines the two leaf siblings into this.  If the potential 
    children are not leaves or if only one potential child is provided, grafts 
    the potential child/children onto this as its child/children.
      
    \return a minimal subpaving from two sibling subpavings
    \param lChild a pointer to the leftChild node to be reunited
    \param rChild a pointer to the rightChild node to be reunited
    \param x is the box of the new subpaving to be returned
    */
    friend AIASubPaving ReUnite(AIASubPaving lChild, AIASubPaving rChild, 
                                ivector x);

    /*! \brief Mince up a subpaving.
   
    \param A a pointer to the node whose box is to be minced
    \param eps the maximum diameter any box in the subpaving should be
    Mince minces recursively until each leaf has maximum diameter smaller than 
    eps. 
    */
    friend void Mince(AIASubPaving A, double eps);

    /*! \brief Evaluate the image.
    
    \param f a pointer to an interval vector function which returns the image 
    box under f of some interval vector
    \param A the node of the subpaving to be evaluated
    \param evalImages a container of image interval vectors
    \param hull the interval hull of the image interval vectors
    Fills in the list of images using f of the subpaving boxes and update the 
    hull of all the images. 
    */
    friend void Evaluate(AIA_PIVF, AIASubPaving A, list<ivector>& evalImages, 
                         ivector& hull);

    /*! \brief Forms a minimal image subpaving from a list of interval vector 
      images.
    
    Uses Reunite() and recursive calls to Regularize() to work up from leaf 
    nodes to form a minimal subpaving
    \return a regular mimimal subpaving with root box hull
    \param hull the interval hull of all the interval vectors in the image list
    \param ivectorList a collection of possibly overlapping interval vectors to 
    be made into a regular subpaving
    \param eps the precision with which the returned subpaving should be 
    formed
    There are some minor changes in this function compared to the 
    implementation in AIA2001, p. 346.  We use the std::list to store the 
    interval vector images instead of a [constantly sorted] set.  We sort the 
    list, using volCompare(), only when we need to.  We sort smallest to 
    largest (the natural std::list sort order) (the Jaulin et al images are 
    sorted largest to smallest) and hence we find the largest image box at the 
    end of the the list rather than at the start.
    */
    friend AIASubPaving Regularize(ivector& hull, list<ivector>& ivectorList, 
                                   double eps);

    /*! \brief Creation of an image subpaving (ImageSp) with 
      Interval Analysis from AIA2001.
      
    \return a minimal regular subpaving coveringing the image of A under some 
    function f
    \param f a (*AIA_PIVF) which specifies the inclusion function of f and 
    returns the inclusion function image [f][x] of an interval vector [x]
    \param A the subpaving for which we wish to find a subpaving covering the 
    image under f
    \param eps the precision with which the returned subpaving should be formed
    ImageSp uses Mince() to chop up A and then Evaluate() and Regularize() to 
    find a regular minimal subpaving covering the set of images of the boxes of 
    the minced A.
    */
    friend AIASubPaving ImageSp(AIA_PIVF, AIASubPaving A, double eps);
    //@}

};

/*! \relates AIASPnode \brief Compare volumes of two boxes

A function for comparing ivectors based on volume.  
Used in sorting a list of ivectors ordered by volume.  Will abort if the index 
sets of the two ivectors are different sizes.  Uses Volume() as defined in 
toolz.hpp
\return TRUE if the Volume of a is strictly less than Volume of b
*/
bool volCompare(const ivector &a, const ivector &b);

//----------------------------------------------------------

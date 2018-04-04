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

#ifndef ___SPALG_HPP__
#define ___SPALG_HPP__

#include "sptypes.hpp"

/*! \file spalgorithms.hpp
\brief SPnode (SubPaving) algorithm function declarations

*/

/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    class SPnode;

    /*! \brief Mince up a subpaving.

    Mince minces recursively until each leaf has maximum diameter
    smaller than eps.

    Mince is now a non-friend non-member function
    (compare to the AIASPnode::Mince())
    which now uses public member functions in the SPnode class.
    \param spn a pointer to the node whose box is to be minced
    \param eps the maximum diameter any box in the subpaving should be
    */
    void Mince(SPnode * const spn, double eps);

    /*! \brief Evalutes images of subpaving.

    Fills in the list of images using f of the subpaving boxes
    and update the hull of all the images.

    Evaluate is now a non-friend non-member function
    (compare to AIASPnode::Evaluate())
    which now uses public member functions in the SPnode class.

    \param f a pointer to an interval vector function which returns
    the image box under f of some an interval vector
    \param spn the node of the subpaving to be evaluated
    \param evalImages a container of image interval vectors
    \param hull the interval hull of the image interval vectors
    */
    void Evaluate(PIVF f, const SPnode * const spn,
                ImageList& evalImages, ivector& hull);

    /*! \brief Expand a leaf node to have two child nodes.

    Uses nodeExpand() method for the type of spn (base/derived) at runtime
    and this nodeExpand() method will find the dimension to split on.

    Expand is now a non-friend non-member function
    (compare to AIASPnode::Expand())
    which now uses public member functions in the SPnode class.
    \param spn a pointer to the node to be expanded.
    */
    void Expand(SPnode * const spn);

    /*! \brief Expand a leaf node to have two child nodes.

    Uses nodeExpand() method for the type of spn (base/derived) at runtime.

    \param spn a pointer to the node to be expanded.
    \param comp is the dimension to split on.
    */
    void Expand(SPnode * const spn, int comp);
    
    /*! \brief Expand a leaf node to have two child nodes.
    
    \param spn a pointer to the node to be expanded.
    */
    void ExpandWithValid(SPnode * const spn, bool boolVal);

    /*! \relates SPnode
    \brief Set Inversion Via Interval Analysis.

    SIVIA progressively subdivides the boxes of the initial search
    subpaving and calls itself recursively to select or reject or retest
    the resulting subpavings until the desired level of precision,
    specified by eps, in forming the subpaving to be returned has been
    achieved.

    Sivia is now a non-friend non-member function
    (compare to AIASPnode::Sivia())
    which now uses public member functions in the SPnode class.

    \param BoolTest a (*PIBT)() which specifies the inclusion function
    of f and tests the inclusion function image [f][x] of an interval
    vector [x] for containment in an SPnode.
    \param toInvert the subpaving we are attempting to find the reciprocal
    image of, the SPnode passed to BoolTest
    \param search a subpaving whose box forms the interval vector
    passed to BoolTest.
    \param eps the precision with which the returned subpaving
    should be formed.
    \return a minimal regular subpaving covering the reciprocal image
    of toInvert under some function f.
    */
    SPnode* Sivia (PIBT BoolTest, const SPnode * const toInvert,
                SPnode * const search, const double eps);

    /*! \relates SPnode
    \brief  Creation of image subpaving with Interval Analysis.

    ImageSp uses Mince() to chop up spn and then Evaluate() and
    Regularize() to find a regular minimal subpaving covering the
    set of images of the boxes of the minced spn.

    ImageSp is now a non-friend non-member function
    (compare to AIASPnode::ImageSp())
    which now uses public member functions in the SPnode class.
    \param f a (*PIVF) which specifies the inclusion function of f and
    returns the inclusion function image [f][x] of an interval vector [x].
    \param spn the subpaving for which we wish to find a subpaving
    covering the image under f.
    \param eps the precision with which the returned subpaving should
    be formed.
    \return a minimal regular subpaving covering the image of a
    subpaving under some function f.
    */
    SPnode* ImageSp(PIVF f, SPnode* spn, double eps);


    /*! \relates SPnode
    \brief Creation of non-minimal image subpaving with Interval Analysis.

    ImageSpNonMinimal uses Mince() to chop up spn and then Evaluate() and
    RegularizeNonMinimal() to find a regular non-minimal subpaving
    covering the set of images of the boxes of the minced spn.
    \param f a (*PIVF) which specifies the inclusion function of f and
    returns the inclusion function image [f][x] of an interval vector [x].
    \param spn the subpaving for which we wish to find a subpaving
    covering the image under f.
    \param eps the precision with which the returned subpaving
    should be formed.
    \return a non-minimal regular subpaving covering the image of a
    subpaving under some function f.
    */
    SPnode* ImageSpNonMinimal(PIVF f, SPnode *spn, double eps);



} // end namespace subpavings

#endif

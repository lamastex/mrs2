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

#ifndef ___SPTEMPLATES_HPP__
#define ___SPTEMPLATES_HPP__

#include "toolz.hpp"
#include "sptools.hpp"

using namespace std;

/*! \file sptemplates.hpp
\brief Templatised functions using node type concepts.

*/


/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {


/*! \note Node types used with these templates should implement
default constructor and constructor taking ivector argument and also implement
nodeExpand(), nodeReunite(), nodeAdoptLeft() and nodeAdoptRight() functions.
*/

/*!\brief Tries to reunite two nodes into to form a single leaf.

Note that the nodes provided, lChild and rChild, are potential children for
the node to be created and returned in this function.  nodeReunite is used
in building a tree upwards (rather than in pruning leaves of formed tree from
the bottom up).

If two potential children are provided and they are both leaves, combines
the two leaf siblings into the returned node.  If the potential children are
not leaves or if only one potential child is provided, graft the potential
child/children onto the node to be created and returned as its child/children.

Calls nodeReunite method for the type of spn (base or derived) at runtime.

Reunite is now a templatised non-friend non-member function
(compare to AIASPnode::Reunite())
which now uses public member functions in the SPnode class.
\param lChild a pointer to the leftChild node to be reunited
\param rChild a pointer to the rightChild node to be reunited
\param x is the box of the new subpaving to be returned
\return a minimal subpaving from two sibling subpavings
*/
template<typename T>
T *Reunite(T *lChild, T *rChild, ivector x)
{
    T* newNode = NULL; // pointer to new node to be returned

    try
    {
        newNode = new T(x);

        // both proposed children are empty, return null
        if(isEmpty(lChild) && isEmpty(rChild)) {
            newNode = NULL;
        }

        // only given a right child, left child is NULL
        if (isEmpty(lChild) && !isEmpty(rChild)) {
            //graft right child on
            newNode->nodeAdoptRight(rChild);
        }

        // only given a left child, right child is NULL
        if (!isEmpty(lChild) && isEmpty(rChild)) {
            //graft left child on
            newNode->nodeAdoptLeft(lChild);
        }

        // both children exist
        if (!isEmpty(lChild) && !isEmpty(rChild)) {

            // otherwise reunite the proposed children on this node
            newNode->nodeReunite(lChild, rChild);
        }
    }

    catch (bad_alloc&)
    {
        std::cout << "Error allocating memory in Reunite()"
            << std::endl;
        throw;
    }

    return newNode;
}


/*! \brief Adopt nodes to build a non-minimal subpaving.

Make a new node and graft the two proposed children on as children.
\param lChild a pointer to the leftChild node to be adopted
\param rChild a pointer to the rightChild node to be adopted
\param x is the box of the new subpaving to be returned
\return a non-minimal subpaving from two sibling subpavings
*/
template<typename T>
T *Adopt(T *lChild, T *rChild, ivector x)
{
    T* newNode = NULL; // pointer to new node to be returned

    try
    {
        newNode = new T(x);


        // both proposed children are empty, return null
        if(isEmpty(lChild) && isEmpty(rChild)) {
            newNode = NULL;
        }

        // add the new right child if there is one
        if (!isEmpty(rChild)) {
            //graft right child on
            newNode->nodeAdoptRight(rChild);
        }

        // add the new left child if there is one
        if (!isEmpty(lChild)) {
            //graft left child on
            newNode->nodeAdoptLeft(lChild);
        }

    }

    catch (bad_alloc&)
    {
        std::cout << "Error allocating memory in Adopt()" << std::endl;
        throw;
    }

    return newNode;
}

/*! \brief Forms a minimal image subpaving

Make a minimal subpaving from a list of interval vector images.
The root of the subpaving will have Box = hull, where hull has already been
formed from the union of all the ivectors in the ivectorList.
Regularize is applied recursively on bisected half of hull and new lists
until either there are no images in the list or the diameter of the hull is
below eps.

Uses Reunite() and recursive calls to Regularize() to work upwards
to form a minimal subpaving

Regularize is now a templatised non-friend non-member function
(compare to AIASPnode::Regularize())
which now uses public member functions in the SPnode class.
\param hull the interval hull of all the interval vectors in the image list
\param ivectorList a collection of possibly overlapping interval vectors
to be made into a regular subpaving
\param eps the precision with which the returned subpaving should be formed
\return a regular minimal subpaving with root box hull
*/
template<typename T>
T *Regularize(ivector& hull, ImageList& ivectorList, double eps)
{
    T* newNode = NULL;  // for return value

    /*sort the list: volCompare makes the sort largest to smallest
    Jaulin et al do not have this step because they have their own
    IMAGELIST class which acts like a set and keeps contents in order
    But we are using the stl std::list and so it is unsorted when
    it is passed to Regularize.  It is more effient to sort it once
    per call to Regularise than to keep it sorted as it is
    being built because the sorted order is only needed when
    the entire list has been built.
    */
    try {

        //sort using the volCompare function
        ivectorList.sort(volCompare);   // sorts smallest to largest

        // test if hull is contained in the first (largest) box in list

        int maxdiamcomp = 0;  // to take value calculated from MaxDiam

        // find the maximum diameter and
        double maxDiamHull = MaxDiam(hull, maxdiamcomp);

        // test if hull is equal to the largest image element, ie the last one
        bool isHullEqual = (hull == (*ivectorList.rbegin()));

        // test if hull is smaller than eps
        bool isHullSmall = (maxDiamHull < eps);

        // if the list has some images in it
        // and either if the hull is equal to the largest box in the list
        // or if the hull max diameter is < eps
        // return a new node based on hull
        if (!(ivectorList.empty()) && (isHullEqual || isHullSmall)) {
                newNode = new T(hull);
        }

        // if the list has some images in it
        // and the hull is not equal to the largest box in the list
        // and the hull max diameter is not < eps
        // return look at the left and right boxes
        if (!(ivectorList.empty()) && !isHullEqual && !isHullSmall) {

            // new ivectors from splitting hull along its biggest dimension
            ivector lefthull = Lower(hull, maxdiamcomp);
            ivector righthull = Upper(hull, maxdiamcomp);

            // create two empty lists for the left and right side
            ImageList leftlist, rightlist;

            ImageListItr it; // iterator to for the list

            // iterate through the current list and put the intersection of any
            // element with the lefthull into new left list, and the intersection
            // of any element with the new right hull into the new right list
            for (it=ivectorList.begin(); it!=ivectorList.end(); it++) {
                ivector interLeft;  // intersection with left hull
                ivector interRight;  // intersection with right hull

                if (Intersection(interLeft, *it, lefthull)) {
                    leftlist.push_back(interLeft);
                }

                if (Intersection(interRight, *it, righthull)) {
                    rightlist.push_back(interRight);
                }

            } // end of iteration through list elements

            // recursively call Regularize with lefthull, leftlist
            // and righthull, rightlist
            // reunite the results using hull as the box for parent node
            // Regularize creates a minimal subpaving
            // (no sibling child nodes) on the hull

            newNode = Reunite<T>(Regularize<T>(lefthull, leftlist, eps),
                                Regularize<T>(righthull,
                                            rightlist, eps), hull);

        } // end of is list has elements and first box does not contain hull
            // and hull is large enough to warrent further splitting

        // if there is nothing in the list we return the default
            // initialisation value of NULL
    }
    catch (bad_alloc& ba)
    {
        string msg(ba.what());
        std::cout << "Error allocating memory in Regularize" << std::endl;
        std::cout << msg << std::endl;
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std::cout << "SPnodeExcepton in Regularize: original error "
                                            << msg << endl;
    }
    catch (exception& e) {
        string msg(e.what());
        std::cout << "Error in Regularize: original error " << msg << endl;
    }

    return newNode;

}

/*! \brief Forms a non-minimal image subpaving

Make a non-minimal subpaving from a list of interval vector images
The root of the subpaving will have Box = hull, where hull has
already been formed from the union of all the ivectors in the ivectorList.
RegularizeNonMinimal is applied recursively on bisected half of hull and
new lists until either there are no images in the list or the diameter of
the hull is below eps.

Uses Reunite() and recursive calls to RegularizeNonMinimal()
to work upwards to form a non-minimal subpaving.
\param hull the interval hull of all the interval vectors in the image list
\param ivectorList a collection of possibly overlapping interval vectors to
be made into a regular subpaving
\param eps the precision with which the returned subpaving should be formed
\return a regular non-minimal subpaving with root box hull
*/
template<typename T>
T *RegularizeNonMinimal(ivector& hull, ImageList& ivectorList, double eps)
{
    T* newNode = NULL;  // for return value

    try {

        // sort the list: volCompare makes the sort largest to smallest
        // Jaulin et al do not have this step because they have their own
        // IMAGELIST class which acts like a set and keeps contents in order
        // But we are using the stl std::list and so it is unsorted when
        // it is passed to Regularize.  It is more efficient to sort it once
        // per call to Regularise than to keep it sorted as it is
        // being built because the sorted order is only needed when
        // the entire list has been built.

        // sort using the volCompare function
        ivectorList.sort(volCompare);   // sorts smallest to largest

        // test if hull is contained in the first (largest) box in the list

        int maxdiamcomp = 0;  // to take value from MaxDiam

        // find the maximum diameter
        double maxDiamHull = MaxDiam(hull, maxdiamcomp);

        // test if hull is equal to the
            // largest image element, ie the last one
        bool isHullEqual = (hull==(*ivectorList.rbegin()));
        bool isHullSmall = (maxDiamHull < eps);

        // if the list has some images in it
        // and either if the hull is equal to the largest box in the list
        // or if the hull max diameter is < eps
        // return a new node based on hull
        if (!(ivectorList.empty()) && (isHullEqual || isHullSmall)) {
            newNode = new T(hull);
        }

        // if the list has some images in it
        // and the hull is not contained in the first box in the list
        // and the hull max diameter is not < eps
        // return look at the left and right boxes
        if (!(ivectorList.empty()) && !isHullEqual && !isHullSmall) {

            // new ivectors from splitting hull along its biggest dimension
            ivector lefthull = Lower(hull, maxdiamcomp);
            ivector righthull = Upper(hull, maxdiamcomp);

            // create two empty lists for the left and right side
            ImageList leftlist, rightlist;

            ImageListItr it; // iterator to for the list

            // iterate through the current list and put the intersection of any
            // element with the lefthull into new left list, and the intersection
            // of any element with the new right hull into the new right list
            for (it=ivectorList.begin(); it!=ivectorList.end(); it++) {
                ivector interLeft;  // intersection with left hull
                ivector interRight;  // intersection with right hull

                if (Intersection(interLeft, *it, lefthull)) {
                    leftlist.push_back(interLeft);
                }

                if (Intersection(interRight, *it, righthull)) {
                    rightlist.push_back(interRight);
                }

            } // end of iteration through list elements

            // recursively call RegularizeNonMinimal with lefthull,
            // leftlist and righthull, rightlist
            // adopt the results using hull as the box for parent node
            // RegularizeNonMinimal creates a non-minimal subpaving
            // (ie has sibling child nodes) on the hull

            newNode = Adopt<T>(RegularizeNonMinimal<T>(lefthull, leftlist,
                                                    eps),
                            RegularizeNonMinimal<T>(righthull,
                                                    rightlist,
                                                    eps),
                                            hull);

        } // end of is list has elements and first box does not contain hull
        // and hull is large enough to warrent further splitting

        // if there is nothing in the list we return the default
            // initialisation value of NULL
    }
    catch (bad_alloc& ba)
    {
        string msg(ba.what());
        std::cout << "Error allocating memory in SRegularizeNonMinimal"
                                            << std::endl;
        std::cout << msg << std::endl;
    }
    catch (SPnodeException& spe) {
        string msg(spe.what());
        std:: cout << "SPnodeExcepton in RegularizeNonMinimal: original error "
                                            << msg << endl;
    }
    catch (exception& e) {
        string msg(e.what());
        std:: cout << "Error in RegularizeNonMinimal: original error "
                                            << msg << endl;
    }

    return newNode;

}

} // end namespace subpavings

#endif

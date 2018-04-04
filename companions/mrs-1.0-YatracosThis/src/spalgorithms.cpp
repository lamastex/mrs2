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

/*!/ \file:     spalgorithms.cpp
\brief SPnode (SubPaving) and algorithm function definitions
*/

#include "spalgorithms.hpp"

// include fstream so as to be able to output a file
#include <fstream>

// to use toolz methods
#include "toolz.hpp"

// spnode headers
#include "sptools.hpp"
#include "spnode.hpp"
#include "sptemplates.hpp"

//src_trunk_0701
//#include "spminimalnode.hpp"

using namespace std;

namespace subpavings {


    // Mince
    // Transforms a minimal subpaving into a non-minimal subpaving
    // ie may have sibling leaves
    // Any leaf subpaving with box with diameter > eps will be expanded
    // Mince will keep mincing until every leaf subpaving
    // has a box with diameter < eps
    void Mince(SPnode * const spn, double eps)
    {
        try {
            // only try to expand if we have a non-empty leaf
            if (!isEmpty(spn) && isLeaf(spn)) {

                int maxdiamcomp; // valued with MaxDiam function below

                // if the leaf's box is >= to eps, keep expanding
                if(MaxDiam(spn->getBox(), maxdiamcomp) >= eps) {
                    // comp is the dimension to split on

                    // if leaf and box not < eps then expand
                    Expand(spn, maxdiamcomp);
                }

            } // end !isEmpty() && isLeaf()

            // not a leaf, so Mince the children
            if (!isEmpty(spn) && !isLeaf(spn)) {
                Mince(spn->getRightChild(), eps);
                Mince(spn->getLeftChild(), eps);
            }

            // if spn points to an empty subpaving, nothing happens
        }
        catch (bad_alloc&)
        {
            std::cout << "Error allocating memory "
                << "in Mince(...)" << std::endl;
            throw;
        }

    }



    // --------------------- implementing SIVIA---------------------
    /* Set Inversion Via Interval Analysis method taken from
    Jaulin et al., pringer 2001.

    SIVIA is not templatised because it only makes sense for basic
    SPnodes, but it uses the templatised function Reunite
    */
    SPnode* Sivia (PIBT BoolTest, const SPnode * const toInvert,
                SPnode * const search, const double eps)
    {
        SPnode* newNode = NULL;  // for return value

        try {
            BOOL_INTERVAL test = BI_FALSE;

            if (!isEmpty(search)) { // if search is not null or empty

                // test the box of the given searchsubpaving
                // using given test (gives f) and given image to invert
                test = BoolTest(search->getBox(), toInvert);
            }

            if (!isEmpty(search) && test!=BI_FALSE) {

                int maxdiamcomp; // for MaxDiam() below

                double boxMaxDiam = MaxDiam(search->getBox(), maxdiamcomp);

                // we know that the test was not BI_FALSE,
                // so it could be BI_TRUE or BI_INDET
                // if it is BI_INDET and the box maximum diameter
                // is >= eps then we keep trying to expand
                if (test==BI_TRUE || boxMaxDiam < eps) {

                    newNode = new SPnode(*search);
                }

                // if test is BI_INDET and the box maximum diameter
                // is >= eps then we keep trying to expand
                else  {

                    // expand search if search is a leaf
                    if (isLeaf(search)) Expand(search,
                                            maxdiamcomp);

                    // ReUnite is used to get a minimal subpaving
                    // from merging two subpavings
                    newNode = Reunite<SPnode>(
                                Sivia(BoolTest, toInvert,
                                search->getLeftChild(), eps),
                            Sivia(BoolTest, toInvert,
                            search->getRightChild(), eps),
                                            search->getBox());
                }


            } // end !isEmpty(search) && test!=BI_FALSE

            // if isEmpty(search) or test==BI_FALSE,
            // newNode will be the initialisation value of NULL
        }
        catch (bad_alloc& ba)
        {
            string msg(ba.what());
            std::cout << "Error allocating memory in Sivia" << std::endl;
            std::cout << msg << std::endl;
        }
        catch (SPnodeException& spe) {
            string msg(spe.what());
            std:: cout << "SPnodeExcepton in Sivia: original error "
                                                << msg << endl;
        }
        catch (exception& e) {
            string msg(e.what());
            std:: cout << "Error in Sivia: original error " << msg << endl;
        }

        return newNode;
    }


    // Function to evaluate all the boxes in subpaving
    // builds an image list and
    // convex hull of images of parts of the subpaving
    void Evaluate(PIVF f, const SPnode * const spn,
                ImageList& evalImages, ivector& hull)
    {
        // later on if we want to generalise this,
        // we might have a base Box class and then
        // use list<Box>& evalImages and Box& hull

        if (spn!=NULL && isLeaf(spn)) {
            // get image using PIVF function f on box of this node
            ivector image = f(spn->getBox());

            // if no images in image set yet, make hull the image
            // if are images in image set, hull is convex hull of
            // the current hull and ivector image from f(Box(A))
            if (evalImages.size() == 0) hull = image;
            else hull = (hull | image);

            // add the image to the list of images
            evalImages.push_back(image);
        } // end of is a leaf

        // recurse on children
        if (spn!=NULL && !isLeaf(spn)) {

            Evaluate(f, spn->getLeftChild(), evalImages, hull);
            Evaluate(f, spn->getRightChild(), evalImages, hull);

        } // end of if is not a leaf

        // case where A == NULL does nothing, just returns

        return;

    }

    // Get a outer paving for the image of spn under f
    // Image evaluation
    // ImageSp is not templatised because it only makes sense for
    // basic SPnodes, but it uses the templatised function Regularize

    SPnode* ImageSp(PIVF f, SPnode *spn, double eps)
    {
        ImageList images;
        ivector hull;

        try {

            Mince(spn, eps);

            Evaluate(f, spn, images, hull);

            /* the output of eval is not included in the AIA examples,
            but it makes an interesting comparison to final subpaving */
            ofstream os2("eval.txt");            // Filename
            list<ivector>::iterator it;
            for (it=images.begin(); it!=images.end(); it++) {
                ivector box = *it;
                os2 << "[ " << Inf(box[1]) << " , " << Sup(box[1])
                    << " ] , [ " << Inf(box[2]) << " , "
                    << Sup(box[2]) << " ]" <<  endl;
            }
            // end of difference from AIA examples
        }
        catch (bad_alloc& ba)
        {
            string msg(ba.what());
            std::cout << "Error allocating memory in ImageSp" << std::endl;
            std::cout << msg << std::endl;
        }
        catch (SPnodeException& spe) {
            string msg(spe.what());
            std:: cout << "SPnodeExcepton in ImageSp: original error "
                                                << msg << endl;
        }
        catch (exception& e) {
            string msg(e.what());
            std:: cout << "Error in ImageSp: original error " << msg << endl;
        }

        // make a minimal subpaving out of images, with root box hull
        return (Regularize<SPnode>(hull, images, eps));

    }

    // Get a non-minimal outer paving for the image of spn under f
    // Image evaluation
    // ImageSpNonMinimal is not templatised because it only makes sense
    // for basic SPnodes, but it uses the templatised function
    // RegularizeNonMinimal
    SPnode* ImageSpNonMinimal(PIVF f, SPnode *spn, double eps)
    {
        ImageList images;
        ivector hull;

        try {

            Mince(spn, eps);

            Evaluate(f, spn, images, hull);
        }
        catch (bad_alloc& ba)
        {
            string msg(ba.what());
            std::cout << "Error allocating memory in ImageSpNonMinimal"
                                                << std::endl;
            std::cout << msg << std::endl;
        }
        catch (SPnodeException& spe) {
            string msg(spe.what());
            std::cout << "SPnodeExcepton in ImageSpNonMinimal: original error "
                                                << msg << endl;
        }
        catch (exception& e) {
            string msg(e.what());
            std:: cout << "Error in ImageSpNonMinimal: original error "
                                                << msg << endl;
        }

        // make a non-minimal subpaving from images,
        // with root box hull
        return (RegularizeNonMinimal<SPnode>(hull, images, eps));

    }

    // -----these functions use inheritance polymorphism

    // graft two sibling nodes onto a leaf SubPaving
    // comp is the dimension to split on
    void Expand(SPnode * const spn, int comp)
    {
        try {

            // uses nodeExpand for type of object pointed to by spn
            if (spn!=NULL) spn->nodeExpand(comp);
        }
        catch (bad_alloc&)
        {
            std::cout << "Error allocating memory in Expand"
                                    << std::endl;
            throw;
        }
    }

    // graft two sibling nodes onto a leaf SubPaving
    // finds its own dimension to split on
    void Expand(SPnode * const spn)
    {
        try {
            // uses nodeExpand() for type of object pointed to by spn
            if (spn!=NULL) spn->nodeExpand();
        }
        catch (bad_alloc&)
        {
            std::cout << "Error allocating memory in Expand" << std::endl;
            throw;
        }
    }

   // graft two sibling nodes onto a leaf SubPaving
    // finds its own dimension to split on
    // brings validation data along
    void ExpandWithValid(SPnode * const spn, bool boolVal)
    {
        try {
			    
            // uses nodeExpand for type of object pointed to by spn
            if (spn!=NULL) {
              spn->nodeExpand(boolVal);
            }
         }
        catch (bad_alloc&)
        {
            std::cout << "Error allocating memory in ExpandWithValid" << std::endl;
            throw;
        }
    }

} // end namespace subpavings


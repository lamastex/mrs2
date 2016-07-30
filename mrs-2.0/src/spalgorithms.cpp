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

/*!/ \file     spalgorithms.cpp
\brief SPnode (SubPaving) and SPMinimalnode (MinimalSubPaving) 
algorithm function definitions.
*/

#include "spalgorithms.hpp"

// spnode headers
#include "sptools.hpp"
#include "spminimalnode.hpp"
#include "sptemplates.hpp"

// to use toolz methods
#include "toolz.hpp"

// include fstream so as to be able to output a file
#include <fstream>

using namespace std;
using namespace subpavings;

// Mince
// Transforms a subpaving with one box into a non-minimal subpaving
// ie may have sibling leaves
// Any leaf subpaving with box with diameter > eps will be expanded
// Mince will keep mincing until every leaf subpaving
// has a box with diameter < eps
void subpavings::Mince(SPnode * const spn, double eps)
{
	
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



// --------------------- implementing SIVIA---------------------
/* Set Inversion Via Interval Analysis method taken from
Jaulin et al., pringer 2001.

SIVIA is not templatised because it only makes sense for basic
SPMinimalnodes, but it uses the templatised function Reunite
*/
SPMinimalnode* subpavings::Sivia (PIBT BoolTest, 
			const SPMinimalnode * const toInvert,
			SPMinimalnode * const search, const double eps)
{
	SPMinimalnode* newNode = NULL;  // for return value

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

				newNode = new SPMinimalnode(*search);
			}

			// if test is BI_INDET and the box maximum diameter
			// is >= eps then we keep trying to expand
			else  {

				// expand search if search is a leaf
				if (isLeaf(search)) Expand(search,
										maxdiamcomp);

				// ReUnite is used to get a minimal subpaving
				// from merging two subpavings
				newNode = Reunite<SPMinimalnode>(
							Sivia(BoolTest, toInvert,
							search->getLeftChild(), eps),
						Sivia(BoolTest, toInvert,
						search->getRightChild(), eps),
										search->getBox());
			}


		} // end !isEmpty(search) && test!=BI_FALSE

		// if isEmpty(search) or test==BI_FALSE,
		// newNode will be the initialisation value of NULL
		
		return newNode;
	}
	catch (std::exception const& e) {
		delete newNode;
		newNode = NULL;
		throw;
	}
}


// Function to evaluate all the boxes in subpaving
// builds an image list and
// convex hull of images of parts of the subpaving
void subpavings::Evaluate(PIVF f, const SPnode * const spn,
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
// basic SPMinimalnodes, but it uses the templatised function Regularize

SPMinimalnode* subpavings::ImageSp(PIVF f, SPMinimalnode *spn, double eps)
{
	SPMinimalnode* newNode = NULL;
	
	try {
		
		ImageList images;
		ivector hull;

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
	
		// make a minimal subpaving out of images, with root box hull
		newNode = (Regularize<SPMinimalnode>(hull, images, eps));
		
		return newNode;
	}
	catch (std::exception const& e) {
		delete newNode;
		newNode = NULL;
		throw;
	}
}


// Get a non-minimal outer paving for the image of spn under f
// Image evaluation
// ImageSpNonMinimal is not templatised because it only makes sense
// for basic SPnodes, but it uses the templatised function
// RegularizeNonMinimal
SPMinimalnode* subpavings::ImageSpNonMinimal(PIVF f, 
										SPMinimalnode *spn, double eps)
{
	SPMinimalnode* newNode = NULL;
	try {
		ImageList images;
		ivector hull;

		Mince(spn, eps);

		Evaluate(f, spn, images, hull);

		// make a non-minimal subpaving from images,
		// with root box hull
		newNode = (RegularizeNonMinimal<SPMinimalnode>(hull, images, eps));
		
		return newNode;
	}
	catch (std::exception const& e) {
		delete newNode;
		newNode = NULL;
		throw;
	}
}

// -----these functions use inheritance polymorphism

// graft two sibling nodes onto a leaf SubPaving
// comp is the dimension to split on
void subpavings::Expand(SPnode * const spn, int comp)
{
	if (spn!=NULL) spn->nodeExpand(comp);
	
}

// graft two sibling nodes onto a leaf SubPaving
// finds its own dimension to split on
void subpavings::Expand(SPnode * const spn)
{
	if (spn!=NULL) spn->nodeExpand();
	
}



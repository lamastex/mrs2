/*
 * Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
 * Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief AIASPnode (AISSubPaving) definitions
*/
#include "AIAsubpaving.hpp"

using namespace std;

// -------------------------- AIASPnode class definitions
//Copy constructor
AIASPnode::AIASPnode(const AIASPnode& n)
{
  theBox=new ivector(*n.theBox);
  //recursion on the children
  if (n.leftChild) leftChild=new AIASPnode(*n.leftChild);
  else leftChild=NULL;
  if (n.rightChild) rightChild=new AIASPnode(*n.rightChild);
  else rightChild=NULL;
}

//Output the leaf boxes in AIASubPaving a
std::ostream & operator<< (std::ostream & os, AIASubPaving a)
{
  if (a==NULL) return os;

  //if (IsLeaf(a)) os << (*(a->theBox)) << std::endl;
  //output for Jenny's matlab format, dimension = 2 max
  if (IsLeaf(a))
  {
    ivector x = *(a->theBox);
    os << "[ " << Inf(x[1]) << " , " << Sup(x[1]) << " ] , "
      << "[ " << Inf(x[2]) << " , " << Sup(x[2]) << " ]"   << std::endl;
  }
  else
    { os << (a->leftChild); os << (a->rightChild); }

    return os;
}

//Number of leaf boxes in AIASubPaving a
//eturns 0 if empty (IsEmpty), 1 if a leaf (IsLeaf), else number of leaves
//Recurses to accumulate number of leaves on child nodes
int NbLeaves (AIASubPaving a )
{
  int nbleaves=0;

  if (IsEmpty(a)) return 0;
  if (IsLeaf(a))  return 1;

  nbleaves += NbLeaves(a->leftChild);
  nbleaves += NbLeaves(a->rightChild);

  return (nbleaves);
}

//Sum of volumes of all the leaf boxes in AIASubPaving a
//Returns 0 if empty (IsEmpty)
//Recurses to accumulate volume of child nodes
double Volume (AIASubPaving a)
{

  if (IsEmpty(a)) return 0.0;
  if (IsLeaf(a))
  {

                    // using Volume taking ivector argument
    return Volume(Box(a));
  }
  double vol=0.0;

  vol += Volume(a->leftChild);
  vol += Volume(a->rightChild);

  return (vol);
}

// check for containment of ivector or box in the AIASubPaving
// returns false if X is empty
// recursively calls <= on child nodes
AIA_BOOL_INTERVAL operator<=(const ivector& z, AIASubPaving X)
{
  // z is assumed not to be empty
  // note that Intersection() gives an error if unequal index sets passed

  if (IsEmpty(X)) return BI_FALSE;

  ivector r;        // ivector to be passed to Intersection()
  // to contain the intersection

  // check if X is a leaf and if so
  if (IsLeaf(X))
  {
    // *************** old code
    //if (z<=Box(X)) return BI_TRUE;
    //if (!Intersection(r,z,Box(X))) return BI_FALSE;
    // what if Intersection() returns 1, ie there is an intersection
    // but z<=Box(X) is not true?
    // at present this goes down to the default (else) case below

    //**************** replacement code
    if (z<=Box(X))
                    // true if z is in the box
        return BI_TRUE;
    else if (!Intersection(r,z,Box(X)))
                    // false if there is no intersection between z and the box
        return BI_FALSE;
    else          // z is not contained in the box but there is an intersection
                  // so result is indeterminate
        return BI_INDET;
  }                 // end isLeaf()

  AIA_BOOL_INTERVAL Ltest=BI_TRUE,Rtest=BI_TRUE;

  ivector Lz,Rz;    // ivectors to be passed to Intersection()
  // will contain intersection after call to Intersection()

                    // both children exist
  if (!IsEmpty(X->leftChild)&&!IsEmpty(X->rightChild))
  {

                    // intersection with left child
    if (Intersection(Lz,z,Box(X->leftChild)))
      // Lz now contains this intersection
    {
                    // compare Lz and left child
      Ltest = (Lz<=(X->leftChild));

                    //check intersection with right child
      if (Intersection(Rz,z,Box(X->rightChild)))
      {
                    // compare Rz and right child
        Rtest = (Rz<=(X->rightChild));
                    // if left and right child tests same
        if (Ltest==Rtest)
                    // return this result from child tests
          return Ltest;
        else
                    // z overlaps the boundary of X
          return BI_INDET;
      }
      else
        return Ltest;
    }               // end if intersection with left child

                    // given no intersection with left child
    else if (Intersection(Rz,z,Box(X->rightChild)))
                    // check intersection with right child
                    // compare Rz and right child
      return (Rz<=(X->rightChild));

    else            // no intersection with either child
                    // used to be BI_TRUE here from the AIA website
        return BI_FALSE;
    // but this cannot possibly be correct!
  }                 // end if both children exist

                    // only left child exists
  else if (!IsEmpty(X->leftChild))
  {
                    // no intersection with left child
    if (!Intersection(Lz,z,Box(X->leftChild)))
      return BI_FALSE;
    // there is an intersection with the left child, this intersection now in Lz
                    // compare Lz and the left child
    Ltest = (Lz<=(X->leftChild));
    // Rtest will still be at default value of TRUE

                    // if left child intersection with z is not z
    if (!(Lz==z)) Rtest = BI_FALSE;
    else            // if left child intersection with z = z
      return Ltest; // ie is whole of z in left child
  }                 // end if  only left child exists

  else              // if is not a leaf, only right child exists
  // BUT if it is a leaf and there is an intersection but not full containment?
  // I have now replaced the code that gave rise to the latter case
  {
                    // no intersection with right child
    if (!Intersection(Rz,z,Box(X->rightChild)))
      return BI_FALSE;
    //there is an intersection with the right child, this intersection now in Rz
                    // compare Rz and right child
    Rtest = (Rz<=(X->rightChild));
    // Ltest will still be at default value of TRUE
                    // if right child intersection with z is not z
    if (!(Rz==z)) Ltest = BI_FALSE;
    else            // if right child intersection with z = z
      return Rtest; // ie is whole of z in right child
  }                 // end if only right child exists

  // only get here if only one child (left or right exists)
  // and the intersection of that child with z exists but is not z
  // in which case the other test will have been set to FALSE?
  // or if is a leaf and both Intersection() returns 1 but z<=Box(X) is not 1
  // in which case Ltest and Rtest will both be TRUE? - that code now replaced

  if (Ltest==Rtest)
    return Ltest;
  else
    // if only one child, and the whole of z is not in that child but
    // the intersection of z with that child is inside the child
    return BI_INDET;// z overlaps the boundary of X
}

// --------------------------  Implementing SIVIA components -------------------

// Expand
//graft two sibling nodes onto an AIASubPaving node provided that node is a leaf
void Expand(AIASubPaving A, int comp)
{
  // comp argument is passed to Upper() and Lower()
  // these functions split the box normal to direction set by comp

  if (!IsLeaf(A)) return;
  ivector lC,rC;
  Lower(Box(A),lC,comp); Upper(Box(A), rC, comp);
  A->leftChild = new AIASPnode(lC);
  A->rightChild = new AIASPnode(rC);
  //A->leftChild = new AIASPnode(Lower(Box(A),comp));
  //A->rightChild = new AIASPnode(Upper(Box(A),comp));
}

//AIASubPaving ReUnite(AIASubPaving lChild, AIASubPaving rChild, ivector&
// computes a minimal AIASubPaving from two sibling AIASubPavings
// an AIASubPaving is minimal if it has no sibling leaves
// a minimal AIASubPaving is created by discarding sibling leaves so that
// their parents become leaves
// passing by value since AIASPnode() requires so...
AIASubPaving ReUnite(AIASubPaving lChild, AIASubPaving rChild, ivector x)
// lChild and rChild are the two AIASubPavings to be reunited
// ivector x is the box for the new AIASubPaving to be created
{
  if(IsEmpty(lChild)&&IsEmpty(rChild)) return NULL;
  AIASubPaving result = new AIASPnode(x);
  // both AIASubPavings are leaves so discard them: new AIASubPaving is a leaf
  if( IsLeaf(lChild)
    &&IsLeaf(rChild)
    &&(  x == ( Box(lChild) | Box(rChild) )  )
    )
    { delete lChild; delete rChild; return result; }

    // if there are not two non-null potential children
    // (so presumably just one child), just graft it on
    // similarly if at least one child is not a leaf, just graft the potential
    // children on
    result->leftChild = lChild; result->rightChild = rChild;

  return result;
}

// ----------------------------------- implementing SIVIA-----------------------
/*! Set Inversion Via Interval Analysis method taken from AIA2001.  Sivia gives
    an AIASubPaving which is based on the reciprocal image X of some function.
    Note that in this implementation, we get the outer AIASubPaving of two,
    (lowerX and UpperX) such that lowerX is in X is in upperX. The AIA_PIBT
    points to a function which specifies the function for which we want the
    reciprocal image and the AIASubPaving Y which should contain the reciprocal
    image.  The AIA_PIBT function evaluates whether the image of some box
    (ivector) is in Y.  The AIA_PIBT can return BI_TRUE, BI_FALSE, or BI_INDET
    (box partly in Y but not totally).  If the AIASubPaving box tests BI_INDET
    and the box diameter is sufficiently small then it is given the benefit of
    the doubt and included in upperX.  Otherwise if the AIA_PIBT function
    returns BI_INDET then the AIASubPaving is expanded, ie two AIASubPaving
    children created, and these AIASubPaving children are then tested.  The
    argument A is an AIASubPaving.  On the first pass through SIVIA, A should
    have a box which is guaranteed to contain upperX; this AIASubPaving is
    progressively refined.  The argument eps specifies the 'sufficiently
    small' width for a box to be included in upperX...
*/
AIASubPaving Sivia (AIA_PIBT BoolTest, AIASubPaving A, double eps)
{
  if (A==NULL) return NULL;

                    // test the box of the given AIASubPaving
  AIA_BOOL_INTERVAL test = BoolTest(Box(A));
  // Test function (may be passed as parameter to Sivia)

  // maxdiamcomp will be given a value by call to MaxDiam() below
  int maxdiamcomp;

  // the box fails the test
  if ( test==BI_FALSE ) return NULL;

  // if the box passes or the result is BI_INDET but the maximum diameter of
  // the box is small enough, then return a new AIASubPaving which is a copy of
  // the current one, ie include this AIASubPaving in outerX, our approximation
  // of X the reciprocal image
  if (test==BI_TRUE || MaxDiam(Box(A),maxdiamcomp)<eps)
    return new AIASPnode(*A);

  // SIVIA will only reach this point if the result of the test was BI_INDET
  // and the maximum diameter of the box is not small enough for the box to be
  // in outerX -- in this case we expand the AIASubPaving by giving it child
  // nodes and test them
  if (IsLeaf(A)) Expand(A,maxdiamcomp);

  // ReUnite is used to get a minimal AIASubPaving from merging two
  // AIASubPavings.  So will ensure that the AIASubPaving we return from Sivia
  // is minimal, ie will not have sibling leaves
  return ReUnite( Sivia(BoolTest,A->leftChild,eps),
    Sivia(BoolTest,A->rightChild,eps),
    Box(A));

}

// ----------------------  implementing ImageSp components ---------------------
// Mince
/*! Mince transforms a minimal AIASubPaving into a non-minimal AIASubPaving,
ie may have sibling leaves.  Any leaf AIASubPaving with a box with diameter
greater than eps will be expanded.  Mince will keep mincing until every leaf
AIASubPaving has a box with diameter < eps
*/
void Mince(AIASubPaving A, double eps)
{
  if (IsEmpty(A)) return;

  if (IsLeaf(A))
  {
    int comp;       // value is given by calling MaxDiam function below

    // if leaf and box smaller than eps then return
    if(MaxDiam(Box(A),comp)<eps) return;

    // if leaf and box not smaller than eps then expand
    Expand(A,comp);

  }                 // end if a leaf

  // not a leaf, recurse Mince() on left and right children
  Mince(A->rightChild,eps);
  Mince(A->leftChild,eps);

}

// Evaluate
/*! Function to evaluate all the boxes in AIASubPaving.  The images of all the
boxes are formed and put into a list and the interval hull of the union of all
the boxes is formed.  Interval hull of a union is implemented under cxsc
with | operator.  Boxes are expected to have been created by Mince().  The
argument f is a pointer to an interval vector function AIA_PIVF (which returns
an interval vector).  AIASubPaving A is the AIASubPaving to be evaluated
(should be const - it is not changed in this function).  evalImages is the
current list of images to be added to in the function call, hull is an ivector
hull of images currently in the list which may be altered in this call if an
image is added to list
*/
void Evaluate(AIA_PIVF f, AIASubPaving A, list<ivector>& evalImages, ivector& hull)
{
  if (A!=NULL && IsLeaf(A))
  {
    // make an ivector image using the AIA_PIVF function f on Box(A)
    ivector image = f(Box(A));

    // if no images in image set yet, make hull the image
    if (evalImages.size() == 0) hull = image;
    // if there are images in the image set, hull is the convex hull of
    else hull = (hull | image);
    // the current hull and the ivector image from f(Box(A))

    // add the image to the list of images
    evalImages.push_back(image);

  }                 // end of is a leaf

  // if not a leaf, recursively call Evaluate on children
  if (A!=NULL && !IsLeaf(A))
  {

    Evaluate(f, A->leftChild, evalImages, hull);
    Evaluate(f, A->rightChild, evalImages, hull);

  }                 // end of if is not a leaf

  // case where A == NULL does nothing, just returns
  return;
}

// return TRUE if volume of a < volume of b
bool volCompare(const ivector &a, const ivector &b)
{
  bool returnValue = 0;

  // Make sure the vectors have the same number of elements and at
  // least one element each
  if( (Ub(a) - Lb(a)) == (Ub(b) - Lb(b)) && (Ub(a) - Lb(a))>=1 )
  {

                    // compare the two volumes
    returnValue = ((Volume(a)<Volume(b)));

  }
  else
  {
    std::cout
      << "Error in volCompare : comparing ivectors of different dimensions"
      << std::endl;
  }

  return returnValue;
}

// Regularize
/*! Regularize creates an AIASubPaving from all the ivectors in the
ivectorList.  The root of the AIASubPaving will have Box = hull where hull has
already been formed from the union of all the ivectors in the ivectorList.
Regularize is applied recursively on bisected half of hull and new lists until
either there are no images in the list or the diameter of the hull is below eps.
*/
AIASubPaving Regularize(ivector& hull, list<ivector>& ivectorList, double eps)
{
                    // return NULL if the list is empty
  if (ivectorList.size()==0) return NULL;

  // sort the list
  // Jaulin et al do not have this step because they have their own IMAGELIST
  // class which acts like a set and keeps the contents in order.  But we are
  // using the stl std::list and so it is unsorted when it is passed to
  // Regularize.  It is more effient to sort it once per call to Regularise
  // than to keep it sorted as it is being built because the sorted order is
  // only needed when the entire list has been built

  // sort using the volCompare function
  ivectorList.sort(volCompare);
  // this sorts smallest to largest (the opposite to Jaulin et al)

  // if the hull is equal to the last (largest) box in the list,
  // this becomes the AIASubPaving
  if (hull==(*ivectorList.rbegin())) return new AIASPnode(hull);

  // un-valued int to take value for larged dimension calculated from MaxDiam
  int maxdiamcomp;

  //if the current maximum diameter is < eps return a new AIASubPaving from hull
  if (MaxDiam(hull,maxdiamcomp)<eps) return new AIASPnode(hull);

  // new ivectors from splitting hull in its biggest dimension
  ivector lefthull = Lower(hull,maxdiamcomp);
  ivector righthull = Upper(hull,maxdiamcomp);

  // create two empty lists
  list<ivector> leftlist,rightlist;

  // iterator to for the list
  list<ivector>::iterator it;

  // iterate through the current list and put the intersection of any element
  // with the lefthull into the new left list, and the intersection of any
  // element with new right hull into the new rightlist.
  for (it=ivectorList.begin(); it!=ivectorList.end(); it++)
  {
    // temporary variables to take the results of call in Intersect
    ivector interLeft, interRight;

    if (Intersection(interLeft, *it, lefthull))
    {
      leftlist.push_back(interLeft);
    }

    if (Intersection(interRight, *it, righthull))
    {
      rightlist.push_back(interRight);
    }

  }  // end of iteration through list elements

  // recursively call Regularize with lefthull, leftlist, righthull, rightlist
  // reunite the results using hull as the box for the parent node
  // Regularize creates a minimal AIASubPaving (no sibling child nodes)
  return ReUnite(Regularize(lefthull,leftlist,eps),
    Regularize(righthull,rightlist,eps),hull);
}

AIASubPaving ImageSp(AIA_PIVF f, AIASubPaving A, double eps)
{
  list<ivector> images;
  ivector hull;

  Mince(A, eps);

  //cout << "After mince " << endl;
  //cout << "A has volume  " << Volume (A) << " and number of leaves "
  //  << NbLeaves(A) << endl;

  Evaluate(f, A, images, hull);

  //cout << "After evaluate " << endl;
  //cout << "Size of image list is : " << images.size()
  //  << "  and hull has volume " << Volume(hull) << endl;

  /* the output of eval is not included in the AIA examples, but it makes
  an interesting comparison to the final subpaving */
  // Filename
  ofstream os2("eval.txt");
  list<ivector>::iterator it;
  for (it=images.begin(); it!=images.end(); it++)
  {
    ivector box = *it;
    os2 << "[ " << Inf(box[1]) << " , " << Sup(box[1]) << " ] , [ "
      << Inf(box[2]) << " , " << Sup(box[2]) << " ]" <<  endl;
  }
  // end of difference from AIA examples

  return (Regularize(hull, images, eps));

}

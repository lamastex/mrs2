/*
MRS: A C++ class library for statistical set processing

Copyright (C) 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin and Thomas York
Copyright (C) 2009, 2010, 2011, 2012, 2013 Jennifer Harlow and Raazesh Sainudiin

Copying Permission Statement:

MRS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
See http://www.gnu.org/copyleft/gpl.html
*/
/*! \file
\brief MRS-1.0 HTML document generator via Doxygen.
*/

/*!
\mainpage

\image html MRS.png "Auto-validating Samples" width=15cm
\image latex MRS.png "Auto-validating Samples" width=15cm

- \ref mainpage_sec_overview
- \ref mainpage_sec_install
- \ref mainpage_sec_examples
- \ref mainpage_sec_history
- \ref mainpage_sec_acknowledgements

<HR>

\section mainpage_sec_overview Overview

MRS is a C++ class library for statistical set processing.

The goal of MRS is to improve the accuracy and reliability of numerical results
in computational statistical problems and provide a cohesive object-oriented
framework for set-valued and set-oriented computational statistics.  This is
generally achieved by extending arithmetic beyond the built-in data types and
applying fixed-point theorems.

MRS builds on
<A HREF="http://www.math.uni-wuppertal.de/wrswt/xsc/cxsc_new.html">C-XSC</A>
 C eXtended for Scientific Computing and
<A HREF="http://www.gnu.org/software/gsl/">GSL</A>
 Gnu Scientific Library.
MRS, C-XSC and GSL are distributed under the terms of GPL, the GNU General Public
License.

<HR>

\section mainpage_sec_install Download and Install
<B> Download the source code of beta version for GNU/Linux, Linux, Unix, Mac OS X (to be released in April 2009)</B>
- <A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/mrs-0.1.2-beta505.tar.gz">
  mrs-0.1.2-beta505.tar.gz (1.1MB)</A>,
  <A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/README">
  README</A>,
  <A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/ChangeLog">
  ChangeLog</A>,
<A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/AUTHORS">
  AUTHORS</A>,
<A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/CONTENTS">
  CONTENTS</A>,
<A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/COPYING">
  COPYING</A>,
<A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/COPYRIGHT">
  COPYRIGHT</A>,
<A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/codes/mrs-0.1.2/docu-mrs-0.1.2-beta505.tar.gz">
API-Documentation (Download-Version, html-format) (18.0 MB)  </A>

<HR>

\section mainpage_sec_examples Examples

MRS library provides the following features for trans-multi-dimensional
inclusion-isotonic target functions:

- \subpage cxscexamples of interval arithmetic for computational statisticians
- Some challenging \subpage target_examples
- \subpage moorerejsam for auto-validating trans-dimensional sampling
- Importance sampler from pseudo and quasi random numbers
- Phylogenetic posterior densities over small tree spaces
- \subpage pavproc are a set of algorithms for set processing in labeled metric spaces
- \subpage AIASubPavings -- [Jaulin et al, Springer 2001]
- \subpage AdaptiveHistograms based on \subpage StatsSubPavings derived from \subpage newsubpavings
- Automatic Differentiation for Gradients, Jacobians and Hessians (via C-XSC)
- Global Optimization (multi-dimensional via C-XSC)

MRS builds upon the CXSC library that provides all features indispensable for modern numerical software development, such as

- Operator concept (user-defined operators)
- Overloading concept
- Module concept
- Dynamic arrays
- Controlled rounding
- Predefined arithmetic data types real, (extended real), complex, interval, complex interval, and corresponding vector and matrix types
- Predefined arithmetic operators of highest accuracy for the arithmetic data types
- Predefined elementary functions of highest accuracy for the arithmetic data types
- Data type dotprecision for the exact representation of dot products
- Library of mathematical problem-solving routines with automatic result verification and high accuracy

<HR>

\section mainpage_sec_history History

\version 1.0
\author 1. Raazesh Sainudiin (r.sainudiin ATADDRESSSYMBOL math.canterbury.ac.nz),
Laboratory for Mathematical Statistical Experiments and
Department of Mathematics and Statistics,
University of Canterbury,
Private Bag 4800,
Christchurch, New Zealand,

\author 2. Jennifer Harlow (jah217 ATADDRESSSYMBOL uclive.ac.nz)

\author 3. Thomas York (tly2 ATADDRESSSYMBOL cornell.edu)

The work on MRS started during Raazesh Sainudiin's PhD dissertation.
Many colleagues have contributed to the realization of MRS. See specific source files
for details.  Special thanks go to:
- Warwick Tucker (warwick ATADDRESSSYMBOL math.uu.se) -- general advise on interval methods since 2005
- Brendan Bycroft (brb44 ATADDRESSSYMBOL student.canterbury.ac.nz) -- GNU Autotools support in 2009
- Gloria Teng (glo.teng ATADDRESSSYMBOL gmail.com) -- Histogram methods 2008 - 2012

\mrscitation

<OL>

<LI>
Posterior expectation of regularly paved random histograms, Raazesh Sainudiin, Gloria Teng, Jennifer Harlow and Dominic Lee, <A HREF="http://dx.doi.org/10.1145/2414416.2414422">ACM Trans. Model. Comput. Simul. 23, 1, Article 6, 20 pages</A>, 2013 (<A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/SubPavingMCMC.pdf">PDF</A> 1864KB)

<LI>
Mapped Regular Pavings, Jennifer Harlow, Raazesh Sainudiin and Warwick Tucker, <A HREF="http://interval.louisiana.edu/reliable-computing-journal/volume-16/reliable-computing-16-pp-252-282.pdf">Reliable Computing, vol. 16, pp. 252-282</A>, 2012 (<A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/MappedRegularPaving.pdf">PDF</A> 972KB)

<LI>
Statistical regular pavings to analyze massive data of aircraft trajectories, Gloria Teng, Kenneth Kuhn and Raazesh Sainudiin, <A HREF="http://arc.aiaa.org/doi/abs/10.2514/1.I010015">Journal of Aerospace Computing, Information, and Communication, Vol. 9, No. 1, pp. 14-25</A>, 2012 (<A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/AAIASubPavingATC.ps">PS</A> 31MB or lossy <A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/AAIASubPavingATC.pdf">PDF</A> 2.9MB or <A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/AAIASubPavingATC_PNG/">26 PNG pages</A>)

<LI>
Auto-validating von Neumann Rejection Sampling from Small Phylogenetic Tree Spaces, Raazesh Sainudiin and Thomas York <A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/AutoValidPhyloSampler.pdf">(PDF 1.1MB)</A>.[older version: December 2006, arXiv/math.ST/0612819 (<A HREF="http://arxiv.org/abs/math.ST/0612819">arXiv link</A>)], Algorithms for Molecular Biology 2009, 4:1
(<A HREF="http://www.almob.org/content/4/1/1">Open Acess -- with sketchy latex support of the online versions as of April 25, 2009.</A>)

<LI>
Applications of interval methods to phylogenetics, Raazesh Sainudiin and Ruriko Yoshida, In L. Pachter and B. Sturmfels (Eds.), Algebraic Statistics for Computational Biology, Cambridge University Press, 2005, DOI: 10.2277/0521857007, <A HREF="http://www.cambridge.org/uk/catalogue/catalogue.asp?isbn=9780521857000">(Cambridge Catalogue)</A>.

<LI> Enclosing the Maximum Likelihood of the Simplest DNA Model Evolving on
Fixed Topologies: Towards a Rigorous Framework for Phylogenetic Inference,
Raazesh Sainudiin, BSCB Dept. Technical Report BU-1653-M, Cornell University,
Ithaca, New York, USA, February 2004 (revised December 2004) <A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/GOPJC69manuscript.pdf">
 (PDF 1.7MB)</A>.

<LI> Machine Interval Experiments, Raazesh Sainudiin, Ph.D. dissertation in
the field of Statistics, Cornell University, Ithaca, New York, USA, May 2005
<A HREF="http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/PHD.pdf">
(PDF 1.8MB)</A>.

<LI> An Auto-validating Rejection Sampler, Raazesh Sainudiin and Thomas York,
BSCB Dept. Technical Report BU-1661-M, Cornell University, Ithaca, New York,
USA, December 2005 <A HREF=
"http://www.math.canterbury.ac.nz/~r.sainudiin/preprints/AutoValidSampler.pdf">
(PDF 1.7MB)</A>.

</OL>

<HR>

\section mainpage_sec_acknowledgements Acknowledgements

The following sources of funding have supported past and current versions of MRS library:
- 2003 by NSF grant DGE-9870631 (United States)
- 2004, 2005 by joint NSF/NIGMS grant DMS-02-01037 (United States)
- 2006, 2007 by Research Fellowship of the Royal Commission for the Exhibition of 1851 (United Kingdom)
- 2007, 2008 by internal research grants from Mathematics and Statistics Department, University of Canterbury (New Zealand)
- 2009 by statistical consulting revenues of Raazesh Sainudiin
- 2012 University of Canterbury Postgraduate Scholarship
<HR>

\subpage GFDL

    Copyright (C)  2008, 2009, 2010, 2011, 2012, 2013 Razesh Sainudiin and Jennifer Harlow.
    Permission is granted to copy, distribute and/or modify this document
    under the terms of the GNU Free Documentation License, Version 1.3
    or any later version published by the Free Software Foundation;
    with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
    A copy of the license is included in the section entitled "GNU
    Free Documentation License".

*/

/*!
\page cxscexamples Basic C-XSC Examples

 To learn how to use the C-XSC class library and the GSL you can read the
 page \ref cxscexamples and follow the examples described there.  For a
 proper introduction to C-XSC and the accompanying Toolbox with more
 sophisticated applications see the <A HREF=
 "http://www.math.uni-wuppertal.de/wrswt/xsc/cxsc/apidoc/html/index.html">
 C-XSC Documentation</A>.

 On this page we will try to show you how to use the C-XSC class library with
 some short and simple examples.

 - \ref cxscexamples_sec_ex1
 - \ref cxscexamples_sec_ex2
 - \ref cxscexamples_sec_ex3
 - \ref cxscexamples_sec_ex4
 - \ref cxscexamples_sec_ex5

 \section cxscexamples_sec_ex1 Example 1 - Intervals in C-XSC

 In the following simple program we use C-XSC intervals to do basic arithmetic
 operations and output the results:

 \dontinclude example.cpp

 \skip include
 \until ---*

 \dontinclude example.cpp

 Let's start examining the code line by line.  The first line:

 \skipline "interval.hpp"

 includes the basic interval class of C-XSC in the program.  The second line:

 \skipline iostream

 includes the standard iostream library for basic input and output operations.
 The next two lines inform the compiler about C-XSC's namespace cxsc and the
 standard library namespace std.  The namespace cxsc is where all of C-XSC's
 classes and methods are stored.  This allows us to use C-XSC classes without
 having to fully qualify their identifiers.

 \skipline using namespace
 \until std

 Next we declare two interval variables and assign adequate values in the
 following lines.

 \skipline interval a
 \until b;

 Finally, we print out the result for our desired subtractions.

 \skipline cout
 \until cout

 To compile the program we edit the Makefile in the examples directory.  First
 we set the 'PROGRAM=example' and the PREFIX to the appropriate directory that
 contains the C-XSC includes and lib directories.  Then we type 'make all' in a
 Unix system to compile the program.

 \section cxscexamples_sec_ex2 Example 2 - Multi-precision Intervals in C-XSC

 \dontinclude lexample.cpp

 \skip include
 \until ---*


 \section cxscexamples_sec_ex3 Example 3 - Mean in GSL and Dot Precision Accumulators in C-XSC

 \section cxscexamples_sec_ex4 Example 4 - Range Enclosure with Automatic Differentiation

 \section cxscexamples_sec_ex5 Example 5 - Global Optimisation for Maximum Likelihood


*/
/*! \page target_examples Example target densities for optimisation and sampling

On this page we show some examples of using Moore Rejection Sampling with selected targets

- Cavender-Farris-Neyman Model
- \ref MRSexamsec_Levy
- \ref MRSexamsec_Rosenbrock
- Witch's Hat
- Trans-dimensional Cavender-Farris-Neyman Model

  \section MRSexamsec_Levy The function of Levy

  The function of Levy is \f$ \Re^2 \rightarrow \Re \f$ with

  \f[
    f(x,y) = \sum_{i=1}^5 i \cos((i - 1)x + i)\sum_{j=1}^5j\cos((j + 1)y + j)
    \f]

  \image html Levy_function.png "Function of Levy on a domain [-2,2]x[-2,2]" width=15cm
  \image latex Levy_function.png "Function of Levy on a domain [-2,2]x[-2,2]" width=15cm

  We convert the above function into a target shape by exponentiating its
negative that is scaled by a temperature parameter \f$T\f$ :
  \f[
    l_{T}(x) =  \exp \left( - \frac{1}{T}\sum_{i=1}^5 i \cos((i - 1)x_1 + i)
                            \sum_{j=1}^5j\cos((j + 1)x_2 + j) \right)
    \f]

                       Levy Target Shape \f$l_{40}(x_1,x_2): [-10,10]\times [-10,10] \to \Re \f$ is
  plotted below.

  \image html LevyT40TargetShape.png "Levy Target Shape on a domain [-10,10]x[-10,10]" width=15cm
  \image latex LevyT40TargetShape.png "Levy Target Shape on a domain [-10,10]x[-10,10]" width=15cm

  \section MRSexamsec_Rosenbrock Rosenbrock's Multi-dimensional function

  Rosenbrock's \f$ D \f$-dimensional function \f$ \Re^2 \rightarrow \Re \f$ with

  \f[
    f(x) = \exp \left( - \sum_{i=2}^D \left( 100(x_i - x_{i-1}^2)^2 +
                                             ( 1 - x_{i-1})^2 \right) \right)
    \f]

  \image html Rosenbrock_function.png "Rosenbrock's D=2 dimensional function on a domain [-2,2]x[-1,3]" width=15cm
  \image latex Rosenbrock_function.png "Rosenbrock's D=2 dimensional function on a domain [-2,2]x[-1,3]" width=15cm

*/

/*! \page moorerejsam Moore Rejection Sampler

The \ref moorerejsam is useful for obtaining independent and identically
distributed (IID) samples from a taget shape using the von Neumann rejection
sampler and interval analysis.


- \ref moorerejsam_sec_RSAlg
- \ref moorerejsam_sec_MRS
- \ref moorerejsam_sec_examples

  \section moorerejsam_sec_RSAlg von Neumann Rejection Sampler

  \image html "RSAlg.png" "von Neumann Rejection Sampling Algorithm" width=15cm
  \image latex "RSAlg.png" "von Neumann Rejection Sampling Algorithm" width=15cm

  \section moorerejsam_sec_MRS Moore Rejection Sampler is an Auto-validating von Neumann Rejection Sampler

  \image html "CFNCin762Cnin133.png" "Idea behind the Moore Rejection Sampling Algorithm" width=15cm
  \image latex "CFNCin762Cnin133.png" "Idea behind the Moore Rejection Sampling Algorithm" width=15cm

  \section moorerejsam_sec_examples Example Targets for Moore Rejection Sampling

- Cavender-Farris-Neyman Model
- \ref MRSexamsec_Levy
- \ref MRSexamsec_Rosenbrock
- Witch's Hat
- Trans-dimensional Cavender-Farris-Neyman Model

*/

/*! \page pavproc Pavings and subpavings

- \ref intro
- \ref setcomputation
- \ref statssetprocessing

<HR>

  \section intro A very short introduction to interval analysis and subpavings
  This is a very concise introduction to some of the concepts involved in
  %subpavings using descriptions given in the introductory text book
  "Applied Interval Analysis", by Jaulin, Kieffer, Didrit and Walter;
  Springer, 2001.  This book is abbreviated by AIA2001.

  For a full description of subpavings see
  AIA2001 http://www.lss.supelec.fr/books/intervals/

  \subsection intervals Intervals
  An interval real [x] (usually referred to simply as an interval [x]) is a
  connected subset of the reals R.  An interval has an upper bound and a lower
  bound. [AIA2001, pp. 18-20]

  \subsection intervalvectors Interval vectors and boxes
  An interval vector [<b>x</b>] is a subset of R<SUP>n</SUP>
  (n-dimensional space) that can be defined as the cartesian product of n
  closed intervals.  [<b>x</b>] is usually referred to as an interval vector
  or \e box.  IR<SUP>n</SUP> is the set of all closed intervals of
  R<SUP>n</SUP>. [AIA2001, p. 23]

  \subsection inclfuncs Inclusion functions
  Given a function <b>x</b> from R<SUP>n</SUP> to R<SUP>m</SUP>, the interval
  function [<b>f</b>] is an inclusion function for <b>f</b> if, for
  all [<b>x</b>] in IR<SUP>n</SUP>, <b>f</b>([<b>x</b>]) is a subset
  of [<b>f</b>]([<b>x</b>]).  The inclusion function <b>f</b> makes it possible
  to compute a box [<b>f</b>]([<b>x</b>]) guaranteed to contain
  f([<b>x</b>]). [AIA2001, p. 27]

  \subsection subpavs Subpavings
  A \e subpaving of a box [<b>x</b>] is a union of non-overlapping subboxes
  of [<b>x</b>]. [AIA2001, p. 48]. To form this subpaving we subdivide
  [<b>x</b>] and then subdivide the subboxes, and then subdivide those
  subboxes . . .  The definition of a subpaving also allows us to not select
  some subbox resulting from any subdivision in the subpaving
  (only the subboxes we keep can be further subdivided, of course).

  Bisection of a box is subdivision of the box in half along its longest
  dimension (or the first longest dimension if it has the same longest length
  in more than one dimension).  In one dimension, bisection is dividing an
  interval in half.  In two dimensions, bisection is drawing a line through the
  middle of a rectangle, normal to its longest dimension.  In three dimensions
  bisection is slicing a plane through a box, again normal to its longest
  dimension and again in order to divide the box in half. In general, in n
  dimensions, bisection is by an (n-1)-dimensional hyperplane normal to the
  longest dimension of the box at the midpoint of the box on that longest
  dimension.

  A \e regular subpaving is obtained from a finite succession of bisections and
  selections.  It is computationally simpler to perform operations, such as
  intersection, on regular %subpavings than on non-regular %subpavings and so
  they are more easily manipulated in computer applications
  [AIA2001, pp. 49-50].

  \subsection regsubpavs Regular subpavings and binary trees
  A regular subpaving can be represented as a binary tree.  Bisection of a box
  into two subboxes corresponds to the node representing that box branching to
  two child nodes.  A decision not to include a subbox in the subpaving means
  that that branch is pruned off the tree.  Thus the growth of the branches is
  defined by how the initial box, which corresponds to the root of the tree,
  is bisected and which boxes are selected.  A node with no children is a
  degenerate node or leaf, and represents a box that is not subdivided.  In
  the tree representation of the subpaving, a leaf indicates that the box it
  represents belongs to the subpaving.\anchor minimal A tree (or the subpaving it represents)
  is \e minimal if it has no sibling leaves; that is, no node has more than one
  leaf-child.  [AIA2001, p. 52]. The presence to two sibling leaves (which
  means that the subpaving is non-minimal) indicates that a box is subdivided
  and both subboxes (child nodes) retained.

  \image html minimalsp1.png "A minimal 2-dimensional subpaving and its representation as a binary tree" width=15cm
  \image latex minimalsp1.png "A minimal 2-dimensional subpaving and its representation as a binary tree" width=15cm

<HR>

  \section setcomputation Regular subpavings and set computation
  Jaulin, Kieffer, Didrit and Walter use regular %subpavings for set
  computation.  In particular, the computation of a reciprocal image
  X = <b>f</b><SUP>-1</SUP>(Y) where Y is regular subpaving of R<SUP>m</SUP>
  (set inversion), and computation of the direct image of a subpaving
  Y  = <b>f</b>(X) where X is a regular subpaving of R<SUP>n</SUP>
  (image evaluation).  [AIA2001] outlines the creation of a SUBPAVINGS class
  and the implementation their algorithms using this class for set inversion
  and image evaluation in C++ given some library for interval analysis: the
  basic requirements of such an interval analysis library are given and their
  own implementation uses the PROFIL/BIAS library.  We have implemented the
  algorithms for set inversion and image evaluation using the C-XSC interval
  library with as little change to the structure and code used in [AIA2001] as
  possible, other than that necessitated by the use of C-XSC.

  See:
- \ref AIASubPavings

<HR>

  \section statssetprocessing Subpavings and statistical set processing
  Our interest in %subpavings is more general than set computation in the sense
  of \ref setcomputation.  Regular
  %subpavings represent a special case of the general kd tree structure: a
  space partitioning data structure for organising points in k-dimensional
  space.  We aim to combine interval analysis and inclusion functions with
  kd trees and manipulation of tree structures.  In particular, we are
  interested in the use of %subpavings for statistical set processing.

  We start with a re-implemention in C++ of the subpaving-as-a-binary-tree class
  based on [AIA2001].  Again, we use C-XSC as the interval analysis library.
  We alter the structure of the class and functions to give us a better basis
  for using the class as a base class for derived classes specifically for
  statistical set processing.  See \ref newsubpavings.

  We extend this base class to a class specifically designed for processing
  statistical sample data.  The subpaving becomes a way to organise and
  summarise the sample data.  The boxes of the subpaving can be thought of a
  'buckets' for the data points contained in its interval vector box.  The
  binary tree representation of the subpaving shows how the subpaving has
  been formed.  Each node in the tree can maintain summaries (such as count,
  and sum) of the data it represents.   See \ref StatsSubPavings.  The
  statistical subpaving class is used in our work on \ref AdaptiveHistograms
  by extending the algoritthms and data structures in \ref setcomputation .

  See:
- \ref newsubpavings
- \ref StatsSubPavings
- \ref AdaptiveHistograms

*/

/*! \page AIASubPavings Implementation of Applied Interval Analysis procedures under C-XSC

- \ref AIAsec_setcomputation
- \ref AIAsec_SIVIA
- \ref AIAsec_imageSp
- \ref AIAsec_examples

<HR>

\section AIAsec_setcomputation Regular subpavings and set computation

[AIA2001] show how %subpavings can approximate compact sets in a guaranteed way
  and how computation on %subpavings allows approximate computation on these
  compact sets.  For any full compact set X we can find a subpaving
  'lower bound' of X (which can be thought of an an 'inner paving' of X) and a
  subpaving 'upper bound' of X (or \anchor outerpaving 'outer paving' of X)
  as 'close' to X as we desire (distance and ordering of sets and hence what
  is meant by 'close' has to be defined of course). X is enclosed by the inner
  and outer paving. [AIA2001, pp. 48-49]

To get 'close' to the set we want to enclose, we can take a large box
  (the root node of the binary tree representing the subpaving) which we know
  encloses the set and then progressively subdivide the box and the subboxes,
  etc, checking to see if we can discard any subbox (prune that branch) to get
  a subpaving/tree representation of the subpaving which meets our criteria
  for 'closeness' to our target compact set.

<HR>

\section AIAsec_SIVIA Set inversion using interval analysis
Set inversion is the computation of a reciprocal image
  X = <b>f</b><SUP>-1</SUP>(Y) where Y is regular subpaving of
  R<SUP>m</SUP>.  [AIA2001] develop an algorithm SIVIA (Set Inverter Via
  Interval Analysis) for finding an outer subpaving of X to a specified level
  of precision.   Starting with a large search box [<b>x</b>](0) to which the
  outer subpaving of X is guaranteed to belong and given a inclusion function
  [<b>f</b>]([<b>x</b>]) for <b>f</b>, we form a subpaving by progressively
  bisecting and then testing the resulting subboxes to see if they should be
  included in our solution subpaving.

The test for any box [<b>x</b>] is a comparison of [<b>f</b>]([<b>x</b>]) to
  Y.  This test has three possible outcomes

<OL>
<LI> If [<b>f</b>]([<b>x</b>]) has an empty intersection with Y then [<b>x</b>]
  does not belong to X and can be cut off from the solution tree
<LI> If [<b>f</b>]([<b>x</b>]) is entirely in Y then [<b>x</b>] belongs to the
  solution and is kept as part of the solution tree
<LI> If [<b>f</b>]([<b>x</b>]) has a non-empty intersection with Y but is not
  entirely in Y then [<b>x</b>] is undetermined.  If [<b>x</b>] has width
  greater than the specified precision parameter it is bisected and each
  subbox (child in the tree structure) is then tested.  If [<b>x</b>] has
  width less than or equal to the specified precision parameter it included in
  the solution, ie the suppaving outer bound of the set X [AIA2001, pp. 55-56].
</OL>

Thus the solution, the outer subpaving of the set X, includes the
  'uncertainty' layer of boxes that cross the boundaries of the set but whose
  width is small enough to satisfy our precision criteria.  We can make the
  uncertainty layer thinner by having a smaller precision and subdividing and
  testing the boxes in that layer further [AIA2001, p. 56].

An implementation of SIVIA using C++ is described in [AIA2001, pp. 339-342].
  This includes an extension of the usual booleans TRUE and FALSE to a set of
  <EM>interval booleans</EM> which includes notion of \e indeterminate.

SIVIA can also be used to evaluate an outer subpaving of the image of set by a
  function provided that the function is invertible in the usual sense.  This
  involves specifying the inclusion function [<b>f</b><SUP>-1</SUP>]
  for <b>f</b><SUP>-1</SUP> the inverse of <b>f</b> and taking as an initial
  search box some box guaranteed to contain the image of the set.

<HR>

\section AIAsec_imageSp Image evaluation

Computation of the direct image of a subpaving Y  = <b>f</b>(X) where X is a
  regular subpaving of R<SUP>n</SUP> (image evaluation) is more complex in the
  case where the function is not invertible in the usual sense.  AIA2001
  develop an algorithm ImageSP for finding an outer subpaving of Y to a
  specified level of precision.

Given a subpaving X for which we seek an outer approximation of the image and
an inclusion function [<b>f</b>]([<b>x</b>]) for <b>f</b>, ImageSP proceeds as
follows:

<OL>
<LI> Mince X (recursively bisect each box) until each box in X has width less
  than the specified level of precision.  (X will no longer be a minimal
  subpaving.)
<LI> For each [<b>x</b>] in the minced X, add [<b>f</b>]([<b>x</b>]) to a list
  of image boxes U and compute the interval hull of the union of all these
  image boxes.
<LI> If [<b>f</b>]([<b>x</b>]) Merge these boxes into a single, minimal,
  subpaving the root of which corresponds to the hull and which only contains
  boxes with width lower than the specified level of precision.
</OL>

An implementation of ImageSP using C++ is described in [AIA2001, pp. 342-347].

We have reimplemented the algorithms for set inversion and image evaluation
  using the C-XSC interval library with as little change to the structure and
  code used in [Applied Interval Analysis, Springer, 2001] as possible, other
  than that necessitated by the use of C-XSC.  Our class is called AIASPnode
  (a pointer to an AIASPnode is aliased as AIASubPaving)

<HR>

\section AIAsec_examples Examples from Applied Interval Analysis using AIASPnodes

AIA2001 show how computation on %subpavings allows approximate computation on
compact sets.  Their SUBPAVINGS class and various supporting procedures are
used to implement their two main algorithms:

<UL>
<LI> SIVIA (Set Inversion Via Interval Analysis).  Set inversion is the
  computation of a reciprocal image X = <b>f</b><SUP>-1</SUP>(Y) where Y is
  regular subpaving of R<SUP>m</SUP>
<LI> ImageSp (Image evaluation). Computation of the direct image of a
  subpaving Y  = <b>f</b>(X) where X is a regular subpaving of R<SUP>n</SUP>
</UL>

The following examples show how we have replicated the SUBPAVING class and the
  implementations of SIVIA and ImageSp created by AIA2001 and available on
their website  http://www.lss.supelec.fr/books/intervals/)

- \ref AIAexamsec_11_33
- \ref AIAexamsec_11_35
- \ref AIAexamsec_3_3
- \ref AIAexamsec_3_4

Our implementation uses the C-XSC library but otherwise makes as little
  alteration to the structure and code used in Applied Interval Analysis as
  possible, other than that necessitated by the use of C-XSC.  Our class is
  called AIASPnode (a pointer to an AIASPnode is aliased as AIASubPaving).

The AIASPnode class declarations and inline definitions are in the header
  file AIAsubpaving.hpp. Other definitions are in the file AIAsubpaving.cpp

\subsection AIAexamsec_11_33 Exercise 11.33

See AIA2001, p. 342.

Exercise 11.33 explores the use of SIVIA and the AIASPnode class for set
  inversion.  In particular, this example demonstrates that it is possible to
  use to evaluate the direct image of a set by a function, <em>provided that
  the function is invertible in the usual sense</em>.

Our implementation of this example is in the file Exr_11_33.cpp, which has
  header file Exr_11_33.hpp

Exr_11_33.hpp shows the includes for this example

\dontinclude Exr_11_33.hpp

\skip include
\until fstream

The file AIAsubpaving.hpp is included in order to be able to use the AIASPnode
  class it declares.  <time.h> is included in order to use the clock method
  and get information on the time take to run the example.  <fstream> is
  included in order to be able to output data to a file.

Turning now to Exm_11_33.cpp,

\dontinclude Exr_11_33.cpp

First we include the header file

\skip include
\until include

Then we declare some AIASubPavings as global (AIASubPaving is declared as an
  alias for a pointer to an AIASPnode in the AIAsubpavings.hpp header file).
  Declaring them as global means that they are available to all the functions
  in the file (they have file scope) and so they do not need to be passed as
  parameters to the functions that use them.

\skip These
\until Sc2

\remarks Use of globals is normally discouraged in good programming,

Then we specify the interval boolean tests we are going to use in the
  example.  The interval boolean test is the key to set inversion with interval
  analysis.  Recall that set inversion is the computation of a reciprocal image
  X = <b>f</b><SUP>-1</SUP>(Y) where Y is regular subpaving of R<SUP>m</SUP>
  (see \ref AIAsec_SIVIA).  Y is the subpaving we want to invert.  The interval
  boolean test takes a box in 'x-space' and tests whether the image of the box
  in 'y-space' ie the image under the inclusion function [<b>f</b>], is in the
  subpaving Y.  The interval boolean test can return one of the special
  interval boolean types BI_TRUE (the image of the box is inside Y), BI_FALSE
  (the image of the box is outside Y), or BI_INDET
  (indeterminate: the image of the box [<b>f</b>]([<b>x</b>] is partly in and
  partly out of the subpaving Y, ie overlaps the boundary of Y).

The first test \anchor IBTAnnular IBTAnnular tests a 2-d box in IR<SUP>2</SUP>
   for inclusion in the set corresponding to the area between circles centred
   at the origin and with radii 1 and 2.  The test is performed on the interval
   image [<b>f</b>]([<b>x</b>]) of the box [<b>x</b>] supplied as the function
   argument.  The test itself specifies the inclusion function f [<b>f</b>]
   (x[1]<SUP>2</SUP> + x[2]<SUP>2</SUP>).  In this case the test also specifies
   the subpaving to invert - in this case this is a very simple subpaving being
   just a single interval [1 2].  The test can return BI_TRUE, BI_FALSE,
   or BI_INDET.

\skip specifications
\until }

The second test \anchor IBTFdirect IBTFdirect is also an interval boolean test
but on first sight looks rather different to IBTAnnular.  Taking it apart, we
can how it does the same thing.  The test takes a ivector argument and
calculates the interval image of this using an inclusion function.  The
inclusion function is again specified within the test.  In this case it is in
fact the inverse of an invertible function: this test will be used to run
SIVIA 'backwards'.  Note that the subpaving we are inverting is one of
the \anchor AIAglobal global AIASubPaving variables, Sc. Because the subpaving
is not a simple interval as in the previous test, the test uses the
AIASPnode::operator<=(const ivector&, AIASubPaving), which returns an
AIA_BOOL_INTERVAL type just like the test above.

\skip IBTFdirect(const
\until }

The third test \anchor IBTFinverse IBTFinverse is similar again.  The function
f is the same function as that for which we used the inverse in the
test \ref IBTFdirect "IBTFdirect" above.  The subpaving we test for inclusion
in is another of the globals, Sc1, and again we use the
AIASPnode::operator<=(const ivector&, AIASubPaving) to returns an
AIA_BOOL_INTERVAL.

\skip IBTFinverse(const
\until specification

Now we say that we are using the std and cxsc namespaces to avoid having to
 type cxsc::ivector or std::cout etc.

\skip cxsc
\until std

In the main process we start by declare some of the variables we will be using

\skip main
\until new

A 2-dimensional interval vector x is set up which is the 2-d box
[-5 5]<SUP>2</SUP>.  This is used to initialise a newed AIASPnode in dynamic
memory, with A as an AIASubPaving or pointer to an AIASPnode.  Creating the
node object in dynamic memory allows us to use pointers to it outside the scope
of the function which created the object

Now we use the SIVIA algorithm (see \ref AIAsec_SIVIA) to find a 'subpaving
characterisation' Sc containing the area between circles centred on the origin
with radii 1 and 2 in 2 dimensional space.  We say 'subpaving characterisation'
because the subpaving is an \ref outerpaving "outer subpaving" of the actual
area between circles centred on the origin with radii 1 and 2.

\skip SIVIA
\until prec

A description of what we are doing is sent to standard output, The program asks
the user to input a precision between 1 and 0.001.

The clock is started and Sivia(AIA_PIBT BoolTest, AIASubPaving A, double eps)
is called.  The AIA_PIBT (pointer to an interval boolean test)
is \ref IBTAnnular, as described above.  A is provided as an initial search
box, and prec is given as the argument for the eps parameter. The clock is
stopped when Sivia has returned a value for Sc.

\skip start
\until end

Sc now points to the root node of a tree representing the subpaving constructed
with Sivia.

The SIVIA algorithm works by taking an initial search box, testing it,
rejecting those returning BI_FALSE, including those returning BI_TRUE in the
subpaving it is building, and bisecting again if the interval boolean test
returns BI_INDET, and recursively sending each subbox to SIVIA to be tested
similarly.  Each subbox continues to be bisected until either the test returns
a clear BI_TRUE or BI_FALSE  or the test is BI_INDET but the width of the box
is the value given for eps, which means that it is thin enough to be included
in the subpaving to be returned.  The eps parameter provides a 'stopping rule'
prevents SIVIA recursing endlessly.  A larger value for eps will mean a thicker
uncertainty layer in the outer subpaving of the area we seek to characterise.

We stop the clock after calling Sivia and report the computing time and the
volume and number of leaf boxes of the subpaving.

\skip cout
\until Volume(Sc)

Then we create an output file of giving the subpaving as a list of interval
vectors.  The file name is specified in the example as AIAannular.txt and the
first three lines of the file give, respectively, the dimensions we are
working in, the initial search box, and the precision used.  The subpaving
itself is then output.  The format used for outputting the subpaving is
specified in operator<<(std::ostream &, AIASubPaving)

\skip realize
\until determined

This example now moves on to use this subpaving in further demonstrations of
the SIVIA algorithm.  First, we delete the AIASPnode that A currently points
to and replace it with a new one, again providing an initial search
box [-5.0 5.0]<SUP>2</SUP>

\skip make
\until new

Now we prepare to use SIVIA to find the direct image of an invertible function,
which can be thought of as using SIVIA in reverse: normally SIVIA characterises
the reciprocal image X = <b>f</b><SUP>-1</SUP>(Y) where Y, the subpaving we
want to invert, is regular subpaving of R<SUP>m</SUP>.  However, if we have X
but want to find Y (or a characterisation for Y) and f is invertible so that we
can specify an inclusion function [<b>f</b><SUP>-1</SUP>] for
<b>f</b><SUP>-1</SUP> then we can use this inclusion function to characterise
Y.  SIVIA can only be used to find a direct image under <b>f</b> when <b>f</b>
is invertible in the normal sense because we need to be able to specify
[<b>f</b><SUP>-1</SUP>] in our interval boolean test.

A description of what we are doing is printed to standard output and the user
again asked to enter a precision, and the clock is started

\skip testing
\until prec

Now we run Sivia as above, but the argument supplied for the interval boolean
test parameter is \ref IBTFdirect as given at the top of the file.  This
specifies [<b>f</b><SUP>-1</SUP>] for
f<SUB>1</SUB>(x) = 2x<SUB>1</SUB> - x<SUB>2</SUB>,
f<SUB>2</SUB>(x) = -x<SUB>1</SUB> + 2x<SUB>2</SUB>.
Recall that \ref IBTFdirect specifies that the subpaving to invert is Sc, ie
the subpaving we found using Sivia above.

\skip start
\until end

Sc1 now points to the root node of a tree representing the subpaving
constructed with Sivia.

The computing time, volume and number of leaves are reported and an output
file produced.

\skip cout
\until example

We have now found a subpaving Sc in 'x-space' and used Sivia on an function
inverse to get a characterisation Sc1 of the image of Sc in 'y-space'.  The
final step is to use Sivia again to go back to 'x-space' and find the
reciprocal image Sc2 of Sc1!

Start by resetting A again, describe what we are doing, and get the precision

\skip initial
\until prec

Now the argument supplied for the interval boolean test parameter
is \ref IBTFinverse from the top of the file.  This specifies [<b>f</b>]
for f<SUB>1</SUB>(x) = 2x<SUB>1</SUB> - x<SUB>2</SUB>,f<SUB>2</SUB>(x)
= -x<SUB>1</SUB> + 2x<SUB>2</SUB>.  \ref IBTFinverse also specifies that the
subpaving to invert is Sc1, ie the subpaving in 'y-space' we have just found
using Sivia above.

\skip start
\until end

Sc2 now points to the root node of a tree representing the subpaving
constructed with Sivia.

The computing time, volume and number of leaves are reported and an output
file produced.

\skip cout
\until precision

Finally we need to delete all the AIASPnode trees we created in dynamic
memory, then the program can return.

\skip delete
\until }

The example can be run with different values supplied for the precision in
each case, to examine the effect that this has on the volume and number of
leaves of the %subpavings.

One of the interesting points about this example is that we can see the effect
of pessimism in inclusion functions (AIA2001, pp. 15-17, p. 342).  We start by
making a subpaving Sc in 'x-space' and then find Sc1, a subpaving
characterisation for the image of Sc under a function <b>f</b> (ie, a
subpaving in 'y-space'), and then go back to 'x-space' again and with Sc2 as a
subpaving characterisation for the reciprocal image of Sc1 under <b>f</b>.

Under the output format currently specified, the the first part of the output
file AIAdirect.txt looks like this:

\verbatim
2
[ -5.000000,  5.000000] [ -5.000000,  5.000000]
Precision is 0.05
[  -1.953125 ,  -1.875000 ] , [  -0.078125 ,   0.000000 ]
[  -1.289062 ,  -1.250000 ] , [  -0.820312 ,  -0.781250 ]
[  -1.406250 ,  -1.328125 ] , [  -0.703125 ,  -0.625000 ]
[  -1.328125 ,  -1.289062 ] , [  -0.742188 ,  -0.703125 ]
\endverbatim

All the output from this example rendered graphically looks like this.

\image html AIAexample11_33.png "Results for Exercise 11.33 using precision 0.05" width=15cm
\image latex AIAexample11_33.png "Results for Exercise 11.33 using precision 0.05" width=15cm

\subsection AIAexamsec_11_35 Exercise 11.35

See examples/AIA/Exr_11_35

\subsection AIAexamsec_3_3 Example 3.3

See AIA2001, p. 61.

Example 3.3 explores the use of ImageSp and the AIASPnode class for image
evaluation.  In \ref AIAexamsec_3_3 used Sivia to find the image of a
subpaving under a function <b>f</b>, but this is only possible when <b>f</b>
is invertible in the usual sense.

Computation of the direct image of a subpaving Y  = <b>f</b>(X) where X is a
regular subpaving of R<SUP>n</SUP> (image evaluation) is more complex in the
case where the function is not invertible in the usual sense.

Our implementation of this example is in the file Exm_3_3.cpp, which has header
file Exm_3_3.hpp.

Exm_3_3.hpp shows is similar to the header for Example 11.33 so we will move
straight on the describing the file Exm_3_3.cpp.

\dontinclude Exm_3_3.cpp

First we include the header file and declare more global %subpavings.

\skip include
\until Sc4

Now we look at the functions used in our main program.

The first is a boolean interval test \ref IBTAnnular is exactly the same as
that used for Exercise 11.33.

\skip specifications
\until end

The next function is an interval vector function, \anchor IVFex3_3 IVFex3_3.
An interval vector function returns the interval vector image of an interval
vector x under <b>f</b> for <b>f</b> as specified in the function and x
supplied as the function argument.

In this case <b>f</b>: R<SUP>2</SUP> -> R<SUP>2</SUP>,
<b>f</b><SUB>1</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) =
x<SUB>1</SUB>x<SUB>2</SUB>, <b>f</b><SUB>2</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>)
= x<SUB>1</SUB> + x<SUB>2</SUB>

\skip specification
\until end

Then we say what namespaces we are using, start the main program, and set up
an initial search box [-3.0 3.0]<SUP>2</SUP> and a pointer A to an AIASPnode
based on this search box.

\skip using
\until new

The first subpaving we create is again Sc, the subpaving covering the set in
R<SUP>2</SUP> between circles centred at the origin and with radii 1 and 2.

\skip SIVIA
\until determined

Now we move on to using ImageSp.  We are going to find Sc4, a characterisation
of the set f(Sc) using f as defined in \ref IVFex3_3 "IVF_ex3_3".

First we say what we are doing and ask for a precision.

\skip cout
\until prec

Then start the clock and run ImageSp(AIA_PIVF, AIASubPaving A, double eps)
giving as arguments our interval vector function, our subpaving node A, and
our precision.

Sc4 now points to the root node of a tree representing the subpaving
constructed with ImageSp.

The ImageSp algorithm works by taking an initial subpaving (in this case, the
box of A), mincing it up into a fine subpaving where every box has width less
than precision, and finding the image of each of these boxes with the
specified interval vector test.  It then forms a minimal regular subpaving
which covers the union of all these image boxes, again with precision as
specified.

\skip start
\until end

We report the running time, volume and number of leaves in the image subpaving
and send the output to a txt file

\skip cout
\until written

Finally, we delete the %subpavings and end the program.

\skip delete
\until }

The %subpavings produced by this program run using precision 0.05 for both
%subpavings is shown graphically below.  As well as showing the initial
subpaving (represented by Sc) and the image subpaving (represented by Sc4),
we capture an intermediate step within ImageSp.  This is the evaluation step,
where we have a large set of (possibly overlapping) image boxes formed from
all the minced up subboxes of the initial box ([-3.0 3.0]<SUP>2</SUP> chopped
up so that each one is less than 0.05 wide).  This set can be compared to the
final regular minimal subpaving characterisation of the image which is formed
by the function Regularize.

\image html AIAexample3_3.png "Results for Example 3.3 using precision 0.05" width=15cm
\image latex AIAexample3_3.png "Results for Example 3.3 using precision 0.05" width=15cm

\subsection AIAexamsec_3_4 Example 3.4

See AIA2001, pp. 61-63.

\anchor AIAexample3_4 Finally we show another example using both Sivia and
ImageSp.  This example shows why we need ImageSp when to evaluate the images
of functions that are only invertible in a set-theoretic sense, and also
illustrates the interesting effects that can occur in these circumstances.

Our implementation of this example is in the file Exm_3_4.cpp, which has
header file Exm_3_4.hpp.

Exm_3_4.hpp shows is similar to the other header files so we will move
straight on the describing the file Exm_3_4.cpp.

\dontinclude Exm_3_4.cpp

First we include the header file and declare more global %subpavings.

\skip include
\until Sc7

Now we look at the functions used in our main program.

The first is a boolean interval test, \anchor IBT_ex3_4 IBT_ex3_4.  The first
specifies a function <b>f</b>: R<SUP>2</SUP> -> R such that
<b>f</b>(x<SUB>1</SUB>, x<SUB>2</SUB>) = x<SUB>1</SUB><SUP>4</SUP> -
x<SUB>1</SUB><SUP>2</SUP> + 4x<SUB>2</SUB><SUP>2</SUP> and tests the image of
a box x under this function for inclusion in the subpaving to be inverted.  As
with \ref IBTAnnular "IBTAnnular" above, the subpaving to be inverted is also
specifed in the test and here is just the interval [-0.1, 0.1].

\skip specifications
\until }

Then we have another boolean interval test,
\anchor IBTinverse_ex3_4 IBTFinverse_ex3_4.  Again this specifies a
function <b>f</b>, this time <b>f</b>: R<SUP>2</SUP> -> R<SUP>2</SUP>,
f<SUB>1</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) = (x<SUB>1</SUB> - 1)<SUP>2</SUP>
-1 + x<SUB>2</SUB>, f<SUB>2</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) =
-x<SUB>1</SUB><SUP>2</SUP> + (x<SUB>2</SUB> - 1)<SUP>2</SUP>.
The subpaving to invert is one of the globally defined subpaving Sc6

\skip IBTFinverse_ex3_4(const
\until end

The next function is an interval vector function, \anchor IVF_ex3_4 IVF_ex3_4.
Recall that an interval vector function returns the interval vector image of an
interval vector x under <b>f</b> for <b>f</b> as specified in the function and
x supplied as the function argument.  In this case <b>f</b> is exactly the same
as in the boolean interval test \ref IBTinverse_ex3_4 "IBTFinverse_ex3_4"
above.

\skip specification
\until end

Then we say what namespaces we are using, start the main program, and set up an
initial search box [-3.0 3.0]<SUP>2</SUP> and a pointer A to an AIASPnode based
on this search box.

\skip using
\until new

First we use SIVIA to get a subpaving Sc5 in 'x-space' using the function f and
subpaving to invert specified in \ref IBT_ex3_4 "IBT_ex3_4", ie Sc5 represents
a subpaving characterisation of the set (x<SUB>1</SUB>, x<SUB>2</SUB>) such
that  x<SUB>1</SUB><SUP>4</SUP> - x<SUB>1</SUB><SUP>2</SUP> +
4x<SUB>2</SUB><SUP>2</SUP> is in the interval [-0.1, 0.1].

\skip SIVIA
\until Number

And send the output to a txt file.

\skip realize
\until written

Now we want to find the image of this set using the function specified in
\ref IVF_ex3_4 "IVF_ex3_4".  <b>f</b>: R<SUP>2</SUP> -> R<SUP>2</SUP>,
f<SUB>1</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) = (x<SUB>1</SUB> - 1)<SUP>2</SUP>
-1 + x<SUB>2</SUB>, f<SUB>2</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) =
-x<SUB>1</SUB><SUP>2</SUP> + (x<SUB>2</SUB> - 1)<SUP>2</SUP>.

\skip ImageSp
\until prec

This function is certainly not invertible in the usual sense and we have to
use ImageSp to find the image. We supply the interval vector function
IVF_ex3_4, the subpaving Sc5 and the precision input by the user as arguments
in the call to ImageSp.

\skip start
\until end

Sc6 now points to a subpaving representation in 'y-space', the image of Sc5
under the function <b>f</b> in \ref IVF_ex3_4 "IVF_ex3_4".

We report computing time, volume and number of leaves and send the suppaving
output to a txt file

\skip cout
\until written

Subpaving A, which will have been minced and mangled in the process of forming
Sc6, is now deleted and remade to provide another initial search box.

\skip remake
\until new

What happens if we now use Sivia on our ImageSp-created image Sc6 in 'y-space',
inverting our image to get back to 'x-space'?

\skip SIVIA
\until prec

The boolean interval test \ref IBTinverse_ex3_4 "IBTFinverse_ex3_4" specifies
the function and also specifies Sc6 as the subpaving to test for inclusion in.
We give the new search box in the subpaving A, and the user-supplied precision.

\skip start
\until end

Sc7 is now points to 'x-space' subpaving characterisation of the reciprocal
image of Sc6, which was in turn a subpaving characterisation of Sc5.

We report computing time, volume and number of leaves and output the subpaving
to a txt file.

\skip cout
\until written

and then delete our %subpavings and end the program

\skip delete
\until }

The %subpavings produced by this program run using precision 0.05 in each case
is shown graphically below.

As well as showing the initial subpaving (represented by Sc5), the image
subpaving (represented by Sc6), and the reciprocal image of the image (Sc7),
we capture an intermediate step in the process of creating Sc6 from Sc5 with
ImageSp.  This is the evaluation step, where we have a large set of (possibly
overlapping) image boxes formed from all the minced up subboxes of the initial
box ([-3.0 3.0]<SUP>2</SUP> chopped up so that each one is less than 0.05
wide).  This set can be compared to the Sc6, the regular minimal subpaving
characterisation of the image.

The most interesting comparison though is between the initial subpaving (Sc5)
and the subpaving for a reciprocal image of the image subpaving (Sc7).  The
initial set (characterised by Sc5) is in there but the result is fatter due to
error accumulation in the process of going to 'y-space' and back, and the final
subpaving has additional parts which appear because <b>f</b>(.) is only
invertible in the set-theoretic sense (AIA2001, pp. 63)

\image html AIAexample3_4.png "Results for Example 3.4 using precision 0.05" width=15cm
\image latex AIAexample3_4.png "Results for Example 3.4 using precision 0.05" width=15cm

*/

/*! \page newsubpavings Subpavings as a basis for statistical set processing

- \ref newsec_subpavingsKd
- \ref newsec_reimplementation
- \ref newsec_minimalsp
- \ref newsec_statisticalsp

<HR>

\section newsec_subpavingsKd Subpavings and kd-trees
Regular %subpavings represent a special case of the general kd tree structure:
  a space partitioning data structure for organising points in k-dimensional
  space.  We aim to combine interval analysis and inclusion functions with kd
  trees and manipulation of tree structures in order to use %subpavings for
  statistical set processing.

<HR>

\section newsec_reimplementation A revised class for subpavings
We start with a re-implemention in C++ of a subpaving-as-a-binary-tree object
based on [AIA2001].  Our new class is based on a class called
\link subpavings::SPnode SPnode\endlink, and a pointer to an SPnode is
aliased as \link subpavings::SubPaving SubPaving\endlink.  Again, we use
C-XSC as the interval analysis library.  This base class and functions 
are designed to also give us a basis for derived classes specifically 
designed for various forms of statistical set
processing.  Our new tree has links from child to parent as well as from parent
to child (ie, each child knows who its parent is); this may be used in 'shrink'
the tree, absorbing children back up into the parent.

Other new members of the SPnode class include \link subpavings::SPnode::nodeName nodeName\endlink.
The nodeName is a used as a convenient way to assist human interpretation of
a tree of nodes.  A \anchor SPnodename nodeName is based on the parent node's name suffixed with
'L' if the node is the left child of that parent and by 'R' if the node
is the right child of that parent.  A root node is usually named 'X'.

The SPnode class declarations and inline definitions are in the header file
spnode.hpp. Other definitions are in the file spnode.cpp.  Type definitions
relevant to SPnodes are in the file sptypes.hpp.  General tool functions
useful for spnodes are in sptools.hpp and sptools.cpp. Algorithms used for set
computation are declared in spalgorithms.hpp and defined in spalgorithms.hpp.
Template functions requiring node concepts are in sptemplates.hpp.


<HR>

\section newsec_minimalsp Extending the SPnode class for set computation

A derived class based on the SPnode class can be used
to reimplement the set computations described
in AIA2001.

	See: 
- \subpage minimalsubpavings

<HR>

\section newsec_statisticalsp Extending the SPnode class for statistical set processing

  We create new classes derived from the 
  \link subpavings::SPnode SPnode\endlink  base class, each class a form of
  SPnode specialised for a particular kind of statistical set processing.

  See:
- \ref StatsSubPavings

*/

/*! \page minimalsubpavings Minimal subpavings for set computation

The derived class (based on the SPnode class) used
to reimplement the set computations described
in AIA2001 is \link subpavings::SPMinimalnode SPMinimalnode\endlink.  
A pointer to an SPMinimalnode is
aliased as \link subpavings::SubPaving MinimalSubPaving\endlink.  This
class is specifically for minimal trees of nodes (representing minimal
subpavings).  As discussed \ref minimal "above", a tree (or the subpaving it represents)
is \e minimal if it has no sibling leaves; that is, no node has more than one
leaf-child. 

<HR>

\section minsec_equivalence Equivalence to the AIASPnode class
In order to demonstrate that this implementation is equivalent to that of
AIA2001, we reimplement the set inversion and image evaluation functions and
using our new SPMinimalnode class.  

For examples using the new SPMinimalnode class for set computation and demonstrating
the equivalence of the results to those in \ref AIAsec_examples
see \ref minexamples


\subsection minexamples Example of set computation with SPMinimalnodes

SPMinimalnodes are re-implemention in C++ of a subpaving-as-a-binary-tree object based
on [AIA2001].  SPMinimalnodes behave exactly as AIASPnodes
do but we highlight here some of the minor changes that we have made the way
that this class is coded and used.



\subsubsection minexamsec_3_4 Example 3.4

In this example show how example 3.4 from AIA2001, pp. 61-63, can be performed
using the SPMinimalnode class (for the same example using AIASPnodes, code
closely based on AIA2001, see \ref AIAexample3_4 "Example 3.4 using AIAPSnodes").
This example uses our implementations of both Sivia
(see \ref AIAsec_SIVIA) for set inversion and ImageSp (see \ref AIAsec_imageSp)
for image evaluation.  The functions we use here are
Sivia (PIBT BoolTest, const SPMinimalnode * const toInvert,
                SPMinimalnode * const search, const double eps)
and
ImageSp(PIVF f, SPMinimalnode* spn, double eps).

Our implementation of this example is in the file NewExm_3_4.cpp, which has
header file NewExm_3_4.hpp. The header file is

\dontinclude NewExm_3_4.hpp

\skip #ifndef
\until #endif

We include headers and libraries we want to be able to use:  <time.h> for a
clock, <fstream> for file output, and suppavings.hpp for our class
declarations.  We also specify the namespaces we want to be able to use
without qualification.  Note that subpavings is now a namespace.  If we did not
specify that we are using this namespace we would have to qualify all
references to entities declared within it, for example
subpavings::MinimalSubPaving (MinimalSubPaving is an alias for a 
pointer to an SPMinimalnode).

In NewExm_3_4.cpp we include the header file

\dontinclude NewExm_3_4.cpp

\skip implementation
\until include

Then we define some functions used in our main program.

The first is a boolean interval test, \anchor NewIBT_ex3_4 IBT_ex3_4.

The interval boolean test is the key to set inversion with interval analysis.
Recall that set inversion is the computation of a reciprocal image
X = <b>f</b><SUP>-1</SUP>(Y) where Y is regular subpaving of R<SUP>m</SUP>
(see \ref AIAsec_SIVIA).  Y is the subpaving we want to invert.  The interval
boolean test takes a box in 'x-space' and tests whether the image of the box
in 'y-space' ie the image under the inclusion function [<b>f</b>], is in the
subpaving Y.  The interval boolean test can return one of the special interval
boolean types BI_TRUE (the image of the box is inside Y), BI_FALSE (the image
of the box is outside Y), or BI_INDET (indeterminate: the image of the box
[<b>f</b>]([<b>x</b>] is partly in and partly out of the subpaving Y, ie
overlaps the boundary of Y).

This test \anchor NewIBT_ex3_4 IBT_ex3_4 tests the image of a box x for
inclusion in a subpaving represented by the SPMinimalnode pointer spn.  Values for x
and spn are given by the caller as the function arguments.  <b>f</b> is
specified in the test as <b>f</b>: R<SUP>2</SUP> -> R such that
<b>f</b>(x<SUB>1</SUB>, x<SUB>2</SUB>) = x<SUB>1</SUB><SUP>4</SUP> - x<SUB>1</SUB><SUP>2</SUP> + 4x<SUB>2</SUB><SUP>2</SUP>.

\skip specifications
\until }

\note Interval boolean tests in the subpavings namespace take a pointer to
subpaving we want to invert as a parameter.  This avoids the use of global
subpaving variables in these tests as in
\ref AIAglobal "example 3.4 with AIASPnodes" and makes 
the functions somewhat clearer to interpret.

Then we have another boolean interval test, \anchor NewIBTFinverse_ex3_4
IBTFinverse_ex3_4.  Again this specifies a function <b>f</b>, this time

<b>f</b>: R<SUP>2</SUP> -> R<SUP>2</SUP>,
f<SUB>1</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) =
(x<SUB>1</SUB> - 1)<SUP>2</SUP> -1 + x<SUB>2</SUB>,
f<SUB>2</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) = -x<SUB>1</SUB><SUP>2</SUP> +
(x<SUB>2</SUB> - 1)<SUP>2</SUP>.
Values for x, the box whose image we test, and spn, a pointer to the subpaving
to invert, are given by the caller as the function arguments.

\skip SIVIA
\until end

The next function is an interval vector function,
\anchor NewIVF_ex3_4 IVF_ex3_4.  An interval vector function returns the
interval vector image of an interval vector x under <b>f</b> for <b>f</b> as
specified in the function and x supplied as the function argument.  In this
case <b>f</b> is exactly the same as in the boolean interval
test \ref NewIBTFinverse_ex3_4 "IBTFinverse_ex3_4" above.

\skip specification
\until end

Then we start the main program, declare some variables, and set up an initial
search box [-3.0 3.0]<SUP>2</SUP>. This is used to initialise a newed
SPMinimalnode in dynamic memory, with A as a MinimalSubPaving or pointer to
a SPMinimalnode.  Creating the node object in dynamic memory allows us to use
pointers to it outside the scope of the function which created the object

\skip main
\until new

We then set up the subpaving we want to invert (again, in dynamic memory).

\skip ivector
\until MinimalSubPaving

And then prepare to use SIVIA to do the set inversion.  We make some standard
output to say what we are doing and ask the user to enter a precision.

\skip SIVIA
\until prec

Now we use the SIVIA algorithm (see \ref AIAsec_SIVIA) to find a
'subpaving characterisation' Sc5 of the set we want with
Sivia(PIBT BoolTest, const SPMinimalnode * const toInvert, SPMinimalnode * const search,
      const double eps).
A PIBT is a pointer to a boolean interval test, in this
case \ref NewIBT_ex3_4 "IBT_ex3_4".  The set we want is defined by the
combination of the SubPaving ToInvert and the <b>f</b> specifed in the
boolean interval test \ref NewIBT_ex3_4 "IBT_ex3_4".  We are finding a
characterisation for the set (x<SUB>1</SUB>, x<SUB>2</SUB>) such that
<b>f</b>(x<SUB>1</SUB>, x<SUB>2</SUB>) = x<SUB>1</SUB><SUP>4</SUP> - x<SUB>1</SUB><SUP>2</SUP> + 4x<SUB>2</SUB><SUP>2</SUP> is in ToInvert,
which is the interval [-0.1 0.1].

\skip start
\until end

The SIVIA algorithm works by taking an initial search box, testing it,
rejecting those returning BI_FALSE, including those returning BI_TRUE in the
subpaving it is building, and bisecting again if the interval boolean test
returns BI_INDET, and recursively sending each subbox to SIVIA to be tested
similarly.  Each subbox continues to be bisected until either the test returns
a clear BI_TRUE or BI_FALSE  or the test is BI_INDET but the width of the box
is the value given for eps, which means that it is thin enough to be included
in the subpaving to be returned.  The eps parameter (precision) provides
a 'stopping rule' prevents SIVIA recursing endlessly.  A larger value for eps
will mean a thicker uncertainty layer in the outer subpaving of the area we
seek to characterise.

\note This version of Sivia takes a parameter for the SubPaving to invert
because it must pass this on to the boolean interval test.

Sc5 now points to the subpaving characterisation we have created.  We
say 'subpaving characterisation' because the subpaving is
an \ref outerpaving "outer subpaving" of the actual set we want.

We print the computing time, volume, and number of leaves of the subpaving we
create to standard output.

\skip time
\until number

And send the output to a txt file.

\skip realize
\until written

Now we want to find the image of this set using the function specified
in \ref NewIVF_ex3_4 "IVF_ex3_4".  <b>f</b>: R<SUP>2</SUP> -> R<SUP>2</SUP>,
f<SUB>1</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) = (x<SUB>1</SUB> - 1)<SUP>2</SUP>
-1 + x<SUB>2</SUB>, f<SUB>2</SUB>(x<SUB>1</SUB>, x<SUB>2</SUB>) =
-x<SUB>1</SUB><SUP>2</SUP> + (x<SUB>2</SUB> - 1)<SUP>2</SUP>.

\skip ImageSp
\until prec

This function is not invertible in the usual sense and we have to use ImageSp
to find the image. The version of ImageSp used is
ImageSp(PIVF f, SPMinimalnode* spn, double eps).  We supply the interval vector
function \ref NewIVF_ex3_4 "IVF_ex3_4", the subpaving Sc5 and the precision input by the user as
arguments in the call to ImageSp.

\skip start
\until end

Sc6 now points to the root node of a tree representing the subpaving
constructed with ImageSp.

The ImageSp algorithm works by taking an initial subpaving (represented by Sc5
in this case), mincing it up into a fine subpaving where every box has width
less than the precision specified, and finding the image of each of these boxes
with the specified interval vector test.  It then forms a minimal regular
subpaving which covers the union of all these image boxes, again with precision
as specified.

Sc6 now points to a subpaving representation in 'y-space', the image of Sc5
under the function <b>f</b> in \ref NewIVF_ex3_4 "IVF_ex3_4".

We report computing time, volume and number of leaves and send the suppaving
output to a txt file.

\skip cout
\until written

Subpaving A, which will have been minced and mangled in the process of
forming Sc6, is now deleted and remade to provide another initial search box.

\skip remake
\until new

What happens if we now use Sivia on our ImageSp-created image Sc6 in 'y-space',
inverting our image to get back to 'x-space'?

\skip Sc6
\until prec

The boolean interval test \ref NewIBTFinverse_ex3_4 "IBTFinverse_ex3_4"
specifies the function. We pass a pointer to this function to Sivia by
supplying it as the value for the PIBT parameter.  We give Sc6 as representing
the MinimalSubPaving to invert and A as reprsenting the initial search box, and we
also give the precision supplied by user.

\skip start
\until end

Sc7 is now points to 'x-space' subpaving characterisation of the reciprocal
image of Sc6, which was in turn a subpaving characterisation of Sc5.

We report computing time, volume and number of leaves and output the subpaving
to a txt file.

\skip cout
\until written

We delete our %subpavings created in dynamic memory before we end the program

\skip delete
\until }

The %subpavings produced by this program run using precision 0.05 in each case
are shown below.

As well as showing the initial subpaving, the image subpaving, and the
reciprocal image of the image, we capture an intermediate step in the process
of creating Sc6 from Sc5 with ImageSp.  This is the evaluation step, where we
have a large set of (possibly overlapping) image boxes formed from all the
minced up subboxes of the initial box [-3.0 3.0]<SUP>2</SUP> chopped up so
that each one is less than 0.05 wide.  This set can be compared to the Sc6,
the regular minimal subpaving characterisation of the image.

The most interesting comparison though is between the initial subpaving and
the subpaving for a reciprocal image of the image subpaving.  The initial set
(characterised by Sc5) is in there but the result is fatter due to error
accumulation in the process of going to 'y-space' and back, and the final
subpaving has additional parts which appear because <b>f</b>(.) is only
invertible in the set-theoretic sense (AIA2001, pp. 63)

\image html NEWexample3_4Coarse.png "Results for Example 3.4 using precision 0.05" width=15cm
\image latex NEWexample3_4Coarse.png "Results for Example 3.4 using precision 0.05" width=15cm

We can rerun the program to explore the effect of the precision.  The image
below was created using a precision of 0.001 for the initial subpaving and
then precisions of 0.01 to create the image of this in 'y-space' and and 0.01
to go back again to the reciprocal image of the image in 'x-space'.  The
evaluated images of the initial subpaving after it has been minced up very
finely (precision 0.001) are also shown.  The reciprocal image of the image is
a much better comparison to the initial subpaving as the errors accumulated are
smaller, but the additional parts due to are again apparent - they are caused
by the process, not by the precision to which we carry it out.

\image html NEWexample3_4Fine.png "Results for Example 3.4 using precision 0.001 for initial paving, then 0.01" width=15cm
\image latex NEWexample3_4Fine.png "Results for Example 3.4 using precision 0.001 for initial paving, then 0.01" width=15cm



*/


/*! \page StatsSubPavings Subpavings for processing statistical sample data

- \ref statssec_intro
- \ref statssec_organising
- \ref statssec_statsclass
- \ref statssec_usingStatsClass
- \ref statssec_technicalities

<HR>

\section statssec_intro Introduction

  %Subpavings for processing statistical sample data are a refinement of regular
  %subpavings where each box in a subpaving can be associated with data points in
  a sample.

<HR>

  \section statssec_organising Organising data
  The subpaving becomes a way to organise the sample data.  The boxes of the
  subpaving can be thought of a 'buckets' for the data points, each box
  containing some of the data points.  A data point x belongs in a box or
  interval vector [<b>x</b>] if [<b>x</b>] contains x.

  The binary tree representation of the subpaving shows how the subpaving has
  been formed.

  \image html nonminimalsp1.png "A 2-dimensional subpaving and its representation as a binary tree" width=15cm
  \image latex nonminimalsp1.png "A 2-dimensional subpaving and its representation as a binary tree" width=15cm

  Note that although we are still using \e regular %subpavings (ie, we are
  bisecting boxes along their longest dimension), our %subpavings are no
  longer \e minimal - that is, if a box is split, both sub-boxes are
  always retained in the paving.

  Using the binary tree representation of a subpaving we can think of the
  process of organising data as adding data starting at the root and having the
  data make its way through the tree to the leaf node whose box contains the
  data, going through the following steps at each node it gets to from the root
down:

  <OL>
  <LI> If the box of the node contains the data, continue <em>[if the box does
    not contain the data, return something indicating "not this node"]</em>
  <LI> If the box has been subdivided (node has children), see which subbox the
  data is in, ie recurse the addition process from step 1. on the children
  <LI> If the box has not been subdivided then the data can be associated with
  this node
  </OL>

  Thus the data moves through the tree until it finds a leaf node whose box
  contains it.

  Because a subpaving comprises  non-overlapping boxes, a data point should
  belong to at most only one box.  If the data point is outside the  box of the
  root node of a subpaving, it will clearly not be in any of the boxes of the
  leaves.

<HR>

  \section statssec_statsclass A class for statistical subpavings

We extend the \link subpavings::SPnode SPnode\endlink base class to a class
  specifically designed for processing statistical sample data.  This class is
called \link subpavings::SPSnode SPSnode\endlink and a pointer to an SPSnode
is aliased as \link subpavings::StatsSubPaving StatsSubPaving\endlink.

SPSnodes summarise as well as organise data:  each node can maintain summaries
  of the data it is associated with, such as a count (number of data points in
the box), sum, etc.  Because of the tree structure we can
  maintain a summary of the data at each level of the tree, not just in relation
  to the final subpaving (the leaf nodes).

  As the data moves down the tree, each node it passes through updates its
  summary statistics, such as a count and sum of data points covered by the box
  of that node.

  \image html statsspaddingdata.png "Summarising data with a binary tree" width=15cm
  \image latex statsspaddingdata.png "Summarising data with a binary tree" width=15cm

\link subpavings::SPSnode SPSnodes\endlink form non-minimal trees: an
SPSnode will either exactly two child nodes or no child nodes.  If a tree
of SPSnodes is representing a statistical
summary of data and the data is imagined plotted as coordinate points
on the root box of the paving, then an SPSnode representing a part of 
the subpaving where no data has fallen will have a
\link subpavings::SPSnode::countsOnly SPSnode::counter\endlink of zero.

By default, all recursively computable statistics provided are maintained
in each SPSnode.  However, since this uses memory and is not always needed,
an SPSnode can be constructed to \anchor SPScountsonly only maintain count statistics.  
See \link subpavings::SPSnode::countsOnly SPSnode::countsOnly\endlink.


<HR>

  \section statssec_usingStatsClass Using the statistical subpavings class
The statistical %subpavings class \link subpavings::SPSnode SPSnode\endlink is
  designed to be a flexible tool.

An SPSnode object represents only part of a subpaving.  The entire statistical subpaving
is represented by a tree of SPSnodes.  The behaviour of the tree of nodes will usually be
controlled by a wrapper or manager object.  In particular, there are no rules
  specified within the SPSnode class itself to determine when a box in the subpaving represented
by a tree of nodes is split. The wrapper or manager class will specify the methods which
determine the development of the subpaving/its tree representation.

An example of such a wrapper or manager object is the 
\link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink class, see \ref AdaptiveHistograms
  and \ref adhsec_examples

These manager classes will usually determine the rules for constructing
the subpaving and simply pass on 'instructions' to the node in the tree representing
a particular box in the subpaving.  This allows decisions to be made on the basis of the
characteristics of the statistical subpaving as a whole and not just an individual node.




<HR>

  \section statssec_technicalities Implementing the statistical subpavings class
For a full description of this class, see the 
\link subpavings::SPSnode SPSnode\endlink class documentation.  
Here we discuss only some of the important
  points of this class.

- \ref statssubsec_splitting
- \ref statssubsec_data
- \ref statssubsec_accumulating

  \subsection statssubsec_splitting Bisecting a box of a statistical subpaving
  When we bisect a box of a statistical subpaving to form two subboxes, we have
  to make a decision about which of the new boxes is considered to include the
  boundary hyperplane. This is so that, if we have a data point which falls
  exactly on that boundary, it is clear which of the subboxes it is associated
  with.

  We consider that the boundary to belong to the \b right child.  That is, the
  interval element of the \b right child's box on the splitting dimension is a
  left-closed interval whereas the interval element of the left child's box on
  the splitting dimension is a right-open interval.  For example, if the parent
  node box interval on the splitting dimension is [1 5] then the split takes place at
  the midpoint (3) and the right child is considered to have a box whose
  interval, on that dimesion, is [3 5] whereas the left child has [1 3).  [1 5]
is split into [1 3)[3 5].  A datapoint which is contained by the parent node's box and
which has value 3 on the splitting dimension will be considered to belong to
the right child's box and not to the left child's box.

\subsection statssubsec_data Associating data with a node
A node of a tree representation of a subpaving is associated with a particular
box of the subpaving.  The box may have been subsequently subdivided, in which
case the node will have children.  A leaf node represents an undivided box in
the subpaving.  As well as keeping summaries of the data associated with a
node (the data points covered by its box) we want to keep the full dataset
itself as well, but without duplication.  Thus, we only associate data with
leaf nodes.  When a box in the paving is subdivided the data associated with
the corresponding node is passed down to the new children.

In order to keep our implementation as flexible as possible, we maintain a
container for the sample data outside the SPSnode class.  A data point in
this 'big collection' of data is associated with a node through the iterator
to that data point in the big collection.  The node does not 'have' the data,
it simply knows where the data it is associated with is stored.  In the
implementation of SPSnode, a node has a data member which is a container for
iterators into the big data collection to store references to the data that it
is associated with.

\subsection statssubsec_accumulating C-XSC and processing statistical sample data

The \link subpavings::SPSnode SPSnode class\endlink is implemented using the
C-XSC library.  A data point is a cxsc::rvector and a box or interval vector
is a cxsc::ivector.  For example, a point in 3-d which would be represented in
cartesian coordinates as (1.0, 2.0, 3.0) is a 3-dimensional rvector whose
elements are, in order, the reals 1.0, 2.0, and 3.0.  This point would be
contained in (in cxsc, the operator <= meaning "is an element of") a 3
dimensional ivector.  The ivector whose elements are [0.0 2.0], [1.0 3.0],
[2.0, 3.0] defines a box which contains this point.

Summaries of data values are stored in csxc::dotprecision accumulators to
mitigate the risk of inaccuracy in summing floating point numbers (this risk
is higher when the numbers being summed have large differences in order of
magnitude, as will be the case when a relatively small new data point value
is added to a much larger accumulated sum of data points so far).
*/

/*! \page AdaptiveHistograms Adaptive histograms

- \ref adhsec_intro
- \ref adhsec_adhclass
- \ref adhsec_tree
- \ref adhsec_inputdata
- \ref adhsec_output
- \ref adhsec_averaging
- \ref adhsec_examples
- \ref adhsec_furtherdev

<HR>

\section adhsec_intro Introduction
An adaptive histogram is a histogram where bin widths and bin centres of
the partition of the data space adapt in some way to the data to be represented in
the histogram.  The basic ideas here are from classical statistical principles.

<HR>

\section adhsec_adhclass The AdaptiveHistogram class
The \link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink class organises
statistical sample data for the purpose
of creating adaptive histograms.  The class holds the container of sample data
and uses a tree of \link subpavings::SPSnode SPSnodes\endlink
(\link subpavings::StatsSubPaving StatsSubPaving\endlink) to represent
a regular non-minimal subpaving or partition of the dataspace
which can develop adaptively according to the data.

The AdaptiveHistogram object manages the tree through the tree's root node.  This
represents the box being \e covered by the subpaving, which is the data sample
space to be partitioned into histogram bins.  The boxes of the statistical
subpaving can be considered as histogram bins and are represented
by the leaf nodes of the tree.  Data is associated with boxes in the subpaving:
the count of data associated with a leaf node of the tree is the number of data
points falling into the box (bin) represented by that node.

\anchor ADHholdallstats The AdaptiveHistogram data member
\link subpavings::AdaptiveHistogram::holdAllStats holdAllStats\endlink controls
whether the tree of SPSnodes managed by the histogram object will maintain
\ref SPScountsonly "all summary statistics for the data in each node, or just count data".
If a value for holdAllStats is not specified in the AdaptiveHistogram
constructor, holdAllStats will be set to false so that, by default, only
count data is maintained.

The image belows shows a partition of a [0, 1]x[0, 1] data space into a
subpaving and the data points in each box of the subpaving, together with
the tree representation with the counts of data in the boxes represented
by each node of the tree.  In the image, the nodes of the tree are
identified with their \ref SPnodename "nodeNames".

\image html ADHExampleSimpleSPandTree.png "A simple 2-dimensional histogram partition with data as subpaving and tree" width=10cm
\image latex ADHExampleSimpleSPandTree.png "A simple 2-dimensional histogram partition with data as subpaving and tree" width=10cm

The equivalent normalised histogram is shown below.  The height \f$ h_j \f$ 
associated with bin \f$ j \f$ is \f$ \frac{n_j}{Nvol_j} \f$ where \f$ n_j \f$ 
is the number of data points associated with bin j, \f$ vol_j \f$ is the
volume of bin j, and \f$ N \f$ is the total number of data points in the
histogram, \f$ N = \sum_{bins j} n_j \f$. Thus
 
\f[
sum_{bins j} h_jvol_j = 1 \f]

\image html ADHExampleSimpleHist.png "The equivalent histogram" width=10cm
\image latex ADHExampleSimpleHist.png "The equivalent histogram" width=10cm

The \link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink class declarations and definitions can be found in
adaptivehistogram.hpp and adaptivehistogram.cpp.

<HR>

\section adhsec_tree Controlling the formation of the histogram partition

The purpose of the \link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink class is to allow the histogram bin widths and
bin centres to adapt in some way to the data to be represented in the histogram.

There are three ways in which a partition can be formed:

- \ref adhsubsec_pq
- \ref adhsubsec_onebyone
- \ref adhsubsec_mcmc

\subsection adhsubsec_pq Priority queue-based partitions

We can use a priority queue to form the partition starting from the point where
all the data is initially associated with a single root box (the sample data
space).  The subpaving is refined by progressively bisecting boxes using a
priority queue to determine which box, of all the boxes in the subpaving at
that point, to bisect next.  Using the tree representation the priority queue
operates on leaf nodes, prioritising leaf nodes for splitting.

The priority queue operation must be able compare SPSnodes in order to
prioritise which to act on first.  Priority queue methods of the
AdaptiveHistogram class use a function object for comparing SPSnodes.
\link subpavings::NodeCompObj NodeCompObj\endlink is an abstract base class
from which concrete classes for these function objects must be derived.  Some
examples of useful node comparison function objects can be found in the file
nodecompobj.hpp.

The priority queue operation must also have a 'stopping rule' to stop further
changes in the subpaving/its tree representation. Typically this rule uses
characteristics of the histogram as a whole (for example, the number of
leaves in the tree).  Priority queue methods of the AdaptiveHistogram class
use a function object to provide the stopping rule. HistEvalObj
is an abstract base class from which concrete classes for these function
objects must be derived.  Some examples of useful function objects to
determine when adaptive changes in the data partition should cease can be
found in the file histevalobj.hpp.

The operation of the priority queue can be further controlled by
specifiying a minimum number of data points which can be associated
with any node in the tree which will be created by the method.  Setting
this parameter to be a value > 0 will effectively remove from the queue
any node which cannot be split because that would result in a child leaf
with less than the required minimum number of data points associated with it.
The priority queue operation will cease when there are no nodes in the queue
if this occurs before the stopping rule is satisfied.  Thus 'large'
nodes (on the basis of the node comparison function) may not be split and the
tree may not yet satisfy the stopping rule when the priority queue method
has finished because of the effect of the value supplied for the minimum
number of data points to be associated with a node.

The prioritySplit method also takes a parameter minVolB which controls the
minimum volume of the boxes represented by the nodes in the tree which
will be create by the method.  If the total number of data points
associated with the whole tree is N, then the minimum volume of the box
represented by any node will be
\f$ \frac{1}{2} minVolB \frac{(logN)^2}{N} \f$.  Setting minVolB to 0.0
will mean that there is no minimm volume restriction on the boxes. The
minimum volume parameter is useful when priority queue splitting to
optimise some histogram fit scoring formula such as
Akaike's Information Criteria (AIC).

Some care must be taken to specify a compatible pairing of node comparision
and stopping rule, to ensure that the number of nodes in the queue will
not just continue to increase without the stopping rule being satisfied.

For splitting using a priority queue method documentation see
\link subpavings::AdaptiveHistogram::prioritySplit AdaptiveHistogram::prioritySplit\endlink.

Priority queues can also be used to coarsen a partition by reabsorbing pairs
of sibling leaf nodes back into their parent node.  This is known as 'merging'
a subpaving.  For merging using a priority queue method documentation see
\link subpavings::AdaptiveHistogram::priorityMerge AdaptiveHistogram::priorityMerge\endlink.

For examples using priority queue approaches see \ref adhsubsec_exambivg
and \ref adhsubsec_examLevy


\subsection adhsubsec_onebyone Partition refinement during data insertion

The histogram partition can be progressively refined as each data point is
inserted if the rules for forming the partition can be applied directly to
each node in the tree representing the subpaving.  For example, a
statistically equivalent blocks (SEB) partition which has a maximum of some
specified number of data points associated with each bin can be formed in
this way.  Each leaf node representing a box in the subpaving will be split
when the number of data points associated with that node exceeds the maximum,
irrespective of the attributes of other nodes or of tree as a whole.

Since each SPSnode in the tree can behave autonomously in this situation,
the decision on whether to split is 'delegated' to the SPSnode objects.
The SPSnode class provides
a method for associating data with a node which takes as a parameter a
function object which will direct whether the node should split after the data
is added.  \link subpavings::SplitDecisionObj SplitDecisionObj\endlink is an
abstract base class from which concrete classes for these function objects
must be derived.  Specific splitting rules can be encoded using these
function objects and then used by the method for inserting data into an SPSnode.
The total regular subpaving/its tree representation will develop accordingly.
Some examples of useful function objects that can be used by nodes to
determine when to split can be found in the file splitdecisionobj.hpp.

See \ref adhsubsec_exambivg for an example using this kind of process for partition formation.

\subsection adhsubsec_mcmc Partition changes in a Markov Chain using a Monte Carlo algorithm (MCMC)

The AdaptiveHistogram class can probabilistically change its partition to form different
histogram states in a MCMC process.  For more details see 
\link subpavings::AdaptiveHistogram::MCMC AdaptiveHistogram::MCMC\endlink
and \link subpavings::AdaptiveHistogram::MCMCsamples AdaptiveHistogram::MCMCsamples\endlink
documentation.

<HR>

\section adhsec_inputdata  Data input with AdaptiveHistogram objects

The AdaptiveHistogram class is designed to deal with multi-dimensional data.
Data associated with each AdaptiveHistogram object is stored as the cxsc::rvector type.

The AdaptiveHistogram class provides methods to take data from appropriately
formatted txt files, from a collection of doubles, from a container of rvectors,
or from an RSSample object.  If necessary, direct data insertion can also be coded in the
user program.  Some of these methods are discussed in more detail below.  Refer
to the \link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink class
documentation for more detail.

When data is input to a AdaptiveHistogram object an attempt will be made to
associate it to the \link subpavings::AdaptiveHistogram::rootPaving StatsSubPaving\endlink
managed by the histogram.   If the data does not fit in the root box of this
StatsSubPaving, it cannot go into the StatsSubPaving.  Any data successfully
associated with the StatsSubPaving managed by the histogram will be put
into the AdaptiveHistogram object's 
\link subpavings::AdaptiveHistogram::dataCollection data collection\endlink.

If the AdaptiveHistogram does not have a StatsSubPaving to manage
when the data is inserted, one will be created with a root box
of the appropriate size and dimensions to fit all the input data
matching the dimensions expected by the input method used.

\subsection adhsubsec_txtinput Data input from a txt file

\subsubsection adhsubsubsec_inputrvectors Input one-dimensinal or multi-dimensional data from txt file

This method reads in lines of data representing rvectors from a txt file.
The dimensions of the rvector are either specified by the user or
deduced from the input
data format.  All the data
is then expected to be in the same dimensions.  Any data line read which is less
than the expected dimensions will be rejected (if a data line with more
values than the expected dimensions is read, the excess values will be
silently ignored).

See \ref adhsubsec_exambivg for an example creating a txt file of data and then
using this input method.

The method expects one line per rvector with the elements separated by white
space (space or tabs) or commas (","), with no non-numeric characters ('e' and 'E' are accepted as part
of a floating point format).

The method can read one-dimensional and  multidimensional data.

The method carries out some basic data checking.  Input lines which do not
pass (because they contain illegal characters, or where the data is not 
of the expected dimension) are printed to
standard error output with an error message but the entire file will continue to be
processed until the end of the file is reached.

For example
<UL>
<LI> A line "12.04 1.00005e-10 -30.0006" will be read as a 3-dimensional rvector
<LI> A line "12,1.0E-10,-30" will be read as a 3-dimensional rvector
<LI> A line "12.04a 1.00005e-10 -30.0006" will be rejected because it contains illegal characters
<LI> A line "12.04, 1.00005e-10, -30 will be rejected
<LI> A line "-30.0006" will be read as a 1-dimensional rvector
<LI> A line "30" will be read as a 1-dimensional rvector.
</UL>

For method documentation see \link subpavings::AdaptiveHistogram::insertRvectorsFromTxtOrd AdaptiveHistogram::insertRvectorsFromTxtOrd\endlink.  
Note that different
levels of checking are available with 
\link subpavings::AdaptiveHistogram::insertRvectorsFromTxtParanoid insertRvectorsFromTxtParanoid\endlink and
\link subpavings::AdaptiveHistogram::insertRvectorsFromTxtFast insertRvectorsFromTxtFast\endlink.


\subsection adhsubsec_dblvecinput Data input from a container of rvectors

Data can be input directly from a collection of containers of type
\link subpavings::VecDbl VecDbl\endlink. 

As with data input from txt files, the dimensions of the rvector are
specified by the user or deduced from the input data and all the data
is then expected to be in the same dimensions.  Any VecDbl containing
less than the expected dimensions will be rejected (if a VecDbl with more
values than the expected dimensions is found, the excess values will be
silently ignored).

For method documentation see
\link subpavings::AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls AdaptiveHistogram::insertRvectorsFromVectorOfVecDbls\endlink.

This method may be useful when reading in output from other programs which
do not use the cxsc::rvector type.

\subsection adhsubsec_rvecinput Data input from a container of rvectors

Data as rvectors can be input directly from a container of type
\link subpavings::RVecData RVecData\endlink. For method documentation see
\link subpavings::AdaptiveHistogram::insertFromRVec AdaptiveHistogram::insertFromRVec\endlink.

For an example using this data input method see \ref adhsubsec_examaveraging

A sub-sample of the data in a container of type
RVecData can also be input to the histogram. See
\link subpavings::AdaptiveHistogram::insertSampleFromRVec AdaptiveHistogram::insertSampleFromRVec\endlink.
This method may be useful when averaging over bootstrapped sub-samples
from a data sample.

\subsection adhsubsec_rssampleinput Data input from a RSSample object

This method takes data from the samples held by an RSSample object.  The
RSSample object will have been created using rejection sampling (see \ref moorerejsam ).

For method documentation see \link subpavings::AdaptiveHistogram::insertFromRSSample AdaptiveHistogram::insertFromRSSample\endlink.

A sub-sample of the data in %RSSample can also be input to the
histogram.  See \link subpavings::AdaptiveHistogram::insertSampleFromRSSample AdaptiveHistogram::insertSampleFromRSSample\endlink.
This method may be useful when averaging over bootstrapped sub-samples
from a data sample.

See \ref adhsubsec_examLevy for an example using averaging over bootstrapped samples from a %RSSample.

\subsection adhsubsec_directinput Direct data input in the user code

The user can customise data input using the AdaptiveHistogram::insertOne() method.
See \ref adhsubsec_exambyhand for an example.

<HR>


\section adhsec_output AdaptiveHistogram output

A summary of a histogram can be output to a txt as either a tab-delimited file
of numeric data or a space-delimited mixture of alphanumeric data.

\subsection adhsubsec_outputtabs Tabbed numeric output to a txt file

Data for the leaf nodes of the tree managed by the AdaptiveHistogram object
is output to a .txt file using the methods 
\link subpavings::AdaptiveHistogram::outputToTxtTabs AdaptiveHistogram::outputToTxtTabs\endlink
or
\link subpavings::AdaptiveHistogram::outputToTxtTabsWithEMPs AdaptiveHistogram::outputToTxtTabsWithEMPs\endlink.
The leaf nodes of the tree
represent boxes (bins) in the subpaving (partition) of the data space.

\link subpavings::AdaptiveHistogram::outputToTxtTabs outputToTxtTabs\endlink produces a tab-delimited file of numeric data
starting with nodeName, then the node box volume, then the node counter,
then the description of the node box as a tab-delimited list of interval upper
and lower bounds, for each leaf node in the tree representing the subpaving
managed by the AdaptiveHistogram object.

The \ref SPnodename "nodeName" is a name for the node which summarises
its local place in the tree (i.e. its immediate parent and whether it is the
left of right child of that parent).

The \link subpavings::SPSnode::counter counter\endlink is a summary of the
number of datapoints associated with the box represented by that node.

The description of the box represented by the node gives the upper and lower
bound of each interval in the \ref intervalvectors "interval vector"
describing that box.

The output format format when the node represents an n-dimensional box
formed of intervals interval_1, interval_2, ...  interval_n is

nodeName   counter   volume  Inf(interval_1) Sup(interval_1) ...  Inf(interval_n)
Sup(interval_n)

The volume of the box is the product of the widths of the intervals comprising
the interval vector defining the box.

Each line of the output file will contain data relating to one leaf node.

This format is designed for flexible subsequent processing, either for further
manipulation or for graphical display.

For example, a section of .txt file output for 2-dimensional data is

\verbatim
XLLLL   2.48672 3    -3.295435   -1.644448   -2.686375   -1.180173
XLLLRL  1.24336 2    -3.295435   -2.469941   -1.180173    0.326030
XLLLRR  1.24336 16   -2.469941   -1.644448   -1.180173    0.326030
XLLRLL  1.24336 13   -1.644448   -0.818954   -2.686375   -1.180173
XLLRLRL 0.62168 9    -0.818954    0.006540   -2.686375   -1.933274
\endverbatim

For method documentation see 
\link subpavings::AdaptiveHistogram::outputToTxtTabs AdaptiveHistogram::outputToTxtTabs\endlink,

Information about the contribution of each node to the empirical risk of the
histogram as a non-parametric density estimate under scoring methods
AIC and COPERR can also be output using 
\link subpavings::AdaptiveHistogram::outputToTxtTabsWithEMPs AdaptiveHistogram::outputToTxtTabsWithEMPs\endlink.

\subsection adhsubsec_outputroot Information on the whole sample to txt file

Summary information about the whole <b> data sample </b> can be output to a txt file using
the method \link subpavings::AdaptiveHistogram::outputRootToTxt AdaptiveHistogram::outputRootToTxt\endlink.
This summarises all the data associated with
the root box of the subpaving (the domain of all
the bins in the histogram), its volume, the total number of datapoints
associated with the histogram, the mean, and the sample
variance-covariance matrix (the mean and variance-covariance matrix can only
be output if the tree managed by the histogram is maintaining all statistics
in its nodes, not just counts).  The information is a mixture of alpha and
numeric characters and is separated by spaces, not tabs.


\subsection adhsubsec_outputroot Information on the whole histogram to console output

Information about the <b> whole histogram </b> can be printed to console output.  This
shows, for each node in the binary tree representing the subpaving whose
leaves form the bins of the histogram, the subpaving box and its volume,
the count, mean and variance-covariance matrix for the data associated
with that node, and the - for leaf nodes only - the actual data.  (The
mean and variance-covariance matrix can only be output if the tree managed
by the histogram is maintaining all statistics in its nodes, not just counts.)
It is best to use this output method only for small histograms with small
samples, since the amount of output is considerable.

For method documentation see operator<<(std::ostream &os, const AdaptiveHistogram& adh).


<HR>

\section adhsec_averaging Averaging AdaptiveHistograms

Summary information from a number of different AdaptiveHistogram objects
can be collated and averaged using the
\link subpavings::AdaptiveHistogramCollator AdaptiveHistogramCollator\endlink
class <B> provided that the statistical %subpavings represented by 
\link subpavings::SPSnode SPSnode\endlink trees managed by each
AdaptiveHistogram collated all have identical root boxes</B>.

The AdaptiveHistogramCollator
creates a subpaving which is the union of the %subpavings of all the 
AdaptiveHistograms collated and keeps a record of the height of each 
collated histogram over each box in the new, unioned, subpaving.  The subpaving is
represented by a tree of \link subpavings::CollatorSPnode CollatorSPnodes\endlink.

Recall that for any individual histogram the height \f$ h_j \f$ associated
with bin j is \f$ \frac{n_j}{Nvol_j} \f$ where \f$ n_j \f$ is the number
of data points associated with bin j, \f$ vol_j \f$ is the volume of bin j,
and \f$ N \f$ is the total number of data points in the
histogram, \f$ N = \sum_{bins j} n_j \f$.

The average histogram is formed by averaging the heights over all the
the histograms collated for each box in the union of the %subpavings of the histograms collated.

The image below shows two histograms, each on the data space
[-5,5]x[-5,5].


\image html DemoHists.png "Two histograms to demonstrate averaging" width=15cm
\image latex DemoHists.png "Two histograms to demonstrate averaging" width=15cm

We can look at the subpaving or parition of the data space for each histogram
and the union of these two %subpavings, which is also a regular %subpaving.

\image html DemoHistsPartitions.png "The union of two subpavings" width=15cm
\image latex DemoHistsPartitions.png "The union of two subpaving" width=15cm

For each original histogram, a histogram of exactly the same shape but formed over the union %subpaving can be found.  Some of the bins have been subdivided  but the height over the subdivided bin is the same as the
height over the original larger bin and so the shape is the same (and the
area of the histogram still integrates to 1).  The images below show the result of representing the histograms to get the same shape but over the
union %subpaving.

\image html DemoHistsOverUnionPartition.png "Maintaining histogram shape over the union partition" width=15cm
\image latex DemoHistsOverUnionPartition.png "Maintaining histogram shape over the union partition" width=15cm

The average histogram is found by averaging the heights of the two histograms for each bin that is part of the union partion.

\image html DemoHistsAv.png "The average histogram" width=15cm
\image latex DemoHistsAv.png "The average histogram" width=15cm

For examples using averaging and \link subpavings::AdaptiveHistogramCollator :AdaptiveHistogramCollator\endlink objects, 
see \ref adhsubsec_examaveraging and \ref adhsubsec_examLevy.

<HR>


\section adhsec_examples Examples using adaptive histograms

- \ref adhsubsec_exambivg
- \ref adhsubsec_examaveraging
- \ref adhsubsec_exambyhand

\subsection adhsubsec_exambivg An AdaptiveHistogram with Bivariate Gaussian data

In this example we show how an AdaptiveHistogram object can be used to process
sample data from a Bivariate Gaussian distribution, ie two-dimensional data.

Our implementation of this example is in the file BiGTest.cpp, which includes the
header file dataprep.hpp.

Dataprep.hpp specifies libaries for generating simulated random
samples from statistical distributions

\anchor ADHexamdataprep
\dontinclude dataprep.hpp

\skip #ifndef
\until #endif

In BiGTest.cpp we include the other header files and libraries we want to be able to use.
The file header histall.hpp includes the headers usually needed for programming
using %AdaptiveHistograms.  We also specify the namespaces we want to be able to use
without qualification.

\dontinclude BiGTest.cpp

\skip #include
\until {

The first part of the file is concerned with setting the up the objects to generate simulated random samples from a bivariate gaussian distribution.

\skip generate
\until gsl_rng_alloc

In this example, we put the samples generated into a txt file which will then be
used to input data into the histogram.  Note that when we create an ofstream object for the file to send the samples to, we specify the formatting for this to be scientific with precision.  This ensures that each value output will include a decimal point.


\skip samplesFileName
\until gsl_rng_free

These are then output into a txt file.  Note the format used in the section where we output the data to a txt file:  Each pair of numbers is on the same line, separated by whitespace (in this case, just space, but tabs would also work).



\skip itx = &x[0];
\until cout

We also set up a way to time the processes, and some boolean variables used later in the process


\skip clock_t
\until successfulPQSplit


The first example in this file will show how to create a histogram with a rule to adjust the partition (split a node in the tree representing the subpaving partition) when the number of data points associated with the bin is greater than some number k.  This is a simple implementation of an upper-bounded version of the statistically equivalent blocks (SEB) principle or test or rule.


\skip example
\until user

Pretending that we did not know how many data points were in the sample, we can find this out using the countLinesInTxt method and get an appropriate input value for k.


\skip count
\until input


Now we make the AdaptiveHistogram object, in this case rather unimaginatively called myHistFirst.


\skip make
\until myHistFirst

myHist is given no initial values. This means that we have not actually specified a root box for the subpaving represented by the tree of SPSnodes that myHistFirst manages. This will be dealt with later by myHistFirst which will create a box tailor-made to fit all the data.

\note This is the safest approach when dealing with data if you are not sure of its magnitude since data which does not fit in the root box of the StatsSubPaving cannot be processed by that StatsSubPaving.

By default, the histogram data member \ref ADHholdallstats "holdAllStats" is false.

Then we start a clock, to time how long processing the data takes.

\skip start=clock()
\until running

Now we create the function object to direct whether or not a node is split after a data point is associated with it.


\skip function
\until SplitOnK

See the class documentation for \link subpavings::SplitOnK SplitOnK\endlink
for more information on this object.  Note that the minimum number of points
to be associated with a node defaults to 0 in this example.  If a minimum 
number of points to be associated with a node is specified in the
SplitOnK constructor this would effectively override the specfified 
maximum k where there was a conflict between the two requirements.  
All nodes would have a least the minimum number of points but might 
have more than k data points associated with them because splitting 
the node would violate the minimum.

Next we input the data from the txt file of samples to myHistFirst, using the
\link subpavings::AdaptiveHistogram::insertRvectorsFromTxtOrd insertRvectorsFromTxtOrd\endlink
method and supplying the function object splitK we 
have just created.  We also specify that there will be no logging of the insertion process
(see the enumeration \link subpavings::LOGGING_LEVEL LOGGING_LEVEL\endlink
for possible values for the logging parameter).

\skip insert
\until splitK

This method returns a boolean (true or false) value indicating whether or not
some data has been successfully read in, put in the
the \link subpavings::AdaptiveHistogram::dataCollection data collection\endlink, and presented to
the \link subpavings::SPSnode StatsSubPaving \endlink managed by the 
histogram. In this case, because we did not initialise myHistFirst with a 
root box for the StatsSubPaving, a root box specifically tailored to the 
data will have been created as part of the data insertion process and so 
we will not have any data left out in the cold.

We now stop the clock and print out the processing time.



\skip end
\until cout

If the data insertion process was successful, we can output the subpaving to a
txt file.  The name of the file to send the output to is given in the program.


\skip only
\until outputToTxtTabs


Then we end the if statement checking if data insertion was successful 
and have finished the first example.


\skip }
\until end


The next part of the program demonstrates the priority queue method of 
forming the histogram:  all the data is initially associated with the 
root box of the subpaving and a priority queue is then used to refine 
the bins in the partition.



\skip successfulInsertion
\until bins

Again we make an AdaptiveHistogram object with no root box.


\skip make
\until myHistSecond

We start the clock and insert the data, but this time we do not give any
rule for splitting the subpaving as the data is inserted.  The 
AdaptiveHistogram will simply make a root box big enough for all the 
data (because we have not pre-specified the root box) and will associate 
all the data with that box.

\skip start
\until samplesFileName

Now, provided that insertion of the data was successful, we do the actual
formation of the histogram by the splitting the boxes most in need of
attention first (ie, a priority queue).

In this example we specify that nodes should be compared on the basis
of the number of data points associated with them by using a node
comparision function object of type
\link subpavings::CompCount CompCount\endlink.  This means that
the node most in need of attention in the priority queue is the
one with the largest number of data points associated with it.

We specifiy a stopping rule with a histogram evaluation function
object of type CritLeaves_GTE.  This object will stop the priority queue
when the number of leaves in the tree representing the subpaving managed
by the histogram is greater than or equal to 50.

The example also specifies a minimum number of data points that can be
associated with any node in the tree.  This requirement will effectively
remove from the priority queue all nodes whose children would have less than
1 data point associated with them if the node were split.  See
\ref adhsubsec_pq and the AdaptiveHistogram::prioritySplit() method
documentation for more information about the operation of this parameter.

\skip if
\until end

Then we show the time taken and output the histogram.



\skip timeTaken
\until }

Free the memory allocated for the data points


\skip delete
\until delete y;

And end the program


\skip return
\until }


One of many possible displays of the result is shown below.  This draws the
histograms as 3-dimensional shapes, with the 'floor' being the boxes of
the subpaving and the height of each histogram bar being the proportion
of the total number of datapoints which were associated with that box
divided by the floor area (volume of the box).  The image shows
both histograms, illustrating the smoothing effect achieved by priority
queue splitting to a fixed number of leaves (bins in the partition).


\image html bivgaussian.png "Graphical representations of AdaptiveHistograms" width=8cm
\image latex bivgaussian.png "Graphical representations of AdaptiveHistograms" width=8cm


\subsection adhsubsec_examaveraging Example of averaging histograms

In this example we show how an \link subpavings::AdaptiveHistogramCollator AdaptiveHistogramCollator\endlink
object can be used
find the average histogram from ten different histograms (AdaptiveHistogram
objects), each managing a sample of bivariate gaussian data.

Our implementation of this example is in the file Averaging.cpp, which includes the
header file dataprep.hpp.

\ref ADHexamdataprep "Dataprep.hpp" specifies libaries for generating simulated random
samples from statistical distributions

In Averaging.cpp we also include the other header files and libraries
we want to be able to use. The file header histall.hpp includes the
headers usually needed for programming
using AdaptiveHistograms.  We also specify the namespaces we want to be able to use
without qualification.

\dontinclude Averaging.cpp

\skip #include
\until {

The first part of the file is concerned with setting the file up to generate
and store simulated random samples from a bivariate gaussian distribution.

\skip generate
\until gsl_rng_alloc

The example is described

\skip example
\until multiple

It is vital that each %AdaptiveHistogram to be collated has the same root box, so we make a large root box which will be used by all the histograms we subsequently make.

\skip make
\until pavingBox[2]

Then we make the AdaptiveHistogramCollator object

\skip make
\until coll;

We have a loop for each histogram to be averaged, and within the loop we
start by generating sample data for the histogram.  Note that the sample
data is stored in a \link subpavings::RVecData RVecData\endlink container.

\skip histograms
\until }


We make an AdaptiveHistogram object to be used with this sample.
Each histogram will use the a copy of the same root box for its
subpaving.

\skip Histogram
\until myHist(pavingBox)

Each AdaptiveHistogram is to be formed under the rule where each bin
should have a maximum number of data points k associated with it, but
that maximum number of data points is a function of the total number
of data points n and the histogram indexer j.  This applies SEB
heuristics for k to satisfy 
\f$ \frac{k}{n} \rightarrow 0 \f$ as \f$ n \rightarrow \infty \f$.
The resulting value for k is used in the constructor for an object of
type SplitOnK.  See the class documentation for 
\link subpavings::SplitOnK SplitOnK\endlink for more information on
this object.  Note that the minimum number of points to be associated
with a node defaults to 0.

\skip find
\until splitK

Next we input the sample data from the container into the histogram, using the
insertFromRVec method and supplying the function object splitK we have just
created.  We also specify that there will be no logging of the insertion process
(see the enumeration \link subpavings::LOGGING_LEVEL LOGGING_LEVEL\endlink
for possible values for the logging parameter).

Note that we choose to split with each data insertion here.  A priority
queue with a stopping rule based on the maximum number of points
associated with any any leaf node of the tree (CritLargestCount_LTE) could also be used.

\skip insert
\until successfulInsertion

Assuming that the data insertion has been successful, we output each histogram to a txt file and add it to the collation.

Adding the histogram to the collation does not change the AdaptiveHistogram object itself, it just adds information from that object to the collation.

A dot file of the tree representation of the histogram can also be created but works best for smaller histograms.

\skip if
\until }

This ends the loop for creating separate histogram objects.  The random number generator can be freed.

\skip }
\until gsl_rng_free

And the collated information output to a txt file.  A dot graph can also be created if desired.

\skip collfileName
\until outputGraphDot

Finally, the average is also output to a txt file.

\skip Make
\until outputToTxtTabs

The example ends.

\skip end
\until program

The ten histograms constructed with different values of k based on the ten
independent simulations of 10000 bivariate Gaussian samples are diplayed next
along with the prabability density of the standard bivariate Gaussian.

\image html BivGaussHists10.png "Graphical representations of the 10 AdaptiveHistograms" width=8cm
\image latex BivGaussHists10.png "Graphical representations of the 10 AdaptiveHistograms" width=8cm

The average histogram or the sample mean histogram of the above 10 histograms
are shown below with the probability density of the standard bivariate Gaussian.

\image html BGAvgHist.png "Graphical representations of the average of the 10 AdaptiveHistograms" width=8cm
\image latex BGAvgHist.png "Graphical representations of the average of the 10 AdaptiveHistograms" width=8cm

\subsection adhsubsec_examLevy Example using data from a RSSample object

This example first generates data in an RSSample object and then inputs samples from this to an %AdaptiveHistogram, again averaging the resulting histograms.

The example must include the header files for the rejection sampling process

\dontinclude LevyTest.cpp

\skip #include
\until example

We start by using the Moore Rejection Sampler to generate the RSSample of data.  See \ref moorerejsam for more information.

\skip sync_with_stdio
\until RejectionSampleMany

We make a root box which will be used by all the histograms to be averaged.

\skip make
\until pavingBox[2]

We specify how many data points are to go into each sample and how many histograms to average over, then make an AdaptiveHistogramCollator to collate the separate histograms together.

\skip number
\until AdaptiveHistogramCollator

We will use a single method of the %AdaptiveHistogramCollator object to both create each histogram using priority queue splitting and collate the histograms.

For this we need to set up function objects and other arguments to be used in the priority queue operations.

In this example we specify that nodes should be compared on the basis of the number of data points associated with them by using a node comparision function object of type
\link subpavings::CompCount CompCount\endlink.  This means that the node most in need of attention in the priority queue is the one with the largest number of data points associated with it.

We specifiy a stopping rule with a histogram evaluation function object of type CritSmallestVol_LTE.  This object is constructed to stop a priority queue when the volume of the smallest box in the subpaving (partition) is less than 0.05.

We could specify a minimum number of data points that can be associated with any node in a tree.  In this example we set this to 0 which will have no effect on the priority queue operations.  Similarly we use minVolB = 0.0 to specify no minimum volume for the boxes represented by nodes in the trees. See \ref adhsubsec_pq and the AdaptiveHistogram::prioritySplit() method documentation for more information about the operation of these arguements.

\skip priority
\until double minVolB = 0.0;

The actual creation of the 'sample' histograms and their collation into the %AdaptiveHistogramCollator is carried out using the
\link subpavings::AdaptiveHistogramCollator::collateFromRSSampleSplitPQ collateFromRSSampleSplitPQ\endlink
method.  If the method is successful, it will output a separate txt file for each 'sample' histogram.
The collated results will be in the AdaptiveHistogramCollator we have 
named "coll" (this is output to the txt file CollatorHistogram.txt).  The
average is also created and output to AverageHistogram.txt.

\skip success
\until cout

The example program then ends.

\skip return
\until program

\subsection adhsubsec_exambyhand Example of insertion of data defined in the program

This is an example of a section of code which sets up an \link subpavings::AdaptiveHistogram AdaptiveHistogram\endlink
object by specifying the root box (the interval [-40.0 40.0]) as well as the
splitting value (2) and then defines the data within the program code and
feeds each value individually to the AdaptiveHistogram.  It is unlikely that
this method would be useful in most practical circumstances but the example is
included to illustrate the use of the pre-specified root box and the
insertOne() method.

\code
  ivector x(1);
  x[1] = interval(-40.0, 40.0);

  AdaptiveHistogram myHist(x);

  size_t k = 2;
  // function object to split when > 2 data points associated with a node
  SplitOnK splitDecision(k);

  real newreal = -3.0;
  // one dimensional rvector, only element is newreal
  myHist.insertOne(newdata, splitDecision);

  newreal = 14.0;
  newdata = _rvector(newreal); // another one
  myHist.insertOne(newdata, splitDecision);

  newreal = 6.0;
  newdata = _rvector(newreal); // another one
  myHist.insertOne(newdata, splitDecision);

  newreal = -10.0;
  newdata = _rvector(newreal); // another one
  myHist.insertOne(newdata, splitDecision);

  newreal = -3.0;
  newdata = _rvector(newreal); // another one
  myHist.insertOne(newdata);

\endcode




<HR>

\section adhsec_furtherdev Further development
  The AdaptiveHistogram class has not yet been fully developed.  More methods
  and functionality can be usefully added.

\todo More flexible user friendly input and output
\todo Generate graphical output for histograms

*/

/*! \page GFDL GNU Free Documentation License

\htmlonly

<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html><head>
 <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
 <title>GNU Free Documentation License - GNU Project - Free Software Foundation (FSF)</title>
</head>
<body>

<h3 style="text-align: center;">GNU Free Documentation License</h3>

<p style="text-align: center;">Version 1.3, 3 November 2008</p>

<p> Copyright &copy; 2000, 2001, 2002, 2007, 2008 Free Software Foundation, Inc.
     &lt;http://fsf.org/&gt;
 </p><p>Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.</p>

<h4><a name="section0"></a>0. PREAMBLE</h4>

<p>The purpose of this License is to make a manual, textbook, or other
functional and useful document &quot;free&quot; in the sense of freedom: to
assure everyone the effective freedom to copy and redistribute it,
with or without modifying it, either commercially or noncommercially.
Secondarily, this License preserves for the author and publisher a way
to get credit for their work, while not being considered responsible
for modifications made by others.</p>

<p>This License is a kind of &quot;copyleft&quot;, which means that derivative
works of the document must themselves be free in the same sense.  It
complements the GNU General Public License, which is a copyleft
license designed for free software.</p>

<p>We have designed this License in order to use it for manuals for free
software, because free software needs free documentation: a free
program should come with manuals providing the same freedoms that the
software does.  But this License is not limited to software manuals;
it can be used for any textual work, regardless of subject matter or
whether it is published as a printed book.  We recommend this License
principally for works whose purpose is instruction or reference.</p>

<h4><a name="section1"></a>1. APPLICABILITY AND DEFINITIONS</h4>

<p>This License applies to any manual or other work, in any medium, that
contains a notice placed by the copyright holder saying it can be
distributed under the terms of this License.  Such a notice grants a
world-wide, royalty-free license, unlimited in duration, to use that
work under the conditions stated herein.  The &quot;Document&quot;, below,
refers to any such manual or work.  Any member of the public is a
licensee, and is addressed as &quot;you&quot;.  You accept the license if you
copy, modify or distribute the work in a way requiring permission
under copyright law.</p>

<p>A &quot;Modified Version&quot; of the Document means any work containing the
Document or a portion of it, either copied verbatim, or with
modifications and/or translated into another language.</p>

<p>A &quot;Secondary Section&quot; is a named appendix or a front-matter section of
the Document that deals exclusively with the relationship of the
publishers or authors of the Document to the Document's overall
subject (or to related matters) and contains nothing that could fall
directly within that overall subject.  (Thus, if the Document is in
part a textbook of mathematics, a Secondary Section may not explain
any mathematics.)  The relationship could be a matter of historical
connection with the subject or with related matters, or of legal,
commercial, philosophical, ethical or political position regarding
them.</p>

<p>The &quot;Invariant Sections&quot; are certain Secondary Sections whose titles
are designated, as being those of Invariant Sections, in the notice
that says that the Document is released under this License.  If a
section does not fit the above definition of Secondary then it is not
allowed to be designated as Invariant.  The Document may contain zero
Invariant Sections.  If the Document does not identify any Invariant
Sections then there are none.</p>

<p>The &quot;Cover Texts&quot; are certain short passages of text that are listed,
as Front-Cover Texts or Back-Cover Texts, in the notice that says that
the Document is released under this License.  A Front-Cover Text may
be at most 5 words, and a Back-Cover Text may be at most 25 words.</p>

<p>A &quot;Transparent&quot; copy of the Document means a machine-readable copy,
represented in a format whose specification is available to the
general public, that is suitable for revising the document
straightforwardly with generic text editors or (for images composed of
pixels) generic paint programs or (for drawings) some widely available
drawing editor, and that is suitable for input to text formatters or
for automatic translation to a variety of formats suitable for input
to text formatters.  A copy made in an otherwise Transparent file
format whose markup, or absence of markup, has been arranged to thwart
or discourage subsequent modification by readers is not Transparent.
An image format is not Transparent if used for any substantial amount
of text.  A copy that is not &quot;Transparent&quot; is called &quot;Opaque&quot;.</p>

<p>Examples of suitable formats for Transparent copies include plain
ASCII without markup, Texinfo input format, LaTeX input format, SGML
or XML using a publicly available DTD, and standard-conforming simple
HTML, PostScript or PDF designed for human modification.  Examples of
transparent image formats include PNG, XCF and JPG.  Opaque formats
include proprietary formats that can be read and edited only by
proprietary word processors, SGML or XML for which the DTD and/or
processing tools are not generally available, and the
machine-generated HTML, PostScript or PDF produced by some word
processors for output purposes only.</p>

<p>The &quot;Title Page&quot; means, for a printed book, the title page itself,
plus such following pages as are needed to hold, legibly, the material
this License requires to appear in the title page.  For works in
formats which do not have any title page as such, &quot;Title Page&quot; means
the text near the most prominent appearance of the work's title,
preceding the beginning of the body of the text.</p>

<p>The &quot;publisher&quot; means any person or entity that distributes copies of
the Document to the public.</p>

<p>A section &quot;Entitled XYZ&quot; means a named subunit of the Document whose
title either is precisely XYZ or contains XYZ in parentheses following
text that translates XYZ in another language.  (Here XYZ stands for a
specific section name mentioned below, such as &quot;Acknowledgements&quot;,
&quot;Dedications&quot;, &quot;Endorsements&quot;, or &quot;History&quot;.)  To &quot;Preserve the Title&quot;
of such a section when you modify the Document means that it remains a
section &quot;Entitled XYZ&quot; according to this definition.</p>

<p>The Document may include Warranty Disclaimers next to the notice which
states that this License applies to the Document.  These Warranty
Disclaimers are considered to be included by reference in this
License, but only as regards disclaiming warranties: any other
implication that these Warranty Disclaimers may have is void and has
no effect on the meaning of this License.</p>

<h4><a name="section2"></a>2. VERBATIM COPYING</h4>

<p>You may copy and distribute the Document in any medium, either
commercially or noncommercially, provided that this License, the
copyright notices, and the license notice saying this License applies
to the Document are reproduced in all copies, and that you add no
other conditions whatsoever to those of this License.  You may not use
technical measures to obstruct or control the reading or further
copying of the copies you make or distribute.  However, you may accept
compensation in exchange for copies.  If you distribute a large enough
number of copies you must also follow the conditions in section 3.</p>

<p>You may also lend copies, under the same conditions stated above, and
you may publicly display copies.</p>

<h4><a name="section3"></a>3. COPYING IN QUANTITY</h4>

<p>If you publish printed copies (or copies in media that commonly have
printed covers) of the Document, numbering more than 100, and the
Document's license notice requires Cover Texts, you must enclose the
copies in covers that carry, clearly and legibly, all these Cover
Texts: Front-Cover Texts on the front cover, and Back-Cover Texts on
the back cover.  Both covers must also clearly and legibly identify
you as the publisher of these copies.  The front cover must present
the full title with all words of the title equally prominent and
visible.  You may add other material on the covers in addition.
Copying with changes limited to the covers, as long as they preserve
the title of the Document and satisfy these conditions, can be treated
as verbatim copying in other respects.</p>

<p>If the required texts for either cover are too voluminous to fit
legibly, you should put the first ones listed (as many as fit
reasonably) on the actual cover, and continue the rest onto adjacent
pages.</p>

<p>If you publish or distribute Opaque copies of the Document numbering
more than 100, you must either include a machine-readable Transparent
copy along with each Opaque copy, or state in or with each Opaque copy
a computer-network location from which the general network-using
public has access to download using public-standard network protocols
a complete Transparent copy of the Document, free of added material.
If you use the latter option, you must take reasonably prudent steps,
when you begin distribution of Opaque copies in quantity, to ensure
that this Transparent copy will remain thus accessible at the stated
location until at least one year after the last time you distribute an
Opaque copy (directly or through your agents or retailers) of that
edition to the public.</p>

<p>It is requested, but not required, that you contact the authors of the
Document well before redistributing any large number of copies, to
give them a chance to provide you with an updated version of the
Document.</p>

<h4><a name="section4"></a>4. MODIFICATIONS</h4>

<p>You may copy and distribute a Modified Version of the Document under
the conditions of sections 2 and 3 above, provided that you release
the Modified Version under precisely this License, with the Modified
Version filling the role of the Document, thus licensing distribution
and modification of the Modified Version to whoever possesses a copy
of it.  In addition, you must do these things in the Modified Version:</p>

<ul>


<li>A. Use in the Title Page (and on the covers, if any) a title distinct
   from that of the Document, and from those of previous versions
   (which should, if there were any, be listed in the History section
   of the Document).  You may use the same title as a previous version
   if the original publisher of that version gives permission.
</li>

<li>B. List on the Title Page, as authors, one or more persons or entities
   responsible for authorship of the modifications in the Modified
   Version, together with at least five of the principal authors of the
   Document (all of its principal authors, if it has fewer than five),
   unless they release you from this requirement.
</li>

<li>C. State on the Title page the name of the publisher of the
   Modified Version, as the publisher.
</li>

<li>D. Preserve all the copyright notices of the Document.
</li>

<li>E. Add an appropriate copyright notice for your modifications
   adjacent to the other copyright notices.
</li>

<li>F. Include, immediately after the copyright notices, a license notice
   giving the public permission to use the Modified Version under the
   terms of this License, in the form shown in the Addendum below.
</li>

<li>G. Preserve in that license notice the full lists of Invariant Sections
   and required Cover Texts given in the Document's license notice.
</li>

<li>H. Include an unaltered copy of this License.
</li>

<li>I. Preserve the section Entitled &quot;History&quot;, Preserve its Title, and add
   to it an item stating at least the title, year, new authors, and
   publisher of the Modified Version as given on the Title Page.  If
   there is no section Entitled &quot;History&quot; in the Document, create one
   stating the title, year, authors, and publisher of the Document as
   given on its Title Page, then add an item describing the Modified
   Version as stated in the previous sentence.
</li>

<li>J. Preserve the network location, if any, given in the Document for
   public access to a Transparent copy of the Document, and likewise
   the network locations given in the Document for previous versions
   it was based on.  These may be placed in the &quot;History&quot; section.
   You may omit a network location for a work that was published at
   least four years before the Document itself, or if the original
   publisher of the version it refers to gives permission.
</li>

<li>K. For any section Entitled &quot;Acknowledgements&quot; or &quot;Dedications&quot;,
   Preserve the Title of the section, and preserve in the section all
   the substance and tone of each of the contributor acknowledgements
   and/or dedications given therein.
</li>

<li>L. Preserve all the Invariant Sections of the Document,
   unaltered in their text and in their titles.  Section numbers
   or the equivalent are not considered part of the section titles.
</li>

<li>M. Delete any section Entitled &quot;Endorsements&quot;.  Such a section
   may not be included in the Modified Version.
</li>

<li>N. Do not retitle any existing section to be Entitled &quot;Endorsements&quot;
   or to conflict in title with any Invariant Section.
</li>

<li>O. Preserve any Warranty Disclaimers.</li>

</ul>

<p>If the Modified Version includes new front-matter sections or
appendices that qualify as Secondary Sections and contain no material
copied from the Document, you may at your option designate some or all
of these sections as invariant.  To do this, add their titles to the
list of Invariant Sections in the Modified Version's license notice.
These titles must be distinct from any other section titles.</p>

<p>You may add a section Entitled &quot;Endorsements&quot;, provided it contains
nothing but endorsements of your Modified Version by various
parties&mdash;for example, statements of peer review or that the text has
been approved by an organization as the authoritative definition of a
standard.</p>

<p>You may add a passage of up to five words as a Front-Cover Text, and a
passage of up to 25 words as a Back-Cover Text, to the end of the list
of Cover Texts in the Modified Version.  Only one passage of
Front-Cover Text and one of Back-Cover Text may be added by (or
through arrangements made by) any one entity.  If the Document already
includes a cover text for the same cover, previously added by you or
by arrangement made by the same entity you are acting on behalf of,
you may not add another; but you may replace the old one, on explicit
permission from the previous publisher that added the old one.</p>

<p>The author(s) and publisher(s) of the Document do not by this License
give permission to use their names for publicity for or to assert or
imply endorsement of any Modified Version.</p>

<h4><a name="section5"></a>5. COMBINING DOCUMENTS</h4>

<p>You may combine the Document with other documents released under this
License, under the terms defined in section 4 above for modified
versions, provided that you include in the combination all of the
Invariant Sections of all of the original documents, unmodified, and
list them all as Invariant Sections of your combined work in its
license notice, and that you preserve all their Warranty Disclaimers.</p>

<p>The combined work need only contain one copy of this License, and
multiple identical Invariant Sections may be replaced with a single
copy.  If there are multiple Invariant Sections with the same name but
different contents, make the title of each such section unique by
adding at the end of it, in parentheses, the name of the original
author or publisher of that section if known, or else a unique number.
Make the same adjustment to the section titles in the list of
Invariant Sections in the license notice of the combined work.</p>

<p>In the combination, you must combine any sections Entitled &quot;History&quot;
in the various original documents, forming one section Entitled
&quot;History&quot;; likewise combine any sections Entitled &quot;Acknowledgements&quot;,
and any sections Entitled &quot;Dedications&quot;.  You must delete all sections
Entitled &quot;Endorsements&quot;.</p>

<h4><a name="section6"></a>6. COLLECTIONS OF DOCUMENTS</h4>

<p>You may make a collection consisting of the Document and other
documents released under this License, and replace the individual
copies of this License in the various documents with a single copy
that is included in the collection, provided that you follow the rules
of this License for verbatim copying of each of the documents in all
other respects.</p>

<p>You may extract a single document from such a collection, and
distribute it individually under this License, provided you insert a
copy of this License into the extracted document, and follow this
License in all other respects regarding verbatim copying of that
document.</p>

<h4><a name="section7"></a>7. AGGREGATION WITH INDEPENDENT WORKS</h4>

<p>A compilation of the Document or its derivatives with other separate
and independent documents or works, in or on a volume of a storage or
distribution medium, is called an &quot;aggregate&quot; if the copyright
resulting from the compilation is not used to limit the legal rights
of the compilation's users beyond what the individual works permit.
When the Document is included in an aggregate, this License does not
apply to the other works in the aggregate which are not themselves
derivative works of the Document.</p>

<p>If the Cover Text requirement of section 3 is applicable to these
copies of the Document, then if the Document is less than one half of
the entire aggregate, the Document's Cover Texts may be placed on
covers that bracket the Document within the aggregate, or the
electronic equivalent of covers if the Document is in electronic form.
Otherwise they must appear on printed covers that bracket the whole
aggregate.</p>

<h4><a name="section8"></a>8. TRANSLATION</h4>

<p>Translation is considered a kind of modification, so you may
distribute translations of the Document under the terms of section 4.
Replacing Invariant Sections with translations requires special
permission from their copyright holders, but you may include
translations of some or all Invariant Sections in addition to the
original versions of these Invariant Sections.  You may include a
translation of this License, and all the license notices in the
Document, and any Warranty Disclaimers, provided that you also include
the original English version of this License and the original versions
of those notices and disclaimers.  In case of a disagreement between
the translation and the original version of this License or a notice
or disclaimer, the original version will prevail.</p>

<p>If a section in the Document is Entitled &quot;Acknowledgements&quot;,
&quot;Dedications&quot;, or &quot;History&quot;, the requirement (section 4) to Preserve
its Title (section 1) will typically require changing the actual
title.</p>

<h4><a name="section9"></a>9. TERMINATION</h4>

<p>You may not copy, modify, sublicense, or distribute the Document
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense, or distribute it is void, and
will automatically terminate your rights under this License.</p>

<p>However, if you cease all violation of this License, then your license
from a particular copyright holder is reinstated (a) provisionally,
unless and until the copyright holder explicitly and finally
terminates your license, and (b) permanently, if the copyright holder
fails to notify you of the violation by some reasonable means prior to
60 days after the cessation.</p>

<p>Moreover, your license from a particular copyright holder is
reinstated permanently if the copyright holder notifies you of the
violation by some reasonable means, this is the first time you have
received notice of violation of this License (for any work) from that
copyright holder, and you cure the violation prior to 30 days after
your receipt of the notice.</p>

<p>Termination of your rights under this section does not terminate the
licenses of parties who have received copies or rights from you under
this License.  If your rights have been terminated and not permanently
reinstated, receipt of a copy of some or all of the same material does
not give you any rights to use it.</p>

<h4><a name="section10"></a>10. FUTURE REVISIONS OF THIS LICENSE</h4>

<p>The Free Software Foundation may publish new, revised versions of the
GNU Free Documentation License from time to time.  Such new versions
will be similar in spirit to the present version, but may differ in
detail to address new problems or concerns.  See
http://www.gnu.org/copyleft/.</p>

<p>Each version of the License is given a distinguishing version number.
If the Document specifies that a particular numbered version of this
License &quot;or any later version&quot; applies to it, you have the option of
following the terms and conditions either of that specified version or
of any later version that has been published (not as a draft) by the
Free Software Foundation.  If the Document does not specify a version
number of this License, you may choose any version ever published (not
as a draft) by the Free Software Foundation.  If the Document
specifies that a proxy can decide which future versions of this
License can be used, that proxy's public statement of acceptance of a
version permanently authorizes you to choose that version for the
Document.</p>

<h4><a name="section11"></a>11. RELICENSING</h4>

<p>&quot;Massive Multiauthor Collaboration Site&quot; (or &quot;MMC Site&quot;) means any
World Wide Web server that publishes copyrightable works and also
provides prominent facilities for anybody to edit those works.  A
public wiki that anybody can edit is an example of such a server.  A
&quot;Massive Multiauthor Collaboration&quot; (or &quot;MMC&quot;) contained in the site
means any set of copyrightable works thus published on the MMC site.</p>

<p>&quot;CC-BY-SA&quot; means the Creative Commons Attribution-Share Alike 3.0
license published by Creative Commons Corporation, a not-for-profit
corporation with a principal place of business in San Francisco,
California, as well as future copyleft versions of that license
published by that same organization.</p>

<p>&quot;Incorporate&quot; means to publish or republish a Document, in whole or in
part, as part of another Document.</p>

<p>An MMC is &quot;eligible for relicensing&quot; if it is licensed under this
License, and if all works that were first published under this License
somewhere other than this MMC, and subsequently incorporated in whole or
in part into the MMC, (1) had no cover texts or invariant sections, and
(2) were thus incorporated prior to November 1, 2008.</p>

<p>The operator of an MMC Site may republish an MMC contained in the site
under CC-BY-SA on the same site at any time before August 1, 2009,
provided the MMC is eligible for relicensing.</p>

<h3><a name="addendum"></a>ADDENDUM: How to use this License for your documents</h3>

<p>To use this License in a document you have written, include a copy of
the License in the document and put the following copyright and
license notices just after the title page:</p>

<pre>    Copyright (C)  YEAR  YOUR NAME.
    Permission is granted to copy, distribute and/or modify this document
    under the terms of the GNU Free Documentation License, Version 1.3
    or any later version published by the Free Software Foundation;
    with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
    A copy of the license is included in the section entitled &quot;GNU
    Free Documentation License&quot;.
</pre>

<p>If you have Invariant Sections, Front-Cover Texts and Back-Cover Texts,
replace the &quot;with &hellip; Texts.&quot; line with this:</p>

<pre>    with the Invariant Sections being LIST THEIR TITLES, with the
    Front-Cover Texts being LIST, and with the Back-Cover Texts being LIST.
</pre>

<p>If you have Invariant Sections without Cover Texts, or some other
combination of the three, merge those two alternatives to suit the
situation.</p>

<p>If your document contains nontrivial examples of program code, we
recommend releasing these examples in parallel under your choice of
free software license, such as the GNU General Public License,
to permit their use in free software.
</p>

</body></html>

\endhtmlonly

*/

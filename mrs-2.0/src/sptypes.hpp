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
\brief Includes and general typedefs for subpavings

*/

#ifndef ___SPTYPES_HPP__
#define ___SPTYPES_HPP__

// to use ivectors and reals
#include "cxsc.hpp"

// to use std::list
#include <list>

// to use std::vector
#include <vector>





/*! \brief The namespace subpavings.

The namespace is used for all classes and non-member methods related to
subpavings.
*/
namespace subpavings {

    //! Forward class declarations
    class SPnode;
    class SPSnode;
    class CollatorSPnode;
	class SPMinimalnode;
	
	/*! \brief SubPaving is an alias for a pointer to an SPnode.    */
    typedef SPnode* SubPaving;

	/*! \brief StatsSubPaving is an alias for a pointer to an SPSnode.*/
    typedef SPSnode* StatsSubPaving;
	
	/*! \brief MinimalSubPaving is an alias for a pointer to a 
	SPMinimalnode.    */
    typedef SPMinimalnode* MinimalSubPaving;
	
    


    /** @name Typedefs for enums.
    */
    //@{
    /*! \brief Define type "Interval Booleans".

    This is used to extend the usual booleans TRUE and FALSE for use with
    intervals.  With intervals the result of a test may be
    indeterminate at a particular level of accuracy rather than clearly
    true or false.
    */
    typedef enum {BI_TRUE,BI_FALSE,BI_INDET} BOOL_INTERVAL;

    /*! \brief Define a type OPERATION_ON.

    This is used in SPSnode to indicate whether an operation is on a
    parent node or on the left or right child of a parent.
    */

    typedef enum{ON_LEFT = -1, ON_RIGHT = 1, ON_PARENT = 0} OPERATIONS_ON;

	//UPDATED JUNE 2012 for logposteriors
    /*! \brief Define a type LOGGING_LEVEL.

    This is used to determine the type of logging output.

	<ul>
	<li>NOLOG should be obvious.</li>
	<li>TXT for logging quite extensive txt files (pq and mcmc).</li>
	<li>GRAPHSAMPLES to graph mcmc samples.</li>
	<li>LOGSAMPLES to send each mcmc sample to a log file.</li>
	<li>LOGANDGRAPHSAMPLES for both of the above.</li>
	<li>LOGSTATETRACE to log leaf trace and log posterior trace for mcmc states.</li>
	<li>LOGMCMCTRACE to log leaf trace and log posterior trace for mcmc states
	\b and also the leaf trace and log posterior for the sample average.</li>
	</ul>
    */

    typedef enum{NOLOG = 0, TXT = 1, TXTANDGRAPH = 2, GRAPHSAMPLES = 3,
                LOGSAMPLES = 4, LOGANDGRAPHSAMPLES = 5, 
				LOGSTATETRACE = 6, LOGMCMCTRACE = 7} LOGGING_LEVEL;


    //@}

    /** @name Typedefs for function pointers.
    */
    //@{
    /*! \brief Define type "Pointer to an interval boolean test".

    The test is for containment of the inclusion function image of an
    ivector in a minimal subpaving.  The inclusion function is specified in the
    body of the interval boolean test function.  The test returns
    BI_INDET if the inclusion function image overlaps the subpaving border.

    \param x the interval vector whose image is found and compared to the
    subpaving.
    \param spn a pointer to the node representing a subpaving we want to
    test for containment in.  This parameter is added to the
    AIASPnode::(*AIA_PIBT) to replace the use of globals there.
    */
    typedef BOOL_INTERVAL (*PIBT)(const ivector& x,
                                const SPMinimalnode * const spn);

    /*! \brief Define type "Pointer to an interval vector function".

    PIVF is an interval vector inclusion function, ie a function which
    returns the interval vector which encloses f(x) for f as specified
    in the function and x given in the function parameter.
    */
    typedef ivector (*PIVF)(const ivector& x);

    //@}


    /** @name Typedefs for containers and iterators.
    */
    //@{

    /*!\brief Define type IntVec as a container for ints.
    */
    typedef std::vector<int> IntVec;

    /*!\brief Define type IntVecItr as iterator over IntVec.
    */
    typedef IntVec::iterator IntVecItr;

    /*!\brief Define type Size_tVec as a container for size_t.
    */
    typedef std::vector<size_t> Size_tVec;

    /*!\brief Define type Size_tVecItr as iterator over Size_tVec.
    */
    typedef Size_tVec::iterator Size_tVecItr;

    /*!\brief Define type RealVec as a container for reals.
    */
    typedef std::vector<real> RealVec;

    /*!\brief Define type RealVecItr as iterator over RealVec.
    */
    typedef RealVec::iterator RealVecItr;

    /*! \brief Define type VecDbl as a container of doubles.
    */
    typedef std::vector<double> VecDbl;

    /*! \brief Define type VecDblIt as an iterator over VecDbl.
    */
    typedef VecDbl::iterator VecDblIt;

    /*!\brief Define type RVecData as a container for rvectors.
    */
    typedef std::vector<rvector> RVecData;

    /*!\brief Define type RVecDataCItr as const_iterator over RVecData.
    */
    typedef RVecData::const_iterator RVecDataCItr;

    /*!\brief Define type RVecDataItr as iterator over RVecData.
    */
    typedef RVecData::iterator RVecDataItr;

    /*! \brief Define type SPnodePtrs as container of pointers to SPnodes.
    */
    typedef std::vector<SPnode*> SPnodePtrs;

    /*! \brief Define type SPnodePtrsItr as an iterator over SPnodePtrs.
    */
    typedef SPnodePtrs::iterator SPnodePtrsItr ;

	/*! \brief Define type SPnodeConstPtrs as container of pointers to const SPnodes.
    */
    typedef std::vector<const SPnode*> SPnodeConstPtrs;

    /*! \brief Define type SPnodeConstPtrsItr as an iterator over SPnodeConstPtrs.
    */
    typedef SPnodeConstPtrs::const_iterator SPnodeConstPtrsItr ;


	/*! \brief Define type SPMinimalnodePtrs as container of pointers to SPMinimalnodes.
    */
    typedef std::vector<SPMinimalnode*> SPMinimalnodePtrs;

    /*! \brief Define type SPMinimalnodePtrsItr as an iterator over SPMinimalnodePtrs.
    */
    typedef SPMinimalnodePtrs::iterator SPMinimalnodePtrsItr ;

	/*!\brief Define type BoxVec as a container of boxes.
    */
    typedef std::vector<ivector> BoxVec;

    /*!\brief Define type BoxVecItr as iterator over BoxVec.
    */
    typedef BoxVec::iterator BoxVecItr;

    /*! \brief Define type ImageList as a container for images of boxes.

    Used in Evaluate() for list of function images.
    */
    typedef std::list<ivector> ImageList;

    /*! \brief Define type iterator over ImageList.
    */
    typedef ImageList::iterator ImageListItr;

    /*! \brief Define type SPSnodePtrs as container of pointers to SPSnodes.
    */
    typedef std::vector<SPSnode*> SPSnodePtrs;

    /*! \brief Define type SPSnodePtrsItr as an iterator over SPSnodePtrs.
    */
    typedef SPSnodePtrs::iterator SPSnodePtrsItr ;

	/*! \brief Define type SPSnodeConstPtrs as container of pointers to const SPSnodes.
    */
    typedef std::vector<const SPSnode*> SPSnodeConstPtrs;

    /*! \brief Define type SPSnodeConstPtrsItr as an iterator over SPSnodeConstPtrs.
    */
    typedef SPSnodeConstPtrs::const_iterator SPSnodeConstPtrsItr ;


    /*! \brief Define type SPSnodeList as a list of pointers to SPSnodes.
    */
    typedef std::list<SPSnode*> SPSnodeList;

    /*! \brief Define type SPSnodeListItr as an iterator over SPSnodeList.
    */
    typedef SPSnodeList::iterator SPSnodeListItr ;

	
	/*! \brief Define type CollatorSPnodePtrs as container of pointers. 
	 to CollatorSPnodes.
    */
    typedef std::vector<CollatorSPnode*> CollatorSPnodePtrs;

    /*! \brief Define type CollatorSPnodePtrsItr as an iterator 
     over CollatorSPnodePtrs.
    */
    typedef CollatorSPnodePtrs::iterator CollatorSPnodePtrsItr ;

    
    /*! \brief Define type BigData Collection as a container for data.

    Used in HistogramWrappers and SPSnode as the container for sample data.
    This container must not be vulnerable to iterator, pointer and
    reference invalidation.  std::list has been chosen because it is a
    node-based container.  The disadvantages of std::list is that it is
    heavy on memory because of the need to maintain pointers to other
    members, and it does not offer a random access iterator.
    */
    typedef std::list<rvector> BigDataCollection;

    /*! \brief Define type BigDataItr as an iterator BigData.
    */
    typedef BigDataCollection::iterator BigDataItr;

	/*! \brief Define type BigDataItr as a BigData const iterator.
    */
    typedef BigDataCollection::const_iterator BigDataConstItr;

    /*! \brief Define type NodeData as a container for iterators a BigData.

    Used by SPSnode to hold iterators to the data associated with node.
    */
    typedef std::vector<BigDataItr> NodeData;

    /*! \brief Define type NodeDataItr as a NodeData iterator.
    */
    typedef NodeData::iterator NodeDataItr;

	/*! \brief Define type NodeDataConstItr as a NodeData const iterator.
    */
    typedef NodeData::const_iterator NodeDataConstItr;

    /*! \brief Define type VecDotPrec as a container of cxsc dotprecision
    variables.

    Can be used to hold accumulators for rvector elements.  Used to hold
    representations of the sum of data points on each of the relevant
    dimensions.
    */
    typedef std::vector<dotprecision> VecDotPrec;


    /*! \brief Define type VecDotPrecIt as an iterator over VecDotPrec container.

    */
    typedef VecDotPrec::iterator VecDotPrecIt;


    //@}

} // end namespace subpavings

#endif

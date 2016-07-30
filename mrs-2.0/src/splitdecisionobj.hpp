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
\brief classes for determining whether to split an SPSnode.
*/

#ifndef ___SPSNODESPLIT_HPP__
#define ___SPSNODESPLIT_HPP__


namespace subpavings {

    //! Forward class declarations
    class SPSnode;


    /*! \brief A Virtual class providing decisions on whether to split spsnodes.

    Used when splitting occurs as data is inserted (not priority queue changes).
    */
    class SplitDecisionObj
    {
        public:

        virtual bool operator()() const = 0;

        virtual bool operator()(const SPSnode * const spn) const = 0;

    };



    /** \brief Classes derived from SplitDecisionObj

    These class provide functions provide a test for whether a node should be
    split when the histogram is being formed as each data point is inserted.
    The function will return true if the node should be split.

    These classes are used with one-at-a-time insertion of data points into
    a subpaving tree.
    */
    //@{

    /*! \brief Class for splitting based on number of points in a node.

    This class will try to split a node with > maxPoints of data points
    associated with it, but will <b>never</b> split any node if the result would
    give at least one child with < minChildPoints of data associated with it.
    i.e. the class will only split a node if the node has > maxPoints of data
    and both children after the split will have >= minChildPoints of data.

    Thus this class may leave a node with > maxPoints of data associated with
    it unsplit because splitting the node would create a child node
    with < minChildPoints of data associated with it.
    */
    class SplitOnK : public SplitDecisionObj
    {
        size_t maxPoints;
        size_t minChildPoints;

        SplitOnK(); // private default constructor

        public:

        explicit  SplitOnK(size_t k, size_t minCP = 0) : maxPoints(k),
                                        minChildPoints(minCP)  {};

        // function used to test if the object provides a real split test
        bool operator()() const
        {
            return true;
        }

        bool operator()(const SPSnode * const spn) const
        {
            bool retvalue = (spn->getCounter() > maxPoints);
            // check whether splitting would give children with enough points
            if (retvalue && (minChildPoints > 0)) {
				
				retvalue = ((spn->getCounter() >= minChildPoints) &&
                    ((spn->getMinChildCountIfSplit() == 0)
                    ||
                    (spn->getMinChildCountIfSplit() >= minChildPoints))
                );
            }
            return retvalue;
        }
    };

    /*! \brief Class for splitting based on average volume per point of a node.

    This class will try to split a node with average volume per point of data
    > avgVol, but will <b>never</b> split any node if the result would
    give at least one child with < minChildPoints of data associated with it.
    i.e. the class will only split a node if the node has > maxPoints of data
    and both children after the split will have >= minChildPoints of data.

     Thus this class may leave a node with average volume per point of data
    > avgVol unsplit because splitting the node would create a child node
    with < minChildPoints of data associated with it.

    */
    class SplitOnVolDivK : public SplitDecisionObj
    {
        double avgVol;
        size_t minChildPoints;

        SplitOnVolDivK(); // private default constructor

        public:

        explicit SplitOnVolDivK(double avg, size_t minCP = 0) : avgVol(avg),
                                                    minChildPoints(minCP) {};

        bool operator()() const
        {
            return true;
        }

        bool operator()(const SPSnode * const spn) const
        {
            bool retvalue = ((spn->nodeVolume())/(spn->getCounter()) > avgVol);
            // check whether splitting would give children with enough points
            if (retvalue && (minChildPoints > 0)) {
				
				retvalue = ((spn->getCounter() >= minChildPoints) &&
                    ((spn->getMinChildCountIfSplit() == 0)
                    ||
                    (spn->getMinChildCountIfSplit() >= minChildPoints))
                );
            }
            return retvalue;
        }
    };

    /*! \brief Class for never splitting.
    */
    class SplitNever : public SplitDecisionObj
    {
        public:

        bool operator()() const
        {
            return false; // this is just a dummy split test object
        }

        bool operator()(const SPSnode * const /*spn*/) const
        {
            return false;
        }
    };

    //@}
}


#endif



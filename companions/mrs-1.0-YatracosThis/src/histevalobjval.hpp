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

/*! \file      HistEvalObjVal.hpp
\brief Declarations for classes for evaluating when to stop changing histograms.
*/

#ifndef ___HISTEVALOBJVAL_HPP__
#define ___HISTEVALOBJVAL_HPP__

namespace subpavings {

//! Forward class declarations
class AdaptiveHistogramValidation;

/*! \brief A Virtual class providing a way to stop histogram changes.
*/
class HistEvalObjVal {

    public:

    /*! return true when splitting should stop. */
    virtual bool operator()(const AdaptiveHistogramValidation * const adh) const = 0;
};


/** @name Concrete classes derived from HistEvalObjVal

These classes provide a stopping test for priority queue splitting or merging.
The function operator () returns true when splitting or merging should stop.

*/
//@{


/*! \brief Class for testing the number of bins of a histogram.
*/
class CritLeaves_GTEV : public HistEvalObjVal
{
    size_t test;

    CritLeaves_GTEV(); // private default constructor

    public:

    CritLeaves_GTEV(size_t t) : test(t) {}

    /*! True if the number of leaves is >= test. */
    bool operator()
        (const AdaptiveHistogramValidation * const adh) const
    {
        return (spLeaves(adh->getSubPaving()) >= test);
    }
};


/*! \brief Class for testing the number of bins of a histogram.
*/
class CritLeaves_LTEV : public HistEvalObjVal
{
    size_t test;

    CritLeaves_LTEV(); // private default constructor

    public:

    CritLeaves_LTEV(size_t t) : test(t) {}

    /*! True if the number of leaves is <= test. */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    {
        return (spLeaves(adh->getSubPaving()) <= test);
    }
};

/*! \brief Class for testing the count of the node with the smallest
count in histogram's subpaving.
*/
class CritSmallestCount_LTEV : public HistEvalObjVal
{
    size_t test;

    CritSmallestCount_LTEV(); // private default constructor

    public:

    CritSmallestCount_LTEV(size_t t) : test(t) {}

    /*! True if the smallest leaf count is <= test. */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    {
        return (((adh->getSubPaving())->getSmallestLeafCount()) <= test);
    }

};

/*! \brief Class for testing the count of the node with the largest
count in histogram's subpaving.
*/
class CritLargestCount_LTEV : public HistEvalObjVal
{
    size_t test;

    CritLargestCount_LTEV(); // private default constructor

    public:

    CritLargestCount_LTEV(size_t t) : test(t) {}

    /*! True if the largest leaf count is <= test. */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    {
        return (((adh->getSubPaving())->getLargestLeafCount()) <= test);
    }
};

/*! \brief Class for testing the volume of the box with the smallest
volume in the histogram's subpaving.
*/
class CritSmallestVol_LTEV : public HistEvalObjVal
{
    double test;

    CritSmallestVol_LTEV(); // private default constructor

    public:

    CritSmallestVol_LTEV(double t) : test(t) {}

    /*! True if the volume of the smallest leaf (by vol) is <= test. */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    {
        return ((adh->getSubPaving())->getSmallestLeafVol() <= test);
    }
};

/*! \brief Class for testing the volume of the box with the largest
volume in the histogram's subpaving.
*/
class CritLargestVol_LTEV : public HistEvalObjVal
{
    double test;

    CritLargestVol_LTEV(); // private default constructor

    public:

    CritLargestVol_LTEV(double t) : test(t) {}

    /*! True if the volume of the largest leaf (by vol) is <= test. */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    {
        return ((adh->getSubPaving())->getLargestLeafVol() <= test);
    }
};

/*! \brief Class to bale out of priority queue splitting.
*/
class CritStopAllV: public HistEvalObjVal
{
    public:

    // use default constructor

    /*! True always */
    bool operator()(const AdaptiveHistogramValidation * const adh) const
    { return true; }
};

//@}

} // end of namespace subpavings

#endif



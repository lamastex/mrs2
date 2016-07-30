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
\brief Tools for voxels

*/

#ifndef ___VOXEL_TOOLS_HPP__
#define ___VOXEL_TOOLS_HPP__

#include "sptypes.hpp"


namespace subpavings {

    /*! \brief parse a .vtk header line for spacings data.

    Expects a line of the form "DIMENSIONS 32 32 32"
    \param line is the line to be parsed.
    \param spacings is an empty container for the spacings.
    \return the container of spacings filled up from the data in the line.
    */
    IntVec parseSpacings(string line, IntVec& spacings);

	/*! \brief Look for a black marker in a line from a .vtk file.

    \param line is the line to be searched.
    \param seek is the string to look for.
    \return true if found.
    */
	bool findBlack(string line, string seek);

    /*! \brief Get coordinates from a .vtk file.

    Expects structured point format data.
    @param Xs a container of integers to put x-coordinate data into.
    @param Ys a container of integers to put y-coordinate data into.
    @param Zs a container of integers to put z-coordinate data into.
    @param s the filename to get the data from.
    @pre empty containers for Xs, Ys, Zs and a properly formated file s.
    @post containers Xs, Ys, Zs filled with coordinate data.
    @return a vector of the coordinate spacing in the x, y, z direction.
    */
    IntVec getCoordinatesFromVtk(IntVec& Xs, IntVec& Ys, IntVec& Zs,
                                            const string& s);


} // end namespace subpavings

#endif

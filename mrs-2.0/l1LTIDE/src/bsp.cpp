/*
* Copyright (C) 2016 Dillon George
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
\brief
*/

#include "DensityTree/bsp.hpp"


void BSP_node::generate_hyperplane(std::vector<Point *> points) {
  // Assume points is a vector of length d, containing the points
  // defining the splitting hyperplane
  std::vector<Point> sample = sample_N_points(points, points[0]->dimension());
  Point orientation = sample[0];

  hyperplane = Plane::Hyperplane_d(sample.begin(), sample.end(), orientation, CGAL::ON_ORIENTED_BOUNDARY);
}


void BSP_node::print_tree(int max_depth, int indent_level) {
  if (  (hyperplane.dimension() == 2)) {
    std::cout << hyperplane << std::endl;
  }

  if (left != nullptr) {
    left->print_tree(max_depth, indent_level);
  }
  if (right != nullptr) {
    right->print_tree(max_depth, indent_level);
  }
}


int BSP_node::Max_Num() {
  // Returns the count of the leaf node containing the most number of points.
  if (left && right) {
    return std::max(left->Max_Num(), right->Max_Num());
  }
  if (left) {
    return left->Max_Num();
  }
  if (right) {
    return right->Max_Num();
  }

  return static_cast<int>(enclosed_points.size());
}


int BSP_node::Min_Num() {
  // Returns the count of the leaf node containing the most number of points.
  if (left && right) {
    return std::min(left->Min_Num(), right->Min_Num());
  }
  if (left) {
    return left->Min_Num();
  }
  if (right) {
    return right->Min_Num();
  }

  return static_cast<int>(enclosed_points.size());

}

std::tuple<BSP_node *, BSP_node *> BSP_node::split() {

  // Generate splitting hyperplane
  generate_hyperplane(enclosed_points);

  // Group points into left and right sides of plane
  left = new BSP_node();
  right = new BSP_node();

  // Reserve memory for points
  left->enclosed_points.reserve(num_points / 2 + 1);
  right->enclosed_points.reserve(num_points / 2 + 1);

  for (auto&& point: enclosed_points) {
    if (hyperplane.has_on_negative_side(*point)) {
      left->enclosed_points.push_back(point);
      left->num_points += 1;
    }
    else if (hyperplane.has_on_positive_side(*point)) {
      right->enclosed_points.push_back(point);
      right->num_points += 1;
    }
    else {
      // Point falls on hyper plane
      // TODO: Deal with this case properly
//      std::cerr << "Point lies on hyperplane. P: " << *point << std::endl;
    }
  }
  enclosed_points.clear();

  return std::make_tuple(left, right);
}

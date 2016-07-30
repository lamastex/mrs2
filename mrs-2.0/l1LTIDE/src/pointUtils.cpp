/*
* copyright (c) 2016 dillon george
*
* this file is part of mrs, a c++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the gnu general public license as published by
* the free software foundation; either version 3 of the license, or (at
* your option) any later version.
*
* this program is distributed in the hope that it will be useful, but
* without any warranty; without even the implied warranty of
* merchantability or fitness for a particular purpose.  see the gnu
* general public license for more details.
*
* you should have received a copy of the gnu general public license
* along with this program; if not, write to the free software
* foundation, inc., 675 mass ave, cambridge, ma 02139, usa.
*/

/*! \file
\brief
*/

#include "DensityTree/PointUtils.hpp"

std::vector<Point> randUnif(int dim, int num_points) {
  double size = 1.0;
  CGAL::Random_points_in_cube_d<Point> gen(dim, size);
  std::vector<Point> points;
  points.reserve(num_points);


  for (int i = 0; i < num_points; i++) {
    points.push_back(*gen++);
  }

  return points;

}


Point randNormPoint(unsigned int dim,
                    std::normal_distribution<double> &dist,
                    std::mt19937_64 &mt) {
  // Return a d-dimensional point with co-ordinates normally distributed
  std::vector<double> coords;
  coords.reserve(dim);
  for (int i = 0; i < dim; i++) {
    coords.push_back(dist(mt));
  }

  // Construct from iterator
  // TODO: Use custom iterator? If this is a bottleneck
  return Point(dim, coords.begin(), coords.end());
}


std::vector<Point> randNorm(unsigned int dim, size_t num_points, double mean, double stddev) {
  std::vector<Point> points;
  points.reserve(num_points);


  // Normally distributed rng
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::normal_distribution<double> dist(mean, stddev);

  for (int i = 0; i < num_points; ++i) {
    points.push_back(randNormPoint(dim, dist, mt));
  }

  return points;
}

std::vector<Point> randPoints(int dim, int num_points, GeneratorType gen) {

  std::vector<Point> points;

  switch (gen) {
    case GeneratorType::Uniform:
      points = randUnif(dim, num_points);
      break;
    case GeneratorType::Normal:
      points = randNorm(dim, num_points);
      break; // Default to mean = 0, stddev = 1
  }
  return points;

}


// Temporary create object that wraps this to avoid creating a new device on every
// call. Probably make a rng class
int uni_rng(int num_points) {
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0, num_points - 1);

  return uni(rng);
}


std::vector<Point> sample_N_points(std::vector<Point *> points, int n) {
  std::vector<Point> sample;
  sample.reserve(n);

  for (int i = 0; i < n; i++) {
    int index = uni_rng(points.size());
    sample.push_back(*points[index]);
  }

  return sample;

}


std::vector<Point> sample_N_points(std::vector<Point> points, int n) {
  std::vector<Point> sample;
  sample.reserve(n);

  for (int i = 0; i < n; i++) {
    int index = uni_rng(points.size());
    sample.push_back(points[index]);
  }


  return sample;
}

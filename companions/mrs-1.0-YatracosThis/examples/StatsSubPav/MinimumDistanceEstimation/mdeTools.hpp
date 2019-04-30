/*! \file
\brief Common routines for MDE

 */

#ifndef __MDE_TOOLS_HPP__
#define __TEST_TOOLS_HPP__

#include <vector>
#include <cstddef>

// return a vector of the top k indices of ...
void topk(std::vector<double> a, std::vector<int> & indtop, size_t k);

#endif


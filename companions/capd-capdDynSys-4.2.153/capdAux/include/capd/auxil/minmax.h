/// @addtogroup capd
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file minmax.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

/* min, max and abs definitions */

#ifndef _CAPD_CAPD_MINMAX_H_
#define _CAPD_CAPD_MINMAX_H_

#undef max
#undef min

namespace capd{

//
// The following lines was prone to errors
//
//template<typename T>
//inline T min(const T& x, const T& y) {
//  return (x<y ? x : y);
//}
//
//template<typename T>
//inline T max(const T& x, const T& y) {
//  return (x<y ? y : x);
//}


template<typename T>
inline T min(const T& /*x*/, const T & /*y*/) {
  return T::Specialization_of_min_function_not_defined();
}

template<>
inline long double min(const long double & x, const long double & y) {
  return (x<y ? x : y);
}

template<>
inline double min(const double & x, const double & y) {
  return (x<y ? x : y);
}

template<>
inline int min(const int & x, const int & y) {
  return (x<y ? x : y);
}



template<typename T>
inline T max(const T& /*x*/, const T & /*y*/) {
  return T::Specialization_of_max_function_not_defined();
}

template<>
inline long double max(const long double & x, const long double & y) {
  return (x<y ? y : x);
}

template<>
inline double max(const double & x, const double & y) {
  return (x<y ? y : x);
}

template<>
inline int max(const int & x, const int & y) {
  return (x<y ? y : x);
}

template<typename T>
inline T abs(const T& /*x*/) {
  return T::Specialization_of_abs_function_not_defined();
}

template<>
inline long double abs(const long double & x) {
  return (x<0.) ? -x : x;
}

template<>
inline double abs(const double & x) {
  return (x<0.) ? -x : x;
}

template<>
inline int abs(const int & x) {
  return (x<0.) ? -x : x;
}

} // end of namespace capd



#endif // _CAPD_CAPD_MINMAX_H_

/// @}

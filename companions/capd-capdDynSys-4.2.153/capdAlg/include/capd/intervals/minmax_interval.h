/////////////////////////////////////////////////////////////////////////////
/// @file minmax_interval
///
/// @author Mateusz Juda <mateusz.juda@{ii.uj.edu.pl,gmail.com}>
///
/// @date 2014-07-23
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H
#define CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H

#include <capd/auxil/minmax.h>

namespace capd{
namespace intervals{
template < typename T_Bound, typename T_Rnd >
class Interval;
}  //namespace intervals

   // definitions in interval.h and interval.cpp
   template < typename T_Bound, typename T_Rnd>
   intervals::Interval< T_Bound, T_Rnd > min(const intervals::Interval< T_Bound, T_Rnd >& x, const intervals::Interval< T_Bound, T_Rnd >& y);
   template < typename T_Bound, typename T_Rnd>
   intervals::Interval< T_Bound, T_Rnd > max(const intervals::Interval< T_Bound, T_Rnd >& x, const intervals::Interval< T_Bound, T_Rnd >& y);
   template < typename T_Bound, typename T_Rnd>
   intervals::Interval< T_Bound, T_Rnd >  abs(const intervals::Interval< T_Bound, T_Rnd >& x);

} // namespace capd

//#ifdef  __USE_CXSC__

namespace capd{
namespace cxsc{
class Interval;
} // end of namespace cxsc

// definitions in capdAlg/capd/filib/Interval.h
template<>
inline cxsc::Interval abs (const cxsc::Interval & A_inter);
template<>
inline cxsc::Interval max(const cxsc::Interval& A_iv1, const cxsc::Interval& A_iv2);
template<>
inline cxsc::Interval min (const cxsc::Interval& A_iv1, const cxsc::Interval& A_iv2);

} // namespace capd
//#endif // __USE_FILIB__


#ifdef  __USE_FILIB__

#include <interval/interval.hpp>

namespace capd{
namespace filib{

template <typename T, ::filib::rounding_strategy R, ::filib::interval_mode M>
class Interval;

} // end of namespace filib

// definitions in capdAlg/capd/filib/Interval.h
template <typename T, ::filib::rounding_strategy R, ::filib::interval_mode M>
filib::Interval<T, R, M> abs (const filib::Interval<T, R, M> & A_inter);
template <typename T, ::filib::rounding_strategy R, ::filib::interval_mode M>
inline filib::Interval<T, R, M> max(const filib::Interval<T, R, M>& A_iv1, const filib::Interval<T, R, M>& A_iv2);
template <typename T, ::filib::rounding_strategy R, ::filib::interval_mode M>
inline filib::Interval<T, R, M> min (const filib::Interval<T, R, M>& A_iv1, const filib::Interval<T, R, M>& A_iv2);

} // namespace capd
#endif // __USE_FILIB__

#endif // CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H

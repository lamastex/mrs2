//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file IComplex.h
///
/// @author Tomasz Kapela   @date 2010-03-15
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_INTERVAL_ICOMPLEX_H_
#define _CAPD_INTERVAL_ICOMPLEX_H_
#include <complex>
#include <stdexcept>
#include "capd/basicalg/power.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/basicalg/doubleFun.h"
#include "capd/intervals/Interval.h"

namespace capd{
namespace intervals{

// TODO: it should be implemented in the CAPD
//template <typename T>
//class IComplex : public std::complex<Interval<T> >{
//}

template< typename T>
bool intersection(const std::complex<T> & a, const std::complex<T> & b, std::complex<T> & result ){
  T re, im;
  bool exist = intersection(real(a),real(b),re) &&  intersection(imag(a),imag(b), im);
  if(exist)
    result = std::complex<T>(re, im);
  return exist;
}

template< typename T>
const std::complex<T> intervalHull(const std::complex<T> & a, const std::complex<T> & b ){
  return std::complex<T>(intervalHull(a.real(),b.real()),intervalHull(a.imag(),b.imag()));
}

} // end of intervals namespace

template <typename T>
class TypeTraits<std::complex<T> > {
public:
  typedef typename capd::TypeTraits<T>::Real Real;
  static inline std::complex<T> zero() {
    return std::complex<T>(static_cast<T> (0.0), static_cast<T> (0.0));
  }
  static inline std::complex<T> one() {
    return std::complex<T>(static_cast<T> (1.0), static_cast<T> (0.0));
  }
  static inline int numberOfDigits(){
    return TypeTraits<T>::numberOfDigits();
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline T epsilon() throw(){
    return TypeTraits<T>::epsilon();
  }
  template <typename S>
  static inline std::complex<T> convert(const S & obj){
    return static_cast<std::complex<T> >(TypeTraits<T>::convert(obj));
  }
  static inline std::complex<T> convert(const T & obj){
    return static_cast<std::complex<T> >(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = TypeTraits<T>::isInterval;

};

template <typename T>
T abs(const T & x);

template <typename T>
T abs(const std::complex<T> & x){
   return  sqrt(sqr(capd::abs(x.real()))+sqr(capd::abs(x.imag())));
}

namespace vectalg {
template<typename T>
// Computes euclidean norm of any vector
T euclNorm(std::complex<T> & x){
  typedef T Scalar;
  Scalar sum = TypeTraits<Scalar>::zero();
  typename std::complex<T>::const_iterator b=x.begin(), e=x.end();
  while(b!=e)
  {
    sum += norm(*b);
    ++b;
  }
  return Scalar(sqrt(nonnegativePart(sum)));
}

}
} // end of namespace capd

namespace std{
// sometimes it is needed for some T
/// TODO: should be corrected
template <typename T>
inline std::complex<T> power(const std::complex<T> & x, int n){
  if(n<0)
    throw std::runtime_error("power implemented only for n>=0");
  std::complex<T> res = capd::TypeTraits<std::complex<T> >::one();
  for(int i=0; i<n; i++){
    res *= x;
  }
  return res;
}


template <typename T>
bool isSingular(const std::complex<T> & x) {
  return ((isSingular(x.real())) && (isSingular(x.imag())));
}

template <typename T>
std::complex<T> mid(const std::complex<T> & x){
   return  std::complex<T>(mid(x.real()),mid(x.imag()));
}

template< typename T>
bool intersection(const std::complex<T> & a, const std::complex<T> & b, std::complex<T> & result ){
  T re, im;
  bool exist = intersection(real(a),real(b),re) &&  intersection(imag(a),imag(b), im);
  if(exist)
    result = std::complex<T>(re, im);
  return exist;
}
}


#endif // _CAPD_INTERVAL_ICOMPLEX_H_

// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <complex>

namespace himalaya {

// some templates to perform operations between int's and complex<double>
template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator* ( const std::complex<T>& c, SCALAR n ) noexcept { return c * T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator* ( SCALAR n, const std::complex<T>& c ) noexcept { return T(n) * c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator/ ( const std::complex<T>& c, SCALAR n ) noexcept { return c / T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator/ ( SCALAR n, const std::complex<T>& c ) noexcept { return T(n) / c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator+ ( const std::complex<T>& c, SCALAR n ) noexcept { return c + T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator+ ( SCALAR n, const std::complex<T>& c ) noexcept { return T(n) + c ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator- ( const std::complex<T>& c, SCALAR n ) noexcept { return c - T(n) ; }

template< typename T, typename SCALAR > inline
typename std::enable_if< !std::is_same<T,SCALAR>::value, std::complex<T> >::type
operator- ( SCALAR n, const std::complex<T>& c ) noexcept { return T(n) - c ; }

template <typename T> T pow2(T x)    noexcept { return x*x; }
template <typename T> T pow3(T x)    noexcept { return x*x*x; }
template <typename T> T pow4(T x)    noexcept { return x*x*x*x; }
template <typename T> T pow5(T x)    noexcept { return x*x*x*x*x; }
template <typename T> T pow6(T x)    noexcept { return x*x*x*x*x*x; }
template <typename T> T pow7(T x)    noexcept { return x*x*x*x*x*x*x; }
template <typename T> T pow8(T x)    noexcept { return x*x*x*x*x*x*x*x; }
template <typename T> T pow9(T x)    noexcept { return x*x*x*x*x*x*x*x*x; }
template <typename T> T power10(T x) noexcept { return x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow11(T x)   noexcept { return x*x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow12(T x)   noexcept { return x*x*x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow13(T x)   noexcept { return x*x*x*x*x*x*x*x*x*x*x*x*x; }
template <typename T> T pow14(T x)   noexcept { return x*x*x*x*x*x*x*x*x*x*x*x*x*x; }

} // namespace himalaya

// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/mh2_fo/pv.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

/**
 * @file pv.cpp
 *
 * @brief Implementation of real Passarino-Veltman loop functions with
 * squared arguments.
 */

namespace himalaya {
namespace mh2_fo {

namespace {

constexpr double EPSTOL = 1.0e-11; ///< underflow accuracy

constexpr double dabs(double a) noexcept { return a >= 0. ? a : -a; }
constexpr double sqr(double a) noexcept { return a*a; }
constexpr double pow3(double a) noexcept { return a*a*a; }
double log_abs(double a) noexcept { return std::log(std::abs(a)); }


/// fast implementation of complex logarithm
template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   const T rz = std::real(z);
   const T iz = std::imag(z);

   return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
}


constexpr bool is_close(double m1, double m2, double tol) noexcept
{
   const double mmax = std::max(dabs(m1), dabs(m2));
   const double mmin = std::min(dabs(m1), dabs(m2));
   const double max_tol = tol * mmax;

   if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
      return mmax - mmin <= tol;

   return mmax - mmin <= max_tol;
}


/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
T divide_finite(T a, T b) noexcept
{
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}


double fB(const std::complex<double>& a) noexcept
{
   const double re = std::real(a);
   const double im = std::imag(a);

   if ((std::abs(re) == 0.0 || std::abs(re) == 1.0) && im == 0.0) {
      return -1.0;
   }

   return std::real(-1.0 + fast_log(1.0 - a) - a*fast_log(1.0 - 1.0/a));
}


} // anonymous namespace


double a0(double m2, double q2) noexcept
{
   if (std::abs(m2) < 1e-8)
      return 0.;

   return m2 * (1.0 - log_abs(m2 / q2));
}


/**
 * B0 function with squared arguments, from hep-ph/9606211.
 * Note it returns the REAL PART ONLY.
 */
double b0(double p2, double m12, double m22, double q2) noexcept
{
   // protect against infrared divergence
   if (is_close(p2, 0., EPSTOL) && is_close(m12, 0., EPSTOL)
       && is_close(m22, 0., EPSTOL))
      return 0.;

   double ans = 0.;
   const double mMax2 = std::max(std::abs(m12), std::abs(m22));
   const double pTest = std::abs(divide_finite(p2, mMax2));

   if (pTest > 1e-10) {
      const double s = p2 - m22 + m12;
      const std::complex<double> iEpsilon(0., EPSTOL * mMax2);
      const std::complex<double> rt = std::sqrt(sqr(s) - 4.*p2*(m12 - iEpsilon));
      const std::complex<double> xPlus = 0.5 * (s + rt) / p2;
      const std::complex<double> xMinus = 0.5 * (s - rt) / p2;

      ans = -log_abs(p2 / q2) - fB(xPlus) - fB(xMinus);
   } else {
      if (is_close(m12, m22, EPSTOL)) {
         ans = -log_abs(m12 / q2);
      } else {
         const double mMin2 = std::min(std::abs(m12), std::abs(m22));

         if (mMin2 < 1.e-30) {
            ans = 1. - log_abs(mMax2 / q2);
         } else {
            ans = (m12*(1. - log_abs(m12/q2)) - m22*(1. - log_abs(m22/q2)))
               / (m12 - m22);
         }
      }
   }

   return ans;
}


/**
 * Derivative of B0(p^2, m1^2, m2^2, Q^2) w.r.t. p^2, for p^2 = 0.
 *
 * @note Implemented only in the p^2 = 0 limit.
 *
 * @param m12 squared mass
 * @param m22 squared mass
 *
 * @return derivative of B0 w.r.t. p^2 at p^2 = 0
 */
double d1_b0(double m12, double m22) noexcept
{
   const double m14 = m12 * m12;
   const double m24 = m22 * m22;

   if ((std::abs(m12) < 0.0001) != (std::abs(m22) < 0.0001)) {
      return (m14 - m24) / (2. * pow3(m12 - m22));
   } else if (std::abs(m12) < 0.0001 && std::abs(m22) < 0.0001) {
      return 0.;
   } else if (std::abs(m22 - m12) < 0.001) {
      return 1./(6. * m12) + (m12 - m22)/(12.* m14);
   }

   return (m14 - m24 + 2. * m12 * m22 * std::log(m22/m12))
      /(2. * pow3(m12 - m22));
}

} // namespace mh2_fo
} // namespace himalaya

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


constexpr double sign(double x) noexcept
{
   return x >= 0.0 ? 1.0 : -1.0;
}


/// fast implementation of complex logarithm
template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   const T rz = std::real(z);
   const T iz = std::imag(z);

   return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
}


/// compares a number for being close to zero
constexpr bool is_zero(double a, double prec) noexcept
{
   return dabs(a) <= prec;
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


double fB(const std::complex<double>& x) noexcept
{
   const double re = std::real(x);
   const double im = std::imag(x);

   if ((std::abs(re) == 0.0 || std::abs(re) == 1.0) && im == 0.0) {
      return -1.0;
   }

   return std::real(-1.0 + fast_log(1.0 - x) - x*fast_log(1.0 - 1.0/x));
}


/// fB(xp) + fB(xm)
double fB(const std::complex<double>& xp, const std::complex<double>& xm) noexcept
{
   const double rep = std::real(xp);
   const double imp = std::imag(xp);

   if ((std::abs(rep) == 0.0 || std::abs(rep) == 1.0) && imp == 0.0) {
      return -1.0 + fB(xm);
   }

   const double rem = std::real(xm);
   const double imm = std::imag(xm);

   if ((std::abs(rem) == 0.0 || std::abs(rem) == 1.0) && imm == 0.0) {
      return -1.0 + fB(xp);
   }

   return std::real(-2.0 + fast_log((1.0 - xp)*(1.0 - xm))
      - xp*fast_log(1.0 - 1.0/xp) - xm*fast_log(1.0 - 1.0/xm));
}


} // anonymous namespace


double a0(double m2, double q2) noexcept
{
   m2 = std::abs(m2);

   if (m2 < 1e-8)
      return 0.;

   return m2 * (1 - std::log(m2 / q2));
}


/**
 * B0 function with squared arguments, from hep-ph/9606211.
 * @note returns only the real part of B0
 */
double b0(double p2, double m12, double m22, double q2) noexcept
{
   // protect against infrared divergence
   if (is_zero(p2, EPSTOL) && is_zero(m12, EPSTOL) && is_zero(m22, EPSTOL)) {
      return 0.;
   }

   p2  = std::abs(p2);
   m12 = std::abs(m12);
   m22 = std::abs(m22);
   q2  = std::abs(q2);

   if (m12 > m22) {
      std::swap(m12, m22);
   }

   // p2 is no 0
   if (p2 > 1e-10*m22) {
      const double s = p2 - m22 + m12;
      const std::complex<double> imin(m12, -EPSTOL);
      const std::complex<double> x = std::sqrt(sqr(s) - 4.0 * p2 * imin);
      const std::complex<double> xp = (s + sign(s)*x) / (2*p2);
      const std::complex<double> xm = imin / (xp*p2);

      return -std::log(p2 / q2) - fB(xp, xm);
   }

   if (is_close(m12, m22, EPSTOL)) {
      return -std::log(m12 / q2);
   }

   if (m12 < 1e-30) {
      return 1 - std::log(m22 / q2);
   }

   return 1.0 - std::log(m22/q2)
        + m12 * std::log(m22/m12) / (m12 - m22);
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
   m12 = std::abs(m12);
   m22 = std::abs(m22);
   const double m14 = m12 * m12;
   const double m24 = m22 * m22;

   if ((m12 < 0.0001) != (m22 < 0.0001)) {
      return (m14 - m24) / (2. * pow3(m12 - m22));
   } else if (m12 < 0.0001 && m22 < 0.0001) {
      return 0.;
   } else if (std::abs(m22 - m12) < 0.001) {
      return 1./(6. * m12) + (m12 - m22)/(12.* m14);
   }

   return (m14 - m24 + 2 * m12 * m22 * std::log(m22/m12))
      /(2 * pow3(m12 - m22));
}

} // namespace mh2_fo
} // namespace himalaya

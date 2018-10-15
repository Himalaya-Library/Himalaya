// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "pv.hpp"

#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>

namespace himalaya {
namespace mh1l {

namespace {

constexpr double EPSTOL = 1.0e-11; ///< underflow accuracy
constexpr double TOL = 1e-4;

constexpr double sqr(double a) noexcept { return a*a; }
constexpr double pow3(double a) noexcept { return a*a*a; }
constexpr double pow6(double a) noexcept { return a*a*a*a*a*a; }
constexpr double log_abs(double a) noexcept { return std::log(std::abs(a)); }

template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   return std::complex<T>(std::log(std::abs(z)), std::arg(z));
}

bool is_close(double m1, double m2, double tol) noexcept
{
   using std::abs;
   const double mmax = abs(std::max(abs(m1), abs(m2)));
   const double max_tol = tol * mmax;

   if (max_tol == 0. && mmax != 0. && tol != 0.)
      return abs(m1 - m2) <= tol;

   return abs(m1 - m2) <= max_tol;
}

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
T divide_finite(T a, T b) noexcept {
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

double fB(const std::complex<double>& a) noexcept
{
   using std::abs;
   const double x = a.real();

   if (abs(x) < EPSTOL)
      return -1. - x + sqr(x) * 0.5;

   if (is_close(x, 1., EPSTOL))
      return -1.;

   return std::real(fast_log(1. - a) - 1. - a * fast_log(1.0 - 1.0 / a));
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
   using std::abs;

   // protect against infrared divergence
   if (is_close(p2, 0., EPSTOL) && is_close(m12, 0., EPSTOL)
       && is_close(m22, 0., EPSTOL))
      return 0.;

   double ans = 0.;
   const double mMax2 = std::max(abs(m12), abs(m22));
   const double pTest = abs(divide_finite(p2, mMax2));

   if (pTest > 1e-10) {
      const double s = p2 - m22 + m12;
      const std::complex<double> iEpsilon(0., EPSTOL * mMax2);
      const std::complex<double> rt = sqrt(sqr(s) - 4.*p2*(m12 - iEpsilon));
      const std::complex<double> xPlus = 0.5 * (s + rt) / p2;
      const std::complex<double> xMinus = 0.5 * (s - rt) / p2;

      ans = -log_abs(p2 / q2) - fB(xPlus) - fB(xMinus);
   } else {
      if (is_close(m12, m22, EPSTOL)) {
         ans = -log_abs(m12 / q2);
      } else {
         const double mMin2 = std::min(abs(m12), abs(m22));

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
 * @param m2a squared mass
 * @param m2b squared mass
 *
 * @return derivative of B0 w.r.t. p^2 at p^2 = 0
 */
double d1_b0(double m2a, double m2b) noexcept
{
   using std::abs;

   const double m4a = m2a * m2a;
   const double m4b = m2b * m2b;

   if ((std::abs(m2a) < 0.0001) != (std::abs(m2b) < 0.0001)) {
      return (m4a - m4b) / (2. * pow3(m2a - m2b));
   } else if (std::abs(m2a) < 0.0001 && std::abs(m2b) < 0.0001) {
      return 0.;
   } else if (std::abs(m2b - m2a) < 0.001) {
      return 1./(6. * m2a) + (m2a - m2b)/(12.* m4a);
   }

   return (m4a - m4b + 2. * m2a * m2b * std::log(m2b/m2a))
      /(2. * pow3(m2a - m2b));
}

} // namespace mh1l
} // namespace himalaya

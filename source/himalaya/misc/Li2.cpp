// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/misc/Li2.hpp"
#include "himalaya/misc/complex.hpp"
#include <cmath>
#include <limits>

/**
 * @file Li2.cpp
 * @brief Implementation of the dilogarithm function
 * @note The implementation has been taken from the polylogarithm package.
 */

namespace himalaya {

namespace {

   template <int Nstart, typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(0 <= Nstart && Nstart < N && N >= 2, "invalid array bounds");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= Nstart; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a
 * maximum error of 4.16e-18.
 */
double dilog(double x) noexcept
{
   const double PI = 3.1415926535897932;
   const double P[] = {
      1.0706105563309304277e+0,
     -4.5353562730201404017e+0,
      7.4819657596286408905e+0,
     -6.0516124315132409155e+0,
      2.4733515209909815443e+0,
     -4.6937565143754629578e-1,
      3.1608910440687221695e-2,
     -2.4630612614645039828e-4
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -4.5355682121856044935e+0,
      8.1790029773247428573e+0,
     -7.4634190853767468810e+0,
      3.6245392503925290187e+0,
     -8.9936784740041174897e-1,
      9.8554565816757007266e-2,
     -3.2116618742475189569e-3
   };

   double y = 0, r = 0, s = 1;

   // transform to [0, 1/2)
   if (x < -1) {
      const double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - std::log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5*l);
      s = 1;
   } else {
      const double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double z = y - 0.25;
   const double z2 = z*z;
   const double z4 = z2*z2;

   const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
                    z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7]));
   const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
                    z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7]));

   return r + s*y*p/q;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
std::complex<double> dilog(const std::complex<double>& z_) noexcept
{
   const double PI = 3.1415926535897932;
   const Complex<double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.0647616451442255e-11,
      + 8.9216910204564526e-13,
      - 1.9939295860721076e-14,
      + 4.5189800296199182e-16
   };

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return dilog(z.re);
      }
      // z.re > 1
      return { dilog(z.re), -PI*std::log(z.re) };
   }

   const double nz = norm_sqr(z);

   if (nz < std::numeric_limits<double>::epsilon()) {
      return z*(1.0 + 0.25*z);
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);
   double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5) {
      if (nz > 1) {
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      } else { // nz <= 1
         u = -log(1.0 - z);
         rest = 0;
         sgn = 1;
      }
   } else { // z.re > 0.5
      if (nz <= 2*z.re) {
         u = -log(z);
         rest = u*log(1.0 - z) + PI*PI/6;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const Complex<double> u2(u*u);

   return sgn*(u + u2*(bf[0] + u*horner<1>(u2, bf))) + rest;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
double clausen_2(double x) noexcept
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon()) {
      return 0.0;
   }

   return std::imag(dilog(std::exp(i*x)));
}

} // namespace himalaya

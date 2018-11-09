// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "DSZHiggs.hpp"
#include "DSZHiggs.h"
#include "dilog.hpp"
#include <cmath>
#include <limits>
#include <mutex>
#include <utility>

/**
 * @file DSZHiggs.cpp
 * @brief Implementation of C++ wrappers for 2-loop FORTRAN routines.
 */

namespace himalaya {
namespace mssm_twoloophiggs {

namespace {

/// locks MSSM fortran functions that are not threadsafe
std::mutex mtx;

template <typename T> T constexpr sqr(T a) { return a * a; }
template <typename T> T sqrtabs(T a) { return std::sqrt(std::abs(a)); }
template <typename T> T logabs(T x) { return std::log(std::abs(x)); }

/// limit st -> 0 and mst1 -> mst2
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_as_st_0_mst1_eq_mst2(
   double mt2, double mg, double mst12, double /* mst22 */,
   double /* sxt */, double /* cxt */, double scale2, double mu,
   double tanb, double vev2, double gs, int /* scheme */)
{
   using std::fabs;
   using std::sqrt;
   using std::atan;
   using std::log;
   using std::sin;
   using std::pow;
   using himalaya::dilog;

   constexpr double Pi2 = M_PI * M_PI;
   constexpr double Pi4 = M_PI * M_PI * M_PI * M_PI;
   const double g = sqr(mg);
   const double g2 = sqr(g);
   const double q = scale2;
   const double t = mt2;
   const double tsq = sqr(t);
   const double T = mst12;
   const double Tsq = sqr(mst12);
   const double del = g2 + tsq + Tsq - 2*(g*t + g*T + t*T);
   const double rdel = sqrtabs(del);
   const double sb = sin(atan(tanb));
   const double ht2 = 2./vev2*mt2/sqr(sb);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = 0.;

   result(0,1) = (sqr(gs)*ht2*mg*mt2*mu* (-1 + logabs(t/q) -
      (T*((2*del* (-(logabs(t/g)/T) - (2*(g + t - T +
      g*sqrtabs(del/g2))* logabs((g + t - T - g*sqrtabs(del/g2))/
      (2.*g)))/ (g*sqrtabs(del/g2)* (t - T + g*(-1 +
      sqrtabs(del/g2)))) + (2*logabs((g - t + T - g*sqrtabs(del/g2))/
      (2.*g)))/(g*sqrtabs(del/g2)) - (2*((g + t - T +
      g*sqrtabs(del/g2))* logabs((g + t - T + g*sqrtabs(del/g2))/
      (2.*g)) + (g - t + T - g*sqrtabs(del/g2))* logabs((g - t + T +
      g*sqrtabs(del/g2))/ (2.*g))))/ (g*sqrtabs(del/g2)* (t - T +
      g*(-1 + sqrtabs(del/g2))))))/ g2 + (2*(g + t - T)*(Pi2 -
      3*logabs(t/g)*logabs(T/g) + 6*logabs((g + t - T -
      g*sqrtabs(del/g2))/(2.*g))*logabs((g - t + T -
      g*sqrtabs(del/g2))/(2.*g)) - 6*dilog((g + t - T -
      g*sqrtabs(del/g2))/(2.*g)) - 6*dilog((g - t + T -
      g*sqrtabs(del/g2))/(2.*g))))/ (3.*g2)))/(2.*pow(fabs(del)/g2,1.5))))/
      (8.*Pi4*T);

   result(1,0) = result(0,1);

   result(1,1) = (sqr(gs)*ht2*mt2*(-2 - (8*(g + t)*(-1 +
      logabs(g/q)))/T + 8*logabs(t/g) + 6*(-1 + logabs(t/q)) +
      8*sqr(logabs(t/q)) - 4*logabs(T/g) + (4*((g2*T - sqr(t - T)*(2*t
      + T) + 2*g*t*(t + 5*T))*logabs(t/g) +
      4*g2*T*logabs(T/g)))/(del*T) + 2*logabs(T/q) -
      8*sqr(logabs(T/q)) + 5*logabs(Tsq/tsq) + sqr(logabs(Tsq/tsq)) +
      (4*g2*(g + t - T)*(Pi2 - 6*dilog((g + t - T -
      g*sqrtabs(del/g2))/ (2.*g)) - 6*dilog((g - t + T -
      g*sqrtabs(del/g2))/(2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((-rdel + g + t - T)/(2.*g))* logabs((-rdel + g - t +
      T)/(2.*g))))/(3.*pow(fabs(del),1.5)) + (4*g*(g + t - T)*(Pi2 -
      6*dilog((g + t - T - g*sqrtabs(del/g2))/(2.*g)) - 6*dilog((g - t
      + T - g*sqrtabs(del/g2))/ (2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((g + t - T - g*sqrtabs(del/g2))/(2.*g))* logabs((g - t
      + T - g*sqrtabs(del/g2))/(2.*g))))/ (3.*del*sqrtabs(del/g2)) +
      (8*mg*mu*(1/T - logabs(t/q)/T + (g*(g + t - T)* (Pi2 -
      6*dilog((g + t - T - g*sqrtabs(del/g2))/(2.*g)) - 6*dilog((g - t
      + T - g*sqrtabs(del/g2))/ (2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((-rdel + g + t - T)/(2.*g))*logabs((-rdel + g - t +
      T)/(2.*g))))/ (3.*pow(fabs(del),1.5)) + (g*(-(logabs(t/g)/T) -
      (2*(-(g*logabs(4.)) + (rdel + g + t - T)*logabs((-rdel + g + t -
      T)/g) + (-rdel + g - t + T)*logabs((-rdel + g - t + T)/g)))/
      (rdel*(rdel - g + t - T)) - 2*(((rdel + g + t - T)* logabs((g +
      t - T + g*sqrtabs(del/g2))/ (2.*g)))/ (rdel*(t - T + g*(-1 +
      sqrtabs(del/g2)))) - ((rdel - g - t + T)* logabs((g - t + T +
      g*sqrtabs(del/g2))/ (2.*g)))/ (rdel*(-t + T + g*(-1 +
      sqrtabs(del/g2)))))))/ rdel))/tanb))/(32.*Pi4);

   return result;
}

/// Pietro Slavich implementation
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_as_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
             &tanb, &vev2, &gs, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

} // anonymous namespace

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_as(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   if (std::abs((mst12 - mst22)/mst12) < 1e-8)
      return delta_mh2_2loop_at_as_st_0_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);

   return delta_mh2_2loop_at_as_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_at(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 2> result;

   {
      std::lock_guard<std::mutex> lg(mtx);

      ddshiggs_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
                &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
                &result(0,0), &result(0,1), &result(1,1));
   }

   result(1,0) = result(0,1);

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_as(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result(delta_mh2_2loop_at_as(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs, scheme));

   std::swap(result(0,0), result(1,1));

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_atau_atau(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   tausqhiggs_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
               &costau, &scale2, &mu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_atau(
   double mtau2, double mb2,
   double mstau12, double mstau22, double msb12, double msb22,
   double sintau, double costau, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 2> result;

   taubot_(&mtau2, &mb2,
           &mstau12, &mstau22, &msb12, &msb22,
           &sintau, &costau, &sxb, &cxb,
           &scale2, &mu, &tanb, &vev2,
           &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

} // namespace mssm_twoloophiggs
} // namespace himalaya

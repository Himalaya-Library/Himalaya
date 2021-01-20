// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "./DSZHiggs.hpp"
#include "./DSZHiggs.h"
#include "./Li2f.hpp"
#include "misc/CouplingOrders.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <mutex>
#include <utility>

/**
 * @file DSZHiggs.cpp
 * @brief Implementation of C++ wrappers for 2-loop FORTRAN routines.
 */

namespace himalaya {
namespace mh2l {

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
   double tanb, double vev2, double g3, int /* scheme */)
{
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
   const double sb = std::sin(std::atan(tanb));
   const double ht2 = 2./vev2*mt2/sqr(sb);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = 0.;

   result(0,1) = (sqr(g3)*ht2*mg*mt2*mu* (-1 + logabs(t/q) -
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
      g*sqrtabs(del/g2))/(2.*g)) - 6*li2((g + t - T -
      g*sqrtabs(del/g2))/(2.*g)) - 6*li2((g - t + T -
      g*sqrtabs(del/g2))/(2.*g))))/ (3.*g2)))/(2.*std::pow(std::abs(del)/g2,1.5))))/
      (8.*Pi4*T);

   result(1,0) = result(0,1);

   result(1,1) = (sqr(g3)*ht2*mt2*(-2 - (8*(g + t)*(-1 +
      logabs(g/q)))/T + 8*logabs(t/g) + 6*(-1 + logabs(t/q)) +
      8*sqr(logabs(t/q)) - 4*logabs(T/g) + (4*((g2*T - sqr(t - T)*(2*t
      + T) + 2*g*t*(t + 5*T))*logabs(t/g) +
      4*g2*T*logabs(T/g)))/(del*T) + 2*logabs(T/q) -
      8*sqr(logabs(T/q)) + 5*logabs(Tsq/tsq) + sqr(logabs(Tsq/tsq)) +
      (4*g2*(g + t - T)*(Pi2 - 6*li2((g + t - T -
      g*sqrtabs(del/g2))/ (2.*g)) - 6*li2((g - t + T -
      g*sqrtabs(del/g2))/(2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((-rdel + g + t - T)/(2.*g))* logabs((-rdel + g - t +
      T)/(2.*g))))/(3.*std::pow(std::abs(del),1.5)) + (4*g*(g + t - T)*(Pi2 -
      6*li2((g + t - T - g*sqrtabs(del/g2))/(2.*g)) - 6*li2((g - t
      + T - g*sqrtabs(del/g2))/ (2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((g + t - T - g*sqrtabs(del/g2))/(2.*g))* logabs((g - t
      + T - g*sqrtabs(del/g2))/(2.*g))))/ (3.*del*sqrtabs(del/g2)) +
      (8*mg*mu*(1/T - logabs(t/q)/T + (g*(g + t - T)* (Pi2 -
      6*li2((g + t - T - g*sqrtabs(del/g2))/(2.*g)) - 6*li2((g - t
      + T - g*sqrtabs(del/g2))/ (2.*g)) - 3*logabs(t/g)*logabs(T/g) +
      6*logabs((-rdel + g + t - T)/(2.*g))*logabs((-rdel + g - t +
      T)/(2.*g))))/ (3.*std::pow(std::abs(del),1.5)) + (g*(-(logabs(t/g)/T) -
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
   double tanb, double vev2, double g3, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
             &tanb, &vev2, &g3, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

Eigen::Matrix<double, 2, 2> rotate_by(double dm2, double tanb)
{
   const double tanb2 = sqr(tanb);
   const double sinb = tanb / sqrtabs(1. + tanb2);
   const double cosb = 1. / sqrtabs(1. + tanb2);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = dm2 * sqr(sinb);
   result(0,1) = - dm2 * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = dm2 * sqr(cosb);

   return result;
}

double calc_At(double mt2, double mst12, double mst22,
   double sxt, double cxt, double mu, double tanb)
{
   const double s2t = 2*cxt*sxt;
   const double Xt = (mst12 - mst22)*s2t/2./sqrtabs(mt2);
   const double At = Xt - mu/tanb;

   return At;
}

double phi(double x, double y, double z)
{
   const double u = x/z, v = y/z;
   const double lambda = sqrtabs(sqr(1 - u - v) - 4*u*v);
   const double xp = 0.5 * (1 + (u - v) - lambda);
   const double xm = 0.5 * (1 - (u - v) - lambda);

   return 1./lambda * (2*logabs(xp)*logabs(xm) - logabs(u)*logabs(v) -
                       2*(li2(xp) + li2(xm)) + M_PI*M_PI/3.);
}

/// First derivative of phi[t,T,g] w.r.t. T
double dphi_010(double t, double T, double g)
{
   constexpr double Pi2 = M_PI * M_PI;
   const double g2 = sqr(g);
   const double abbr = (-4*t*T)/g2 + sqr(1 - t/g - T/g);
   const double rabbr = sqrtabs(abbr);

   return ((g + t - T)*(Pi2 - 6*li2((g - rabbr*g + t - T)/(2.*g)) -
      6*li2((g - rabbr*g - t + T)/(2.*g)) -
      3*logabs(t/g)*logabs(T/g) + 6*logabs((g - rabbr*g + t -
      T)/(2.*g))*logabs((g - rabbr*g - t + T)/(2.*g))) + (3*rabbr*g* (
      rabbr*g*((-1 + rabbr)*g + t - T)*logabs(t/g) +
      2*T*(-2*g*logabs(4.) + (g + rabbr*g + t - T)*logabs((g - rabbr*g
      + t - T)/g) + (g + rabbr*g + t - T)*logabs((g + rabbr*g + t -
      T)/g) + g*logabs((g - rabbr*g - t + T)/g) - rabbr*g*logabs((g -
      rabbr*g - t + T)/g) - t*logabs((g - rabbr*g - t + T)/g) +
      T*logabs((g - rabbr*g - t + T)/g) + g*logabs((g + rabbr*g - t +
      T)/g) - rabbr*g*logabs((g + rabbr*g - t + T)/g) - t*logabs((g +
      rabbr*g - t + T)/g) + T*logabs((g + rabbr*g - t + T)/g)) ) ) /
      (T*(g - rabbr*g - t + T)))/(3.*std::pow(std::abs(abbr),1.5)*g2);
}

} // anonymous namespace

double delta_ma2_2loop_at_as_mst1_eq_mst2(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double g3)
{
   constexpr double Pi2 = M_PI * M_PI;
   const double g = sqr(mg);
   const double g2 = sqr(g);
   const double q = scale2;
   const double q2 = sqr(scale2);
   const double t = mt2;
   const double T = mst12;
   const double sb = std::sin(std::atan(tanb));
   const double ht2 = 2./vev2*mt2/sqr(sb);
   const double At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

   const double result = (-2*(g*(2*At*g + 2*At*t - At*T + mg*T + mg*(g
      - t)*logabs(g/t) - At*T*logabs(g/q)*logabs(t/q) -
      mg*T*logabs(g/q)*logabs(t/q) - 4*mg*T*logabs(T/q) -
      2*At*T*sqr(logabs(T/q)) + logabs((g*t)/q2)*(-(At*(g + t - T)) +
      mg*T + (At + mg)*T*logabs(T/q))) - 2*(At + mg)*(g + t -
      T)*T*phi(t,T,g) + T*(At*(g2 + sqr(t - T) - 2*g*T) + mg*(g2 +
      sqr(t - T) - 2*g*(t + T)))*dphi_010(t,T,g)))/ (g*T);

   const double pref = 4*sqr(g3)/sqr(16*Pi2) * ht2*mu*(1./tanb + tanb);

   return pref * result;
}

double delta_ma2_2loop_at_as_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double g3)
{
   double result = 0.;

   dszodd_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
           &tanb, &vev2, &g3, &result);

   return result;
}

/// 2-loop contribution to CP-odd Higgs O(at*as)
double delta_ma2_2loop_at_as(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double g3)
{
   if (std::abs((mst12 - mst22)/mst12) < 1e-8) {
      const double At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

      // if At = 0 => mu = 0 => dMA(2L) = 0
      if (std::abs(At) < std::numeric_limits<double>::epsilon())
         return 0.;

      return delta_ma2_2loop_at_as_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3);
   }

   return delta_ma2_2loop_at_as_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3);
}

double delta_ma2_2loop_at_at(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2, int atat)
{
   double result = 0.;

   {
      std::lock_guard<std::mutex> lg(mtx);

      ddsodd_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2, &result, &atat);
   }

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_as(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double g3,
   int include_heavy_higgs)
{
   Eigen::Matrix<double, 2, 2> result;

   if (std::abs((mst12 - mst22)/mst12) < 1e-8) {
      result = delta_mh2_2loop_at_as_st_0_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3, 0);
   } else {
      result = delta_mh2_2loop_at_as_general(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3, 0);
   }

   const double dMA = include_heavy_higgs
      ? delta_ma2_2loop_at_as(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, g3)
      : 0.;

   return result + rotate_by(dMA, tanb);
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_at(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2,
   int include_heavy_higgs, int atat)
{
   Eigen::Matrix<double, 2, 2> result;

   {
      std::lock_guard<std::mutex> lg(mtx);

      ddshiggs_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
                &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
                &result(0,0), &result(0,1), &result(1,1), &atat);
   }

   result(1,0) = result(0,1);

   const double dMA = include_heavy_higgs
      ? delta_ma2_2loop_at_at(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2, atat)
      : 0.;

   return result + rotate_by(dMA, tanb);
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_as(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double g3,
   int include_heavy_higgs)
{
   Eigen::Matrix<double, 2, 2> result(delta_mh2_2loop_at_as(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, g3, include_heavy_higgs));

   std::swap(result(0,0), result(1,1));

   return result;
}

double delta_ma2_2loop_atau_atau(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   double result = 0.;

   tausqodd_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result);

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_atau_atau(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2,
   int include_heavy_higgs)
{
   int scheme = 0;
   Eigen::Matrix<double, 2, 2> result;

   tausqhiggs_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
               &costau, &scale2, &mu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   const double dMA = include_heavy_higgs
      ? delta_ma2_2loop_atau_atau(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2)
      : 0.;

   return result + rotate_by(dMA, tanb);
}

double delta_ma2_2loop_ab_atau(
   double mtau2, double mb2,
   double mstau12, double mstau22, double msb12, double msb22,
   double sintau, double costau, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   double result = 0.;

   taubotodd_(&mtau2, &mb2,
              &mstau12, &mstau22, &msb12, &msb22,
              &sintau, &costau, &sxb, &cxb,
              &scale2, &mu, &tanb, &vev2, &result);

   return result;
}

Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_atau(
   double mtau2, double mb2,
   double mstau12, double mstau22, double msb12, double msb22,
   double sintau, double costau, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2,
   int include_heavy_higgs)
{
   Eigen::Matrix<double, 2, 2> result;

   taubot_(&mtau2, &mb2,
           &mstau12, &mstau22, &msb12, &msb22,
           &sintau, &costau, &sxb, &cxb,
           &scale2, &mu, &tanb, &vev2,
           &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   const double dMA = include_heavy_higgs
      ? delta_ma2_2loop_ab_atau(
         mtau2, mb2, mstau12, mstau22, msb12, msb22, sintau, costau, sxb, cxb,
         scale2, mu, tanb, vev2)
      : 0.;

   return result + rotate_by(dMA, tanb);
}

} // namespace mh2l
} // namespace himalaya
